module evaluate_fields
using StaticArrays, JLD2, Ferrite, NearestNeighbors
export setup_fast_eval, fast_evaluate_field
#=
----Asi se lee el archivo de fluido-----

file = jldopen("data/fluid_data.jld2", "r")
f_name =  file["f_name"]
f_grid = file["f_mesh"]
f_dh = file["f_dh"]
f_field = file["f_u"] 
close(file)

-----Asi se lee el archivo de campo magnetico----
file = jldopen("data/magnetic_data.jld2", "r")
mag_name =  file["mag_name"]
mag_grid = file["mag_mesh"]
mag_dh = file["mag_projector"]
mag_field = file["mag_H2_field"] 
close(file)
=#


# Obtener arbol de datos
function setup_fast_eval(grid, dh,field_name::Symbol=:_)
    #Obtener coordenadas de centros de celdas
    centroids = map(1:getncells(grid)) do i
        coords = getcoordinates(grid, i)
        center = sum(coords)/length(coords)
        SVector{2, Float64}(center[1], center[2])
    end
    #generar arbol de datos, interpolacion y pointValues
    tree = KDTree(centroids)
    field_idx = Ferrite.find_field(dh, field_name)
    ip = Ferrite.getfieldinterpolation(dh, field_idx)
    pv = PointValues(Float64, ip)
    return [tree, pv]
end

function fast_evaluate_field(pos::SVector{2, Float64}, tree::KDTree, grid, dh, field_data, pv, field_name::Symbol = :_)
    #Se obtienen las 3 celdas mas cercanas al punto/posicion
    nearest_cell_ids, _ = knn(tree, pos, 3)

    #Para cada celda se revisa si la particula esta dentro de esa celda, para eso>
    for cell_id in nearest_cell_ids
        #Se obtienen las coordenadas de los nodos de la celda
        cell_coords = getcoordinates(grid, cell_id)
        #Las coordenadas se guardan en x1, x2, x3
        x1,x2,x3 = cell_coords[1], cell_coords[2], cell_coords[3]
        #=
        Se tienen los nodos de los triangulos en las coordenadas globales,
         pero se quieren transformar al sistema de coordenadas locales del triangulo
         (0,0),(0,1),(1,0)

         (0,1)
           |  \
           |   \
           |     \
         (0,0)---(1,0)

         En donde un punto (x,y) se puede expresar en estas coordenadas al realizar una transormacion
         En donde:
         x(ξ, η) = x1 + (x2-x1)ξ + (x3-x1)η
         En donde para resolver por ξ y η
         x-x1 = (x2-x1)ξ + (x3-x1)η
         => 
         [(x2-x1) (x3-x1) ; (y2-y1) (y3-y1)] [ξ ; η] = [x-x1 ; y - y1]
        En donde se usa la regla de Cramer para resolver por  [ξ ; η], siendo>
        ξ = (J₂₂(x-x1)-J₁₂(y-y1)) / (det(J))
        η = (-J₂₁(x-x1)+J₁₁(y-y1)) / (det(J))

        Luego se revisa si el punto (ξ, η) existe dentro de la celda, en donde se debe cumplir:
        1. ξ ≥ 0    -> Esta a la derecha del punto vertical
        2. η ≥ 0    -> Esta sobre el punto horizontal
        3. η+ξ ≤ 1  -> Estan debajo de la hipotenusa del triangulo

        Si se cumple, se realiza una interpolacion en el punto
        Lo explicado se observa en el siguiente segmento ====>
        =#
        J11 = x2[1] - x1[1] ; J12 = x3[1] - x1[1]
        J21 = x2[2] - x1[2] ; J22 = x3[2] - x1[2]
        detJ = J11*J22 - J12 * J21
        b1 = pos[1] - x1[1]
        b2 = pos[2] - x1[2]

        xi_1 = (J22*b1 - J12 * b2) / detJ
        xi_2 = (-J21*b1 + J11 * b2) / detJ
        xi = Ferrite.Vec{2, Float64}((xi_1, xi_2))

        tol = 1e-3
        if xi_1 >= -tol && xi_2 >= -tol && (xi_1 + xi_2) <= 1.0 + tol
            reinit!(pv, cell_coords, xi)
            cell_dofs = celldofs(dh, cell_id)
            u_local = field_data[cell_dofs]
            dr = dof_range(dh, field_name)
            value_field = function_value(pv, u_local, dr)
            SFIELD_VALUE = SVector{2}(value_field[1], value_field[2])

            return SFIELD_VALUE
        end
    end
    #println("PELIGRO:La particula en pos = $pos no se encuentra dentro de ninguna celda !")
    return SVector{2,Float64}(0.0,0.0)
end

end