module evaluate_fields
using StaticArrays, JLD2, Ferrite
export get_field_value
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


#Funcion para obtener el valor del campo en un punto en especifico
function get_field_value(data,point, isFluid)
    #Toma un vector SVector y lo convierte a Ferrite.Vec
    fvec = Vec{2}(Tuple(point)) 
    #Calcula el PointHandler para el punto en especifico
    ph = PointEvalHandler(data[1], [fvec])
    #Obtiene el valor del vector dependiendo del tipo de campo a evaluar (fluido o magnetico)
    value = isFluid ? evaluate_at_points(ph, data[2], data[3], :u ) :  evaluate_at_points(ph, data[2], data[3])
    #Convierte el valor obtenido el cual es Ferrite.Vec a un SVector
    svec = SVector{2}(value[1])
    return svec
end

#= EJEMPLO
        point = SVector{2, Float64}(0.0,1.0)
        mag_data = [mag_grid, mag_dh, mag_field]
        field_value = get_field_value(mag_data, point)

        println(field_value, "is a type:", typeof(field_value))
=#

end