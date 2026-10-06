module setup
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays

export setup_fevalues, setup_dofs, setup_constraints



function setup_fevalues(ipu, ipp, ipg)
    qr = QuadratureRule{RefTriangle}(2)
    cvu = CellValues(qr, ipu, ipg)
    cvp = CellValues(qr, ipp, ipg; update_gradients = false, update_detJdV = false)
    return cvu, cvp
end

function setup_dofs(grid, ipu, ipp)
    dh = DofHandler(grid)
    add!(dh, :u, ipu)
    add!(dh, :p, ipp)
    close!(dh)
    return dh
end
#Busca el centro del facetset (para buscar el centro del eje.)
function boundary_center(grid::Ferrite.AbstractGrid, setname::String)
    ids = Set{Int}()
    for fi in getfacetset(grid, setname)
        cellid, localfacet = fi
        cell = getcells(grid, cellid)
        union!(ids, Ferrite.facets(cell)[localfacet])
    end
    coords = [grid.nodes[i].x for i in ids]
    return sum(coords) / length(coords)
end



function setup_constraints(dh, angle, height, vel_buffer)
    ch = ConstraintHandler(dh)
    #Inlet Constraints
    θ = deg2rad(angle)
    U_max1 = 7.5 #u/s
    U_max2 = vel_buffer #u/s
    c_in1 = boundary_center(dh.grid, "IN1")
    c_in2 = boundary_center(dh.grid, "IN2")

    function parabolic_velocity(x,t, U_max, height, θ, center)
        dx = x[1] - center[1]
        dy = x[2] - center[2]
        V_x = U_max * (1-(-dx*sin(θ)+dy*cos(θ))^2/((height/2)^2))*cos(θ)
        V_y = U_max * (1-(-dx*sin(θ)+dy*cos(θ))^2/((height/2)^2))*sin(θ)
        return Ferrite.Vec{2}((V_x, V_y))
    end
    input1 = Dirichlet(:u, getfacetset(dh.grid, "IN1"), (x,t) -> parabolic_velocity(x,t,U_max1, height,-θ, c_in1), [1,2])
    input2 = Dirichlet(:u, getfacetset(dh.grid, "IN2"), (x,t) -> parabolic_velocity(x,t,U_max2, height, θ, c_in2), [1,2])

    add!(ch, input1)
    add!(ch, input2)
    #Dirichlet Wall Constraints
    CH_global = union(getfacetset(dh.grid, "CH1"),getfacetset(dh.grid, "CH2"), getfacetset(dh.grid, "CH3"), getfacetset(dh.grid, "CH4"))
    wall_dbc = Dirichlet(:u, CH_global, (x,t)->[0.0,0.0], [1,2])
    add!(ch, wall_dbc)
    #Output Constraints
    OUT_global = union(getfacetset(dh.grid, "OUT1"), getfacetset(dh.grid, "OUT2"))
    output_dbc = Dirichlet(:p, OUT_global, (x,t)-> 0.0)
    add!(ch, output_dbc)

    #finalize
    close!(ch)
    Ferrite.update!(ch,0)

    return ch

end

end