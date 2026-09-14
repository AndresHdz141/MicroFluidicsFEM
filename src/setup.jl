module setup
export setup_fevalues



function setup_fevalues(ipu, ipp, ipg)
    qr = QuadratureRule{RefTriangle}(2)
    cvu = CellValues(qr, ipu, ipg)
    cvp = CellValues(qr, ipp, ipg; update_gradients = false, update_detJdV = false)
    qr_facet = FacetQuadratureRule{RefTriangle}(2)
    fvp = FacetValues(qr_facet, ipp, ipg; update_gradients = false)
    return cvu, cvp, fvp
end

function setup_dofs(grid, ipu, ipp)
    dh = DofHandler(grid)
    add!(dh, :u, ipu)
    add!(dh, :p, ipp)
    close!(dh)
    return dh
end
#Busca el centro del facetset (para buscar el centro del eje.)
function get_facetset_center(grid, set_name::String)
    facet_set = getfacetset(grid, set_name)
    total_coord = Vec{2}((0.0,0.0))
    count = 0 
    for (cell_idx, facet_idx) in facet_set
        facet_coords = getcoordinates(grid, cell_idx)[getnodes(grid.cells[cell_idx][facet_idx])]
        for pt in facet_coords
            total_coord += pt 
            count+=1
        end
    end
    return total_coord/count
end


function setup_constraints(dh, angle, height)
    ch = ConstraintHandler(dh)
    #Inlet Constraints
    θ = deg2rad(angle)
    U_max1 = 75.0 #μm/s
    U_max2 = 120.0 #μm/s
    c_in1 = get_facetset_center(dh.grid, "IN1")
    c_in2 = get_facetset_center(dh.grid, "IN2")


    function parabolic_velocity(x,t, U_max, height, θ, xc, yc)
        dx = x[1] - xc
        dy = x[2] - yc
        V_x = U_max * (1-(-dx*sin(θ)+dy*cos(θ))^2/((height/2)^2))*cos(θ)
        V_y = U_max * (1-(-dx*sin(θ)+dy*cos(θ))^2/((height/2)^2))*sin(θ)
        return Vec{2}((V_x, V_y))
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
    update!(ch,0)

    return ch

end

end