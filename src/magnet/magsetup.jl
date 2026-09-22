module magsetup
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays

export magsetup_fevalues, magsetup_dofs, magsetup_constraints


function magsetup_fevalues(ipu, ipg)
    qr = QuadratureRule{RefTriangle}(2)
    cvu = CellValues(qr, ipu, ipg)
    
    return qr, cvu
end

#degrees of fredom
function magsetup_dofs(grid, ipu, ipp)
    dh = DofHandler(grid)
    add!(dh, :u, ipu)
    close!(dh)
    return dh
end


#restriction

function magsetup_const(dh)
    ch = ConstraintHandler(dh)

    #bordes del aire osea que el potencial es 0en ese punto
    AirLm_C = Dirichlet(:u, getfacetset(dh.grid, "AIR_Lm"), (x,t) -> 0.0)
    add!(ch, AirLm_C)

        #aire en todo lo demas creoque  esto node

    #=Air_C = Dirichlet(:u, getfacet(dh.grid, "Air"),(x,t) -> [])


    #condición de magnetización dentro del iman en dirección de y
    Magnet_c = Dirichlet(:u, getfaceset(dh.gdid, "magnets"), (x,t) -> [0,Mᵣ], [1,2]) # creo que Mᵣ
    add!(ch, Magnet_c)=#

    close!(ch)
    return ch
end


end