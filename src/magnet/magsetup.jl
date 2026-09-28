module magsetup
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays

export magsetup_fevalues, magsetup_dofs, magsetup_const


function magsetup_fevalues(ipu, ipg)
    qr = QuadratureRule{RefTriangle}(2)
    cvu = CellValues(qr, ipu, ipg)
    return cvu, qr
end

#degrees of fredom
function magsetup_dofs(grid, ipu)
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

    close!(ch)
    return ch
end


end