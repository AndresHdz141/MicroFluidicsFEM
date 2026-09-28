module PostProcess

using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays
export get_H_field

function get_H_field(grid, dh, cvu,ipu, qr, u)
    dim = Ferrite.getspatialdim(grid)
    grad_u = [[zero(Vec{dim, Float64}) for _ in 1:getnquadpoints(cvu)] for _ in 1:getncells(grid)]

    #grad_H² = [[zero(Vec{dim, Float64}) for _ in 1:getnquadpoints(cvu)] for _ in 1:getncells(grid)]

    for (idx, cell) in enumerate(CellIterator(dh))
        reinit!(cvu, cell)

        potential = grad_u[cellid(cell)]

        for qp in 1:getnquadpoints(cvu)
            #dΩ = getdetJdV(cvu, qp)
            ∇u = function_gradient(cvu, qp, u, celldofs(cell))
            potential[qp] = -∇u     

        end

    end

    projector = L2Projector(grid)
    add!(projector, 1:getncells(grid), ipu, qr_rhs=qr)
    close!(projector)
    H_field = project(projector, grad_u)
    return H_field, projector
end

end