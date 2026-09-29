module PostProcess

using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays
export get_H_field, get_grad_magH

function get_H_field(grid, dh, cvu,ipu, qr, u)
    dim = Ferrite.getspatialdim(grid)
    grad_u = [[zero(Vec{dim, Float64}) for _ in 1:getnquadpoints(cvu)] for _ in 1:getncells(grid)]
    H_mag = [[zero(Float64) for _ in 1:getnquadpoints(cvu)] for _ in 1:getncells(grid)]

    #grad_H² = [[zero(Vec{dim, Float64}) for _ in 1:getnquadpoints(cvu)] for _ in 1:getncells(grid)]

    for (idx, cell) in enumerate(CellIterator(dh))
        reinit!(cvu, cell)

        potential = grad_u[cellid(cell)]
        H² = H_mag[cellid(cell)]
        #∇H² = grad_H_mag[cellid(cell)]

        for qp in 1:getnquadpoints(cvu)
            #dΩ = getdetJdV(cvu, qp)
            ∇u = function_gradient(cvu, qp, u, celldofs(cell))
            potential[qp] = -∇u     
            H²[qp] = potential[qp] ⋅ potential[qp]
        end

    end

    projector = L2Projector(grid)
    add!(projector, 1:getncells(grid), ipu, qr_rhs=qr)
    close!(projector)
    H_field = project(projector, grad_u)
    mag_H_field = project(projector, H_mag)
    return H_field, mag_H_field, projector
end

function get_grad_magH(grid,cvu, ipu,qr,H)
    dim = Ferrite.getspatialdim(grid)
    dh_H² = DofHandler(grid)
    add!(dh_H², :H2, ipu)
    close!(dh_H²)

    ∇H2_grid = [[zero(Vec{dim, Float64}) for _ in 1:getnquadpoints(cvu)] for _ in 1:getncells(grid)]

    for cell in CellIterator(dh_H²)
        reinit!(cvu, cell)
        gradient= ∇H2_grid[cellid(cell)] 
        for qp in 1:getnquadpoints(cvu)
            ∇H2_gradient_evaluator = function_gradient(cvu, qp, H, celldofs(cell))
            gradient[qp] = ∇H2_gradient_evaluator
        end

    end
    projector = L2Projector(grid)
    add!(projector, 1:getncells(grid), ipu, qr_rhs = qr)
    close!(projector)
    ∇H² = project(projector, ∇H2_grid)
    return ∇H², projector
end

end