module assemble
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays

export assemble_system!

function assemble_system!(K,f,dh,cvu,cvp, height)

    assembler = start_assemble(K,f)
    ke = zeros(ndofs_per_cell(dh),ndofs_per_cell(dh) )
    fe = zeros(ndofs_per_cell(dh))

    range_u = dof_range(dh, :u)
    ndofs_u = length(range_u)
    range_p = dof_range(dh,:p)
    ndofs_p = length(range_p)

    Φᵤ=Vector{Ferrite.Vec{2,Float64}}(undef, ndofs_u)
    ∇Φᵤ = Vector{Tensor{2,2,Float64,4}}(undef, ndofs_u)
    divΦᵤ=Vector{Float64}(undef, ndofs_u)
    Φₚ = Vector{Float64}(undef, ndofs_p)


    for cell in CellIterator(dh)
        reinit!(cvu, cell)
        reinit!(cvp, cell)
        coords = getcoordinates(cell)
        
        y_coords = [c[2] for c in coords]
        y_center = sum(y_coords) / length(y_coords)

        if y_center >=height/2
            μ = 0.005
        else
            μ = 0.007
        end
        ke .= 0
        fe .= 0
        for qp in 1:getnquadpoints(cvu)
            dΩ = getdetJdV(cvu, qp)
            for i in 1:ndofs_u
                Φᵤ[i] = shape_value(cvu, qp, i)
                ∇Φᵤ[i] = shape_gradient(cvu, qp, i)
                divΦᵤ[i] = shape_divergence(cvu, qp, i)
            end
            for i in 1:ndofs_p
                Φₚ[i] = shape_value(cvp, qp, i)
            end
            #u-u
            for(i, I) in pairs(range_u), (j,J) in pairs(range_u)
                ke[I,J] += μ*(∇Φᵤ[i] ⊡ ∇Φᵤ[j]) * dΩ
            end
            #u-p
            for(i, I) in pairs(range_u), (j,J) in pairs(range_p)
            ke[I,J] += (-divΦᵤ[i] * Φₚ[j]) * dΩ
            end
            #p-u
            for(i, I) in pairs(range_p), (j,J) in pairs(range_u)
            ke[I,J] += (-divΦᵤ[j] * Φₚ[i]) * dΩ
            end
        end
        assemble!(assembler, celldofs(cell), ke, fe)
    end
    return K,f
end





end