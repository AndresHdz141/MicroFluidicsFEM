module magssemble
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays

export magssemble_system!

    #esto es lo chido hay que reafirmar que las ecuaciones son correctas
    #ahorita lo checamos
    function magssemble_system!(k,f,dh,cvu, array, magrid)


        μ₀ = (4*pi)*1e-7
        μᵣ = 1.05 #Permeabilidad relativa neodimio

        Bᵣ₁x =  Ferrite.Vec{2,Float64}((1.2,0.0)) 
        Bᵣ₁y =  Ferrite.Vec{2,Float64}((0.0,1.2)) 
        Bᵣ₂ = Ferrite.Vec{2,Float64}((0.0, 0.0))

        magassembler = start_assemble(k,f)
        ke = zeros(ndofs_per_cell(dh),ndofs_per_cell(dh))
        fe = zeros(ndofs_per_cell(dh))

        range_u = dof_range(dh,:u)
        ndofs_u = length(range_u)

        total_cells = getncells(magrid)
        println("TOTAL CELLS> ", total_cells)
        ∇ϕ = Vector{Ferrite.Vec{2,Float64}}(undef, ndofs_u)
        
        if array == 1
        magnets_imp = getcellset(magrid, "magnets_impar")
        magnets_par = getcellset(magrid, "magnets_par")
        magnet_domain = union(magnets_imp, magnets_par)

        elseif array ==2
        magnets_1 = getcellset(magrid, "magnets_1")
        magnets_2 = getcellset(magrid, "magnets_2")
        magnets_3 = getcellset(magrid, "magnets_3")
        magnets_4 = getcellset(magrid, "magnets_4")
        magnet_domain = union(magnets_1, magnets_2, magnets_3, magnets_4)
        else
            magnets = getcellset(magrid, "magnets")
            
        end
    

        for (i,cell) in enumerate(CellIterator(dh))
            reinit!(cvu,cell)

            if (i % (total_cells/10) == 0) || (i == total_cells)
                println("Cell n> ", i)
                println("Progreso de assembly: ", i/total_cells * 100, "%")
            end
            
            ke .= 0
            fe .= 0
            for qp ∈ 1:getnquadpoints(cvu)

                if array == 1
                    if i in magnet_domain
                        μ = μᵣ*μ₀
                        Bᵣ = i in magnets_par ? Bᵣ₁y : -Bᵣ₁y
                    else
                        μ= μ₀
                        Bᵣ= Bᵣ₂
                    end #:)
                elseif array == 2
                    if i in magnet_domain
                        μ = μᵣ*μ₀
                        Bᵣ = i in magnets_1 ? Bᵣ₁y : i in magnets_2 ? -Bᵣ₁x : i in magnets_3 ? -Bᵣ₁y : Bᵣ₁x
                    else
                        μ= μ₀
                        Bᵣ= Bᵣ₂
                    end
                else
                    μ = i in magnets ? μᵣ*μ₀ : μ₀
                    Bᵣ = i in magnets ? Bᵣ₁y : Bᵣ₂
                end

            Mᵣ = Bᵣ / μ₀

                dΩ = getdetJdV(cvu, qp)
                for i in 1:ndofs_u
                ∇ϕ[i] = shape_gradient(cvu, qp, i)
                end
                for (i,I) in pairs(range_u), (j,J) in pairs(range_u)
                ke[I,J] += μ * (∇ϕ[i] ⋅ ∇ϕ[j]) * dΩ
                end
                for (i,I) in pairs(range_u)
                fe[I] += μ₀ * (∇ϕ[i] ⋅ Mᵣ) * dΩ
                end

            end
            assemble!(magassembler,celldofs(cell), ke, fe)
        end
        return k, f
    end


end
