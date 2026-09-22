module magssemble
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays

export magssemble_system!

    #esto es lo chido hay que reafirmar que las ecuaciones son correctas
    #ahorita lo checamos
    function magssemble_system!(k,f,dh,cvu, magrid)


        μ₀ = (4*pi)*1e-7
        μᵣ = 1.05 #Permeabilidad relativa neodimio

        Bᵣ₁ =  Ferrite.Vec{2,Float64}((0, 1.2)) 
        Bᵣ₂ = Ferrite.Vec{2,Float64}((0, 0))

        magassembler = start_assemble(k,f)
        ke = zeros(ndofs_per_cell(dh),ndofs_per_cell(dh))
        fe = zeros(ndofs_per_cell(dh))

        range_u = dof_range(dh,:u)
        ndofs_u = length(range_u)

        total_cells = getncells(magrid)
        println("TOTAL CELLS> ", total_cells)
        ∇ϕ = Vector{Ferrite.Vec{2,Float64}}(undef, ndofs_u)

        for (i,cell) in enumerate(CellIterator(dh))
            reinit!(cvu,cell)
            domaincell = getcells(magrid, "AIR")

            if (i % 50000 == 0) || (i == total_cells)
                println("Cell n> ", i)
                println("Progreso de assembly: ", i/total_cells * 100, "%")
            end

            if cell ∈ domaincell

                μ = μ₀
                Bᵣ = Bᵣ₁  # teslas dominio  
                
            else
                μ = μᵣ*μ₀
                Bᵣ = Bᵣ₂  # teslas dominio  
            end

            ke .= 0
            fe .= 0
            for qp ∈ 1:getnquadpoints(cvu)
                dΩ = getdetJdV(cvu, qp)
                for i in 1:ndofs_u
                ∇ϕ[i] = shape_gradient(cvu, qp, i)
                end
                for (i,I) in pairs(range_u), (j,J) in pairs(range_u)
                ke[I,J] += μ * (∇ϕ[i] ⋅ ∇ϕ[j]) * dΩ
                end
                for (i,I) in pairs(range_u)
                fe[I] += (μ/μ₀) * ∇ϕ[i] ⋅ Bᵣ * dΩ
                end

            end
            assemble!(magassembler,celldofs(cell), ke, fe)
        end
        return k, f
    end


end
