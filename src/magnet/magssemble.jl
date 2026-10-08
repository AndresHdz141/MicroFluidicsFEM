module magssemble
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays

export magssemble_system!

function magssemble_system!(k, f, dh, cvu, array, magrid)

    μ₀ = 4π * 1e-7
    μᵣ = 1.05

    Bᵣ₁x = Ferrite.Vec{2,Float64}((1.2, 0.0))
    Bᵣ₁y = Ferrite.Vec{2,Float64}((0.0, 1.2))
    Bᵣ₀  = Ferrite.Vec{2,Float64}((0.0, 0.0))

    # Remanencia por grupo. Un grupo ausente (p. ej. magnets_3 con 2 imanes) se ignora.
    getset(name) = haskey(magrid.cellsets, name) ? getcellset(magrid, name) : Set{Int}()

    groups = if array == 1
        ["magnets_par" => Bᵣ₁y, "magnets_impar" => -Bᵣ₁y]
    elseif array == 2
        ["magnets_1" => Bᵣ₁y, "magnets_2" => -Bᵣ₁x,
         "magnets_3" => -Bᵣ₁y, "magnets_4" => Bᵣ₁x]
    else
        ["magnets" => Bᵣ₁y]
    end

    # cell id => Bᵣ
    Br_of_cell = Dict{Int,Ferrite.Vec{2,Float64}}()
    for (name, Br) in groups
        for c in getset(name)
            Br_of_cell[c] = Br
        end
    end

    total_cells = getncells(magrid)
    println("TOTAL CELLS> ", total_cells)
    println("CELDAS DE IMAN> ", length(Br_of_cell), " de ", total_cells)
    @assert !isempty(Br_of_cell) "Sin celdas de iman. cellsets disponibles: $(collect(keys(magrid.cellsets)))"

    magassembler = start_assemble(k, f)
    nd = ndofs_per_cell(dh)
    ke = zeros(nd, nd)
    fe = zeros(nd)

    range_u = dof_range(dh, :u)
    ndofs_u = length(range_u)
    ∇ϕ = Vector{Ferrite.Vec{2,Float64}}(undef, ndofs_u)

    step = max(1, total_cells ÷ 10)

    for cell in CellIterator(dh)
        cid = cellid(cell)
        reinit!(cvu, cell)

        if cid % step == 0 || cid == total_cells
            println("Cell n> ", cid)
            println("Progreso de assembly: ", round(cid / total_cells * 100, digits=1), "%")
        end

        Bᵣ = get(Br_of_cell, cid, Bᵣ₀)
        μ  = haskey(Br_of_cell, cid) ? μᵣ * μ₀ : μ₀
        Mᵣ = Bᵣ / μ₀

        fill!(ke, 0.0)
        fill!(fe, 0.0)

        for qp in 1:getnquadpoints(cvu)
            dΩ = getdetJdV(cvu, qp)
            for a in 1:ndofs_u
                ∇ϕ[a] = shape_gradient(cvu, qp, a)
            end
            for (a, I) in pairs(range_u), (b, J) in pairs(range_u)
                ke[I, J] += μ * (∇ϕ[a] ⋅ ∇ϕ[b]) * dΩ
            end
            for (a, I) in pairs(range_u)
                fe[I] += μ₀ * (∇ϕ[a] ⋅ Mᵣ) * dΩ
            end
        end

        assemble!(magassembler, celldofs(cell), ke, fe)
    end
    return k, f
end

end