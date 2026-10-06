module simulate_fluid
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays, JLD2


include("../src/fluid/grid.jl")
include("../src/fluid/setup.jl")
include("../src/fluid/assemble.jl")

using .grid
using .assemble
using .setup

export fluid_sim

    function fluid_sim(length, height,vel_buffer,n_ranuras,nombre, angle_d,h=0.05)

        #setup grid
        grid = setup_grid(length, height, angle_d,n_ranuras, h)
        println("Grid listo")
        #Interpolation
        ipu = Lagrange{RefTriangle, 2}()^2
        ipp = Lagrange{RefTriangle, 1}()
        println("Interpolacion listo")

        #dofs
        dh = setup_dofs(grid, ipu, ipp)
        println("DOFS listo")

        #FE Values
        ipg = Lagrange{RefTriangle, 1}()
        cvu, cvp = setup_fevalues(ipu, ipp, ipg)
        println("FE Values listo")

        #Boundary conditions
        ch = setup_constraints(dh, angle_d, height, vel_buffer)
        println("Constraints listo")

        coupling = [true true; true false]
        K = allocate_matrix(dh, ch;coupling=coupling)
        f = zeros(ndofs(dh))
        println("Coupling listo")

        assemble_system!(K,f,dh,cvu,cvp, height)
        println("Assemble listo")

        apply!(K,f,ch)
        u = K \ f
        apply!(u,ch)
        println("Solucion lista")

        nombre_data = nombre * ".jld2"
        location_vis = joinpath("data", "visual_fluido", nombre)
        location_data = joinpath("data", "campo_fluido", nombre_data)
        VTKGridFile(location_vis, grid) do vtk
            write_solution(vtk, dh, u)

        jldsave(location_data;
            f_name = "Fluid data",
            f_mesh = grid, 
            f_dh = dh,
            f_u = u
        )
        end
        return
    end

end