module simulate_magnet
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays, GLMakie, JLD2

include("../src/magnet/magnet_grid.jl")
include("../src/magnet/magsetup.jl")
include("../src/magnet/magssemble.jl")
include("../src/magnet/postprocess.jl")


using .grid_mag
using .magssemble
using .magsetup
using .PostProcess

export mag_sim

    function mag_sim(length, height,n, arreglo,nombre,dist_iman, h = 0.1)

        grid = magnet_setup_grid(length, height,n,dist_iman, arreglo, h)
        println("Grid listo")
        #Interpolation
        ipu = Lagrange{RefTriangle, 2}()
        println("Interpolacion listo")

        #dofs
        dh = magsetup_dofs(grid, ipu)
        println("DOFS listo")

        #FE Values
        ipg = Lagrange{RefTriangle, 1}()
        cvu, qr= magsetup_fevalues(ipu, ipg)
        println("FE Values listo")

        #Boundary conditions
        ch = magsetup_const(dh)
        println("Constraints listo")

        K = allocate_matrix(dh, ch)
        f = zeros(ndofs(dh))
        println("Coupling listo")

        magssemble_system!(K,f,dh,cvu,arreglo, grid)
        println("Assemble listo")

        apply!(K,f,ch)
        u = K \ f
        apply!(u,ch)
        println("Solucion lista")

        H,H_mag, proj = get_H_field(grid, dh, cvu,ipu, qr, u)
        println("H field calculado")

        ∇H², proj2 = get_grad_magH(grid,cvu, ipu,qr,H_mag)
        println("∇H² calculado")

        nombre_data = nombre * ".jld2"
        location_vis = joinpath("data", "visual_imanes", nombre)
        location_data = joinpath("data", "arreglos_imanes",nombre_data)
        VTKGridFile(location_vis, grid) do vtk
        write_solution(vtk, dh, u, :Potential)
        write_projection(vtk, proj, H, "H_Field")   
        write_projection(vtk, proj, H_mag, "Mag_H^2")
        write_projection(vtk, proj2,∇H², "grad_H2" )

        jldsave(location_data;
            mag_name = "Magnetic Data",
            mag_mesh = grid, 
            mag_projector = proj2,
            mag_H2_field = ∇H²
        )

        end
        return
    end

end