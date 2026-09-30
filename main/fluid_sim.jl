using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays, JLD2


include("../src/fluid/grid.jl")
include("../src/fluid/setup.jl")
include("../src/fluid/assemble.jl")

using .grid
using .assemble
using .setup



function main(length, height,n_ranuras, angle_d,h=0.05)

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
    ch = setup_constraints(dh, angle_d, height)
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

    VTKGridFile("data/stokes-flow", grid) do vtk
        write_solution(vtk, dh, u)

    jldsave("data/fluid_data.jld2";
        f_name = "Fluid data",
        f_mesh = grid, 
        f_dh = dh,
        f_u = u
    )
    end
    return
end

main(24.0,4.0,5,45.0)
