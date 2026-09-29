using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays, GLMakie


include("../src/fluid/grid.jl")
include("../src/fluid/setup.jl")
include("../src/fluid/assemble.jl")

using .grid
using .assemble
using .setup



function main(length, height, angle_d,h=0.01)

    #setup grid
    grid = setup_grid(length, height, angle_d, h)
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

    VTKGridFile("data/stokes-flow3", grid) do vtk
        write_solution(vtk, dh, u)
    end
    return
end

main(24.0,4.0,45.0)
