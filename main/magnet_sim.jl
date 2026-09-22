using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays, GLMakie

include("../src/magnet/magnet_grid.jl")
include("../src/magnet/magsetup.jl")
include("../src/magnet/magssemble.jl")

using .grid_mag
using .magssemble
using .magsetup


function main(length, height, h = 0.05)

    grid = magnet_setup_grid(length, height, h)
    println("Grid listo")
    #Interpolation
    ipu = Lagrange{RefTriangle, 2}()
    println("Interpolacion listo")

    #dofs
    dh = magsetup_dofs(grid, ipu)
    println("DOFS listo")

    #FE Values
    ipg = Lagrange{RefTriangle, 1}()
    cvu= magsetup_fevalues(ipu, ipg)
    println("FE Values listo")

    #Boundary conditions
    ch = magsetup_const(dh)
    println("Constraints listo")

    K = allocate_matrix(dh, ch)
    f = zeros(ndofs(dh))
    println("Coupling listo")

    magssemble_system!(K,f,dh,cvu, grid)
    println("Assemble listo")

    apply!(K,f,ch)
    u = K \ f
    apply!(u,ch)
    println("Solucion lista")

    VTKGridFile("Scalar_field2", grid) do vtk
    write_solution(vtk, dh, u)
    end
    return
end

main(24.0, 4.0, 0.4)