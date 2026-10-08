module grid_mag
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays

export magnet_setup_grid

function magnet_setup_grid(L, height, n_imanes, dist_iman, array, h=0.05)

    a     = L / 2
    paso  = L / n_imanes
    lado  = 0.6 * paso
    x1    = 2.2 * a
    y1    = a * (3 / 2)
    x_init = -a + (paso - lado) / 2
    y_init = -dist_iman - lado

    Gmsh.initialize()
    gmsh.clear()
    try
        gmsh.option.set_number("General.Verbosity", 2)

        # ---- puntos de los imanes ----
        imanes = Vector{Vector{Int}}()
        for i in 1:n_imanes
            x = x_init + (i - 1) * paso
            y = y_init
            push!(imanes, [
                gmsh.model.geo.add_point(x,        y,        0.0, h),
                gmsh.model.geo.add_point(x + lado, y,        0.0, h),
                gmsh.model.geo.add_point(x + lado, y + lado, 0.0, h),
                gmsh.model.geo.add_point(x,        y + lado, 0.0, h),
            ])
        end

        # ---- dominio de aire ----
        p1 = gmsh.model.geo.add_point( x1,  y1, 0.0, h)
        p2 = gmsh.model.geo.add_point(-x1,  y1, 0.0, h)
        p3 = gmsh.model.geo.add_point(-x1, -y1, 0.0, h)
        p4 = gmsh.model.geo.add_point( x1, -y1, 0.0, h)

        l1 = gmsh.model.geo.add_line(p1, p2)
        l2 = gmsh.model.geo.add_line(p2, p3)
        l3 = gmsh.model.geo.add_line(p3, p4)
        l4 = gmsh.model.geo.add_line(p4, p1)
        println("Check 1 - Grid")

        # ---- lineas y loops de imanes ----
        lopps = Int[]
        for i in 1:n_imanes
            pts = imanes[i]
            lin = [
                gmsh.model.geo.add_line(pts[1], pts[2]),
                gmsh.model.geo.add_line(pts[2], pts[3]),
                gmsh.model.geo.add_line(pts[3], pts[4]),
                gmsh.model.geo.add_line(pts[4], pts[1]),
            ]
            push!(lopps, gmsh.model.geo.add_curve_loop(lin))
        end
        loop = gmsh.model.geo.add_curve_loop([l1, l2, l3, l4])
        println("Check 2 - Grid")

        # ---- superficies ----
        surf = gmsh.model.geo.add_plane_surface([loop, lopps...])

        # grupos de superficies de imanes: nombre => vector de superficies
        magnet_groups = Pair{String,Vector{Int}}[]
        if array == 1
            push!(magnet_groups, "magnets_impar" => Int[])
            push!(magnet_groups, "magnets_par"   => Int[])
            for i in 1:n_imanes
                srf = gmsh.model.geo.add_plane_surface([lopps[i]])
                push!(magnet_groups[iseven(i) ? 2 : 1].second, srf)
            end
        elseif array == 2
            for k in 1:4
                push!(magnet_groups, "magnets_$k" => Int[])
            end
            for i in 1:n_imanes
                srf = gmsh.model.geo.add_plane_surface([lopps[i]])
                push!(magnet_groups[mod1(i, 4)].second, srf)
            end
        else
            push!(magnet_groups, "magnets" => Int[])
            for i in 1:n_imanes
                srf = gmsh.model.geo.add_plane_surface([lopps[i]])
                push!(magnet_groups[1].second, srf)
            end
        end

        # ---- sincronizar ANTES de crear grupos fisicos ----
        gmsh.model.geo.synchronize()

        # ---- grupos fisicos (una sola vez) ----
        gmsh.model.add_physical_group(1, [l1, l2, l3, l4], -1, "AIR_Lm")
        gmsh.model.add_physical_group(2, [surf], -1, "AIR")
        for (name, surfs) in magnet_groups
            isempty(surfs) && continue
            gmsh.model.add_physical_group(2, surfs, -1, name)
        end

        gmsh.model.mesh.generate(2)
        println("Check 3 - Grid")

        grid = mktempdir() do dir
            path = joinpath(dir, "meshmag.msh")
            gmsh.write(path)
            togrid(path)
        end

        println("cellsets del grid magnetico: ", sort(collect(keys(grid.cellsets))))
        return grid
    finally
        Gmsh.finalize()
    end
end

end