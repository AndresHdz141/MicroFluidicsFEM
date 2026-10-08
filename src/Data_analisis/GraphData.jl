module DataGraph

using GLMakie, Ferrite, FerriteGmsh, Random
export PlotTray,AnimateParticles
#Se realizo con ayuda de Claude - no sabia como graficar los meshes

function PlotTray(data, grid, grid2, env, a)
    fig = Figure(size = (1920, 1080))
    ax  = Axis(fig[1, 1]; title = "Trayectorias de particulas",
               xlabel = "X (10 μm)", ylabel = "Y (10 μm)",
               aspect = DataAspect())

    draw_boundaries!(ax, grid)
    draw_magnets!(ax, grid2)
    #hlines!(ax, [a + env.gap_cap]; color = :blue, linestyle = :dash)   # capture line
    for (i,part) in enumerate(data)
        lines!(ax, part.pos, color = :red, alpha = 0.4)
    end
    axislegend(ax)
    display(fig)
end

function facet_segments(grid, setname)
    pts = Point2f[]
    for (cid, lf) in getfacetset(grid, setname)
        nds = Ferrite.facets(getcells(grid, cid))[lf]
        push!(pts, Point2f(grid.nodes[nds[1]].x))   # first two nodes = the edge ends
        push!(pts, Point2f(grid.nodes[nds[2]].x))
    end
    return pts          # consecutive pairs = one segment each
end

function draw_boundaries!(ax, grid; scale = 1.0)
    for (name, _) in grid.facetsets                  # every set in the grid
        pts = facet_segments(grid, name) .* scale
        isempty(pts) && continue
        linesegments!(ax, pts; color = :black,linewidth = 2, label = "")
    end
end

function cellset_outline(grid, setname)
    cnt = Dict{Tuple{Int,Int},Int}()
    for cid in getcellset(grid, setname)
        for f in Ferrite.facets(getcells(grid, cid))
            k = f[1] < f[2] ? (f[1], f[2]) : (f[2], f[1])
            cnt[k] = get(cnt, k, 0) + 1
        end
    end
    pts = Point2f[]
    for (k, c) in cnt
        c == 1 || continue                      # interior edges are shared by 2 cells
        push!(pts, Point2f(grid.nodes[k[1]].x), Point2f(grid.nodes[k[2]].x))
    end
    return pts
end

function draw_magnets!(ax, mag_grid; scale = 1.0)
    for name in keys(mag_grid.cellsets)
        startswith(name, "magnets") || continue
        pts = cellset_outline(mag_grid, name) .* scale
        isempty(pts) && continue
        linesegments!(ax, pts; color = :black, linewidth = 2)
    end
end

function AnimateParticles(data, grid, mag_grid, a, dt; every = 3, fps = 30,
                          file = "particles.mp4", max_delay = 3.0, seed = 1,
                          hold = true)
    fig   = Figure(size = (1600, 800))
    title = Observable("t = 0.00 s")
    ax    = Axis(fig[1, 1]; title = title, xlabel = "X (10 μm)", ylabel = "Y (10 μm)",
                 aspect = DataAspect())
    draw_boundaries!(ax, grid)
    draw_magnets!(ax, mag_grid)

    rng    = Random.MersenneTwister(seed)
    delays = rand(rng, 0:round(Int, max_delay / dt), length(data))   # start delay in steps

    nsteps = maximum(length(p.pos) + d for (p, d) in zip(data, delays))
    frames = 1:every:nsteps

    # before the start: hidden (NaN). after the end: hold at the last point, or hide
    function posat(s)
        map(zip(data, delays)) do (p, d)
            k = s - d
            k < 1 && return Point2f(NaN, NaN)
            k > length(p.pos) && return hold ? Point2f(p.pos[end]) : Point2f(NaN, NaN)
            return Point2f(p.pos[k])
        end
    end

    pos = Observable(posat(1))
    scatter!(ax, pos; marker = Circle, markersize = 2a, markerspace = :data,
             color = (:red, 0.6), strokecolor = :black, strokewidth = 0.5)

    record(fig, file, frames; framerate = fps) do s
        pos[]   = posat(s)
        title[] = "t = $(round((s - 1) * dt, digits = 2)) s"
    end
    return file
end

end