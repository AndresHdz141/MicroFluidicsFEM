module solvers

using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays, JLD2, StaticArrays, Random
using ..Physics

export RK4!, pos_inicial



function RK4!(obj::Any, Env::Any,f_data, mag_data, dt::Float64) 
    pos0 = obj.pos
    vel0 = obj.vel
    #k1
    k1x = copy(obj.vel)
    k1v = aceleracion(obj, Env, pos0, vel0, f_data, mag_data)
    #k2
    k2x = obj.vel + 1/2 * k1v * dt
    k2v = aceleracion(obj, Env, pos0 + dt/2 * k1x, vel0 + dt/2 * k1v, f_data, mag_data)
    #k3
    k3x = obj.vel + 1/2 * k2v * dt
    k3v = aceleracion(obj, Env, pos0 + dt/2 * k2x, vel0 + dt/2 * k2v, f_data, mag_data)
    #k4
    k4x = obj.vel + k3v * dt
    k4v = aceleracion(obj, Env, pos0 + k3x * dt, vel0 + k3v * dt, f_data, mag_data)
    #Update
    obj.vel = obj.vel + dt/6 * (k1v + 2*k2v + 2k3v + k4v)
    obj.pos = obj.pos + dt/6 * (k1x + 2*k2x + 2k3x + k4x)
    obj.acc = aceleracion(obj, Env, obj.pos, obj.vel, f_data, mag_data)
end

function get_grid_points(grid::Ferrite.AbstractGrid, setname::String)
    ids = Set{Int}()
    for fi in getfacetset(grid, setname)
        cellid, localfacet = fi
        cell = getcells(grid, cellid)
        union!(ids, Ferrite.facets(cell)[localfacet])
    end
    coords = [grid.nodes[i].x for i in ids]
    points = fill(zero(SVector{2, Float64}), 2)

    points[1] = coords[argmax(v[2] for v in coords)]
    points[2] = coords[argmin(v[2] for v in coords)]
    return points
end 

function pos_inicial(grid::Ferrite.AbstractGrid, setname::String, diametro; sigma = 0.15)
    points = get_grid_points(grid, setname)
    margen_borde = diametro/2
    longitud = norm(points[2]-points[1])

    s_min = margen_borde/longitud
    s_max = 1.0 - s_min

    s = 0.5 + sigma*randn()

    while !(s_min < s < s_max)
        s =0.5 + sigma * randn()
    end
    ini_point = points[1] + s * (points[2]-points[1])
    return SVector{2}(ini_point)

end



    
end