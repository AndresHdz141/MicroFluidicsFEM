module solvers

using Random, StaticArrays, LinearAlgebra
using ..Physics

export RK4!, dist_P



function RK4!(obj::Any, Env::Any, dt::Float64) 
    pos0 = obj.pos
    vel0 = obj.vel
    #k1
    k1x = copy(obj.vel)
    k1v = aceleracion(obj, Env, pos0, vel0)
    #k2
    k2x = obj.vel + 1/2 * k1v * dt
    k2v = aceleracion(obj, Env, pos0 + dt/2 * k1x, vel0 + dt/2 * k1v)
    #k3
    k3x = obj.vel + 1/2 * k2v * dt
    k3v = aceleracion(obj, Env, pos0 + dt/2 * k2x, vel0 + dt/2 * k2v)
    #k4
    k4x = obj.vel + k3v * dt
    k4v = aceleracion(obj, Env, pos0 + k3x * dt, vel0 + k3v * dt)
    #Update
    obj.vel = obj.vel + dt/6 * (k1v + 2*k2v + 2k3v + k4v)
    obj.pos = obj.pos + dt/6 * (k1x + 2*k2x + 2k3x + k4x)
    obj.acc = aceleracion(obj, Env, obj.pos, obj.vel)
end

function pos_inicial(p10 , p11, diametro; sigma = 0.15)

    margen_borde = diametro/2
    longitud = norm(p11-p10)

    s_min = margen_borde/longitud
    s_max = 1.0 - s_min

    s = 0.5 + sigma*randn()

    while !(s_min < s < s_max)
        s =0.5 + sigma * randn()
    end

    return p10 + s * (p11-p10)

end



    
end