module solvers
using ..Physics

export RK4!


function RK4!(Hemozoin::Any, Env::Any, dt::Float64) 
    pos0 = obj.pos
    vel0 = obj.vel
    #k1
    k1x = copy(obj.vel)
    k1v = aceleracion(obj, env, pos0, vel0)
    #k2
    k2x = obj.vel + 1/2 * k1v * dt
    k2v = aceleracion(obj, env, pos0 + dt/2 * k1x, vel0 + dt/2 * k1v)
    #k3
    k3x = obj.vel + 1/2 * k2v * dt
    k3v = aceleracion(obj, env, pos0 + dt/2 * k2x, vel0 + dt/2 * k2v)
    #k4
    k4x = obj.vel + k3v * dt
    k4v = aceleracion(obj, env, pos0 + k3x * dt, vel0 + k3v * dt)
    #Update
    obj.vel = obj.vel + dt/6 * (k1v + 2*k2v + 2k3v + k4v)
    obj.pos = obj.pos + dt/6 * (k1x + 2*k2x + 2k3x + k4x)
    obj.acc = aceleracion(obj,env, obj.pos, obj.vel)
end


end