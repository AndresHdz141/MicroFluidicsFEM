module Physics
using StaticArrays
using LinearAlgebra

export aceleracion
function aceleracion(obj::Hemozoin, env::Env,)
    F_m = magforce(obj, env, position, velocity)
    F_d = dragforce(obj, env, position, velocity)
    F_total = F_m + F_d
return F_total / obj.mass
end
function magforce()
grad_H2 = getHvalue(obj.pos)

end
function dragforce(obj::Hemozoin, env::Env, )
    relative_velocity = velocityfield - obj.vel
end
function gbforce(obj::Hemozoin, env::Env)
    density_diff = obj.density - env.density
end
end