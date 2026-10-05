module Physics
    using ..structdata
    using ..evaluate_fields
    using StaticArrays,LinearAlgebra

    export aceleracion,get_velocity
    
    function aceleracion(obj::Hemozoin, env::EnvironmentData, c_pos, c_vel, f_data, fluid_eval, mag_data, mag_eval)
        F_m = magforce(obj, env, mag_data, mag_eval, c_pos)
        #println("MagneticForce:", F_m)
        F_d = dragforce(obj, env, f_data, fluid_eval,c_pos, c_vel)
        #println("DragForce:", F_d)
        F_total = F_m + F_d
    return F_total / obj.mass
    end
    
    #La velocidad toma valores de arrastre distintos, cambian segun la distancia a la pared (solo se asume la pared en y=0)
    function get_velocity(obj::Hemozoin, env::EnvironmentData, c_pos, f_data,  fluid_eval, mag_data, mag_eval)
        a = obj.diam / 2
        F_m = magforce(obj, env, mag_data, mag_eval, c_pos)
        τ = 6 * pi * env.η * a
        u_f = fast_evaluate_field(c_pos, fluid_eval[1], f_data[1], f_data[2], f_data[3], fluid_eval[2], :u)
        #Fuera de donde esta la pared inferior (saldas y entradas)
        (env.x_min_c < c_pos[1] < env.x_max_c) || return u_f + F_m / τ
        #Distancia del centro de la particula a la pared dividida entre radio de particula
        ξ = max(c_pos[2] - a, env.δc) / a
        #Ajuste en donde la correcion solo toma lugar cuando la distancia a la pared/radio es menor a 3
        ξ > 3 && return u_f + F_m / τ
        #Arrastre con correccion de Faxén cerca de las pareces.
        λp, λn = wall_drag(ξ)
        #vel_fluido veces el shear se le suma la fuerza magnetica dividido entre el arrastre corregido
        v_x = u_f[1] + F_m[1] / (τ * λp)  #f_shear(ξ)*
        #velocidad de fluido se considera 0 (La particula no sube ni baja), Fuerza magnetica con arrastre corregido
        v_y = u_f[2]+ F_m[2] / (τ * λn)
        #Regresa ambos valores en sus componentes como vector
        return SVector{2, Float64}(v_x, v_y)
    end
    #Arrastre de la pared con correcion de Faxen
    function wall_drag(ξ)
        r = 1 / (1 + ξ)                                   # a / z
        den_p = 1 - 9/16*r + 1/8*r^3          # Faxén al termino ^3 en paralelo a superficie
        den_n = 1 - 9/8*r + 1/2*r^3           # Faxen al termino ^3 en normal a superficie
        return 1 / max(den_p, 1e-3), 1 / max(den_n, 1e-3)    #parallel and normal parts
    end
    #
    function f_shear(ξ)
        r = 1/(1+ξ)
        far = 1 - 5/16 * r^3
        near = 0.7431 /(0.6376 - 0.2 * log(ξ))
        w = clamp((ξ - 0.5) / 1.5, 0.0, 1.0)           # blend ξ = 0.5 → 2
        return (1 - w) * near + w * far
    end

    function magforce(obj::Hemozoin, env::EnvironmentData, mag_data, mag_eval, c_pos)
        μ₀ = env.μ₀
        Volumen_p = obj.volumen
        ΔXᵣ = obj.ΔXᵣ
        ∇H² = fast_evaluate_field(c_pos,mag_eval[1], mag_data[1], mag_data[2], mag_data[3], mag_eval[2], :_) * 1e-5
        return (1/2) * μ₀ * Volumen_p * ΔXᵣ * ∇H²
    end
    function dragforce(obj::Hemozoin, env::EnvironmentData, f_data, fluid_eval, c_pos, c_vel)
        fluid_vel = fast_evaluate_field(c_pos, fluid_eval[1], f_data[1], f_data[2], f_data[3], fluid_eval[2], :u)
        relative_velocity = fluid_vel - c_vel
        return 3*π*env.η*obj.diam*relative_velocity 
    end
    function gbforce(obj::Hemozoin, env::EnvironmentData)
        density_diff = obj.dens - env.dens
        Volumen_p = obj.volumen
        g = env.g
        return density_diff * Volumen_p * g
    end
end

