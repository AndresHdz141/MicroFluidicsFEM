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
    
    function get_velocity(obj::Hemozoin, env::EnvironmentData, c_pos, f_data,  fluid_eval, mag_data, mag_eval)
        F_m = magforce(obj, env, mag_data, mag_eval, c_pos)
        τ = 3*pi*env.η*obj.diam
        F_m_effective = F_m ./τ
        println("MAGNETIC DRIVING FORCE>", F_m_effective)
        fluid_velocity = fast_evaluate_field(c_pos, fluid_eval[1], f_data[1], f_data[2], f_data[3], fluid_eval[2], :u)
        println("FLUID VEL>" ,fluid_velocity)

        v_particle = fluid_velocity+F_m_effective
        return v_particle
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

