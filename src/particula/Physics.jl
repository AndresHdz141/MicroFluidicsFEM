module Physics
    using ..structdata
    using ..evaluate_fields
    using StaticArrays,LinearAlgebra

    export aceleracion
    function aceleracion(obj::Hemozoin, env::Env, c_pos, c_vel, f_data, mag_data)
        F_m = magforce(obj, env, mag_data, c_pos)

        F_d = dragforce(obj, env, f_data,c_pos, c_vel )

        F_total = F_m + F_d
    return F_total / obj.mass
    end
    function magforce(obj::Hemozoin, env::Env, mag_data, c_pos)
        μ₀ = env.μ₀
        Volumen_p = obj.volumen
        ΔXᵣ = obj.ΔXᵣ
        ∇H² = get_field_value(mag_data, c_pos)
        return (1/2) * μ₀ * Volumen_p * ΔXᵣ * ∇H²
    end
    function dragforce(obj::Hemozoin, env::Env, f_data, c_pos, c_vel)
        relative_velocity = get_field_value(f_data, c_pos) - c_vel
        return 3*π*env.η*obj.diam*relative_velocity 
    end
    function gbforce(obj::Hemozoin, env::Env)
        density_diff = obj.dens - env.dens
        Volumen_p = obj.volumen
        g = env.g
        return density_diff * Volumen_p * g
    end
end