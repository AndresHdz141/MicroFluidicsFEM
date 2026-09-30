module structs
export Hemozoin, Env, ParticleHandler
base.@kwdef multable struct Hemozoin


    mass::Float64 
    volumen::Float64
    diametro::Float64
    dens::Float64
    pos::SVector{2, Float64}
    vel::SVector{2, Float64}
    ΔXᵣ::Float64                    #susceptibilidad
    
end

base.@kwdef multable struct Env

    μ₀::Float64 = 4π * 1e-7  

    η::Float64 = 1.20     #mPa    tal vez  lo cambie pero ne


end

struct ParticleHandler
    ΔXᵣ::Float64
    pos::Vector{SVector{2, Float64}}
    vel::Vector{SVector{2, Float64}}
    acc::Vector{SVector{2, Float64}}
end


end