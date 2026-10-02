module structdata
export Hemozoin, Env, ParticleHandler
using StaticArrays
Base.@kwdef mutable struct Hemozoin


    mass::Float64 
    volumen::Float64
    diam::Float64
    dens::Float64
    pos::SVector{2, Float64}
    vel::SVector{2, Float64}
    acc::SVector{2, Float64}
    ΔXᵣ::Float64                    #susceptibilidad
    
end

Base.@kwdef mutable struct environment

    μ₀::Float64 = (4π * 1e-7) * 1e-5 #H/m -> H/u
    η::Float64 = (3.5*1e-3) * 1e-5         #kg/ms -> kg/us
    dens::Float64 = 1056 * 1e-15   #kg/m^2 -> kg/u³ - Densidad agua/plasma
    g::SVector{2, Float64} = (0.0,0.0) .* 1e5 #m/s^2 -> u/s^2
end

struct ParticleHandler
    ΔXᵣ::Float64
    pos::Vector{SVector{2, Float64}}
    vel::Vector{SVector{2, Float64}}
    acc::Vector{SVector{2, Float64}}
end


end