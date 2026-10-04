module structdata
export Hemozoin, EnvironmentData, ParticleHandler
using StaticArrays
Base.@kwdef mutable struct Hemozoin


    mass::Float64                   #masa
    volumen::Float64                #Volumen
    diam::Float64                   #diametro
    dens::Float64                   #Densidad
    pos::SVector{2, Float64}        #posicion  
    vel::SVector{2, Float64}        #velocidad
    acc::SVector{2, Float64}        #aceleracion
    ΔXᵣ::Float64                    #susceptibilidad
    
end

Base.@kwdef mutable struct cases

    

end

Base.@kwdef mutable struct EnvironmentData

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