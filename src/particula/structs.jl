module structs

base.@kwdef multable struct Hemozoin


    mass::Float64 
    volumen::Float64
    dens::Float64
    pos::SVector{3, Float64}
    vel::SVector{3, Float64}
    ΔXᵣ::Float64                    #susceptibilidad



end

base.@kwdef multable struct Env

    μ₀::Float64 = 4π * 1e-7  

    η::Float64 = 1.20     #mPa    tal vez  lo cambie pero ne


end


end