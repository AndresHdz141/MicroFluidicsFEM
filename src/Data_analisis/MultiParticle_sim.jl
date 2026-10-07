module MultiParticle

using GLMakie, Ferrite, FerriteGmsh
import ..DataGraph: draw_boundaries!, draw_magnets!
export plot_particles

function plot_particles(data_infectada, data_sana, grid, mag_grid,
                        location = joinpath("data", "particle_sim_data"),
                        name = "Trayectorias_Particulas.png")

    fig = Figure(size = (1920, 1080))
    ax = Axis(fig[1, 1]; title = "Trayectorias Partículas susceptibles vs no susceptibles", xlabel = "X (10 μm)", ylabel = "Y (10 μm)",
                aspect = DataAspect())

    draw_boundaries!(ax,grid)
    draw_magnets!(ax,mag_grid)
    
    function draw_trajectories!(part, color, name)
        pos = Point2f.(part.pos)
        idx = 1:70:length(pos)
        n = length(idx)

        scatter!(ax, pos[idx], marker = Circle, color = color, 
                strokecolor = :black, strokewidth = 0.5, label = name)

    end

    draw_trajectories!(data_sana[5], (:crimson, 0.7), "Healthy RBC")
    draw_trajectories!(data_infectada[5], (:royalblue4, 0.7), "Infected RBC")

    axislegend(ax)

    mkpath(location)
    save(joinpath(location, name), fig)
    display(fig)
    return fig
end
end