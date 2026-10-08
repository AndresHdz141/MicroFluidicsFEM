module MultiParticle

using GLMakie, Ferrite, FerriteGmsh
import ..DataGraph: draw_boundaries!, draw_magnets!
export plot_particles

function plot_particles(data_infectada, data_sana, grid, mag_grid,
                        location = joinpath("data", "particle_sim_data"),
                        name = "Trayectorias_Particulas.png")

    fig = Figure(size = (1920, 1080))
       ax = Axis(
        fig[1, 1]; 
        title = "Trayectorias Partículas susceptibles vs no susceptibles", 
        xlabel = "X (10 μm)", 
        ylabel = "Y (10 μm)",
        aspect = DataAspect(),
        titlesize = 32,          # Scale up the main title
        xlabelsize = 24,         # Scale up X axis label
        ylabelsize = 24,         # Scale up Y axis label
        xticklabelsize = 18,     # Scale up X tick numbers
        yticklabelsize = 18      # Scale up Y tick numbers
    )

    draw_boundaries!(ax,grid)
    draw_magnets!(ax,mag_grid)
    
    function draw_trajectories!(part, color, name)
        pos = Point2f.(part.pos)
        idx = 1:60:length(pos)
        n = length(idx)

        scatter!(ax, pos[idx], marker = Circle,markersize = 0.8,markerspace = :data, color = color, 
                strokecolor = :black, strokewidth = 0.5, label = name)

    end

    draw_trajectories!(data_sana[5], (:crimson, 0.7), "Healthy RBC")
    draw_trajectories!(data_infectada[5], (:royalblue4, 0.7), "Infected RBC")

        axislegend(
        ax, 
        position = :rt,
        labelsize = 22,          # Size of the text inside the legend
        titlesize = 24,          # Size of the legend title (if added)
        markersize = 30,         # Size of the colored line/marker indicators
        patchsize = (50, 30),    # Width and height of the symbol preview box
        padding = (15, 15, 15, 15) # Inner spacing around the legend content
    )


    mkpath(location)
    save(joinpath(location, name), fig)
    display(fig)
    return fig
end
end