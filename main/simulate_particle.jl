using Pkg
Pkg.activate(".")
Pkg.instantiate()

include("../src/particula/structs.jl")
include("../src/particula/fields_eval.jl")
include("../src/particula/Physics.jl")
include("../src/particula/solvers.jl")

using .evaluate_fields
using .Physics
using .solvers
using .structdata
using Ferrite,FerriteViz, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays, JLD2, StaticArrays, GLMakie


function main()
    #=
        /\_/\
    = ( • . • ) =
      /      \     
    =#  
    #Leer campo de fluido y de campo magnetico
    file = jldopen("data/fluid_data.jld2", "r")
    f_name =  file["f_name"]
    f_grid = file["f_mesh"]
    f_dh = file["f_dh"]
    f_field = file["f_u"] 
    close(file)

    file = jldopen("data/magnetic_data.jld2", "r")
    mag_name =  file["mag_name"]
    mag_grid = file["mag_mesh"]
    mag_dh = file["mag_projector"]
    mag_field = file["mag_H2_field"] 
    close(file)

    #Vector de campos de velocidades y magneticos
    f_data = [f_grid, f_dh, f_field]
    mag_data = [mag_grid, mag_dh, mag_field]
    #generar environment

    env = Env()

    #N de particulas y max_steps
    n_particles = 1
    max_steps = 10
    dt = 0.01
    #Generar almacenamiento para info de las particulas
    all_particle_data = Vector{ParticleHandler}(undef, n_particles)
    #Seleccion de Caso o proporcion de hemozoina (HZ) en particulas

    #caso = Evento() - Falta crear un struct

    #loop para cada n particula
    for i ∈ 1:n_particles

        #Generar particula del tipo struct ::Hemozoin
        hz = Hemozoin(
            mass = 5.0,
            volumen = 3.0,
            diam = 0.8, #μm x10⁻¹
            dens = 0.4,
            pos = SVector{2,Float64}(0.0,0.0),
            vel = SVector{2,Float64}(0.0,0.0),
            acc = SVector{2, Float64}(0.0, 0.0),
            ΔXᵣ = 0.0004
            )

        hz.pos = pos_inicial(f_data[1], "IN1",hz.diam )
        hz.vel = get_field_value(f_data, hz.pos)
        println("Pos0:", hz.pos)
        println("Vel0:", hz.vel)

        #Generar almacenamiento de historia de particula, pos, vel, acc, etc (lo necesario unicamente)
        pos_hist = Vector{SVector{2, Float64}}(undef, max_steps) 
        vel_hist = Vector{SVector{2, Float64}}(undef, max_steps)
        acc_hist = Vector{SVector{2, Float64}}(undef, max_steps)
        #Asignar valores iniciales a la particula, pos inicial, vel inicial (depende de pos inicial)
        #y guardar en el espacio [1] del almacenamiento
        pos_hist[1] = hz.pos
        vel_hist[1] = hz.vel
        acc_hist[1] = hz.acc
        #step = 1
        step = 1
        time = 0.0
        #while step < max_steps
        while step < max_steps
            #posicion anterior = obj.pos (Puede funcionar)
            pos_ant = hz.pos
            #aumentar timpo t por un dt, RK4! y aumentar un step 
            time += dt
            RK4!(hz, env, f_data, mag_data, dt)
            #Revisar si hay colision o si llega al final del circuito
            #Guardar historia con valores nuevos
            step +=1
            pos_hist[step] = hz.pos
            vel_hist[step] = hz.vel
            acc_hist[step] = hz.acc
            println("POS:", hz.pos)

        end
        #Ajustar tamaño de vectores de historia (elimina ceros no utilizados)
        resize!(pos_hist, step)
        resize!(vel_hist, step)
        resize!(acc_hist, step)
        #agregar datos a la historia de todas las particuals en i
        all_particle_data[i] = ParticleHandler(hz.ΔXᵣ, pos_hist, vel_hist, acc_hist)
        println(all_particle_data[i])
    end
    #Guardar la variable de toda la informacion en un JLD2, para posprocesamiento despues
    println("SIMULACION FINALIZADA")
    #Ejecutar codigo
    PlotTray(all_particle_data, f_data[1])
    readline()
    #Guardar valores de las n particulas

    return 
end

function PlotTray(data, grid)
    fig = Figure(size = (1920, 1080))
    ax = Axis(
            fig[1, 1],
            title = "Trayectorias de particulas",
            xlabel = "X (10⁻⁶ m)",
            ylabel = "Y (10⁻⁶ m)",
        )
    for part in data
        lines!(ax, part.pos,color = :red, alpha = 0.5)
    end
    FerriteViz.wireframe
    display(fig)
end

main()























function test()

    file = jldopen("data/fluid_data.jld2", "r")
    f_name =  file["f_name"]
    f_grid = file["f_mesh"]
    f_dh = file["f_dh"]
    f_field = file["f_u"] 
    close(file)

    file = jldopen("data/magnetic_data.jld2", "r")
    mag_name =  file["mag_name"]
    mag_grid = file["mag_mesh"]
    mag_dh = file["mag_projector"]
    mag_field = file["mag_H2_field"] 
    close(file)

    #Vector de campos de velocidades y magneticos
    f_data = [f_grid, f_dh, f_field]
    mag_data = [mag_grid, mag_dh, mag_field]

    points = get_grid_points(f_data[1], "IN1")
    ps = [pos_inicial(points, 0.7) for i in 1:200]
    all_points = vcat([points[1], points[2]], ps)

    xs = getindex.(all_points, 1)
    ys = getindex.(all_points, 2)
    p = scatter(xs, ys; aspect_ratio = :equal, legend = false)
    display(p)
    readline()
end
