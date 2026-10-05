using Pkg

Pkg.activate(".")
Pkg.instantiate()

include("../src/particula/structs.jl")
include("../src/particula/fields_eval.jl")
include("../src/particula/Physics.jl")
include("../src/particula/solvers.jl")
include("../src/Data_analisis/data.jl")
include("../src/Data_analisis/GraphData.jl")

using .structdata
using .evaluate_fields
using .Physics
using .solvers
using .Data
using .DataGraph
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays, JLD2, StaticArrays, GLMakie


function main(name)
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

    println("READ FIELD DATA")
    #Vector de campos de velocidades y magneticos
    f_data = (f_grid, f_dh, f_field)
    mag_data = (mag_grid, mag_dh.dh, mag_field)
    #Obtener arbles del campo de fluidos y magnetico
    fluid_eval = setup_fast_eval(f_data[1], f_data[2], :u)
    mag_eval = setup_fast_eval(mag_data[1], mag_data[2], :_)
    println("GENERATED DATA TREES")
    #generar Environment
    length = 28.0
    env = structdata.EnvironmentData()
    env.x_min_c = -length/2;  env.x_max_c = length/2                 

    println("GENERATED ENVIRONMENT DATA")
    #N de particulas y max_steps
    n_particles = 500
    dt = 0.01

    sim_seconds = 10
    max_steps = round(Int, sim_seconds / dt) 
    save_interval = 1
    hist_size = cld(max_steps, save_interval)
    all_particle_data = Vector{structdata.ParticleHandler}(undef, n_particles)
    #Seleccion de Caso o proporcion de hemozoina (HZ) en particulas
    #caso = Evento() - Falta crear un struct
    #loop para cada n particula
    for i ∈ 1:n_particles
        #Generar particula del tipo struct ::Hemozoin
        hz = Hemozoin(
            mass = 3.3* 1e-14, #kg
            volumen = (4/3)*pi*(4)^3 * 1e-3,  #μm^3 -> u^-3
            diam = 8 * 1e-1,       #μm -> 10^-1 u (e-6-1)
            dens = 1.15* 1e-12,     #g/m^3 -> 10^-12 kg/u^3
            pos = SVector{2,Float64}(0.0,0.0), #u -> 10μm
            vel = SVector{2,Float64}(0.0,0.0), #u/2 -> 10μm/s
            acc = SVector{2, Float64}(0.0, 0.0), #u/s^2 -> 10μm/s^2
            ΔXᵣ = 1.3   #Adimensional
            )

        hz.pos = pos_inicial(f_data[1], "IN1",hz.diam )
        hz.vel = fast_evaluate_field(hz.pos, fluid_eval[1], f_data[1], f_data[2], f_data[3], fluid_eval[2], :u)
        #println("Pos0:", hz.pos)
        #println("Vel0:", hz.vel)

        #Generar almacenamiento de historia de particula, pos, vel, acc, etc (lo necesario unicamente)
        pos_hist = Vector{SVector{2, Float64}}(undef, hist_size) 
        vel_hist = Vector{SVector{2, Float64}}(undef, hist_size)
        acc_hist = Vector{SVector{2, Float64}}(undef, hist_size)
        #Asignar valores iniciales a la particula, pos inicial, vel inicial (depende de pos inicial)
        #y guardar en el espacio [1] del almacenamiento
        pos_hist[1] = hz.pos
        vel_hist[1] = hz.vel
        acc_hist[1] = hz.acc
        #step = 1
        step = 1
        time = 0.0
        colision = false
        #while step < max_steps
        while step < max_steps
            #posicion anterior = obj.pos (Puede funcionar)
            pos_ant = hz.pos
            #aumentar timpo t por un dt, RK4! y aumentar un step 
            time += dt
            if !colision
            RK4!_vel(hz, env, f_data, fluid_eval, mag_data, mag_eval, dt)
            #println("RK4[",step, "]")
                in_zone = env.x_min_c < hz.pos[1] < env.x_max_c
                if in_zone && hz.pos[2] < 0.0 + hz.diam/2+env.gap_cap
                    fluid_vel =  fast_evaluate_field(pos_ant, fluid_eval[1], f_data[1], f_data[2], f_data[3], fluid_eval[2], :u)
                    hz.vel = fluid_vel * 0.0
                    hz.pos = pos_ant + hz.vel*dt
                    hz.acc = SVector{2,Float64}(0.0,0.0)
                    step +=1
                    if step % save_interval ==0
                        pos_hist[step] = hz.pos
                        vel_hist[step] = hz.vel
                        acc_hist[step] = hz.acc
                    end
                    colision = true
                    break
                end
            end
            #Si la posicion[2] part <=
            #Revisar si hay colision o si llega al final del circuito
            #Guardar historia con valores nuevos
            step +=1
            if step % save_interval ==0
            pos_hist[step] = hz.pos
            vel_hist[step] = hz.vel
            acc_hist[step] = hz.acc
            end
            #println("POS:", hz.pos)

        end
        #Ajustar tamaño de vectores de historia (elimina ceros no utilizados)
        resize!(pos_hist, step)
        resize!(vel_hist, step)
        resize!(acc_hist, step)
        #agregar datos a la historia de todas las particuals en i
        all_particle_data[i] = structdata.ParticleHandler(hz.ΔXᵣ, pos_hist, vel_hist, acc_hist)
    end
    #Guardar la variable de toda la informacion en un JLD2, para posprocesamiento despues
    println("SIMULACION FINALIZADA")
    join_name = name * ".jld2"
    path = joinpath("data", "particulas", join_name)
    jldsave(path;
        data = all_particle_data
    )

    part_capt, part_nocapt = inenout(all_particle_data)

    println("total = ", n_particles)
    println("capturadas = ", part_capt)
    println("No capturadas = ", part_nocapt)
    println("Eficiencia: ", part_capt/n_particles*100, "%")
    #Ejecutar codigo
    AnimateParticles(all_particle_data, f_data[1], mag_data[1], 0.4, dt)
    PlotTray(all_particle_data, f_data[1], mag_data[1], env, 0.8/2)
    readline()
    #Guardar valores de las n particulas

    return 
end



main("RBC_TROPHOZOITE")

