using Pkg

Pkg.activate(".")
Pkg.instantiate()

include("../src/particula/structs.jl")
include("../src/ReadLineFunctions/interact.jl")
include("magnet_sim.jl")
include("fluid_sim.jl")
include("simulate_particle.jl")

using .structdata
using .ask_data
using .simulate_fluid
using .simulate_magnet
using .simulate_cell

using CSV, DataFrames, JLD2

function main_run(longitud_micro,altura_micro, angulo_canal,dist_imanes_micro, vel_buffer_micro,n_ran, n_iman, n_particulas)

    #Leer valores iniciales. #micro esta en micrometros, u esta en unidad del codigo 
    #10 micrometros = 1 unidad de longitud
    longitud_u = longitud_micro/10
    altura_u= altura_micro/10
    angulo = angulo_canal
    vel_buffer_u = vel_buffer_micro/10
    num_ranuras = n_ran
    cant_imanes = n_iman
    dist_imanes_u = dist_imanes_micro/10
    env = structdata.EnvironmentData()
    env.x_min_c = -longitud_u/2;  env.x_max_c = longitud_u/2                 

    simular_nuevo = ask_simulate_new_field()
    #Generar fluido
    if simular_nuevo
    println("SIMULANDO FLUIDO")
    fluid_sim(longitud_u, altura_u, vel_buffer_u, num_ranuras, "Campo_Fluido", angulo)
    #Generar tres arreglos de imanes
    println("_________________________")
    grid = 
    #Unipolar
    println("SIMULANDO ARREGLO UNIPOLAR")
    mag_sim(longitud_u, altura_u, cant_imanes, 0, "Unipolar", dist_imanes_u)
    println("_________________________")

    #Alterno
    println("SIMULANDO ARREGLO ALTERNO")
    mag_sim(longitud_u, altura_u, cant_imanes, 1, "Alterno", dist_imanes_u)
    println("_________________________")

    #Halbach
    println("SIMULANDO ARREGLO HALBACH")
    mag_sim(longitud_u, altura_u, cant_imanes, 2, "Halbach", dist_imanes_u)
    println("_________________________")
    end
    #Preguntar tipo de particula

    simulate_again = true
    while simulate_again

        
        array_type = ask_array_type()
        fluid_location_data = joinpath(@__DIR__, "..","data", "campo_fluido", "Campo_Fluido.jld2")
        magnet_location_data = joinpath(@__DIR__, "..","data", "arreglos_imanes",array_type)
        #Leer campo de fluido y de campo magnetico
            file = jldopen(fluid_location_data, "r")
            f_name =  file["f_name"]
            f_grid = file["f_mesh"]
            f_dh = file["f_dh"]
            f_field = file["f_u"] 
            close(file)

            file = jldopen(magnet_location_data, "r")
            mag_name =  file["mag_name"]
            mag_grid = file["mag_mesh"]
            mag_dh = file["mag_projector"]
            mag_field = file["mag_H2_field"] 
            close(file)

            println("READ FIELD DATA")
            #Vector de campos de velocidades y magneticos
            f_data = (f_grid, f_dh, f_field)
            mag_data = (mag_grid, mag_dh.dh, mag_field)

        #Extraes info de particula
        particle_path = joinpath(@__DIR__, "..", "data", "particle_data", "ParticleData.csv")
        df_particles = CSV.read(particle_path, DataFrame;header = 1, skipto = 3)
        #Simular particula
        particle_data = select_particle(df_particles)
        println(particle_data)

        efficiency = particle_simulation(
            particle_data[2],
            longitud_u,
            f_data,
            mag_data,
            env,
            particle_data,
            n_particulas
        )
        #Generar Graficos


        simulate_again = ask_simulate_again()
    end
end

#altura canal μm
#longitud canal μm
#angulo entre canales grados

#cantidad de 

#=
    

=#

main_run(
    280,          #longitud Canal   | μm
    40,           #Altura Canal     | μm
    35,           #Angulo del Canal | °
    10,           #Dist iman-canal  | μm
    80,           #Velocidad buffer | μm/s
    5,            #Num. ranuras     | n
    4,            #Num. imanes      | n
    100           #Num. Particulas  | n
    )