using Pkg

Pkg.activate(".")
Pkg.instantiate()

include("../particula/structs.jl")
include("../ReadLineFunctions/interact.jl")
include("../../main/magnet_sim.jl")
include("../../main/fluid_sim.jl")
include("../../main/simulate_particle.jl")

using .structdata
using .ask_data
using .simulate_fluid
using .simulate_magnet
using .simulate_cell

using BayesianOptimization, GaussianProcesses, CSV, DataFrames, JLD2, Random

function optimize_function()

 best_run_tracker = (
        max_pseudo = Ref(0.0),       # Lo que guía a la IA
        max_real_eff = Ref(0.0),     # Lo que te interesa a ti
        best_individuals = zeros(4)  # El porcentaje físico de cada tag
    )

    particle_path = joinpath(@__DIR__, "..", "..", "data", "particle_data", "ParticleData.csv")
    df_particles = CSV.read(particle_path, DataFrame; header = 1, skipto = 3)

    function blackbox_function(u)
        p = decode(u)

        longitud_micro    = p.longitud
        altura_micro      = p.altura
        angulo            = p.angulo
        dist_imanes_micro = p.dist_imanes
        vel_buffer_micro  = p.vel_buffer
        num_ranuras       = p.num_ranuras
        cant_imanes       = p.cant_imanes
        println("longitud: ", p.longitud)
        println("altura_micro: ", p.altura)
        println("angulo: ", p.angulo)
        println("dist_imanes_micro: ", p.dist_imanes)
        println("vel_buffer_micro: ", p.vel_buffer)
        println("num_ranuras: ", p.num_ranuras)
        println("cant_imanes: ", p.cant_imanes)


        n_particulas = 500
        longitud_u = longitud_micro / 10
        altura_u = altura_micro / 10
        vel_buffer_u = vel_buffer_micro / 10
        dist_imanes_u = dist_imanes_micro / 10

        env = structdata.EnvironmentData()
        env.x_min_c = -longitud_u/2
        env.x_max_c = longitud_u/2 

        names = ("Unipolar", "Alterno", "Halbach")

        # "0 - Unipolar | 1 - Alterno | 2 - Halbach"
        select_array = 2
        array_type = names[select_array + 1] # Simplified selection logic

        # Run simulations
        fluid_sim(longitud_u, altura_u, vel_buffer_u, num_ranuras, "Campo_Fluido", angulo)
        mag_sim(longitud_u, altura_u, cant_imanes, select_array, array_type, dist_imanes_u)

        mag_rad_name = array_type * ".jld2"
        fluid_location_data = joinpath(@__DIR__, "..", "..", "data", "campo_fluido", "Campo_Fluido.jld2")
        magnet_location_data = joinpath(@__DIR__, "..", "..", "data", "arreglos_imanes", mag_rad_name)
            
        file = jldopen(fluid_location_data, "r")
        f_grid = file["f_mesh"]
        f_dh = file["f_dh"]
        f_field = file["f_u"] 
        close(file)

        file = jldopen(magnet_location_data, "r")
        mag_grid = file["mag_mesh"]
        mag_dh = file["mag_projector"]
        mag_field = file["mag_H2_field"] 
        close(file)

        f_data = (f_grid, f_dh, f_field)
        mag_data = (mag_grid, mag_dh.dh, mag_field)

        total_pseudo = 0.0
        total_real = 0.0
        current_real_individuals = zeros(4)
        for i in 2:4
            matching_rows = filter(row -> row.tag == i, df_particles)
            particle_data = Vector(matching_rows[1, :])

            pseudo_val,real_eff, data = particle_simulation(
                particle_data[2],
                altura_u,
                f_data,
                mag_data,
                env,
                particle_data,
                n_particulas
            )
                current_real_individuals[i] = real_eff
                total_pseudo += pseudo_val
                total_real += real_eff
        end
        pseudo_promedio = total_pseudo / 3.0
        
        if pseudo_promedio > best_run_tracker.max_pseudo[]
            best_run_tracker.max_pseudo[] = pseudo_promedio
            best_run_tracker.max_real_eff[] = total_real
            best_run_tracker.best_individuals .= current_real_individuals
        end
        
        return pseudo_promedio
    end

lower_bounds = [270.0, 35.0, 30.0, 5.0, 80.0, 2.0, 2.0]
upper_bounds = [300.0, 45.0, 45.0, 25.0, 120.0, 6.0, 8.0]

to_physical(u) = lower_bounds .+ u .* (upper_bounds .- lower_bounds)
known_good_var = [280.0, 40.0, 35.0, 10.0, 80.0, 5.0, 4.0]
known_good_u = (known_good_var .- lower_bounds) ./ (upper_bounds .- lower_bounds)

function decode(u)
    x = to_physical(u)
    return (
        longitud     = x[1],
        altura       = x[2],
        angulo       = x[3],
        dist_imanes  = x[4],
        vel_buffer   = x[5],
        num_ranuras  = clamp(round(Int, x[6]), 2, 6),
        cant_imanes  = clamp(round(Int, x[7]), 2, 8),
    )
end
# inside blackbox_function(u):  x = to_physical(u)  and use x[1], x[2], ...

X_init = rand(7, 20)                                  # unit cube
X_init[:,1] = known_good_u

y_init = [blackbox_function(X_init[:, j]) for j in 1:20]

model = ElasticGPE(7,
    mean = MeanConst(0.0),
    kernel = SEArd(fill(log(0.3), 7), log(50.0)),     # 7 log-lengthscales + log σ
    logNoise = -2.0,
    capacity = 3000)
append!(model, X_init, y_init)

modeloptimizer = MAPGPOptimizer(
    every = 10,
    noisebounds = [-4.0, 1.0],
    kernbounds = [[fill(-3.0, 7); -1.0], [fill(1.0, 7); 6.0]],
    maxeval = 100)

opt = BOpt(blackbox_function, model, ExpectedImprovement(), modeloptimizer,
           zeros(7), ones(7);
           maxiterations = 50,
           initializer_iterations = 5)

result = boptimize!(opt)
best = decode(result.observed_optimizer)
    println("===========================")
    println("Best efficiency: ", result.observed_optimum)
    println("Best parameters: ", best)
    println("Per-particle breakdown: ", best_run_tracker.best_individuals)
    println("\nIndividual Efficiency Breakdown for this Optimal Set:")
    for i in 1:4
        println("  • Particle Tag $i: ", best_run_tracker.best_individuals[i])
    end
    println("=============================================")

    # save it
    CSV.write(joinpath(@__DIR__, "best_params.csv"), DataFrame([best]))
    jldsave(joinpath(@__DIR__, "best_params.jld2"); best, efficiency = result.observed_optimum)
    return result
end

optimize_function()