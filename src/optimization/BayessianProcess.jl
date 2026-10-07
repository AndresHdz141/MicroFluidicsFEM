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

using BayesianOptimization, GaussianProcesses,CSV, DataFrames, JLD2, Random

function optimize_function()

best_run_tracker = (
    max_efficiency = Ref(0.0),
    best_individuals = zeros(4)
)

function blackbox_function(x)
    
    longitud_micro = x[1]
    altura_micro = x[2]
    angulo = x[3]
    dist_imanes_micro = x[4]
    vel_buffer_micro = x[5]
    num_ranuras = round(Int, x[6])
    cant_imanes = round(Int, x[7])
    
    
    n_particulas = 100
    longitud_u = longitud_micro/10
    altura_u= altura_micro/10
    vel_buffer_u = vel_buffer_micro/10
    dist_imanes_u = dist_imanes_micro/10


    env = structdata.EnvironmentData()
    env.x_min_c = -longitud_u/2;  env.x_max_c = longitud_u/2 

    names = ("Unipolar", "Alterno", "Halbach")

    # "0 - Unipolar | 1 - Alterno | 2 - Halbach"
    select_array = 0

    if select_array == 0
    array_type = names[1]

    elseif select_array == 1
    array_type = names[2]

    elseif select_array == 2
    array_type = names[3]

    end


    fluid_sim(longitud_u, altura_u, vel_buffer_u, num_ranuras, "Campo_Fluido", angulo)

    mag_sim(longitud_u, altura_u, cant_imanes, select_array, array_type, dist_imanes_u)


    mag_rad_name = array_type *".jld2"
    fluid_location_data = joinpath(@__DIR__,"..", "..","data", "campo_fluido", "Campo_Fluido.jld2")
    magnet_location_data = joinpath(@__DIR__, "..","..","data", "arreglos_imanes",mag_rad_name)
        
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

    f_data = (f_grid, f_dh, f_field)
    mag_data = (mag_grid, mag_dh.dh, mag_field)

    particle_path = joinpath(@__DIR__, "..", "..","data", "particle_data", "ParticleData.csv")
    df_particles = CSV.read(particle_path, DataFrame;header = 1, skipto = 3)
    total_eficiencia = 0
    current_individuals = zeros(4)

    for i in 1:4
        matching_rows = filter(row -> row.tag == i, df_particles)
        particle_data = Vector(matching_rows[1, :])

        efficiency = particle_simulation(
            particle_data[2],
            longitud_u,
            f_data,
            mag_data,
            env,
            particle_data,
            n_particulas
        )
    total_eficiencia += efficiency
    end

    if total_eficiencia > best_run_tracker.max_efficiency[]
        best_run_tracker.max_efficiency[] = total_eficiencia
        best_run_tracker.best_individuals .= current_individuals # Copy values in-place
    end

    return -(total_eficiencia)

end

#Definir bounds del espacio de busqueda
lower_bounds = [200, 35, 30, 5, 80,2, 2]
upper_bounds = [300, 45, 45, 40, 120,6, 8]
println("Bounds")

num_dimensions = length(lower_bounds) # 7
num_initial_points = 20  
# Creates a 7 × 20 matrix
X_init = zeros(num_dimensions, num_initial_points)
for i in 1:num_dimensions
    # Scale random numbers (0 to 1) into your specific parameter ranges
    X_init[i, :] = lower_bounds[i] .+ rand(num_initial_points) .* (upper_bounds[i] - lower_bounds[i])
end

y_init = [blackbox_function(X_init[:, j]) for j in 1:num_initial_points]

gp_mean = MeanZero()
gp_kernel = SEArd(zeros(7), 0.0) 

gp_surrogate = GP(X_init, y_init, gp_mean, gp_kernel, -2.0)
println("Surrogate")

opt = BOpt(blackbox_function, 
           gp_surrogate, 
           ExpectedImprovement(), 
           NoModelOptimizer(), # O MAPGPOptimizer() si deseas actualizar las escalas del GP
           lower_bounds, 
           upper_bounds; 
           maxiterations = 3)

println("Beginning core Bayesian Optimization loops...")

result = boptimize!(opt)

println(result)
    println("\n=============================================")
    println("OPTIMIZATION BREAKDOWN COMPLETE")
    println("=============================================")
    println("Best Global Efficiency Found: ", best_run_tracker.max_efficiency[])
    println("Optimal Parameter Set: ", result.observed_optimizer)
    println("\nIndividual Efficiency Breakdown for this Optimal Set:")
    for i in 1:4
        println("  • Particle Tag $i: ", best_run_tracker.best_individuals[i])
    end
    println("=============================================")

    return result
end

optimize_function()
