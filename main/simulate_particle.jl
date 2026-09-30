using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays, JLD2, StaticArrays

include("../src/particula/fields_eval.jl")
include("../src/particula/Physics.jl")
include("../src/particula/solvers.jl")
include("../src/particula/structs.jl")

using .evaluate_fields
using .Physics
using .solvers
using .structs


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
    n_particles = 1000
    max_steps = 10000
    dt = 0.1
    #Generar almacenamiento para info de las particulas
    all_particle_data = Vector{ParticleHandler}(undef, n_particles)
    #Seleccion de Caso o proporcion de hemozoina (HZ) en particulas

#caso = Evento() - Falta crear un struct

    #loop para cada n particula
    for i ∈ 1:n_particles

        #Generar particula del tipo struct ::Hemozoin

        #Generar almacenamiento de historia de particula, pos, vel, acc, etc (lo necesario unicamente)

        #Asignar valores iniciales a la particula, pos inicial, vel inicial (depende de pos inicial)
        #y guardar en el espacio [1] del almacenamiento

        #step = 1
        #while step < max_steps

            #posicion anterior = obj.pos (Puede funcionar)

            #aumentar timpo t por un dt, RK4! y aumentar un step 
            #Revisar si hay colision o si llega al final del circuito
            #Guardar historia con valores nuevos

        #Ajustar tamaño de vectores de historia (elimina ceros no utilizados)

        #agregar datos a la historia de todas las particuals en i
    end
    #Guardar la variable de toda la informacion en un JLD2, para posprocesamiento despues
    

    #Ejecutar codigo











    #Guardar valores de las n particulas

    return 
end

