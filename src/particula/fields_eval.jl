using Pkg
Pkg.activate(".")
Pkg.instantiate()

using StaticArrays, JLD2, Ferrite

#=

point = SVector{2, Float64}(0.0,0.0)
mag_path = "magnetic_data.jld2"
mag_grid = load(mag_path, "mag_mesh")
mag_field = load(mag_path, "mag_H2_field")

PointEvalHandler(mag_grid, point)

case = CaseOfEval("ring")
for

part = Hemozoin(ΔXᵣ = case.ΔXᵣ)

pos0 = setposition, ya la sabemos, distribucion

vel0 = getVelocity(field, :u)

part.pos = pos0
part.vel = vel0

pos_history = zeros(max_iter)
vel_history = zeros(max_iter)

pos_history[1] = part.pos
pos0

for  in maxIter, if end_colision
    rk4!
    if check_wall_collision(part.pos) -> velocity = part.pos * 0.1

end

end

=#

file = jldopen("data/fluid_data.jld2", "r")
f_name =  file["f_name"]
f_grid = file["f_mesh"]
f_dh = file["f_dh"]
f_u = file["f_u"] 
close(file)

file = jldopen("data/magnetic_data.jld2", "r")
mag_name =  file["mag_name"]
mag_grid = file["mag_mesh"]
mag_proj = file["mag_projector"]
mag_field = file["mag_H2_field"] 
close(file)

println("Data type ",typeof(f_u))
println("Data read ",typeof(mag_field))


function get_field_value(grid, point, dh,field )
    ph = PointEvalHandler(grid, [point])
    value = typeof(field) == Vector{Float64} ? "fluid" : "magnetic" #evaluate_at_points(ph, dh, field, :u ) :  evaluate_at_points(ph, dh, field)
    println(value)
end

point = Vec((0.0,0.00))

get_field_value(f_grid, point, f_dh, f_u)
get_field_value(mag_grid, point, mag_proj, mag_field)