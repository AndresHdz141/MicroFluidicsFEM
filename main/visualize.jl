using ReadVTK, GLMakie

vtk = VTKFile("../stokes-flow.vtu")
points = get_points(vtk)                  # matrix of node coordinates
velocity = get_data(get_point_data(vtk)["u"])  # your velocity field, however you named it

# Quick scatter colored by velocity magnitude
mag = [norm(velocity[:, i]) for i in axes(velocity, 2)]
scatter(points[1, :], points[2, :], color = mag, colormap = :viridis)

readline()