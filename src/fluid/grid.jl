module grid
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays

 export setup_grid

    function setup_grid(L, height, phi,n_ranuras, h=0.5)
        a = L/2
        Gmsh.initialize()
        gmsh.clear()

        gmsh.option.set_number("General.Verbosity", 2)

        #Coordenadas de puntos
        n = n_ranuras     # Cantidad de ranuras
        s = L/4n  # 2s = grosor de la ranura
        #X
        x1 = a
        x2 = a+a*cosd(phi)
        x3 =a+(a+height)*cosd(phi)
        x4 = a + height*cosd(phi)
        x5 = a + (a+height)*cosd(phi)
        x6 = a + a*cosd(phi)
        x7 = a
        # Coordenadas de los puntos inferiores de las ranuras
        x_ranuras = Float64[]                 # Matriz para almacenar las coordenadas x de las ranuras
        x_0 = x7 - s                          # Coordenada x de la primera ranura
        for i in 1:2n-1                       # Generar las coordenadas x de las ranuras
            push!(x_ranuras, x_0 - i * 2 * s) 
        end
        
        #Y
        y1 = 0
        y2 = -a*sind(phi)
        y3 = -a*sind(phi)+height/2
        y4 = height/2
        y5 = a*sind(phi)+height/2
        y6 = height+a*sind(phi)
        y7 = height
        # Coordenadas de los puntos superiores de las ranuras
        y_ranuras = height + height/3 # Altura de las ranuras

        #Add points
        #Origin
        o = gmsh.model.geo.add_point(0.0, 0.0,0.0, h)

        #Right symmetry
        p1 = gmsh.model.geo.add_point(x1, y1,0.0, h)
        p2 = gmsh.model.geo.add_point(x2, y2,0.0, h)
        p3 = gmsh.model.geo.add_point(x3, y3,0.0, h)
        p4 = gmsh.model.geo.add_point(x4, y4,0.0, h)
        p5 = gmsh.model.geo.add_point(x5, y5,0.0, h)
        p6 = gmsh.model.geo.add_point(x6, y6,0.0, h)
        p7 = gmsh.model.geo.add_point(x7, y7,0.0, h)

        #Groove points
        p8 = gmsh.model.geo.add_point(x_0, height,0.0, h)       # Punto de la primera ranura
        p8_up = gmsh.model.geo.add_point(x_0, y_ranuras,0.0, h) # Punto superior de la primera ranura
        points_ranuras_inf = Int[] # Matriz para almacenar los puntos inferiores de las ranuras
        points_ranuras_sup = Int[] # Matriz para almacenar los puntos superiores de las ranuras
        for i in 1:2n-1            # Generar los puntos inferiores y superiores de las ranuras
            push!(points_ranuras_inf, gmsh.model.geo.add_point(x_ranuras[i], height,0.0, h))
            push!(points_ranuras_sup, gmsh.model.geo.add_point(x_ranuras[i], y_ranuras,0.0, h))
        end

        #Left symmetry
        p9 = gmsh.model.geo.add_point(-x7, y7,0.0, h)
        p10 = gmsh.model.geo.add_point(-x6, y6,0.0, h)
        p11 = gmsh.model.geo.add_point(-x5, y5,0.0, h)
        p12 = gmsh.model.geo.add_point(-x4, y4,0.0, h)
        p13 = gmsh.model.geo.add_point(-x3, y3,0.0, h)
        p14 = gmsh.model.geo.add_point(-x2, y2,0.0, h)
        p15 = gmsh.model.geo.add_point(-x1, y1,0.0, h)

        #Generate lines
        #Pared inferior 14, 15, o, 1, 2
        #Pared derecha 3, 4, 5
        #Pared superior 6,7,ranuras,9,10
        #Pared izquieda 11,12,13

        #Inlets 9-10 y 12-13
        #Outlets 2-3 y 5-6

        l1 = gmsh.model.geo.add_line(p14,p15)
        l2 = gmsh.model.geo.add_line(p15,o)
        l3 = gmsh.model.geo.add_line(o,p1)
        l4 = gmsh.model.geo.add_line(p1,p2) 
        l5 = gmsh.model.geo.add_line(p2,p3) #OUTLET
        l6 = gmsh.model.geo.add_line(p3,p4)
        l7 = gmsh.model.geo.add_line(p4,p5)
        l8 = gmsh.model.geo.add_line(p5,p6) #OUTLET
        l9 = gmsh.model.geo.add_line(p6,p7)
        l10 = gmsh.model.geo.add_line(p7,p8)
        groove_points = Int[p8,p8_up] # Matriz para almacenar el orden de los puntos de las ranuras

        i = 1
        while i + 1 <= length(points_ranuras_inf)
        push!(groove_points, points_ranuras_sup[i])   # tramo plano arriba
        push!(groove_points, points_ranuras_inf[i])   # baja a la pared
        push!(groove_points, points_ranuras_inf[i+1]) # tramo plano abajo   
        push!(groove_points, points_ranuras_sup[i+1]) # sube a la ranura
        i += 2
        end
        #println("groove_points: ", groove_points)
        # Si 2n-1 es impar, sobra un índice
        if i <= length(points_ranuras_inf)
            push!(groove_points, points_ranuras_sup[i])
            push!(groove_points, points_ranuras_inf[i])
        end

        lines_ranuras = Int[]              # Matriz para almacenar las líneas de las ranuras
        for k in 1:length(groove_points)-1 # Generar las líneas de las ranuras
            push!(lines_ranuras, gmsh.model.geo.add_line(groove_points[k], groove_points[k+1]))
        end

        l11 = gmsh.model.geo.add_line(groove_points[end], p9)
        l12 = gmsh.model.geo.add_line(p9,p10)
        l13 = gmsh.model.geo.add_line(p10,p11) #INLET
        l14 = gmsh.model.geo.add_line(p11,p12)
        l15 = gmsh.model.geo.add_line(p12,p13)
        l16= gmsh.model.geo.add_line(p13,p14) #INLET

        #Create loop and surface
        loop = gmsh.model.geo.add_curve_loop(
            [l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,lines_ranuras...,l11,l12,l13,l14,l15,l16])

        surf = gmsh.model.geo.add_plane_surface([loop])

        #Synchronize model
        gmsh.model.geo.synchronize()

        #Create physical domains (for constrains)

        #Paredes
        gmsh.model.add_physical_group(1,[l1,l2,l3,l4], -1, "CH1")       #CH1 es la pared 1 (Inferior)
        gmsh.model.add_physical_group(1,[l6,l7], -1, "CH2")             #CH2 es la pared 2 (Derecha)
        gmsh.model.add_physical_group(1,[l9,l10,lines_ranuras...,l11,l12], -1, "CH3")    #CH3 es la pared 3 (Superior)
        gmsh.model.add_physical_group(1,[l14,l15], -1, "CH4")           #CH4 es la pared 4 (Izquierda)

        #Inlets
        gmsh.model.add_physical_group(1,[l13], -1, "IN1") #IN1 es Inlet 1 (superior) - inlet
        gmsh.model.add_physical_group(1,[l16], -1, "IN2") #IN2 es Inlet 2 (inferior) - buffer

        #Outlets
        gmsh.model.add_physical_group(1,[l5], -1, "OUT1") #OUT1 es Outlet 1 (superior) - waste
        gmsh.model.add_physical_group(1,[l8], -1, "OUT2") #OUT2 es Outlet 2 (inferior) - outlet

        gmsh.model.add_physical_group(2,[surf])

        #Generate 2d Mesh
        gmsh.model.mesh.generate(2)

        #Save mesh and read back in Ferrite Grid

        grid = mktempdir() do dir
            path = joinpath(dir, "mesh.msh")
            gmsh.write(path)
            togrid(path)
        end
        #= Eliminar Comment para imprimir ventana de GMSH
        if !("-nopopup" in ARGS)
            gmsh.fltk.run()
        end
        
        Gmsh.finalize()
        =#
        return grid

    end



end