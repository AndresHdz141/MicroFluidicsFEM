module grid_mag
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays



export magnet_setup_grid

    function magnet_setup_grid(L, height, h=0.05)

        a = L/2
        n_imanes =  9
        paso = L / n_imanes
        lado = 0.6 * paso

        Gmsh.initialize()
        gmsh.option.set_number("General.Verbosity", 2)

        #Coordenadas de puntos para el aire
        #X
        x1 = 2*a
        #Y 
        y1 = a*(3/2)

        #coordenadas para los imanes

        x_init = -a  + (paso-lado)/2
        y_init = -height/2 -lado

        #crea los puntos de los imanes partiendo de unas cordenadas iniciales
        imanes= Vector{Vector{Int}}()
        for i in 1:n_imanes
            x = x_init + (i-1) * paso
            y = y_init

            puntos = [
                gmsh.model.geo.add_point(x,y,0.0,h),
                gmsh.model.geo.add_point(x+lado,y,0.0,h),
                gmsh.model.geo.add_point(x+lado,y+lado,0.0,h),
                gmsh.model.geo.add_point(x,y+lado,0.0,h)
            ]
            push!(imanes, puntos)

        end
       
        #Add points
        #Origin
        #o = gmsh.model.geo.add_point(0.0, 0.0,0.0, h)

        #superficie del aire
        p1 = gmsh.model.geo.add_point(x1, y1,0.0, h)
        p2 = gmsh.model.geo.add_point(-x1, y1,0.0, h)
        p3 = gmsh.model.geo.add_point(-x1, -y1,0.0, h)
        p4 = gmsh.model.geo.add_point(x1, -y1,0.0, h)

       
        #tecnicamente aquí los inlets y outlets son los imanes pero luego hago eso

        l1 = gmsh.model.geo.add_line(p1,p2)
        l2 = gmsh.model.geo.add_line(p2,p3)
        l3 = gmsh.model.geo.add_line(p3,p4)
        l4 = gmsh.model.geo.add_line(p4,p1) 

     
        #un for que crea las linas y las cuarda en vectores por grupos
        lineas_imanes = Vector{Vector{Int}}()

        for i in 1:n_imanes 
            puntos_act = imanes[i]
            lineas = Int[]
            for j in 1:3
                push!(lineas,gmsh.model.geo.add_line(puntos_act[j],puntos_act[j+1]))
            end

            #la cuarta linea
            push!(lineas, gmsh.model.geo.add_line(puntos_act[4],puntos_act[1]))

            push!(lineas_imanes, lineas)
        end


        #Create loop and surface
        loop = gmsh.model.geo.add_curve_loop([l1,l2,l3,l4]) #el del aire

        lopps = Int[]   #los loops de los imanes

        #Aqui es un for que crea los lopps
        for i in 1:n_imanes
            lin = lineas_imanes[i]
            
            loopp = gmsh.model.geo.add_curve_loop([lin[1], 
                                                   lin[2],
                                                   lin[3],
                                                   lin[4]])      
            push!(lopps, loopp)
        end

        #surfs

        surf = gmsh.model.geo.add_plane_surface([loop,lopps...])

        surfs = Int[]   #los surfs de los imanes
        for i in 1:n_imanes
            srf = gmsh.model.geo.add_plane_surface([lopps[i]])
            push!(surfs, srf)
        end

        #Synchronize model
        gmsh.model.geo.synchronize()

        #Create physical domains (for constrains)

        #Paredes
        gmsh.model.add_physical_group(1,[l1,l2,l3,l4], -1, "AIR_Lm")       #Aire
        
        #Imanes

        gmsh.model.add_physical_group(2,[surf], -1, "AIR")
        gmsh.model.add_physical_group(2, surfs, -1, "magnets")


        #Generate 2d Mesh
        gmsh.model.mesh.generate(2)

        #Save mesh and read back in Ferrite Grid

        grid = mktempdir() do dir
            path = joinpath(dir, "meshmag.msh")
            gmsh.write(path)
            togrid(path)
        end
            
        if !("-nopopup" in ARGS)
            gmsh.fltk.run()
        end
        
        Gmsh.finalize()

        return grid

    end    

end