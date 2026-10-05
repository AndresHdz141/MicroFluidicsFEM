module grid_mag
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays



export magnet_setup_grid

    function magnet_setup_grid(L, height,n_imanes, array, h=0.05)

        a = L/2 
        paso = L / n_imanes
        lado = 0.6 * paso

        Gmsh.initialize()
        gmsh.option.set_number("General.Verbosity", 2)

        #Coordenadas de puntos para el aire
        #X
        x1 = 2.2*a
        #Y 
        y1 = a*(3/2)

        #coordenadas para los imanes

        x_init = -a  + (paso-lado)/2
        y_init = -height/8 -lado

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

        println("Check 1 - Grid")
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
        println("Check 2 - Grid")

        #surfs

        surf = gmsh.model.geo.add_plane_surface([loop,lopps...])


        if array == 1     
            surfs_imp = Int[]   #los surfs de los imanes
            surfs_par = Int[]

            for i in 1:n_imanes
            srf = gmsh.model.geo.add_plane_surface([lopps[i]])
            if i%2 == 0
            push!(surfs_par, srf)
            else
            push!(surfs_imp, srf)
            end
        end

        else
            surfs_1 = Int[]
            surfs_2 = Int[]
            surfs_3 = Int[]
            surfs_4 = Int[]

            for num in 1:n_imanes
                srf = gmsh.model.geo.add_plane_surface([lopps[num]])
                surf_action = [
                (x) -> push!(surfs_1, x),
                (x) -> push!(surfs_2, x),
                (x) -> push!(surfs_3, x),
                (x) -> push!(surfs_4, x)
                ]
                function push_surf(i)
                    surf_action[mod1(i, 4)](srf)
                end
                push_surf(num)
            end     
        end

        #Synchronize model
        gmsh.model.geo.synchronize()

        #Create physical domains (for constrains)

        #Paredes
        gmsh.model.add_physical_group(1,[l1,l2,l3,l4], -1, "AIR_Lm")       #Aire
        
        #area de Aire

        gmsh.model.add_physical_group(2,[surf], -1, "AIR")

        #Imanes

        if array == 1

            gmsh.model.add_physical_group(2, surfs_imp, -1, "magnets_impar")
            gmsh.model.add_physical_group(2, surfs_par, -1, "magnets_par")

        else
            gmsh.model.add_physical_group(2, surfs_1, -1, "magnets_1")
            gmsh.model.add_physical_group(2, surfs_2, -1, "magnets_2")
            gmsh.model.add_physical_group(2, surfs_3, -1, "magnets_3")
            gmsh.model.add_physical_group(2, surfs_4, -1, "magnets_4")
        end


            


        #Generate 2d Mesh
        gmsh.model.mesh.generate(2)
        println("Check 3 - Grid")

        #Save mesh and read back in Ferrite Grid

        grid = mktempdir() do dir
            path = joinpath(dir, "meshmag.msh")
            gmsh.write(path)
            togrid(path)
        end
           
        #=
        if !("-nopopup" in ARGS)
            gmsh.fltk.run()
        end
        =#
        Gmsh.finalize()

        return  grid

    end    

end