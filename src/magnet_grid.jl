module grid_mag
using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays



export magnet_setup_grid

    function magnet_setup_grid(length, height, h=0.05)
        a = length/2
        b = height/2
        g = length/9
        Gmsh.initialize()
        gmsh.option.set_number("General.Verbosity", 2)

        #Coordenadas de puntos
        #X
        x1 = 2*a
        x2 = a
        x3 = a
        x4 = a-g
        x5 = a-g
        #Y 
        y1 = a*(3/2)
        y2 = -b
        y3 = -b-g
        y4 = -b-g
        y5 = -b
        

       
        #Add points
        #Origin
        #o = gmsh.model.geo.add_point(0.0, 0.0,0.0, h)

        #superficie del aire
        p1 = gmsh.model.geo.add_point(x1, y1,0.0, h)
        p2 = gmsh.model.geo.add_point(-x1, y1,0.0, h)
        p3 = gmsh.model.geo.add_point(-x1, -y1,0.0, h)
        p4 = gmsh.model.geo.add_point(x1, -y1,0.0, h)

        p5 = gmsh.model.geo.add_point(x2, y2,0.0, h)
        p6 = gmsh.model.geo.add_point(x3, y3,0.0, h)
        p7 = gmsh.model.geo.add_point(x4, y4,0.0, h)
        p8 = gmsh.model.geo.add_point(x5, y5,0.0, h)
        

        #tecnicamente aquí los inlets y outlets son los imanes pero luego hago eso

        l1 = gmsh.model.geo.add_line(p1,p2)
        l2 = gmsh.model.geo.add_line(p2,p3)
        l3 = gmsh.model.geo.add_line(p3,p4)
        l4 = gmsh.model.geo.add_line(p4,p1) 

        l5 = gmsh.model.geo.add_line(p5,p6) 
        l6 = gmsh.model.geo.add_line(p6,p7)
        l7 = gmsh.model.geo.add_line(p7,p8)
        l8 = gmsh.model.geo.add_line(p8,p5)

        #Create loop and surface
        loop = gmsh.model.geo.add_curve_loop([l1,l2,l3,l4])
        loop2 =gmsh.model.geo.addCurveLoop([l5,l6,l7,l8])

        surf = gmsh.model.geo.add_plane_surface([loop,loop2])
        surf2 = gmsh.model.geo.add_plane_surface([loop2])

        #Synchronize model
        gmsh.model.geo.synchronize()

        #Create physical domains (for constrains)

        #Paredes
        gmsh.model.add_physical_group(1,[l1,l2,l3,l4], -1, "Air")       #CH1 es la pared
        gmsh.model.add_physical_group(1,[l5,l6,l7,l8],-1, "magnet")
        #Inlets
    
        #Outlets
       

        gmsh.model.add_physical_group(2,[surf])
        gmsh.model.addPhysicalGroup(2,[surf2])

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