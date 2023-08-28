from HelpFile import read_input
def createMesh():
    data = read_input()
    if data['FEM']["elementtype"] == "Spherical":
        mesh1D(0,data['Geometry']['radius'])
    elif data['FEM']["elementtype"] == "Axisym2D":
        mesh2D(0, data['Geometry']['radius'])
    else:
        raise KeyError('Not implemented in mesh')

def mesh1D(r1, r2):
    import gmsh

    meshsize = (r2-r1)/100
    gmsh.initialize()
    model = gmsh.model()

    factory = model.geo
    factory.addPoint(0, 0, 0, meshsize, 1)
    factory.addPoint(r2, 0, 0, meshsize, 2)
    factory.addLine(1, 2, 1)
    factory.synchronize()

    gmsh.model.mesh.generate(1)

    gmsh.write("Resultfiles/Mesh.msh")
    gmsh.write("Resultfiles/Mesh.nas")
    gmsh.finalize()

def mesh2D(r1, r2):
    import gmsh

    gmsh.initialize()
    gmsh.model.add("QuarterCirc")
    gdim = 2
    gmsh.model.occ.addDisk(0, 0, 0, r2, 1)
    gmsh.model.occ.addRectangle(0, 0, 0, r2, r2, 2)
    gmsh.model.occ.intersect([(gdim, 1)], [(gdim, 2)], 3)
    gmsh.model.occ.synchronize()

    gmsh.model.addPhysicalGroup(gdim, [3], 1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.02)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.03)
    gmsh.model.mesh.generate(gdim)
    # ----------------------
    gmsh.write("Resultfiles/Mesh.msh")
    gmsh.write("Resultfiles/Mesh.nas")
    gmsh.finalize()

    #gmsh_model_rank = 0
    #mesh_comm = MPI.COMM_WORLD
    #domain, cell_markers, facet_markers = gmshio.model_to_mesh(gmsh.model, mesh_comm, gmsh_model_rank, gdim=gdim)