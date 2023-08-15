import numpy as np


def createMesh(modelvar):
    if modelvar["elementtype"] == "Spherical":
        mesh1D(0, modelvar['radius'])
    elif modelvar["elementtype"] == "Axisym2D":
        mesh2D(0, modelvar['radius'])
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
    gmsh.write("Resultfiles/Mesh.vtk")
    gmsh.finalize()

def mesh2D(r1, r2):
    import gmsh

    gmsh.initialize()
    model = gmsh.model()
    model.add("circle")
    model.set_current('circle')

    factory = model.geo

    lcent = r2 / 50
    lcirc = r2 / 50

    # membrane = model.occ.addDisk(0, 0, 0, 1, 1)

    factory.addPoint(0, 0, 0, lcent, 1)
    factory.addPoint(r2, 0, 0, lcirc, 2)
    factory.addPoint(0, r2, 0, lcirc, 3)
    factory.addLine(3, 1, 1)
    factory.addLine(1, 2, 2)
    factory.addCircleArc(2, 1, 3, 3)

    factory.addCurveLoop([1, 2, 3], 4)
    factory.addPlaneSurface([4], 1)
    factory.synchronize()
    model.addPhysicalGroup(2, [1], name="My surface")
    #surface = []
    #for surface in gmsh.model.getEntities(dim=2):
    #    surfaces = surface[1]
    #gmsh.model.addPhysicalGroup(2, [surfaces], 1)

    #outside = []
    #for line in gmsh.model.getEntities(dim=1):
    #    com = gmsh.model.occ.getCenterOfMass(line[0], line[1])
    #    if np.isclose(np.sqrt(com[0]**2+com[1]**2), r2):
    #        outside.append(line[1])
    #gmsh.model.addPhysicalGroup(1, outside, 2)

    # model.addPhysicalGroup(2, [membrane], 1)
    # gmsh.option.setNumber("Mesh.CharacteristicLengthMin",0.05)
    # gmsh.option.setNumber("Mesh.CharacteristicLengthMax",0.05)
    model.mesh.generate(2)
    gmsh.write("Resultfiles/Mesh.msh")
    gmsh.write("Resultfiles/Mesh.vtk")
    gmsh.finalize()
    # domain, cell_markers, facet_markers = gmshio.model_to_mesh(model, MPI.COMM_WORLD, 0)
