

def meshnew():
    import gmsh

    r = 1

    gmsh.initialize()
    model = gmsh.model()
    model.add("circle")

    factory = model.geo

    lcent = r/100
    lcirc = r/100

    #membrane = model.occ.addDisk(0, 0, 0, 1, 1)

    factory.addPoint(0,0,0,lcent,1)
    factory.addPoint(r,0,0,lcirc,2)
    factory.addPoint(0,r,0,lcirc,3)

    factory.addLine(3,1,1)
    factory.addLine(1,2,2)
    factory.addCircleArc(2,1,3,3)

    factory.addCurveLoop([1,2,3],4)
    factory.addPlaneSurface([4],5)
    factory.synchronize()
    #model.addPhysicalGroup(2, [membrane], 1)
    #gmsh.option.setNumber("Mesh.CharacteristicLengthMin",0.05)
    #gmsh.option.setNumber("Mesh.CharacteristicLengthMax",0.05)
    model.mesh.generate(2)
    gmsh.write("circle.msh")
    #domain, cell_markers, facet_markers = gmshio.model_to_mesh(model, MPI.COMM_WORLD, 0)