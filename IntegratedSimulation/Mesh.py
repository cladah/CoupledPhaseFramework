from HelpFile import read_input
import numpy as np
def createMesh():
    print('Mesh module')
    data = read_input()
    if data['FEM']["elementtype"] == "Spherical":
        print('Creating 1D mesh...')
        mesh1D()
    elif data['FEM']["elementtype"] == "Axisym2D":
        print('Creating 2D mesh...')
        mesh2D()
    else:
        raise KeyError('Meshtype implemented in mesh module')

def mesh1D():
    import gmsh
    data = read_input()
    lc = data['Geometry']['radius'] * data['FEM']['Meshscaling'] ** (data['FEM']['nodes'] - 1) / np.sum(
        data['FEM']['Meshscaling'] **
        np.linspace(0, data['FEM']['nodes'] - 1, data['FEM']['nodes']))
    gmsh.initialize()
    model = gmsh.model()

    factory = model.geo
    factory.addPoint(0, 0, 0, lc, 1)
    factory.addPoint(data['Geometry']['radius'], 0, 0, lc, 2)
    factory.addLine(1, 2, 1)
    factory.synchronize()
    gmsh.model.addPhysicalGroup(1, [1], 1, 'Line')
    gmsh.model.mesh.set_transfinite_curve(1, data['FEM']['nodes'], 'Progression', data['FEM']['Meshscaling'])
    gmsh.model.mesh.generate(1)

    gmsh.write("Resultfiles/Mesh.msh")
    gmsh.write("Resultfiles/Mesh.nas")
    gmsh.finalize()

def mesh2D():
    import gmsh
    data = read_input()
    lc = data['Geometry']['radius'] * data['FEM']['Meshscaling'] ** (data['FEM']['nodes'] - 1) / np.sum(data['FEM']['Meshscaling'] **
        np.linspace(0, data['FEM']['nodes'] - 1, data['FEM']['nodes']))

    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("QuarterCirc")
    gdim = 2

    gmsh.model.occ.addPoint(0, 0, 0, lc, 1)
    gmsh.model.occ.addPoint(data['Geometry']['radius'], 0, 0, lc, 2)
    gmsh.model.occ.addPoint(0, data['Geometry']['radius'], 0, lc, 3)
    gmsh.model.occ.addLine(1, 2, 1)
    gmsh.model.occ.addLine(1, 3, 2)
    gmsh.model.occ.addCircleArc(2, 1, 3, 3)
    gmsh.model.occ.addCurveLoop([1, 2, 3], 4)
    gmsh.model.occ.addPlaneSurface([4], 5)
    gmsh.model.occ.synchronize()

    gmsh.model.addPhysicalGroup(1, [1], 1, 'Bottom')
    gmsh.model.addPhysicalGroup(1, [2], 2, 'Side')
    gmsh.model.addPhysicalGroup(1, [3], 3, 'Circumference')
    gmsh.model.addPhysicalGroup(gdim, [5], 4, 'Sphere')
    gmsh.model.mesh.set_transfinite_curve(1, data['FEM']['nodes'], 'Progression', data['FEM']['Meshscaling'])
    gmsh.model.mesh.set_transfinite_curve(2, data['FEM']['nodes'], 'Progression', data['FEM']['Meshscaling'])
    #gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 2*data['Geometry']['radius']/100)
    #gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 3*data['Geometry']['radius']/100)
    gmsh.model.mesh.generate(gdim)
    # ----------------------
    gmsh.write("Resultfiles/Mesh.msh")
    gmsh.write("Resultfiles/Mesh.nas")
    gmsh.finalize()