import mph
import os
def runComsol():
    directory = os.getcwd()
    directory = directory.replace('\\Solvers','\\Resultfiles')

    client = mph.start()
    pymodel = client.create('Model')
    model = pymodel.java

    model.component().create("comp1")
    g = model.component('comp1').geom().create('geom1', 3)
    m = model.component('comp1').mesh().create('Import')
    imp1 = m.create('imp1', 'Import')
    imp1.set('type', 'mesh')
    imp1.set('mesh', 'mesh1')
    m.run()
    #model.modelNode().create("sphere")
    #model.geom().create("Sphere", 2)
    #model.modelNode().mesh().create("Import")
    #model.component().mesh('Import').feature().set()
    #model.component().mesh('Import').feature().getType()
    #model.component().mesh('Import').feature().importData(directory)
    #model.geom().mesh().create().import_(directory + 'Resultfiles/Mesh.nas')
    #model.geom().feature().importData()

    #model.mesh().run()
    #model.solve('static')
    model.save('Resultfiles/model')