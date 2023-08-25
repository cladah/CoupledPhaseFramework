import mph
import os
def runComsol():
    directory = os.getcwd()
    directory = directory.replace('\\Solvers','\\Resultfiles')
    meshdirec = directory + '\\Mesh.nas'

    client = mph.start()
    model = client.create('Sphere')
    model.component('Sphere').geom().create('geo1', 2)
    model.component('Sphere').geom('geo1').isaxisymmetric()
    #model.mesh().create("mesh1", 2)
    #model.geom("geom1").feature().create("blk1", "Block");
    #model.geom("geom1").feature("blk1").set("size", ["0.1", "0.2", "0.5"]);
    #model.geom("geom1").run("fin");

    #model.mesh().run()
    #model.solve('static')
    model.save('Resultfiles/model')