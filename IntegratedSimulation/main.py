from Inputfile import createmodelinput
from Mesh import createMesh
#from Carbonitriding import *
from RunDocker import *


def start():

    modelvar = createmodelinput()
    # --------------- Input ------------------#
    #inputvariable = input('What do you want to do? ')
    inputvariable = 'Run'
    if inputvariable == 'Run':
        createMesh(modelvar)
        rundocker(modelvar)
        # runcarbonitriding(modelvar)
        # runquenching(modelvar)
        # runtempering(modelvar)
        # runfatiguetest(modelvar)
        # runplot(modelvar)

    elif inputvariable == "Test":
        rundocker(modelvar)

if __name__ == "__main__":
    start()