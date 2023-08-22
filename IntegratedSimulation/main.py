from Inputfile import createmodelinput
from Mesh import createMesh
#from Carbonitriding import *
from RunDocker import *
import h5py

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

        with h5py.File('C:/Users/ClasD/Documents/GitHub/CoupledPhaseFramework/IntegratedSimulation/Resultfiles/displacement.h5', "r") as f:
            disp = f['Function']['Displacement']['0'][...]
            print(disp)
            temp = f['Function']['Temperature']['0'][...]
            print(temp)

    elif inputvariable == "Test":
        rundocker(modelvar)



if __name__ == "__main__":
    start()