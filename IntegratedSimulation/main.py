from Inputfile import createmodelinput
from Quenching import runquenching
#from Carbonitriding import runcarbonitriding
from Mesh import createMesh
from HelpFile import read_input

# Test imports
from Solvers.ComsolSolver import runComsol
from Solvers.RunDocker import rundocker
import json
import h5py


def start():
    # --------------- Input ------------------#
    #inputvariable = input('What do you want to do? ')
    inputvariable = 'Test'
    if inputvariable == 'Run':
        createMesh()
        #runcarbonitriding(modelvar)
        runquenching()
        # runtempering(modelvar)
        # runfatiguetest(modelvar)
        # runplot(modelvar)



    elif inputvariable == "Test":
        #createMesh()
        #runComsol()
        rundocker()
        #with h5py.File('C:/Users/ClasD/Documents/GitHub/CoupledPhaseFramework/IntegratedSimulation/Resultfiles/displacement.h5', "r") as f:
        #    disp = f['Function']['Displacement']['0'][...]
        #    print(disp)
        #    temp = f['Function']['Temperature']['0'][...]
        #    print(temp)

if __name__ == "__main__":
    start()