from Inputfile import createmodelinput
from Quenching import runquenching
from Carbonitriding import runcarbonitriding
#from CCT import runCCT
from Solvers.TTTmodelfit import TTTfit
from Mesh import createMesh
from Post_processing import runplot
from HelpFile import *

# Test imports
#from Solvers.ComsolSolver import runComsol
#from Solvers.RunDocker import rundocker
#import json
#import h5py
#from Post_processing import *


def start():
    # --------------- Input ------------------#
    #inputvariable = input('What do you want to do? ')
    inputvariable = 'Run'
    if inputvariable == 'Run':
        createresultfile()
        createMesh()
        runcarbonitriding()
        # runCCT()
        TTTfit()
        # runquenching()
        # runtempering()
        # runfatiguetest()
        # runplot()
        print(readresultfile("JMAK_perlite"))
        print('Simulation done')
        print('Caching data')
        #createinputcache()


    elif inputvariable == "Test":
        pass
        #data = readCNfile()
        #createMesh()
        #runComsol()
        #rundocker()
        #runpostprocess()
        #plotPyVista()
        #with h5py.File('C:/Users/ClasD/Documents/GitHub/CoupledPhaseFramework/IntegratedSimulation/Resultfiles/displacement.h5', "r") as f:
        #    disp = f['Function']['Displacement']['0'][...]
        #    print(disp)
        #    temp = f['Function']['Temperature']['0'][...]
        #    print(temp)

if __name__ == "__main__":
    start()