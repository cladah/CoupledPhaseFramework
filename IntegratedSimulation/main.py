from Inputfile import createmodelinput
from Quenching import runquenching
#from Carbonitriding import runcarbonitriding
from Mesh import createMesh

# Test imports
from Solvers.ComsolSolver import runComsol
from Solvers.RunDocker import rundocker


def start():

    modelvar = createmodelinput()
    # --------------- Input ------------------#
    #inputvariable = input('What do you want to do? ')
    inputvariable = 'Test'
    if inputvariable == 'Run':
        createMesh(modelvar)
        #runcarbonitriding(modelvar)
        runquenching(modelvar)
        # runtempering(modelvar)
        # runfatiguetest(modelvar)
        # runplot(modelvar)



    elif inputvariable == "Test":
        createMesh(modelvar)
        runComsol()
        #rundocker(modelvar)
        #with h5py.File('C:/Users/ClasD/Documents/GitHub/CoupledPhaseFramework/IntegratedSimulation/Resultfiles/displacement.h5', "r") as f:
        #    disp = f['Function']['Displacement']['0'][...]
        #    print(disp)
            #temp = f['Function']['Temperature']['0'][...]
            #print(temp)


if __name__ == "__main__":
    start()