from Inputfile import createmodelinput
from Mesh import meshnew
from Quenching import *
#from Carbonitriding import *
from Post_processing import *
import subprocess
import docker
import h5py
from Solvers.RunDocker import *
import time

def start():
    modelvar = createmodelinput()
    # --------------- Input ------------------#
    inputvariable = input('What do you want to do? ')

    if inputvariable == 'Run':
        meshnew()
        # createmesh(modelvar)
        # runcarbonitriding(modelvar)
        # runquenching(modelvar)
        # runtempering(modelvar)
        # runfatiguetest(modelvar)
        # runplot(modelvar)

    elif inputvariable == '1D':
        if not testchacheinfo(modelvar):
            print("Modelinput changed recalculation matrixes")

        RHS = 1
        #sphericalsolver(RHS, modelvar)

        # --------------- Adding new info to cache ------------------#
        savetocache("info",
                    [modelvar["radius"], modelvar["nodesFEM"], modelvar["E"], modelvar["nu"], modelvar["shapef"],
                     modelvar["Meshscaling"]])

    elif inputvariable == "Test":
        createresultfile(modelvar)
        createcachfile(modelvar)
        createmesh(modelvar)

        #[displacement, err] = sphericalsolver(1,modelvar)
        saveresult("displacement",displacement)
        plotstrain(modelvar)
        plotdisplacement(modelvar)
        #strains = getstrains(modelvar, results['displacement'])
        #sigma = getstresses(modelvar, strains)
        #createResults(modelvar, results['displacement'])
        #plotstrain(modelvar)

    elif inputvariable == 'y':
        rundocker(modelvar)

if __name__ == "__main__":
    start()