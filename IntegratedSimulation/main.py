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

if __name__ == "__main__":
    start()