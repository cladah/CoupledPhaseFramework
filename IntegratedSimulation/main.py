from Inputfile import createmodelinput
from Quenching import runquenching
from Carbonitriding import runcarbonitriding
from CCT import runCCT
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
    import matplotlib.pyplot as plt
    # --------------- Input ------------------#
    #inputvariable = input('What do you want to do? ')
    inputvariable = 'Run'
    if inputvariable == 'Run':
        # createresultfile()
        createMesh()
        runcarbonitriding()
        runCCT()
        # TTTfit()
        # runquenching()
        # runtempering()
        # runfatiguetest()
        # runplot()
        print(readresultfile("TTT/Surface/Bainite/Start"))
        print(readresultfile("TTT/Surface/Bainite/Half"))
        print(readresultfile("TTT/Surface/Bainite/Finish"))
        plt.plot(readresultfile("TTT/Surface/Bainite/Start"),np.linspace(300,1000,10))
        plt.xscale("log")
        plt.show()
        return


        print(readresultfile("JMAK_perlite"))
        Tbai,nbai,taubai = readresultfile("JMAK_bainite")
        Tper, nper, tauper = readresultfile("JMAK_perlite")
        xpoints = -tauper * (-np.log(0.98)) ** (1 / nper)
        ypoints = Tper
        if 1==1:
            plt.plot(xpoints, ypoints)
            xpoints = -tauper * (-np.log(0.02)) ** (1 / nper)
            ypoints = Tper
            plt.plot(xpoints, ypoints)
            xpoints = -taubai * (-np.log(0.98)) ** (1 / nbai)
            xpoints = np.append(xpoints, 5000 * (-np.log(0.98)) ** (1 / 4))
            ypoints = np.append(Tbai, 1000)
            plt.plot(xpoints, ypoints)
            #xpoints = -taubai * (-np.log(0.02)) ** (1 / nbai)
            #ypoints = Tbai
            #plt.plot(xpoints, ypoints)
            plt.xscale("log")
            plt.show()
        plt.plot(Tper,tauper)
        plt.show()
        print('Simulation done')
        print('Caching data')
        #createinputcache()


    elif inputvariable == "Test":
        createresultfile()
        createMesh()
        runcarbonitriding()
        # runCCT()
        TTTfit()
        # runquenching()
        # runtempering()
        # runfatiguetest()
        # runplot()
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