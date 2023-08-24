from IntegratedSimulation.Solvers.Thermocalc import *
from HelpFile import *


def runcarbonitriding(modelvar):
    if modelvar["programs"]["Diff"] == "TC":
        activityenv = TCequalibrium(modelvar,"env")
        #activitymat = TCequalibrium(modelvar,"steel")
        CN = TCcarbonitriding(modelvar,activityenv)
        concentrationC = CN[0]
        concentrationN = CN[1]
    else:
        raise KeyError('Program not implemented for carbonitriding')

    savetocache("CNcurves/C", np.array(concentrationC))
    savetocache("CNcurves/N", np.array(concentrationN))

