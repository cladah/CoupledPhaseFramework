from IntegratedSimulation.Solvers.Thermocalc import *
from HelpFile import *
import h5py

def runcarbonitriding():

    data = read_input()
    if data["Programs"]["CNDiffusion"] == "TC":
        print('Running carbon-nitriding module with ThermoCalc')
        activityenv = TCequalibrium("env")
        CN = TCcarbonitriding(activityenv)
    else:
        raise KeyError('Program not implemented for carbonitriding')

    with h5py.File("Resultfiles/ThermoCache.hdf5", "w") as f:
        for element in CN[1].keys():
            f.create_dataset("CNcurves/"+element, data=np.array(CN[1][element]))
        f.create_dataset("CNcurves/Position", data=np.array(CN[0]))

