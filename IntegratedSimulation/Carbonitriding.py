from IntegratedSimulation.Solvers.Thermocalc import *
from HelpFile import *
import h5py

def runcarbonitriding():
    print('Carbonitriding module')
    data = read_input()
    if data["Programs"]["CNDiffusion"] == "TC":
        print('Running carbon-nitriding module with ThermoCalc')
        activityenv = TCequalibrium("env")
        CN = TCcarbonitriding(activityenv)
    else:
        raise KeyError(str(data["Programs"]["CNDiffusion"]) + ' not implemented for carbonitriding')

    with h5py.File("Resultfiles/ThermoResult.hdf5", "w") as f:
        for element in CN[1].keys():
            f.create_dataset("CNcurves/"+element, data=np.array(CN[1][element]))
        f.create_dataset("CNcurves/Position", data=np.array(CN[0]))

