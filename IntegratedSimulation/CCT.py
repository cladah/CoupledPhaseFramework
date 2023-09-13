from IntegratedSimulation.Solvers.Thermocalc import *
def modelfitting(model, x, y):
    import scipy
    import numpy as np
    from scipy.interpolate import BSpline, make_interp_spline
    if model == "KM":
        def objective(x, Ms, beta):
            return 1-np.exp(-beta*(Ms-x))
    elif model == "JMAK":
        def objective(x, tau, n):
            return 1 - np.exp(-(x/tau)**n)
    elif model == "Devol":
        def objective(x, a, b, c):
            return a * x + b
    elif model == "Devol":
        def objective(x, K, L):
            return K * x + L * (1-x)
    else:
        raise KeyError("Fitting model wrong or not implemented")
    fitvalues = scipy.optimize.curve_fit(objective, x, y)
    spline = make_interp_spline(y, x, k=3, t=None, bc_type=None, axis=0, check_finite=True)

    return spline


def runCCT():
    import h5py
    import matplotlib.pyplot as plt
    import numpy as np
    from HelpFile import read_input
    print('CCT module')
    data = read_input()
    with h5py.File("Resultfiles/ThermoResult.hdf5", "r") as f:
        x = np.array(f.get("CNcurves/Position"))
        a = dict()
        for element in data['Material']['Composition'].keys():
            #plt.plot(x, 100 * np.array(f.get("CNcurves/"+element)))
            a[element] = 100 * np.array(f.get("CNcurves/" + element)) # Composition curves at all points along x

    #Ccurve = np.array(f.get("CNcurves/C"))
    #Ncurve = np.array(f.get("CNcurves/N"))
    #plt.plot(Ccurve)
    #plt.plot(Ncurve)
    #plt.legend(data['Material']['Composition'].keys(), loc='upper left')
    #plt.show()
    x = [239.7, 126.4, 73.5, 53.0, 55.6, 67.6, 68.1, 90.6, 185.0, 838.0]
    y = [780, 800, 820, 840, 860, 880, 900, 920, 940, 960]
    Perlite = modelfitting('JMAK', x, y)
    print(Perlite)
    T_new = np.linspace(700,1000,50)
    t_new = Perlite(T_new).T

    #plt.plot(x, y)
    #plt.plot(t_new,T_new)
    #plt.show()
    #Martensite = modelfitting(data['Material']['Martensite']['model'], [1, 1], 2)
    #Perlite = modelfitting(data['Material']['Perlite']['model'], [1, 1], [1, 1])
    #Bainite = modelfitting(data['Material']['Bainite']['model'], [1, 1], [1, 1])

    with h5py.File("Resultfiles/ThermoResult.hdf5", "r+") as f:
        try:
            f.create_dataset('JMAK/Tau', data=1)
            f.create_dataset('JMAK/n', data=1)
        except ValueError:
            del f['JMAK/Tau']
            del f['JMAK/n']
            f.create_dataset('JMAK/Tau', data=1)
            f.create_dataset('JMAK/n', data=1)