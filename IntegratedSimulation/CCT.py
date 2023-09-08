def runCCT():
    import h5py
    import matplotlib.pyplot as plt
    import numpy as np
    from HelpFile import read_input
    print('CCT module')
    data = read_input()
    with h5py.File("Resultfiles/ThermoCache.hdf5", "r") as f:
        x = np.array(f.get("CNcurves/Position"))
        a = dict()
        for element in data['Material']['Composition'].keys():
            plt.plot(x, 100 * np.array(f.get("CNcurves/"+element)))
            a[element] = 100 * np.array(f.get("CNcurves/" + element))[-1]

    #Ccurve = np.array(f.get("CNcurves/C"))
    #Ncurve = np.array(f.get("CNcurves/N"))
    #plt.plot(Ccurve)
    #plt.plot(Ncurve)
    plt.legend(data['Material']['Composition'].keys(), loc='upper left')
    plt.show()
    print(data['Material']['Composition'])
    print(a)
    #print(Ccurve.min())
    #print(Ncurve.min())
    #print(Ccurve.max())
    #print(Ncurve.max())