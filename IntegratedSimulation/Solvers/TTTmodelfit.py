def TTTfit():
    print("TTTFitting module")
    from HelpFile import saveresult
    import csv
    import scipy
    import numpy as np
    from scipy.interpolate import BSpline, make_interp_spline, interpn
    import matplotlib.pyplot as plt

    # 'Ferrite start (2%)', 'Pearlite start (2%)', 'Bainite start (2%)', 'Total ferrite start (2%)',
    # 'Total ferrite+cementite start (2%)', 'Austenite transformed 2%', 'Austenite transformed 50%',
    # 'Austenite transformed 98%', 'Martensite start', 'Martensite 50%', 'Martensite 98%', 'Temperature [K]'

    with open("Resultfiles/TTT_Si0,4_Mn1,2_P0,025_S0,014_Cr1,15_C0,2_N0,01_all.txt") as f:
        freader = csv.reader(f, delimiter='\t')
        header = freader.__next__()
        numlines = len(header)
        TTTdata = list()
        for line in freader:
            TTTdata.append(line)
        TTTdata = np.array(TTTdata)

        for i in range(0,len(header)):
            if header[i]=='Temperature [K]':
                Tnum = i
        Tdata = TTTdata[:,Tnum]
        data = dict()
        for x in header:
            listnum = header.index(x)
            tmpdata = np.array([TTTdata[:,listnum],Tdata])
            deleteindex = [i for i, x in enumerate(TTTdata[:,listnum]) if x == '']
            tmpdata = np.delete(tmpdata,deleteindex,1)
            data[x] = list(np.float_(tmpdata))
            #[i for i, x in enumerate(TTdata) if x == "whatever"]
            if x =='Temperature [K]':
                break
    plt.xscale("log")
    legend = list()
    print(header)
    for x in header:
        if x == 'Temperature [K]':
            break
        legend.append(x)
        xpoints = data[x][0]
        ypoints = data[x][1]
        #plt.plot(xpoints, ypoints)
    #plt.legend(legend)
    #plt.show()
    perlite, bainite = dict(), dict()
    perlite["T"], perlite["n"], perlite["tau"] = JMAKfit(data['Start time (2% pearlite)'], data['Half time (50% pearlite)'], data['Finish time (98% pearlite)'])
    bainite["T"], bainite["n"], bainite["tau"] = JMAKfit(data['Start time (2% bainite)'], data['Half time (50% bainite)'], data['Finish time (98% bainite)'])
    martensite = KMfit(400,300,200)
    # Tfer, nfer, taufer = JMAKfit(data['Ferrite start'], data['Ferrite half'], data['Ferrite finish'])
    # print(Tfer)

    #xpoints = -tauper * (-np.log(0.98)) ** (1 / nper)
    #ypoints = Tper
    #plt.plot(xpoints, ypoints)
    #xpoints = -taubai * (-np.log(0.98)) ** (1 / nbai)
    #ypoints = Tbai
    #plt.plot(xpoints, ypoints)
    #plt.xscale("log")
    #plt.show()

    saveresult("JMAK_perlite", [perlite["T"], perlite["n"], perlite["tau"]])
    saveresult("JMAK_bainite", [bainite["T"], bainite["n"], bainite["tau"]])
def JMAKfit(data1,data2,data3):
    import numpy as np
    tau = np.array([])
    n = np.array([])
    Tlist = np.array([])
    for x in data1[1]:
        i = np.where(data3[1] == x)[0]
        j = np.where(data1[1] == x)[0]
        if i.size==0 or j.size==0:
            pass
        elif data1[0][j[0]] == data3[0][i[0]]:
            pass
        else:
            i = i[0]
            j = j[0]
            tmpn = np.log(np.log(0.98)/np.log(0.02))/np.log(data1[0][j]/data3[0][i])
            tmptau = - data1[0][j]/(-np.log(0.98))**(1/tmpn)
            n = np.append(n, tmpn)
            tau = np.append(tau, tmptau)
            Tlist = np.append(Tlist, x)
    return Tlist, n, tau

def KMfit(data1,data2,data3):
    import numpy as np
    def objective(x, Ms, beta):
        return 1 - np.exp(-beta * (Ms - x))
    # 0.02 = 1- exp(-beta * (Ms - data1))
    pass