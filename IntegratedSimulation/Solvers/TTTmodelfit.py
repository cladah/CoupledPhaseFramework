def TTTfit():
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
            if x =='Martensite start':
                break
    plt.xscale("log")
    legend = list()
    for x in header:
        if x == 'Temperature [K]':
            break
        legend.append(x)
        xpoints = data[x][0]
        ypoints = data[x][1]
        #plt.plot(xpoints, ypoints)
    #plt.legend(legend)
    #plt.show()
    #JMAKfit(data['Start time (2% pearlite)'], data['Half time (50% pearlite)'], data['Finish time (98% pearlite)'])
def JMAKfit(data1,data2,data3):
    import numpy as np
    tau = list()
    n = list()
    Tlist = list()
    tmpdata = list()
    for x in data1[1]:
       if np.where(data3[1]==x)[0].size!=0:

           i = np.where(data3[1] == x)[0][0]
           print(data3[0][i])
           print(data1[0][i])
           print(np.log(0.98-0.02))
           tmpn = np.log(np.log(0.98-0.02)-data1[1][i]/data3[1][i])
           print(tmpn)
           print(x)
           tmptau = - x/np.log(0.02)^(1/tmpn)
           n.append(tmpn)
           tau.append(tmptau)
           Tlist.append(x)
    print(Tlist)
    #0.2 = (1-np.exp(-(data1[0]/tau)^(n)))
    #0

    pass