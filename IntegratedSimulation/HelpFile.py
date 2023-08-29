import numpy as np

class Material:
    def __init__(self, name, E, nu, alpha, f):
        self.name = name
        self.E = E              # Youngs modulus [Pa]
        self.nu = nu            # Poissons ratio [-]
        self.alpha = alpha      # Thermal expansion coefficient
        self.f = f              # Material fraction
def read_input():
    import json
    f = open('Cachefiles/Input.json', 'r')
    data = json.load(f)
    f.close()
    return data



def createinputcach():
    import json
    f = open('Cachefiles/InputCache.json', 'w')
    data = read_input()
    json.dump(data)
    f.close()

def checkinput():
    import json
    import numpy as np
    f = open('Cachefiles/Input.json', 'r')
    indata = json.load(f)
    f.close()
    f = open('Cachefiles/InputCache.json', 'r')
    cachedata = json.load(f)
    f.close()
    x = np.zeros(5)
    if indata["Geometry"] == cachedata["Geometry"]:
        x[0] = 1
    if indata["Thermo"] == cachedata["Thermo"]:
        x[1] = 1
    return

def savetocache(dataname ,data):
    import h5py
    with h5py.File("Cache.hdf5", "r+") as f:
        del f[dataname]
        f.create_dataset(dataname, data=data)


def retrievecache(dataname):
    import h5py
    with h5py.File("Cache.hdf5", "r") as f:
        data = np.array(f.get(dataname))
    return data

def comparecache(dataname,data):
    import h5py
    with h5py.File("Cache.hdf5", "r") as f:
        testdata = np.array(f.get(dataname))
    if testdata == data:
        return True
    else:
        return False

def testchacheinfo(modelvar):
    fileinfo = retrievecache("Info")
    info = np.array([modelvar["radius"],modelvar["nodesFEM"],modelvar["E"],modelvar["nu"],modelvar["shapef"],modelvar["Meshscaling"]])
    if np.equal(fileinfo,info).all():
        return True
    else:
        return False

def checkcalculation(calcname):
    try:
        #precalc = ["K","B","Mesh"]
        #calc = 0
        #for x in precalc:
        #    dataname = "Cache/"+x
        #    calc += retrievecache(dataname)
        if retrievecache(calcname)==0:
            return False
        else:
            return True
    except:
        print("Cache file unfunctioning, resetting")
        resetcalculations()


def resetcalculations():
    print("Resetting the cache file")
    createcachfile()


def createresultfile(modelvar):
    import h5py
    cachename = "Result.hdf5"
    tmpinfo = [modelvar["E"]]
    with h5py.File(cachename, "w") as f:
        info = f.create_dataset("Info",data=tmpinfo)
        CNcurves = f.create_dataset("CNcurves",data=0)
        displacements = f.create_dataset("displacement",data= 0)

def saveresult(dataname,data):
    import h5py
    with h5py.File("Result.hdf5", "r+") as f:
        del f[dataname]
        f.create_dataset(dataname, data=data)

def readresultfile(dataname):
    import h5py
    try:
        with h5py.File("Result.hdf5", "r") as f:
            data = np.array(f.get(dataname))
        return data
    except:
        raise KeyError("Result "+str(dataname)+" doesn't exist in result file")


def renameResultfile():
    import h5py
    from datetime import datetime
    now = datetime.now().strftime("%Y%m%d%H%M%S")
    cachename = "Result_" + now + ".hdf5"
    with h5py.File(cachename, "w") as f:
        CNcurves = f.create_dataset("CNcurves",data=[1,2,3])

def h5_tree(val, pre=''):
    import h5py
    items = len(val)
    for key, val in val.items():
        items -= 1
        if items == 0:
            # the last item
            if type(val) == h5py._hl.group.Group:
                print(pre + '└── ' + key)
                h5_tree(val, pre+'    ')
            else:
                print(pre + '└── ' + key + ' (%d)' % len(val))
        else:
            if type(val) == h5py._hl.group.Group:
                print(pre + '├── ' + key)
                h5_tree(val, pre+'│   ')
            else:
                print(pre + '├── ' + key + ' (%d)' % len(val))