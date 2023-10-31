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



def createinputcache():
    import json
    f = open('Cachefiles/InputCache.json', 'w')
    data = read_input()
    json.dump(data, f, indent=2)
    f.close()

def adjustinputcache(model):
    import json
    f = open('Cachefiles/Input.json', 'r')
    indata = json.load(f)
    f.close()
    f = open('Cachefiles/InputCache.json', 'w')
    data = read_input()
    if model == 'Mesh':
        data['Geometry'] = indata['Geometry']
    elif model == 'Quenching':
        for x in ['Geometry', 'Material', 'Thermo', 'FEM', 'Programs']:
            data[x] = indata[x]
    else:
        raise KeyError('Adjust input cache not implemented')
    json.dump(data, f, indent=2)
    f.close()

def checkinput(model):
    import json
    import numpy as np
    f = open('Cachefiles/Input.json', 'r')
    indata = json.load(f)
    f.close()
    f = open('Cachefiles/InputCache.json', 'r')
    cachedata = json.load(f)
    f.close()
    # Check rerun criteria
    if indata['Rerun'][model] == 1:
        return False
    if model == 'Mesh':
        for x in ['Geometry']:
            if indata[x] != cachedata[x]:
                return False
    elif model == 'Quenching':
        for x in ['Geometry', 'Material', 'Thermo', 'FEM', 'Programs']:
            if indata[x] != cachedata[x]:
                return False
    return True



def readCNfile():
    import queue
    with open('Cachefiles/50min_CarbonNitration.txt') as f:
        lines = f.readlines()
        header = lines.pop(0).split()
        material = dict()
        for i in range(3, len(header), 2):
            material[header[i].strip("W()")] = []
        for line in lines:
            line = line.split()
            for x in material.keys():
                material[x].append(float(line.pop(0)))
    return material




def savetocache(dataname ,data):
    import h5py
    try:
        with h5py.File("ThermoCache.hdf5", "r+") as f:
            pass
    except:
        with h5py.File("ThermoCache.hdf5", "w") as f:
            pass

    with h5py.File("ThermoCache.hdf5", "r+") as f:
        del f[dataname]
        f.create_dataset(dataname, data=data)


def retrievecache(dataname):
    import h5py
    with h5py.File("ThermoCache.hdf5", "r") as f:
        data = np.array(f.get(dataname))
    return data

def comparecache(dataname,data):
    import h5py
    with h5py.File("ThermoCache.hdf5", "r") as f:
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


def createresultfile():
    data = read_input()
    # Composition / point
    geo = [1 * 0.9**i for i in range(data['Geometry']['nodes'])]
    x = [1 * data['Geometry']['meshscaling']**i for i in range(data['Geometry']['nodes'])]*np.linspace(0,1,data['Geometry']['nodes'])
    x = data['Geometry']['radius']*x/x[-1]
    # Mesh
    import h5py
    cachename = "Result.hdf5"
    with h5py.File(cachename, "w") as f:
        CNcurves = f.create_dataset("CNcurves",data=0)
        displacements = f.create_dataset("displacement",data= 0)

def saveresult(dataname,data):
    import h5py
    with h5py.File("Result.hdf5", "r+") as f:
        try:
            del f[dataname]
        except:
            pass
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