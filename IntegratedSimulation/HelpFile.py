import numpy as np
import h5py
from datetime import datetime

def createcachfile(modelvar):
    info = np.array([modelvar["radius"], modelvar["nodesFEM"], modelvar["E"], modelvar["nu"], modelvar["shapef"],
                     modelvar["Meshscaling"]])
    with h5py.File("Cache.hdf5", "w") as f:
        CNcurves = f.create_group("CNcurves")
        CNcurves.create_dataset("C", data=0)
        CNcurves.create_dataset("N", data=0)
        Cache = f.create_group("Cache")
        Cache.create_dataset("K", data=0)
        Cache.create_dataset("B", data=0)
        Cache.create_dataset("Mesh", data=0)
        matrixes = f.create_group("Matrixes")
        matrixes.create_dataset("K_FEM", data=0)
        matrixes.create_dataset("Mesh", data=0)
        matrixes.create_dataset("Bmatrix", data=0)
        f.create_dataset("Info", data=info)


def savetocache(dataname ,data):
    with h5py.File("Cache.hdf5", "r+") as f:
        del f[dataname]
        f.create_dataset(dataname, data=data)


def retrievecache(dataname):
    with h5py.File("Cache.hdf5", "r") as f:
        data = np.array(f.get(dataname))
    return data

def comparecache(dataname,data):
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
    cachename = "Result.hdf5"
    tmpinfo = [modelvar["E"]]
    with h5py.File(cachename, "w") as f:
        info = f.create_dataset("Info",data=tmpinfo)
        CNcurves = f.create_dataset("CNcurves",data=0)
        displacements = f.create_dataset("displacement",data= 0)

def saveresult(dataname,data):
    with h5py.File("Result.hdf5", "r+") as f:
        del f[dataname]
        f.create_dataset(dataname, data=data)

def readresultfile(dataname):
    try:
        with h5py.File("Result.hdf5", "r") as f:
            data = np.array(f.get(dataname))
        return data
    except:
        raise KeyError("Result "+str(dataname)+" doesn't exist in result file")


def renameResultfile():
    now = datetime.now().strftime("%Y%m%d%H%M%S")
    cachename = "Result_" + now + ".hdf5"
    with h5py.File(cachename, "w") as f:
        CNcurves = f.create_dataset("CNcurves",data=[1,2,3])

def h5_tree(val, pre=''):
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