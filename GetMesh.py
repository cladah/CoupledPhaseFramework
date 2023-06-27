import pandas as pd
import numpy as np
import h5py
from HelpFile import *

def createmesh(modelvar):
    if modelvar["elementtype"] in ["Spherical", "AxisymPstrain", "AxisymPstress"]:
        mesh1D(modelvar)
    elif modelvar["elementtype"] == "Axisym2D":
        mesh2D(modelvar)
    else:
        raise KeyError("Can't create mesh due to wrong element type")

def mesh1D(modelvar):
    columns = ['coordinates', 'dispx', 'loadx']
    data = list()
    data.append([0, 0, np.nan])                         # Adding 0 displacement for node at center
    for x in range(1, modelvar["nodesFEM"]):
        data.append([modelvar["radius"] * (x / (modelvar["nodesFEM"] - 1)), np.nan, 0])
    nodes = pd.DataFrame(data=data, columns=columns)

    columns = ['start', 'end', 'area', 'material']
    data = list()
    for x in range(modelvar["nodesFEM"] - 1):
        data.append([x, x + 1, 5, 1])
    elements = pd.DataFrame(data=data, columns=columns)

    with h5py.File("Cache_Mesh.hdf5", "w") as f:
        f.create_dataset("Nodes", data=nodes)
        f.create_dataset("Elements", data=elements)

def mesh2D(modelvar):
    columns = ["cordx", "cordy", "dispx", "dispy", "loadx", "loady"]
    data = list()
    data.append([0, 0, 0, 0, np.nan, np.nan])  # Adding 0 displacement for node at center
    for x in range(1, modelvar["nodesFEM"]):
        for y in range(1, modelvar["nodesFEM"]):
            if y < modelvar["nodesFEM"]/2 and x < modelvar["nodesFEM"]/2:
                data.append([modelvar["radius"] * (x / (modelvar["nodesFEM"] - 1)), modelvar["radius"] * (y / (modelvar["nodesFEM"] - 1)), np.nan, np.nan, 0, 0])
            elif y > modelvar["nodesFEM"]/2:
                data.append([modelvar["radius"] * (x / (modelvar["nodesFEM"] - 1)), modelvar["radius"]/2 + modelvar["radius"]/2 * (y / (modelvar["nodesFEM"] - 1)), np.nan, np.nan, 0, 0])
            elif x > modelvar["nodesFEM"]/2:
                pass
    nodes = pd.DataFrame(data=data, columns=columns)

    columns = ['n1', 'n2', 'n3', 'n4']
    data = list()
    for x in range(modelvar["nodesFEM"] - 1):
        data.append([x, x + 1, 5, 1])
    elements = pd.DataFrame(data=data, columns=columns)
    with h5py.File("Cache_Mesh.hdf5", "w") as f:
        f.create_dataset("Nodes", data=nodes)
        f.create_dataset("Elements", data=elements)
    pass

def readMesh(modelvar):
    #if not testchacheinfo(modelvar):
        #print('Change in input, remaking mesh')
    createmesh(modelvar)

    with h5py.File("Cache_Mesh.hdf5", "r") as f:
        nodes = pd.DataFrame(f.get('Nodes'))
        elements = pd.DataFrame(f.get('Elements'))
    nodes.columns = ['coordinates', 'displacement', 'load']
    elements.columns = ['start', 'end', 'area', 'material']

    return nodes, elements

