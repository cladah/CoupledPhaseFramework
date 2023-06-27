import pandas as pd
import numpy as np
import sympy as sym
from sympy import *
import matplotlib.pyplot as plt
import scipy.linalg as la
from GetMesh import *
from GetMatrixes import *
from HelpFile import *

def calculatelengths(modelvar):
    Mesh = readMesh(modelvar)
    nodes = Mesh[0]
    elements = Mesh[1]
    lengths = list()
    for e in range(len(elements)):
        lengths.append(
            abs(nodes.loc[elements.loc[e]['start'], 'coordinates'] - nodes.loc[elements.loc[e]['end'], 'coordinates']))
    return lengths

def sphericalsolver(RHS, modelvar):
    if not testchacheinfo(modelvar):
        print("Resetting all matrixes")
        resetcalculations()
    else:
        pass

    Kfull = readK(modelvar)
    #print(Kfull)
    Mesh = readMesh(modelvar)
    nodes = Mesh[0]
    elements = Mesh[1]

    # thstrains = readstrains()
    # tfstrains = readstrains()
    nodes.loc[modelvar["nodesFEM"]-1, 'load'] = 1E9 # Adding a force at the last node

    A = nodes['displacement'].isna()            # Known loads
    B = nodes['load'].isna()                        # Unknown displacements

    KAA = Kfull[np.ix_(A, A)]
    KAB = Kfull[np.ix_(A, B)]
    KBA = Kfull[np.ix_(B, A)]
    KBB = Kfull[np.ix_(B, B)]

    U = nodes['displacement']
    P = nodes['load']
    Uk = U[B]
    Pk = P[A]
    # The known displacements are Uk The unknown are Uuk
    # The known loads are Pk The unknown are Puk

    Uuk = np.dot(np.linalg.inv(KAA), (Pk - np.dot(KAB, Uk)))
    Puk = np.dot(KBA,Uuk)+np.dot(KBB,Uk)
    results = nodes.copy()
    results.loc[A,'displacement'] = Uuk
    results.loc[B,'load'] = Puk


    displacement = results['displacement']
    err = 0
    return displacement, err

