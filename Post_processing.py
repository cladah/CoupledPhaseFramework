from GetMesh import *
from GetMatrixes import *
import numpy as np
import pandas as pd
import sympy as sym
from sympy import *
import h5py
from HelpFile import *
import matplotlib.pyplot as plt

def getstrains(modelvar,U):
    Bfull = readB(modelvar)
    eps = np.dot(Bfull,U)
    #print(U.shape)
    #print(eps.shape)
    return eps

def getstresses(modelvar, eps):
    nu = modelvar["nu"]
    E = modelvar["E"]
    sigma = zeros(len(eps),1)
    for i in range(int(len(eps)/2)):
        sigma[2*i] = E / ((1 + nu) * (1 - 2 * nu)) * np.array([1 - nu, nu]) * eps[2*i]
        sigma[2 * i] = E / ((1 + nu) * (1 - 2 * nu)) * np.array([nu,1 - nu]) * eps[2 * i + 1]
    return sigma

def plotstrain(modelvar):
    Mesh = readMesh(modelvar)
    elements = Mesh[1]

    displacement = readresultfile("displacement")
    strains = getstrains(modelvar, displacement)

    plt.figure(2, figsize=(10, 6))
    plt.rcParams.update({'font.size': 16})
    plt.plot(elements['start'], strains[range(1, 2 * len(elements.index), 2)])
    plt.plot(elements['start'], strains[range(0, 2 * len(elements.index), 2)])
    # plt.plot(elements['start'],np.dot(B,results['displacement']))
    plt.legend(['Circumferential [-]', 'Radial [-]'])
    plt.xlabel('Radius [mm]')
    plt.ylabel('Strain [-]')
    plt.show()

def plotdisplacement(modelvar):     # Plotting the strain at the gausspoints of the mesh
    Mesh = readMesh(modelvar)
    nodes = Mesh[0]

    displacement = readresultfile("displacement")

    plt.figure(1, figsize=(10, 6))
    plt.rcParams.update({'font.size': 16})
    plt.plot(nodes['coordinates'], displacement)
    # plt.plot(elements['start'],np.dot(B,results['displacement']))
    plt.legend(['Displacements'])
    plt.xlabel('Radius [m]')
    plt.ylabel('Displacement [m]')
    plt.show()

def plotstress(modelvar):       # Plotting the stress at the gausspoints of the mesh
    Mesh = readMesh(modelvar)
    elements = Mesh[1]
    displacement = readresultfile("displacement")
    strains = getstrains(modelvar,displacement)
    sigma = getstresses(modelvar,strains)

    plt.figure(3, figsize=(10, 6))
    plt.rcParams.update({'font.size': 16})
    plt.plot(elements['start'], sigma[range(1, 2 * len(elements.index), 2),0]/1E6)
    plt.plot(elements['start'], sigma[range(0, 2 * len(elements.index), 2),0]/1E6)
    # plt.plot(elements['start'],np.dot(B,results['displacement']))
    plt.legend(['Circumferential [-]', 'Radial [-]'])
    plt.xlabel('Radius [mm]')
    plt.ylabel('Stress [MPa]')
    plt.show()