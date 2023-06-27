from HelpFile import *
from GetMesh import *
from sympy import pi
import sympy as sym
import numpy as np

def getJacobian(r1,r2):
    J = 1 / (r2 - r1)  # dxi/dr
    return J

def gaussintegration(modelvar):
    xi = sym.symbols("xi")
    if modelvar["shapef"] == "Linear":
        #xi = 0.577350269189625764509148780502
        g = [np.sqrt(1/3), -np.sqrt(1/3)]
        w = [1, 1]
        N = [1-xi, xi]
        dN = [-1, 1]

    elif modelvar["shapef"]=="Quad":
        g = [0, np.sqrt(3/5), -np.sqrt(3/5)]
        w = [8/9, 5/9, 5/9]
        N = [2*xi**2-3*xi+1, -4*xi**2+4*xi, 2*xi**2-xi]
        dN = [4*xi-3, -8*xi+4, 4*xi-1]
    else:
        raise KeyError("Shapefunction unfunctioning in gaussintegration")
    return N, dN, g, w

def getBmatrix(modelvar,r1,r2):
    xi = sym.symbols("xi")
    nodes = readMesh(modelvar)[0]
    if modelvar["shapef"] == "Linear":
        Bfull = np.zeros((2 * modelvar["nodesFEM"] - 2, modelvar["nodesFEM"]))
    else:
        raise KeyError("Shapefunction unfunctioning in Bmatrix")
    N, dN = gaussintegration(modelvar)[0:2]

    r = xi * (r2 - r1) + r1
    J = getJacobian(r1,r2)
    B = sym.Matrix(np.zeros((len(N),2)))
    for j in range(0,len(N)-1):
        tmpN = N[j]
        tmpdN = dN[j]
        tmpB = np.array([tmpdN * J, tmpN / r])
        B[j, 0] = tmpB[0]
        B[j, 1] = tmpB[1]
    return B

def getKmatrix(modelvar):
    xi = sym.symbols("xi")

    Kfull = np.zeros((modelvar["nodesFEM"], modelvar["nodesFEM"]))
    nodes = readMesh(modelvar)[0]

    g, w = gaussintegration(modelvar)[2:4]
    nu = modelvar["nu"]
    E = modelvar["E"]
    C = E / ((1 + nu) * (1 - 2 * nu)) * np.array(([1 - nu, nu], [nu, 1 - nu]))

    if modelvar["elementtype"] == "Spherical":
        C = C[0:2,0:2]
        for i in range(modelvar["nodesFEM"] - 1):
            r1 = nodes.loc[i, 'coordinates']
            r2 = nodes.loc[i + 1, 'coordinates']
            B = getBmatrix(modelvar,r1,r2)
            B.subs(xi,)
            r = xi * (r2 - r1) + r1
            k = np.zeros((len(g), modelvar["nodesFEM"]))
            for j in range(0,len(g)-1): # Gauss integration
                tmpB = B.subs(xi,g[j])
                tmpr = r.sbs(xi,g[j])
                ktmp = 4 * pi * np.dot(np.dot(tmpB.transpose(), C), tmpB) * tmpr ** 2
                k += w[j]*ktmp

            # Adding K matrix for element to full matrix
            Kfull[i, i] = Kfull[i, i] + k[0, 0]
            Kfull[i, i + 1] = k[0, 1]
            Kfull[i + 1, i] = Kfull[i, i + 1]
            Kfull[i + 1, i + 1] = k[1, 1]
        savetocache("Matrixes/K_FEM", Kfull)
    elif modelvar["elementtype"] == "AxisymPstrain":
        pass
    else:
        raise KeyError("Wrong element type in Kmatrix function")


def readN(modelvar):
    xi = 0.577350269189625764509148780502
    N1 = 1 - xi
    N2 = xi
    N = [N1,N2]
    return N

def createB(modelvar):
    Mesh = readMesh(modelvar)
    nodes = Mesh[0]

    if modelvar["elementtype"]=="Spherical"or modelvar["elementtype"]=="AxisymPstrain":
        Bfull = np.zeros((2 * modelvar["nodesFEM"] - 2, modelvar["nodesFEM"]))
        if modelvar["shapef"] == 1:
            xi = 0.577350269189625764509148780502
            N1 = 1-xi
            N2 = xi

            # Differentiating the shape functions
            dN1 = -1  # dN1/dxi
            dN2 = 1  # dN2/dxi
            for i in range(0, modelvar["nodesFEM"] - 1):
                r1 = nodes.loc[i, 'coordinates']
                r2 = nodes.loc[i + 1, 'coordinates']
                J = 1 / (r2 - r1)  # dxi/dr
                r = xi * (r2 - r1) + r1
                B = np.array([[dN1 * J, dN2 * J], [N1 / r, N2 / r]])
                Bfull[2 * i, i] = B[0, 0]
                Bfull[2 * i, i + 1] = B[0, 1]
                Bfull[2 * i + 1, i] = B[1, 0]
                Bfull[2 * i + 1, i + 1] = B[1, 1]
        else:
            raise KeyError
        savetocache("Matrixes/Bmatrix",Bfull)
    elif modelvar["elementtype"]=="AxisymPstress":
        # rr
        # pp
        # rp

        Bfull = np.zeros((2 * modelvar["nodesFEM"] - 2, modelvar["nodesFEM"]))
        if modelvar["shapef"] == 1:
            # Gausintegration of linear element
            xi = 0.577350269189625764509148780502

            # Shapefinction
            N1 = 1 - xi
            N2 = xi

            # Differentiating the shape functions
            dN1 = -1  # dN1/dxi
            dN2 = 1  # dN2/dxi

            for i in range(0, modelvar["nodesFEM"] - 1):
                r1 = nodes.loc[i, 'coordinates']
                r2 = nodes.loc[i + 1, 'coordinates']
                J = 1 / (r2 - r1)  # dxi/dr
                r = xi * (r2 - r1) + r1
                B = np.array([[dN1 * J, dN2 * J], [N1 / r, N2 / r]])
                Bfull[2 * i, i] = B[0, 0]
                Bfull[2 * i, i + 1] = B[0, 1]
                Bfull[2 * i + 1, i] = B[1, 0]
                Bfull[2 * i + 1, i + 1] = B[1, 1]
        else:
            raise KeyError("Shapefunction is not implemented")
        savetocache("Matrixes/Bmatrix", Bfull)
    else:
        raise KeyError("Wrong element type")

def readB(modelvar):
    if not checkcalculation("Cache/B"):
        print('Change in input, recalculating strain matrix')
        createB(modelvar)
        savetocache("Cache/B", 1)
    Bfull = retrievecache("Matrixes/Bmatrix")
    return Bfull

def createK(modelvar):
    Kfull = np.zeros((modelvar["nodesFEM"], modelvar["nodesFEM"]))

    Mesh = readMesh(modelvar)
    nodes = Mesh[0]
    elements = Mesh[1]

    nu = modelvar["nu"]
    E = modelvar["E"]
    if modelvar["elementtype"]=="Spherical":
        C = E / ((1 + nu) * (1 - 2 * nu)) * np.array(([1 - nu, nu], [nu, 1 - nu]))

        xi = 0.577350269189625764509148780502

        #
        N1 = 1 - xi
        N2 = xi

        # Differentiating the shape functions
        dN1 = -1  # dN1/dxi
        dN2 = 1  # dN2/dxi

        for i in range(modelvar["nodesFEM"] - 1):
            r1 = nodes.loc[i, 'coordinates']
            r2 = nodes.loc[i + 1, 'coordinates']
            J = 1 / (r2 - r1)
            r = xi*(r2-r1)+r1
            B = np.array([[dN1 * J, dN2 * J], [N1 / r, N2 / r]])
            k = 4 * pi * np.dot(np.dot(B.transpose(),C),B)*r**2

            # Adding K matrix for element to full matrix
            Kfull[i, i] = Kfull[i, i] + k[0, 0]
            Kfull[i, i + 1] = k[0, 1]
            Kfull[i+1, i] = Kfull[i, i + 1]
            Kfull[i+1, i + 1] = k[1, 1]
        savetocache("Matrixes/K_FEM", Kfull)
    elif modelvar["elementtype"]=="AxisymPstrain":
        pass
    else:
        raise KeyError("Wrong element type")

def readK(modelvar):
    if not checkcalculation("Cache/K"):
        print('Change in input, recalculating stiffness matrix')
        createK(modelvar)
        savetocache("Cache/K", 1)
    else:
        print("Getting stiffness matrix")
    KFEM = retrievecache("Matrixes/K_FEM")

    return KFEM