def heatequation(T,r):
    y = 1/r**2*dr*(r**2*dT/dr)
    return y

def heatsolver(modelvar,mesh,ti,dt,Ti,bcth):
    from BackwardEulerClas import *
    yi = eulerbackwardstep(heatequation(), modelvar, mesh, ti, dt, Ti, bcth)
    return yi

