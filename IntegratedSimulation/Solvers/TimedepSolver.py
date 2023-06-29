from IntegratedSimulation.Solvers.Heatsolver import *
def coupled_solver(modelvar):
    heat = 1
    phase = 1
    FEM = 1

    mesh = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    dt = 1
    ti = 0
    tstop = 10
    y0 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    T0 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    bcth = [0]
    yi = y0
    Ti = T0
    while ti < tstop:
        if heat == 1:
            Ti = heatsolver(modelvar)
        if phase == 1:
            pass
        if FEM == 1:
            pass
        ti += dt



