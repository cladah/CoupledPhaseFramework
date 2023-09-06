from dolfinx import fem, io
from mpi4py import MPI
import ufl

from ..HelpFile import read_input
#import pyvista as pv
import ufl
import h5py
import os
def post_stress(uh, sigma):
    directory = os.getcwd() + '/Resultfiles'
    directory = directory.replace('\\', '/')
    indata = read_input()

    # Material
    E = indata["material"]["Austenite"]["E"]
    nu = indata["material"]["Austenite"]["nu"]
    rho_g = indata["material"]["rho"] * 9.82
    Cp_rho = indata["material"]["Cp"] * indata["material"]["rho"]
    k = indata["material"]["k"]
    mu = E / 2 / (1 + nu)
    lmbda = E * nu / (1 + nu) / (1 - 2 * nu)
    alpha = indata["material"]["alpha_k"]

    import ufl
    def eps(v):
        return ufl.sym(ufl.grad(v))

    def sig(v, T, T0):
        return 2.0 * mu * eps(v) + lmbda * ufl.tr(eps(v)) * ufl.Identity(len(v)) - (3 * lmbda + 2 * mu) * alpha * (
                T - 800.) * ufl.Identity(len(v))

    s = sigma - 1. / 3 * ufl.tr(sigma) * ufl.Identity(len(uh))
    von_Mises = ufl.sqrt(3. / 2 * ufl.inner(s, s))
    return (von_Mises)