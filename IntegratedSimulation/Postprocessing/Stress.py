from dolfinx import fem, io
from mpi4py import MPI
import ufl
def post_stress(uh, sigma):
    s = sigma - 1. / 3 * ufl.tr(sigma) * ufl.Identity(len(uh))
    von_Mises = ufl.sqrt(3. / 2 * ufl.inner(s, s))
    return (von_Mises)