from mpi4py import MPI
from dolfinx import mesh, fem, io
from dolfinx.fem import FunctionSpace
from dolfinx.io import gmshio
import ufl
import numpy as np
#import gmsh
from petsc4py.PETSc import ScalarType

from UC_FEMsolver2D import runFEM
from UC_Heatsolver import *
from UC_Phasesolver import runPhase


runFEM()
runHeat()
#runPhase()