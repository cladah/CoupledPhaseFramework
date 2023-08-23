from mpi4py import MPI
from dolfinx import mesh, fem, io
from dolfinx.fem import FunctionSpace
from dolfinx.io import gmshio
import ufl
import numpy as np
#import gmsh
from petsc4py.PETSc import ScalarType


from MeshConverter import convert_mesh
from UC_FEMsolver import runFEM
from UC_Heatsolver import runHeat
from demo_stokes import runDEMO

from ThermoMech import runThermoMech
from UC_Phasesolver import runPhase
from Linearelast import runelast

convert_mesh()
#runFEM()
runHeat()
#runDEMO()
#runThermoMech()
#runPhase()