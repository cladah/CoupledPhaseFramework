from mpi4py import MPI
from dolfinx import mesh, fem, io
from dolfinx.fem import FunctionSpace
from dolfinx.io import gmshio
import ufl
import numpy as np
from petsc4py.PETSc import ScalarType
import gmsh

from CreateCircleMesh import *
r = 1

gmsh.initialize()
model = gmsh.model()
model.add("circle")

factory = model.geo

lcent = r/100
lcirc = r/100

#membrane = model.occ.addDisk(0, 0, 0, 1, 1)

factory.addPoint(0,0,0,lcent,1)
factory.addPoint(r,0,0,lcirc,2)
factory.addPoint(0,r,0,lcirc,3)

factory.addLine(3,1,1)
factory.addLine(1,2,2)
factory.addCircleArc(2,1,3,3)

factory.addCurveLoop([1,2,3],4)
factory.addPlaneSurface([4],5)
factory.synchronize()
#model.addPhysicalGroup(2, [membrane], 1)
#gmsh.option.setNumber("Mesh.CharacteristicLengthMin",0.05)
#gmsh.option.setNumber("Mesh.CharacteristicLengthMax",0.05)
model.mesh.generate(2)

#gmsh.write("circle.msh")


#domain = mesh.create_interval(MPI.COMM_WORLD, nx=50, points=(0.0, 1.0))

domain, cell_markers, facet_markers = gmshio.model_to_mesh(model, MPI.COMM_WORLD, 0)
V = fem.FunctionSpace(domain, ("Lagrange", 1))
#P1 = ufl.element("Lagrange", domain.basix_cell(), 1)
#P2v = ufl.element("Lagrange", domain.basix_cell(), 3)

#P2v = VectorFunctionSpace(mesh, "Lagrange", 2)  # Solid mechanics space
#P1  = FunctionSpace(mesh, "CG", 1)              # Heat and
#ME  = MixedFunctionSpace([P2v, P1, P1])
with io.XDMFFile(domain.comm, "output.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)




#(w, q, r, z) = TestFunctions(ME)
#dU           = TrialFunction(ME)