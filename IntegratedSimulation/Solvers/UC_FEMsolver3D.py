from mpi4py import MPI
from dolfinx import mesh, fem, io
from dolfinx.fem import FunctionSpace
from dolfinx.io import gmshio
import ufl
import numpy as np
#import gmsh
from petsc4py.PETSc import ScalarType

def runFEM():
    # --------------- Loading mesh ------------------#

    msh, cell_markers, facet_markers = gmshio.read_from_msh("Resultfiles/Mesh.msh", MPI.COMM_WORLD, 0, gdim=3)
    domain = msh

    # --------------- Input assignment ------------------#

    # Material parameters #--------------------------------- Read from file saved under /Cachefiles
    E = fem.Constant(domain, 1e5)
    nu = fem.Constant(domain, (0.3))
    rho_g = 1e-3
    mu = E / 2 / (1 + nu)
    lmbda = E * nu / (1 + nu) / (1 - 2 * nu)

    # --------------- Creating functionspaces ------------------#

    V = fem.VectorFunctionSpace(domain, ("Lagrange", 1))

    # --------------- Finding boundary conditions ------------------#

    def xsym_boundary(x):
        return np.isclose(x[0], 0)
    def ysym_boundary(x):
        return np.isclose(x[1], 0)

    def outer(x):
        return np.isclose(np.sqrt(x[0]**2+x[1]**2), 1)

    fdim = domain.topology.dim - 1
    xsym_facets = mesh.locate_entities_boundary(domain, fdim, xsym_boundary)
    ysym_facets = mesh.locate_entities_boundary(domain, fdim, ysym_boundary)

    u_D = np.array([0, 0, 0], dtype=ScalarType)
    xsym_bc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, xsym_facets), V)
    ysym_bc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, ysym_facets), V)
    bc = [xsym_bc, ysym_bc]

    # --------------- Tractions on all remaining edges ------------------#

    T = fem.Constant(domain, ScalarType((0, 100, 0)))
    ds = ufl.Measure("ds", domain=domain)

    # --------------- Strain and stress definitions ------------------#
    def epsilon(u):
        return ufl.sym(ufl.grad(u))  # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)

    def sigma(u):
        return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
    L = ufl.dot(T, v) * ds

    # --------------- Problem formulation ------------------#

    problem = fem.petsc.LinearProblem(a, L, bcs=bc, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    uh = problem.solve()
    uh.name = "displacement"


    with io.XDMFFile(domain.comm, "Resultfiles/displacement.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(uh)
        #xdmf.write_function(stresses)