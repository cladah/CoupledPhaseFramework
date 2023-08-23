import numpy as np
from mpi4py import MPI
from petsc4py.PETSc import ScalarType
from petsc4py import PETSc
from dolfinx import mesh, fem, io, nls
from dolfinx.io import gmshio
import ufl
import gmsh



def runFEM():

    # --------------- Loading mesh ------------------#
    domain, cell_markers, facet_markers = gmshio.read_from_msh("Resultfiles/Mesh.msh", MPI.COMM_WORLD, 0, gdim=2)

    # --------------- Loading input ------------------#
    E = 1.0e9
    nu = 0.3
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))

    # --------------- Formulating mixed element ------------------#
    V = fem.VectorFunctionSpace(domain, ("Lagrange", 1), dim=2)

    # --------------- Boundary conditions ------------------#
    def yaxis(x):
        return np.isclose(x[0], 0.)

    def xaxis(x):
        return np.isclose(x[1], 0.)

    fdim = 1
    ndim = 0
    yaxis_facets = mesh.locate_entities_boundary(domain, ndim, yaxis)
    xaxis_facets = mesh.locate_entities_boundary(domain, ndim, xaxis)

    # --------------- boundary condition value
    u_D = np.array(0.0, dtype=ScalarType)

    # --------------- Displacement field
    xbc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V.sub(1), ndim, xaxis_facets), V.sub(1))  # sub defines the component, y on xaxis
    ybc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V.sub(0), ndim, yaxis_facets), V.sub(0))  # sub defines the component x on yaxis
    # --------------- Total
    bcs = [xbc, ybc]

    # --------------- Loading ------------------#
    def outside(x):
        return np.isclose(np.sqrt(x[0]**2+x[1]**2), 1.)

    outside_f = mesh.locate_entities(domain, 1, outside)
    marked_values = np.hstack([np.full_like(outside_f, 1)])
    outside_ft = mesh.meshtags(domain, fdim, outside_f, marked_values)

    f = fem.Constant(domain, 1E9)
    ds = ufl.Measure("ds", domain=domain, subdomain_data=outside_ft)
    n = ufl.FacetNormal(domain)

    # --------------- Variational formulation ------------------#
    u = ufl.TrialFunction(V)
    du = ufl.TestFunction(V)

    def eps(v):
        return ufl.sym(ufl.grad(v))
    def sig(v):
        return 2.0 * mu * eps(v) + lmbda * ufl.nabla_div(v) * ufl.Identity(len(v))

    # --------------- Problem formulation ------------------#

    a = ufl.inner(sig(u), eps(du)) * ufl.dx
    L = ufl.dot(f*n, du) * ds

    # Setting up problem
    problem = fem.petsc.LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    disp = problem.solve()
    disp.name = "Displacement"


    with io.XDMFFile(domain.comm, "Resultfiles/displacement.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(disp)