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

    msh, cell_markers, facet_markers = gmshio.read_from_msh("Resultfiles/Mesh.msh", MPI.COMM_WORLD, 0, gdim=2)
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

    def ysym_boundary(x):
        return np.isclose(x[0], 0)

    def xsym_boundary(x):
        return np.isclose(x[1], 0)

    fdim = domain.topology.dim - 1
    ysym = mesh.locate_entities_boundary(domain, fdim, ysym_boundary)
    xsym = mesh.locate_entities_boundary(domain, fdim, xsym_boundary)

    u_D = np.array([0, 0], dtype=ScalarType)
    xbc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, xsym), V)
    ybc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, ysym), V)
    bc = [xbc, ybc]

    # --------------- Tractions on all remaining edges ------------------#

    T = fem.Constant(domain, ScalarType((0, 0)))
    ds = ufl.Measure("ds", domain=domain) # No idea
    p = fem.Constant(domain, ScalarType(-10))

    f = fem.Function(V)
    dofs = fem.locate_dofs_geometrical(V, lambda x: np.isclose(x.T, 1.0))
    f.x.array[dofs] = 100

    x = ufl.SpatialCoordinate(domain)

    # --------------- Strain and stress definitions ------------------#
    def eps(u):
        return ufl.sym(ufl.grad(u))

    def axisym_eps(v):
        return ufl.sym(ufl.as_tensor([[v[0].dx(0), 0, v[0].dx(1)],
                              [0, v[0]/x[0], 0],
                              [v[1].dx(0), 0, v[1].dx(1)]]))

    # Stress given a displacement field
    def sigma(u):
        return lmbda*ufl.tr(axisym_eps(u))*ufl.Identity(2) + 2.0*mu*axisym_eps(u)

    # --------------- Variational formulation ------------------#

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    # Defining the problem
    #a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
    #L = p * v * ufl.dx
    a = (ufl.inner(sigma(u), axisym_eps(v))*x[0]**2)*ufl.dx # axisymmetric definitions
    L = (ufl.inner(f, v)*x[0]**2)*ufl.dx # axisymmetric definitions

    # --------------- Problem formulation ------------------#

    problem = fem.petsc.LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    uh = problem.solve()
    uh.name = "displacement"

    s = sigma(uh) -1/3*ufl.tr(sigma(uh))*ufl.Identity(2)
    von_Mises = ufl.sqrt(3./2*ufl.inner(s, s))
    V_vm = fem.FunctionSpace(domain, ("DG", 1))
    stress_expr = fem.Expression(von_Mises, V_vm.element.interpolation_points())
    vm_stresses = fem.Function(V_vm, name="vm_stress")
    vm_stresses.interpolate(stress_expr)
    V_stress = fem.VectorFunctionSpace(domain, ("DG", 1))
    stress_expr = fem.Expression(s, V_stress.element.interpolation_points())
    stresses = fem.Function(V_stress, name="stress")
    #stresses.interpolate(stress_expr.sub(0))
    #.sub(0)

    with io.XDMFFile(domain.comm, "Resultfiles/displacement.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(uh)
        xdmf.write_function(vm_stresses)
        #xdmf.write_function(stresses)