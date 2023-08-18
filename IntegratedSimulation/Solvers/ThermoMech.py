from mpi4py import MPI
from dolfinx import mesh, fem, io
from dolfinx.fem import FunctionSpace, Constant
from dolfinx.io import gmshio
import ufl
import numpy as np
#import gmsh
from petsc4py.PETSc import ScalarType

def runThermoMech():
    # --------------- Loading mesh ------------------#

    msh, cell_markers, facet_markers = gmshio.read_from_msh("Resultfiles/Mesh.msh", MPI.COMM_WORLD, 0, gdim=2)
    domain = msh

    # --------------- Loading inputs ------------------#

    T0 = Constant(domain,293.)
    E = Constant(domain,70e3)
    nu = 0.3
    lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / 2 / (1 + nu)
    rho = Constant(domain,2700.)
    alpha = Constant(domain,2.31e-5)  # thermal expansion coefficient
    kappa = alpha * (2 * mu + 3 * lmbda)
    cV = 910e-6 * rho  # specific heat per unit volume at constant strain
    k = Constant(domain, 237e-6)  # thermal conductivity

    # --------------- Defining elements ------------------#

    Vue = ufl.VectorElement('CG', domain.ufl_cell(), 2)  # displacement finite element
    Vte = ufl.FiniteElement('CG', domain.ufl_cell(), 1)  # temperature finite element


    mixed_el = ufl.MixedElement([Vue, Vte])
    V = FunctionSpace(domain, mixed_el)
    num_subs = V.num_sub_spaces
    spaces = []
    maps = []
    for i in range(num_subs):
        space_i, map_i = V.sub(i).collapse()
        spaces.append(space_i)
        maps.append(map_i)
    # --------------- Defining boundary conditions ------------------#
    def y_axis(x):
        return np.isclose(x[0], 0)

    def x_axis(x):
        return np.isclose(x[1], 0)

    fdim = domain.topology.dim - 1
    u_D = np.array([0.0, 0.0], dtype=ScalarType)
    bc1 = fem.dirichletbc(u_D, fem.locate_dofs_topological(V.sub(0), fdim, x_axis), V.sub(0).sub(1))
    bc2 = fem.dirichletbc(u_D, fem.locate_dofs_topological(V.sub(0), fdim, y_axis), V.sub(0).sub(0))
    bcs = [bc1, bc2]

    # --------------- Variational formulation ------------------#

    U = ufl.TestFunction(V)
    (u, temp) = ufl.split(U)
    dU = ufl.TrialFunction(V)
    (du, dtemp) = ufl.split(dU)
    Uold = ufl.Function(V)
    (uold, Thetaold) = ufl.split(Uold)

    def eps(v):
        return ufl.sym(ufl.grad(v))

    def sigma(v, Theta):
        return (lmbda * ufl.tr(eps(v)) - kappa * Theta) * ufl.Identity(2) + 2 * mu * eps(v)

    dt = Constant(domain,0.)
    mech_form = ufl.inner(sigma(du, dtemp), eps(du)) * ufl.dx
    therm_form = (cV * (dtemp - Thetaold) / dt * temp +
                  kappa * T0 * ufl.tr(eps(du - uold)) / dt * temp +
                  ufl.dot(k * ufl.grad(dtemp), ufl.grad(temp))) * ufl.dx
    form = mech_form + therm_form

    # --------------- Solving ------------------#

    Nincr = 100
    t = np.logspace(1, 4, Nincr + 1)
    Nx = 100
    x = np.linspace(R, L, Nx)
    T_res = np.zeros((Nx, Nincr + 1))
    U = fem.Function(V)
    for (i, dti) in enumerate(np.diff(t)):
        print("Increment " + str(i + 1))
        dt.assign(dti)
        ufl.solve(ufl.lhs(form) == ufl.rhs(form), U, bcs)
        Uold.assign(U)
        T_res[:, i + 1] = [U(xi, 0.)[2] for xi in x]
    u, Theta = ufl.split(U)