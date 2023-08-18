from mpi4py import MPI
from dolfinx import mesh, fem, io, nls
import ufl
import numpy as np
from petsc4py import PETSc

def runHeat():
    # --------------- Loading mesh ------------------#

    domain, cell_markers, facet_markers = io.gmshio.read_from_msh("Resultfiles/Mesh.msh", MPI.COMM_WORLD, 0, gdim=2)
    domain.name = 'Sphere'

    # --------------- Loading inputs ------------------#
    T0 = fem.Constant(domain,293.)
    DThole = fem.Constant(domain,10.)
    E = fem.Constant(domain,70e3)
    nu = fem.Constant(domain,0.3)
    lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / 2 / (1 + nu)
    rho = fem.Constant(domain,2700.)  # density
    alpha = fem.Constant(domain,2.31e-5)  # thermal expansion coefficient
    kappa = alpha * (2 * mu + 3 * lmbda)
    cV = 910e-6 * rho  # specific heat per unit volume at constant strain
    k = fem.Constant(domain, 237e-6)  # thermal conductivity

    # --------------- Creating function space ------------------#

    V = fem.FunctionSpace(domain, ("CG", 1))

    # --------------- Assigning initial conditions ------------------#
    def initial_condition(x, a=500):
        return np.exp(-a * (x[0] ** 2 + x[1] ** 2))+200

    t0 = fem.Function(V)
    t0.name = "Temperature"
    t0.interpolate(initial_condition)

    # --------------- Finding boundary conditions ------------------#
    def ysym_boundary(x):
        return np.isclose(x[0], 0)

    def xsym_boundary(x):
        return np.isclose(x[1], 0)

    def circumf_boundary(x):
        return np.isclose(np.sqrt(x[0]**2+x[1]**2), 1) # Insert radius

    bccirc = fem.dirichletbc(PETSc.ScalarType(200), fem.locate_dofs_geometrical(V, circumf_boundary), V)
    bc = [bccirc]

    # --------------- Setting up transient inputs ------------------#

    t0 = 0  # Start time
    tf = 1.0  # Final time

    num_steps = 50

    # --------------- Setting formulation ------------------#

    Uold = fem.Function(V) #

    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)

    dt = fem.Constant(domain,0.)
    F = (cV * (u - Uold) / dt * v + ufl.dot(k * ufl.grad(u), ufl.grad(v))) * ufl.dx
    a = fem.form(ufl.lhs(F))
    L = fem.form(ufl.rhs(F))
    # a = rho*c*T*v*dx + theta*dt*kappa*inner(grad(T), grad(v))*dx
    # L = (rho*c*T_prev*v + dt*f*v - (1-theta)*dt*kappa*inner(grad(T), grad(v)))*dx
    t = np.logspace(t0, tf, num_steps)

    xdmf = io.XDMFFile(domain.comm, "Resultfiles/diffusion.xdmf", "w")
    xdmf.write_mesh(domain)
    for (dti) in t:
        dt = fem.Constant(domain, dti)
        problem = fem.petsc.LinearProblem(a, L, bcs=bc, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        uh = problem.solve()
        Uold.assign(uh)

        # Write solution to file
        xdmf.write_function(uh, t)
    xdmf.close()