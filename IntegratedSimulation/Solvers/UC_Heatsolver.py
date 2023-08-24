from mpi4py import MPI
from dolfinx import mesh, fem, io
from dolfinx.io import gmshio
from dolfinx.fem import FunctionSpace
from dolfinx.io import gmshio
import ufl
import numpy as np
#import gmsh
from petsc4py.PETSc import ScalarType
from petsc4py import PETSc


def runHeat():
    # --------------- Loading mesh ------------------#

    domain, cell_markers, facet_markers = gmshio.read_from_msh("Resultfiles/Mesh.msh", MPI.COMM_WORLD, 0, gdim=2)

    # --------------- Creating function space ------------------#
    V = fem.FunctionSpace(domain, ("CG", 1))

    # --------------- Loading inputs ------------------#
    tstart = 0  # Start time
    tstop = 100.  # Final time

    num_steps = 50
    dt = tstop / num_steps  # time step size

    # Material input
    rho = fem.Constant(domain, 2700.)
    cV = fem.Constant(domain, 910e-6) * rho
    k = fem.Constant(domain, 237e-6)

    # --------------- Assigning initial conditions ------------------#
    def initial_condition(x):
        return np.isreal(x[0])*800.

    T0 = fem.Function(V)
    T0.name = "Initial temperature"
    T0.interpolate(initial_condition)

    # --------------- Finding boundary conditions ------------------#
    def yaxis(x):
        return np.isclose(x[0], 0.)

    def xaxis(x):
        return np.isclose(x[1], 0.)

    fdim = 1
    ndim = 0
    yaxis_facets = mesh.locate_entities_boundary(domain, fdim, yaxis)
    xaxis_facets = mesh.locate_entities_boundary(domain, fdim, xaxis)

    # --------------- boundary condition value
    u_D = np.array(0.0, dtype=ScalarType)

    # --------------- Temperature field
    xbc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, xaxis_facets), V)  # su§b defines the component, y on xaxis
    ybc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, yaxis_facets), V)  # sub defines the component x on yaxis
    # --------------- Total
    bcs = [xbc, ybc]
    bcs = [] # No boundary conditions

    # --------------- Outside surface loading ------------------#
    def outside(x):
        return np.isclose(x[0], 0)
        #return np.isclose(np.sqrt(x[0] ** 2 + x[1] ** 2), 1.)

    outside_f = mesh.locate_entities_boundary(domain, fdim, outside)
    marked_facets = np.hstack([outside_f])
    marked_values = np.hstack([np.full_like(outside_f, 1)])
    sorted_facets = np.argsort(marked_facets)
    outside_ft = mesh.meshtags(domain, fdim, marked_facets[sorted_facets], marked_values[sorted_facets])

    f = fem.Constant(domain, ScalarType(200.))
    ds = ufl.Measure("ds", domain=domain, subdomain_data=outside_ft)
    n = ufl.FacetNormal(domain)

    # --------------- Setting up  ------------------#

    xdmf = io.XDMFFile(domain.comm, "Resultfiles/Temperature.xdmf", "w")
    xdmf.write_mesh(domain)

    temp = fem.Function(V)
    temp.name = "Temperature"
    temp.interpolate(initial_condition)
    xdmf.write_function(temp, tstart)

    # --------------- Variational formulation ------------------#
    T, dT = ufl.TrialFunction(V), ufl.TestFunction(V)
    Told = fem.Function(V)
    Told.name = "Incremental temperature"
    Told.interpolate(initial_condition)

    # --------------- Problem formulation ------------------#

    # a = cV * T * (dT) * ufl.dx + k * dt * ufl.dot(ufl.grad(T), ufl.grad(dT)) * ufl.dx
    # L = u_n * dT * ufl.dx + f * dT * ds

    # therm_form = (cV * (dT - Thetaold) / dt * T +
    #                   kappa * T0 * ufl.tr(eps(du - uold)) / dt * T +
    #                   ufl.dot(k * ufl.grad(dT), ufl.grad(T))) * ufl.dx
    # therm_form = (cV * (dT - Thetaold) / dt * T + ufl.dot(k * ufl.grad(dT), ufl.grad(T))) * ufl.dx

    #bilinear_form = fem.form(a)
    #linear_form = fem.form(L)

    t = tstart
    for i in range(num_steps):
        # Updating solution time
        t += dt

        a = cV * T * dT * ufl.dx + k * dt * ufl.dot(ufl.grad(T), ufl.grad(dT)) * ufl.dx
        L = T0 * dT * ufl.dx + f * dT * ds

        a = cV * T * dT * ufl.dx + k * dt * ufl.dot(ufl.grad(T), ufl.grad(dT)) * ufl.dx
        L = cV * T0 * dT * ufl.dx + f * dT * ds

        a = cV * T * dT * ufl.dx + k * dt * ufl.dot(ufl.grad(T), ufl.grad(dT)) * ufl.dx
        L = (dt * f) * dT * ds

        # setting up problem and solver
        problem = fem.petsc.LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        temp = problem.solve()
        temp.name = "Temperature"
        # Updating the temperature vector with the ith solution before next iteration
        #Told.x.array[:] = temp.x.array

        # Writing into result file
        xdmf.write_function(temp, t)

    return
    A = fem.petsc.assemble_matrix(bilinear_form, bcs=bcs)
    A.assemble()
    b = fem.petsc.create_vector(linear_form)

    solver = PETSc.KSP().create(domain.comm)
    solver.setOperators(A)
    solver.setType(PETSc.KSP.Type.PREONLY)
    solver.getPC().setType(PETSc.PC.Type.LU) # Linear solver

    # Loops over timesteps
    for i in range(num_steps):
        t += dt

        # Update the right hand side reusing the initial vector
        with b.localForm() as loc_b:
            loc_b.set(0)
        fem.petsc.assemble_vector(b, linear_form)

        # Apply Dirichlet boundary condition to the vector
        fem.petsc.apply_lifting(b, [bilinear_form], [bcs])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        fem.petsc.set_bc(b, bcs)

        # Solve linear problem
        solver.solve(b, uh.vector)
        uh.x.scatter_forward()

        # Update solution at previous time step (u_n)
        u_n.x.array[:] = uh.x.array

        # Write solution to file
        xdmf.write_function(uh, t)
    xdmf.close()