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
def runPhase():
    # --------------- Loading mesh ------------------#

    msh, cell_markers, facet_markers = gmshio.read_from_msh("Resultfiles/Mesh.msh", MPI.COMM_WORLD, 0, gdim=2)
    domain = msh

    # --------------- Loading inputs ------------------#

    t = 0  # Start time
    T = 1.0  # Final time

    num_steps = 50
    dt = T / num_steps  # time step size

    # --------------- Creating function space ------------------#

    V = fem.FunctionSpace(domain, ("CG", 1))

    # --------------- Assigning initial conditions ------------------#
    def initial_condition(x, a=5):
        return np.exp(-a * (x[0] ** 2 + x[1] ** 2)) + 800

    u_n = fem.Function(V)
    u_n.name = "u_n"
    u_n.interpolate(initial_condition)

    u_bc = fem.Function(V)
    fdim = domain.topology.dim - 1
    boundary_facets = mesh.locate_entities_boundary(
        domain, fdim, lambda x: np.full(x.shape[1], True, dtype=bool))
    outside_facets = facet_markers.find(2)
    outside_dofs = fem.locate_dofs_topological(V, domain.topology.dim - 1, outside_facets)

    # --------------- Finding boundary conditions ------------------#
    def ysym_boundary(x):
        return np.isclose(x[0], 0)

    def xsym_boundary(x):
        return np.isclose(x[1], 0)

    def circumf_boundary(x):
        return np.isclose(np.sqrt(x[0] ** 2 + x[1] ** 2), 1)  # Insert radius

    # u_zero = np.array((0,) * mesh.geometry.dim, dtype=ScalarType)
    # bcxsym = fem.dirichletbc(u_zero, fem.locate_dofs_geometrical(V, xsym_boundary), V)
    # bcysym = fem.dirichletbc(u_zero, fem.locate_dofs_geometrical(V, ysym_boundary), V)
    bccirc = fem.dirichletbc(PETSc.ScalarType(200), fem.locate_dofs_geometrical(V, circumf_boundary), V)
    # bcs = [bcxsym, bcysym, bccirc]

    # bc = fem.dirichletbc(PETSc.ScalarType(0), outside_dofs, V)
    # bc = fem.dirichletbc(PETSc.ScalarType(200), fem.locate_dofs_topological(V, fdim, circumf_boundary), V)
    bc = bccirc

    # --------------- Setting up  ------------------#

    xdmf = io.XDMFFile(domain.comm, "Resultfiles/diffusion.xdmf", "w")
    xdmf.write_mesh(domain)

    # Define solution variable, and interpolate initial solution for visualization in Paraview
    uh = fem.Function(V)
    uh.name = "Temperature"
    uh.interpolate(initial_condition)
    xdmf.write_function(uh, t)

    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    f = fem.Constant(domain, PETSc.ScalarType(0))

    a = u * v * ufl.dx + dt * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
    L = (u_n + dt * f) * v * ufl.dx

    bilinear_form = fem.form(a)
    linear_form = fem.form(L)

    A = fem.petsc.assemble_matrix(bilinear_form, bcs=[bc])
    A.assemble()
    b = fem.petsc.create_vector(linear_form)

    solver = PETSc.KSP().create(domain.comm)
    solver.setOperators(A)
    solver.setType(PETSc.KSP.Type.PREONLY)
    solver.getPC().setType(PETSc.PC.Type.LU)

    for i in range(num_steps):
        t += dt

        # Update the right hand side reusing the initial vector
        with b.localForm() as loc_b:
            loc_b.set(0)
        fem.petsc.assemble_vector(b, linear_form)

        # Apply Dirichlet boundary condition to the vector
        fem.petsc.apply_lifting(b, [bilinear_form], [[bc]])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        fem.petsc.set_bc(b, [bc])

        # Solve linear problem
        solver.solve(b, uh.vector)
        uh.x.scatter_forward()

        # Update solution at previous time step (u_n)
        u_n.x.array[:] = uh.x.array

        # Write solution to file
        xdmf.write_function(uh, t)
    xdmf.close()