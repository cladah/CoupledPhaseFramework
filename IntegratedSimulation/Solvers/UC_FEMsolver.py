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
    gmsh.initialize()
    gmsh.model.add("QuarterCirc")
    gdim = 2
    gmsh.model.occ.addDisk(0, 0, 0, 1, 1)
    gmsh.model.occ.addRectangle(0, 0, 0, 1, 1, 2)
    gmsh.model.occ.intersect([(gdim, 1)], [(gdim, 2)], 3)
    gmsh.model.occ.synchronize()

    gmsh.model.addPhysicalGroup(gdim, [3], 1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.02)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.03)
    gmsh.model.mesh.generate(gdim)
    gmsh_model_rank = 0
    mesh_comm = MPI.COMM_WORLD
    domain, cell_markers, facet_markers = gmshio.model_to_mesh(gmsh.model, mesh_comm, gmsh_model_rank, gdim=gdim)

    # --------------- Loading input ------------------#
    E = 1.0e9
    nu = 0.3
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))

    # --------------- Formulating mixed element ------------------#
    V = fem.VectorFunctionSpace(domain, ("CG", 2), dim=2)

    # --------------- Boundary conditions ------------------#
    def yaxis(x):
        return np.isclose(x[0], 0.)

    def xaxis(x):
        return np.isclose(x[1], 0.)
    def corner(x):
        return np.isclose(x[0]+x[1],0.)

    fdim = domain.topology.dim - 1
    yaxis_facets = mesh.locate_entities_boundary(domain, fdim, yaxis)
    xaxis_facets = mesh.locate_entities_boundary(domain, fdim, xaxis)
    corner_f = mesh.locate_entities(domain, fdim, corner)

    u_D = np.array(0.0, dtype=ScalarType)
    # --------------- Displacement field
    xbc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, xaxis_facets), V.sub(1))  # sub defines the component
    cornbc = fem.dirichletbc(np.array((0.0,0.0), dtype=ScalarType), fem.locate_dofs_topological(V, fdim, corner_f), V)  # sub defines the component
    #ybc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, yaxis_facets), V.sub(0))  # sub defines the component
    # --------------- Total
    #bcs = [xbc, ybc]
    bcs = [cornbc]

    # --------------- Loading ------------------#
    f = fem.Constant(domain, (1000., 1000.))

    # --------------- Variational formulation ------------------#
    u = ufl.TrialFunction(V)
    du = ufl.TestFunction(V)

    def eps(v):
        return ufl.sym(ufl.grad(v))
    def sig(v):
        return 2.0 * mu * eps(v) + lmbda * ufl.nabla_div(v) * ufl.Identity(len(v))

    # --------------- Problem formulation ------------------#

    #F = ufl.inner(sig(du), eps(u)) * ufl.dx - ufl.inner(f, du) * ufl.dx
    #a = fem.form(ufl.lhs(F))
    #L = fem.form(ufl.rhs(F))

    a = ufl.inner(sig(u), eps(du)) * ufl.dx
    L = ufl.dot(f, du) * ufl.dx

    # Left hand side
    problem = fem.petsc.LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    disp = problem.solve()
    disp.name = "Displacement"


    with io.XDMFFile(domain.comm, "Resultfiles/displacement.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(disp)