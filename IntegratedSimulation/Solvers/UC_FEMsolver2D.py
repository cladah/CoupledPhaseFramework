from dolfinx import mesh, fem, io
from dolfinx.fem import FunctionSpace, VectorFunctionSpace
from dolfinx.io import gmshio
import ufl
import numpy as np
from petsc4py.PETSc import ScalarType
from mpi4py import MPI

def runFEM():
    # --------------- Loading mesh ------------------#
    #with io.XDMFFile(MPI.COMM_WORLD, "Resultfiles/Mesh.xdmf", "r") as file:
    #    msher = file.read_mesh(name="Sphere")
    #    cell_markers = file.read_meshtags(msher, name="Sphere")
    #    facet_markers = file.read_meshtags(msher, name="Sphere")
    msh, cell_markers, facet_markers = gmshio.read_from_msh("Resultfiles/Mesh.msh", MPI.COMM_WORLD, 0, gdim=2)
    domain = msh
    cord = ufl.SpatialCoordinate(domain)
    dx = ufl.Measure('dx',domain=domain)
    n = ufl.FacetNormal(domain)

    # --------------- Input assignment ------------------#

    # Material parameters #--------------------------------- Read from file saved under /Cachefiles
    E = fem.Constant(domain, 1e5)
    nu = fem.Constant(domain, 0.3)
    rho_g = 1e-3
    mu = E / 2 / (1 + nu)
    lmbda = E * nu / (1 + nu) / (1 - 2 * nu)

    # --------------- Creating functionspaces ------------------#

    V = VectorFunctionSpace(domain, ("Lagrange", 1), dim=2)

    # --------------- Finding boundary conditions ------------------#

    def ysym_boundary(x):
        return np.isclose(x[0], 0)

    def xsym_boundary(x):
        return np.isclose(x[1], 0)

    fdim = domain.topology.dim - 1
    ysym = mesh.locate_entities_boundary(domain, fdim, ysym_boundary)
    xsym = mesh.locate_entities_boundary(domain, fdim, xsym_boundary)

    u_D = np.array(0.0, dtype=ScalarType)
    xbc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, xsym), V.sub(0))  # sub defines the component
    ybc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, ysym), V.sub(1))  # sub defines the component
    bc = [xbc, ybc]

    # --------------- Tractions on outside edge ------------------#
    def outer_bound(x):
        return np.isclose(np.sqrt(x[0]**2+x[1]**2), 1)

    #outer_facets = mesh.locate_entities_boundary(domain, 1, outer_bound)
    #ds = ufl.Measure("ds", subdomain_data=ysym)
    p = fem.Constant(domain, ScalarType((0, 1000)))

    bottom_facets = mesh.locate_entities_boundary(domain, 1, lambda x: np.isclose(x[1], 0))
    ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_markers.find(0))

    # --------------- Strain and stress definitions ------------------#
    def eps(u):
        return ufl.sym(ufl.grad(u))

    def axisym_eps(v):
        return ufl.sym(ufl.as_tensor([[v[0].dx(0), 0, v[0].dx(1)],
                              [0, v[0]/x[0], 0],
                              [v[1].dx(0), 0, v[1].dx(1)]]))

    # Stress given a displacement field
    def sigma(u):
        return lmbda*ufl.tr(eps(u))*ufl.Identity(2) + 2.0*mu*eps(u)

    # --------------- Variational formulation ------------------#
    n = ufl.FacetNormal(domain)
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    #print(thermal_initial(u).ufl_shape)
    #print(eps(v).ufl_shape)

    a = (ufl.inner(sigma(u), eps(v)))*ufl.dx # axisymmetric definitions
    L = (ufl.inner(p, v))*ds # axisymmetric definitions


    # --------------- Problem formulation ------------------#

    problem = fem.petsc.LinearProblem(a, L, bcs=bc, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    uh = problem.solve()
    uh.name = "displacement"
    print(uh.ufl_shape)
    ux,uy = uh.split()
    print(ux.ufl_shape)
    # --------------- Postprocessing ------------------#

    #s = sigma(uh) -1/3*ufl.tr(sigma(uh))*ufl.Identity(2)
    #von_Mises = ufl.sqrt(3./2*ufl.inner(s, s))
    #V_vm = fem.FunctionSpace(domain, ("DG", 1))
    #stress_expr = fem.Expression(von_Mises, V_vm.element.interpolation_points())
    #vm_stresses = fem.Function(V_vm, name="vm_stress")
    #vm_stresses.interpolate(stress_expr)
    #V_stress = fem.VectorFunctionSpace(domain, ("DG", 1))
    #stress_expr = fem.Expression(s, V_stress.element.interpolation_points())
    #stresses = fem.Function(V_stress, name="stress")
    #stresses.interpolate(stress_expr.sub(0))
    #.sub(0)

    with io.XDMFFile(domain.comm, "Resultfiles/displacement.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(uh)
        #xdmf.write_function(ux)
        #xdmf.write_function(vm_stresses)
        #xdmf.write_function(stresses)