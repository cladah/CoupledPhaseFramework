from mpi4py import MPI
from dolfinx import mesh, fem, io
from dolfinx.fem import FunctionSpace, Constant
from dolfinx.io import gmshio
import ufl
import numpy as np
#import gmsh
from petsc4py.PETSc import ScalarType

def runelast():
    msh, cell_markers, facet_markers = gmshio.read_from_msh("Resultfiles/Mesh.msh", MPI.COMM_WORLD, 0, gdim=2)
    msh.name = 'Quarter Circle'

    element = ufl.VectorElement('Lagrange',mesh.ufl_cell(),degree=1,dim=2)
    V = FunctionSpace(msh, element)

    V_x = V.sub(0).collapse() # Collapse space to x and y
    V_y = V.sub(1).collapse()


    def left(x):
        return np.isclose(x[0], 0)

    yaxis = fem.locate_dofs_geometrical((V.sub(0), V_x), left)

    zero_uy = fem.Function(V_y)
    with zero_uy.vector.localForm() as bc_local:
        bc_local.set(0.0)

    zero_ux = fem.Function(V_x)
    with zero_ux.vector.localForm() as bc_local:
        bc_local.set(0.0)

    bc = fem.DirichletBC(zero_ux, yaxis, V.sub(0))


    dx = ufl.Measure("dx",domain=msh)
    Circumf = msh.locate_entities_boundary(msh, 1, lambda x : np.isclose(np.sqrt(x[0]**2+x[1]**2), 1))
    mt = msh.MeshTags(msh, 1, Circumf, 1)
    ds = ufl.Measure("ds", subdomain_data=mt)




    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    E = 1.
    nu = 0.3
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    # this is for plane-stress
    lmbda = 2*mu*lmbda/(lmbda+2*mu)

    def eps(u):
        """Strain"""
        return ufl.sym(ufl.grad(u))

    def sigma(eps):
        """Stress"""
        return 2.0 * mu * eps + lmbda * ufl.tr(eps) * ufl.Identity(2)

    def a(u,v):
        """The bilinear form of the weak formulation"""
        k = 1.e+6
        return ufl.inner(sigma(eps(u)), eps(v)) * dx

    def L(v):
        """The linear form of the weak formulation"""
        # Volume force
        b = ufl.Constant(msh,ufl.as_vector((0,0)))

        # Surface force on the top
        f = ufl.Constant(msh, ufl.as_vector((0,0.1)))
        return ufl.dot(b, v) * dx + ufl.dot(f, v) * ds(1)

    problem = fem.LinearProblem(a(u,v), L(v), bcs=bcs,
                                        petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    uh = problem.solve()
    uh.name = "displacement"

    with io.XDMFFile(msh.comm, "Resultfiles/displacement.xdmf", "w") as xdmf:
        xdmf.write_mesh(msh)
        xdmf.write_function(uh)