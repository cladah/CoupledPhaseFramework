from mpi4py import MPI
from dolfinx import mesh, fem, io
from dolfinx.fem import FunctionSpace, Function
from dolfinx.io import gmshio
import ufl
import numpy as np
from petsc4py.PETSc import ScalarType
import gmsh
from ufl import dx, ds, grad, inner, div

# --------------- Defining domain/mesh ------------------#
# domain, cell_markers, facet_markers = gmshio.read_from_msh("Resultfiles/Mesh.msh", MPI.COMM_WORLD, 0, gdim=3)
domain = mesh.create_interval(MPI.COMM_WORLD, nx=50, points=(0.0, 1.0))
el = ufl.FiniteElement("CG", 1)
mixed_el = ufl.MixedElement([el, el, el])
V = FunctionSpace(domain, mixed_el)

# --------------- Material parameters ------------------#

E = fem.Constant(domain, 1e5)
nu = fem.Constant(domain, 0.3)
rho_g = 1e-3
mu = E/2/(1+nu)
lmbda = E*nu/(1+nu)/(1-2*nu)
alpha = 2.31e-5

# --------------- Boundary conditions ------------------#

# Thermal
g_expr = '1 + x[0]*x[0] + beta *t'
g = ufl.Expression(g_expr, alpha=3.0, beta=1.2, t=0, degree=2)
u_n = fem.Function(V)
u_n.name = "u_n"
u_n.interpolate(g)

#Inner bc
center_ent = mesh.locate_entities_boundary(domain, dim=0,marker=lambda x: np.isclose(x[0], 0.0))
center_dofs = fem.locate_dofs_topological(V=V, entity_dim=0, entities=center_ent)
bc = fem.dirichletbc(value=ScalarType(0), dofs=center_dofs, V=V)
circ_ent = mesh.locate_entities_boundary(domain, dim=0,marker=lambda x: np.isclose(x[0], 1.0))
circ_dofs = fem.locate_dofs_topological(V=V, entity_dim=0, entities=center_ent)
bcT = fem.dirichletbc(value=g, dofs=circ_dofs, V=V) # Thermal bc

# Force at end node
f = fem.Function(V)
dofs = fem.locate_dofs_geometrical(V, lambda x: np.isclose(x.T, 1.0))
f.x.array[dofs] = 100

# --------------- Trail- and test-functions ------------------#

# Displacement (u,v)
# Heat (p,q)
# Phase Martensite (f1,e1)
(u, p, f1) = ufl.TrialFunctions(V)
(v, q, e1) = ufl.TestFunctions(V)

# --------------- Defining problem ------------------#

# Domain coordinates
x = ufl.SpatialCoordinate(domain)


# Strain given displacement field
def sphere_eps(v):
    # return ufl.grad(v)
    return ufl.as_tensor([[v.dx(0), 0], [0, v / x[0]]])


# Stress given a displacement field
def sigma(u):
    return lmbda * ufl.tr(sphere_eps(u)) * ufl.Identity(2) + 2.0 * mu * sphere_eps(u)

# Displacement
au = (ufl.inner(sigma(u), sphere_eps(v))*x[0]**2)*ufl.dx
Lu = (ufl.inner(f, v)*x[0]**2)*ufl.dx


# Heat
dt = fem.Constant(domain, 0.3)
t = float(dt)
T = 1.8

fT = fem.Constant(domain, 1.2 - 2.0 - 2*3.0)

# Previous and current solution
u0 = ufl.interpolate(g, V)
u1 = Function(V)

aT = p * q * dx + dt * ufl.dot(grad(p), grad(q)) * dx
LT = (u0 + dt * fT) * q * dx


# --------------- Solving problem (Displacement) ------------------#

solver = dolfinx.fem.petsc.LinearProblem(
    au, Lu, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu",
                              "pc_factor_mat_solver_type": "mumps"})
uh = solver.solve()
uh.name = "displacement"


# --------------- Solving problem (Heat equation) ------------------#
t_i = 0
t_f = 100
Nsteps = 50

time = np.linspace(t_i, t_f, Nsteps+1)

bilinear_form = fem.form(aT)
linear_form = fem.form(LT)
A = fem.petsc.assemble_matrix(bilinear_form, bcs=[bcT])
A.assemble()
b = fem.petsc.create_vector(linear_form)

solver = PETSc.KSP().create(domain.comm)
solver.setOperators(A)
solver.setType(PETSc.KSP.Type.PREONLY)
solver.getPC().setType(PETSc.PC.Type.LU)

for i in range(Nsteps):
    t += dt

    # Update the right hand side reusing the initial vector
    with b.localForm() as loc_b:
        loc_b.set(0)
    fem.petsc.assemble_vector(b, linear_form)

    # Apply Dirichlet boundary condition to the vector
    fem.petsc.apply_lifting(b, [bilinear_form], [[bcT]])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
    fem.petsc.set_bc(b, [bcT])

    # Solve linear problem
    solver.solve(b, uh.vector)
    uh.x.scatter_forward()

    # Update solution at previous time step (u_n)
    u_n.x.array[:] = uh.x.array

    xdmf.write_function(uh, t)
xdmf.close()