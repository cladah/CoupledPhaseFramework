from mpi4py import MPI
from dolfinx import mesh, fem, io
from dolfinx.fem import FunctionSpace
from dolfinx.io import gmshio
import ufl
import numpy as np
#import gmsh
from petsc4py.PETSc import ScalarType


domain = mesh.create_interval(MPI.COMM_WORLD, nx=50, points=(0.0, 1.0))

# Function space over domain
V = fem.FunctionSpace(domain, ("Lagrange", 1)) #Lagrange

material = dict()
i = 0
with open('Material.txt') as f:
    lines = f.readlines()
    for line in lines:
        line = line.split()
        material[line[0]] = [float(line[1]), float(line[2]), float(line[3]), float(line[4])]
    i += 1
    f.close()

# Material parameters
E = fem.Constant(domain, material['Austenite'][0])
nu = fem.Constant(domain, material['Austenite'][1])
rho_g = 1e-3
mu = E/2/(1+nu)
lmbda = E*nu/(1+nu)/(1-2*nu)

#center = mesh.locate_entities_boundary(domain, dim=0, marker=lambda x: np.isclose(x[0], 0.0))
#circumference = mesh.locate_entities_boundary(domain, dim=0, marker=lambda x: np.isclose(x[0], 1.0))
#dofs = fem.locate_dofs_topological(V=V, entity_dim=0, entities=center)


def in_center(x):
    return np.isclose(x[0], 0)

def on_circumference(x):
    return np.isclose(x[0], 1)

#Inner bc
center_ent = mesh.locate_entities_boundary(domain, dim=0,marker=lambda x: np.isclose(x[0], 0.0))
center_dofs = fem.locate_dofs_topological(V=V, entity_dim=0, entities=center_ent)
bc = fem.dirichletbc(value=ScalarType(0), dofs=center_dofs, V=V)

#Circumference bc
circumference_ent = mesh.locate_entities_boundary(domain, dim=0, marker=lambda x: np.isclose(x[0], 1.0))
circumferencedofs = fem.locate_dofs_topological(V=V, entity_dim=0, entities=circumference_ent)
#circumference_ds = ufl.Measure("dp", subdomain_data=circumferencedofs)

p = fem.Constant(domain, ScalarType(-10))


f = fem.Function(V)
dofs = fem.locate_dofs_geometrical(V, lambda x: np.isclose(x.T, 1.0))
f.x.array[dofs] = 100

x = ufl.SpatialCoordinate(domain)

def eps(u):
    return ufl.sym(ufl.grad(u))


def axisym_eps(v):
    return ufl.sym(ufl.as_tensor([[v[0].dx(0), 0, v[0].dx(1)],
                          [0, v[0]/x[0], 0],
                          [v[1].dx(0), 0, v[1].dx(1)]]))
def sphere_eps(v):
    #return ufl.grad(v)
    return ufl.as_tensor([[v.dx(0), 0],[0, v/x[0]]])

# Stress given a displacement field
def sigma(u):
    return lmbda*ufl.tr(sphere_eps(u))*ufl.Identity(2) + 2.0*mu*sphere_eps(u)


# Gravity on the full domain
#f = fem.Constant(domain, ScalarType((0, 0, -rho_g)))


# Variational problem
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
# Defining the problem
#a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
#L = p * v * ufl.dx
a = (ufl.inner(sigma(u), sphere_eps(v))*x[0]**2)*ufl.dx
L = (ufl.inner(f, v)*x[0]**2)*ufl.dx

# Setting up problem
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

# Heat solver
alpha = 3
beta = 1.2
u_D = ufl.Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t',
                 degree=2, alpha=alpha, beta=beta, t=0)
VT = fem.FunctionSpace(domain, ('DG', 1))
bcT = fem.DirichletBC(VT, u_D, boundary)




with io.XDMFFile(domain.comm, "output.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_function(uh)
    xdmf.write_function(vm_stresses)
    #xdmf.write_function(stresses)