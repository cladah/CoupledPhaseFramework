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

    rho = fem.Constant(domain, 2700.)
    alpha = fem.Constant(domain, 2.31e-5)
    k = fem.Constant(domain, 237e-6)
    kappa = alpha * (2 * mu + 3 * lmbda)
    cV = fem.Constant(domain, 910e-6) * rho
    T0 = fem.Constant(domain, 293.)

    # --------------- Formulating mixed element ------------------#
    V_el = ufl.VectorElement("CG", domain.ufl_cell(), 2) # Displacement
    W_el = ufl.FiniteElement("CG", domain.ufl_cell(), 1) # Heat
    P_el = ufl.FiniteElement("CG", domain.ufl_cell(), 1) # Phases
    Mix_el = ufl.MixedElement([V_el, W_el, P_el])
    V = fem.FunctionSpace(domain, Mix_el)

    num_subs = V.num_sub_spaces
    spaces = []
    maps = []
    for i in range(num_subs):
        space_i, map_i = V.sub(i).collapse()
        spaces.append(space_i)
        maps.append(map_i)

    # --------------- Boundary conditions ------------------#
    def yaxis(x):
        return np.isclose(x[0], 0)

    def xaxis(x):
        return np.isclose(x[1], 0)

    fdim = domain.topology.dim - 1
    yaxis_facets = mesh.locate_entities_boundary(domain, fdim, yaxis)
    xaxis_facets = mesh.locate_entities_boundary(domain, fdim, xaxis)

    u_D = np.array(0.0, dtype=ScalarType)
    T_D = np.array(100.0, dtype=ScalarType)
    # --------------- Displacement field
    xbc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, xaxis_facets), V.sub(0).sub(0))  # sub defines the component
    ybc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, yaxis_facets), V.sub(0).sub(1))  # sub defines the component
    # --------------- Temperature field
    Txbc = fem.dirichletbc(T_D, fem.locate_dofs_topological(V, fdim, xaxis_facets), V.sub(1))  # sub defines the component
    # --------------- Phase field
    Pxbc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, xaxis_facets), V.sub(2))  # sub defines the component
    # --------------- Total
    bcs = [xbc, ybc, Txbc, Pxbc]

    # --------------- Loading ------------------#
    f = fem.Constant(domain, (PETSc.ScalarType(100), PETSc.ScalarType(100)))
    Toutside = fem.Constant(domain, (PETSc.ScalarType(10), PETSc.ScalarType(10)))

    # --------------- Variational formulation ------------------#
    u, T, psi = ufl.TrialFunctions(V)
    du, dT, dpsi = ufl.TestFunctions(V)

    # Defining loading function
    Uold = fem.Function(V)
    uold, Thetaold, psiold = ufl.split(Uold)
    #w_n = fem.Function(V)
    #u_n, p_n, r_n = ufl.split(w_n)

    def eps(u):
        return ufl.sym(ufl.grad(u))
    def sig(v):
        """Return an expression for the stress σ given a displacement field"""
        return 2.0 * mu * eps(v) + lmbda * ufl.nabla_div(v) * ufl.Identity(len(v))

    def eps_th(u):
        return k * ufl.inner(ufl.grad(T), ufl.grad(dT)) * ufl.dx

    def sigtest(u,Theta):
        return (lmbda * ufl.tr(eps(u)) - kappa * Theta) * ufl.Identity(len(u)) + 2 * mu * eps(u)

    # --------------- Problem formulation ------------------#
    #a = ufl.inner(sig(u), eps(du)) * ufl.dx
    #L = ufl.dot(f, du) * ufl.dx
    #F = ufl.inner(sig(u), eps(du)) * ufl.dx - ufl.dot(f, du) * ufl.dx
    #a = fem.form(ufl.lhs(F))
    #L = fem.form(ufl.rhs(F))
    dt = fem.Constant(domain, 0.1)
    mech_form = ufl.inner(sigtest(du, dT), eps(u)) * ufl.dx + ufl.inner(f, du) * ufl.dx
    therm_form = (cV * (dT - Thetaold) / dt * T + kappa * T0 * ufl.tr(eps(du - uold)) / dt * T + ufl.dot(k * ufl.grad(dT), ufl.grad(T))) * ufl.dx
    F = mech_form + therm_form
    a = fem.form(ufl.lhs(F))
    L = fem.form(ufl.rhs(F))

    a = fem.form(ufl.inner(sig(u), eps(du)) * ufl.dx)
    L = fem.form(ufl.dot(f, du) * ufl.dx)

    # Left hand side
    uh = fem.Function(V)
    A = fem.petsc.assemble_matrix(a, bcs=bcs)
    A.assemble()

    # Right hand side
    b = fem.petsc.assemble_vector(L)
    fem.petsc.apply_lifting(b, [a], [bcs])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
    fem.petsc.set_bc(b, bcs)

    # Solver
    solver = PETSc.KSP().create(domain.comm)
    solver.setOperators(A)
    #solver.setType(PETSc.KSP.Type.PREONLY)
    solver.setType(PETSc.KSP.Type.BCGS)
    solver.getPC().setType(PETSc.PC.Type.LU)

    Problem = fem.petsc.NonlinearProblem(F, uh, bcs=bcs)
    solver = nls.petsc.NewtonSolver(mesh.comm, Problem)
    # Absolute tolerance
    solver.atol = 5E-10
    # relative tolerance
    solver.rtol = 1E-11
    solver.convergence_criterion = " incremental "
    solver.solve(uh)

    a = ufl.inner(sig(u), eps(du)) * ufl.dx
    L = ufl.dot(f, du) * ufl.dx
    #solver.solve(b, uh.vector)
    #problem = fem.petsc.LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    #uh = problem.solve()
    #problem = fem.petsc.LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    #uh = problem.solve()
    disp, temp, phase = uh.split() # Splitting solution into field arrays
    disp.name = "Displacement"
    temp.name = "Temperature"
    phase.name = "Phasefraction"

    #therm_form = (cV * (dT - Thetaold) / dt * T + kappa * T0 * ufl.tr(eps(du - uold)) / dt * T + ufl.dot(k * ufl.grad(dT), ufl.grad(T))) * ufl.dx


    # F = fem.form(ufl.inner(sig(u), ufl.grad(v)) * ufl.dx - ufl.inner(f, v) * ufl.dx)
    # a = fem.form(ufl.inner(sig(u), ufl.grad(v)) * ufl.dx)
    # L = fem.form(ufl.inner(f, v) * ufl.dx)
    #A = fem.petsc.assemble_matrix(a, bcs=bcs)
    #A.assemble()

    #b = fem.petsc.assemble_vector(L)
    #fem.petsc.apply_lifting(b, [a], bcs=[[bcs]])
    #fem.petsc.set_bc(b, [bcs])
    #uh = fem.Function(V)
    #problem = fem.petsc.NonlinearProblem(F, u, bcs)
    #solver = PETSc.KSP().create(domain.comm)
    #solver.setFromOptions()
    #solver.setOperators(A)
    #solver = nls.petsc.NewtonSolver(domain.comm, problem)



    # F_heat = k * inner(grad(T), grad(q)) * dx - alpha * T * div(u) * q * dx
    # F_disp = inner(sym(nabla_grad(u)), grad(v)) * dx
    # F =F_disp + F_heat





    with io.XDMFFile(domain.comm, "Resultfiles/displacement.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(disp)
        xdmf.write_function(temp)
    #    xdmf.write_function(phase)