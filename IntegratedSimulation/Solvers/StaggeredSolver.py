from dolfinx import mesh, fem, io
from dolfinx.fem import FunctionSpace, VectorFunctionSpace
from dolfinx.io import gmshio
import ufl
import numpy as np
from petsc4py.PETSc import ScalarType
from mpi4py import MPI
from HelpFile import read_input
from Postprocessing.Stress import post_stress





def runFEM():
    indata = read_input()
    # --------------- Loading mesh ------------------#
    domain, cell_markers, facet_markers = gmshio.read_from_msh("Resultfiles/Mesh.msh", MPI.COMM_WORLD, 0, gdim=2)
    cord = ufl.SpatialCoordinate(domain)
    dx = ufl.Measure('dx',domain=domain)
    n = ufl.FacetNormal(domain)

    # --------------- Creating functionspaces ------------------#
    # Displacement, Heat, Phase
    Vu = VectorFunctionSpace(domain, ("Lagrange", indata["FEM"]["element_f"]), dim=indata["FEM"]["elementdim"])
    VT = FunctionSpace(domain, ("Lagrange", indata["FEM"]["element_f"]))
    Vpsi = FunctionSpace(domain, ("Lagrange", indata["FEM"]["element_f"]))

    # --------------- Input assignment ------------------#
    tstart = 0  # Start time
    tstop = indata["Thermo"]["quenchtime"]  # Final time
    num_steps = indata["Thermo"]["quench_steps"]
    dt = tstop / num_steps  # time step size

    E = fem.Constant(domain, indata["material"]["Austenite"]["E"])
    nu = fem.Constant(domain, indata["material"]["Austenite"]["E"])
    rho_g = fem.Constant(domain, ScalarType(indata["material"]["rho"])) * 9.82
    Cp_rho = fem.Constant(domain, indata["material"]["Cp"]) * fem.Constant(domain, ScalarType(indata["material"]["rho"]))
    k = fem.Constant(domain, indata["material"]["k"])
    mu = E / 2 / (1 + nu)
    lmbda = E * nu / (1 + nu) / (1 - 2 * nu)
    alpha = fem.Constant(domain, indata["material"]["alpha_k"])
    Tstart = indata["Thermo"]["CNtemp"]

    kM = fem.Constant(domain, 1e-5)
    Ms = fem.Constant(domain, 400.)
    # --------------- Assigning initial conditions ------------------#
    def initial_condition(x):
        return np.logical_or(np.isreal(x[0]), np.isreal(x[1])) * Tstart

    T0 = fem.Function(VT)
    T0.name = "Initial temperature"
    T0.interpolate(initial_condition)

    # --------------- Finding boundary conditions ------------------#

    def yaxis(x):
        return np.isclose(x[0], 0.)

    def xaxis(x):
        return np.isclose(x[1], 0.)

    def outside(x):
        return np.isclose(np.sqrt(x[0] ** 2 + x[1] ** 2), 1.)

    fdim = 1
    ndim = 0
    yaxis_facets = mesh.locate_entities_boundary(domain, fdim, yaxis)
    xaxis_facets = mesh.locate_entities_boundary(domain, fdim, xaxis)
    outside_facets = mesh.locate_entities_boundary(domain, fdim, outside)

    # --------------- boundary condition value
    u_D = np.array(0.0, dtype=ScalarType)
    T_D = np.array(200., dtype=ScalarType)

    # --------------- Displacement field
    xbc = fem.dirichletbc(u_D, fem.locate_dofs_topological(Vu.sub(1), fdim, xaxis_facets),
                          Vu.sub(1))  # sub defines the component, y on xaxis
    ybc = fem.dirichletbc(u_D, fem.locate_dofs_topological(Vu.sub(0), fdim, yaxis_facets),
                          Vu.sub(0))  # sub defines the component x on yaxis

    # --------------- Heat field
    bcT = fem.dirichletbc(T_D, fem.locate_dofs_topological(VT, fdim, outside_facets), VT)  # sub defines the component, y on xaxis

    # --------------- Total
    bcu = [xbc, ybc]
    bcT = [bcT]

    # --------------- Loading ------------------#

    outside_f = mesh.locate_entities(domain, 1, outside)
    marked_values = np.hstack([np.full_like(outside_f, 1)])
    outside_ft = mesh.meshtags(domain, fdim, outside_f, marked_values)

    fu = fem.Constant(domain, 0.)
    fT = fem.Constant(domain, 1.)
    ds = ufl.Measure("ds", domain=domain, subdomain_data=outside_ft)
    n = ufl.FacetNormal(domain)

    # --------------- Strain and stress definitions ------------------#
    def eps(v):
        return ufl.sym(ufl.grad(v))

    def sig(v, T):
        return 2.0 * mu * eps(v) + (lmbda * ufl.nabla_div(v) - alpha*(3*lmbda+2*mu)*T) * ufl.Identity(len(v))

    # --------------- Variational formulation ------------------#
    n = ufl.FacetNormal(domain)
    u, du = ufl.TrialFunction(Vu), ufl.TestFunction(Vu)
    T, dT = ufl.TrialFunction(VT), ufl.TestFunction(VT)
    psi, dpsi = ufl.TrialFunction(Vpsi), ufl.TestFunction(Vpsi)

    Delta_T = fem.Function(VT, name="Temperature increase")
    Told = fem.Function(VT)
    Told.x.array[:] = T0.x.array
    uh = fem.Function(Vu)
    uh.x.array[:] = np.zeros(len(Told.x.array[:])*2)
    uh.name = "Displacement"
    Told.name = "Temperature"
    xdmf = io.XDMFFile(domain.comm, "Resultfiles/Result.xdmf", "w")
    time = 0
    xdmf.write_mesh(domain)
    xdmf.write_function(uh, time)
    xdmf.write_function(Told, time)
    for i in range(num_steps):
        time += dt
        # --------------- Thermal problem ------------------#
        FT = (Cp_rho/dt * T * dT + k * ufl.dot(ufl.grad(T), ufl.grad(dT)) - Cp_rho/dt * Told * dT)*ufl.dx
        aT, LT = ufl.lhs(FT), ufl.rhs(FT)
        problem_T = fem.petsc.LinearProblem(aT, LT, bcs=bcT, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        Th = problem_T.solve()
        Told.x.array[:] = Th.x.array # Assigning the function values to the Ti-1 function

        # --------------- Mechanical problem ------------------#
        F = ufl.inner(sig(u, Th), eps(du)) * ufl.dx - ufl.dot(fu * n, du) * ds
        au, Lu = ufl.lhs(F), ufl.rhs(F)
        problem_u = fem.petsc.LinearProblem(au, Lu, bcs=bcu, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        uh = problem_u.solve()

        # --------------- Phase problem ------------------#
        def Koistinen(Th):
            return 1 - ufl.exp(-kM * (Ms - Th))

        Fpsi = Koistinen(Th)
        apsi, Lpsi = ufl.lhs(Fpsi), ufl.rhs(Fpsi)
        problem_psi = fem.petsc.LinearProblem(apsi, Lpsi, bcs=[], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        psih = problem_psi.solve()

        # --------------- Postprocessing and saving ------------------#
        von_Mises = post_stress(uh, sig(uh, Th))
        V_von_mises = fem.FunctionSpace(domain, ("DG", 1))
        stress_expr = fem.Expression(von_Mises, V_von_mises.element.interpolation_points())
        stresses = fem.Function(V_von_mises)
        stresses.interpolate(stress_expr)
        # Write solution to file
        uh.name = "Displacement"
        Th.name = "Temperature"
        stresses.name = "von_Mises"
        xdmf.write_function(uh, time)
        xdmf.write_function(Th, time)
        xdmf.write_function(stresses, time)

    xdmf.close()

runFEM()
