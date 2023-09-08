from dolfinx import mesh, fem, io
from dolfinx.fem import FunctionSpace, VectorFunctionSpace
from dolfinx.io import gmshio
import ufl
import numpy as np
from petsc4py.PETSc import ScalarType
from mpi4py import MPI
from HelpFile import read_input
#from Postprocessing.Stress import post_stress
from Postprocessing.FeniCSx_post import runpost_FCSx
from Thermodynamic import Koistinen
import os

def runsolver():
    #import h5py
    indata = read_input()
    # --------------- Loading mesh ------------------#
    domain, cell_markers, facet_markers = gmshio.read_from_msh("Resultfiles/Mesh.msh", MPI.COMM_WORLD, 0, gdim=2)
    cord = ufl.SpatialCoordinate(domain)
    dx = ufl.Measure('dx', domain=domain)
    n = ufl.FacetNormal(domain)

    # --------------- Creating functionspaces ------------------#
    # Displacement, Heat, Phase
    Vu = VectorFunctionSpace(domain, ("Lagrange", indata["FEM"]["element_f"]), dim=indata["FEM"]["elementdim"])
    VT = FunctionSpace(domain, ("Lagrange", indata["FEM"]["element_f"]))
    Vpsi = FunctionSpace(domain, ("Lagrange", indata["FEM"]["element_f"]))

    # --------------- Input assignment ------------------#
    # Time
    tstop = indata["Thermo"]["quenchtime"]  # Final time
    num_steps = indata["Thermo"]["quench_steps"]
    dt = tstop / num_steps  # time step size

    # Material
    E = fem.Constant(domain, indata["material"]["Austenite"]["E"])
    nu = fem.Constant(domain, indata["material"]["Austenite"]["nu"])
    rho_g = fem.Constant(domain, ScalarType(indata["material"]["rho"])) * 9.82
    Cp_rho = fem.Constant(domain, indata["material"]["Cp"]) * fem.Constant(domain, ScalarType(indata["material"]["rho"]))
    k = fem.Constant(domain, indata["material"]["k"])
    mu = E / 2 / (1 + nu)
    lmbda = E * nu / (1 + nu) / (1 - 2 * nu)
    alpha = fem.Constant(domain, indata["material"]["alpha_k"])
    alpha_psiM = fem.Constant(domain, indata["material"]["alpha_psiM"])
    Tstart = indata["Thermo"]["CNtemp"]
    Et = E / 100.  # tangent modulus
    EH = E * Et / (E - Et)  # hardening modulus



    # Models
    beta = fem.Constant(domain, indata["Models"]["KM"]["beta"])
    Ms = fem.Constant(domain, indata["Models"]["KM"]["Ms"])
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
        return np.isclose(np.sqrt(x[0] ** 2 + x[1] ** 2), indata["Geometry"]["radius"])

    fdim = 1
    ndim = 0
    yaxis_facets = mesh.locate_entities_boundary(domain, fdim, yaxis)
    xaxis_facets = mesh.locate_entities_boundary(domain, fdim, xaxis)
    outside_facets = mesh.locate_entities_boundary(domain, fdim, outside)

    # --------------- boundary condition value
    u_D = np.array(0.0, dtype=ScalarType)
    T_D = np.array(293.15, dtype=ScalarType)

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

    # --------------- Strain and stress definitions ------------------#
    def eps(v):
        return ufl.sym(ufl.grad(v))

    def sig(v, T, T0, psi):
        return 2.0 * mu * eps(v) + lmbda * ufl.tr(eps(v)) * ufl.Identity(len(v)) - (3*lmbda+2*mu) * (alpha * (T - 800.) - alpha_psiM*psi) * ufl.Identity(len(v)) # alpha*(3*lmbda+2*mu)

    def eps_th(T, T0):
        return alpha * (T-T0)
    # --------------- Variational formulation ------------------#
    n = ufl.FacetNormal(domain)
    u, du = ufl.TrialFunction(Vu), ufl.TestFunction(Vu)
    T, dT = ufl.TrialFunction(VT), ufl.TestFunction(VT)
    psi, dpsi = ufl.TrialFunction(Vpsi), ufl.TestFunction(Vpsi)

    # Setting up solution functions
    uh = fem.Function(Vu)
    uh.x.array[:] = np.zeros(len(uh.x.array[:]))
    Told = fem.Function(VT)
    Told.x.array[:] = T0.x.array
    psiMh = fem.Function(Vpsi)
    psiMh.x.array[:] = np.zeros(len(psiMh.x.array[:]))

    # saving 0 point
    uh.name = "Displacement"
    Told.name = "Temperature"
    psiMh.name = "Martensite fraction"
    xdmf = io.XDMFFile(domain.comm, "Resultfiles/Result.xdmf", "w")
    time = 0
    xdmf.write_mesh(domain)
    xdmf.write_function(uh, time)
    xdmf.write_function(Told, time)
    xdmf.write_function(psiMh, time)

    for i in range(num_steps):
        time += dt

        # --------------- Thermal problem ------------------#
        FT = (Cp_rho/dt * T * dT + k * ufl.dot(ufl.grad(T), ufl.grad(dT)) - Cp_rho/dt * Told * dT)*ufl.dx

        aT, LT = ufl.lhs(FT), ufl.rhs(FT)
        problem_T = fem.petsc.LinearProblem(aT, LT, bcs=bcT, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        Th = problem_T.solve()
        Told.x.array[:] = Th.x.array  # Assigning the function values to the Ti-1 function

        # --------------- Phase problem ------------------#

        psiMh = fem.Function(Vpsi)
        psiMh.interpolate(Koistinen(Th, Vpsi, Ms, beta))

        # --------------- Mechanical problem ------------------#
        F = ufl.inner(sig(u, Th, T0, psiMh), eps(du)) * ufl.dx  # - ufl.dot(fu * n, du) * ds, DT is Th-Told

        au, Lu = ufl.lhs(F), ufl.rhs(F)
        problem_u = fem.petsc.LinearProblem(au, Lu, bcs=bcu, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        uh = problem_u.solve()



        # --------------- Postprocessing and saving ------------------#
        stress = sig(uh, Th, T0, psiMh)
        s = sig(uh, Th, T0, psiMh) - 1. / 3 * ufl.tr(sig(uh, Th, T0, psiMh)) * ufl.Identity(len(uh))
        von_Mises = ufl.sqrt(3. / 2 * ufl.inner(s, s))
        #s = sih(uh) - 1. / 3 * ufl.tr(sig(uh)) * ufl.Identity(len(uh))
        #von_Mises = ufl.sqrt(3. / 2 * ufl.inner(s, s))

        Stress_V = fem.TensorFunctionSpace(domain, ("DG", 1))
        vM_V = fem.FunctionSpace(domain, ("DG", 1))
        stress_expr = fem.Expression(stress, Stress_V.element.interpolation_points())
        vM_expr = fem.Expression(von_Mises, vM_V.element.interpolation_points())
        stresses = fem.Function(Stress_V)
        stresses.interpolate(stress_expr)
        vM_f = fem.Function(vM_V)
        vM_f.interpolate(vM_expr)
        #s = runpost_FCSx(uh, Th, T0)
        #s.name = "Stress"
        vM_f.name = "von Mises"

        # Write solution to file
        uh.name = "Displacement"
        Th.name = "Temperature"
        #stresses.name = "von_Mises"
        psiMh.name = "Martensite fraction"
        xdmf.write_function(uh, time)
        xdmf.write_function(Th, time)
        xdmf.write_function(psiMh, time)
        xdmf.write_function(stresses, time)
        xdmf.write_function(vM_f, time)



    xdmf.close()
#os.system("export PYTHONPATH=/usr/local/lib/python3/dist-packages:$PYTHONPATH")
#os.system('HDF5_MPI="ON" HDF5_DIR="/usr/local" pip3 install --no-cache-dir --no-binary=h5py h5py --upgrade')
runsolver()

