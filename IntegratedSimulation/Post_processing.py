from dolfinx import mesh, fem, io
from mpi4py import MPI
import ufl
from HelpFile import read_input
def runpostprocess():
    indata = read_input()

    domain, cell_markers, facet_markers = io.gmshio.read_from_msh("Resultfiles/Mesh.msh", MPI.COMM_WORLD, 0, gdim=2)
    xdmf = io.XDMFFile(domain.comm, "Resultfiles/Result.xdmf", "w")

    # Material
    E = fem.Constant(domain, indata["material"]["Austenite"]["E"])
    nu = fem.Constant(domain, indata["material"]["Austenite"]["nu"])
    rho_g = fem.Constant(domain, indata["material"]["rho"]) * 9.82
    Cp_rho = fem.Constant(domain, indata["material"]["Cp"]) * fem.Constant(domain, indata["material"]["rho"])
    k = fem.Constant(domain, indata["material"]["k"])
    mu = E / 2 / (1 + nu)
    lmbda = E * nu / (1 + nu) / (1 - 2 * nu)
    alpha = fem.Constant(domain, indata["material"]["alpha_k"])




    def post_stress(uh, Th, T0, time):
        def eps(v):
            return ufl.sym(ufl.grad(v))
        def sig(v, T, T0):
            return 2.0 * mu * eps(v) + lmbda * ufl.tr(eps(v)) * ufl.Identity(len(v)) - (3 * lmbda + 2 * mu) * alpha * (
                        T - 800.) * ufl.Identity(len(v))

        s = sig(uh, Th, T0)
        V_stress = fem.TensorFunctionSpace(domain, ("DG", 1))
        # von_Mises = post_stress(uh, sig(uh, Th))
        # V_von_mises = fem.FunctionSpace(domain, ("DG", 1))
        stress_expr = fem.Expression(s, V_stress.element.interpolation_points())
        stresses = fem.Function(V_stress)
        stresses.interpolate(stress_expr)
        stresses.name = "Stress"
        stresses.sub(0).name = "xx"
        stresses.sub(1).name = "xy"
        xdmf.write_function(stresses, time)
        return
    def post_strain(uh):
        def eps(v):
            return ufl.sym(ufl.grad(v))

        epsilon = eps(uh)
        V_eps = fem.TensorFunctionSpace(domain, ("DG", 1))
        strain_expr = fem.Expression(epsilon, V_eps.element.interpolation_points())
        strains = fem.Function(V_eps)
        strains.interpolate(strain_expr)
        strains.name = 'Strain'
        return