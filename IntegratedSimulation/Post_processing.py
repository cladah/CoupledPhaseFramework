from HelpFile import read_input
# import paraview
import h5py
import os
def runpostprocess():
    directory = os.getcwd() + '/Resultfiles'
    directory = directory.replace('\\', '/')
    indata = read_input()

    # Material
    E = indata["material"]["Austenite"]["E"]
    nu = indata["material"]["Austenite"]["nu"]
    rho_g = indata["material"]["rho"] * 9.82
    Cp_rho = indata["material"]["Cp"] * indata["material"]["rho"]
    k = indata["material"]["k"]
    mu = E / 2 / (1 + nu)
    lmbda = E * nu / (1 + nu) / (1 - 2 * nu)
    alpha = indata["material"]["alpha_k"]
    print(directory)
    with h5py.File(directory + '/Result.h5', "r") as f:
        geometry = f['Mesh']['mesh']['geometry'][...]
        disp = f['Function']['Displacement']
        for x in disp.keys():
            print(disp[x])



    def post_stress(uh, Th, T0, time):
        import ufl
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
        names = ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"]
        for i in range(len(names)):
            stresses.sub(i).name = names[i]
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