from HelpFile import read_input
import pyvista as pv

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

    def post_stress(disp, Th, T0, time):
        import numpy as np
        def eps(v):
            return np.gradient(v)
        def sig(v, T, T0):
            return 2.0 * mu * eps(v) + lmbda * np.trace(eps(v)) * np.identity(len(v)) - (3 * lmbda + 2 * mu) * alpha * (
                        T - 800.) * np.identity(len(v))
        print(np.array(disp).shape)
        #s = sig(disp, T, T0)
        #print(sig(disp, Th, T0))


    with h5py.File(directory + '/Result.h5', "r") as f:
        geometry = f['Mesh']['mesh']['geometry'][...]
        disp = f['Function']['Displacement']
        T = f['Function']['Temperature']
        post_stress(disp['1000'][...], T['1000'][...], T['0'], float('1000'))
        #for x in disp.keys():





    #post_stress(disp['100'][...], T['100'][...], T['0'][...], 100)
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

def plotPyVista():
    def load_time(value):
        """ Load solution at value specified using the slider """
        reader.set_active_time_value(value)
        grid = reader.read()[2]
        p.add_mesh(grid, cmap=cmap)

    #directory = os.getcwd() + '/Resultfiles/Result.xdmf'
    #directory = directory.replace('\\', '/')
    reader = pv.get_reader('Resultfiles/Result.xdmf')
    reader.set_active_time_value(100.0)
    mesh = reader.read()[2]
    #mesh.set_active_scalars("Displacement")

    cmap = "plasma"

    p = pv.Plotter()
    p.add_mesh(mesh, cmap=cmap)
    p.view_xy()
    p.add_slider_widget(load_time, [0.0, 1000.0], value=0.0, title="Time",
                        interaction_event="always", pointa=(0.25, 0.93),
                        pointb=(0.75, 0.93))
    p.show()