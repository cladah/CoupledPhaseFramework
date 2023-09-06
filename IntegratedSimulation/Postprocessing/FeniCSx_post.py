def runpost_FCSx(uh,Th,T0):
    import ufl
    import os
    from HelpFile import read_input
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

    def post_stress(uh, Th, T0):
        import numpy as np
        def eps(v):
            return ufl.sym(ufl.grad(v))

        def sig(v, T, T0):
            return 2.0 * mu * eps(v) + lmbda * ufl.tr(eps(v)) * ufl.Identity(len(v)) - (3 * lmbda + 2 * mu) * alpha * (
                    T - 800.) * ufl.Identity(len(v))  # alpha*(3*lmbda+2*mu)
        return sig(uh, Th, T0)
    print(post_stress(uh, Th, T0))
    return post_stress(uh, Th, T0)