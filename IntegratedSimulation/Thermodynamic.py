def modelfitting(model,x,y):
    if model == "KM":
        return 1
    elif model == "JMAK":
        return 1
    elif model == "Devol":
        return 1
    else:
        raise KeyError("Fitting model wrong or not implemented")
def interpolateThermMod():
    pass


def Koistinen(Th, V, Ms, beta):
    from dolfinx import fem
    import ufl
    return fem.Expression(1 - ufl.exp(-beta * (Ms - Th)), V.element.interpolation_points())

def JMAK(Th, V, K, A):
    pass

