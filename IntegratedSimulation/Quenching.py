#from IntegratedSimulation.Solvers.SphericalFEM import *

def modelfitting(model,x,y):
    if model == "KM":
        return 1
    elif model == "JMAK":
        return 1
    elif model == "Devol":
        return 1
    else:
        raise KeyError("Fitting model wrong or not implemented")


def runquenching(modelvar):



    a = 1
    if a==1:
        pass
    else:
        if modelvar["programs"]["FEM"] == "1Dsphere":
            #sphericalsolver(1,modelvar)
            pass

