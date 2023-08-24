from Solvers.RunDocker import rundocker
from Solvers.ComsolSolver import runComsol
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
    if modelvar['programs']['FEM'] == 'FCSx':
        rundocker(modelvar)
    elif modelvar['programs']['FEM'] == 'Comsol':
        runComsol()
        pass
    else:
        raise KeyError('Program not implemented')
