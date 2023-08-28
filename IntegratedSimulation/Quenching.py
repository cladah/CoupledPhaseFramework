from Solvers.RunDocker import rundocker
from Solvers.ComsolSolver import runComsol
from HelpFile import read_input
def modelfitting(model,x,y):
    if model == "KM":
        return 1
    elif model == "JMAK":
        return 1
    elif model == "Devol":
        return 1
    else:
        raise KeyError("Fitting model wrong or not implemented")


def runquenching():
    data = read_input()
    if data['Programs']['FEM'] == 'FCSx':
        rundocker()
    elif data['Programs']['FEM'] == 'Comsol':
        runComsol()
        pass
    else:
        raise KeyError('Program not implemented')
