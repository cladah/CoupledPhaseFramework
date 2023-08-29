from Solvers.RunDocker import rundocker
from Solvers.ComsolSolver import runComsol
from HelpFile import read_input

def runquenching():
    data = read_input()
    if data['Programs']['FEM'] == 'FCSx':
        rundocker()
    elif data['Programs']['FEM'] == 'Comsol':
        runComsol()
        pass
    else:
        raise KeyError('Program not implemented')
