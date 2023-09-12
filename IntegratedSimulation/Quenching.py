from Solvers.RunDocker import rundocker
from Solvers.ComsolSolver import runComsol
from HelpFile import read_input

def runquenching():
    print('Quenching module')
    data = read_input()
    if data['Programs']['FEM'] == 'FCSx':
        print('Using FeniCSx for FEM calculation')
        rundocker()
    elif data['Programs']['FEM'] == 'Comsol':
        print('Using COMSOL for FEM calculation')
        runComsol()
        pass
    else:
        raise KeyError('Program not implemented')
