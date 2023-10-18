from Solvers.RunDocker import rundocker
from Solvers.ComsolSolver import runComsol
from HelpFile import read_input, checkinput, adjustinputcache

def runquenching():
    if checkinput('Quenching'):
        print('Using old quenching simulation')
        return
    print('Quenching module')
    data = read_input()
    if data['Programs']['FEM'] == 'FCSx':
        print('Using FeniCSx for FEM calculation')
        rundocker()
    elif data['Programs']['FEM'] == 'Comsol':
        print('Using COMSOL for FEM calculation')
        runComsol()
    else:
        raise KeyError('Program not implemented')
    adjustinputcache('Mesh')