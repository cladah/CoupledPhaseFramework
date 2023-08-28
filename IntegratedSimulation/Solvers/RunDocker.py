import docker
import os
from HelpFile import read_input
def rundocker():
    directory = os.getcwd()
    directory = directory.strip('\\Solvers')
    dockervolume = directory + ':/root/shared'
    dockervolume = dockervolume.replace('\\', '/')

    data = read_input()

    if data['Programs']['Coupling'] == 0:
        client = docker.from_env()
        container = client.containers.run('dolfinx/dolfinx:stable', ["python3", "Solvers/StaggeredSolver.py"],
                                          detach=True,
                                          auto_remove=True,
                                          volumes=[dockervolume],
                                          working_dir='/root/shared',
                                          name='fenicscxcont')
        for log in container.logs(stream=True, stdout=True, stderr=True):
            print(log)
        print('Fenicsx calculation done!')
    else:
        raise KeyError('Solver not implemented')
        client = docker.from_env()
        container = client.containers.run('dolfinx/dolfinx:stable', ["python3", "Solvers/CoupledSolver.py"],
                                          detach=True,
                                          auto_remove=True,
                                          volumes=[dockervolume],
                                          working_dir='/root/shared',
                                          name='fenicscxcont')
        for log in container.logs(stream=True, stdout=True, stderr=True):
            print(log)
        print('Fenicsx calculation done!')