import docker
import os
from HelpFile import read_input
import time

import subprocess
def rundocker():
    directory = os.getcwd()
    directory = directory.strip('\\Solvers')
    dockervolume = directory + ':/root/shared'
    dockervolume = dockervolume.replace('\\', '/')

    data = read_input()
    if data['Programs']['Coupling'] == "Stagg":
        client = docker.from_env()
        print('Creating docker container')
        container = client.containers.run('dolfinx/dolfinx:stable', ["python3", "Solvers/StaggeredSolver.py"],
                                          detach=True,
                                          auto_remove=True,
                                          #tty=True,
                                          #stdin_open=True,
                                          volumes=[dockervolume],
                                          working_dir='/root/shared',
                                          environment=['PYTHONPATH=/usr/local/lib/python3/dist-packages:/usr/local/dolfinx-real/lib/python3.10/dist-packages:/usr/local/lib:'],
                                          name='fenicscxcont')

        print('Running FeniCSx')
        for log in container.logs(stream=True, stdout=True, stderr=True):
            print(log)
        #container.stop()
        print('FeniCSx calculation done!')
        print('Removing docker container')

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