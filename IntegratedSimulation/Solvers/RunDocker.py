import docker
def rundocker(modelvar):
    if modelvar['Coupling'] == 0:
        if modelvar['geotype'] == 'Spherical':
            FCSxcmd = 'Spherical_Uncoupled.py'
        elif modelvar['geotype'] == 'Axisym':
            FCSxcmd = 'Axisym_uncoupled.py'
        elif modelvar['geotype'] == 'Axisym2D':
            FCSxcmd = 'Axisym2D_uncoupled.py'
        else:
            raise KeyError('Docker has not implemented geometry')

        client = docker.from_env()
        container = client.containers.run('dolfinx/dolfinx:stable', ["python3", "Testsolver.py"],
                                          detach=True,
                                          auto_remove=True,
                                          volumes=['C:/Users/ClasD/Documents/GitHub/CoupledPhaseFramework/IntegratedSimulation/Solvers:/root/shared'],
                                          working_dir='/root/shared',
                                          name='fenicscxcont')
        for log in container.logs(stream=True, stdout=True, stderr=True):
            print(log)
        print('Fenicsx calculation done!')
    else:
        raise KeyError('Coupled solvers not implemented')

        client = docker.from_env()
        container = client.containers.run('dolfinx/dolfinx:stable', ["python3", "CoupledSolver.py"],
                                          detach=True,
                                          auto_remove=True,
                                          volumes=[
                                              'C:/Users/ClasD/Documents/GitHub/CoupledPhaseFramework/IntegratedSimulation/Solvers:/root/shared'],
                                          working_dir='/root/shared',
                                          name='fenicscxcont')
        for log in container.logs(stream=True, stdout=True, stderr=True):
            print(log)
        print('Fenicsx calculation done!')