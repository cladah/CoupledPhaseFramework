import docker
def rundocker(modelvar):
    if modelvar['Coupling'] == 0:
        client = docker.from_env()
        #container = client.containers.run('dolfinx/dolfinx:stable', ["python3", "FEMFCSx.py"], detach=True,
                                            #auto_remove=True,
                                            #volumes=['C:/Users/ClasD/Documents/GitHub/coupledfenicsx:/root/shared'],
                                            #working_dir='/root/shared')
        container = client.containers.run('dolfinx/dolfinx:stable', ["python3", "Testsolver.py"],
                                          detach=True,
                                          auto_remove=True,
                                          volumes=['C:/Users/ClasD/Documents/GitHub/CoupledPhaseFramework/IntegratedSimulation/Solvers:/root/shared'],
                                          working_dir='/root/shared',
                                          name='fenicscxcont')
        for log in container.logs(stream=True, stdout=True, stderr=True):
            print(log)
    else:
        raise KeyError('Coupled solvers not implemented')