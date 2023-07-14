import docker
def rundocker(modelvar):
    austenite = modelvar['Austenite']
    martensite = modelvar['Martensite']
    perlite = modelvar['Perlite']
    bainite = modelvar['Bainite']
    f = open("Material.txt", "w")
    f.write('Austenite ' + str(austenite.E)+' ' + str(austenite.nu)+' ' + str(austenite.alpha) + ' ' + str(austenite.f) + '\n')
    f.write('Martensite ' + str(martensite.E)+' ' + str(martensite.nu)+' ' + str(martensite.alpha) + ' ' + str(martensite.f) + '\n')
    f.write('Perlite ' + str(perlite.E)+' ' + str(perlite.nu)+' ' + str(perlite.alpha) + ' ' + str(perlite.f) + '\n')
    f.write('Bainite ' + str(bainite.E)+' ' + str(bainite.nu)+' ' + str(bainite.alpha) + ' ' + str(bainite.f) + '\n')
    f.close()

    material = dict()
    i = 0
    with open('Material.txt') as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            print(line[0])
            material[line[0]] = [float(line[1]), float(line[2]), float(line[3]), float(line[4])]
        i += 1
        f.close()

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