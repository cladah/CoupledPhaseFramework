import docker
def rundocker():
    client = docker.from_env()
    for container in client.containers.list():
      container.stop()

    fenicscont = client.containers.run('dolfinx/dolfinx:stable',
        #name='test',
        auto_remove=True,
        detach=True,
        entrypoint='bash',
        volumes=['C:/Users/ClasD/Documents/Docker/CoupledPhaseFramework:/root/'],
        working_dir='/root/',
        command=['echo FEMFCS.py'])

    #fenicscont.exec_run(cmd=['python3', 'my-code.py'])

    for line in fenicscont.logs(stream=True):
        print(str(line))