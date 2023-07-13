from Quenching import *
from Carbonitriding import *
from Post_processing import *
import subprocess
import docker
import h5py
from HelpFile import Material
from Solvers.RunDocker import *
import time

def start():
    # Model variables stored in a dictionary
    modelvar = dict()
    # --------------- Material ------------------#
    modelvar['Coupling'] = 0                        # 0 - Uncoupled solvers

    # --------------- Material ------------------#
    modelvar["dependentmat"] = "Fe"
    modelvar["composition"] = {"C": 0.3, "Si": 0.2, "Mn": 0.3, "Cr": 1.4, "P": 0.025, "N": 0.025}

    # --------------- Geometry ------------------#
    modelvar["geotype"] = "Spherical"  # Type of geometry Spherical, Axisymmetric,
    modelvar["radius"] = 0.01  # Radius of sphere [m]
    modelvar["thickness"] = 0.01  # Thickness of cylinder [m]
    modelvar["innerradius"] = 0.005  # inner radius of cylinder [m]

    # --------------- Material ------------------#
    modelvar["lambda"] = 1E-12              #
    modelvar["c_k"] = 1                     #
    austenite = Material(E=200E9, nu=0.3, alpha=10E-6, f=1)
    martensite = Material(200E9, 0.3, 10E-6, 0)
    perlite = Material(200E9, 0.3, 10E-6, 0)
    bainite = Material(200E9, 0.3, 10E-6, 0)
    materialmix = Material(200E9, 0.3, 10E-6, 0)
    modelvar['Austenite'] = austenite
    modelvar['Martensite'] = martensite
    modelvar['Perlite'] = perlite
    modelvar['Bainite'] = bainite
    modelvar['Materialmix'] = materialmix       # This is recalculated in FEM simulation

    # --------------- Thermodynamics ------------------#
    modelvar["nodesThermo"] = 3                 # Nodes in thermodynamic calculation
    modelvar["CNtemp"] = 800 + 273.15           # Temperature for carbonitriding/annealing degC
    modelvar["CNenv"] = 'Methane'               # Type of gas environment for carbonitriding
    modelvar["temperingtemp"] = 400 + 273.15    # Temperature of tempering
    modelvar["quenchtime"] = 1                  # Time in quenching medium [h]
    modelvar["quench_dt"] = 1/500               # Time step for quenching simulation
    modelvar["quenchmedium"] = 'oil'            # Type of quenching medium

    # --------------- FEM ------------------#
    modelvar["elementtype"] = "Spherical"   # Spherical, AxisymPstrain, AxisymPstress, Axisym2D
    modelvar["nodesFEM"] = 4                # Nodes radial FEM calculation
    modelvar["nodesFEM_phi"] = 4            # Nodes circumferential FEM calculation
    modelvar["nodesFEM_z"] = 4              # Nodes circumferential FEM calculation
    modelvar["Meshscaling"] = 0.95          # scaling towards circumference
    modelvar["shapef"] = "Linear"           # FEM shapefunction Linear, Quad
    modelvar["peclet"] = abs(2E-9 * modelvar["radius"] / modelvar["nodesFEM"] / 2 / modelvar["lambda"])

    # --------------- Material models ------------------#
    models = dict()
    models['Martensite'] = 'KM'             # Koitinen-Marburger
    models['Perlite'] = 'JMAK'              # JMAK, Devold
    models['Banite'] = 'JMAK'               # JMAK, Devold
    models['Ferrite'] = 'JMAK'              # JMAK, Devold
    models['CombinedMat'] = 'RuleofMix'     # RuleofMix,
    models['Mattype'] = 'Elastic'           # Elastic, Elasto-plastic,
    modelvar["models"] = models

    # --------------- Program choices ------------------#
    programs = dict()
    programs['Diff'] = 'TC'                 # TC, Comsol, FCSx_1Dsph, FCSx_1Daxi
    programs['CCT'] = 'TC'                  # TC, Comsol, FCSx_1Dsph, FCSx_1Daxi
    programs['Therm'] = 'FCSx'              # Comsol, FCSx
    programs['FEM'] = 'FCSx'                # Comsol, FCSx
    programs['Phasetr'] = 'FCSx'            # Comsol, FCSx
    modelvar["programs"] = programs

    # --------------- Time ------------------#
    modelvar["CNtime"] = 1                  # Time for carbonitration [h]
    modelvar["quenchingtime"] = 60 / 3600   # Quenching calculation time [h]
    modelvar["temperingtime"] = 1           # time at tempering temperature [h]

    # --------------- Info ------------------#
    modelvar["info"] = [modelvar["dependentmat"], modelvar["composition"], modelvar["nodesFEM"]]

    # --------------- Input ------------------#
    inputvariable = input('What do you want to do? ')

    if inputvariable == 'Run':
        checkcalculation()
        createmesh(modelvar)
        runcarbonitriding(modelvar)
        #runquenching(modelvar)
        # runtempering(modelvar)
        # runfatiguetest(modelvar)
        # runplot(modelvar)

    elif inputvariable == '1D':
        if not testchacheinfo(modelvar):
            print("Modelinput changed recalculation matrixes")

        RHS = 1
        #sphericalsolver(RHS, modelvar)

        # --------------- Adding new info to cache ------------------#
        savetocache("info",
                    [modelvar["radius"], modelvar["nodesFEM"], modelvar["E"], modelvar["nu"], modelvar["shapef"],
                     modelvar["Meshscaling"]])

    elif inputvariable == "Test":
        createresultfile(modelvar)
        createcachfile(modelvar)
        createmesh(modelvar)

        #[displacement, err] = sphericalsolver(1,modelvar)
        saveresult("displacement",displacement)
        plotstrain(modelvar)
        plotdisplacement(modelvar)
        #strains = getstrains(modelvar, results['displacement'])
        #sigma = getstresses(modelvar, strains)
        #createResults(modelvar, results['displacement'])
        #plotstrain(modelvar)

    elif inputvariable == 'y':
        rundocker(modelvar)

if __name__ == "__main__":
    start()