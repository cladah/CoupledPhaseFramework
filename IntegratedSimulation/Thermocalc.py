# Thermomecanical
# Phases, Composition
import numpy as np
from tc_python import *

def TCequalibrium(modelvar,type):
    if type == "env":
        database = "SSUB6"
        material = setmaterial(modelvar,type)
        dependentmat = material[0]
        composition = material[1]
        phases = ["GAS", "C_S"]
        dormantphases = ["C_S"]
        referencestates = {"C": "C_S", "N": "GAS"}
    elif type == "mat":
        database = "TCFE12"
        dependentmat = modelvar["dependentmat"]
        composition = modelvar["composition"]
        phases = ["FCC_A1","FCC_A1#2","GAS","GRAPHITE_A9"]
        dormantphases = ["GAS","GRAPHITE_A9"]
        referencestates = {"C":"Graphite_A9","N":"GAS"}
    else:
        raise KeyError

    with TCPython() as start:
        # create and configure a single equilibrium calculation
        calculation = (
            start
            .select_database_and_elements(database, [dependentmat] + list(composition))
            .get_system()
            .with_single_equilibrium_calculation()
            .set_condition(ThermodynamicQuantity.temperature(), modelvar["CNtemp"])
            .set_phase_to_suspended('*')
            #.disable_global_minimization()
        )

        for element in composition:
            calculation.set_condition(ThermodynamicQuantity.mass_fraction_of_a_component(element), composition[element]/100)
        for phase in phases:
            calculation.set_phase_to_entered(phase)
        for phase in dormantphases:
            calculation.set_phase_to_dormant(phase)
        for element in referencestates:
            calculation.with_reference_state(element,referencestates[element])
        calc_result = (calculation
                       .calculate()  # Aktiverar beräkningen
                       )
        activityC = calc_result.get_value_of(ThermodynamicQuantity.activity_of_component('C'))
        activityN = calc_result.get_value_of(ThermodynamicQuantity.activity_of_component('N'))
        return activityC, activityN

def TCcarbonitriding(modelvar,activityair):
    with TCPython() as session:
        system = (session
                  .select_thermodynamic_and_kinetic_databases_with_elements("TCFE12", "MOBFE7", [modelvar["dependentmat"]] + list(modelvar["composition"]))
                  .without_default_phases().select_phase("FCC_A1").select_phase("GAS").select_phase("FCC_A1#2").select_phase("CEMENTITE_D011").select_phase("GRAPHITE_A9")
                  .get_system())
        austenite = Region("Austenite")
        austenite.set_width(modelvar["radius"])
        austenite.with_grid(CalculatedGrid.geometric()
                       .set_no_of_points(modelvar["nodesThermo"])
                       .set_geometrical_factor(0.95))
        austenite.add_phase("FCC_A1")
        austenite.add_phase("FCC_A1#2")
        austeniteprofile = CompositionProfile()
        for element in modelvar["composition"]:
            austeniteprofile.add(element, ElementProfile.constant(modelvar["composition"][element]))
        austenite.with_composition_profile(austeniteprofile)
        calculation = (system
                      .with_isothermal_diffusion_calculation()
                      .with_reference_state("N", "GAS").with_reference_state("C", "GRAPHITE_A9")
                      .set_temperature(modelvar["CNtemp"])
                      .set_simulation_time(modelvar["CNtime"] * 3600)
                      .with_right_boundary_condition(BoundaryCondition.mixed_zero_flux_and_activity()
                                                     .set_activity_for_element('C', str(activityair[0]))
                                                     .set_activity_for_element('N', str(activityair[1])))
                      .with_spherical_geometry().remove_all_regions()
                      .add_region(austenite))


        result = calculation.calculate()
        mass_frac = list()
        for element in modelvar["composition"]:
            distance = []
            distance, mass_frac_temp = result.get_mass_fraction_of_component_at_time(element, SimulationTime.LAST)
            mass_frac.append(mass_frac_temp)

        return distance, mass_frac

def setmaterial(modelvar,type):
    if type == "mat":
        return modelvar["dependentmat"],modelvar["compositon"]
    elif type == "env":
        pass
    else:
        raise KeyError("TCcalculation error")

    if modelvar["CNenv"] =="Argon":
        print('Using argon as atmosphere')
        env_dep = 'N'
        env_comp = {'H': 1, 'C': 1, 'O': 1}
    elif modelvar["CNenv"] == "Methane":
        print('Using methane as atmosphere')
        env_dep = 'N'
        env_comp = {'H': 4.5788, 'C': 13.6344, 'O': 18.1614}
        activityair = [1.471, 0.639]
    elif modelvar["CNenv"] =='Propane':
        print('Using propane as atmosphere')
        env_dep = 'N'
        env_comp = {'H': 3.2762, 'C': 14.0238, 'O': 18.6801}
    else:
        raise KeyError("Wrong carbonitiding environment, use Argon, Methane, or Propane")

    return env_dep, env_comp