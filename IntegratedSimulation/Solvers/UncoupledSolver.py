from MeshConverter import convert_mesh
from UC_FEMsolver import runFEM
from UC_Heatsolver import runHeat
from ThermoMech import runThermoMech
from UC_Phasesolver import runPhase
from Linearelast import runelast
from HelpFile import read_input


#convert_mesh()
runFEM()
#runHeat()
#runDEMO()
#runThermoMech()
#runPhase()