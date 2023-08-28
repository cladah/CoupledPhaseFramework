# CoupledPhaseFramework
Version 0.1

Simulation of a quenching process implemented with Thermocalc and FeniCSx/COMSOL.

**Python 3 modules to run simulation**
- numpy -
- MPh - 1.2.3
- tc_python - 2023b
- gmsh - 4.1.1

**Running FeniCSx through Docker**

Docker image 95b47fc536a5 - dolfinx/dolfinx:stable - FeniCSx 0.6

**Running COMSOL through MPh for JAVA compatibility**

MPh - 1.2.3

**Running ThermoCalc through TC-Python**

ThermoCalc 2023b

# Calculation map
- Read input
- Create Mesh
- CCT calculation (Composition)->(Psi(r))
- Diffusion of C and N (Mesh, Composition)->(Composition(r))
