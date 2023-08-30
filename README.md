# CoupledPhaseFramework
Version 0.1

Simulation of a quenching process implemented with Gmsh, Thermocalc, and FeniCSx/COMSOL.

**Python 3 modules to run simulation**
- gmsh - 4.1.1
- numpy -
- MPh - 1.2.3
- tc_python - 2023b


**Running FeniCSx through Docker**

Docker image 95b47fc536a5 - dolfinx/dolfinx:stable - FeniCSx 0.6

**Running COMSOL through MPh for JAVA compatibility**

MPh - 1.2.3

**Running ThermoCalc through TC-Python**

ThermoCalc 2023b

# Calculation map
- Read input
- Create Mesh
- Activity of C and N calculation (Composition, Temperature) -> (aC, aN)
- Diffusion of C and N (Mesh, Composition)->(Composition(r))
- CCT calculation (Composition(r))->(CCT(T,t,r,dT))
- Estimate model parmeters (CCT(T,t,r,dT)) -> (Modelparameters(r))
- Time dependent solver, t = range(0,tmax,numstep):
  - Convergance analysis, err < maxerr:
  - Heat solver (Ti(r), dt) -> (Ti+1(r), eps_th(r))
  - Phase solver (Ti+1(r), dt, psi_ji(r) sigi+1(r)) -> (psi_ji+1(r), eps_psi(r))
    - CCT interpolation(Ti+1(r), dt, psi_ji(r) sigi+1(r)) -> (psi_ji+1(r))
  - Solid mech solver (eps_th(r), eps_psi(r), Ti+1(r)) -> (sig(r), eps(r), eps_pl(r))
