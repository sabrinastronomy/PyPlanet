# About PyPlanet 

PyPlanet is a rocky exoplanet modeling software with realistic temperature profiles and numerical thermal evolution. It was developed by Sabrina Berger while an undergraduate at UC Berkeley and into her time as a graduate student at McGill University alongsider her supervisor Professor Leslie Rogers. Jisheng Zhang and Ellen Price also contributed parts of this code. 

# PyPlanet Science Goals (abstract)

ocky planets are diverse in structure and composition compared to the Earth. Their temperature profiles could also differ greatly from Earth depending on their mass and distance from their host stars. Interior structure models of rocky exoplanets have not yet studied the full range of possible temperature profiles. We develop an open source Python module, \texttt{\textit{PyPlanet}}, for rocky planets with an arbitrary number of layers and equations of state. We apply this model to explore many possible temperature profiles and quantify the thermal effects on the mass-radius relations of rocky planets with compositions similar to other rocky planets in our solar system, such as Mars and Mercury. This detailed modeling will be crucial for making robust inferences about rocky planet structure and composition from transit and radial velocity observations. Our planetary grid explores central pressures around that of the Earth (approximately 360 \si{\giga \pascal}) between \e{9} to \e{12} \si{\pascal}. We find that thermal profile considerations are significant for super-Earths but have an even larger effect on sub-Earth planets. A rocky planet's interior temperature profile can affect its radius by several percent and can exceed the observational uncertainties on measured exoplanet radii.

# Requirements and Installation
Dependencies
scipy==1.3.1
numpy==1.17.0
burnman==0.9.0
matplotlib==3.3.2
astropy==3.2.1
imageio==2.6.1

# What can I do with PyPlanet?

PyPlanet gives the following many overall planetary properties of the planet:
-radius
-mass
-core mass fraction
-core radius fraction
-total energy
-simulated transition pressures and surface pressures to determine over or underestimations of a layer's boundary.

You can also interpolate to other masses and radii beyond your simulation by using the interpolation feature. Here's how the interpolation works:
1. flatten all grids
2. dictionary of central pressures is created where keys are all unique p_c and values are element indices
3. gathers all values of parameters for p_c dictionary values
4. interpolate between p_cmb/p_c ratio and all parameters, interpolate between CMFs and p_cmb/p_c ratio
5. get value of each parameter at p_cmb/p ratio value for desired CMF
6. collect u_value for given anchor temperature into separate array
Now you have function for U(T_s) for fixed CMF value 
surface temp for given CMF and mass for a bunch 
for each anchor temperature
get planet radius and U (grav binding energy)
numerical integration of dT/dt

# Description of Individual .py files
-run.py
-planet.py
-eos.py
-planet_grid.py
-therm_ev.py
-mass_radius_relations.py
-interpCMF.py
## \Plotting
-T_equilibrium.py
-plotting.py

# Examples
