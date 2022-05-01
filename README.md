# About *PyPlanet* 

*PyPlanet* is a rocky exoplanet modeling software with realistic temperature profiles and numerical thermal evolution. It was developed by Sabrina Berger while she was an undergraduate at UC Berkeley and into her time as a graduate student at McGill University. This code was developed alongside her supervisor Professor Leslie Rogers. Jisheng Zhang and Ellen Price also contributed parts of this code. 

# *PyPlanet* Science Goals

Rocky planets are diverse in structure and composition compared to the Earth. Their temperature profiles could also differ greatly from Earth depending on their mass and distance from their host stars. Interior structure models of rocky exoplanets have not yet studied the full range of possible temperature profiles. We develop the following open source Python module, *PyPlanet*, for rocky planets with an arbitrary number of layers and equations of state. We apply this model to explore many possible temperature profiles and quantify the thermal effects on the mass-radius relations of rocky planets with compositions similar to other rocky planets in our solar system, such as Mars and Mercury. This detailed modeling will be crucial for making robust inferences about rocky planet structure and composition from transit and radial velocity observations. Our planetary grid explores central pressures around that of the Earth (approximately 360 GPa) between 10^9 and 10^12 Pa. We find that thermal profile considerations are significant for super-Earths but have an even larger effect on sub-Earth planets. A rocky planet's interior temperature profile can affect its radius by several percent and can exceed the observational uncertainties on measured exoplanet radii.

# Requirements and Installation
Dependencies:
- scipy==1.3.1
- numpy==1.17.0
- burnman==0.9.0
- matplotlib==3.3.2
- astropy==3.2.1
- imageio==2.6.1

# What can I do with *PyPlanet*?

*PyPlanet* gives the following many overall planetary properties of the planet:
- radius
- mass
- core mass fraction
- core radius fraction
- total energy
- simulated transition pressures and surface pressures to determine over or underestimations of a layer's boundary.

You can also interpolate to other masses and radii beyond your simulation by using the interpolation feature.

# Description of Individual .py files
## run.py
This wraps planet_grid.py such that you can integrate multiple planetary grids at an array of temperatures. These profiles can be either constant or adiabatic.

## planet_grid.py
This generates planetary grids of arbitrary size according to the user's specifications. 

## eos.py
This generates the equations of state for the planet. The equations of state are generated with Burnman except for the molten silicate upper mantle. This region pulls from a pre-generated table of pressures, temperatures, heat capacities, and densities.


## planet.py
This does the main integration of the planet given the parameters set in planet_grid.py and the equation of state generated in eos.py.

## therm_ev.py 
This computes the thermal evolution of a planet.
## mass_radius_relations.py
This generates the mass radius relationship of planets given a particular grid generated in planet_grid.py.
## interpCMF.py
This interpolates to various core mass fractions that the user specifies from a planetary grid.
# Plotting
## T_equilibrium.py
This generates an equilibrium temperature plot as a function of mass and radius.
## plotting.py
This provides code to plot a grid of planets. An example of what's plotted for a grid of 300K constant planets is shown below:
(I can't figure out how to put an inline figure hold up ugh)