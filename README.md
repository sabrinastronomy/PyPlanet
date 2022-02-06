# PyPlanet 
## A rocky exoplanet modeling software with realistic temperature profiles and numerical thermal evolution
Rocky planets are diverse in structure and composition compared to the Earth. Their temperature profiles could also differ greatly from Earth depending on their mass and distance from their host stars. Interior structure models of rocky exoplanets have not yet studied the full range of possible temperature profiles. We develop an open source Python module, \texttt{\textit{PyPlanet}}, for rocky planets with an arbitrary number of layers and equations of state. We apply this model to explore many possible temperature profiles and quantify the thermal effects on the mass-radius relations of rocky planets with compositions similar to other rocky planets in our solar system, such as Mars and Mercury. This detailed modeling will be crucial for making robust inferences about rocky planet structure and composition from transit and radial velocity observations. Our planetary grid explores central pressures around that of the Earth (approximately 360 \si{\giga \pascal}) between \e{9} to \e{12} \si{\pascal}. We find that thermal profile considerations are significant for super-Earths but have an even larger effect on sub-Earth planets. A rocky planet's interior temperature profile can affect its radius by several percent and can exceed the observational uncertainties on measured exoplanet radii. 
## Dependencies
* scipy==1.3.1
* numpy==1.17.0
* burnman==0.9.0
* matplotlib==3.3.2
* astropy==3.2.1
* imageio==2.6.1
