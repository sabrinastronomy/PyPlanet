from planet_grid import PlanetGrid
from eos import *
import subprocess
from mass_radius_relations_in_progress import *
from plotting import *

testing = True

location = "/Users/sabrinaberger/All Research/RockyPlanets/without_phases/"
default_central_pressures = [9, 12]
default_grid_size = [2, 2]

all_data_location = location + "/PyPlanetData"
grid_planets_location = all_data_location + "/gridPlanets"

if not testing:
    subprocess.call(["mkdir", all_data_location])
    subprocess.call(["mkdir", grid_planets_location])

    def varying_temp(type_eos, temperatures, central_pressures, grid_size):
        planetary_grids = []

        for temp in temperatures:
            temp_plan_grid = PlanetGrid(temp, central_pressures, grid_size, str(type_eos), all_data_location)
            temp_plan_grid.integrateGrid()
            planetary_grids.append(temp_plan_grid)


        # TODO currently unable to save all grids in pickle file
            #  np.save(grid_planets_location + "/" + type_eos + str(temp) + "_.pyc", planetary_grids)


    if __name__ == "__main__":
        temp_range = [300, 3000]

        # adiabatic
        # varying_temp("_adiabatic_", temp_range, default_central_pressures, default_grid_size)
        # constant temperature
        # varying_temp("_constant_", temp_range, default_central_pressures, default_grid_size)

        planet_interp(location + "/DataFiles/", "_adiabatic_", temp_range)


location = location + "/DataFiles/"

for temp in [300]:
    for type in ['Constant']:
        p_c = location + "p_c_" + type + '_' + str(temp) +".pyc.npy"
        p_cmb_percentage = location + "p_cmb_percentage_" + type + '_' + str(temp) +".pyc.npy"
        radius = location + "radius_grid_" + type + '_' + str(temp) +".pyc.npy"
        mass = location + "mass_grid_" + type + '_' + str(temp) +".pyc.npy"
        surf_press = location + "press_grid_" + type + '_' + str(temp) +".pyc.npy"
        core_mass = location + "core_mass_grid_" + type + '_' + str(temp) +".pyc.npy"
        core_rad = location + "core_rad_grid_" + type + '_' + str(temp) +".pyc.npy"
        p_cmb_simulated = location + "p_cmb_simulated_" + type + '_' + str(temp) +".pyc.npy"

        radius_higher = location + "radius_grid_" + type + '_' + str(3000) +".pyc.npy"
        PlotBasic(temp, type, p_c, p_cmb_percentage, radius, mass, surf_press, core_mass, core_rad, p_cmb_simulated, radius_higher)

# Testing MT EoS - in progress

# import os
# import sys
# # hack to allow scripts to be placed in subdirectories next to burnman:
# if not os.path.exists('burnman') and os.path.exists('../burnman'):
#     sys.path.insert(1, os.path.abspath('..'))
#
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.integrate import odeint
# from scipy.integrate import quad
# from scipy.interpolate import interp1d
# from scipy.interpolate import UnivariateSpline
#
# import burnman
# import burnman.minerals as minerals
#
#
# class iron(burnman.Mineral):
#     def __init__(self):
#         # Parameters for gamma - Fe are from Tsujino et al. 2013
#         # G_0 G_Prime0 from Mao et al. 2001 (fig. 3)
#         self.params = {
#             'equation_of_state': 'mt',
#             'T_0': 1273.,
#             'V_0': 7.381e-06,
#             'K_0': 111.5e9,
#             'Kdprime_0': 5.2,
#             'Kprime_0': 5.2,
#             'G_0': 83.2e9,  # Shear modulus and derivative from Gleason and Mao, 2013
#             'Gprime_0': 2.04,
#             'molar_mass': 55.845 / 1000.,
#             'n': 1,
#             'Debye_0': 340.,
#             'grueneisen_0': 2.28,
#             'q_0': 0.21,
#             'eta_s_0': 2.0  # Wholly invented value
#         }
#         burnman.Mineral.__init__(self)
#
# p_cmb = 10e8
# p_c = 10e12
# pressures = np.linspace(0, p_c, 20).tolist()
#
# # temperatures_mantle = geotherm.adiabatic(pressures_mantle, anchor_temperature,
# #
# #
# # c                                        self.mantle.composite)  # first pressure is 0
#
# core = iron()
# # temperatures_core = geotherm.adiabatic(pressures_core, 300, core)  # first pressure is p_cmb
# temperatures_core = np.full(1000, 200)
# print(temperatures_core)
# core = core.evaluate(['density'], pressures, temperatures_core)
