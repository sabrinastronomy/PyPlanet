from planet_grid import PlanetGrid
from plotting import *
import numpy as np

testing = True

location = "/Users/sabrinaberger/RockyPlanets"


thermal_location = location + "/thermalData"
test_eos_location = location + "/testEoSData"

# try:
#     subprocess.call(["mkdir", all_data_location])
#     subprocess.call(["mkdir", grid_planets_location])

def varying_temp(type_eos, temperatures, central_pressures, grid_size, loc):
    """
    :param type_eos:
    :param temperatures:
    :param central_pressures:
    :param grid_size:
    :param loc:
    :return:

    Creates planets grids with varying anchor temperatures
    """
    planetary_grids = []

    for temp in temperatures:
        temp_plan_grid = PlanetGrid(temp, central_pressures, grid_size, str(type_eos), loc)
        temp_plan_grid.integrateGrid()
        planetary_grids.append(temp_plan_grid)


if __name__ == "__main__":

    # TEST EOS
    # default_central_pressures = [8, 14]
    # default_grid_size = [1, 1]
    # temp_range = [300]
    # # adiabatic
    # varying_temp("_adiabatic_", temp_range, default_central_pressures, default_grid_size, test_eos_location)
    # # constant temperature
    # varying_temp("_constant_", temp_range, default_central_pressures, default_grid_size, test_eos_location)
    #

    # THERMAL EVOLUTION
    default_central_pressures = [9, 12]
    default_grid_size = [4, 4]
    temp_range = np.linspace(300, 3000, 100)
    # adiabatic
    varying_temp("_adiabatic_", temp_range, default_central_pressures, default_grid_size, thermal_location)
    # constant temperature
    # varying_temp("_constant_", temp_range, default_central_pressures, default_grid_size, thermal_location)


