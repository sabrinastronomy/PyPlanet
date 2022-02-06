from planet_grid import PlanetGrid
from plotting import *
import numpy as np
import burnman.minerals

testing = True

location = "/Users/sabrinaberger/RockyPlanets"

minfractions_normal = [[1], [[0.9, 0.1, 0.0], [0.8, 0.2]]]
minfractions_hot = [[1], [0.9, 0.1, 0.0]]

# core molar fractions corresponding to:
# formula = 'FeSiO3', i.e., Murakami_2013.fe_perovskite

# (lower) mantle molar fractions corresponding to
# [mg_perovskite(), '[Mg][Si]O3'], [fe_perovskite(), '[Fe][Si]O3'], [al_perovskite(), '[Al][Al]O3']
# i.e., SLB_2011.mg_fe_silicate_perovskite

# (upper) mantle molar fractions corresponding to
# 'Mg2Si2O6', i.e., SLB_2011.enstatite

# BURNMAN MATERIALS USED WITHIN PLANET (not including liquid silicates)
core_material = burnman.minerals.Murakami_2013.fe_perovskite()
mantle_material = burnman.minerals.SLB_2011.mg_fe_silicate_perovskite()
upper_mantle_material = burnman.minerals.SLB_2011.enstatite()

layers_normal = [[core_material], [mantle_material, upper_mantle_material]]
layers_hot = [[core_material], [mantle_material]]

thermal_location = location + "/complete_data_with_silicate_mantle"
test_eos_location = location + "/testEoSData"

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
        if temp > 2500:
            temp_plan_grid = PlanetGrid(temp, central_pressures, grid_size, str(type_eos), loc, layers_types=layers_hot, minfractions=minfractions_hot)
        else:
            temp_plan_grid = PlanetGrid(temp, central_pressures, grid_size, str(type_eos), loc, layers_types=layers_normal, minfractions=minfractions_normal)
        xx, yy, radius_grid, mass_grid, press_grid, core_mass_grid, core_rad_grid, u_grid = temp_plan_grid.integrateGrid()
        planetary_grids.append(temp_plan_grid)
    np.save("samp_planet_rad.npy", radius_grid)
    np.save("samp_planet_mass.npy", mass_grid)
    np.save("samp_planet_press.npy", press_grid)
    np.save("samp_planet_u.npy", u_grid)


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
    default_central_pressures = [12, 12]
    default_grid_size = [2, 2]
    temp_range = np.linspace(300, 3000, 10)
    # print(temp_range)
    temp_range = [2600]
    # adiabatic
    varying_temp("_adiabatic_", temp_range, default_central_pressures, default_grid_size, thermal_location)
    # constant temperature
    # varying_temp("_constant_", temp_range, default_central_pressures, default_grid_size, thermal_location)


