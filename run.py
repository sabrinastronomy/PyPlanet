"""
This file wraps all other scripts such that you can run a full varying temperature rocky planet grid from this script.
Written by Sabrina Berger
"""

# importing packages
from planet_grid import PlanetGrid  # where the planet grid is created and run
import numpy as np  # numpy
import burnman.minerals  # relevant burnman minerals package

# This testing variable is for debugging purposes
testing = False

# This is the location where ALL output planetary grids are stored
# location = "/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/RockyPlanets/paper_data"
location = "/home/sabrinab/scratch"
thermal_location = location + "/complete_data_with_silicate_mantle"
test_eos_location = location + "/testEoSData"

# TODO add olivine information
""" 
OVERVIEW OF BURNMAN MATERIALS USED WITHIN PLANET (not including liquid silicates which are implemented in eos.py)
---------------------------------------------------------
CORE 
one layer and an iron silicate
formula = 'FeSiO3', i.e., Murakami_2013.fe_perovskite
---------------------------------------------------------
MANTLE 
lower and upper mantle 

If less than 2500 K:
LOWER
-----------------------------
If perovskite (higher pressures):
lower mantle molar fractions corresponding to
[90% mg_perovskite(), '[Mg][Si]O3'], [10 % fe_perovskite(), '[Fe][Si]O3'], [0 % al_perovskite(), '[Al][Al]O3']
i.e., SLB_2011.mg_fe_silicate_perovskite

transition between lower and upper found through olivine transition pressure

UPPER
-----------------------------
upper mantle molar fractions corresponding to
forsterite: 'Mg2Si2O6' and fayalite: '[Fe]2SiO4' i.e., SLB_2011.mg_fe_olivine
molar fractions are [0.9, 0.1]


If greater than 2500 K:
LOWER AND UPPER 
-----------------------------
upper mantle is a liquid silicate (see molten_silicates_tables)
"""

core_material = burnman.minerals.Murakami_2013.fe_perovskite()
mantle_material = burnman.minerals.SLB_2011.mg_fe_silicate_perovskite()
upper_mantle_material = burnman.minerals.SLB_2011.mg_fe_olivine(molar_fractions=(0.9, 0.1)) # when T < 2500 K and solid material in all mantle

layers_normal = [[core_material], [mantle_material, upper_mantle_material]]
layers_hot = [[core_material], [mantle_material]] # when T > 2500 K and liquid silicate in upper mantle

# 0th element = core molar fraction, 1st element = mantle molar fractions
minfractions_normal = [[1], [[0.9, 0.1, 0.0], [1]]]  # used when temperatures are less than 2500 K, es
minfractions_hot = [[1], [0.9, 0.1, 0.0]]  # molten upper mantle, liquid silicates added within eos.py


def varying_temp(type_eos, temperatures, central_pressures, grid_size, loc, testing):
    """
    Creates planetary grids with varying anchor or constant temperatures and saves them
    :param type_eos: adiabatic/constant
    :param temperatures: temperature range for grids
    :param central_pressures: central pressures to explore for each grid
    :param grid_size: size of grid, format [row, col], e.g., [4,4]
    :param loc: location where data for each grid be saved
    :return: all planetary grids (note these have already been saved individually in loc)
    """
    planetary_grids = []  # list to hole all planetary grid objects

    for temp in temperatures:
        # iterates over every temperature in range and initializes a planetary grid
        if temp > 2500:
            # using molten silicates for upper mantle
            temp_plan_grid = PlanetGrid(temp, central_pressures, grid_size, str(type_eos), loc, layers_types=layers_hot,
                                        minfractions=minfractions_hot, testing=testing)
        else:
            # using enstatite for upper mantle
            temp_plan_grid = PlanetGrid(temp, central_pressures, grid_size, str(type_eos), loc, layers_types=layers_normal,
                                        minfractions=minfractions_normal, testing=testing)
        # integrates grid
        temp_plan_grid.integrateGrid()
        # add grid to planetary_grids list
        planetary_grids.append(temp_plan_grid)
    return planetary_grids


if __name__ == "__main__":  # only executes if running run.py versus calling a function
    # CREATE GRIDS AND DO THERMAL EVOLUTION

    # done for most
    # default_central_pressures = [8, 12]
    # default_grid_size = [10, 10]
    # temp_range = [2100, 1500, 2000]

    #second step
    default_central_pressures = [8, 12]
    default_grid_size = [100, 100]
    temp_range = [300, 1000, 3000]

    # adiabatic planetary grids
    varying_temp("_adiabatic_", temp_range, default_central_pressures, default_grid_size, thermal_location, testing)
    # constant temperature planetary grids
    varying_temp("_constant_", temp_range, default_central_pressures, default_grid_size, thermal_location, testing)


