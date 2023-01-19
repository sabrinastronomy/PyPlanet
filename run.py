"""
This file wraps all other scripts such that you can run a full varying temperature rocky planet grid from this script.
Written by Sabrina Berger
"""

import sys
# importing packages
from mpi4py import MPI
from planet_grid_mpi import PlanetGrid  # where the planet grid is created and run
import numpy as np  # numpy
import burnman.minerals  # relevant burnman minerals package
# This testing variable is for debugging purposes
testing = False

# This is the location where ALL output planetary grids are stored
location = "/home/scberger/scratch-midway2/PyPlanet"
# location = "/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/Current Research/RockyPlanets/PyPlanet/paper"
# location = "/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/RockyPlanets/paper_data"
thermal_location = location + "/complete_data_with_silicate_mantle" #TO DO change back

# Grid of planets
# initializing MPS comm
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print(f"CURRENT RANK IS {rank}")

#### local directory
# location = "/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/Current Research/RockyPlanets/PyPlanet"
# thermal_location = location
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
enstatite: 'Mg2Si2O6' i.e., SLB_2011.enstatite
molar fractions are [0.9, 0.1]


If greater than 2500 K:
LOWER AND UPPER 
-----------------------------
upper mantle is a liquid silicate (see molten_silicates_tables)
"""

core_material = burnman.minerals.other.Fe_Dewaele()
mantle_material = burnman.minerals.SLB_2011.mg_fe_silicate_perovskite()
upper_mantle_material = burnman.minerals.SLB_2011.enstatite() # when no silicate mantle
#
layers_normal = [[core_material], [mantle_material, upper_mantle_material]]
# layers_hot = [[core_material], [mantle_material]] # when no silicate mantle

# 0th element = core molar fraction, 1st element = mantle molar fractions
minfractions_normal = [[1], [[1, 0, 0.0], [1]]]  # used when there isn't a silicate mantle
minfractions_hot = [[1], [1, 0, 0.0]]  # molten upper mantle, liquid silicates added within eos.py


def varying_temp(type_eos, entropy_range, central_pressures, grid_size, loc, testing, restart=False, last_index=None):
    """
    Creates planetary grids with varying anchor or constant temperatures and saves them
    :param type_eos: adiabatic/constant
    :param temperatures: temperature range for grids
    :param central_pressures: central pressures to explore for each grid
    :param grid_size: size of grid, format [row, col], e.g., [4,4]
    :param loc: location where data for each grid be saved
    :return: all planetary grids (note these have already been saved individually in loc)
    """
    planetary_grids = []  # list to hold all planetary grid objects

    for S in entropy_range:
        # iterates over every entropy in range and initializes a planetary grid

        # molten, _ = eos_new.silicate_melt_or_not(press, temp)
        #
        # if molten:
        #     # using molten silicates for upper mantle
        #     temp_plan_grid = PlanetGrid(temp, central_pressures, grid_size, str(type_eos), loc, layers_types=layers_hot,
        #                                 minfractions=minfractions_hot, testing=testing)
        # else:
            # using enstatite for upper mantle

        temp_plan_grid = PlanetGrid(S, central_pressures, grid_size, str(type_eos), loc, layers_types=layers_normal,
                                    minfractions=minfractions_normal, testing=testing, restart=restart,
                                    last_index=last_index, rank=rank)
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
    default_central_pressures = [11, 12]
    default_grid_size = [3, 3]
    # entropy = sys.argv[1] #TO DO change back
    entropy = 1000
    entropy_range = [int(entropy)] # new in 2022

    # adiabatic planetary grids
    varying_temp("_adiabatic_", entropy_range, default_central_pressures, default_grid_size, thermal_location, testing)
    # constant temperature planetary grids
    # varying_temp("_constant_", entropy_range, default_central_pressures, default_grid_size, thermal_location, testing)


