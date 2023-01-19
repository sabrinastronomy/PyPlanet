"""
This file contains the class PlanetGrid that initializes a planetary grid of varying temperatures
and central pressure/CMB pressure ratios. The planetary grid can be integrated from here.
Written by Sabrina Berger
"""

# importing packages
import subprocess  # for creating directories
import time  # for timing integrations
import eos_newest
from eos_newest import *  # where eoses
from planet import Planet
from numpy import *
import os

# Grid of planets


class PlanetGrid:
    def __init__(self, entropy, central_pressures, grid_size, temp_profile, location, layers_types, minfractions,
                 relative_tolerance=1e-5, testing=False, debug_certain_planet=False, certain_planet=None, restart=False, last_index=None):
        # Initial parameters
        self.location = location
        self.temp_profile = temp_profile
        self.entropy = entropy
        self.central_pressures = central_pressures
        self.grid_size = grid_size  # grid_size[0] = number of central pressures, grid_size[1] = number of ratios of CMB pressure/central pressure
        self.relative_tolerance = relative_tolerance
        save_folder = self.save_folder = location + "/" + str(entropy) + temp_profile + "data/"
        self.debug_certain_planet = debug_certain_planet # boolean if you want to debug a specific planet in the grid
        self.certain_planet = certain_planet # row, col of certain planet that you want to debut
        self.restart = restart
        self.last_index = last_index
        try:
            os.makedirs(self.save_folder)
            print("Directory ", self.save_folder, " Created ")
        except FileExistsError:
            print("Directory ", self.save_folder, " already exists")

            # p_cmb = core mantle boundary pressure
        # p_c = central pressure

        # generating central pressures and ratio of core mantle boundary to central pressures
        self.p_c_list = np.array(
            [10 ** x for x in linspace(self.central_pressures[0], self.central_pressures[1], self.grid_size[0])])
        self.p_cmb_p_c = linspace(0, 1, self.grid_size[1])  # transition pressure when multiplied by p_c
        self.xx, self.yy = meshgrid(self.p_c_list, self.p_cmb_p_c)
        self.num_rows = self.grid_size[0]
        self.num_cols = self.grid_size[1]
        self.t0 = 1.0 # inital radius!
        self.dt = 200
        self.layers_types = layers_types
        self.minfractions = minfractions
        self.testing = testing

        # Aliases
        num_rows = self.num_rows
        num_cols = self.num_cols

        # To be calculated parameters
        if restart:
            self.radius_grid = load(save_folder + "radius_grid" + temp_profile + str(entropy) + ".npy")
            self.mass_grid = load(save_folder + "mass_grid" + temp_profile + str(entropy) + ".npy")
            self.press_grid = load(save_folder + "surf_press_grid" + temp_profile + str(entropy) + ".npy")
            self.core_mass_grid = load(save_folder + "core_mass_grid" + temp_profile + str(entropy) + ".npy")
            self.core_rad_grid = load(save_folder + "core_rad_grid" + temp_profile + str(entropy) + ".npy")
            self.p_cmb_simulated = load(save_folder + "p_cmb_simulated" + temp_profile + str(entropy) + ".npy")
            self.p_cmb_grid=load(save_folder + "p_cmb_grid" + temp_profile + str(entropy) + ".npy")
            self.u_grid = load(save_folder + "u_grid" + temp_profile + str(entropy) + ".npy")
            self.press_mat_actual = load(save_folder + "materials" + temp_profile + str(entropy) + ".npy", allow_pickle=True)
            self.full_radius_grid = load(save_folder + "full_rad_grids" + temp_profile + str(entropy) + ".npy")  # full radius grids

        else:
            self.radius_grid = full(shape=(num_rows, num_cols), fill_value=-1., dtype=float64)  # the total radii
            self.mass_grid = full(shape=(num_rows, num_cols), fill_value=-1, dtype=float64)  # the total masses
            self.press_grid = full(shape=(num_rows, num_cols), fill_value=-1., dtype=float64)  # the total pressures
            self.core_mass_grid = full(shape=(num_rows, num_cols), fill_value=-1., dtype=float64)  # masses at CMB
            self.core_rad_grid = full(shape=(num_rows, num_cols), fill_value=-1., dtype=float64)  # radii at CMB
            self.p_cmb_simulated = full(shape=(num_rows, num_cols), fill_value=-1., dtype=float64)  # transition pressures computed during integration
            self.p_cmb_grid = full(shape=(num_rows, num_cols), fill_value=-1., dtype=float64)  # transition pressures inputted


            # Thermal evolution
            self.u_grid = full(shape=(num_rows, num_cols), fill_value=-1., dtype=float64)  # relative total thermal energy

            # This is getting the materials in the planet at each pressure
            self.full_radius_grid = ndarray(shape=(num_rows, num_cols), dtype=object)  # full radius grids
            self.press_mat_actual = ndarray(shape=(num_rows, num_cols), dtype=object)


    def integrateGrid(self):
        # Time and Count Parameters
        planet_number = 0
        initial_time = time.time()

        # Aliases
        radius_grid = self.radius_grid
        mass_grid = self.mass_grid
        press_grid = self.press_grid
        core_mass_grid = self.core_mass_grid
        core_rad_grid = self.core_rad_grid
        p_cmb_simulated = self.p_cmb_simulated
        p_cmb_grid = self.p_cmb_grid
        entropy = self.entropy
        full_radius_grid = self.full_radius_grid
        press_mat_actual = self.press_mat_actual

        xx = self.xx
        yy = self.yy
        temp_profile = self.temp_profile # whether a constant or adiabatic planet is being integrated
        u_grid = self.u_grid
        save_folder = self.save_folder
        layers_types = self.layers_types
        minfractions = self.minfractions

        if self.restart:
            mask = self.radius_grid == -1  # creating mask where True if the element hasn't been explored yet
            if self.last_index == None:
                self.last_index = np.unravel_index(np.argmax(mask, axis=None), mask.shape) # should get the first index that's -1

            print("STARTING GRID FROM: ")
            print(self.last_index)
            iterator_rows = range(self.last_index[0], self.num_rows)
            iterator_cols = range(self.last_index[1], self.num_cols)
        else:
            if self.last_index != None:
                print("STARTING GRID FROM: ")
                print(self.last_index)
                iterator_rows = range(self.last_index[0], self.num_rows)
                iterator_cols = range(self.last_index[1], self.num_cols)
            else:
                iterator_rows = range(self.num_rows)
                iterator_cols = range(self.num_cols)

        print("Type of EoS used in this grid: " + self.temp_profile)

        # iterating over mesh grid of p_cmb and p_cmb/p_c ratio where each distinct i and j pair correspond to a different planet
        for i in iterator_rows:
            for j in iterator_cols:
                print("i, j")
                print(i)
                print(j)
                planet_number += 1
                self.layers = [] # instances of EoS layer
                layers = self.layers
                self.intermediate_transition_pressures = []
                intermediate_transition_pressures = self.intermediate_transition_pressures
                print("----------------------------------------------------------------------------------------------")
                print("Location in mesh grid: %d %d \n P_c = %e \n P_cmb/P_c =% g" % (i, j, self.xx[i][j], self.yy[i][j]))
                print("Planet: " + str(planet_number))

                # Instantiating all pressures
                p_c = xx[i][j]
                fraction = yy[i][j]
                intermediate_transition_pressures.append(p_c) # central pressure is the first transition pressure
                intermediate_transition_pressures.append(fraction * p_c) # CMB pressure is the 2nd transition pressure
                p_cmb = intermediate_transition_pressures[1]

                # if p_cmb == 0: # skipping planets with no mantles TODO add this to the paper
                #     continue

                for k, (layer, minfrac) in enumerate(zip(layers_types, minfractions)):
                    # TODO take this out
                    # layer or layers_types, contains materials that will be used in EoS.
                    # These are not used for adiabatic EoSes
                    layers.append(BurnmanLayer("layer{}".format(k), layer, minfrac, self.temp_profile))

                # TODO MESA CORES
                intermediate_transition_pressures.append(0.) # surface pressure is always the final transition pressure stopping integration
                print("p_c = {:e}, p_cmb = {:e}, p_surface = {:e}".format(intermediate_transition_pressures[0], intermediate_transition_pressures[1], intermediate_transition_pressures[2]))
                # print("layers " + str(layers))
                planet_eos = EoS(p_c, p_cmb, self.temp_profile,
                                 layers, self.entropy)
                planet = Planet(planet_eos, self.t0, self.dt, self.relative_tolerance,
                                intermediate_transition_pressures)
                if not self.testing:
                    planet.integratePlanet()
                    # Taking last value from each integrated property array which is the estimated value
                    # of the given property for the planet just integrated
                    radius_grid[i][j] = planet.rad[-1]
                    mass_grid[i][j] = planet.mass[-1]
                    press_grid[i][j] = planet.press[-1]

                    core_mass_grid[i][j] = planet.transition_mass_list[1] / planet.mass[-1]
                    core_rad_grid[i][j] = planet.transition_rad_list[1] / planet.rad[-1]

                    p_cmb_simulated[i][j] = planet.transition_press_list[1]
                    p_cmb_grid[i][j] = intermediate_transition_pressures[1] # getting the core mantle boundary pressures

                    # Thermal evolution
                    u_grid[i][j] = planet.u[-1]
                    full_radius_grid[i][j] = planet.rad

                    ### Get amount of material in each layer, e.g., to get magma ocean percentage
                    materials = np.full(np.shape(planet.rad), -1, dtype=int)
                    # iterate over all pressures actually found in planet
                    # for k, actual_press in enumerate(planet.press):
                    #     materials[k] = eos_newest.helper_func_get_closest(actual_press, planet_eos.pressures_concatenate, planet_eos.materials_concatenate)
                    # self.press_mat_actual[i][j] = materials

                    # saves after each planet that's integrated, this means that if your job runs out of time, your grid might be incomplete
                    save(save_folder + "p_c_grid" + temp_profile + str(entropy), xx)
                    save(save_folder + "p_cmb_percentage_grid" + temp_profile + str(entropy), yy)
                    save(save_folder + "radius_grid" + temp_profile + str(entropy), radius_grid)
                    save(save_folder + "mass_grid" + temp_profile + str(entropy), mass_grid)
                    save(save_folder + "surf_press_grid" + temp_profile + str(entropy), press_grid)
                    save(save_folder + "core_mass_grid" + temp_profile + str(entropy), core_mass_grid)
                    save(save_folder + "core_rad_grid" + temp_profile + str(entropy), core_rad_grid)
                    save(save_folder + "p_cmb_simulated" + temp_profile + str(entropy), p_cmb_simulated)
                    save(save_folder + "p_cmb_grid" + temp_profile + str(entropy), p_cmb_grid)
                    save(save_folder + "u_grid" + temp_profile + str(entropy), u_grid)
                    save(save_folder + "materials" + temp_profile + str(entropy), press_mat_actual)
                    save(save_folder + "full_radius_grid" + temp_profile + str(entropy), full_radius_grid)

                    print("The following data files have been created or overwritten: ")
                    print(save_folder + "p_c_grid" + temp_profile + str(entropy))
                    print(save_folder + "p_cmb_percentage_grid" + temp_profile + str(entropy))
                    print(save_folder + "radius_grid" + temp_profile + str(entropy))
                    print(save_folder + "mass_grid" + temp_profile + str(entropy))
                    print(save_folder + "surf_press_grid" + temp_profile + str(entropy))
                    print(save_folder + "core_mass_grid" + temp_profile + str(entropy))
                    print(save_folder + "core_rad_grid" + temp_profile + str(entropy))
                    print(save_folder + "p_cmb_simulated" + temp_profile + str(entropy))
                    print(save_folder + "p_cmb_grid" + temp_profile + str(entropy))
                    print(save_folder + "u_grid" + temp_profile + str(entropy))
                    print(save_folder + "materials" + temp_profile + str(entropy))
                    print(save_folder + "full_radius_grid" + temp_profile + str(entropy))
                    print("Saved materials and full radius profile...")

        print("Done. Integrated " + str(planet_number) + " " + temp_profile[1:-1] + " planets successfully!" + " This took " +
              str((time.time() - initial_time) / 60) + " minutes" + ".")
