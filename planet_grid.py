import time
from planet import Planet
from eos import *
import subprocess

# Grid of planets

class PlanetGrid:
    def __init__(self, anchor_temp, central_pressures, grid_size, temp_profile, location, want_full_profile=False,
                 relative_tolerance=1e-5,
                 mantle_material=minerals.SLB_2011.mg_fe_silicate_perovskite(),
                 core_material=minerals.Murakami_2013.fe_perovskite()):
        # Initial parameters
        self.location = location

        self.mantle_material = mantle_material
        self.core_material = core_material
        self.temp_profile = temp_profile
        self.anchor_temp = anchor_temp
        self.central_pressures = central_pressures
        self.grid_size = grid_size  # grid_size[0] = number of central pressures, grid_size[1] = number of ratios of CMB pressure/central pressure
        self.want_full_profile = want_full_profile
        self.relative_tolerance = relative_tolerance
        self.save_folder = location + "/" + str(anchor_temp) + temp_profile + "data/"

        # p_cmb = core mantle boundary pressure
        # p_c = central pressure

        self.p_c_list = np.array(
            [10 ** x for x in linspace(self.central_pressures[0], self.central_pressures[1], self.grid_size[0])])
        self.p_cmb_p_c = linspace(0, 1, self.grid_size[1])  # transition pressure when multiplied by p_c

        self.xx, self.yy = meshgrid(self.p_c_list, self.p_cmb_p_c)

        self.num_rows = self.grid_size[0]
        self.num_cols = self.grid_size[1]
        self.t0 = float(1.0)
        self.dt = 200

        # aliases
        num_rows = self.num_rows
        num_cols = self.num_cols

        # To be calculated parameters
        self.radius_grid = ndarray(shape=(num_rows, num_cols))  # the total radii
        self.mass_grid = ndarray(shape=(num_rows, num_cols))  # the total masses
        self.press_grid = ndarray(shape=(num_rows, num_cols))  # the total pressures
        self.core_mass_grid = ndarray(shape=(num_rows, num_cols))  # masses at CMB
        self.core_rad_grid = ndarray(shape=(num_rows, num_cols))  # radii at CMB
        self.p_cmb_simulated = ndarray(shape=(num_rows, num_cols))  # transition pressures computed during integration
        self.p_cmb_grid = ndarray(shape=(num_rows, num_cols))  # transition pressures inputted

    def integrateGrid(self):
        # Time and Count Parameters
        planet_number = 1
        initial_time = time.clock()

        # Aliases
        location = self.location
        radius_grid = self.radius_grid
        mass_grid = self.mass_grid
        press_grid = self.press_grid
        core_mass_grid = self.core_mass_grid
        core_rad_grid = self.core_rad_grid
        p_cmb_simulated = self.p_cmb_simulated
        p_cmb_grid = self.p_cmb_grid
        anchor_temp = self.anchor_temp
        want_full_profile = self.want_full_profile
        xx = self.xx
        yy = self.yy
        save_folder = self.save_folder
        temp_profile = self.temp_profile

        print("Type of EoS used in this grid: " + self.temp_profile)

        # iterating over mesh grid of p_cmb and p_cmb/p_c ratio
        for i in range(self.num_rows):
            for j in range(self.num_cols):
                print(
                    "Location in mesh grid: %d %d \n P_c = %e \n P_cmb/P_c =% g" % (i, j, self.xx[i][j], self.yy[i][j]))
                print("Planet: " + str(planet_number) + " out of " + str(self.num_rows * self.num_cols))

                p_c = self.xx[i][j]
                transition_pressure = self.yy[i][j] * p_c

                needed_pressures = [p_c, transition_pressure, 0.]  # TODO extend to arbitrary number of layers
                print(needed_pressures)
                print("here")
                mantle = Layer("mantle", self.mantle_material, [0.1, 0.1, 0.8], self.temp_profile)
                core = Layer("core", self.core_material, [1], self.temp_profile)

                planet_eos = EoS(p_c, transition_pressure, self.temp_profile, core, mantle,
                                 self.anchor_temp)  # TODO MAKE FOR ARBITRARY NUMBER OF LAYERS
                planet = Planet([planet_eos.core_eos, planet_eos.mantle_eos], self.t0, self.dt, self.relative_tolerance,
                                needed_pressures)
                planet.integratePlanet()

                radius_grid[i][j] = planet.rad[-1]
                mass_grid[i][j] = planet.mass[-1]
                press_grid[i][j] = planet.press[-1]

                core_mass_grid[i][j] = planet.transition_mass_list[1] / planet.mass[-1]
                core_rad_grid[i][j] = planet.transition_rad_list[1] / planet.rad[-1]

                p_cmb_simulated[i][j] = planet.transition_press_list[1]
                p_cmb_grid[i][j] = transition_pressure

                planet_number += 1

        subprocess.call(["mkdir", save_folder])
        save(save_folder + "p_c" + temp_profile + str(anchor_temp) + ".pyc", xx)
        save(save_folder + "p_cmb_percentage" + temp_profile + str(anchor_temp) + ".pyc", yy)
        save(save_folder + "radius_grid" + temp_profile + str(anchor_temp) + ".pyc", radius_grid)
        save(save_folder + "mass_grid" + temp_profile + str(anchor_temp) + ".pyc", mass_grid)
        save(save_folder + "press_grid" + temp_profile + str(anchor_temp) + ".pyc", press_grid)
        save(save_folder + "core_mass_grid" + temp_profile + str(anchor_temp) + ".pyc", core_mass_grid)
        save(save_folder + "core_rad_grid" + temp_profile + str(anchor_temp) + ".pyc", core_rad_grid)
        save(save_folder + "p_cmb_simulated" + temp_profile + str(anchor_temp) + ".pyc", p_cmb_simulated)
        save(save_folder + "p_cmb_grid" + temp_profile + str(anchor_temp) + ".pyc", p_cmb_grid)

        print("The following data files have been created: ")
        print(save_folder + "p_c" + temp_profile + str(anchor_temp) + ".pyc")
        print(save_folder + "p_cmb_percentage" + temp_profile + str(anchor_temp) + ".pyc")
        print(save_folder + "radius_grid" + temp_profile + str(anchor_temp) + ".pyc")
        print(save_folder + "mass_grid" + temp_profile + str(anchor_temp) + ".pyc")
        print(save_folder + "press_grid" + temp_profile + str(anchor_temp) + ".pyc")
        print(save_folder + "core_mass_grid" + temp_profile + str(anchor_temp) + ".pyc")
        print(save_folder + "core_rad_grid" + temp_profile + str(anchor_temp) + ".pyc")
        print(save_folder + "p_cmb_simulated" + temp_profile + str(anchor_temp) + ".pyc")
        print(save_folder + "p_cmb_grid" + temp_profile + str(anchor_temp) + ".pyc")

        if want_full_profile:
            return xx, yy, radius_grid, mass_grid, press_grid, core_mass_grid, core_rad_grid

        print("Done. Integrated " + repr(planet_number - 1) + " " + temp_profile[1:-1] + " planets successfully!" + " This took " +
              repr((time.clock() - initial_time) / 60) + " minutes" + ".")
