"""
This file contains the interpolation method to go from full planetary grids to mass-radius relations.
Written by Sabrina Berger
"""

from scipy import interpolate
import numpy as np
import astropy.constants as const

mars_cmf = 0.26
earth_cmf = 0.33
cmfs_of_interest = [earth_cmf, mars_cmf]
# source_location = "/Users/sabrinaberger/RockyPlanets/DataFiles/"
# final_location = "/Users/sabrinaberger/RockyPlanets/MassRadiusDiagramData/"


def grid_to_1Darray(grid):
    compressed = []
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            compressed.append(grid[i][j])
    return compressed


def planet_interp(location, label, anchor_temps, cmf_of_interest):
    mass_plots = []
    radius_plots = []
    for temp in anchor_temps:
        properties = []
        p_c_list = grid_to_1Darray(np.load(location + "p_c_grid" + label + str(temp) + ".pyc" + ".npy"))
        p_cmb_pc_list = grid_to_1Darray(np.load(location + "p_cmb_percentage_grid" + label + str(temp) + ".pyc" + ".npy"))
        mass_list = grid_to_1Darray(np.load(location + "mass_grid" + label + str(temp) + ".pyc" + ".npy"))
        radius_list = grid_to_1Darray(np.load(location + "radius_grid" + label + str(temp) + ".pyc" + ".npy"))
        cmf_list = grid_to_1Darray(np.load(location + "core_mass_grid" + label + str(temp) + ".pyc" + ".npy"))
        cmr_list = grid_to_1Darray(np.load(location + "core_rad_grid" + label + str(temp) + ".pyc" + ".npy"))
        p_c_unique_dict = {}
        mass_values = [] # test plots
        radius_values = [] # test plots
        i = 0
        while i < p_c_list.__len__():
            p_c = p_c_list[i]
            if str(p_c) in p_c_unique_dict:
                p_c_unique_dict.get(str(p_c)).append(i)
            else:
                p_c_unique_dict[str(p_c)] = [i]
            i += 1
        for p_c in p_c_unique_dict.keys():
            desired_indices = p_c_unique_dict.get(str(p_c))
            desired_p_cmb_p_c_list = []
            desired_mass_list = []
            desired_radius_list = []
            desired_cmf_list = []
            desired_cmr_list = []
            for j in desired_indices:
                desired_p_cmb_p_c_list.append(p_cmb_pc_list[j])
                desired_mass_list.append(mass_list[j])
                desired_radius_list.append(radius_list[j])
                desired_cmf_list.append(cmf_list[j])
                desired_cmr_list.append(cmr_list[j])

            p_cmb_p_c = interpolate.interp1d(desired_cmf_list, desired_p_cmb_p_c_list, bounds_error=False, fill_value="extrapolate")
            p_cmb_p_c_cmf0 = p_cmb_p_c(cmf_of_interest)
            mass = interpolate.interp1d(desired_p_cmb_p_c_list, desired_mass_list, bounds_error=False, fill_value="extrapolate")
            radius = interpolate.interp1d(desired_p_cmb_p_c_list, desired_radius_list, bounds_error=False, fill_value="extrapolate")
            cmr = interpolate.interp1d(desired_p_cmb_p_c_list, desired_cmr_list, bounds_error=False, fill_value="extrapolate")

            mass_cmf0 = mass(p_cmb_p_c_cmf0)
            radius_cmf0 = radius(p_cmb_p_c_cmf0)
            cmr_cmf0 = cmr(p_cmb_p_c_cmf0)
            # properties = [temp, cmf_of_interest, p_cmb_p_c_cmf0, cmr, mass_cmf0, radius_cmf0, cmr_cmf0]
            # np.save(final_location + str(temp) + "_" + str(cmf_of_interest) + ".pyc", properties)
            mass_values.append(mass_cmf0)
            radius_values.append(radius_cmf0)
        mass_values = np.array(mass_values)/const.M_earth.value
        radius_values = np.array(radius_values)/const.R_earth.value
        print(f"done with temp: {temp}")
        return np.sort(mass_values), np.sort(radius_values)

        #save interpolated values

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    # temp_range = np.linspace(300, 3000, 10)
    data_files_stored_in = "/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/RockyPlanets/PyPlanet/paper_data/complete_data_with_silicate_mantle/"
    temp_type = "_constant_"
    # getting constant planet
    for i, temp in enumerate([300]):
        mass_plots_300K, radius_plots_300K = planet_interp(data_files_stored_in + "{}{}data".format(temp, temp_type) + "/", temp_type, [temp], 0.33)

    # temp_range = np.linspace(300, 3000, 10)[:8]
    temp_range_test = [2400.0]
    # mass_plots = np.zeros(len(temp_range))
    # radius_plots = np.zeros(len(temp_range))
    temp_type = "_adiabatic_"
    for i, temp in enumerate(temp_range_test):
        mass_plots, radius_plots = planet_interp(data_files_stored_in + "{}{}data".format(temp, temp_type) + "/", temp_type, [temp], 0.33)
        plt.plot(mass_plots, 100*(radius_plots-radius_plots_300K)/radius_plots_300K)
    print("plotted")
    plt.savefig("paper_plots/new_mr.pdf")
