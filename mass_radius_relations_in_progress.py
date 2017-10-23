from scipy import interpolate
import numpy as np
import astropy.constants as const

# mars_cmf = 0.26
# earth_cmf = 0.33
# cmfs_of_interest = [earth_cmf, mars_cmf]
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
        p_c_list = grid_to_1Darray(np.load(location + "p_c" + label + str(temp) + ".pyc" + ".npy"))
        p_cmb_pc_list = grid_to_1Darray(np.load(location + "p_cmb_percentage" + label + str(temp) + ".pyc" + ".npy"))
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
            mass = interpolate.interp1d(desired_p_cmb_p_c_list, desired_mass_list)
            radius = interpolate.interp1d(desired_p_cmb_p_c_list, desired_radius_list)
            cmr = interpolate.interp1d(desired_p_cmb_p_c_list, desired_cmr_list)

            mass_cmf0 = mass(p_cmb_p_c_cmf0)
            radius_cmf0 = radius(p_cmb_p_c_cmf0)
            cmr_cmf0 = cmr(p_cmb_p_c_cmf0)
            # properties = [temp, cmf_of_interest, p_cmb_p_c_cmf0, cmr, mass_cmf0, radius_cmf0, cmr_cmf0]
            # np.save(final_location + str(temp) + "_" + str(cmf_of_interest) + ".pyc", properties)
            mass_values.append(mass_cmf0)
            radius_values.append(radius_cmf0)
        mass_values = np.array(mass_values)/const.M_earth.value
        radius_values = np.array(radius_values)/const.R_earth.value
        print("done")
        return np.sort(mass_values), np.sort(radius_values)

        #save interpolated values



