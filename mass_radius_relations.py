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
        label = "_" + label + "_"
        p_c_list = grid_to_1Darray(np.load(location + "p_c_grid" + label + str(temp) + ".pyc" + ".npy"))
        print("length {}".format(len(p_c_list)))
        print("max {}".format(max(p_c_list)))
        print("min {}".format(min(p_c_list)))

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

    # temp_range_test = [3000]
    # colors = ["blue", "orange"]
    # styles_adiabatic = ["solid", "dashed", "dashdot"]
    # styles_constant = ["dashed", "solid"]
    # style_w_temps = [styles_adiabatic, styles_constant]
    # temp_types = ["adiabatic", "constant"]
    # cmfs = [0.33, 0.67]
    # for cmf, color in zip(cmfs, colors):
    #     for temp_type, styles in zip(temp_types, style_w_temps):
    #         for temp, style in zip(temp_range_test, styles):
    #             # if temp == 1000 and temp_type == "constant":
    #             #     continue
    #             data_files_stored_in = f"/Users/sabrinaberger/complete_data_with_silicate_mantle/{temp}_{temp_type}_data/"
    #             mass_plots, radius_plots = planet_interp(data_files_stored_in, temp_type, [temp], cmf)
    #             plt.plot(mass_plots, radius_plots, c=color, ls=style)
    #             if cmf == cmfs[0]:
    #                 plt.plot(np.nan, np.nan, c="k", ls=style, label=f"{temp} K {temp_type.title()}")
    # # create legend
    # plt.plot(np.nan, np.nan, c="blue", ls="-", label="Earth: CMF = 0.33 ")
    # plt.plot(np.nan, np.nan, c="orange", ls="-", label="Mercury: CMF = 0.66 ")
    # # plt.xlim(1, 4.5)
    #
    # plt.legend()
    # # plt.ylabel(r'$\log(R_{pl})$', r'[$R_{\bigoplus}$]')
    # # plt.xlabel(r'$\log(M_{pl})$', r'[$M_{\bigoplus}$]')
    # plt.title("Linear Mass-Radius Relationships")
    #     # plt.plot(mass_plots, 100*(radius_plots-radius_plots_300K)/radius_plots_300K)
    # # print("plotted")
    # plt.savefig("paper_plots/linear_new_mr.pdf")
    # plt.savefig("paper_plots/linear_new_mr.png", dpi=300)
    # plt.close()

    ### DIFFERENCE
    temp_range_test = [300, 1000, 3000]
    colors = ["blue", "orange"]
    styles_adiabatic = ["dashdot", "dashed", "dotted"]
    styles_constant = ["dashed", "solid"]
    style_w_temps = [styles_adiabatic, styles_constant]
    temp_types = ["adiabatic"]
    cmfs = [0.33, 0.67]

    for cmf, color in zip(cmfs, colors):
        for temp_type, styles in zip(temp_types, style_w_temps):
            for temp, style in zip(temp_range_test, styles):
                compare_mass_plots, compare_radius_plots = planet_interp("/Users/sabrinaberger/complete_data_with_silicate_mantle/300_constant_data/", "constant", [300], cmf)
                data_files_stored_in = f"/Users/sabrinaberger/complete_data_with_silicate_mantle/{temp}_{temp_type}_data/"
                mass_plots, radius_plots = planet_interp(data_files_stored_in, temp_type, [temp], cmf)
                print(min(compare_mass_plots))

                print(min(mass_plots))

                perc_mass = 100*(mass_plots - compare_mass_plots)/compare_mass_plots
                perc_rad = 100*(radius_plots - compare_radius_plots)/compare_radius_plots
                plt.plot(mass_plots, perc_rad, c=color, ls=style)
                if cmf == cmfs[0]:
                    plt.plot(np.nan, np.nan, c="k", ls=style, label=f"{temp} K {temp_type.title()}")
    # create legend
    plt.plot(np.nan, np.nan, c="blue", ls="-", label="Earth: CMF = 0.33 ")
    plt.plot(np.nan, np.nan, c="orange", ls="-", label="Mercury: CMF = 0.66 ")
    # plt.xlim(1, 4.5)
    plt.ylabel("Radial % Difference from 300K Constant Planets")
    plt.legend()
    # plt.ylabel(r'$\log(R_{pl})$', r'[$R_{\bigoplus}$]')
    # plt.xlabel(r'$\log(M_{pl})$', r'[$M_{\bigoplus}$]')
    plt.title("Comparing Mass-Radius Relationships")
        # plt.plot(mass_plots, 100*(radius_plots-radius_plots_300K)/radius_plots_300K)
    # print("plotted")
    plt.savefig("paper_plots/diff_new_mr.pdf")
    plt.savefig("paper_plots/diff_new_mr.png", dpi=300)