from scipy import interpolate
import numpy as np
import astropy.constants as const
import matplotlib.pyplot as plt
import scipy.constants as sp


mars_cmf = 0.26
earth_cmf = 0.33
cmfs_of_interest = [earth_cmf, mars_cmf]
# source_location = "/Users/sabrinaberger/RockyPlanets/DataFiles/"
# final_location = "/Users/sabrinaberger/RockyPlanets/MassRadiusDiagramData/"


location_src = "/Users/sabrinaberger/RockyPlanets"
data = location_src + "/thermalData/"
temps = np.linspace(300, 1000, 100)

def grid_to_1Darray(grid):
    compressed = []
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            compressed.append(grid[i][j])
    return compressed


def planet_interp(top_location, anchor_temps, val_of_interest, type = "adiabatic", label='u'):
    mass_temped = []
    radius_temped = []
    u_temped = []
    p_cmb_p_c_temped = []

    type = "_" + type + "_"
    for temp in anchor_temps:
        location = top_location + "{}{}data/".format(str(temp), type)
        p_c_list = grid_to_1Darray(np.load(location + "p_c_grid" + type + str(temp) + ".pyc" + ".npy"))
        p_cmb_pc_list = grid_to_1Darray(np.load(location + "p_cmb_percentage_grid" + type + str(temp) + ".pyc" + ".npy"))
        mass_list = grid_to_1Darray(np.load(location + "mass_grid" + type + str(temp) + ".pyc" + ".npy"))
        radius_list = grid_to_1Darray(np.load(location + "radius_grid" + type + str(temp) + ".pyc" + ".npy"))
        cmf_list = grid_to_1Darray(np.load(location + "core_mass_grid" + type + str(temp) + ".pyc" + ".npy"))
        cmr_list = grid_to_1Darray(np.load(location + "core_rad_grid" + type + str(temp) + ".pyc" + ".npy"))
        u_list = grid_to_1Darray(np.load(location + "u_grid" + type + str(temp) + ".pyc" + ".npy"))


        p_c_unique_dict = {}
        mass_values = [] # test plots
        radius_values = [] # test plots
        u_values = []

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
            desired_mass_list = []
            desired_radius_list = []
            desired_u_list = []
            desired_p_cmb_p_c_list = []
            desired_cmf_list = []
            desired_cmr_list = []
            for j in desired_indices:
                desired_p_cmb_p_c_list.append(p_cmb_pc_list[j])
                desired_mass_list.append(mass_list[j])
                desired_radius_list.append(radius_list[j])
                desired_cmf_list.append(cmf_list[j])
                desired_cmr_list.append(cmr_list[j])
                desired_u_list.append(u_list[j])

            p_c_unique_dict[p_c] = [desired_mass_list, desired_radius_list, desired_u_list, desired_p_cmb_p_c_list, desired_cmf_list, desired_cmr_list]

            if label is 'mr':
                # Mass & Radius
                p_cmb_p_c = interpolate.interp1d(desired_cmf_list, desired_p_cmb_p_c_list, bounds_error=False, fill_value="extrapolate")
                p_cmb_p_c_cmf0 = p_cmb_p_c(val_of_interest)

                mass = interpolate.interp1d(desired_p_cmb_p_c_list, desired_mass_list)
                radius = interpolate.interp1d(desired_p_cmb_p_c_list, desired_radius_list)
                cmr = interpolate.interp1d(desired_p_cmb_p_c_list, desired_cmr_list)
                u = interpolate.interp1d(desired_p_cmb_p_c_list, desired_u_list)

                mass_cmf0 = mass(p_cmb_p_c_cmf0)
                radius_cmf0 = radius(p_cmb_p_c_cmf0)
                # properties = [temp, cmf_of_interest, p_cmb_p_c_cmf0, cmr, mass_cmf0, radius_cmf0, cmr_cmf0]
                # np.save(final_location + str(temp) + "_" + str(cmf_of_interest) + ".pyc", properties)

                mass_values.append(mass_cmf0)
                radius_values.append(radius_cmf0)


                # Thermal Energy
                # p_cmb_p_c = interpolate.interp1d(desired_mass_list, desired_p_cmb_p_c_list, fill_value="extrapolate")
                # p_cmb_p_c_m0 = p_cmb_p_c(val_of_interest)
                #
                # radius = interpolate.interp1d(desired_p_cmb_p_c_list, desired_radius_list, fill_value="extrapolate")
                # u = interpolate.interp1d(desired_p_cmb_p_c_list, desired_u_list, fill_value="extrapolate")
                #
                # radius_m0 = radius(p_cmb_p_c_m0)
                # u_m0 = u(p_cmb_p_c_m0)
                #
                # radius_values.append(radius_m0)
                # u_values.append(u_m0)

        # mass_values = np.array(mass_values)/const.M_earth.value
        # radius_values = np.array(radius_values)/const.R_earth.value

        if label is "u":
            masses = []
            radii = []
            us = []
            p_cs = []
            for p_c in p_c_unique_dict.keys():
                arrs = p_c_unique_dict[p_c]
                desired_mass_list, desired_radius_list, desired_u_list, desired_p_cmb_p_c_list, desired_cmf_list, desired_cmr_list = arrs[0], arrs[1], arrs[2], arrs[3], arrs[4], arrs[5]
                masses.append(desired_mass_list)
                radii.append(desired_radius_list)
                us.append(desired_u_list)
                # p_cs.extend([p_c]*len(masses))

            masses = np.array(masses).flatten()
            radii = np.array(radii).flatten()
            us = np.array(us).flatten()
            # p_cs = np.array(p_cs).flatten()

            order = masses.argsort()
            masses = masses[order]
            radii = radii[order]
            us = us[order]
            # p_cs = p_cs[order]

            u_func = interpolate.interp1d(masses, us)
            r_func = interpolate.interp1d(masses, radii)

            u_calc = u_func(val_of_interest)
            r_calc = r_func(val_of_interest)

            u_temped.append(u_calc)
            radius_temped.append(r_calc)


    if label is not "u":
        mass_temped.append(mass_values)
        return np.array(u_temped), np.array(mass_temped), np.array(radius_temped), list(p_c_unique_dict.keys())

    return np.array(u_temped).flatten(), np.array(radius_temped).flatten(), list(p_c_unique_dict.keys())


def time(R_p, T_s, T_star, R_star, a, delU):
    # average values for all input params
    const = (-sp.sigma* 4 * sp.pi * R_p**2)*(T_s**4 -T_star**4*(R_star/(2*a))**2)
    return delU/const


