from scipy import interpolate
import numpy as np
import astropy.constants as const
import matplotlib.pyplot as plt
import scipy.constants as sp
from glob import glob
import os
import imageio



mars_cmf = 0.26
earth_cmf = 0.33
cmfs_of_interest = [earth_cmf, mars_cmf]
# source_location = "/Users/sabrinaberger/RockyPlanets/DataFiles/"
# final_location = "/Users/sabrinaberger/RockyPlanets/MassRadiusDiagramData/"
location = "/Users/sabrinaberger/RockyPlanets"

# class PlanetInterp:

#stackoverflow clode
def non_increasing(L):
    return all(x>=y for x, y in zip(L, L[1:]))

def non_decreasing(L):
    return all(x<=y for x, y in zip(L, L[1:]))

def monotonic(L):
    return non_increasing(L) or non_decreasing(L)

location_src = "/Users/sabrinaberger/RockyPlanets"
data = location_src + "/thermalData/"
# temps = np.linspace(300, 1000, 100)

def grid_to_1Darray(grid):
    compressed = []
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            compressed.append(grid[i][j])
    return compressed

def file_load(direc, name, temp):
    # Globs for file in direc with part of name, returns loaded first element of list
    # print(glob(direc + name + "*")[0])
    # print(glob(direc + name + "*")[0])
    return np.load(direc + name + "_adiabatic_" + temp + ".pyc.npy")


def makePropertyPlot(dir, propertyName, temp, units, plotTitle, filename, type = "adiabatic", colorscheme='magma_r', difference=False, top_location= location + "/thermalData/"):
    temp2 = int(temp)
    temp = str(temp)
    print(temp)
    type = "_" + type + "_"

    location = top_location + "{}{}data/".format(temp, type)
    p_c = file_load(location, "p_c_grid", temp)
    p_cmb_percentage = file_load(location, "p_cmb_percentage_grid", temp)
    propertyArray = file_load(location, propertyName, temp)

    if difference:
        property = plt.contourf(np.log10(p_c), p_cmb_percentage, propertyArray,
                                cmap=plt.get_cmap(colorscheme))
    else:
        property = plt.contourf(np.log10(p_c), p_cmb_percentage, np.log10(propertyArray),
                                cmap=plt.get_cmap(colorscheme))
    # property.set_xscale('log')
    # property.set_yscale('log')
    cbar = plt.colorbar(property)
    cbar.ax.set_ylabel("log(U)" + ' ' + units)

    # levs = np.power(10, lev_exp)

    plt.xlabel('log($P_c$) (Pa)')
    plt.ylabel('$P_{CMB}$/$P_c$')
    # plt.xlim(10, 12)
    plt.title(plotTitle + " (T = {})".format(int(temp2)))
    print(dir + filename)
    plt.savefig(dir + filename, dpi=600)
    plt.close()


def animate(name, Tarr, dir='/Users/sabrinaberger/PyPlanet/png/'):
    # for T in Tarr:
    #     makePropertyPlot(dir, "u_grid", T, "[J]", "Thermal Energy", "therm_contour_{}.png".format(T))
    images = []
    filelist = sorted(os.listdir(dir), key=lambda x: int(os.path.splitext(x)[0]))
    print(filelist)
    for file_name in filelist:
        if file_name.endswith('.png'):
            file_path = os.path.join(dir, file_name)
            images.append(imageio.imread(file_path))
    imageio.mimsave(dir + "{}_movie.gif".format(name), images)


def planet_interp(top_location, anchor_temps, cmfs_of_interest, masses_of_interest, type = "adiabatic", label='u'):
    # CMFs of interest
    # Masses of interest
    mark = False
    mass_temped = []
    radius_temped = []
    u_temped = []
    p_cmb_p_c_temped = []
    type = "_" + type + "_"
    for temp in anchor_temps:
        temp = str(temp)
        location = top_location + "{}{}data/".format(temp, type)
        p_c_list = grid_to_1Darray(file_load(location, "p_c_grid", temp))
        p_cmb_pc_list = grid_to_1Darray(file_load(location, "p_cmb_percentage_grid", temp))
        mass_list = grid_to_1Darray(file_load(location, "mass_grid", temp))
        radius_list = grid_to_1Darray(file_load(location, "radius_grid", temp))
        cmf_list = grid_to_1Darray(file_load(location, "core_mass_grid", temp))
        cmr_list = grid_to_1Darray(file_load(location, "core_rad_grid", temp))
        u_list = grid_to_1Darray(file_load(location, "u_grid", temp))
        if temp == 3000:
            u_grid = file_load(location, "u_grid", temp)
            r_grid = file_load(location, "radius_grid", temp)
            p_grid = file_load(location, "p_c_grid", temp)
            p_perc_grid = file_load(location, "p_cmb_percentage_grid", temp)

            plt.contourf(np.log10(p_grid), p_perc_grid, np.log10(u_grid))
            plt.xlabel('log($P_c$) (Pa)')
            plt.ylabel('$P_{CMB}$/$P_c$')
            plt.title("U")
            plt.savefig("u_example_contour.pdf")
            plt.close()

            plt.contourf(np.log10(r_grid), p_perc_grid, np.log10(u_grid))
            plt.xlabel('log($P_c$) (Pa)')
            plt.ylabel('$P_{CMB}$/$P_c$')
            plt.title("R")
            plt.savefig("r_example_contour.pdf")
            plt.close()

        p_c_unique_dict = {}
        mass_values = [] # test plots
        radius_values = [] # test plots
        u_values = []

        i = 0
        while i < p_c_list.__len__():

            p_c = p_c_list[i]
            if str(p_c) in p_c_unique_dict:
                # if p_c is already there, add the index to the values in the dict

                p_c_unique_dict.get(str(p_c)).append(i)
            else:
                # add new p_c and associated index

                p_c_unique_dict[str(p_c)] = [i]
            i += 1
        for p_c in p_c_unique_dict.keys():

            # separate parameters based on the same p_c
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
                p_cmb_p_c_cmf0 = p_cmb_p_c(cmfs_of_interest)

                mass = interpolate.interp1d(desired_p_cmb_p_c_list, desired_mass_list)
                radius = interpolate.interp1d(desired_p_cmb_p_c_list, desired_radius_list)
                cmr = interpolate.interp1d(desired_p_cmb_p_c_list, desired_cmr_list)
                u = interpolate.interp1d(desired_p_cmb_p_c_list, desired_u_list)

                mass_cmf0 = mass(p_cmb_p_c_cmf0)
                radius_cmf0 = radius(p_cmb_p_c_cmf0)
                # properties = [temp, cmf_of_interest, p_cmb_p_c_cmf0, cmr, mass_cmf0, radius_cmf0, cmr_cmf0]
                # np.save(final_location + temp + "_" + str(cmf_of_interest) + ".pyc", properties)

                mass_values.append(mass_cmf0)
                radius_values.append(radius_cmf0)

        # put parameters in terms of Earth units
        # mass_values = np.array(mass_values)/const.M_earth.value
        # radius_values = np.array(radius_values)/const.R_earth.value

        if label is "u":

            masses = []
            radii = []
            us = []
            cmfs = []
            for p_c in p_c_unique_dict.keys():
                arrs = p_c_unique_dict[p_c]
                desired_mass_list, desired_radius_list, desired_u_list, desired_p_cmb_p_c_list, desired_cmf_list, desired_cmr_list = arrs[0], arrs[1], arrs[2], arrs[3], arrs[4], arrs[5]
                masses.append(desired_mass_list)
                radii.append(desired_radius_list)
                us.append(desired_u_list)
                cmfs.append(desired_cmf_list)

            masses = np.array(masses).flatten()
            radii = np.array(radii).flatten()
            us = np.array(us).flatten()
            cmfs = np.array(cmfs).flatten()
            # order = masses.argsort()
            order = radii.argsort()
            masses = masses[order]
            radii = radii[order]
            us = us[order]
            cmfs = cmfs[order]



            mass_cmfs = interpolate.interp1d(cmfs, masses, bounds_error=False,
                                             fill_value="extrapolate")


            u_func = interpolate.interp1d(masses, us, fill_value="extrapolate")
            r_func = interpolate.interp1d(masses, radii, fill_value="extrapolate")

            u_calc = u_func(masses_of_interest)
            r_calc = r_func(masses_of_interest)

            # Checking monotonicity
            # u_calc = u_func(mass_0)
            # r_calc = r_func(mass_0)

            # if not monotonic(masses):
            #     print("mass is not monotonic")
            #
            # if not monotonic(us):
            #     plt.close()
            #     plt.plot(us)
            #     plt.plot(u_func(masses))
            #     plt.show()
            #     print(us)
            #     print("u is not monotonic")
            # if len(u_temped) is not 0:
            #     if u_temped[-1] < u_calc:
            #         mark = True
            #         print("oops")
            #         print(temp)
            # elif mark is True:
            #     print("temp {}",format(temp))


            u_temped.append(u_calc)
            radius_temped.append(r_calc)

    if label is not "u":
        mass_temped.append(mass_values)
        return np.array(u_temped), np.array(mass_temped), np.array(radius_temped), list(p_c_unique_dict.keys())

    print(temp)
    print("U_monotonicity {}".format(monotonic(np.array(u_temped).flatten())))
    print("R_monotonicity {}".format(monotonic(np.array(radius_temped).flatten())))

    return np.array(u_temped).flatten(), np.array(radius_temped).flatten(), list(p_c_unique_dict.keys())


def time(R_p, T_s, T_star, R_star, a, delU):
    # average values for all input params
    const = (-sp.sigma* 4 * sp.pi * R_p**2)*(T_s**4 -T_star**4*(R_star/(2*a))**2)
    return delU/const


