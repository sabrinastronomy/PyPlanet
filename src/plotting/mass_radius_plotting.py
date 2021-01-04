import matplotlib.pyplot as plt
import csv
import numpy as np
# from mass_radius_interpolation import planet_interp

data_files_stored_in = "/Users/sabrinaberger/Desktop/DataFiles/"

temp = np.load("/Users/sabrinaberger/Research with Leslie Rogers - 2016/Code/DataFiles/prefix_adiabat_.pyc.npy")
radius = np.load("/Users/sabrinaberger/Research with Leslie Rogers - 2016/Code/DataFiles/radius_adiabat_.pyc.npy")
mass = np.load("/Users/sabrinaberger/Research with Leslie Rogers - 2016/Code/DataFiles/mass_adiabat_.pyc.npy")

mass = log10(mass)
radius = log10(radius)

plt.plot(temp, mass)

plt.xlabel('Surface Temperature (K)')
plt.ylabel('log($mass$) (kg)')
plt.savefig("adiabat_mass.pdf")
plt.close()

plt.plot(temp, radius)
plt.xlabel('Surface Temperature (K)')
plt.ylabel('log($radius$) (m)')
plt.savefig("adiabat_radius.pdf")
plt.close()

property = plt.contourf(temp, radius, mass, p)
cbar = plt.colorbar(property)
cbar.ax.set_ylabel('log($mass$) (kg)')
plt.xlabel('Surface Temperature (K)')
plt.ylabel('log($radius$) (m)')
plt.title("Final Mass and Radius, p_c = 1.7e11 Pa, p_cmb = 136e9 Pa")
plt.savefig("/Users/sabrinaberger/Research with Leslie Rogers - 2016/Code/Plots_Burnman/adiabat_total.py")
plt.show()
# plt.close()


#
# temp_range = [300, 3000]
#
# mass_plots = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# radius_plots = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#
# mass_plots[0], radius_plots[0] = planet_interp(data_files_stored_in, "_adiabatic_", [300], 0.33)
# mass_plots[1], radius_plots[1] = planet_interp(data_files_stored_in, "_constant_", [300], 0.33)
#
# # mass_plots[2], radius_plots[2] = planet_interp(data_files_stored_in, "_adiabatic_", [1000], 0.33)
# # mass_plots[3], radius_plots[3] = planet_interp(data_files_stored_in, "_constant_", [1000], 0.33)
#
# mass_plots[4], radius_plots[4] = planet_interp(data_files_stored_in, "_adiabatic_", [3000], 0.33)
# mass_plots[5], radius_plots[5] = planet_interp(data_files_stored_in, "_constant_", [3000], 0.33)
#
# mass_plots[6], radius_plots[6] = planet_interp(data_files_stored_in, "_adiabatic_", [300], 0.69)
# mass_plots[7], radius_plots[7] = planet_interp(data_files_stored_in, "_constant_", [300], 0.69)
#
# # mass_plots[8], radius_plots[8] = planet_interp(data_files_stored_in, "_adiabatic_", [1000], 0.69)
# # mass_plots[9], radius_plots[9] = planet_interp(data_files_stored_in, "_constant_", [1000], 0.69)
#
# mass_plots[10], radius_plots[10] = planet_interp(data_files_stored_in, "_adiabatic_", [3000], 0.69)
# mass_plots[11], radius_plots[11] = planet_interp(data_files_stored_in, "_constant_", [3000], 0.69)
#
#
# mass_plots[12], radius_plots[12] = planet_interp(data_files_stored_in, "_adiabatic_", [300], 0.1)
# mass_plots[13], radius_plots[13] = planet_interp(data_files_stored_in, "_constant_", [300], 0.1)
#
# # mass_plots[14], radius_plots[14] = planet_interp(data_files_stored_in, "_adiabatic_", [1000], 0.1)
# # mass_plots[15], radius_plots[15] = planet_interp(data_files_stored_in, "_constant_", [1000], 0.1)
#
# mass_plots[16], radius_plots[16] = planet_interp(data_files_stored_in, "_adiabatic_", [3000], 0.1)
# mass_plots[17], radius_plots[17] = planet_interp(data_files_stored_in, "_constant_", [3000], 0.1)
#
# mass_plots[18], radius_plots[18] = planet_interp(data_files_stored_in, "_adiabatic_", [300], 0.8)
# mass_plots[19], radius_plots[19] = planet_interp(data_files_stored_in, "_constant_", [300], 0.8)
#
# # mass_plots[20], radius_plots[20] = planet_interp(data_files_stored_in, "_adiabatic_", [1000], 0.8)
# # mass_plots[21], radius_plots[21] = planet_interp(data_files_stored_in, "_constant_", [1000], 0.8)
#
# mass_plots[22], radius_plots[22] = planet_interp(data_files_stored_in, "_adiabatic_", [3000], 0.8)
# mass_plots[23], radius_plots[23] = planet_interp(data_files_stored_in, "_constant_", [3000], 0.8)
#


# earth = [0,2,3,4,5]
# mercury = [6,8,9,10,11]
#
# plt.plot(mass_plots[0], radius_plots[0], color="blue", linestyle="dashed")
# plt.plot(mass_plots[1], radius_plots[1], color="blue", linestyle=(0, (0.3, 0.3)))
#
# plt.plot(mass_plots[2], radius_plots[2], color="blue", linestyle="dotted")
# plt.plot(mass_plots[3], radius_plots[3], color="blue", linestyle=(0, (0.3, 0.3)))
#
# plt.plot(mass_plots[4], radius_plots[4], color='blue', linestyle="dashdot")
# plt.plot(mass_plots[5], radius_plots[5], color="blue", label="Earth: CMF = 0.33", linestyle="solid")
#
#
# plt.plot(mass_plots[6], radius_plots[6], color="orange", linestyle="dashed")
# plt.plot(mass_plots[7], radius_plots[7], color="orange", linestyle=(0, (0.3, 0.3)))
#
# plt.plot(mass_plots[8], radius_plots[8], color="orange", linestyle="dotted")
# plt.plot(mass_plots[9], radius_plots[9], color="orange", linestyle=(0, (0.3, 0.3)))
#
# plt.plot(mass_plots[10], radius_plots[10], color='orange', linestyle="dashdot")
# plt.plot(mass_plots[11], radius_plots[11], color="orange", label="Mercury: CMF = 0.69", linestyle="solid")
#
# plt.ylim(0.9, 1.6)
#
# plt.plot([], [], color="black", label="300 K Adiabatic", linestyle="dashed")
# plt.plot([], [], color="black", label="1000 K Adiabatic", linestyle="dotted")
# plt.plot([], [], color="black", label="3000 K Adiabatic", linestyle="dashdot")
#
# plt.plot([], [], color="black", label="1000 K Constant", linestyle=(0, (0.3, 0.3)))
# plt.plot([], [], color="black", label="3000 K Constant", linestyle="solid")

# with open('/Users/sabrinaberger/All Research/RockyPlanets/planets.csv') as csvfile:
#     readCSV = csv.reader(csvfile, delimiter=',')
#     mass_planets = []
#     mass_upper = []
#     mass_lower = []
#
#     radius_planets = []
#     radius_upper = []
#     radius_lower = []
#
#     for row in readCSV:
#         if float(row[22])/float(row[21]) < 0.4 and float(row[27])/float(row[26]) < 0.4:
#             mass_planets.append(float(row[21])*317)
#             print(row[21])
#             mass_upper.append(float(row[22])*317)
#             print(row[22])
#             mass_lower.append(-float(row[23])*317)
#
#             radius_planets.append(float(row[26])*11.2)
#             radius_upper.append(float(row[27])*11.2)
#             radius_lower.append(-float(row[28]) * 11.2)
#
#     mass_planets = np.array(mass_planets)
#     radius_planets = np.array(radius_planets)
#
#
# # This can also be done by calling plt.figure() and not setting it equal to fig
# # fig = plt.figure()
# # ax = fig.add_subplot(111)
# # ax.scatter(mass_planets, radius_planets, marker='o', color='b')
# # # If you call errorbar without linestyle="None", it will connect all the
# # # points and look weird
# # ax.errorbar(mass_planets, radius_planets, xerr=(mass_lower, mass_upper), yerr=(radius_lower, radius_upper))
#
# plt.errorbar(mass_planets, radius_planets, xerr=(mass_lower, mass_upper), yerr=(radius_lower, radius_upper), c = 'k', marker = 'o', linestyle = 'None', capsize=2, markersize=3, zorder = 20)
#
#
# plt.xlabel("Mass ($M_{Earth}$)")
# plt.ylabel("Radius ($R_{Earth}$)")
# # plt.ylim(1, 1.7)
# plt.legend(loc=0)
# plt.title("Mass-Radius Relationships")
# plt.savefig("/Users/sabrinaberger/Desktop/mass_radius_all_with_obs.pdf")
# plt.show()


# for i in earth:
#     radius_plots[i] = 100*(radius_plots[i]-radius_plots[1])/radius_plots[1]
#
#
# for i in mercury:
#     radius_plots[i] = 100*(radius_plots[i]-radius_plots[7])/radius_plots[7]
#


plt.plot(mass_plots[0], radius_plots[0], color="blue", linestyle="dashed")
plt.plot(mass_plots[1], radius_plots[1], color="blue", linestyle=(0, (0.3, 0.3)))

# plt.plot(mass_plots[2], radius_plots[2], color="blue", linestyle="dotted")
# plt.plot(mass_plots[3], radius_plots[3], color="blue", linestyle=(0, (0.3, 0.3)))

plt.plot(mass_plots[4], radius_plots[4], color='blue', linestyle="dashdot")
plt.plot(mass_plots[5], radius_plots[5], color="blue", label="Earth: CMF = 0.33", linestyle="solid")


plt.plot(mass_plots[6], radius_plots[6], color="orange", linestyle="dashed")
plt.plot(mass_plots[7], radius_plots[7], color="orange", linestyle=(0, (0.3, 0.3)))

# plt.plot(mass_plots[8], radius_plots[8], color="orange", linestyle="dotted")
# plt.plot(mass_plots[9], radius_plots[9], color="orange", linestyle=(0, (0.3, 0.3)))

plt.plot(mass_plots[10], radius_plots[10], color='orange', linestyle="dashdot")
plt.plot(mass_plots[11], radius_plots[11], color="orange", label="Mars: CMF = 0.69", linestyle="solid")

plt.plot(mass_plots[12], radius_plots[12], color="purple", linestyle="dashed")
plt.plot(mass_plots[13], radius_plots[13], color="purple", linestyle=(0, (0.3, 0.3)))
#
# plt.plot(mass_plots[14], radius_plots[14], color="purple", linestyle="dotted")
# plt.plot(mass_plots[15], radius_plots[15], color="purple", linestyle=(0, (0.3, 0.3)))

plt.plot(mass_plots[16], radius_plots[16], color='purple', linestyle="dashdot")
plt.plot(mass_plots[17], radius_plots[17], color="purple", label="Silicate Planet: CMF = 0.1", linestyle="solid")

plt.plot(mass_plots[18], radius_plots[18], color="red", linestyle="dashed")
plt.plot(mass_plots[19], radius_plots[19], color="red", linestyle=(0, (0.3, 0.3)))
#
# plt.plot(mass_plots[20], radius_plots[20], color="red", linestyle="dotted")
# plt.plot(mass_plots[21], radius_plots[21], color="red", linestyle=(0, (0.3, 0.3)))

plt.plot(mass_plots[22], radius_plots[22], color='red', linestyle="dashdot")
plt.plot(mass_plots[23], radius_plots[23], color="red", label="Iron Planet: CMF = 0.8", linestyle="solid")

plt.show()

# with open('/Users/sabrinaberger/All Research/RockyPlanets/planets.csv') as csvfile:
#     readCSV = csv.reader(csvfile, delimiter=',')
#     mass_planets = []
#     mass_upper = []
#     mass_lower = []
#
#     radius_planets = []
#     radius_upper = []
#     radius_lower = []
#     planet_names = []
#
#     star_temp = []
#     star_rad = []
#     a = []
#
#     for row in readCSV:
#         if row[54] is not "" and row[64] is not "" and row[22] is not "" and row[21] is not "" and row[26] is not "" and row[27] is not "" and float(row[22])/float(row[21]) < 0.4 and float(row[27])/float(row[26]) < 0.4:
#             mass_planets.append(float(row[21])*317)
#             print("mass of " + row[1] + " is " + row[21])
#             # print(row[21])
#             mass_upper.append(float(row[22])*317)
#             # print(row[22])
#             mass_lower.append(-float(row[23])*317)
#
#             radius_planets.append(float(row[26])*11.2)
#             radius_upper.append(float(row[27])*11.2)
#             radius_lower.append(-float(row[28]) * 11.2)
#             planet_names.append(row[1])
#
#             star_temp.append(float(row[54]))
#
#             star_rad.append(float(row[64]))
#             a.append(float(row[9]))
#
#     mass_planets = np.array(mass_planets)
#     radius_planets = np.array(radius_planets)
#     star_temp = np.array(star_temp)
#     star_rad = np.array(star_rad)
#     star_rad = star_rad*6.957e8
#     a = np.array(a)
#     a_m = a*1.496e11
#
# T_eq = star_temp * (star_rad/(2*a_m))**(0.5)
# print(len(T_eq))
#

# plt.errorbar(mass_planets, radius_planets, xerr=(mass_lower, mass_upper), yerr=(radius_lower, radius_upper), c = 'k', marker = 'o', linestyle = 'None', capsize=2, markersize=3, zorder = 20)
# i=0

#
# bbox_props = dict(boxstyle="square", fc="white", ec="w", lw=2)
#
#
# while i != len(planet_names) - 3:
#     print(i)
#     plt.text(mass_planets[i], radius_planets[i], planet_names[i], ha="left", va="bottom",
#             size=5,
#             bbox=bbox_props)
#     i += 1


print(max(a))
plt.plot([], [], color="black", label="300 K Adiabatic", linestyle="dashed")
plt.plot([], [], color="black", label="1000 K Adiabatic", linestyle="dotted")
plt.plot([], [], color="black", label="3000 K Adiabatic", linestyle="dashdot")

plt.plot([], [], color="black", label="300 K Constant", linestyle=(0, (0.3, 0.3)))
plt.plot([], [], color="black", label="1000 K Constant", linestyle=(0, (0.3, 0.3)))
plt.plot([], [], color="black", label="3000 K Constant", linestyle="solid")


my_cmap_r = 'viridis_r'
a = np.log10(a)

properties = plt.scatter(radius_planets, T_eq, c = a, s=10, cmap = my_cmap_r)
cbar = plt.colorbar(properties)
cbar.set_label('log(Semi-Major Axis) [AU]')

plt.xlabel("Radius [$R_{\oplus}$]", fontsize = 'large')
plt.ylabel("$T_{equilibrium}$ [K]", fontsize = 'large')
plt.title("Distribution of Equilibrium Temperatures of Exoplanets Discovered", fontsize = 'medium')

#
# plt.savefig("/Users/sabrinaberger/Desktop/teq.pdf", dpi=300)
# plt.show()
# plt.close()
#
# plt.ylim(1, 1.7)
# plt.legend(loc=0, markerscale=0.2, fontsize=8)
# plt.title("Mass-Radius Relationships")
# plt.savefig("/Users/sabrinaberger/Desktop/plot_all.pdf")
# plt.show()
