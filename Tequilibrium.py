import matplotlib.pyplot as plt
import csv
import numpy as np

# I pull this csv file from the NASA Exoplanet Archive
with open('/Users/sabrinaberger/All Research/RockyPlanets/cumulative-2.csv') as file:
    # readCSV = csv.reader(file, delimiter=',')

    planet_id = []

    radius_planets = []
    radius_upper = []
    radius_lower = []

    a = []
    a_upper = []
    a_lower = []

    t_eq_given = []

    star_temp = []
    star_temp_upper = []
    star_temp_lower = []

    star_rad = []
    star_rad_upper = []
    star_rad_lower = []


    i = 0
    for line in file:
        if i >= 23 and i <= 100:
            line = line.strip()
            linesplit = line.split(',')

            planet_id.append(float(linesplit[0]))

            radius_planets.append(float(linesplit[1]))
            radius_upper.append(float(linesplit[2]))
            radius_lower.append(float(linesplit[3]))

            a.append(float(linesplit[4]))

            if linesplit[5] is not "":
                a_upper.append(float(linesplit[5]))
            else:
                a_upper.append(float(0))

            if linesplit[6] is not "":
                a_lower.append(float(linesplit[6]))
            else:
                a_lower.append(0)

            t_eq_given.append(float(linesplit[7]))

            star_temp.append(float(linesplit[10]))

            if linesplit[11] is not "":
                star_temp_upper.append(float(linesplit[11]))
            else:
                star_temp_upper.append(0)

            if linesplit[12] is not "":
                star_temp_lower.append(float(linesplit[12]))
            else:
                star_temp_lower.append(0)

            star_rad.append(float(linesplit[13]))

            if linesplit[14] is not "":
                star_rad_upper.append(float(linesplit[14]))
            else:
                star_rad_upper.append(0)

            if linesplit[15] is not "":
                print(i)
                star_rad_lower.append(float(linesplit[15]))
            else:
                star_rad_lower.append(0)

        i += 1

    # for row in readCSV:
    #     #
    #     if row[54] is not "" and row[64] is not "" and row[22] is not "" and row[21] is not "" and row[26] is not "" and row[27] is not "" and row[9] is not ""  and float(row[22]) / float(row[21]) < 0.4 and float(row[27]) / float(row[26]) < 0.4:
    #         mass_planets.append(float(row[21])*317)
    #         mass_upper.append(float(row[22])*317)
    #         mass_lower.append(abs(float(row[23])*317))
    #
    #         radius_planets.append(float(row[26])*11.2)
    #         radius_upper.append(float(row[27])*11.2)
    #         radius_lower.append(abs(float(row[28]) * 11.2))
    #         planet_names.append(row[1])
    #
    #         star_temp.append(float(row[54]))
    #         if row[55] is "":
    #             star_temp_upper_lower.append(0)
    #         else:
    #             star_temp_upper_lower.append(abs(float(row[55])))
    #
    #         star_rad.append(float(row[64]))
    #         if row[65] is "":
    #             star_rad_upper.append(0)
    #         else:
    #             star_rad_upper.append(abs(float(row[65])))
    #
    #         if row[66] is "":
    #             star_rad_lower.append(0)
    #         else:
    #             star_rad_lower.append(abs(float(row[66])))
    #
    #         a.append(float(row[9]))
    #         if row[10] is "":
    #             a_upper.append(0)
    #         else:
    #             a_upper.append(abs(float(row[10])))
    #
    #         if row[11] is "":
    #             a_lower.append(0)
    #         else:
    #             a_lower.append(abs(float(row[11])))

    radius_planets = np.array(radius_planets)
    radius_upper = np.array(radius_upper)
    radius_lower = -np.array(radius_lower)

    star_temp = np.array(star_temp)
    star_temp_upper = np.array(star_temp_upper)
    star_temp_lower = np.array(star_temp_lower)

    star_rad = np.array(star_rad)
    star_rad = star_rad*6.957e8
    star_rad_upper = np.array(star_rad_upper)
    star_rad_upper = star_rad_upper*6.957e8
    star_rad_lower = np.array(star_rad_lower)
    star_rad_lower = star_rad_lower*6.957e8

    a = np.array(a)
    a_m = a*1.496e11 # semi-major axis in meters
    a_upper = np.array(a_upper)
    a_upper = a_upper*1.496e11
    a_lower = np.array(a_lower)
    a_lower = a_lower*1.496e11

##### - see here first

T_eq = star_temp * (star_rad/(2*a_m))**(0.5)
T_eq_error_upper = np.sqrt(((star_rad/(2*a_m))*np.power(star_temp_upper, 2)) + (((star_temp**2)/(8*a_m*star_rad))*np.power(star_rad_upper, 2)) + (a_upper**2 * ((star_temp**2 * star_rad)/(8 * a_m**3))))
T_eq_error_lower = np.sqrt(((star_rad/(2*a_m))*np.power(star_temp_lower, 2)) + (((star_temp**2)/(8*a_m*star_rad))*np.power(star_rad_lower, 2)) + (a_lower**2 * (star_temp**2 * star_rad)/(8 * a_m**3)))


my_cmap_r = 'viridis_r'
a = np.log10(a)

fig = plt.figure()
ax = fig.add_subplot(111)

# im = ax.scatter(radius_planets, T_eq, c = a, cmap = my_cmap_r, zorder=100)
# fig.colorbar(im)
ax.errorbar(radius_planets, T_eq, xerr = [radius_lower, radius_upper], yerr = [T_eq_error_lower, T_eq_error_upper], ecolor = "b", zorder = 0, linestyle = "None", lw=0.5, fmt='o', capsize = 10)


ax.set_label('log(Semi-Major Axis) [AU]')
ax.set_xlabel("Radius [$R_{\oplus}$]", fontsize = 'large')
ax.set_ylabel("$T_{equilibrium}$ [K]", fontsize = 'large')
ax.set_title("Distribution of Equilibrium Temperatures")


fig.savefig("/Users/sabrinaberger/Desktop/teq.pdf", dpi=1200)
fig.show()
#####