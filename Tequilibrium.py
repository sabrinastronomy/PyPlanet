import matplotlib.pyplot as plt
import numpy as np

# I pull this csv file from the NASA Exoplanet Archive
with open('/Users/sabrinaberger/All Research/RockyPlanets/Planets < 2 Re/cumulative-2.csv') as file:

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
        if i >= 24:
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

my_cmap = 'Greys'
my_cmap_r = 'Greys_r'

a = np.log10(a)

fig = plt.figure()
ax = fig.add_subplot(111)

quad = np.ndarray(shape = len(radius_lower))

for i in range(len(radius_lower)):
    quad[i] = np.sqrt(max(radius_lower[i], radius_upper[i])**2 + max(T_eq_error_lower[i], T_eq_error_upper[i]))
    i+=1


im = ax.scatter(radius_planets, T_eq, c = quad, cmap = my_cmap_r, zorder=100)
# cbar = fig.colorbar(im)
# ax.errorbar(radius_planets, T_eq, xerr = [radius_lower, radius_upper], yerr = [T_eq_error_lower, T_eq_error_upper], ecolor = "k", zorder = 0, linestyle = "None", lw=0.5, fmt='o')
# cbar.ax.set_ylabel('$\sqrt{\delta a ^2 + \delta T_{eq}^2}$')

ax.set_xlabel("Radius [$R_{\oplus}$]", fontsize = 'large')
ax.set_ylabel("$T_{equilibrium}$ [K]", fontsize = 'large')


plt.show()
fig.savefig("/Users/sabrinaberger/Desktop/teq.pdf", dpi=1200)
#####