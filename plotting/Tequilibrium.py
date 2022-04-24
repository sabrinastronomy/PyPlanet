import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib
import numpy as np
import pandas as pd
import math

plt.rc('font', family='serif')

# helper functions

def convert_stellar(radius):
    # input radius in stellar radius
    return radius * 6.957e8 # returned value in meters

def convert_AU(semi):
    # input radius in stellar radius
    return semi * 1.496e11 # returned value in meters

def calculate_Teq_and_err(star_temp, star_rad, a_m, star_temp_lim, star_rad_lim, a_m_lim):
    T_eq = star_temp * (star_rad / (2 * a_m)) ** (0.5)
    T_eq_error_lim = np.sqrt(((star_rad / (2 * a_m)) * np.power(star_temp_lim, 2)) + (
                ((star_temp ** 2) / (8 * a_m * star_rad)) * np.power(star_rad_lim, 2)) + (
                                           a_m_lim ** 2 * ((star_temp ** 2 * star_rad) / (8 * a_m ** 3))))
    return T_eq, T_eq_error_lim

# COLUMN pl_name:        Planet Name
# COLUMN pl_orbsmax:     Orbit Semi-Major Axis [au])
# COLUMN pl_orbsmaxerr1: Orbit Semi-Major Axis Upper Unc. [au]
# COLUMN pl_orbsmaxerr2: Orbit Semi-Major Axis Lower Unc. [au]
# COLUMN pl_orbsmaxlim:  Orbit Semi-Major Axis Limit Flag
# COLUMN pl_rade:        Planet Radius [Earth Radius]
# COLUMN pl_radeerr1:    Planet Radius Upper Unc. [Earth Radius]
# COLUMN pl_radeerr2:    Planet Radius Lower Unc. [Earth Radius]
# COLUMN pl_radelim:     Planet Radius Limit Flag
# COLUMN pl_masse:       Planet Mass [Earth Mass]
# COLUMN pl_masseerr1:   Planet Mass [Earth Mass] Upper Unc.
# COLUMN pl_masseerr2:   Planet Mass [Earth Mass] Lower Unc.
# COLUMN pl_masselim:    Planet Mass [Earth Mass] Limit Flag
# COLUMN pl_eqt:         Equilibrium Temperature [K]
# COLUMN pl_eqterr1:     Equilibrium Temperature Upper Unc. [K]
# COLUMN pl_eqterr2:     Equilibrium Temperature Lower Unc. [K]
# COLUMN pl_eqtlim:      Equilibrium Temperature Limit Flag
# COLUMN st_teff:        Stellar Effective Temperature [K]
# COLUMN st_tefferr1:    Stellar Effective Temperature Upper Unc. [K]
# COLUMN st_tefferr2:    Stellar Effective Temperature Lower Unc. [K]
# COLUMN st_tefflim:     Stellar Effective Temperature Limit Flag
# COLUMN st_rad:         Stellar Radius [Solar Radius]
# COLUMN st_raderr1:     Stellar Radius Upper Unc. [Solar Radius]
# COLUMN st_raderr2:     Stellar Radius Lower Unc. [Solar Radius]
# COLUMN st_radlim:      Stellar Radius Limit Flag

file_archive = pd.read_csv("cumulative-1.csv", skiprows=3)

radius_planets = file_archive["pl_rade"].values
radius_upper = file_archive["pl_radeerr1"].values
radius_lower = file_archive["pl_radeerr2"].values

mass_planets = file_archive["pl_masse"].values
mass_upper = file_archive["pl_masseerr1"].values
mass_lower = file_archive["pl_masseerr2"].values

star_temp = file_archive["st_teff"].values
star_temp_upper = file_archive["st_tefferr1"].values
star_temp_lower = file_archive["st_tefferr2"].values

star_rad = convert_stellar(file_archive["st_rad"].values) # converting from solar radius to meters
star_rad_upper = convert_stellar(file_archive["st_raderr1"].values) # converting from solar radius to meters
star_rad_lower = convert_stellar(file_archive["st_raderr2"].values) # converting from solar radius to meters

a = file_archive["pl_orbsmax"].values
a_m = convert_AU(a)
a_upper = convert_AU(file_archive["pl_orbsmaxerr1"].values)
a_lower = convert_AU(file_archive["pl_orbsmaxerr2"].values)

T_eq_arr, T_eq_error_upper = calculate_Teq_and_err(star_temp, star_rad, a_m, star_temp_upper, star_rad_upper, a_upper)
T_eq_arr, T_eq_error_lower = calculate_Teq_and_err(star_temp, star_rad, a_m, star_temp_lower, star_rad_lower, a_lower)

a_m_Earth = convert_AU(1)
star_rad_Earth = convert_stellar(1)
star_temp_Earth = 5800
Teq_Earth, _ = calculate_Teq_and_err(star_temp_Earth, star_rad_Earth, a_m_Earth, 0, 0, 0)


a_m_Venus = convert_AU(1)
Teq_Earth, _ = calculate_Teq_and_err(star_temp_Earth, star_rad_Earth, a_m_Earth, 0, 0, 0)

##### - see here first


my_cmap = 'jet'
my_cmap_r = 'viridis_r'
fig = plt.figure()
ax = fig.add_subplot(111)
T_eq_arr = T_eq_arr

im = ax.scatter(radius_planets, mass_planets, c=T_eq_arr, s=20, cmap=my_cmap_r, norm=matplotlib.colors.LogNorm())
ax.scatter(1, 1, c = Teq_Earth, label="Earth", cmap=my_cmap_r, s=100, alpha = 0.7)
plt.text(1.5, -0.3,'Earth')
plt.arrow(1.5, 0, -0.4, 0.8, width =0.04, color='k')
# formatter = matplotlib.ticker.LogFormatter(10, labelOnlyBase=False)
cbar = fig.colorbar(im, ticks=[100, 200, 300, 400, 500, 700, 1000, 2000, 3000, 4000, 5000, 7000], format='%.0f')

cbar.ax.set_ylabel(r'$T_{equilibrium} [T_{eq, \oplus}]$')

ax.set_xlabel("Radius [$R_{\oplus}$]", fontsize = 'large')
ax.set_ylabel("Mass [$M_{\oplus}$]", fontsize = 'large')

plt.show()
fig.savefig("/Users/sabrinaberger/Desktop/teq.png", dpi=350, bbox_inches='tight')
#####