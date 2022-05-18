import matplotlib.pyplot as plt
from numpy import *
from matplotlib import ticker, colors
import astropy.constants as const

plt.rc('font', family='serif')

# Property Plots

class PlotBasic:

    def __init__(self, temp, typ, p_c, p_cmb_percentage, radius, mass, temp_diff=0, radius_lower="",
                 path="../paper_plots/", difference=True, surf_press="", core_mass="", core_rad="", p_cmb_percentage_simulated=""):

        self.temperature = temp
        self.temp_diff = temp_diff
        self.typ = typ
        self.path = path
        self.difference = difference

        # collecting places where files are stored and loading grid
        self.p_c = load(p_c)
        self.p_cmb_percentage = load(p_cmb_percentage)
        self.radius = load(radius)/const.R_earth.value

        if self.difference:
            self.radius_lower = load(radius_lower)/const.R_earth.value
            self.radius_difference = (self.radius - self.radius_lower)/self.radius_lower

        self.mass = load(mass)/const.M_earth.value
        self.num_planets = len(self.mass.flatten())
        # self.surf_press = load(surf_press)
        # self.core_mass = load(core_mass)
        # self.core_rad = load(core_rad)
        # self.p_cmb_percentage_simulated = load(p_cmb_percentage_simulated)

        self.contourify()  # where x axis is automatically log(P_c) and y axis is p_cmb_percentage/P_c

        # generate all contour plots with contourify

    def makePropertyPlot(self, propertyName, units, propertyArray, plotTitle, filename, difference=False, colorscheme='magma_r'):
        if difference:
            property = plt.contourf(log10(self.p_c), self.p_cmb_percentage, propertyArray, cmap=plt.get_cmap(colorscheme))
        else:
            property = plt.contourf(log10(self.p_c), self.p_cmb_percentage, log10(propertyArray), cmap=plt.get_cmap(colorscheme))

        cbar = plt.colorbar(property)
        cbar.ax.set_ylabel(propertyName + ' ' + units)

        plt.xlabel(r'$\log(P_c)~[{\rm Pa}]$')
        plt.ylabel(r'$P_{\rm CMB}/P_c$')
        plt.title(plotTitle)
        plt.tight_layout()
        print(self.path + filename)
        plt.savefig(self.path + filename, dpi=600)
        plt.close()

    def contourify(self):
        t_anchor_str = r"$T_{\rm anchor}$"
        self.makePropertyPlot(r'$\log(R_{pl})$', r'[$R_{\bigoplus}$]', self.radius, f"Total Radius for {self.num_planets} Planets ({self.typ}, {t_anchor_str} = {self.temperature}K)",
                         f"radius_{self.typ}_{self.temperature}.pdf")
        self.makePropertyPlot(r'$\log(M_{pl})$', r'[$M_{\bigoplus}$]', self.mass, f"Total Mass for {self.num_planets} Planets ({self.typ}, {t_anchor_str} = {self.temperature}K)", f"mass_{self.typ}_{self.temperature}.pdf")
        if self.difference:
            self.makePropertyPlot(r'$(R_{3000} - R_{300}) / R_{300}$', '', self.radius_difference, f'Radial Difference in {self.temperature}K {self.typ}', f'difference_{self.typ}_{self.temperature}.pdf', difference=True)


temp_diff = 300

#
# for temp in [300, 3000]:
#     direc_with_phases_adiabatic = f"/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/RockyPlanets/paper_data/complete_data_with_silicate_mantle/{temp}_adiabatic_data/"
#     p_c = direc_with_phases_adiabatic + f"p_c_grid_adiabatic_{temp}.pyc.npy"
#     p_cmb_percentage = direc_with_phases_adiabatic + f"p_cmb_percentage_grid_adiabatic_{temp}.pyc.npy"
#     radius = direc_with_phases_adiabatic + f"radius_grid_adiabatic_{temp}.pyc.npy"
#     mass = direc_with_phases_adiabatic + f"mass_grid_adiabatic_{temp}.pyc.npy"
#     PlotBasic(temp, "Adiabatic", p_c, p_cmb_percentage, radius, mass, temp_diff=temp_diff, radius_lower="", difference=False)
#
# diff_direc_with_phases_constant = f"/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/RockyPlanets/paper_data/complete_data_with_silicate_mantle/{temp_diff}_constant_data/"
diff_temp = 300
types = ["constant", "adiabatic"]
for typ in types:
    for temp in [3000]:
        direc_with_phases = f"/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/RockyPlanets/paper_data/complete_data_with_silicate_mantle/{temp}_{typ}_data/"
        diff_direc_with_phases = f"/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/RockyPlanets/paper_data/complete_data_with_silicate_mantle/{diff_temp}_constant_data/"
        p_c = direc_with_phases + f"p_c_grid_{typ}_{temp}.pyc.npy"
        p_cmb_percentage = direc_with_phases + f"p_cmb_percentage_grid_{typ}_{temp}.pyc.npy"
        radius = direc_with_phases + f"radius_grid_{typ}_{temp}.pyc.npy"
        mass = direc_with_phases + f"mass_grid_{typ}_{temp}.pyc.npy"
        radius_lower = diff_direc_with_phases + f"radius_grid_constant_{diff_temp}.pyc.npy"
        PlotBasic(temp, typ.title(), p_c, p_cmb_percentage, radius, mass, radius_lower=radius_lower, difference=True)
#

# Optimization Plots
#
# x = load("DataFiles/OPTrelativeTolerance.pyc.npy")
# radii = load("DataFiles/OPTtotalRadii.pyc.npy")
# radii1000 = load("DataFiles/1000totalRadii.pyc.npy")
# masses = load("DataFiles/OPTtotalMasses.pyc.npy")
# masses1000 = load("DataFiles/1000totalMasses.pyc.npy")
# pressures = load("DataFiles/OPTtotalPressures.pyc.npy")
# times = load("DataFiles/OPTtotalTime.pyc.npy")

# relative_tolerance, totalTime, totalRadii, totalMasses, totalPressures = loadtxt("optimize_grid_data.txt", dtyp='float', usecols=(0,1,2,3,4), unpack=True)

# print(masses[0], masses[-1])
# print(masses1000[0], masses1000[-1])
# print(pressures[0], pressures[-1])
# print(times[0], times[-1])

# Mass Plot

# masses = masses/(float("1e24"))
# masses1000 = masses1000/(float("1e24"))
#
# xmin = x[0]
# xmax = x[-1]
#
# ymin = min(min(masses), min(masses1000))
# ymax = max(max(masses), max(masses1000))
#
# masses = plt.semilogx(x, masses, linewidth=2, color="m", label='dt = 100')
# masses1000 = plt.semilogx(x, masses1000, linewidth=2, color="r", label='dt = 1000')
# plt.xlabel("Relative Tolerance")
# plt.ylabel("Mass/(1e24) (kg)")
#
# axes = plt.gca()
# plt.legend(loc=3)
# axes.set_xlim([xmin, xmax])
# axes.set_ylim([ymin-0.0001, ymax+0.0001])
#
# plt.savefig('/Users/sabrinaberger/RockyPlanets/Code/Plots_Optimize/optimize_reltol_mass.pdf')
# plt.close()
#
# # Radius Plot
# radii = radii/float("1e6")
# radii1000 = radii1000/float("1e6")
#
#
# ymin = min(min(radii), min(radii1000))
# ymax = max(max(radii), max(radii1000))
#
# radii = plt.semilogx(x, radii, linewidth=2, color="m", label = 'dt = 100')
# radii1000 = plt.semilogx(x, radii1000, linewidth=2, color = "r", label = 'dt = 1000')
# plt.xlabel("Relative Tolerance")
# plt.ylabel("Radii/(1e6) (m)")
#
# axes = plt.gca()
# plt.legend(loc=3)
# axes.set_xlim([xmin, xmax])
# axes.set_ylim([ymin-0.0001, ymax+0.0001])
#
# plt.savefig('/Users/sabrinaberger/RockyPlanets/Code/Plots_Optimize/optimize_reltol_radii.pdf')
# # plt.show()
# plt.close()
#


# WEIRD TWO PLOT THING BELOW
# host = host_subplot(111, axes_class=AA.Axes)
# plt.subplots_adjust(right=0.75)
#
# par1 = host.twinx()
#
# offset = 60
#
# host.set_xlabel("dt (m)")
# host.set_ylabel("Radius (m)")
# par1.set_ylabel("Mass (kg)")
#
# p1, = host.plot(x, log(radii), label="log(Radius)", linewidth=2, color = "m")
# p2, = par1.plot(x, log(masses), label="log(Mass)", linewidth=2, color = "g")
#
# host.legend()
#
# host.set_ylim(1e5, 3e6)
# par1.set_ylim(min(masses), max(masses))
#
# host.axis["left"].label.set_color(p1.get_color())
# par1.axis["right"].label.set_color(p2.get_color())
#
# plt.draw()
