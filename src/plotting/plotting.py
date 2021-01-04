import matplotlib.pyplot as plt
from numpy import *
from matplotlib import ticker, colors

# Property Plots

class PlotBasic:

    def __init__(self, temp, type, p_c, p_cmb_percentage, radius, mass, surf_press="", core_mass="", core_rad="", p_cmb_percentage_simulated="", radius_higher="",
                 path="/Users/sabrinaberger/All Research/RockyPlanets/2017Plots/"):

        self.temperature = temp
        self.type = type
        self.path = path

        # collecting places where files are stored and loading grid
        self.p_c = load(p_c)
        self.p_cmb_percentage = load(p_cmb_percentage)

        self.radius = load(radius)

        if radius_higher is not "":
            self.radius_higher = load(radius_higher)
            print(max(self.radius_higher.any()))
            print(max(self.radius))
            self.radius_difference = (self.radius_higher - self.radius)/self.radius
        print(self.radius_difference)
        self.mass = load(mass)
        self.surf_press = load(surf_press)
        self.core_mass = load(core_mass)
        self.core_rad = load(core_rad)
        self.p_cmb_percentage_simulated = load(p_cmb_percentage_simulated)

        self.contourify()  # where x axis is automatically log(P_c) and y axis is p_cmb_percentage/P_c

        # generate all contour plots with contourify

    def contourify(self):
        def makePropertyPlot(propertyName, units, propertyArray, plotTitle, filename, colorscheme='magma_r', difference = True):

            if difference:
                property = plt.contourf(log10(self.p_c), self.p_cmb_percentage, propertyArray,
                                        cmap=plt.get_cmap(colorscheme))
            else:
                property = plt.contourf(log10(self.p_c), self.p_cmb_percentage, log10(propertyArray),  cmap=plt.get_cmap(colorscheme))
            # property.set_xscale('log')
            # property.set_yscale('log')
            cbar = plt.colorbar(property)
            cbar.ax.set_ylabel(propertyName + ' ' + units)



            # levs = np.power(10, lev_exp)

            plt.xlabel('log($P_c$) (Pa)')
            plt.ylabel('$P_{CMB}$/$P_c$')
            # plt.xlim(10, 12)
            plt.title(plotTitle)
            plt.savefig(self.path + filename, dpi=600)
            plt.close()

        makePropertyPlot('log(Radius)', '(m)', self.radius, 'Total Radius ' + str(self.type) + " " + str(self.temperature) + ' K',
                         'radius_' + self.type + "_" + str(self.temperature) + '.pdf')
        makePropertyPlot('log(Mass)', '(kg)', self.mass, 'Total Mass ' + str(self.type) + " " + str(self.temperature) + ' K',
                         'mass_' + self.type + "_" + str(self.temperature) + '.pdf')
        makePropertyPlot('(3000K Radius - 300K Radius)/300K Radius', '', self.radius_difference, 'Radial difference - 3000 K and 300 K Constant', 'DIFFERENCE.pdf')
#
# class CompareRadius(PlotBasic):
#
#     def __init__(self, other_radius, other_type, other_temperature):
#         PlotBasic.__init__(self, temperature, type, radius, path)
#         self.other_radius = other_radius
#         self.other_type = other_type
#         self.other_temperature = other_temperature
#         self.difference = self.other_radius - self.radius


#
#
#
#     property = plt.contourf(x_1, y_1, rad_grid_diff, cmap=plt.get_cmap('magma_r'))
#     cbar = plt.colorbar(property)
#     cbar.ax.set_ylabel("(3000K Radius - 300K Radius)/300K Radius")
#     plt.xlabel('log($P_c$) (Pa)')
#     plt.xlim(10, 12)
#     plt.ylabel('$P_{CMB}$/$P_c$')
#     plt.savefig('/Users/sabrinaberger/Desktop/2017_Paper_Plots/difference.pdf')
#     plt.close()
#     plt.show()
#
#
# radius_grid300 = load("/Users/sabrinaberger/RockyPlanets/DataFiles/radius_grid300.pyc.npy")
# radius_grid3000 = load("/Users/sabrinaberger/RockyPlanets/DataFiles/radius_grid3000.pyc.npy")
# rad_grid_diff = (radius_grid3000 - radius_grid300)/radius_grid300
#



# prop_1 = load("/Users/sabrinaberger/RockyPlanets/MassRadiusDiagramData/300_0.26.pyc.npy")
# prop_2 = load("/Users/sabrinaberger/RockyPlanets/MassRadiusDiagramData/300_0.33.pyc.npy")
# print(prop_1)
#
# print([prop_1[7], prop_2[7]], [prop_1[6], prop_2[6]])

# core_mass_1 = load("DataFiles/core_mass_grid_adiabatic_300.0.pyc.npy")
# core_mass_2 = load("DataFiles/core_mass_grid_adiabatic_600.0.pyc.npy")
#
# core_press_1 = load("DataFiles/p_c_adiabatic_300.0.pyc.npy")
#
# print(core_mass_1)
# print(core_press_1)

# Optimization Plots
#
# x = load("DataFiles/OPTrelativeTolerance.pyc.npy")
# radii = load("DataFiles/OPTtotalRadii.pyc.npy")
# radii1000 = load("DataFiles/1000totalRadii.pyc.npy")
# masses = load("DataFiles/OPTtotalMasses.pyc.npy")
# masses1000 = load("DataFiles/1000totalMasses.pyc.npy")
# pressures = load("DataFiles/OPTtotalPressures.pyc.npy")
# times = load("DataFiles/OPTtotalTime.pyc.npy")

# relative_tolerance, totalTime, totalRadii, totalMasses, totalPressures = loadtxt("optimize_grid_data.txt", dtype='float', usecols=(0,1,2,3,4), unpack=True)

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
