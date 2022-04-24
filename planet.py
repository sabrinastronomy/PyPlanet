# VARIABLES
# t is radius
# y[0] = mass, y[1] = pressure

from math import *
import numpy as np
from scipy.integrate import ode

# defining constants
G = float("6.67e-11")  # G = [N *(m/kg)^2]


class Planet:
    def __init__(self, planet_eos, t0, dt, relative_tolerance, transpress, save_all=False):

        # Initial parameters
        self.eoss = planet_eos.layer_eos
        # self.mm = [planet_eos.core_mm, planet_eos.mantle_mm]
        self.t0 = t0
        self.dt = dt
        self.relative_tolerance = relative_tolerance
        self.transpress = transpress
        self.nLayers = len(transpress) - 1

        # if save_all = True, then all radius, mass and pressure calculated values will be saved
        # if save_all = False, then only final radius, mass and pressure values calculated will be saved
        self.save_all = save_all
        # initializing arrays to hold all values of radius, pressure and mass throughout the planet

        # Parameters to be calculated

        self.rad = []
        self.press = []
        self.mass = []
        self.u = []

        # initializing arrays to hold total radius, pressure and mass at the end of each layer of the planet
        self.transition_rad_list = np.zeros(shape=(len(transpress),))
        self.transition_press_list = np.zeros(shape=(len(transpress),))
        self.transition_mass_list = np.zeros(shape=(len(transpress),))
        self.transition_u_list = np.zeros(shape=(len(transpress),))



    def all_rho_C_p_T(self, pressure, n, i): # TODO Check this
        # first element of EoS is density(pressure) i = 0
        # second element of EoS is C_p(pressure) i = 1
        # third element of EoS is T(pressure) i = 2
        eoss = self.eoss
        if n == 1 or n == 0:
            # print("n = {}".format(n))
            # print("i = {}".format(i))

            return eoss[0][i](pressure)

        else:
            # print("n = {}".format(n))
            # print("i = {}".format(i))

            return eoss[1][i](pressure)


    def integratePlanet(self):
        print("Planet initialized.")

        # Aliases
        t0 = self.t0
        dt = self.dt
        relative_tolerance = self.relative_tolerance
        transpress = self.transpress
        save_all= self.save_all
        rad = self.rad
        press = self.press
        mass = self.mass

        # Heat capacity
        u = self.u

        # initializing arrays to hold total radius, pressure and mass at the end of each layer of the planet
        transition_rad_list = self.transition_rad_list
        transition_press_list = self.transition_press_list
        transition_mass_list = self.transition_mass_list
        transition_u_list = self.transition_u_list


        # checking that list is sorted and monotonically decreasing
        if sorted(transpress, reverse=True) != transpress:
            raise ValueError("Transition pressure list entered is not sorted.")
            sys.exit(0)
        elif transpress[0] == transpress[self.nLayers]:
            print("%g %g %g" % (transpress[0], transpress[1], transpress[2]))
            raise ValueError("No layers to integrate.")
            sys.exit(0)

        def f(t, y, eos, n):
            """
            :param t: current radius
            :param y: current pressure
            :param eos: equation of state being used for this layer
            :param n: layer number
            :return: [mass, pressure, total energy]
            """
            pressure = y[1]
            next_rho = eos(pressure, n, 0)
            next_C_p = eos(pressure, n, 1)
            next_T = eos(pressure, n, 2)

            mass = 4. * pi * (t ** 2.) * next_rho
            press = -(G * y[0] * next_rho / (t ** 2.))
            thermal_energy = 4. * pi * (t ** 2.) * next_rho * next_C_p * next_T
            grav_energy = 4 * pi * t * next_rho * G * mass
            # print(f"Thermal energy = {thermal_energy}")
            # print(f"grav energy = {grav_energy}")

            tot_energy = thermal_energy - grav_energy
            return [mass, press, tot_energy]



        def RK4(derivatives, y0, t0, dt, transition_pressure, n):
            z = ode(derivatives).set_integrator('dopri5').set_initial_value(y0, t0).set_f_params(self.all_rho_C_p_T, n)
            # mass - y[0]
            # press - y[1]
            # u - y[2]
            # radius - t
            which_loop = 'first loop'

            while z.successful() and z.y[1] - transition_pressure > 0 and which_loop == 'first loop':
                t_next = z.t + dt
                z.integrate(t_next)
                rad.append(np.array(z.t).tolist())
                mass.append(np.array(z.y[0]).tolist())
                press.append(np.array(z.y[1]).tolist())

                # Heat capacity
                u.append(z.y[2].tolist())

                # Enters if-statement if passed transition pressure
                press_b = press[-1]
                if press_b - transition_pressure < 0:
                    press_a = press[-2]
                    which_loop = 'second loop'
                    a = rad[-2]  # last positive pressure
                    b = rad[-1]  # first negative pressure
                    c = 0.5 * (a + b)

                # if z.y[1]) == transition_pressure, will not enter second while loop
                if z.y[1] == transition_pressure:
                    rad[-1] = z.t
                    mass[-1] = z.y[0]
                    press[-1] = z.y[1]
                    u[-1] = z.y[2]

            # Bisection Method
            counter = 0
            while z.successful() and abs(transition_pressure - z.y[1]) > max(transition_pressure * relative_tolerance, 10000) and which_loop == 'second loop':
                counter += 1
                z.integrate(c)
                press_c = z.y[1]
                if press_a > transition_pressure and press_c > transition_pressure:
                    press_a = press_c
                    a = c
                else:
                    press_b = press_c
                    b = c
                c = 0.5 * (a + b)
            print("This layer went through " + repr(counter) + " bisection method steps.")
            # TODO: MAKE SURE PRESSURE RETURNS A POSITIVE VALUE
            rad[-1] = z.t
            mass[-1] = z.y[0]
            press[-1] = z.y[1]
            u[-1] = z.y[2]
            print("rad {}".format(z.t))
            print("mass {}".format(z.y[0]))
            print("press {}".format(z.y[1]))
            print("u {}".format(z.y[2]))



            if not z.successful():
                raise RuntimeError("Integration failed.")
            if save_all:
                return rad[-1], mass[-1], press[-1], u[-1], rad, mass, press, u
            else:
                return rad[-1], mass[-1], press[-1], u[-1]

        i = 1
        first_layer = True
        mass_init = 0
        rad_init = 0
        u_init = 0
        p_c = transpress[0]
        press_init = p_c
        transition_rad_list[0] = rad_init
        transition_mass_list[0] = mass_init
        transition_press_list[0] = press_init
        transition_u_list[0] = u_init

        while i <= self.nLayers:
            print("layer %g" % (i))
            # First tests whether or not to integrate layer
            # This takes care of [p_c, p_c, 0] and [p_c, 0, 0] cases
            if transpress[i] == transpress[i - 1]:  # Don't integrate layer
                transition_rad_list[i] = transition_rad_list[i - 1]
                transition_mass_list[i] = transition_mass_list[i - 1]
                transition_press_list[i] = transition_press_list[i - 1]
                transition_u_list[i] = transition_u_list[i - 1]


            else:  # need to integrate layer
                # set ICs for integration
                if first_layer:
                    # Need to carefully set ICs in first layer
                    # print("first layer %d" % (i))
                    mass_0 = (4 / 3) * pi * pow(t0, 3) * self.all_rho_C_p_T(p_c, i, 0)
                    C_p_0 = self.all_rho_C_p_T(p_c, i, 1)
                    T_0 = self.all_rho_C_p_T(p_c, i, 2)
                    y0 = [mass_0, p_c, mass_0 * C_p_0 * T_0] # HEAT CAPACITY
                    print("first layer %d Initial Radius: %g Initial Mass: %g Central Pressure: %g U: %g" % (i, rad_init, y0[0], y0[1], y0[2]))
                    rad_init = t0
                    # Update detailed profile arrays
                    rad.append(t0)
                    mass.append(y0[0])
                    press.append(y0[1])

                    # heat capacity
                    u.append(y0[2])

                    first_layer = False
                else:
                    # Take values at end of last layer as ICs
                    y0 = [mass_init, press_init, u_init]

                # integrate layer
                if save_all:
                    print("integrate layer %d" % (i))
                    # full lists of all radii, mass, and pressures calculated for each integration by RK4
                    # properties given t0 and y0
                    rad_init, mass_init, press_init, u_init, more_rad, more_mass, more_press, more_u = RK4(f, y0, rad_init, dt,
                                                                                           transpress[i], i)
                    rad.extend(more_rad)
                    mass.extend(more_mass)
                    press.extend(more_press)

                    # heat capacity
                    u.extend(more_u)
                else:
                    print("integrate layer %d" % (i))
                    rad_init, mass_init, press_init, u_init = RK4(f, y0, rad_init, dt, transpress[i], i)
                # Update arrays containing properties of planets at transition
                transition_rad_list[i] = rad_init
                transition_mass_list[i] = mass_init
                transition_press_list[i] = press_init
                transition_u_list[i] = u_init

            i += 1

        print("Planet complete.")

        if save_all:
            return rad, mass, press, u, transition_rad_list, transition_mass_list, transition_press_list
        else:
            return rad_init, mass_init, press_init, u_init, transition_rad_list, transition_mass_list, transition_press_list, transition_u_list
