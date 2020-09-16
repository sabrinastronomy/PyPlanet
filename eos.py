from burnman import *
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import pandas as pd
import xmeos
from xmeos import models
from xmeos import datamod

CONSTS = models.CONSTS

def create_function(pressures, dependent_var):  # TODO incorporate into object
    # Interpolates from data sets with density in the first column
    # OR C_p in the first column
    # OR other variable
    # and pressures in the second column

    # Extrapolates data outside of range of data

    # random variable names x = press, y = density/C_p

    return interpolate.interp1d(pressures, dependent_var, bounds_error=False, fill_value="extrapolate")

def upper_lower_mantle(temp, P0=25.0, T0=800.0, a=-0.0017): # perovskite phase transition P->GPa
    # Credit: Jisheng from Sotin et al 2007
    # The function defines the linear boundary between perovskite and olivine.
    # Temperature input is in K. If the actual temperature is above the temperature
    # calculated using this function, then the material is in perovskite phase, otherwise olivine.
    # Pressure output in Pa
    return (a*(temp-T0) + P0) * 10**9 # converting to Pascals

def func_transition(self, n):
    if self.transpress[n] == "lower mantle transpress": # == checks VALUE, is checks pointer
        transition_pressure = self.upper_lower_mantle()
        return transition_pressure
    else:
        return self.transpress[n]

def non_decreasing(L):
    return all(x<=y for x, y in zip(L, L[1:]))

# n-layers

class Layer:
    def __init__(self, name, material, composition, temp_profile):
        self.name = name
        print("Initializing {} layer to generate equation of state...".format(name))

        if len(material) > 1:
            self.material = material
            self.composition = composition

            self.temp_profile = temp_profile
            if name == "layer1":  # TO DO FIX THIS: saving for debug?
                # for i, mat in enumerate(self.material): # iterating through materials in mantle
                #     try:
                # Setting lower mantle composition (an Fe/Silicate/Mg perovskite composite)
                self.material[0].set_composition(composition[0])

            self.composite = [Composite([mat], [1]) for mat in self.material]
            self.molar_mass = [comp.molar_mass for comp in self.composite]

        else:
            self.material = material[0]
            self.composition = composition[0]

            self.temp_profile = temp_profile
            if name == "layer1": # TO DO FIX THIS: saving for bug?
                self.material.set_composition(composition)
            self.composite = Composite([self.material], [1])
            self.molar_mass = self.composite.molar_mass
        print("Layer initialized with material: {}, composition: {}, and temperature profile: {}.".format(material, composition, temp_profile))

def get_molten_silicates(press, temp):
    print('Read molten silicate density table')
    rho_liquidpv_file=np.loadtxt('/Users/sabrinaberger/PyPlanet/density_liquidpv.txt')
    rho_liquidpv=rho_liquidpv_file[:,2]
    pressure_liquidpv=rho_liquidpv_file[:,1]
    temperature_liquidpv=rho_liquidpv_file[:,0]

    #reshape density to #pressure*temperature table
    rho_liquidpv_table=np.reshape(rho_liquidpv, (-1, 1501)) ## note that 1501 is the length of pressure array you define above in line 111.

    #remesh pressure temperature to get pt table
    pressure_liquidpv_array=np.linspace(min(pressure_liquidpv), max(pressure_liquidpv), 1501)  ## note that the length of pressure_liquidpv_array and temperature_liquidpv_array should match the pressure and temperature arrays you define above in lines 111 and 112
    temperature_liquidpv_array=np.linspace(min(temperature_liquidpv), max(temperature_liquidpv), 627)
    pv,tv=np.meshgrid(pressure_liquidpv_array, temperature_liquidpv_array, sparse=True)

    #2d interp
    f_liquidpv = interpolate.interp2d(pv,tv,rho_liquidpv_table,kind='cubic')
    return f_liquidpv

def get_molten_heatcap(press, temp):

# density at a specific T and P

class EoS:
    def __init__(self, p_c, p_cmb, temp_profile, layers,
                 anchor_temperature=0):

        # anchor temperature ???
        self.p_c = p_c
        self.p_cmb = p_cmb

        self.pressures = []
        self.temperatures = []

        self.anchor_temperature = anchor_temperature
        self.temp_profile = temp_profile

        self.layers = layers
        self.layer_eos = [] # each element contains the density, heat capacity and temperature as a function of pressure for each corresponding layer
        self.layer_rho_data = []
        self.layer_cp_data = []

        self.get_temp = []
        self.layer_mm = [layer.molar_mass for layer in self.layers]

        self.eos_function_generate()


    def eos_function_generate(self):
        p_cmb = self.p_cmb
        p_c = self.p_c
        anchor_temperature = self.anchor_temperature



        if self.temp_profile == "constant": # constant temperatures

            olivine_trans = upper_lower_mantle(anchor_temperature)

            # Returns the pressure of the phase transition to the
            # lower mantle (olivine/ringwoodite)

            pressures = self.pressures = np.linspace(0, self.p_c, 1e4)

            for i in enumerate(self.layers):
                if i == 0:
                    self.temperatures.append(np.empty(len(pressures)))
                else:
                    self.temperatures.append(np.fill(self.anchor_temperature))

            temperatures = self.temperatures
            # lower mantle pressure
            for mm, layer in zip(self.layer_mm, self.layers):
                print("mm {}".format(mm))
                if type(mm) is not float:
                    pressures_lower = np.linspace(olivine_trans, p_c, 1e4)
                    pressures_upper = np.linspace(0, olivine_trans, 1e4)
                    pressures = np.asarray([pressures_lower, pressures_upper])
                    layer_rho_data = []
                    layer_cp_data = []
                    for i, lay in enumerate(layer):
                        if i == len(layer) - 1: # LIQUID SILICATE
                            self.layer_eos.append(f_liquidpv, lambda x: print("not implemented"))

                        else:
                            layer_rho_data.append(lay.composite.evaluate(['density'], pressures[i], temperatures))
                            layer_cp_data.append(lay.composite.evaluate(['heat_capacity_p'], pressures[i], temperatures))
                        pressures = self.pressures = pressures.flatten()
                        layer_rho_data = np.asarray(layer_rho_data).flatten()
                        layer_cp_data = np.asarray(layer_cp_data).flatten()
                        print("layer_rho_data {}".format(layer_rho_data))
                        print("layer_cp_data {}".format(layer_cp_data))

                        self.layer_eos.append([create_function(pressures, layer_rho_data), create_function(pressures, layer_cp_data/mm)])
                else:
                    pressures = self.pressures = np.linspace(0, self.p_c, 1e4)
                    layer_rho_data = layer.composite.evaluate(['density'], pressures, temperatures)
                    layer_cp_data = layer.composite.evaluate(['heat_capacity_p'], pressures, temperatures)
                    self.layer_eos.append([create_function(pressures, layer_rho_data), create_function(pressures, layer_cp_data/mm)])

        else:  # adiabatic temperature profile assumed
            self.pressures_other = np.linspace(100, p_cmb, 1e4) # pressures in other layers
            self.pressures_core = np.linspace(p_c, p_cmb, 1e4).tolist() # pressures in core
            self.pressures_core.reverse()
            pressures = self.pressures_other # starting at lower pressure = 0

            # CREATING ALL LAYER EOSs
            if p_cmb == 0:
                # NO MANTLE
                # CHECK THIS

                self.temperatures.append([anchor_temperature])
                self.get_temp.append(lambda x: anchor_temperature) # returns anchor temperature for every input pressure
                self.layer_eos.append(["", ""])
            else:
                # lower mantle pressure
                print("self.layer_mm {}".format(self.layer_mm))
                for mm, layer in zip(self.layer_mm, self.layers):
                    if type(mm) is float: # case for one phase mantle, core, etc.
                        print("mm {}".format(mm))
                        print("layer.composite {}".format(layer.composite))
                        temperatures = geotherm.adiabatic(pressures, anchor_temperature, layer.composite)

                        self.temperatures.append(temperatures)  # first pressure is 0, anchor_temperature is surface temperature
                        self.get_temp.append(create_function(pressures, temperatures))
                        if layer.name == "layer1":
                            layer_rho_data = layer.composite.evaluate(['density'], pressures, temperatures)[0]
                            layer_cp_data = layer.composite.evaluate(['heat_capacity_p'], pressures, temperatures)[0] / mm
                            self.layer_eos.append([create_function(pressures, layer_rho_data),
                                                   create_function(pressures, layer_cp_data)])
                        plt.close()
                        plt.scatter(pressures, temperatures)
                        plt.title("Mantle PT")
                        plt.show()

                    else:
                        key = 0
                        pressures_upper = self.pressures_other[(self.pressures_other < 26e9)] # Making sure to stay within max upper mantle transition pressure
                        print("max is {:e}".format(np.max(pressures_upper)))
                        print("anchor temperature {}".format(anchor_temperature))
                        temperatures_upper = geotherm.adiabatic(pressures_upper, anchor_temperature, layer.composite[1]) # we have to do this in reverse
                        for temp, press in zip(temperatures_upper, pressures_upper):
                            # getting mantle upper transition pressure
                            olivine_trans = upper_lower_mantle(temp)
                            if press > olivine_trans:
                                pressures_upper = np.asarray(pressures[:key])
                                temperatures_upper = np.asarray(temperatures[:key])
                                pressures_lower = np.linspace(olivine_trans, p_c, 1e4)
                                temperatures_lower = geotherm.adiabatic(pressures_lower, temperatures_upper[-1], layer.composite[0])
                                # else:
                                #     temperatures_lower = geotherm.adiabatic(pressures_lower, anchor_temperature,
                                #                                             layer.composite[0])

                                break
                            key += 1
                        print("olivine_trans {:e}".format(olivine_trans))
                        pressures = self.pressures = np.asarray([pressures_lower, pressures_upper])
                        temperatures = self.temperatures = np.asarray([temperatures_lower, temperatures_upper])

                        layer_rho_data = []
                        layer_cp_data = []
                        print(temperatures[-1])
                        for i, lay in enumerate(layer.composite):
                            if i == len(layer) - 1:  # LIQUID SILICATE
                                self.layer_eos.append(f_liquidpv, lambda x: print("not implemented"))
                            print("max P {:e}".format(np.max(pressures[i])))
                            print("max T {:e}".format(np.max(temperatures[i])))
                            print(lay)
                            layer_rho_data.append((lay.evaluate(['density'], pressures[i], temperatures[i])[0]))
                            layer_cp_data.append((lay.evaluate(['heat_capacity_p'], pressures[i], temperatures[i])[0])/mm[i])
                        layer_rho_data = np.concatenate((layer_rho_data[0], layer_rho_data[1]))
                        layer_cp_data = np.concatenate((layer_cp_data[0], layer_cp_data[1]))
                        pressures = self.pressures = np.concatenate((pressures[0], pressures[1]))
                        temperatures = self.pressures = np.concatenate((temperatures[0], temperatures[1]))

                        self.layer_eos.append([create_function(pressures, layer_rho_data),
                                               create_function(pressures, layer_cp_data)])
                        self.get_temp.append(create_function(pressures, temperatures))
                        # plt.close()
                        # plt.plot(pressures, temperatures)
                        # plt.title("Mantle PT")
                        # plt.show()

                # for mm, layer in zip(self.layer_mm, self.layers):
                #     temperatures = geotherm.adiabatic(pressures, anchor_temperature,
                #                                                   layer.composite)
                #     self.temperatures.append(temperatures)  # first pressure is 0, anchor_temperature is surface temperature
                #     layer_rho_data = layer.composite.evaluate(['density'], pressures, temperatures)
                #     layer_cp_data = layer.composite.evaluate(['heat_capacity_p'], pressures, temperatures)
                #     self.layer_eos.append([create_function(pressures, layer_rho_data),
                #                            create_function(pressures, layer_cp_data / mm)])
                #     self.get_temp.append(create_function(self.pressures_other, self.temperatures))

                # mantle_rho_data = self.mantle.composite.evaluate(['density'], self.pressures_mantle, self.temperatures_mantle)
                #
                # mantle_cp_data = self.mantle.composite.evaluate(['heat_capacity_p'], self.pressures_mantle, self.temperatures_mantle)
                #
                # self.mantle_eos = [create_function(self.pressures_mantle, mantle_rho_data), create_function(self.pressures_mantle, mantle_cp_data/self.mantle_mm)]

            # CREATING CORE EOSs
            if p_c == p_cmb:
                # NO CORE
                print("self.temperatures " + str(self.temperatures))
                print("self.get_temp " + str(self.get_temp))

                # self.temperatures.insert(0, ([self.temperatures[-1]]))
                self.get_temp.insert(0, self.get_temp[-1])
                self.layer_eos.insert(0, ["", ""])

            else:
                first_layer_temps = self.temperatures[0][-1] # lower mantle top temperature which is anchor temperature
                core = self.layers[0]
                core_mm = self.layer_mm[0]
                self.temperatures.insert(0, geotherm.adiabatic(self.pressures_core, first_layer_temps,
                                                           core.composite))  # first pressure is p_cmb
                core_temperatures = self.temperatures[0]

                core_rho_data = core.composite.evaluate(['density'], self.pressures_core, core_temperatures)
                core_cp_data = core.composite.evaluate(['heat_capacity_p'], self.pressures_core, core_temperatures)

                core_rho_data = core_rho_data.flatten()
                core_cp_data = core_cp_data.flatten()

                # print("layer_rho_data {}".format(core_rho_data))
                # print("layer_cp_data {}".format(core_cp_data))
                #

                self.layer_eos.insert(0, [create_function(self.pressures_core, core_rho_data), create_function(self.pressures_core, core_cp_data/core_mm)])
                self.get_temp.insert(0, create_function(self.pressures_core, core_temperatures))

            # if self.core_eos is not "":
            for i, layer in enumerate(self.layer_eos):     # adding the temperature function as third element in each element in EoS
                self.layer_eos[i].append(self.get_temp[i])
            print("layers" + str(self.layers))
            print("temp" + str(self.get_temp))
            print("eos" + str(self.layer_eos))

# TODO: implement iron class, wasn't working in adiabat eos so try later

