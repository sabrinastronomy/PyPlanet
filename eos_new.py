"""
This file contains the class Layer and EoS. EoSs are generated prior to planet integration.
Written by Sabrina Berger and Jisheng Zhang
"""
import matplotlib.pyplot as plt
import numpy as np
from burnman import *
from scipy import interpolate

##############################################################################################################
# helper functions below
def create_function(pressures, dependent_var):
    """
    Interpolates from data sets containing pressure and pressure dependent variable
    Extrapolates data outside of range of data

    :param pressures: pressure
    :param dependent_var: pressure dependent variable, e.g., density, C_p, etc.
    :return: interpolation function
    """
    return interpolate.interp1d(pressures, dependent_var, bounds_error=False, fill_value="extrapolate")

def upper_lower_mantle(temp, P0=25.0, T0=800.0, a=-0.0017):  # perovskite phase transition P->GPa
    """
    Credit: molten_silicates_tables Zhang (UChicago) from Sotin et al 2007
    The function defines the linear boundary between perovskite and olivine.
    Temperature input is in K. If the actual temperature is above the temperature
    calculated using this function, then the material is in perovskite phase, otherwise olivine.
    :param temp: temperature in kelvins
    :param P0: pressure constant in GPa
    :param T0: presure constant in K
    :param a: other constant a
    :return: pressure in Pascals
    """
    return (a * (temp - T0) + P0) * 10 ** 9  # converting to Pascals

def get_molten_silicates_rho(press, temp):
    """
    Get molten silicates EoS from molten_silicates_tables Zhang's (UChicago) code's output data file
    :param press: pressures to interpolate to (Pascals)
    :param temp: temperatures to interpolate to (K)
    :return: densities for molten silicate layer from given input params

    Note: Part of this code was written by molten_silicates_tables Zhang
    """
    print('Read molten silicate density table')
    silicate_data = np.loadtxt('../molten_silicates_tables/density_liquidpv.txt')
    # Set path of density_liquidpv.txt
    temperatures = silicate_data[:, 0]
    pressures = silicate_data[:, 1]
    rho = silicate_data[:, 2]

    # reshape density to pressure & temperature table
    rho = np.reshape(rho, (
        -1, 1501))  # should match the array lengths found when running testing_pre_adiabat2.py code

    # remesh pressure temperature to get pt table
    pressure_array = np.linspace(min(pressures), max(pressures),
                                 1501)  # should match the array lengths found when running testing_pre_adiabat2.py code
    temperature_array = np.linspace(min(temperatures), max(temperatures), 627)
    pv, tv = np.meshgrid(pressure_array, temperature_array, sparse=True)

    # 2-D interpolation
    p_t_rho_interp = interpolate.interp2d(pv, tv, rho, kind='cubic')
    return p_t_rho_interp(press, temp)

def get_molten_silicates_cv(press, temp):
    """
    Get molten silicates heat capacity from molten_silicates_tables Zhang's (UChicago) code's output data file
    :param press: pressures to interpolate to (Pascals)
    :param temp: temperatures to interpolate to (K)
    :return: heat capacities for molten silicate layer from given input params

    Note: Part of this code was written by molten_silicates_tables Zhang
    """
    CV_file = np.loadtxt('../molten_silicates_tables/CV_silicate_melt.txt')
    CV_P = CV_file[0]
    CV_file = CV_file[1:]
    CV_T = np.linspace(1800.0, 3000.0, 1201)
    pv, tv = np.meshgrid(CV_P, CV_T, sparse=True)
    f_CV = interpolate.interp2d(pv, tv, CV_file, kind='cubic')
    return f_CV(press, temp)
##############################################################################################################


class Layer:
    def __init__(self, name_layer, material, composition, temp_profile):
        """
        Object for layer within planet
        :param name_layer: name of layer, e.g., layer1, etc.
        :param material: type of Burnman material
        :param composition: composition
        :param temp_profile: type of temperature profile
        """
        self.name_layer = name_layer
        print("Initializing {} layer to generate equation of state...".format(self.name_layer))
        if len(material) > 1:
            # Usually mantle currently since material length is greater than 1
            # Layers with more than one type of material
            self.material = material
            self.composition = composition
            self.temp_profile = temp_profile
            if self.name_layer == "layer1":  # mantle is called layer1
                self.material[0].set_composition(composition[0])
            self.composite = [Composite([mat], [1]) for mat in self.material]
            self.molar_mass = [comp.molar_mass for comp in self.composite]

        else:
            # Core currently since material length is usually less than 1, i.e., just iron
            # Layers with just one type of material
            self.material = material[0]
            self.composition = composition[0]

            self.temp_profile = temp_profile
            if self.name_layer == "layer1":  # mantle is called layer1
                self.material.set_composition(composition) # have to set compositions for mantle materials but not core
            self.composite = Composite([self.material], [1])
            self.molar_mass = self.composite.molar_mass
        print("{} initialized with material: {}, composition: {}, and temperature profile: {}.".format(self.name_layer, material,
                                                                                                          composition,
                                                                                                          temp_profile))


class EoS:
    def __init__(self, p_c, p_cmb, temp_profile, layers, anchor_temperature=0):
        """
        Object representing equation of state of layers of planet
        :param p_c: central pressure (Pascals)
        :param p_cmb: core mantle boundary pressure (Pascals)
        :param temp_profile: type of temperature profile0
        :param layers: array of layer objects
        :param anchor_temperature: anchor temperature for adiabat

        ***Only parameter used in integration of planet is layer_eos***
        ***The layer_eos list contains the density, heat capacity and temperature as a function of
        pressure for each corresponding layer.***
        """

        # anchor temperature
        self.p_c = p_c
        self.p_cmb = p_cmb

        self.pressures_core = np.linspace(self.p_cmb, self.p_c, int(1e4))  # pressures in core
        self.pressures_mantle = np.linspace(100, self.p_cmb, int(1e4))  # pressures in mantle layers

        self.temperatures = []  # array of temperatures for each layer

        self.anchor_temperature = anchor_temperature # this is the temperature of the whole planet if the temperature profile is constant
        self.temp_profile = temp_profile

        self.layers = layers # all layers, containing multiple instances of Layer
        self.layer_eos = []  # each element contains the density, heat capacity and temperature as a function of
        # pressure for each corresponding layer
        self.layer_rho_data = []
        self.layer_cp_data = []

        self.get_temp = []  # array of functions of interpolated pressures and temperatures
        self.layer_mm = [layer.molar_mass for layer in self.layers]  # getting molar mass for each layer

        # Create mantle EoS
        if self.temp_profile == "_constant_":  # constant temperatures
            self.temperatures = np.full(int(1e4), self.anchor_temperature) # creating a constant temperature array
            self.mantle_constant_eos_function_generate()
        elif self.temp_profile == "_adiabatic_":  # constant temperatures
            self.mantle_adiabatic_eos_function_generate()
        else:
            print("Invalid temperature profile type.")
            exit()

        # Create core EoS
        self.core_eos_function_generate()

    def mantle_constant_eos_function_generate(self):
        # these are the layers and molar mass for mantle
        mantle_layer_composite = self.layers[1].composite # note composite is a list
        mantle_mm = self.layer_mm[1]

        if self.p_cmb == 0: ### edge case, NO MANTLE, ONLY CORE ###

            self.get_temp.append(lambda x: self.anchor_temperature)  #  returns anchor temperature for every input pressure
            self.layer_eos.append(["", "", ""]) # return empty string place holders in the layer eos list
            return

        # creating placeholder lists for density and heat capacity
        layer_rho_data = []
        layer_cp_data = []

        ### Getting corresponding EoSs below for upper/lower mantle
        if self.anchor_temperature < 2500:  # case when enstatite in upper mantle and iron perovskite in lower mantle

            ### getting upper mantle pressures and temperatures
            pressures_upper = self.pressures_mantle[(
                        self.pressures_mantle < 26e9)]  # Getting upper mantle pressures, making sure to stay within max upper mantle transition pressure
            temperatures_upper = np.full(len(pressures_upper), self.anchor_temperature)

            ### getting upper and lower mantle transition pressure, pressures below this will be olivine
            max_olivine_trans_press = upper_lower_mantle(
                max(temperatures_upper))  # this will just be the constant temperature of the planet

            if max(pressures_upper) < max_olivine_trans_press:  # all olivine upper mantle
                pressures_lower = ""
                temperatures_lower = ""
            else:  # lower/upper mantle
                ### Determining the array element closest to the max_olivine_trans_press
                absolute_val_array = np.abs(pressures_upper - max_olivine_trans_press)
                smallest_difference_index = absolute_val_array.argmin()

                ### separating pressures and temperatures for upper and lower mantle
                pressures_upper = np.asarray(pressures_upper[:smallest_difference_index])
                temperatures_upper = np.full(len(pressures_upper), self.anchor_temperature)
                pressures_lower = np.linspace(max_olivine_trans_press, self.p_cmb, int(1e4))
                temperatures_lower = np.full(len(pressures_lower), self.anchor_temperature)

            ### pulling out lower and upper mantle
            # (note sometimes entire mantle will be olivine if pressures are low enough)
            lower_mantle_perovskite_layer = mantle_layer_composite[0]
            lower_mantle_perovskite_mm = mantle_mm[0]

            upper_mantle_olivine_layer = mantle_layer_composite[1]
            upper_mantle_olivine_mm = mantle_mm[1]


            if pressures_lower == "": # edge case when pressures never high enough to form silicate perovskite mantle
                print("Planet created with one phase olivine mantle so pressures never high enough to form silicate perovskite mantle.")
                layer_rho_data = upper_mantle_olivine_layer.evaluate(['density'], pressures_upper, temperatures_upper)[0]
                layer_cp_data = upper_mantle_olivine_layer.evaluate(['molar_heat_capacity_v'], pressures_upper, temperatures_upper)[0] / upper_mantle_olivine_mm
            else: # creating silicate perovskite and then olivine mantle
                layer_rho_data.append((lower_mantle_perovskite_layer.evaluate(['density'], pressures_lower, temperatures_lower)[0]))
                layer_cp_data.append((lower_mantle_perovskite_layer.evaluate(['molar_heat_capacity_v'], pressures_lower, temperatures_lower)[0]) / lower_mantle_perovskite_mm)

                layer_rho_data.append((upper_mantle_olivine_layer.evaluate(['density'], pressures_upper, temperatures_upper)[0]))
                layer_cp_data.append((upper_mantle_olivine_layer.evaluate(['molar_heat_capacity_v'], pressures_upper, temperatures_upper)[0]) / upper_mantle_olivine_mm)

                ### Makes this one dimensional array through concatenation
                layer_rho_data = np.concatenate((layer_rho_data[0], layer_rho_data[1]))
                layer_cp_data = np.concatenate((layer_cp_data[0], layer_cp_data[1]))

            ### if edge case (all olivine so no lower mantle), deleting the string
            if pressures_lower == "":
                self.pressures_mantle = pressures_upper
                self.temperatures_mantle = temperatures_upper
            ### generating array of mantle pressures and temperatures
            else:
                self.pressures_mantle = np.concatenate((pressures_lower, pressures_upper))
                self.temperatures_mantle = np.concatenate((temperatures_lower, temperatures_upper))

        else:  # hotter planets, case when only molten silicates in mantle
            self.temperatures_mantle = np.full(len(self.pressures_mantle), self.anchor_temperature) # Getting upper mantle pressures, making sure to stay within max upper mantle transition pressure

            # returns molten silicates 2-D interpolated functions
            densities = get_molten_silicates_rho(self.pressures_mantle, self.temperatures_mantle)[0]
            heat_capacities = get_molten_silicates_cv(self.pressures_mantle, self.temperatures_mantle)[0] / mantle_mm # converting to C_P, note mantle_mm is that of silicate perovskite
            layer_rho_data = densities
            layer_cp_data = heat_capacities

        ### add to eos layers for entire mantle, get_temp array for layers other than single phase (unless transition pressure never reached)
        # only reaches here if there IS a mantle
        self.layer_eos.append([create_function(self.pressures_mantle, layer_rho_data), create_function(self.pressures_mantle, layer_cp_data), create_function(self.pressures_mantle, self.temperatures_mantle)])

    def mantle_adiabatic_eos_function_generate(self):
        # these are the layers and molar mass for mantle
        mantle_layer_composite = self.layers[1].composite
        mantle_mm = self.layer_mm[1]

        if self.p_cmb == 0: ### edge case, NO MANTLE, ONLY CORE ##
            self.temperatures_mantle = [] # empty placeholder for temperature
            self.pressures_mantle = [] # empty placeholder for pressure
            self.get_temp.append(lambda x: self.anchor_temperature)  #  returns anchor temperature for every input pressure
            self.layer_eos.append(["", "", ""]) # return empty string place holders in the layer eos list
            return

        # creating placeholder lists for density and heat capacity
        layer_rho_data = []
        layer_cp_data = []

        ### Getting corresponding EoSs below for upper/lower mantle
        if self.anchor_temperature < 2500:  # case when enstatite in upper mantle and iron perovskite in lower mantle
            #     #     adiabatic_temperatures_mantle = geotherm.adiabatic(pressures, anchor_temperature, self.layers[1].composite[0]) #TODO does it make sense to have it checking for 0

            ### getting upper mantle pressures and temperatures
            pressures_upper = self.pressures_mantle[(self.pressures_mantle < 26e9)]  # Getting upper mantle pressures, making sure to stay within max upper mantle transition pressure
            temperatures_upper = geotherm.adiabatic(pressures_upper, self.anchor_temperature, mantle_layer_composite[1])  # adiabatic temperatures for olivine upper mantle

            ### getting upper and lower mantle transition pressure, pressures below this will be olivine
            max_olivine_trans_press = upper_lower_mantle(max(temperatures_upper))  # this will just be the constant temperature of the planet

            if max(pressures_upper) < max_olivine_trans_press:  # all olivine upper mantle
                pressures_lower = ""
                temperatures_lower = ""
            else:  # lower/upper mantle
                ### Determining the array element closest to the max_olivine_trans_press
                absolute_val_array = np.abs(pressures_upper - max_olivine_trans_press)
                smallest_difference_index = absolute_val_array.argmin()

                ### separating pressures and temperatures for upper and lower mantle
                pressures_upper = np.asarray(pressures_upper[:smallest_difference_index])
                temperatures_upper = np.asarray(temperatures_upper[:smallest_difference_index])
                pressures_lower = np.linspace(max_olivine_trans_press, self.p_cmb, int(1e4))
                temperatures_lower = geotherm.adiabatic(pressures_lower, self.anchor_temperature,
                                                        mantle_layer_composite[
                                                            0])  # adiabatic temperatures for silicate lower mantle

            ### pulling out lower and upper mantle
            # (note sometimes entire mantle will be olivine if pressures are low enough)
            lower_mantle_perovskite_layer = mantle_layer_composite[0]
            lower_mantle_perovskite_mm = mantle_mm[0]

            upper_mantle_olivine_layer = mantle_layer_composite[1]
            upper_mantle_olivine_mm = mantle_mm[1]


            if pressures_lower == "": # edge case when pressures never high enough to form silicate perovskite mantle
                print("Planet created with one phase olivine mantle so pressures never high enough to form silicate perovskite mantle.")
                layer_rho_data = upper_mantle_olivine_layer.evaluate(['density'], pressures_upper, temperatures_upper)[0]
                layer_cp_data = upper_mantle_olivine_layer.evaluate(['molar_heat_capacity_v'], pressures_upper, temperatures_upper)[0] / upper_mantle_olivine_mm
            else: # creating silicate perovskite and then olivine mantle
                layer_rho_data.append((lower_mantle_perovskite_layer.evaluate(['density'], pressures_lower, temperatures_lower)[0]))
                layer_cp_data.append((lower_mantle_perovskite_layer.evaluate(['molar_heat_capacity_v'], pressures_lower, temperatures_lower)[0]) / lower_mantle_perovskite_mm)

                layer_rho_data.append((upper_mantle_olivine_layer.evaluate(['density'], pressures_upper, temperatures_upper)[0]))
                layer_cp_data.append((upper_mantle_olivine_layer.evaluate(['molar_heat_capacity_v'], pressures_upper, temperatures_upper)[0]) / upper_mantle_olivine_mm)

                ### Makes this one dimensional array through concatenation
                layer_rho_data = np.concatenate((layer_rho_data[0], layer_rho_data[1]))
                layer_cp_data = np.concatenate((layer_cp_data[0], layer_cp_data[1]))

            ### if edge case (all olivine so no lower mantle), deleting the string
            if pressures_lower == "":
                self.pressures_mantle = pressures_upper
                self.temperatures_mantle = temperatures_upper
            ### generating array of mantle pressures and temperatures
            else:
                self.pressures_mantle = np.concatenate((pressures_lower, pressures_upper))
                self.temperatures_mantle = np.concatenate((temperatures_lower, temperatures_upper))

        else:  # hotter planets, case when only molten silicates in mantle
            self.temperatures_mantle = np.add(geotherm.adiabatic(self.pressures_mantle, self.anchor_temperature, mantle_layer_composite), np.full(len(self.pressures_mantle),400)) # using solid silicate to generate adiabatic temperature profile
            # returns molten silicates 2-D interpolated functions
            densities = get_molten_silicates_rho(self.pressures_mantle, self.temperatures_mantle)[0]
            heat_capacities = get_molten_silicates_cv(self.pressures_mantle, self.temperatures_mantle)[0] / mantle_mm # converting to C_P
            layer_rho_data = densities
            layer_cp_data = heat_capacities
        ### add to eos layers for entire mantle, get_temp array for layers other than single phase (unless transition pressure never reached)
        # only reaches here if there IS a mantle
        self.layer_eos.append([create_function(self.pressures_mantle, layer_rho_data), create_function(self.pressures_mantle, layer_cp_data), create_function(self.pressures_mantle, self.temperatures_mantle)])


    def core_eos_function_generate(self):
        ### this is the core layer and molar mass we'll be setting
        core = self.layers[0]
        core_mm = self.layer_mm[0]

        ### CREATING CORE EOSs ###
        if self.p_c == self.p_cmb: # NO CORE
            ### inserting empty array for core layer_eos
            self.layer_eos.insert(0, ["", "", ""])
            return

        else:
            if self.temp_profile == "_adiabatic_":
                if len(self.temperatures_mantle) > 1:  # only enters if there is a mantle
                    first_layer_temps = max(self.temperatures_mantle)  # lower mantle top temperature which is anchor temperature
                else: # enters here if there isn't a mantle (temperatures in mantle list is empty)
                    first_layer_temps = self.anchor_temperature

                core_temp_adiabatic_pre = geotherm.adiabatic(self.pressures_core, first_layer_temps, core.composite)  # core temps before putting them into the list of all core tempeartures
                self.core_temperatures = core_temp_adiabatic_pre  # first temperature corresponds to cmb

            elif self.temp_profile == "_constant_":
                self.core_temperatures = np.full(len(self.pressures_core), self.anchor_temperature)
            else:
                print("Invalid temperature profile type.")
                exit()

            ### Getting core density and heat capacity data
            core_rho_data = core.composite.evaluate(['density'], self.pressures_core, self.core_temperatures)
            core_cp_data = core.composite.evaluate(['molar_heat_capacity_v'], self.pressures_core, self.core_temperatures)

            ### Flattening arrays
            core_rho_data = core_rho_data.flatten()
            core_cp_data = core_cp_data.flatten()

            ### Inserting empty array for core layer_eos
            if self.temp_profile == "_constant_":

                self.layer_eos.insert(0, [create_function(self.pressures_core, core_rho_data),
                                          create_function(self.pressures_core, core_cp_data / core_mm),
                                          lambda x: self.anchor_temperature])
            elif self.temp_profile == "_adiabatic_":
                self.layer_eos.insert(0, [create_function(self.pressures_core, core_rho_data),
                                          create_function(self.pressures_core, core_cp_data / core_mm),
                                          create_function(self.pressures_core, self.core_temperatures)])


