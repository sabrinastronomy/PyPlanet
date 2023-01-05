"""
This file contains the class Layer and EoS. EoSs are generated prior to planet integration.
Written by Sabrina Berger and Jisheng Zhang
"""
import matplotlib.pyplot as plt
import numpy as np
from burnman import *
from scipy import interpolate
import eos_phase_mantle
import eos_phase_core

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

##############################################################################################################

def helper_func_get_closest(press, pressures, materials):
    ind_closest = min(range(len(pressures)), key=lambda i: abs(pressures[i] - press)) # stole from https://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value
    return materials[ind_closest]

class BurnmanLayer:
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
    def __init__(self, p_c, p_cmb, temp_profile, layers, entropy):
        """
        Object representing equation of state of layers of planet
        :param p_c: central pressure (Pascals)
        :param p_cmb: core mantle boundary pressure (Pascals)
        :param temp_profile: type of temperature profile0
        :param layers: array of layer objects
        :param entropy: entropy corresponding to adiabatic profile

        ***Only parameter used in integration of planet is layer_eos***
        ***The layer_eos list contains the density, heat capacity and temperature as a function of
        pressure for each corresponding layer.***
        """
        assert temp_profile == "_constant" or temp_profile == "_adiabatic_"

        # anchor temperature
        self.p_c = p_c
        self.p_cmb = p_cmb
        self.anchor_temperature = "" # not yet defined until mantle has been generated

        self.pressures_core = np.linspace(self.p_cmb, self.p_c, int(1e4))  # pressures in core
        self.pressures_mantle = np.linspace(100, self.p_cmb, int(1e4))  # pressures in mantle layers

        self.temperatures = []  # array of temperatures for each layer

        self.entropy = entropy # this is the entropy of the planet (constant and adiabatic)
        self.temp_profile = temp_profile

        self.layers = layers # all layers, containing multiple instances of Layer
        self.layer_eos = []  # each element contains the density, heat capacity and temperature as a function of
        # pressure for each corresponding layer
        self.layer_rho_data = []
        self.layer_cp_data = []

        self.get_temp = []  # array of functions of interpolated pressures and temperatures
        self.layer_mm = [layer.molar_mass for layer in self.layers]  # getting molar mass for each layer
        self.layers_mat_dict = {"mantle": None, "core": None}

        # Create mantle EoS
        if self.temp_profile == "_constant_":  # constant temperatures
            print("Not working")
            exit()
            # self.temperatures = np.full(int(1e4), self.anchor_temperature) # creating a constant temperature array
            # self.mantle_constant_eos_function_generate()
        elif self.temp_profile == "_adiabatic_":  # adiabatic temperatures
            self.mantle_adiabatic_eos_function_generate()

        # Create core EoS
        self.core_eos_function_generate()
        if self.layers_mat_dict["core"] != None and self.layers_mat_dict["mantle"] != None:
            self.pressures_concatenate = np.concatenate((self.layers_mat_dict["core"][0], self.layers_mat_dict["mantle"][0]))
            self.materials_concatenate = np.concatenate((self.layers_mat_dict["core"][1], self.layers_mat_dict["mantle"][1]))
        elif self.layers_mat_dict["core"] == None:
            self.pressures_concatenate = self.layers_mat_dict["mantle"][0]
            self.materials_concatenate = self.layers_mat_dict["mantle"][1]
        elif self.layers_mat_dict["mantle"] == None:
            self.pressures_concatenate = self.layers_mat_dict["core"][0]
            self.materials_concatenate = self.layers_mat_dict["core"][1]
    # def mantle_constant_eos_function_generate(self):
    #     molten = self.anchor_temperature > 2500 # checking whether planet is fully molten (True) or differentiated into enstate and perovskite
    #     # these are the layers and molar mass for mantle
    #     mantle_layer_composite = self.layers[1].composite # note composite is a list
    #     mantle_mm = self.layer_mm[1]
    #
    #     if self.p_cmb == 0: ### edge case, NO MANTLE, ONLY CORE ###
    #
    #         self.get_temp.append(lambda x: self.anchor_temperature)  #  returns anchor temperature for every input pressure
    #         self.layer_eos.append(["", "", ""]) # return empty string place holders in the layer eos list
    #         return
    #
    #     # creating placeholder lists for density and heat capacity
    #     layer_rho_data = []
    #     layer_cp_data = []
    #
    #     ### Getting corresponding EoSs below for upper/lower mantle
    #     if not molten:  # case when enstatite in upper mantle and iron perovskite in lower mantle
    #
    #         ### getting upper mantle pressures and temperatures
    #         pressures_upper = self.pressures_mantle[(
    #                     self.pressures_mantle < 26e9)]  # Getting upper mantle pressures, making sure to stay within max upper mantle transition pressure
    #         temperatures_upper = np.full(len(pressures_upper), self.anchor_temperature)
    #
    #         ### getting upper and lower mantle transition pressure, pressures below this will be olivine
    #         max_olivine_trans_press = upper_lower_mantle(
    #             max(temperatures_upper))  # this will just be the constant temperature of the planet
    #
    #         if max(pressures_upper) < max_olivine_trans_press:  # all olivine upper mantle
    #             pressures_lower = ""
    #             temperatures_lower = ""
    #         else:  # lower/upper mantle
    #             ### Determining the array element closest to the max_olivine_trans_press
    #             absolute_val_array = np.abs(pressures_upper - max_olivine_trans_press)
    #             smallest_difference_index = absolute_val_array.argmin()
    #
    #             ### separating pressures and temperatures for upper and lower mantle
    #             pressures_upper = np.asarray(pressures_upper[:smallest_difference_index])
    #             temperatures_upper = np.full(len(pressures_upper), self.anchor_temperature)
    #             pressures_lower = np.linspace(max_olivine_trans_press, self.p_cmb, int(1e4))
    #             temperatures_lower = np.full(len(pressures_lower), self.anchor_temperature)
    #
    #         ### pulling out lower and upper mantle
    #         # (note sometimes entire mantle will be olivine if pressures are low enough)
    #         lower_mantle_perovskite_layer = mantle_layer_composite[0]
    #         lower_mantle_perovskite_mm = mantle_mm[0]
    #
    #         upper_mantle_olivine_layer = mantle_layer_composite[1]
    #         upper_mantle_olivine_mm = mantle_mm[1]
    #
    #         if pressures_lower == "": # edge case when pressures never high enough to form silicate perovskite mantle
    #             print("Planet created with one phase olivine mantle so pressures never high enough to form silicate perovskite mantle.")
    #             layer_rho_data = upper_mantle_olivine_layer.evaluate(['density'], pressures_upper, temperatures_upper)[0]
    #             layer_cp_data = upper_mantle_olivine_layer.evaluate(['molar_heat_capacity_v'], pressures_upper, temperatures_upper)[0] / upper_mantle_olivine_mm
    #         else: # creating silicate perovskite and then olivine mantle
    #             layer_rho_data.append((lower_mantle_perovskite_layer.evaluate(['density'], pressures_lower, temperatures_lower)[0]))
    #             layer_cp_data.append((lower_mantle_perovskite_layer.evaluate(['molar_heat_capacity_v'], pressures_lower, temperatures_lower)[0]) / lower_mantle_perovskite_mm)
    #
    #             layer_rho_data.append((upper_mantle_olivine_layer.evaluate(['density'], pressures_upper, temperatures_upper)[0]))
    #             layer_cp_data.append((upper_mantle_olivine_layer.evaluate(['molar_heat_capacity_v'], pressures_upper, temperatures_upper)[0]) / upper_mantle_olivine_mm)
    #
    #             ### Makes this one dimensional array through concatenation
    #             layer_rho_data = np.concatenate((layer_rho_data[0], layer_rho_data[1]))
    #             layer_cp_data = np.concatenate((layer_cp_data[0], layer_cp_data[1]))
    #
    #         ### if edge case (all enstatite so no lower mantle), deleting the string
    #         if pressures_lower == "":
    #             self.pressures_mantle = pressures_upper
    #             self.temperatures_mantle = temperatures_upper
    #         ### generating array of mantle pressures and temperatures
    #         else:
    #             self.pressures_mantle = np.concatenate((pressures_lower, pressures_upper))
    #             self.temperatures_mantle = np.concatenate((temperatures_lower, temperatures_upper))
    #
    #     else:  # hotter planets, case when only molten silicates in mantle
    #         self.temperatures_mantle = np.full(len(self.pressures_mantle), self.anchor_temperature) # Getting upper mantle pressures, making sure to stay within max upper mantle transition pressure
    #
    #         # returns molten silicates 2-D interpolated functions
    #         densities = get_molten_silicates_rho(self.pressures_mantle, self.temperatures_mantle)[0]
    #         heat_capacities = get_molten_silicates_cv(self.pressures_mantle, self.temperatures_mantle)[0] / mantle_mm # converting to C_P, note mantle_mm is that of silicate perovskite
    #         layer_rho_data = densities
    #         layer_cp_data = heat_capacities
    #
    #     ### add to eos layers for entire mantle, get_temp array for layers other than single phase (unless transition pressure never reached)
    #     # only reaches here if there IS a mantle
    #     self.layer_eos.append([create_function(self.pressures_mantle, layer_rho_data), create_function(self.pressures_mantle, layer_cp_data), create_function(self.pressures_mantle, self.temperatures_mantle)])

    def mantle_adiabatic_eos_function_generate(self):
        if self.p_cmb == 0: ### edge case, NO MANTLE, ONLY CORE ##
            self.temperatures_mantle = [] # empty placeholder for temperature
            self.pressures_mantle = [] # empty placeholder for pressure
            self.get_temp.append(lambda x: self.anchor_temperature)  #  returns anchor temperature for every input pressure
            self.layer_eos.append(["", "", ""]) # return empty string place holders in the layer eos list
            return
        # creating placeholder lists for density and heat capacity
        mantle_holder = eos_phase_mantle.AdiabaticMantleEOS(S_setting=self.entropy, P_mantle_planet=self.pressures_mantle, need_plot=False)
        layer_rho_data = mantle_holder.mantle_global_rho
        layer_cp_data = mantle_holder.mantle_global_Cp
        self.temperatures_mantle = mantle_holder.mantle_global_T
        self.mantle_mat = mantle_holder.mantle_global_mat
        self.anchor_temperature = np.max(self.temperatures_mantle)



        ### add to eos layers for entire mantle, get_temp array for layers other than single phase (unless transition pressure never reached)
        # only reaches here if there IS a mantle
        self.layer_eos.append([create_function(self.pressures_mantle, layer_rho_data),
                               create_function(self.pressures_mantle, layer_cp_data),
                               create_function(self.pressures_mantle, self.temperatures_mantle)])

        self.layers_mat_dict["mantle"] = (self.pressures_mantle, self.mantle_mat)

    def core_eos_function_generate(self):
        """
        Breakdown of where we're getting information for the core
        Adiabatic:
        pressures: assigned above starting from central pressure and going to CMB pressure
        temperatures:
        solid iron 1) We get the first round from Burnman using an adiabatic temperature profile for a solid iron core.
        liq iron 2) In eos_phase_core, we replace what should be the cooler temperatures with adiabatic temperatures from Jisheng's code
        densities: All from Jisheng's code
        heat capacity: constant set at 840

        Constant: the same except we set pressures and temperatures explicitly
        """
        # Getting the solid Fe core layer from Burnamn
        core_sol = self.layers[0]
        ### CREATING CORE EOSs ###
        if self.p_c == self.p_cmb: # NO CORE
            ### inserting empty array for core layer_eos
            self.layer_eos.insert(0, ["", "", ""])
            return
        else:
            if self.temp_profile == "_adiabatic_":
                # grabbing max temperature from mantle to use as anchor temperature for the core
                if len(self.temperatures_mantle) > 1:  # only enters if there is a mantle
                    first_layer_temps = np.max(self.temperatures_mantle)  # lower mantle max temperature is the first pressure in the list of core pressures
                    print(f"last anchor temperature of core is {first_layer_temps}.")
                    np.save("anchor_temp.npy", first_layer_temps)
                else:
                    ### NO MANTLE
                    first_layer_temps = np.load("anchor_temp.npy")
                core_temp_adiabatic_pre = geotherm.adiabatic(self.pressures_core, first_layer_temps, core_sol.composite)  # core temps before putting them into the list of all core tempeartures
                self.core_temperatures = core_temp_adiabatic_pre  # first temperature corresponds to cmb
            elif self.temp_profile == "_constant_":
                self.core_temperatures = np.full(len(self.pressures_core), self.anchor_temperature)

            ### Below, the liquid iron core portion is found using Jisheng's tables
            core_phases = eos_phase_core.CoreEos(self.temp_profile, self.core_temperatures, pressures_core=self.pressures_core)

            ### Getting core density and heat capacity data
            core_rho_data = core_phases.densities
            core_cp_data = core_phases.Cp
            # we need to update the pressures with those from eos_phase_core (only different for adiabatic temp. profiles)
            self.pressures_core = core_phases.pressures
            self.temperatures = core_phases.temperatures
            self.core_mat = core_phases.core_mat

            ### Inserting empty array for core layer_eos
            if self.temp_profile == "_constant_":
                self.layer_eos.insert(0, [create_function(self.pressures_core, core_rho_data),
                                          create_function(self.pressures_core, core_cp_data),
                                          lambda x: self.anchor_temperature])
            elif self.temp_profile == "_adiabatic_":
                self.layer_eos.insert(0, [create_function(self.pressures_core, core_rho_data),
                                          create_function(self.pressures_core, core_cp_data),
                                          create_function(self.pressures_core, self.core_temperatures)])


            self.layers_mat_dict["core"] = (self.pressures_core, self.core_mat)
