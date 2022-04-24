"""
Written by Sabrina Berger and Jisheng Zhang
"""
import matplotlib.pyplot as plt
import numpy as np
from burnman import *
from scipy import interpolate


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
    rho = np.reshape(rho, (-1, 1501))  # should match the array lengths found when running testing_pre_adiabat2.py code

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
    print("molten")
    print(temp)
    print(press)
    print("____")
    return f_CV(press, temp)


class Layer:
    def __init__(self, name, material, composition, temp_profile):
        """
        Object for layer within planet
        :param name: name of layer, e.g., layer1, etc.
        :param material: type of Burnman material
        :param composition: composition
        :param temp_profile: type of temperature profile
        """
        self.name = name
        print("Initializing {} layer to generate equation of state...".format(name))
        print(f"material {material}")
        if len(material) > 1:
            # Layers with more than one type of material (most frequent)
            self.material = material
            self.composition = composition
            self.temp_profile = temp_profile
            if name == "layer1": # mantle is called layer1
                self.material[0].set_composition(composition[0])

            self.composite = [Composite([mat], [1]) for mat in self.material]
            self.molar_mass = [comp.molar_mass for comp in self.composite]

        else:
            # Layers with just one type of material
            self.material = material[0]
            self.composition = composition[0]

            self.temp_profile = temp_profile
            if name == "layer1": # mantle is called layer1
                self.material.set_composition(composition)
            self.composite = Composite([self.material], [1])
            self.molar_mass = self.composite.molar_mass
        print("Layer initialized with material: {}, composition: {}, and temperature profile: {}.".format(material,
                                                                                                          composition,
                                                                                                          temp_profile))


class EoS:
    def __init__(self, p_c, p_cmb, temp_profile, layers,
                 anchor_temperature=0):
        """
        Object representing equation of state of layers of planet
        :param p_c: central pressure (Pascals)
        :param p_cmb: core mantle boundary pressure (Pascals)
        :param temp_profile: type of temperature profile0
        :param layers: array of layer objects
        :param anchor_temperature: anchor temperature for adiabat
        """

        # anchor temperature
        self.p_c = p_c
        self.p_cmb = p_cmb

        self.pressures = [] # array of pressures for each layer
        self.temperatures = [] # array of temperatures for each layer

        self.anchor_temperature = anchor_temperature
        self.temp_profile = temp_profile

        self.layers = layers
        self.layer_eos = []  # each element contains the density, heat capacity and temperature as a function of pressure for each corresponding layer
        self.layer_rho_data = []
        self.layer_cp_data = []

        self.get_temp = []  # array of functions of interpolated pressures and temperatures
        self.layer_mm = [layer.molar_mass for layer in self.layers] # getting molar mass for each layer

        if self.temp_profile == "_constant_":  # constant temperatures
            self.constant_eos_function_generate()
        elif self.temp_profile == "_adiabatic_":  # constant temperatures
            self.mantle_adiabatic_eos_function_generate()
        else:
            print("Invalid temperature profile type.")
            exit()
        self.core_eos_function_generate()

    def constant_eos_function_generate(self):
        p_cmb = self.p_cmb
        p_c = self.p_c
        anchor_temperature = self.anchor_temperature
        # adiabatic temperature profile assumed
        self.pressures_other = np.linspace(100, p_cmb, int(1e4))  # pressures in other layers
        self.pressures_core = np.linspace(p_c, p_cmb, int(1e4)).tolist()  # pressures in core
        self.pressures_core.reverse()
        pressures = self.pressures_other  # starting at lower pressure = 0
        for i in enumerate(self.layers):
            if i == 0:
                self.temperatures.append(np.empty(len(pressures)))
            else:
                self.temperatures.append(np.full(len(pressures), self.anchor_temperature))

        temperatures = self.temperatures

        # CREATING ALL LAYER EOS ARRAYS
        if p_cmb == 0:
            ### NO MANTLE, ONLY CORE ###
            self.get_temp.append(lambda x: anchor_temperature)  # returns anchor temperature for every input pressure
            self.layer_eos.append(["", ""])
            return

        layer_rho_data = []
        layer_cp_data = []

        for mm, layer in zip(self.layer_mm, self.layers):
            if type(mm) != list:  # case for one phase mantle, core, etc.
                # print("layer.composite {}".format(layer.composite))
                self.get_temp.append(create_function(pressures, self.temperatures))
                if layer.name == "layer1": # mantle
                    layer_rho_data = layer.composite.evaluate(['density'], pressures, temperatures)[0]
                    layer_cp_data = layer.composite.evaluate(['molar_heat_capacity_v'], pressures, temperatures)[0] / mm
                    self.layer_eos.append([create_function(pressures, layer_rho_data),
                                           create_function(pressures, layer_cp_data),
                                           self.get_temp[-1]]) # all needed interpolations from pressure to density, heat capacity, and temperature
            else:
                """
                case when there is both an upper and lower mantle
                """
                key = 0
                pressures_upper = self.pressures_other[(
                            self.pressures_other < 26e9)]  # Making sure to stay within max upper mantle transition pressure so calculating adiabat doesn't error
                temperatures_upper = np.full(len(pressures_upper), self.anchor_temperature)
                max_olivine_trans_press = upper_lower_mantle(max(temperatures_upper))
                if max(pressures_upper) < max_olivine_trans_press: # all olivine upper mantle
                    pressures_holder = np.asarray(["", pressures_upper], dtype=object)
                    temperatures_holder = self.temperatures = np.asarray(["", temperatures_upper], dtype=object)
                else: # lower/upper mantle
                    for temp, press in zip(temperatures_upper, pressures_upper):
                        # getting mantle upper transition pressure
                        olivine_trans = upper_lower_mantle(temp)
                        if press > olivine_trans:
                            pressures_upper = np.asarray(pressures[:key])
                            temperatures_upper = np.asarray(temperatures_upper[:key])
                            pressures_lower = np.linspace(olivine_trans, p_c, int(1e4))
                            temperatures_lower = geotherm.adiabatic(pressures_lower, temperatures_upper[-1],
                                                                    layer.composite[0])
                            break
                        key += 1
                    pressures_holder = np.asarray([pressures_lower, pressures_upper], dtype=object)
                    temperatures_holder = np.asarray([temperatures_lower, temperatures_upper], dtype=object)

                ### Getting corresponding EoSs below for upper/lower mantle
                if anchor_temperature < 2500:  # case when enstatite in upper mantle and iron perovskite in lower mantle
                    self.pressures = []
                    self.temperatures = []
                    for i, lay in enumerate(layer.composite):
                        if pressures_holder[i] == "":
                            # edge case
                            print("Planet created with one phase olivine mantle.")
                            layer_rho_data = lay.evaluate(['density'], pressures_holder[1], temperatures_holder[1])[0]
                            layer_cp_data = lay.evaluate(['molar_heat_capacity_v'], pressures_holder[1], temperatures_holder[1])[0] / mm[1]
                            self.pressures = pressures_holder[1]
                            self.temperatures = temperatures_holder[1]
                            break
                        layer_rho_data.extend((lay.evaluate(['density'], pressures_holder[i], temperatures_holder[i])[0]))
                        layer_cp_data.extend((lay.evaluate(['molar_heat_capacity_v'], pressures_holder[i], temperatures_holder[i])[0]) / mm[i])
                        self.pressures.extend(pressures_holder[i])
                        self.temperatures.extend(temperatures_holder[i])

                else:  # case when molten silicates in upper mantle and iron perovskite in lower mantle
                    # generate lower mantle equations of state (0th element)
                    layer_rho_data.append((layer[0].evaluate(['density'], pressures_holder[0], temperatures_holder[0])[0]))
                    layer_cp_data.append((layer[0].evaluate(['molar_heat_capacity_v'], pressures_holder[0], temperatures_holder[0])[0]) / mm[i])
                    # generate upper mantle equation of state (1st element)
                    # returns molten silicates 2-D interpolated functions
                    densities = get_molten_silicates_rho(pressures_holder[1], temperatures_holder[1])[0]
                    heat_capacities = get_molten_silicates_cv(pressures_holder[1], temperatures_holder[1])
                    layer_rho_data.append(densities)
                    layer_cp_data.append(heat_capacities)

                    layer_rho_data = np.concatenate((layer_rho_data[0], layer_rho_data[1])) # makes this one dimensional array through concatenation
                    layer_cp_data = np.concatenate((layer_cp_data[0], layer_cp_data[1]))
                    self.pressures = np.concatenate((pressures_holder[0], pressures_holder[1]))
                    self.temperatures = np.concatenate((temperatures_holder[0], temperatures_holder[1]))

                # add to eos layers, get_temp array for layers other than single phase (unless transition pressure never reached)
                self.get_temp.append(create_function(self.pressures, self.temperatures))
                self.layer_eos.append([create_function(self.pressures, layer_rho_data),
                                       create_function(self.pressures, layer_cp_data),
                                       self.get_temp[-1]])

                # else:
                #     pressures = self.pressures = np.linspace(0.01, self.p_c, int(1e4))
                #     print("mm {}".format(mm))
                #     print("layer.composite {}".format(layer.composite))
                #     if layer.name == "layer1":
                #         print("self.p_c")
                #         print("trouble PT")
                #         print(len(pressures))
                #         print(pressures)
                #
                #         print(len(temperatures))
                #         print(temperatures)
                #
                #         print("---------------")
                #         layer_1_temperatures = temperatures[0]
                #         layer_rho_data = layer.composite.evaluate(['density'], pressures, layer_1_temperatures)[0]
                #         layer_cp_data = layer.composite.evaluate(['molar_heat_capacity_v'], pressures, layer_1_temperatures)[
                #                             0] / mm
                #         self.layer_eos.append([create_function(pressures, layer_rho_data),
                #                                create_function(pressures, layer_cp_data)])
                #     plt.close()
                #     plt.scatter(pressures, temperatures[0])
                #     plt.title("Mantle PT")
                #     plt.show()


    def mantle_adiabatic_eos_function_generate(self):
        p_cmb = self.p_cmb
        p_c = self.p_c
        anchor_temperature = self.anchor_temperature
        # adiabatic temperature profile assumed
        self.pressures_other = np.linspace(100, p_cmb, int(1e4))  # pressures in other layers
        self.pressures_core = np.linspace(p_c, p_cmb, int(1e4)).tolist()  # pressures in core
        self.pressures_core.reverse()
        pressures = self.pressures_other  # starting at lower pressure = 0

        # CREATING ALL LAYER EOS ARRAYS
        if p_cmb == 0:
            ### NO MANTLE, ONLY CORE ###
            self.temperatures.append(anchor_temperature)
            self.get_temp.append(lambda x: anchor_temperature)  # returns anchor temperature for every input pressure
            self.layer_eos.append(["", ""])
            return

        layer_rho_data = []
        layer_cp_data = []

        # if anchor_temperature < 2500: # Only entered if it's cool enough to be. If there's a silicate mantle, there's an upper mantle.
        #     adiabatic_temperatures_mantle = geotherm.adiabatic(pressures, anchor_temperature, self.layers[1].composite[0]) #TODO does it make sense to have it checking for 0
        #     olivine_trans_max_pressure = upper_lower_mantle(max(adiabatic_temperatures_mantle))
        #     if max(pressures) < olivine_trans_max_pressure:  # NO upper mantle
        #         lay_perovskite = self.layers[1].composite[0]
        #         """
        #         case when there is only a lower mantle, perovskite
        #         """
        #         print(f"self.layer_mm[1] {self.layer_mm[1]}")
        #         layer_rho_data.append((lay_perovskite.evaluate(['density'], pressures, adiabatic_temperatures_mantle)[0]))
        #         layer_cp_data.append((lay_perovskite.evaluate(['molar_heat_capacity_v'], pressures, adiabatic_temperatures_mantle)[0]) / self.layer_mm[1][0])
        #         self.get_temp.append(create_function(pressures, adiabatic_temperatures_mantle))
        #         self.layer_eos.append([create_function(pressures, layer_rho_data),
        #                                create_function(pressures, layer_cp_data),
        #                                self.get_temp[-1]])
        #         return

        for mm, layer in zip(self.layer_mm, self.layers):
            if type(mm) != list:  # case for one phase mantle, core, etc.
                # print("layer.composite {}".format(layer.composite))
                adiabatic_temperatures = geotherm.adiabatic(pressures, anchor_temperature, layer.composite)
                self.temperatures.append(adiabatic_temperatures)  # first pressure is 0, anchor_temperature is surface temperature
                self.get_temp.append(create_function(pressures, self.temperatures))
                if layer.name == "layer1": # mantle
                    layer_rho_data = layer.composite.evaluate(['density'], pressures, adiabatic_temperatures)[0]
                    layer_cp_data = layer.composite.evaluate(['molar_heat_capacity_v'], pressures, adiabatic_temperatures)[0] / mm
                    self.layer_eos.append([create_function(pressures, layer_rho_data),
                                           create_function(pressures, layer_cp_data),
                                           self.get_temp[-1]]) # all needed interpolations from pressure to density, heat capacity, and temperature
            else:
                """
                case when there is both an upper and lower mantle
                """
                key = 0
                pressures_upper = self.pressures_other[(
                            self.pressures_other < 26e9)]  # Making sure to stay within max upper mantle transition pressure so calculating adiabat doesn't error
                print(f"length of pressures {len(pressures_upper)}")
                print(f"anchor_temperature {anchor_temperature}")

                temperatures_upper = geotherm.adiabatic(pressures_upper, anchor_temperature,
                                                        layer.composite[1])  # we have to do this in reverse so using enstatite (1)
                # temperatures_other = geotherm.adiabatic(self.pressures_other, anchor_temperature,
                #                                         layer.composite[1])  # we have to do this in reverse
                max_olivine_trans_press = upper_lower_mantle(max(temperatures_upper))
                if max(pressures_upper) < max_olivine_trans_press: # all olivine upper mantle
                    pressures_holder = np.asarray(["", pressures_upper], dtype=object)
                    temperatures_holder = self.temperatures = np.asarray(["", temperatures_upper], dtype=object)
                else: # lower/upper mantle
                    for temp, press in zip(temperatures_upper, pressures_upper):
                        # getting mantle upper transition pressure
                        olivine_trans = upper_lower_mantle(temp)
                        if press > olivine_trans:
                            pressures_upper = np.asarray(pressures[:key])
                            temperatures_upper = np.asarray(temperatures_upper[:key])
                            pressures_lower = np.linspace(olivine_trans, p_c, int(1e4))
                            temperatures_lower = geotherm.adiabatic(pressures_lower, temperatures_upper[-1],
                                                                    layer.composite[0])
                            break
                        key += 1
                    pressures_holder = np.asarray([pressures_lower, pressures_upper], dtype=object)
                    temperatures_holder = np.asarray([temperatures_lower, temperatures_upper], dtype=object)

                ### Getting corresponding EoSs below for upper/lower mantle
                if anchor_temperature < 2500:  # case when  enstatite in upper mantle and iron perovskite in lower mantle
                    longer_length = max(len(pressures_holder[0]), len(pressures_holder[1]))
                    self.pressures = []
                    self.temperatures = []
                    for i, lay in enumerate(layer.composite):
                        if pressures_holder[i] == "":
                            # edge case
                            print("Planet created with one phase olivine mantle.")
                            layer_rho_data = lay.evaluate(['density'], pressures_holder[1], temperatures_holder[1])[0]
                            layer_cp_data = lay.evaluate(['molar_heat_capacity_v'], pressures_holder[1], temperatures_holder[1])[0] / mm[1]
                            self.pressures = pressures_holder[1]
                            self.temperatures = temperatures_holder[1]
                            break
                        layer_rho_data.extend((lay.evaluate(['density'], pressures_holder[i], temperatures_holder[i])[0]))
                        layer_cp_data.extend((lay.evaluate(['molar_heat_capacity_v'], pressures_holder[i], temperatures_holder[i])[0]) / mm[i])
                        self.pressures.extend(pressures_holder[i])
                        self.temperatures.extend(temperatures_holder[i])


                else:  # case when molten silicates in upper mantle and iron perovskite in lower mantle
                    # generate lower mantle equations of state (0th element)
                    layer_rho_data.append((layer[0].evaluate(['density'], pressures_holder[0], temperatures_holder[0])[0]))
                    layer_cp_data.append((layer[0].evaluate(['molar_heat_capacity_v'], pressures_holder[0], temperatures_holder[0])[0]) / mm[i])
                    # generate upper mantle equation of state (1st element)
                    # returns molten silicates 2-D interpolated functions
                    densities = get_molten_silicates_rho(pressures_holder[1], temperatures_holder[1])[0]
                    heat_capacities = get_molten_silicates_cv(pressures_holder[1], temperatures_holder[1])
                    layer_rho_data.append(densities)
                    layer_cp_data.append(heat_capacities)

                    layer_rho_data = np.concatenate((layer_rho_data[0], layer_rho_data[1])) # makes this one dimensional array through concatenation
                    layer_cp_data = np.concatenate((layer_cp_data[0], layer_cp_data[1]))
                    self.pressures = np.concatenate((pressures_holder[0], pressures_holder[1]))
                    self.temperatures = np.concatenate((temperatures_holder[0], temperatures_holder[1]))

                # add to eos layers, get_temp array for layers other than single phase (unless transition pressure never reached)
                self.get_temp.append(create_function(self.pressures, self.temperatures))
                self.layer_eos.append([create_function(self.pressures, layer_rho_data),
                                       create_function(self.pressures, layer_cp_data),
                                       self.get_temp[-1]])

    def core_eos_function_generate(self):
        p_cmb = self.p_cmb
        p_c = self.p_c
        anchor_temperature = self.anchor_temperature
        ### CREATING CORE EOSs ###
        if p_c == p_cmb:
            # NO CORE
            if self.temp_profile == "_adiabatic_":
                self.get_temp.insert(0, self.get_temp[-1])
            self.layer_eos.insert(0, ["", ""])

        else:
            if self.temp_profile == "_adiabatic_":
                if len(self.temperatures) > 1:  # there is a mantle
                    first_layer_temps = self.temperatures[-1]  # lower mantle top temperature which is anchor temperature
                else:
                    first_layer_temps = anchor_temperature
                core = self.layers[0]
                core_mm = self.layer_mm[0]
                print("first layer temperatures core adiabatic {}".format(first_layer_temps))
                print("first composite core adiabatic {}".format(core.composite))

                core_temp_adiabatic_pre = geotherm.adiabatic(self.pressures_core, first_layer_temps,
                                                               core.composite) # core temps before putting them into the list of all core tempeartures
                # self.temperatures.insert(0, core_temp_adiabatic_pre)  # first pressure is p_cmb
                core_temperatures = core_temp_adiabatic_pre  # first pressure is p_cmb
            elif self.temp_profile == "_constant_":
                core = self.layers[0]
                core_mm = self.layer_mm[0]
                # self.temperatures = self.temperatures.flatten()
                core_temperatures = np.full(len(self.pressures_core), self.anchor_temperature)
                # self.pressures_core = self.pressures
            # print(f"core_temperatures: {core_temperatures}")
            # print(f"core_pressures: {self.pressures_core}")


            core_rho_data = core.composite.evaluate(['density'], self.pressures_core, core_temperatures)
            core_cp_data = core.composite.evaluate(['molar_heat_capacity_v'], self.pressures_core, core_temperatures)

            core_rho_data = core_rho_data.flatten()
            core_cp_data = core_cp_data.flatten()

            self.layer_eos.insert(0, [create_function(self.pressures_core, core_rho_data),
                                      create_function(self.pressures_core, core_cp_data / core_mm)])
            if self.temp_profile == "_adiabatic_":
                self.get_temp.insert(0, create_function(self.pressures_core, core_temperatures))
            for i, layer in enumerate(
                self.layer_eos):
                # adding the temperature function as third element in each element in EoS
                if self.temp_profile == "_adiabatic_":
                    self.layer_eos[i].append(self.get_temp[i])
                if self.temp_profile == "_constant_":
                    self.layer_eos[i].append(lambda x: anchor_temperature)
        # print("self.layer_eos" + str(self.layer_eos))
        # print("layers" + str(self.layers))
        # print("temp" + str(self.get_temp))
        # print("eos" + str(self.layer_eos))