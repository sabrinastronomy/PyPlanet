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
    Credit: Jisheng Zhang (UChicago) from Sotin et al 2007
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
    Get molten silicates EoS from Jisheng Zhang's (UChicago) code's output data file
    :param press: pressures to interpolate to (Pascals)
    :param temp: temperatures to interpolate to (K)
    :return: densities for molten silicate layer from given input params
    """
    print('Read molten silicate density table')
    silicate_data = np.loadtxt('/Users/sabrinaberger/RockyPlanets/Jisheng/density_liquidpv.txt')
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

        if len(material) > 1:
            # Layers with more than one type of material
            self.material = material
            self.composition = composition

            self.temp_profile = temp_profile
            if name == "layer1":
                self.material[0].set_composition(composition[0])

            self.composite = [Composite([mat], [1]) for mat in self.material]
            self.molar_mass = [comp.molar_mass for comp in self.composite]

        else:
            # Layers with just one type of material
            self.material = material[0]
            self.composition = composition[0]

            self.temp_profile = temp_profile
            if name == "layer1":
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

        self.pressures = []
        self.temperatures = []

        self.anchor_temperature = anchor_temperature
        self.temp_profile = temp_profile

        self.layers = layers
        self.layer_eos = []  # each element contains the density, heat capacity and temperature as a function of pressure for each corresponding layer
        self.layer_rho_data = []
        self.layer_cp_data = []

        self.get_temp = []
        self.layer_mm = [layer.molar_mass for layer in self.layers]

        if self.temp_profile == "_constant_":  # constant temperatures
            self.constant_eos_function_generate()
        elif self.temp_profile == "_adiabatic_":  # constant temperatures
            self.adiabatic_eos_function_generate()
        else:
            print("Invalid temperature profile type.")
            exit()
        self.core_eos_function_generate()

    def constant_eos_function_generate(self):
        p_cmb = self.p_cmb
        p_c = self.p_c
        anchor_temperature = self.anchor_temperature
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
                    if i == len(layer) - 1:  # LIQUID SILICATE
                        # returns molten silicates 2-D interpolated functions
                        densities = get_molten_silicates_rho(pressures[i], temperatures[i])[0]
                        # TODO: heat_capacities = get_molten_silicates_cv(pressures[i], temperatures[i])
                        heat_capacities = np.full(len(densities), 1)  # TODO: wrong!
                        layer_rho_data.append(densities)
                        layer_cp_data.append(heat_capacities)
                    else:
                        layer_rho_data.append(lay.composite.evaluate(['density'], pressures[i], temperatures))
                        layer_cp_data.append(
                            lay.composite.evaluate(['heat_capacity_v'], pressures[i], temperatures))
                    pressures = self.pressures = pressures.flatten()
                    layer_rho_data = np.asarray(layer_rho_data).flatten()
                    layer_cp_data = np.asarray(layer_cp_data).flatten()
                    self.layer_eos.append([create_function(pressures, layer_rho_data),
                                           create_function(pressures, layer_cp_data / mm)])
            else:
                pressures = self.pressures = np.linspace(0, self.p_c, 1e4)
                layer_rho_data = layer.composite.evaluate(['density'], pressures, temperatures)
                layer_cp_data = layer.composite.evaluate(['heat_capacity_v'], pressures, temperatures)
                self.layer_eos.append(
                    [create_function(pressures, layer_rho_data), create_function(pressures, layer_cp_data / mm)])

    def adiabatic_eos_function_generate(self):
        p_cmb = self.p_cmb
        p_c = self.p_c
        anchor_temperature = self.anchor_temperature
        # adiabatic temperature profile assumed
        self.pressures_other = np.linspace(100, p_cmb, 1e4)  # pressures in other layers
        self.pressures_core = np.linspace(p_c, p_cmb, 1e4).tolist()  # pressures in core
        self.pressures_core.reverse()
        pressures = self.pressures_other  # starting at lower pressure = 0
        # CREATING ALL LAYER EOSs
        if p_cmb == 0:
            # NO MANTLE
            # TODO: CHECK THIS
            self.temperatures.append([anchor_temperature])
            self.get_temp.append(
                lambda x: anchor_temperature)  # returns anchor temperature for every input pressure
            self.layer_eos.append(["", ""])
        else:
            # lower mantle pressure
            print("self.layer_mm {}".format(self.layer_mm))
            for mm, layer in zip(self.layer_mm, self.layers):
                if type(mm) is float:  # case for one phase mantle, core, etc.
                    print("mm {}".format(mm))
                    print("layer.composite {}".format(layer.composite))
                    temperatures = geotherm.adiabatic(pressures, anchor_temperature, layer.composite)

                    self.temperatures.append(
                        temperatures)  # first pressure is 0, anchor_temperature is surface temperature
                    self.get_temp.append(create_function(pressures, temperatures))
                    if layer.name == "layer1":
                        layer_rho_data = layer.composite.evaluate(['density'], pressures, temperatures)[0]
                        layer_cp_data = layer.composite.evaluate(['heat_capacity_v'], pressures, temperatures)[
                                            0] / mm
                        self.layer_eos.append([create_function(pressures, layer_rho_data),
                                               create_function(pressures, layer_cp_data)])
                    plt.close()
                    plt.scatter(pressures, temperatures)
                    plt.title("Mantle PT")
                    plt.show()

                else:
                    key = 0
                    pressures_upper = self.pressures_other[(
                                self.pressures_other < 26e9)]  # Making sure to stay within max upper mantle transition pressure
                    print("max is {:e}".format(np.max(pressures_upper)))
                    print("anchor temperature {}".format(anchor_temperature))
                    temperatures_upper = geotherm.adiabatic(pressures_upper, anchor_temperature,
                                                            layer.composite[1])  # we have to do this in reverse

                    ### Getting correct temperatures and pressures below
                    for temp, press in zip(temperatures_upper, pressures_upper):
                        # getting mantle upper transition pressure
                        olivine_trans = upper_lower_mantle(temp)
                        if press > olivine_trans:
                            pressures_upper = np.asarray(pressures[:key])
                            temperatures_upper = np.asarray(temperatures[:key])
                            pressures_lower = np.linspace(olivine_trans, p_c, 1e4)
                            temperatures_lower = geotherm.adiabatic(pressures_lower, temperatures_upper[-1],
                                                                    layer.composite[0])
                            break
                        key += 1

                    print("olivine_trans {:e}".format(olivine_trans))
                    pressures = self.pressures = np.asarray([pressures_lower, pressures_upper])
                    temperatures = self.temperatures = np.asarray([temperatures_lower, temperatures_upper])

                    layer_rho_data = []
                    layer_cp_data = []
                    print(temperatures[-1])

                    ### Getting corresponding EoSs below
                    for i, lay in enumerate(layer.composite):
                        print(layer)
                        if i == len(layer.composite) - 1:  # LIQUID SILICATE
                            # returns molten silicates 2-D interpolated functions
                            densities = get_molten_silicates_rho(pressures[i], temperatures[i])[0]
                            print("densities {}".format(densities))
                            # TODO: heat_capacities = get_molten_silicates_cv(pressures[i], temperatures[i])
                            heat_capacities = np.full(len(densities), 1) # TODO: wrong!
                            layer_rho_data.append(densities)
                            layer_cp_data.append(heat_capacities)
                        else:
                            layer_rho_data.append((lay.evaluate(['density'], pressures[i], temperatures[i])[0]))
                            layer_cp_data.append(
                                (lay.evaluate(['heat_capacity_v'], pressures[i], temperatures[i])[0]) / mm[i])
                    layer_rho_data = np.concatenate((layer_rho_data[0], layer_rho_data[1]))
                    layer_cp_data = np.concatenate((layer_cp_data[0], layer_cp_data[1]))
                    pressures = self.pressures = np.concatenate((pressures[0], pressures[1]))
                    temperatures = self.pressures = np.concatenate((temperatures[0], temperatures[1]))

                    self.layer_eos.append([create_function(pressures, layer_rho_data),
                                           create_function(pressures, layer_cp_data)])
                    self.get_temp.append(create_function(pressures, temperatures))

    def core_eos_function_generate(self):
        ### CREATING CORE EOSs ###
        p_cmb = self.p_cmb
        p_c = self.p_c
        anchor_temperature = self.anchor_temperature
        if p_c == p_cmb:
            # NO CORE
            print("self.temperatures " + str(self.temperatures))
            print("self.get_temp " + str(self.get_temp))

            # self.temperatures.insert(0, ([self.temperatures[-1]]))
            self.get_temp.insert(0, self.get_temp[-1])
            self.layer_eos.insert(0, ["", ""])

        else:
            first_layer_temps = self.temperatures[0][-1]  # lower mantle top temperature which is anchor temperature
            core = self.layers[0]
            core_mm = self.layer_mm[0]
            self.temperatures.insert(0, geotherm.adiabatic(self.pressures_core, first_layer_temps,
                                                           core.composite))  # first pressure is p_cmb
            core_temperatures = self.temperatures[0]

            core_rho_data = core.composite.evaluate(['density'], self.pressures_core, core_temperatures)
            core_cp_data = core.composite.evaluate(['heat_capacity_v'], self.pressures_core, core_temperatures)

            core_rho_data = core_rho_data.flatten()
            core_cp_data = core_cp_data.flatten()

            self.layer_eos.insert(0, [create_function(self.pressures_core, core_rho_data),
                                      create_function(self.pressures_core, core_cp_data / core_mm)])
            self.get_temp.insert(0, create_function(self.pressures_core, core_temperatures))

        for i, layer in enumerate(
                self.layer_eos):  # adding the temperature function as third element in each element in EoS
            self.layer_eos[i].append(self.get_temp[i])
        print("layers" + str(self.layers))
        print("temp" + str(self.get_temp))
        print("eos" + str(self.layer_eos))