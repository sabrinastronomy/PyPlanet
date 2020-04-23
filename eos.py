from burnman import *
import numpy as np
from scipy import interpolate

def create_function(pressures, dependent_var):  # TODO incorporate into object
    # Interpolates from data sets with density in the first column
    # OR C_p in the first column
    # OR other variable
    # and pressures in the second column

    # Extrapolates data outside of range of data

    # random variable names x = press, y = density/C_p

    return interpolate.interp1d(pressures, dependent_var, bounds_error=False, fill_value="extrapolate")


# n-layers

class Layer:
    def __init__(self, name, material, composition, temp_profile):
        self.name = name
        self.material = material
        self.composition = composition
        print(material)
        print(name)
        self.temp_profile = temp_profile
        if name == "layer1" or name == "layer2": # TO DO FIX THIS: saving for bug?
            print("material " + str(material))
            material.set_composition(composition)
        self.composite = Composite([self.material], [1])
        self.molar_mass = self.composite.molar_mass


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
        print(self.layers[1].molar_mass)
        self.layer_mm = [layer.molar_mass for layer in self.layers]

        self.eos_function_generate()


    def eos_function_generate(self):
        p_cmb = self.p_cmb
        p_c = self.p_c
        anchor_temperature = self.anchor_temperature



        if self.temp_profile == "constant":
            pressures = self.pressures = np.linspace(0, self.p_c, 1e4)

            for i in enumerate(self.layers):
                if i == 0:
                    self.temperatures.append(np.empty(len(pressures)))
                else:
                    self.temperatures.append(np.fill(self.anchor_temperature))

            temperatures = self.temperatures
            # lower mantel pressure
            for mm, layer in zip(self.layer_mm, self.layers):
                layer_rho_data = layer.composite.evaluate(['density'], pressures, temperatures)
                layer_cp_data = layer.composite.evaluate(['heat_capacity_p'], pressures, temperatures)
                self.layer_eos.append([create_function(pressures, layer_rho_data), create_function(pressures, layer_cp_data/mm)])

        else:  # adiabatic temperature profile assumed
            self.pressures_other = np.linspace(0, p_cmb, 20).tolist() # pressures in other layers
            self.pressures_core = np.linspace(p_c, p_cmb, 20).tolist() # pressures in core
            self.pressures_core.reverse()
            pressures = self.pressures_other

            # CREATING ALL LAYER EOSs
            if p_cmb == 0:
                # NO MANTLE
                # CHECK THIS
                print("self.temperatures " + str(self.temperatures))

                self.temperatures.append([anchor_temperature])
                self.get_temp.append(lambda x: anchor_temperature) # returns anchor temperature for every input pressure
                self.layer_eos.append(["", ""])
            else:
                for mm, layer in zip(self.layer_mm, self.layers):
                    temperatures = geotherm.adiabatic(pressures, anchor_temperature,
                                                                  layer.composite)
                    self.temperatures.append(temperatures)  # first pressure is 0, anchor_temperature is surface temperature
                    layer_rho_data = layer.composite.evaluate(['density'], pressures, temperatures)
                    layer_cp_data = layer.composite.evaluate(['heat_capacity_p'], pressures, temperatures)
                    self.layer_eos.append([create_function(pressures, layer_rho_data),
                                           create_function(pressures, layer_cp_data / mm)])
                    self.get_temp.append(create_function(self.pressures_other, self.temperatures))

                # mantle_rho_data = self.mantle.composite.evaluate(['density'], self.pressures_mantle, self.temperatures_mantle)
                #
                # mantle_cp_data = self.mantle.composite.evaluate(['heat_capacity_p'], self.pressures_mantle, self.temperatures_mantle)
                #
                # self.mantle_eos = [create_function(self.pressures_mantle, mantle_rho_data), create_function(self.pressures_mantle, mantle_cp_data/self.mantle_mm)]
            print(self.temperatures)
            # CREATING CORE EOSs
            if p_c == p_cmb:
                # NO CORE
                print("self.temperatures " + str(self.temperatures))
                # self.temperatures.insert(0, ([self.temperatures[-1]]))
                self.get_temp.insert(0, self.get_temp[-1])

            else:
                first_layer_temps = self.temperatures[0][-1] # lower mantle top temperature which is anchor temperature
                core = self.layers[0]
                core_mm = self.layer_mm[0]
                self.temperatures.insert(0, geotherm.adiabatic(self.pressures_core, first_layer_temps,
                                                           core.composite))  # first pressure is p_cmb
                core_temperatures = self.temperatures[0]

                core_rho_data = core.composite.evaluate(['density'], self.pressures_core, core_temperatures)
                core_cp_data = core.composite.evaluate(['heat_capacity_p'], self.pressures_core, core_temperatures)


                self.layer_eos.insert(0, [create_function(self.pressures_core, core_rho_data), create_function(self.pressures_core, core_cp_data/core_mm)])
                self.get_temp.insert(0, create_function(self.pressures_core, core_temperatures))

            # if self.core_eos is not "":
            for i, layer in enumerate(self.layer_eos):     # adding the temperature function as third element in each element in EoS

                self.layer_eos[i].append(self.get_temp[i])
            print("layers" + str(self.layers))
            print("temp" + str(self.get_temp))
            print("eos" + str(self.layer_eos))

# TODO: implement iron class, wasn't working in adiabat eos so try later

