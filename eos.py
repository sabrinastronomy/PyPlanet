from burnman import *
from numpy import *
from scipy import interpolate

def create_function(pressures, dependent_var):  # TODO incorporate into object
    # Interpolates from data sets with density in the first column
    # OR C_p in the first column
    # OR other variable
    # and pressures in the second column

    # Extrapolates data outside of range of data

    # random variable names x = press, y = density/C_p

    return interpolate.interp1d(pressures, dependent_var, bounds_error=False, fill_value="extrapolate")



# Two layer

class Layer:
    def __init__(self, name, material, composition, temp_profile):
        self.name = name
        self.material = material
        self.composition = composition
        self.temp_profile = temp_profile
        if name == "mantle":
            material.set_composition(composition)
        self.composite = Composite([self.material], [1])


class EoS:  # TODO generalize so you can specify temp profile for each layer
    def __init__(self, p_c, p_cmb, temp_profile, layer1, layer2,
                 anchor_temperature=0):  # TODO MAKE FOR ARBITRARY NUMBER OF LAYERS
        self.p_c = p_c
        self.p_cmb = p_cmb

        self.pressures_mantle = []
        self.pressures_core = []

        self.temperatures_mantle = []
        self.temperatures_core = []

        self.anchor_temperature = anchor_temperature
        self.temp_profile = temp_profile
        self.core = layer1
        self.mantle = layer2
        self.core_eos = ""
        self.mantle_eos = ""
        self.get_temp = ""

        self.core_mm = self.core.composite.molar_mass
        self.mantle_mm = self.mantle.composite.molar_mass


        self.eos_function_generate()


    def eos_function_generate(self):
        p_cmb = self.p_cmb
        p_c = self.p_c
        anchor_temperature = self.anchor_temperature

        if self.temp_profile == "constant":
            self.pressures_mantle = self.pressures_core = linspace(0, self.p_c, 1e4)
            pressures = self.pressures_mantle

            self.temperatures_core = empty(len(pressures))
            self.temperatures_core.fill(self.anchor_temperature)

            self.temperatures_mantle = self.temperatures_core

            temperatures = self.temperatures_mantle

            mantle_rho_data = self.mantle.composite.evaluate(['density'], pressures, temperatures)
            core_rho_data = self.core.composite.evaluate(['density'], pressures, temperatures)

            mantle_cp_data = self.mantle.composite.evaluate(['heat_capacity_p'], pressures, temperatures)
            core_cp_data = self.core.composite.evaluate(['heat_capacity_p'], pressures, temperatures)

            self.mantle_eos = [create_function(pressures, mantle_rho_data), create_function(pressures, self.mantle_mm*mantle_cp_data)]
            self.core_eos = [create_function(pressures, core_rho_data), create_function(pressures, self.core_mm*core_cp_data)]


        else:  # TODO make more generic - currently adiabatic
            self.pressures_mantle = np.linspace(0, p_cmb, 20).tolist()
            self.pressures_core = np.linspace(p_c, p_cmb, 20).tolist()
            self.pressures_core.reverse()


            if p_cmb == 0:
                self.temperatures_mantle = [anchor_temperature]
                self.get_temp = lambda x: self.temperatures_mantle

            else:
                self.temperatures_mantle = geotherm.adiabatic(self.pressures_mantle, anchor_temperature,
                                                         self.mantle.composite)  # first pressure is 0

                mantle_rho_data = self.mantle.composite.evaluate(['density'], self.pressures_mantle, self.temperatures_mantle)

                mantle_cp_data = self.mantle.composite.evaluate(['heat_capacity_p'], self.pressures_mantle, self.temperatures_mantle)

                self.mantle_eos = [create_function(self.pressures_mantle, mantle_rho_data), create_function(self.pressures_mantle, self.mantle_mm*mantle_cp_data)]
                self.get_temp = create_function(self.pressures_mantle, self.temperatures_mantle)

            if p_c == p_cmb:
                self.temperatures_core = [self.temperatures_mantle[-1]]
                self.get_temp = lambda x: self.temperatures_core

            else:
                self.temperatures_core = geotherm.adiabatic(self.pressures_core, self.temperatures_mantle[-1],
                                                           self.core.composite)  # first pressure is p_cmb
                core_rho_data = self.core.composite.evaluate(['density'], self.pressures_core, self.temperatures_core)
                core_cp_data = self.core.composite.evaluate(['heat_capacity_p'], self.pressures_core, self.temperatures_core)


                self.core_eos = [create_function(self.pressures_core, core_rho_data), create_function(self.pressures_core, self.core_mm*core_cp_data)]

                self.get_temp = create_function(self.pressures_core, self.temperatures_core)

            if self.core_eos is not "":
                self.core_eos.append(self.get_temp)

            if self.mantle_eos is not "":
                self.mantle_eos.append(self.get_temp)



# TODO: implement iron class, wasn't working in adiabat eos so try later

