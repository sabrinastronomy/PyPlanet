from burnman import *
from numpy import *
from scipy import interpolate



def create_function(pressures, density):  # TODO incorporate into object
    # Interpolates from data sets with density in the first column and pressure in the second column
    # Extrapolates data outside of range of data

    # random variable names x = press, y = density

    press = np.array(pressures)

    return interpolate.interp1d(press, density, bounds_error=False, fill_value="extrapolate")


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

        self.pressures = []

        self.anchor_temperature = anchor_temperature
        self.temp_profile = temp_profile
        self.core = layer1
        self.mantle = layer2

        if temp_profile == "constant":
            self.pressures = linspace(0, self.p_c, 1e4)
            pressures = self.pressures

            temperatures = empty(len(pressures))
            temperatures.fill(self.anchor_temperature)

            mantle_eos_data = self.mantle.composite.evaluate(['density'], pressures, temperatures)
            core_eos_data = self.core.composite.evaluate(['density'], pressures, temperatures)


        else:  # TODO make more generic - currently adiabatic
            pressures_mantle = np.linspace(0, p_cmb, 20).tolist()
            pressures_core = np.linspace(p_c, p_cmb, 20).tolist()
            pressures_core.reverse()

            temperatures_mantle = geotherm.adiabatic(pressures_mantle, anchor_temperature,
                                                     self.mantle.composite)  # first pressure is 0
            temperatures_core = geotherm.adiabatic(pressures_core, temperatures_mantle[-1],
                                                   self.core.composite)  # first pressure is p_cmb

            mantle_eos_data = self.mantle.composite.evaluate(['density'], pressures_mantle, temperatures_mantle)
            core_eos_data = self.core.composite.evaluate(['density'], pressures_core, temperatures_core)

            core_eos = create_function(pressures_core, core_eos_data)
            mantle_eos = create_function(pressures_mantle, mantle_eos_data)

        self.core_eos = core_eos
        self.mantle_eos = mantle_eos

# TODO: implement iron class, wasn't working in adiabat eos so try later

