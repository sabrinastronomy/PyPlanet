from burnman import *
import numpy as np
from eos import *
from scipy import interpolate

data = "/Users/sabrinaberger/All Research/RockyPlanets/DataFiles"

temps = [1000, 300, 3000]

# TO DO make part of eos class instead

    pcmb = np.load(data + 'p_cmb_grid_adiabatic_{}.pyc.npy'.format(temp))
    mass = np.load(data + 'core_mass_grid_adiabatic_{}.pyc.npy'.format(temp))
    pc = np.load(data + 'p_c_adiabatic_{}.pyc.npy'.format(temp))

    pressures_mantle = np.linspace(0, p_cmb, 20).tolist()
    pressures_core = np.linspace(p_c, p_cmb, 20).tolist()
    pressures_core.reverse()

    mantle_material = minerals.SLB_2011.mg_fe_silicate_perovskite()
    mantle_material.set_composition([0.1, 0.1, 0.8])
    mantle_composite = composite([mantle_material], [1])

    temperatures_mantle = geotherm.adiabatic(pressures_mantle, 3000, mantle_composite)
    mantle_heats = mantle_composite.evaluate(['C_v'], pressures_mantle, temperatures_mantle)

    core_material = minerals.Murakami_2013.fe_perovskite()


    temperatures_core = geotherm.adiabatic(pressures_core, temperatures_mantle[-1],
                                       self.core.composite)  # first pressure is p_cmb

    core_heats = core_composite.evaluate(['C_v'], pressures_core, temperatures_core)
