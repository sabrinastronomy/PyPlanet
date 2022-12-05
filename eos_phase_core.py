import os

import numpy as np
from scipy import interpolate
from pynbody.analysis.interpolate import interpolate3d

class CoreEos:
    def __init__(self, type_prof, temperatures, pressures_core=None):
        assert type_prof == "_constant" or type_prof == "_adiabatic_"
        self.eos_loc = "/Users/sabrinaberger/Desktop/EoS/"
        os.chdir(self.eos_loc)  # loading in files from EoS directory
        self.type_prof = type_prof
        self.read_in_tables_interpolate()
        self.pressures = pressures_core
        self.temperatures = temperatures
        self.core_mat = np.empty(np.shape(self.temperatures))
        self.get_densities()

    def read_in_tables_interpolate(self):
        print('Read liquid Fe table')
        rho_Fel = np.loadtxt('rho_Fel_2000GPa_9e3K.txt')
        T_Fel = np.loadtxt('T_Fel_2000GPa_9e3K.txt')
        P_Fel = np.loadtxt('P_Fel_2000GPa_9e3K.txt')

        ### Mixture table is not used but could be implemented later on
        print('Read FeSi mixture table')
        rho_Fea = np.loadtxt('rho_Fe16Si_2000GPa_9e3K.txt')
        T_Fea = np.loadtxt('T_Fe16Si_2000GPa_9e3K.txt')
        P_Fea = np.loadtxt('P_Fe16Si_2000GPa_9e3K.txt')

        print('Read solid Fe table')
        rho_Fes = np.loadtxt('rho_Fes_2000GPa_9e3K.txt')
        T_Fes = np.loadtxt('T_Fes_2000GPa_9e3K.txt')
        P_Fes = np.loadtxt('P_Fes_2000GPa_9e3K.txt')

        Tlg, Plg = np.meshgrid(T_Fel, P_Fel, sparse=True)
        Tag, Pag = np.meshgrid(T_Fea, P_Fea, sparse=True)
        Tsg, Psg = np.meshgrid(T_Fes, P_Fes, sparse=True)

        # looks like these return a density for a given (P, T)
        self.f_rho_Fel = interpolate.interp2d(Tlg, Plg, rho_Fel)
        self.f_rho_Fea = interpolate.interp2d(Tag, Pag, rho_Fea)
        self.f_rho_Fes = interpolate.interp2d(Tsg, Psg, rho_Fes)

    def melting_curve_iron(self, press, curr_T):
        # This is equation (4) in Zhang & Rogers 2022
        T_Fe_o = 1900
        T_melted = T_Fe_o * (press/31.3 + 1)**(1/1.99)
        molten = curr_T > T_melted # will be molten if current T is greater than T_melted
        return molten

    def get_Fel_adiabatic_temp_profile(self, pressures_upper, anchor_temperature_liq_iron):
        # This is the adiabat table for the liquid iron core.
        # The interpolation function requires three input, x (mass fraction of Fe16Si alloy),
        # Tref or T0 (the anchor temperature, antropy temperature, see my paper section 2.3), and pressure.
        loaded_T = np.loadtxt('Fe_adiabat_2000GPa_7e3K.txt')
        load_original_T = loaded_T.reshape(loaded_T.shape[0], loaded_T.shape[1] // 391, 391)

        # Grid for x, Tref/T0 and pressure
        # check these grid to find the range for interpolation.
        x_core_grid = np.loadtxt('Fe_adiabat_xgrid_2000GPa_7e3K.txt')
        Tref_core_grid = np.loadtxt('Fe_adiabat_Tgrid_2000GPa_7e3K.txt')
        P_core_grid = np.loadtxt('Fe_adiabat_Pgrid_2000GPa_7e3K.txt')
        # Tgrid_core_grid = np.loadtxt('Fe_adiabat_P_Tgridgrid_2000GPa_7e3K.txt')

        x_alloy = 0.5  # x_alloy is calculated a x_core/0.16, where x_core is the mass fraction of Si
        T_en = anchor_temperature_liq_iron
        adiabat_array = interpolate3d(np.ones(np.shape(pressures_upper)) * x_alloy,
                                      np.ones(np.shape(pressures_upper)) * T_en,
                                      pressures_upper,
                                      x_core_grid,
                                      Tref_core_grid,
                                      P_core_grid,
                                      load_original_T)
        return adiabat_array

    def get_densities(self):
        vect_melting_curve_iron = np.vectorize(self.melting_curve_iron) # vectorizing melting curve function
        binary_array = vect_melting_curve_iron(self.pressures, self.temperatures) # gets a binary array where True means liquid Fe and False means solid Fe
        # getting molten and non molten masks
        molten_mask = binary_array
        non_molten_mask = np.invert(binary_array)
        if self.type_prof == "_adiabatic_":
            # getting LIQUID iron core adiabatic temperatures for given pressures
            anchor_temperature = np.min(self.temperatures[non_molten_mask]) # getting the anchor temperature at the boundary between liq/sol core
            self.temperatures[molten_mask] = self.get_Fel_adiabatic_temp_profile(self.pressures[molten_mask], anchor_temperature)

        self.densities = np.empty(np.shape(self.pressures))
        molten_densities = self.f_rho_Fel(self.temperatures[molten_mask], self.pressures[molten_mask])
        non_molten_densities = self.f_rho_Fes(self.temperatures[non_molten_mask], self.pressures[non_molten_mask])[0, :]
        if molten_densities.ndim == 2:
            self.densities[molten_mask] = molten_densities[0, :]
        else:
            self.densities[molten_mask] = molten_densities
        if non_molten_densities.ndim == 2:
            self.densities[non_molten_mask] = non_molten_densities[0, :]
        else:
            self.densities[non_molten_mask] = non_molten_densities

        self.Cp = np.full(np.shape(self.densities), 840) # 840 J/K/kg - constant heat capacity in core
        self.core_mat = ["liq_iron" if bin_core else "sol_iron" for bin_core in molten_mask]

