import sys

import numpy as np
import run_perplex


def get_percents(*args):
    FeMg = args[1]
    SiMg = args[2]
    CaMg = args[3]
    AlMg = args[4]
    mol_frac_Fe_mantle = args[5]
    wt_frac_Si_core = args[6]
    wt_frac_O_core = args[7]
    wt_frac_S_core = args[8]

    MgSi = 1. / SiMg
    FeSi = FeMg * MgSi
    CaSi = CaMg * MgSi
    AlSi = AlMg * MgSi

    # constants, atomic masses
    mFe = 55.845
    mMg = 24.306
    mSi = 28.0867
    mO = 15.9994
    mS = 32.0650
    mCa = 40.078
    mAl = 26.981

    # system of equation which when solved gives molar values of elements in mantle and core
    # x = [nFec,nSic,nOc, nSc | ,nFem,nMgm,nSim,nOm, nCam, nAlm]

    # Sum of all masses = 100 g
    b = np.array([0., 0., 0., 0., 0., 0., 0, 0., 0., 100.])

    ############################################################################################

    # Comment here
    # Double checks
    # Return only moles of elements, no molecules
    A = np.array([[0., 0., 0., 0., -1., -1., -2., 1., -1., -1.5],
                  [1., -1. * FeSi, 0., 0., 1., 0., -1. * FeSi, 0., 0., 0.],
                  [0., -1. * MgSi, 0., 0., 0., 1., -1. * MgSi, 0., 0., 0.],
                  [0., -1. * CaSi, 0., 0., 0., 0., -1. * CaSi, 0., 1., 0.],
                  [0., -1. * AlSi, 0., 0., 0., 0., -1. * AlSi, 0., 0., 1.],
                  [mol_frac_Fe_mantle, 0., 0., 0., (mol_frac_Fe_mantle - 1.), 0., 0., 0., 0., 0.],
                  [wt_frac_Si_core * mFe, (wt_frac_Si_core - 1.) * mSi, wt_frac_Si_core * mO, wt_frac_Si_core * mS \
                      , 0., 0., 0., 0., 0., 0.],
                  [wt_frac_O_core * mFe, (wt_frac_O_core) * mSi, (wt_frac_O_core - 1.) * mO, wt_frac_O_core * mS
                      , 0., 0., 0., 0., 0., 0.],
                  [wt_frac_S_core * mFe, (wt_frac_S_core) * mSi, (wt_frac_S_core) * mO, (wt_frac_S_core - 1.) * mS
                      , 0., 0., 0., 0., 0., 0.],
                  [mFe, mSi, mO, mS, mFe, mMg, mSi, mO, mCa, mAl]])

    # returns number of moles of each element between core and mantle.

    # returns number of moles assuming 100 g planet
    #

    Num_moles = np.linalg.solve(A, b)

    # To do: Adjust so oxygen is last in list
    ## find masses and wt% below for perplex ##

    # Splitting up into lists
    # in order Fe, Si, O, S in core
    Core_moles = Num_moles[:4]

    # in order Fe, Mg, Si, O, Ca, Al in mantle
    Mantle_moles = Num_moles[4:]

    tot_moles_core = sum(Core_moles)

    mass_of_Core = (mFe * (Core_moles[0]) + (mSi * Core_moles[1]) \
                    + ((mO) * Core_moles[2]) + (mS * Core_moles[3]))
    mass_of_Mantle = (mFe * Mantle_moles[0]) + (mMg * Mantle_moles[1]) \
                     + (mSi * Mantle_moles[2]) + (mO * Mantle_moles[3]) \
                     + (mCa * Mantle_moles[4]) + (mAl * Mantle_moles[5])

    Mtot = mass_of_Core + mass_of_Mantle  # in g

    core_mass_frac = mass_of_Core / Mtot
    # Weight percents of mantle oxides
    # Weight percents assuming FeO, MgO, SiO2, CaO, Al2O3

    FeO_mant_wt = Mantle_moles[0] * (mFe + mO) / mass_of_Mantle
    MgO_mant_wt = Mantle_moles[1] * (mMg + mO) / mass_of_Mantle
    SiO2_mant_wt = Mantle_moles[2] * (mSi + (2. * mO)) / mass_of_Mantle
    CaO_mant_wt = Mantle_moles[4] * (mCa + mO) / mass_of_Mantle
    Al2O3_mant_wt = (Mantle_moles[5] / 2.) * (2. * mAl + 3. * mO) / mass_of_Mantle

    # Throw exception, not if statement
    # make inequality not, absolute if. Use machine precision
    total = float(FeO_mant_wt + MgO_mant_wt + SiO2_mant_wt + CaO_mant_wt + Al2O3_mant_wt)
    if total > 1. + (5. * np.finfo(float).eps) or total < 1. - (5. * np.finfo(float).eps):
        print(total)
        print('\n\n Mantle wt% don\'t add to 1')
        sys.exit()

    # same as above edited for printing
    # Fix for PerPlex to input and round to the 8th decimal point
    # converts to percentages

    FeO_mant_wt = abs(round(FeO_mant_wt * 100., 8))
    SiO2_mant_wt = abs(round(SiO2_mant_wt * 100., 8))
    MgO_mant_wt = abs(round(MgO_mant_wt * 100., 8))
    CaO_mant_wt = abs(round(CaO_mant_wt * 100., 8))
    Al2O3_mant_wt = abs(round(Al2O3_mant_wt * 100., 8))

    Mantle_wt_per = {'FeO': FeO_mant_wt, 'SiO2': SiO2_mant_wt, 'MgO': MgO_mant_wt, \
                     'CaO': CaO_mant_wt, 'Al2O3': Al2O3_mant_wt}

    # decimal fraction of materials in CORE by mass, these are perplex inputs (hence the rounding)
    # this is the wt% of Fe and the Fe in FeSi *** so total Fe in core

    Fe_core_wt = (Core_moles[0]) * (mFe) / mass_of_Core
    Si_core_wt = Core_moles[1] * (mSi) / mass_of_Core
    O_core_wt = Core_moles[2] * (mO) / mass_of_Core
    S_core_wt = Core_moles[3] * (mS) / mass_of_Core

    # Throw exception, not if statement
    # make inequality not, absolute if. Use machine precision
    corwt_tot = round(S_core_wt + O_core_wt + Si_core_wt + Fe_core_wt, 8)
    if corwt_tot != 1. and corwt_tot > 1:
        print('\n\n*****Exiting program*****')
        print('Core wt%% don\'t add up')
        print('S_core_wt + O_core_wt + Si_core_wt + Fe_core_wt = %r' % (corwt_tot))
        print('Siwt = %.5f' % (Si_core_wt))
        print('Siwt_input = %.5f' % wt_frac_Si_core)
        print('*************************')
        sys.exit()

    Fe_core_wt = abs(round(Fe_core_wt * 100., 8))
    Si_core_wt = abs(round(Si_core_wt * 100., 8))
    O_core_wt = abs(round(O_core_wt * 100., 8))
    S_core_wt = abs(round(S_core_wt * 100., 8))

    Core_wt_per = {'Fe': Fe_core_wt, 'Si': Si_core_wt, 'O': O_core_wt, 'S': S_core_wt}
    Core_mol_per = {'Fe': Core_moles[0] / tot_moles_core, 'Si': Core_moles[1] / tot_moles_core, \
                    'O': Core_moles[2] / tot_moles_core, 'S': Core_moles[3] / tot_moles_core}

    return (Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac)


# First user must input the ratios
# Radius_planet = 1.



wt_frac_Si_core = 0.
wt_frac_water = 0.
mol_frac_Fe_mantle = 0.0

Pressure_range_mantle_UM = '3000 1400000'
Temperature_range_mantle_UM = '1400 3500'
resolution_UM = '80 140'  # 11200

Pressure_range_mantle_LM = '1250000 6200000'
Temperature_range_mantle_LM = '2200 5000'
resolution_LM = '80 80'  # 6400

wt_frac_O_core = 0.
wt_frac_S_core = 0.

num_mantle_layers = 1500
num_core_layers = 1000
number_h2o_layers = 0

Core_rad_frac_guess = .54

layers = [num_mantle_layers, num_core_layers, number_h2o_layers]

Radius_planet = [0.5, .6, .7, .8, .9, 1., 1.1, 1.2, 1.3, 1.4, 1.5]

Mantle_potential_temp = 1700.

Star = 'Sun'
CaMg = 0.0616595
SiMg = 0.954993
AlMg = 0.0851138
FeMg = 0.812831

compositional_params = [wt_frac_water, FeMg, SiMg, CaMg, AlMg, mol_frac_Fe_mantle, wt_frac_Si_core, \
                        wt_frac_O_core, wt_frac_S_core]

structure_params = [Pressure_range_mantle_UM, Temperature_range_mantle_UM, resolution_UM,
                    Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                    Core_rad_frac_guess, Mantle_potential_temp]

Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = get_percents(*compositional_params)

Mantle_filename = run_perplex.run_perplex(
    *[Mantle_wt_per, compositional_params, [structure_params[0], structure_params[1], structure_params[2]], "test",
      True])
