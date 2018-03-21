import os
import pexpect as pe
from functions import *

PerPlex_path = "/Users/sabrinaberger/PerPlex"

def run_perplex(*args):

    Mantle_wt_per = args[0]

    FeMg = args[1][1]
    SiMg = args[1][2]
    CaMg = args[1][3]
    AlMg = args[1][4]
    mol_frac_Fe_mantle = args[1][5]
    wt_frac_Si_core = args[1][6]

    Pressure_range_mantle = args[2][0]
    Temperature_range_mantle = args[2][1]
    resolution = args[2][2]

    filename = args[3]

    #RENAME FILENAME
    #TODO: come up with a better way to handle the filenames
    #
    solutionFileName = "test"


    # define perplex inputs in terms of components, this is for legibility

    component1 = 'MGO'
    component2 = 'SIO2'
    component3 = 'FEO'
    component4 = 'CAO'
    component5 = 'AL2O3'
    component6 = 'NA2O'

    plxMan = str(Mantle_wt_per.get('MgO')) + ' ' + str(Mantle_wt_per.get('SiO2')) + ' ' \
             + str(Mantle_wt_per.get('FeO')) + ' ' + str(Mantle_wt_per.get('CaO')) \
             + ' ' + str(Mantle_wt_per.get('Al2O3')) + ' ' + str(0.)  # last value included for Na

    p = pe.spawn(PerPlex_path+"/./build")


    p.sendline(solutionFileName)
    p.sendline(PerPlex_path+'/stx11ver.dat')
    p.sendline(PerPlex_path+'/perplex_option.dat')

    # Transform them (Y/N)?
    p.sendline('N')

    # Calculations with saturated components (Y/N)?
    p.sendline('N')

    # Use chemical potentials, activities or fugacities as independent variables (Y/N)?
    p.sendline('N')

    # Select thermodynamic components from the set:
    p.sendline(component1)  # MGO
    p.sendline(component2)  # SIO2
    p.sendline(component3)  # FEO
    p.sendline(component4)  # CAO
    p.sendline(component5)  # AL2O3
    p.sendline(component6)  # NA2O
    p.sendline('')

    # Specify computational mode:
    p.sendline('2')

    # Make one dependent on the other, e.g., as along a geothermal gradient (y/n)?
    p.sendline('N')

    # Select x-axis variable:
    p.sendline('2')

    # Enter minimum and maximum values, respectively, for: T(K)
    p.sendline(Temperature_range_mantle)

    # Enter minimum and maximum values, respectively, for: P(bar)
    # P(Pa) = P(bar)*100000
    p.sendline(Pressure_range_mantle)

    # Specify component amounts by weight (Y/N)?
    p.sendline('Y')

    # Enter weight amounts of the components:
    # MGO SIO2 FEO CAO AL2O3
    # for the bulk composition of interest:
    # NOTE*** This is a wt%
    p.sendline(plxMan)
    # Output a print file (Y/N)?
    p.sendline('N')
    # Exclude pure and/or endmember phases (Y/N)?
    p.sendline('N')

    # Include solution models (Y/N)?
    p.sendline('Y')
    p.sendline(PerPlex_path+'/stx11_solution_model.dat')
    p.sendline('C2/c') #C2C Phase of clinopyroxene
    p.sendline('Wus')
    p.sendline('Pv')
    p.sendline('Pl')
    p.sendline('Sp')
    p.sendline('O')
    p.sendline('Wad')
    p.sendline('Ring')
    p.sendline('Opx')
    p.sendline('Cpx')
    p.sendline('Aki')
    p.sendline('Gt_maj') #kitchen sink
    p.sendline('Ppv')
    p.sendline('CF')
    p.sendline('')

    # Enter calculation title:
    p.sendline(solutionFileName + 'calc')

    p.logfile = open('build.log', 'wb')
    p.read()
    p.wait()

    print("Done with Build, moving on to Vertex")


    # Spawn Vertex ----------------#
    # Enter the project name (the name assigned in BUILD) [default = my_project]:
    p = pe.spawn(PerPlex_path+"/./vertex",timeout=2400)

    p.sendline(solutionFileName )

    #p.expect('$$$$',timeout=None)
    p.logfile = open('vertex.log', 'wb')
    p.read()
    p.wait()

    print('Finished with Vertex, beginning Werami')



    p = pe.spawn(PerPlex_path+"/./werami",timeout=None)


    p.sendline(solutionFileName )
    # select 2D grid
    p.sendline('2')
    # Below, select parameters density, alpha, cp.
    # Ns for no calculating individual phase properties
    p.sendline('2')
    p.sendline('N')
    p.sendline('12')
    p.sendline('N')
    p.sendline('13')
    p.sendline('N')
    p.sendline('14')
    p.sendline('N')
    p.sendline('4')
    p.sendline('N')
    p.sendline('19')
    p.sendline('N')
    ####### the next lines will pass requests to perplex to print phases and their proportions into the .tab file

    # 21 species, in all for Fe-Si-Mg-O regime
    # p.sendline('7')
    # p.sendline('C2/c')  # 0
    # p.sendline('7')
    # p.sendline('Wus')  # 1
    # p.sendline('7')
    # p.sendline('Pv')  # 2
    # p.sendline('7')
    # p.sendline('an')  # 3
    # p.sendline('7')
    # p.sendline('Sp')  #4
    # p.sendline('7')
    # p.sendline('O')  # 4
    # p.sendline('7')
    # p.sendline('Wad')  # 5
    # p.sendline('7')
    # p.sendline('Ring')  # 6  #if statement about no FeO or some shit
    # p.sendline('7')
    # p.sendline('Opx')  # 7
    # p.sendline('7')
    # p.sendline('Cpx')  # 8
    # p.sendline('7')
    # p.sendline('Aki')  # 9
    # p.sendline('7')
    # p.sendline('Gt_maj')  # 10
    # p.sendline('7')
    # p.sendline('Ppv')  # 11
    # p.sendline('7')
    # p.sendline('CF')   # 12
    # p.sendline('7')
    # p.sendline('st')  # 12
    # p.sendline('7')
    # p.sendline('q')  # 13
    # p.sendline('7')
    # p.sendline('ca-pv')  # 14
    # p.sendline('7')
    # p.sendline('cfs')  # 15
    # p.sendline('7')
    # p.sendline('coe')  # 16
    # p.sendline('7')
    # p.sendline('ky')  # 17
    # p.sendline('7')
    # p.sendline('seif')  # 18

    # exit parameter choosing

    p.sendline('0')

    #  Change default variable range (y/n)?
    p.sendline('N')

    p.sendline('0')
    p.sendline('0')
    p.sendline('0')


    # Enter number of nodes in the T(K)     and P(bar)   directions:

    p.sendline(resolution)
    p.logfile = open('werami.log','wb')

    p.read()
    p.terminate()
    print("Done with PerPlex")

    successful = True

    #
    # os.rename('build.log', 'ERROR_'+solutionFileName +'_build.log')
    # os.rename('vertex.log', 'ERROR_'+solutionFileName +'_vertex.log')
    # os.rename('werami.log', 'ERROR_'+solutionFileName +'_werami.log')



# os.remove(solutionFileName +'.arf')
# os.remove(solutionFileName +'.blk')
# os.remove(solutionFileName +'.dat')
# os.remove(solutionFileName +'.plt')
# os.remove(solutionFileName +'.tof')
# os.remove(solutionFileName +'_VERTEX_options.txt')
# os.remove(solutionFileName +'_WERAMI_options.txt')
# os.remove(solutionFileName +'_auto_refine.txt')

    filename = '../Solutions/'+filename

    return filename


#TEST
Radius_planet = 1.3

Star = 'Sun'
CaMg =0.0616595
SiMg =0.954993
AlMg = 0.0851138
FeMg = 0.812831

wt_frac_Si_core = 0.
wt_frac_water = 0.
mol_frac_Fe_mantle = 0.0

#(Perplex) P&T parameter space definitions for perplex
Pressure_range_mantle_UM = '3000 1400000'
Temperature_range_mantle_UM = '1400 3000'
resolution_UM = '60 60'

Pressure_range_mantle_LM = '1250000 6500000'
Temperature_range_mantle_LM = '2500 5000'
resolution_LM = '50 50'

wt_frac_O_core = 0.
wt_frac_S_core = 0.

#layers, like concentric shells set here in each region: core, mantle, h20 envelope
num_mantle_layers = 2000
num_core_layers = 1000
number_h2o_layers = 0

#temperature at surface if no water layer. Essentially temperature below the crust
Mantle_potential_temp = 1700.

#h2o potential Temp, surface temperature if there exists an h2o layer
T_surface_h2o = 300. # K

#initialize planet with these guesses for radial fraction of core and water layer
Core_rad_frac_guess = .54
h20_radfrac_guess = 0.1


#lists of compositional and structural inputs used to build planet
# lists of compositional and structural inputs used to build planet
compositional_params = [wt_frac_water, FeMg, SiMg, CaMg, AlMg, mol_frac_Fe_mantle, wt_frac_Si_core, wt_frac_O_core, wt_frac_S_core]

structure_params = [Pressure_range_mantle_UM, Temperature_range_mantle_UM, resolution_UM,
                    Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                    Core_rad_frac_guess, Mantle_potential_temp, h20_radfrac_guess, T_surface_h2o]

# Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = functions.get_percents(*compositional_params)
# Mantle_filename = run_perplex(*[Mantle_wt_per,compositional_params,[structure_params[0],structure_params[1],structure_params[2]],"test",True])
#
#
Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = get_percents(*compositional_params)

# Mantle_filename = run_perplex.run_perplex(
#     *[Mantle_wt_per, compositional_params, [structure_params[0], structure_params[1], structure_params[2]], "test",
#       True])
