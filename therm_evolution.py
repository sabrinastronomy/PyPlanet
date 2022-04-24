from eos import *
import scipy.constants as sp
import astropy.constants as ap
from interp_CMF import planet_interp
import matplotlib.pyplot as plt
import astropy.constants as const


def ttt(u_arr, temps_arr, r_arr, T_star, R_star, a, T_eq):
    t_arr = [0]
    for i in range(len(temps_arr)-1):
        j = i + 1

        U_2 = u_arr[i]
        R_2 = r_arr[i]

        U_1 = u_arr[j]
        R_1 = r_arr[j]

        T_2 = temps_arr[i]
        T_1 = temps_arr[j]

        if T_2 <= T_eq:
            t_1 = float("NaN")

        else:
            delt = T_s_t(R_2, R_1, U_2, U_1, T_2, T_1, T_star, R_star, a)
            t_1 = t_arr[i] + delt

        t_arr.append(t_1)

    return np.array(t_arr)


def T_s_t(R_p_2, R_p_1, U_2, U_1, T_s_2, T_s_1, T_star, R_star, a):
    R_p_avg = np.mean([R_p_2, R_p_1])
    T_s_avg = np.mean([T_s_2, T_s_1])
    delU = U_2-U_1
    beta = (-sp.sigma * 4 * sp.pi * R_p_avg**2)
    alpha = (T_s_avg**4) - ((T_star**4) * (R_star/(2*a))**2)

    delt = delU/(beta * alpha)
    return abs(delt)

def T_equil(T_star, R_star, a):
        return T_star * (R_star/(2*a))**(1/2)

def T_s_t_plotter(star_type, cmfs_of_interest, masses_of_interest, names_of_interest, lss=['dotted', '-'], a_colors=['k', 'm', 'g'], afrac_values=[0.001, 0.006], host_temp=6000, host_rs=ap.R_sun.value):
    # Get u and r from mass_radius_relations interpolated values for one CMF
    # plt.hlines(y=700, ls='dashed', colors='r', xmin=0, xmax=10e14, label="~ $T_{eq, max}$ discovered (K2-266c)")
    for mass, name, ls in zip(masses_of_interest, names_of_interest, lss):
        u_arr, r_arr, p_c_arr = planet_interp(data, temp_range, cmfs_of_interest, masses_of_interest)
        plt.plot(r_arr, u_arr)
        plt.savefig("r_u.pdf")
        plt.close()
        for afrac, a_color, r_test in zip(afrac_values, a_colors, r_arr):
            a = afrac*ap.au.value
            T_eq = T_equil(host_temp, host_rs, a)
            t_arr = ttt(u_arr, temp_range, r_arr, host_temp, host_rs, a, T_eq)
            print(t_arr[-1] - t_arr[0])

            plt.hlines(y=T_eq, ls='dashed', colors=a_color, xmin=0, xmax=10e13)
            plt.semilogx(t_arr, temp_range, label="$T_{eq} = $" + str(np.round(T_eq, decimals=-1))[:-2] + "K (for " + name + " at $a = {} AU$".format(afrac) + ")", ls=ls, c=a_color)
            plt.legend(loc=1)
            plt.title("Cooling Rates")
            plt.xlabel("$\Delta t$ from $T_{s, max}$")
            plt.ylabel("$T_s$ [K]")
            plt.savefig("num_TE_{}_{}.pdf".format(star_type, a))
            plt.close()


            plt.loglog(t_arr, r_arr)
            plt.ylabel("r [m]")
            plt.xlabel("t [s]")
            plt.title("Radius and Cooling")
            plt.tight_layout()
            plt.savefig("num_R_{}_{}.pdf".format(star_type, a))
            plt.close()




if __name__ == "__main__":
    location = "/Users/sabrinaberger/RockyPlanets"
    data = location + "/thermalData_upper_mantle/adiabatic/"
    save = location + "/thermalData_upper_mantle/adiabatic/"
    temp_range = np.linspace(300, 3000, 10)
    mars_cmf = 0.26
    earth_cmf = 0.33
    cmfs_of_interest = [earth_cmf]

    earth_mass = ap.M_earth.value
    mars_mass = 0.64171e24
    merc_mass = 3.285e23

    # masses_of_interest = [earth_mass, mars_mass]
    masses_of_interest = [earth_mass, mars_mass, merc_mass]
    labels = [r"$M_{\oplus}$", r"$M_{Mars}$", r"$M_{Mercury}$"]
    colors = ['k', 'b', 'r']
    for mass, label, c in zip(masses_of_interest, labels, colors):
        u_arr, r_arr, p_c_arr = planet_interp(data, temp_range, cmfs_of_interest, mass)
        radius_values = np.array(r_arr)/const.R_earth.value
        plt.loglog(temp_range, radius_values, c=c, label=label)

    plt.xlabel(r"T [K]")
    plt.ylabel(r"$R_{p}$ $[R_{\oplus}]$")
    plt.title("Radius and Surface Temperature")
    plt.legend()
    plt.savefig(save + "rad_T_s_log.png")
    plt.close()

    # animate("therm", temp_range)


    T_s_t_plotter("Sun-like host", cmfs_of_interest, masses_of_interest, ["$M_{Earth}$"], afrac_values = [0.1, 1], lss=["-"])
    # T_s_t_plotter("M dwarf host", cmfs_of_interest, masses_of_interest, ["$M_{Earth}$"], afrac_values = [0.01, 0.1], host_temp=3000, host_rs=0.5*ap.R_sun.value, lss=["-"])



