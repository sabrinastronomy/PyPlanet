from eos import *
import scipy.constants as sp
import astropy.constants as ap
from interpCMF import planet_interp
import matplotlib.pyplot as plt

def ttt(u_arr, temps_arr, r_arr, T_star, R_star, a):
    t_arr = [0]
    for i in range(len(temps_arr)-1):
        j = i + 1

        U_2 = u_arr[i]
        R_2 = r_arr[i]

        U_1 = u_arr[j]
        R_1 = r_arr[j]

        T_2 = temps_arr[i]
        T_1 = temps_arr[j]

        delt = T_s_t(R_2, R_1, U_2, U_1, T_2, T_1, T_star, R_star, a)/22896000
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

def T_s_t_plotter(masses_of_interest, names_of_interest, lss=['dotted', '-'], a_colors=['k', 'm', 'g'], afrac_values=[1, 5, 0.1], host_temp=6000, host_rs=ap.R_sun.value, num_planets=100):
    # Get u and r from mass_radius_relations interpolated values for one CMF
    ts = []
    for mass, name, ls in zip(masses_of_interest, names_of_interest, lss):
        u_arr, r_arr, p_c_arr = planet_interp(data, temp_range, mass)
        for afrac, a_color in zip(afrac_values, a_colors):
            t_arr = ttt(u_arr, temp_range, r_arr, host_temp, host_rs, afrac*ap.au.value)
            print(t_arr[-1]-t_arr[0])
            t_arr = t_arr - t_arr[0]
            plt.semilogx(t_arr, temp_range, label = name + " $a = {} AU$".format(afrac), ls=ls, c=a_color)


    # plt.semilogx(ts[1], temp_range, label="$M_{Mars}$", ls='dotted', c ='k')
    # plt.plot(0, 0, label="$a = 1 AU$",  c ='k')

    plt.legend(loc=3)
    plt.title("Interpolation of {} planets to determine cooling rate".format(num_planets))
    plt.xlabel("time []")
    plt.ylabel("$T_s$ [K]")
    plt.savefig("num_TE.pdf")


if __name__ == "__main__":
    location = "/Users/sabrinaberger/RockyPlanets"
    data = location + "/thermalData/"
    temp_range = np.linspace(300, 3000, 100)[::-1]
    mars_cmf = 0.26
    earth_cmf = 0.33
    cmfs_of_interest = [earth_cmf, mars_cmf]

    earth_mass = ap.M_earth.value
    mars_mass = 0.64171e24

    masses_of_interest = [earth_mass, mars_mass]

    T_s_t_plotter(masses_of_interest, ["$M_{Earth}$", "$M_{Mars}$"] )


