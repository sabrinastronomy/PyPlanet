#!python
#cython: boundscheck=False

import numpy as np
from scipy.optimize import fsolve
from libc cimport math

cpdef double P_vinet(double V, double V0, double K0, double KP0):
    cdef double x=V/V0
    cdef double pressure=3.0*K0*x**(-2.0/3.0)*(1.0-x**(1.0/3.0))*math.exp(1.5*(KP0-1.0)*(1.0-x**(1.0/3.0)))
    return pressure

cpdef double b_V(double V,double b0,double b1,double b2,double b3,double b4,double V0):
    cdef double x=V/V0-1.0
    cdef double bv0=b0
    cdef double bv1=b1*x
    cdef double bv2=b2*x**2.0
    cdef double bv3=b3*x**3.0
    cdef double bv4=b4*x**4.0
    cdef double total=bv0+bv1+bv2+bv3+bv4
    return total

cpdef double b_V_deriv(double V,double b0,double b1,double b2,double b3,double b4,double V0):
    cdef double x=V/V0-1.0
    cdef double bv0=0.0
    cdef double bv1=b1/V0
    cdef double bv2=b2*2.0/V0*x
    cdef double bv3=b3*3.0/V0*x**2.0
    cdef double bv4=b4*4.0/V0*x**3.0
    cdef double total=bv0+bv1+bv2+bv3+bv4
    return total

cpdef double b_V_deriv2(double V,double b0,double b1,double b2,double b3,double b4,double V0):
    cdef double x=V/V0-1.0
    cdef double bv0=0.0
    cdef double bv1=0.0
    cdef double bv2=b2*2.0/V0**2.0
    cdef double bv3=b3*6.0/V0**2.0*x
    cdef double bv4=b4*12.0/V0**2.0*x**2.0
    cdef double total=bv0+bv1+bv2+bv3+bv4
    return total

cpdef double f_T(double T, double T0, double m):
    cdef double f=(T/T0)**m-1.0
    return f

cpdef double f_T_deriv(double T, double T0, double m):
    cdef double f_deriv=m/T0*(T/T0)**(m-1.0)
    return f_deriv

cpdef double P_E(double T,double V,double b0,double b1,double b2,double b3,double b4,double V0,double T0,double m):
    #convert rho to V first
    cdef double pressure=-b_V_deriv(V,b0,b1,b2,b3,b4,V0)*f_T(T,T0,m)
    return pressure

cpdef double P_S(double T,double V,double gamma0,double gamma0P,double T0,double V0,double b0,double b1,double b2,double b3,double b4,double m):
    cdef double a1=6.0*gamma0
    cdef double a2=-12.0*gamma0+36.0*gamma0**2.0-18.0*gamma0P
    cdef double f=0.5*((V0/V)**(2.0/3.0)-1.0)
    cdef double gamma0S=((2.0*f+1.0)*(a1+a2*f))/(6.0*(1.0+a1*f+0.5*a2*f**2.0))
    cdef double T0S=T0*(1.0+a1*f+0.5*a2*f**2.0)**0.5
    cdef double CV0S=b_V(V,b0,b1,b2,b3,b4,V0)*f_T_deriv(T0S, T0, m)+620.86 #620.86->1.5Nkb
    cdef double first=b_V_deriv(V,b0,b1,b2,b3,b4,V0)/(m-1.0)*(T*(f_T_deriv(T,T0,m)-f_T_deriv(T0S,T0,m))-T0*(f_T_deriv(T0,T0,m)-f_T_deriv(T0S,T0,m)))
    cdef double second=gamma0S*CV0S*(T-T0)/V
    cdef double pressure=first+second
    return pressure
    
cpdef double P_liquidpv(double V,double T,double V0,double K0,double KP0,double b0,double b1,double b2,double b3,double b4,double gamma0,double gamma0P,double T0,double m):
    cdef double P_total=P_vinet(V,V0,K0,KP0)+P_E(T,V,b0,b1,b2,b3,b4,V0,T0,m)+P_S(T,V,gamma0,gamma0P,T0,V0,b0,b1,b2,b3,b4,m)
    return P_total

def rho_liquidpv(guess,pressure,T,V0,K0,KP0,b0,b1,b2,b3,b4,gamma0,gamma0P,T0,m):
    guess_V=1.0/guess
    def func(x):
        func=P_liquidpv(x,T,V0,K0,KP0,b0,b1,b2,b3,b4,gamma0,gamma0P,T0,m)-pressure
        return func
    solution_V=fsolve(func, x0=guess_V)[0]
    solution=1.0/solution_V
    return solution

cpdef double CV(double T,double V,double b0,double b1,double b2,double b3,double b4,double V0,double T0,double m):
    cdef double CV=b_V(V,b0,b1,b2,b3,b4,V0)*f_T_deriv(T,T0,m)+620.86
    return CV

cpdef double gamma0S_T0S(double V,double V0,double gamma0,double gamma0P,double T0):
    cdef double a1=6.0*gamma0
    cdef double a2=-12.0*gamma0+36.0*gamma0**2.0-18.0*gamma0P
    cdef double f=0.5*((V0/V)**(2.0/3.0)-1.0)
    cdef double gamma0S=((2.0*f+1.0)*(a1+a2*f))/(6.0*(1.0+a1*f+0.5*a2*f**2.0))
    return gamma0S

cpdef double T0S_gamma0S(double V,double V0,double gamma0,double gamma0P,double T0):
    cdef double a1=6.0*gamma0
    cdef double a2=-12.0*gamma0+36.0*gamma0**2.0-18.0*gamma0P
    cdef double f=0.5*((V0/V)**(2.0/3.0)-1.0)
    cdef double T0S=T0*(1.0+a1*f+0.5*a2*f**2.0)**0.5
    return T0S

cpdef double S_pot(double T,double V,double T0S,double b0,double b1,double b2,double b3,double b4,double V0,double T0,double m):
    #T0S is calculated using gamma0S_T0S
    cdef double first=b_V(V,b0,b1,b2,b3,b4,V0)/(m-1.0)*(f_T_deriv(T,T0,m)-f_T_deriv(T0S,T0,m))
    cdef double second=620.86*math.log(T/T0S)
    return first+second

cpdef double gamma_melt(double T,double V,double T0S,double gamma0S,double delta_S_pot,double b0,double b1,double b2,double b3,double b4,double V0,double T0,double m):
    # T0S gamma0S->gamma0S_T0S delta_S_pot->S_pot
    cdef double CV0S=CV(T0S,V,b0,b1,b2,b3,b4,V0,T0,m)
    cdef double CV_T=CV(T,V,b0,b1,b2,b3,b4,V0,T0,m)
    cdef double bVP=b_V_deriv(V,b0,b1,b2,b3,b4,V0)
    cdef double bV=b_V(V,b0,b1,b2,b3,b4,V0)
    cdef double gamma=gamma0S*CV0S/CV_T+V*bVP/bV*(delta_S_pot/CV_T)
    return gamma

cpdef double K_T_melt(double V,double T,double T0S,double V0,double K0,double KP0,double b0,double b1,double b2,double b3,double b4,double gamma0,double gamma0P,double T0,double m):
    #T0S->calculated using gamma0S_T0S
    #dPv/dV
    cdef double x=V/V0
    cdef double a=x**(-2.0/3.0)
    cdef double b=1.0-x**(1.0/3.0)
    cdef double c=math.exp(1.5*(KP0-1.0)*(1.0-x**(1.0/3.0)))
    cdef double ap=-2.0/3.0*x**(-5.0/3.0)/V0
    cdef double bp=-1.0/3.0*x**(-2.0/3.0)/V0
    cdef double cp=c*1.5*(KP0-1.0)*bp
    cdef double dPvdV=3.0*K0*(ap*b*c+a*bp*c+a*b*cp)
    #dPE/dV
    cdef double bVP2=b_V_deriv2(V,b0,b1,b2,b3,b4,V0)
    cdef double dPEdV=-bVP2*f_T(T,T0,m)
    #dPS/dV
    ##dPS1/dV
    cdef double dPS1dV=(T*(f_T_deriv(T,T0,m)-f_T_deriv(T0S,T0,m))-T0*(f_T_deriv(T0,T0,m)-f_T_deriv(T0S,T0,m)))/(m-1.0)*bVP2
    ##dPS2/dV
    ###PS2=gamma0S*CV0S*(T-T0)/V;
    cdef double a1=6.0*gamma0
    cdef double a2=-12.0*gamma0+36.0*gamma0**2.0-18.0*gamma0P
    cdef double f=0.5*((V0/V)**(2.0/3.0)-1.0)
    ####dgamma0S/dV; gamma0S=gamma0S1*gamma0S2*gamma0S3
    cdef double gamma0S1=(2.0*f+1.0)/6.0
    cdef double gamma0S2=a1+a2*f
    cdef double gamma0S3=(1.0+a1*f+0.5*a2*f**2.0)**(-1.0)
    cdef double gamma0S1P=1.0/3.0
    cdef double gamma0S2P=a2
    cdef double gamma0S3P=-(1.0+a1*f+0.5*a2*f**2.0)**(-2.0)*(a1+a2*f)
    cdef double PS2a=gamma0S1*gamma0S2*gamma0S3
    cdef double dgamma0Sdf=gamma0S1P*gamma0S2*gamma0S3+gamma0S1*gamma0S2P*gamma0S3+gamma0S1*gamma0S2*gamma0S3P
    cdef double dfdV=0.5*((-2.0/3.0)*x**(-5.0/3.0)/V0)
    cdef double dgamma0SdV=dgamma0Sdf*dfdV
    ####dCV0S/dV
    cdef double PS2b=CV(T0S,V,b0,b1,b2,b3,b4,V0,T0,m)
    cdef double dCV0SdV=b_V_deriv(V,b0,b1,b2,b3,b4,V0)*f_T_deriv(T0S,T0,m)
    ####d((T-T0)/V)/dV
    cdef double PS2c=(T-T0)/V
    cdef double dPS2cdV=(T0-T)/V**2.0

    cdef double dPS2dV=dgamma0SdV*PS2b*PS2c+PS2a*dCV0SdV*PS2c+PS2a*PS2b*dPS2cdV

    cdef double KT=-V*(dPvdV+dPEdV+dPS1dV+dPS2dV)
    return KT

cpdef double alpha_melt(double gamma,double rho,double CV_T,double KT):
    # gamma->gamma_melt; CV_T->CV; KT->K_T_melt; rho->1/V
    cdef double alpha=gamma*rho*CV_T/KT
    return alpha

cpdef double CP(double T,double gamma,double alpha,double CV_T):
    #gamma->gamma_melt;alpha->alpha_melt;CV_T->CV
    cdef double CP=CV_T*(1.0+gamma*alpha*T)
    return CP

cpdef double KS_melt(double T,double alpha,double gamma,double KT):
    #gamma->gamma_melt;alpha->alpha_melt;KT->K_T_melt
    cdef double KS=KT*(1.0+alpha*gamma*T)
    return KS

cpdef double dTdP(double T,double gamma,double KS):
    #gamma->gamma_melt; KS->KS_melt
    cdef double dTdP=gamma*T/KS
    return dTdP


cdef double[:] T_ref=np.linspace(1800.0, 10000.0, 1641)
cdef double[:] P=np.linspace(1.0*10.0**5.0, 1500.0*10.0**9.0, 3001)
cdef double[:] T=np.zeros(3001)

cdef double hP=P[1]-P[0]

cdef double gamma0S_value
cdef double T0S_value
cdef double delta_S_pot_value
cdef double gamma_value
cdef double KT_value
cdef double CV_T_value
cdef double alpha_value
cdef double KS_value

cdef double rho_value
cdef double volume
cdef double k1
cdef double k2
cdef double k3
cdef double k4

cdef Py_ssize_t i
cdef Py_ssize_t j

cdef double G=6.674e-11 #SI
cdef double k_b=1.38e-23 #J/K
cdef double N_A=6.022e+23 # mol^-1

#EoS parameter
cdef double rho0=2574.80
cdef double K0=13.2*10.0**9.0
cdef double KP0=8.238
cdef double gamma0=0.1899
cdef double gamma0P=-1.94
cdef double V0=1.0/rho0
cdef double m=0.6
cdef double mass=1.0
cdef double T0=3000.0
cdef double b0=4804882.956*(0.9821)
cdef double b1=4804882.956*(0.615)
cdef double b2=4804882.956*(1.31)
cdef double b3=4804882.956*(-3.0)
cdef double b4=4804882.956*(-4.1)

# Define the pressure and temperature range. 
cdef double[:] pressure=np.linspace(1.0*10.0**5.0, 1500.0*10.0**9.0, 1501)
cdef double[:] temperature=np.linspace(10, 4000, 627)

cdef double guess=5000.0
cdef double density_value

# Calculate and save to file. 
# f=open("density_liquidpv.txt", "w")
# for i in range(len(temperature)):
#     for j in range(len(pressure)):
#         density_value=rho_liquidpv(guess,pressure[j],temperature[i],V0,K0,KP0,b0,b1,b2,b3,b4,gamma0,gamma0P,T0,m)
#         guess=density_value
#         f.write("%s %s %s\n"%(temperature[i], pressure[j], density_value))
#     f.write("\n")
# f.close()

#### heat capacity stuff below:
pre_adiabat=open('density_Cp_melt.txt', 'w')
for j in range(len(T_ref)):
    T[0]=T_ref[j]
    for i in range(len(P)):
        rho_value=rho_liquidpv(guess,P[i],T[i],V0,K0,KP0,b0,b1,b2,b3,b4,gamma0,gamma0P,T0,m)
        volume=1.0/rho_value
        gamma0S_value=gamma0S_T0S(volume,V0,gamma0,gamma0P,T0)
        T0S_value=T0S_gamma0S(volume,V0,gamma0,gamma0P,T0)
        delta_S_pot_value=S_pot(T[i],volume,T0S_value,b0,b1,b2,b3,b4,V0,T0,m)
        gamma_value=gamma_melt(T[i],volume,T0S_value,gamma0S_value,delta_S_pot_value,b0,b1,b2,b3,b4,V0,T0,m)
        KT_value=K_T_melt(volume,T[i],T0S_value,V0,K0,KP0,b0,b1,b2,b3,b4,gamma0,gamma0P,T0,m)
        CV_T_value=CV(T[i],volume,b0,b1,b2,b3,b4,V0,T0,m)
        alpha_value=alpha_melt(gamma_value,rho_value,CV_T_value,KT_value)
        CP_T_value=CP(T[i],gamma_value, alpha_value, CV_T_value)
        guess=rho_value
        pre_adiabat.write("%s %s %s\n"%(temperature[j], pressure[i], CP_T_value))

    pre_adiabat.write("\n")
print("heat capacity written")
pre_adiabat.close()