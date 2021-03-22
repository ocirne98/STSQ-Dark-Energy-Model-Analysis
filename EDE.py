from math import sin, cos, log, exp, sqrt, pi
from numpy import array,arange, linspace, subtract
from pylab import plot,xlabel,ylabel,show,legend,hlines,xlim,savefig
from main_class import model


#choose your fighter!
eta = 0.5
N = int(1e5)
z_i = 500
a_i = 1/(z_i+1)
z_e = z_i
a_e = a_i
omega_de_0 = 0.68
omega_m_0 = 0.32
omega_EDE = 0.0036 
n = 3

#sets z_c to lower EDE constraint
f1 = 2*omega_de_0/(omega_m_0*omega_EDE)
f2 = 2*omega_de_0/omega_m_0
#z_c = ((f1-f2)**(1/3))
z_c = 4
a_c = 1/(z_c+1)


#applies chosen parameters to the STSQ model
stsq_model = model(z_i, z_c, N, eta,True)

#constants
z_dec = 1089
a_dec = 1/(z_dec+1)
a_0=1

#distance STSQ
#distance 2 (a>a_e)
stsq = stsq_model.stsq()
a2_ax = stsq[2]
da2_ax = subtract(a2_ax[1:],a2_ax[0:-1])
HSTSQ2_ax = stsq[7][0:-1]

dSTSQ2 = stsq_model.distance(a2_ax[0:-1],da2_ax,HSTSQ2_ax)

#distance 1 (a<a_e)
#finds rho_de at a_e
V_i = stsq_model.V_i
x_i = stsq_model.x_i
rho_de_ae = stsq_model.rho_de(V_i,x_i) 

def HSTSQ1(a):
    h = 0.67 #0.674 actually
    omega_r_0 = (2.47e-5)/(h**2)
    rho_r0 = omega_r_0*stsq_model.rho_c
    rho_tot = rho_de_ae + stsq_model.rho_m0*(1/a)**3 + rho_r0*(1/a)**4
    H = sqrt(rho_tot/3)/stsq_model.m
    return(H)


# diagnostic function for 
"""
modello = stsq_model.stsq()[7]
def HSTSQ1(a):
    omegade_e = stsq_model.n/(stsq_model.lam**2)
    omegam_e = 1-omegade_e
    H_e = modello[0]
    omegade = omegade_e*(a/a_e)**n
    omegam = omegam_e*(a_e/a)**n
    H = H_e*sqrt(omegade+omegam)
    return(H)
"""

a1_ax = list(linspace(a_dec,a_e,int(1e5)))
da1_ax = subtract(a1_ax[1:],a1_ax[0:-1])
HSTSQ1_ax = [HSTSQ1(a) for a in a1_ax[0:-1]]

dSTSQ1 = stsq_model.distance(a1_ax[0:-1],da1_ax,HSTSQ1_ax)

dSTSQ = dSTSQ1+dSTSQ2
print('distanceSTSQ =',dSTSQ,'Mpc')



#distance CDM
#distance 1 (a<a_e)
def HCDM1(a):
    h = 0.67 #0.674 actually
    omega_r_0 = (2.47e-5)/(h**2)
    rho_r0 = omega_r_0*stsq_model.rho_c
    rho_tot = stsq_model.rho_de0 + stsq_model.rho_m0*(1/a)**3 + rho_r0*(1/a)**4
    H = sqrt(rho_tot/3)/stsq_model.m
    return(H)
HCDM1_ax = [HCDM1(a) for a in a1_ax[0:-1]]
dCDM1 = stsq_model.distance(a1_ax[0:-1],da1_ax,HCDM1_ax)

#distance 2 (a>a_e)
HCDM2_ax = stsq_model.CDM(a2_ax)[3][0:-1]
dCDM2 = stsq_model.distance(a2_ax[0:-1],da2_ax,HCDM2_ax)
dCDM = dCDM1+dCDM2

print('distanceCDM =',dCDM,'Mpc')


#fractional deviation
print((dCDM-dSTSQ)/dCDM)


# just some diagnostic functions
"""
aCDM1_ax = a1_ax
HCDM1_ax = stsq_model.CDM(aCDM1_ax)[3]
daCDM1_ax = da1_ax
dCDM1 = stsq_model.distance(aCDM1_ax,daCDM1_ax,HCDM1_ax)

aCDM2_ax = a2_ax
HCDM2_ax = stsq_model.CDM(aCDM2_ax)[3]
daCDM2_ax = da2_ax
dCDM2 = stsq_model.distance(aCDM2_ax,daCDM2_ax,HCDM2_ax)

H = stsq_model.stsq()[7]
lista1 = aCDM2_ax[1:]
lista2 = aCDM2_ax[0:-1]
da = subtract(lista1,lista2)
d = stsq_model.distance(aCDM2_ax,da,H)

print(dCDM2,d)

"""


