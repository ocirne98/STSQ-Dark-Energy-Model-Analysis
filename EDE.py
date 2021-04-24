from math import sin, cos, log, exp, sqrt, pi
from numpy import array,arange, linspace, subtract
from pylab import plot,xlabel,ylabel,show,legend,hlines,xlim,savefig
from main_class import model

#choose your fighter!
eta = 0.1
N = int(1e5)
z_i = 800

#constants
omega_de_0 = 0.68
omega_m_0 = 0.32
#omega_EDE = 0.0036
omega_EDE = 0.02

#sets z_c according to the EDE constraint
f1 = 2*omega_de_0/(omega_m_0*omega_EDE)
f2 = 2*omega_de_0/omega_m_0
#z_c = ((f1-f2)**(1/3))
z_c = 750

#conversions
z_e = z_i
a_i = 1/(z_i+1)
a_e = a_i
a_c = 1/(z_c+1)
n = 3

#other constants
z_dec = 1089
a_dec = 1/(z_dec+1)
a_0=1

#applies chosen parameters to the STSQ model
stsq_model = model(z_i, z_c, N, eta,True)

#distances STSQ
#distance STSQ (for a>a_e)
stsq = stsq_model.stsq()
a2_ax = stsq[2]
da2_ax = subtract(a2_ax[1:],a2_ax[0:-1])
HSTSQ2_ax = stsq[7][0:-1]

dSTSQ2 = stsq_model.distance(a2_ax[0:-1],da2_ax,HSTSQ2_ax)


#distance STSQ (a<a_e)
#finds rho_de at a_e
V_i = stsq_model.V_i
x_i = stsq_model.x_i
rho_de_ae = stsq_model.rho_de(V_i,x_i) 

def HSTSQ1(a):
    h = 0.674
    omega_r_0 = (2.47e-5)/(h**2)
    rho_r0 = omega_r_0*stsq_model.rho_c
    rho_tot = rho_de_ae + stsq_model.rho_m0*(1/a)**3 + rho_r0*(1/a)**4
    H = sqrt(rho_tot/3)/stsq_model.m
    return(H)

a1_ax = list(linspace(a_dec,a_e,N))
da1_ax = subtract(a1_ax[1:],a1_ax[0:-1])
HSTSQ1_ax = [HSTSQ1(a) for a in a1_ax[0:-1]]

dSTSQ1 = stsq_model.distance(a1_ax[0:-1],da1_ax,HSTSQ1_ax)

dSTSQ = dSTSQ1+dSTSQ2




#distances LCDM
#distance LCDM (a<a_e)
def HCDM1(a):
    h = 0.674 
    omega_r_0 = (2.47e-5)/(h**2)
    rho_r0 = omega_r_0*stsq_model.rho_c
    rho_tot = stsq_model.rho_de0 + stsq_model.rho_m0*(1/a)**3 + rho_r0*(1/a)**4
    H = sqrt(rho_tot/3)/stsq_model.m
    return(H)
HCDM1_ax = [HCDM1(a) for a in a1_ax[0:-1]]
dCDM1 = stsq_model.distance(a1_ax[0:-1],da1_ax,HCDM1_ax)

#distance LCDM (a>a_e)
HCDM2_ax = stsq_model.CDM(a2_ax)[3][0:-1]
dCDM2 = stsq_model.distance(a2_ax[0:-1],da2_ax,HCDM2_ax)
dCDM = dCDM1+dCDM2




#fractional deviation between LCDM and STSQ
print((dCDM-dSTSQ)/dCDM) #has to be <0.0006
print(((dCDM-dSTSQ)/dCDM)<0.0006)