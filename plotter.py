from main_class import model
from numpy import array, sqrt, log, exp, subtract
from pylab import plot,xlabel,ylabel,show,legend,hlines,xlim,savefig,axhspan

#choose your fighter!
eta = 0.5
N = int(1e4)
z_i = 20
omega_de_0 = 0.68
omega_m_0 = 0.32
omega_EDE = 0.0036

#sets z_c to lower EDE constraint
f1 = 2*omega_de_0/(omega_m_0*omega_EDE)
f2 = 2*omega_de_0/omega_m_0
z_c = ((f1-f2)**(1/3))
#z_c = 12

#applies chosen parameters to the STSQ model
stsq_model = model(z_i, z_c, N, eta,True)

#cosmo parameters calculator (w,H,a,...) for STSQ
a_ax = stsq_model.stsq()[2]
z_ax = stsq_model.stsq()[3]

#Hubbles
HSTSQ_ax = stsq_model.stsq()[7]
HCDM_ax = stsq_model.CDM(a_ax)[3]

#HDGP_ax = stsq_model.DGP(a_ax)[2]



#distance STSQ
#distSTSQ_ax = stsq_model.distance(a_ax,HSTSQ_ax)
lista1 = a_ax[1:]
lista2 = a_ax[0:-1]
lista3 = subtract(lista1,lista2)
distSTSQ_ax = []
for k in range((len(a_ax))):
    a = a_ax[k]
    lista = a_ax[k:-1]
    lista_da = lista3[k:]
    distSTSQ_ax.append(stsq_model.distance(lista,lista_da,HSTSQ_ax)/a)


#distance CDM 
#distCDM_ax = stsq_model.distance(a_ax,HCDM_ax)
lista1 = a_ax[1:]
lista2 = a_ax[0:-1]
lista3 = subtract(lista1,lista2)
distCDM_ax = []

for k in range((len(a_ax))):
    a = a_ax[k]
    lista = a_ax[k:-1]
    lista_da = lista3[k:]
    distCDM_ax.append(stsq_model.distance(lista,lista_da,HCDM_ax)/a)



#growth factor STSQ
#gSTSQ_ax = stsq_model.stsq()[4]


#growth factor CDM
#gCDM_ax = stsq_model.CDM(a_ax)[0]


#growth factor DGP
#gDGP_ax = stsq_model.DGP(a_ax)[0]


#omegas

#omegamSTSQ_ax = stsq_model.stsq()[6]
#omegamCDM_ax = stsq_model.CDM(a_ax)[2]
#omegadeSTSQ_ax = stsq_model.stsq()[8]
#omegamDGP_ax = stsq_model.DGP(a_ax)[3]

#gammas
"""

fgSTSQ_ax = stsq_model.stsq()[5]
gammaSTSQ_ax = stsq_model.gammaSTSQ(a_ax, gSTSQ_ax, fgSTSQ_ax, HSTSQ_ax, omegamSTSQ_ax)[0]
fgaSTSQ_ax = stsq_model.gammaSTSQ(a_ax, gSTSQ_ax, fgSTSQ_ax, HSTSQ_ax, omegamSTSQ_ax)[1]

fgCDM_ax = stsq_model.CDM(a_ax)[1]
gammaCDM_ax = stsq_model.gammaCDM(a_ax, gCDM_ax, fgCDM_ax, omegamCDM_ax)[0]
fgaCDM_ax = stsq_model.gammaCDM(a_ax, gCDM_ax, fgCDM_ax, omegamCDM_ax)[1]



fgDGP_ax = stsq_model.DGP(a_ax)[1]
gammaDGP_ax = stsq_model.gammaDGP(a_ax,gDGP_ax,fgDGP_ax,omegamDGP_ax)[0]




#plots

plot(z_ax, distCDM_ax, 'g--', linewidth=1)
#plot(z_ax, distSTSQ_ax, 'r--', linewidth=1)
xlabel("z")
ylabel("Luminosity Distance (Mpc) - " r'$d_{L}$')
legend([r'$Lambda$''CDM model','STSQ model'])
show()
#savefig('myimage.pdf', format='pdf', dpi=1200)


"""

plot(z_ax, subtract(distCDM_ax,distSTSQ_ax)/distCDM_ax,'g--', linewidth=1)
hlines(0.0018, z_ax[0], z_ax[-1],'r', 'dashed',label='', data=None, linewidth=1)
ylabel('LumDis Fractional Deviation')
xlabel('z')
xlim([1/2.7, 1/1.2])
legend(['STSQ/$\Lambda$CDM fractional deviation','WFIRST 1$\sigma$ confidence limit'])
show()
#savefig('myimage.pdf', format='pdf', dpi=1200)

"""

plot(a_ax, gCDM_ax,'r--' ,linewidth=1)
plot(a_ax, g_ax,'g--' ,linewidth=1)
ylabel('Growth Factor g')
xlabel('z')


legend([r'$\Lambda$''CDM model','STSQ model'])
show()


plot(z_ax, omegadeSTSQ_ax)
plot(z_ax, omegamSTSQ_ax)
show()


#plot(a_ax, gammaCDM_ax, 'r--',linewidth=1)
plot(a_ax, gammaSTSQ_ax,'g',linewidth=1)
#plot(a_ax, gammaDGP_ax,'g--',linewidth=1)
axhspan(0.54, 0.56, alpha=0.2, color='green')
left, right = xlim()
xlim(left, right) 
hlines(0.54, left, right,'g', 'dashed',label='', data=None, linewidth=1)
hlines(0.56, left, right,'g', 'dashed',label='', data=None, linewidth=1)
legend(['a','b'])
show()

for i in gammaDGP_ax:
    if i>=0.56:
        print(i)
print('done')
"""


"""
plot(z_ax, fgaCDM_ax, 'r--',linewidth=1)
plot(z_ax, fgaSTSQ_ax, 'g--',linewidth=1)
show()


plot(z_ax,abs(subtract(fgaCDM_ax,fgaSTSQ_ax)/fgaCDM_ax),'g',linewidth=1)
hlines(0.012, 1, 2,'r', 'dashed',label='', data=None, linewidth=1)
xlim([1, 2])
ylabel('Fractional Deviation of $g(a)$')
xlabel('$z$')
legend(['1-(g($\Lambda$$CDM$)/g($STSQ$))','WFIRST Aggregate Error - 0.012'],loc='upper right')
show()
"""



