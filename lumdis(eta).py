from main_class import model
from numpy import array, sqrt, log, exp, subtract
from pylab import plot,xlabel,ylabel,show,legend,hlines,xlim,ylim,savefig,axhspan,figure,tick_params

#choose the class parameters!
N = int(1e3)
z_i = 20

#constants
omega_de_0 = 0.68
omega_m_0 = 0.32
omega_EDE = 0.0036

#sets z_c according to the EDE constraint
f1 = 2*omega_de_0/(omega_m_0*omega_EDE)
f2 = 2*omega_de_0/omega_m_0
z_c = ((f1-f2)**(1/3))

#plots of the fractional deviation for different etas
figure(figsize= ((13.0, 10.8)))
for i in range(3):

    etas = [0.1,0.5,0.9] #different values for eta
    colors = ['g','b','c']

    stsq_model = model(z_i, z_c, N, etas[i], True)

    #cosmo parameters calculator (w,H,a,...) for STSQ
    a_ax = stsq_model.stsq()[2]
    z_ax = stsq_model.stsq()[3]

    #Hubble paramater for STSQ and LCDM models
    HSTSQ_ax = stsq_model.stsq()[7]
    HCDM_ax = stsq_model.CDM(a_ax)[3]

    
    #distance STSQ
    lista1 = a_ax[1:]
    lista2 = a_ax[0:-1]
    lista3 = subtract(lista1,lista2)
    distSTSQ_ax = []
    for k in range((len(a_ax))):
        a = a_ax[k]
        lista = a_ax[k:-1]
        lista_da = lista3[k:]
        distSTSQ_ax.append(stsq_model.distance(lista,lista_da,HSTSQ_ax)/a)
    

    
    #distance LCDM 
    lista1 = a_ax[1:]
    lista2 = a_ax[0:-1]
    lista3 = subtract(lista1,lista2)
    distCDM_ax = []
    for k in range((len(a_ax))):
        a = a_ax[k]
        lista = a_ax[k:-1]
        lista_da = lista3[k:]
        distCDM_ax.append(stsq_model.distance(lista,lista_da,HCDM_ax)/a)   

    plot(z_ax, subtract(distCDM_ax,distSTSQ_ax)/distCDM_ax, colors[i], linewidth=1)

    #this loop finds the z where distinguishability starts
    listone = list(subtract(distCDM_ax,distSTSQ_ax)/distCDM_ax)
    for j in range(len(listone)):
        z = z_ax[j]
        if listone[j]<=0.0018:
            print(z)
            break

hlines(0.0018, z_ax[0], z_ax[-1],'r', 'dashed', data=None, linewidth=1)
tick_params(labelsize='18', direction='inout',length=5)
ylabel(r'$\Delta d_{L}$',fontsize=18)
xlabel('z',fontsize=18)
xlim([0.2, 1.7])
legend([r'$\eta = 0.1$',r'$\eta = 0.5$',r'$\eta = 0.9$','WFIRST 1$\sigma$ aggregate error'],fontsize=15)
show()
#savefig('lumdis1.pdf')