from main_class import model
from numpy import array, sqrt, log, exp, subtract
from pylab import plot,xlabel,ylabel,show,legend,hlines,xlim,savefig,tick_params,figure

#choose your fighter!
N = int(1e5)
z_i = 20
eta = 0.9

figure(figsize= ((15.0, 12.8)))
for i in range(3):
    z_c = [12,15,17]
    colors = ['g','b','c']
    stsq_model = model(z_i, z_c[i], N, eta,True)

    modello = stsq_model.stsq()
    a_ax = modello[2]
    z_ax = modello[3]

    #Hubbles
    HSTSQ_ax = modello[7]
    #growth factor STSQ
    gSTSQ_ax = modello[4]
    #growth factor CDM
    gCDM_ax = stsq_model.CDM(a_ax)[0]
    #omegas
    omegamSTSQ_ax = stsq_model.stsq()[6]
    omegamCDM_ax = stsq_model.CDM(a_ax)[2]
    #gammas
    fgSTSQ_ax = stsq_model.stsq()[5]
    gammaSTSQ_ax = stsq_model.gammaSTSQ(a_ax, gSTSQ_ax, fgSTSQ_ax, HSTSQ_ax, omegamSTSQ_ax)[0]
    fgaSTSQ_ax = stsq_model.gammaSTSQ(a_ax, gSTSQ_ax, fgSTSQ_ax, HSTSQ_ax, omegamSTSQ_ax)[1]

    fgCDM_ax = stsq_model.CDM(a_ax)[1]
    gammaCDM_ax = stsq_model.gammaCDM(a_ax, gCDM_ax, fgCDM_ax, omegamCDM_ax)[0]
    fgaCDM_ax = stsq_model.gammaCDM(a_ax, gCDM_ax, fgCDM_ax, omegamCDM_ax)[1]

    #plots
    plot(z_ax,subtract(fgaCDM_ax,fgaSTSQ_ax)/fgaCDM_ax,colors[i],linewidth=1)

    # this loop finds the z at which the aggregate precision is overtaken
    # if it never does, no value is given
    listone = list(subtract(fgaCDM_ax,fgaSTSQ_ax)/fgaCDM_ax)
    for j in range(len(listone)):
        z = z_ax[j]
        if listone[j]>=0.012:
            print(z)
            break


hlines(0.012, 1, 2,'r', 'dashed',label='', data=None, linewidth=1)
tick_params(labelsize='20', direction='inout',length=5)
xlim([1, 2])
ylabel(r'$\Delta(f(z)\sigma_m(z))$',fontsize=20)
xlabel('$z$',fontsize=20)
legend([r'$z_c = 12$',r'$z_c = 15$',r'$z_c = 17$','WFIRST aggregate error'],fontsize=20,loc='lower right')
savefig('fsigma(z).pdf', format='pdf', dpi=1200)
#show()