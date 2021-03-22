from main_class import model
from numpy import array, sqrt, log, exp, subtract
from pylab import plot,xlabel,ylabel,show,legend,hlines,xlim,savefig

#choose your fighter!
#eta = 0.8 
N = int(1e5)
z_i = 20
z_c = 10.55

for i in range(7):
    etas = [0.2,0.3,0.4,0.5,0.6,0.7,0.8]
    colors = ['g','b','c','m','g--','b--','c--']
    stsq_model = model(z_i, z_c, N, etas[i],True)

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
    #plot
    plot(z_ax,subtract(fgaCDM_ax,fgaSTSQ_ax)/fgaSTSQ_ax,colors[i],linewidth=1)
    #plot(z_ax,fgaCDM_ax,'m--',linewidth=0.8)
    #plot(z_ax,fgaSTSQ_ax,colors[i],linewidth=0.8)

xlim([1, 2])
hlines(0.012, 1, 2,'r', 'dashed',label='', data=None, linewidth=1)
ylabel('Fractional Deviation of $f(z)G(z)$')
#ylabel('$f(z)G(z)$')
xlabel('$z$')
legend(['eta = 0.2','eta = 0.3','eta = 0.4','eta = 0.5','eta = 0.6','eta = 0.7','eta = 0.8','WFIRST Agg Error'])
show()



