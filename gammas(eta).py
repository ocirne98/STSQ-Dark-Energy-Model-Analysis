from main_class import model
from numpy import array, sqrt, log, exp, subtract
from pylab import plot,xlabel,ylabel,show,legend,hlines,xlim,savefig,axhspan,figure,tick_params

#choose your fighter!
N = int(1e5)
z_i = 20
omega_de_0 = 0.68
omega_m_0 = 0.32
omega_EDE = 0.0036

#sets z_c to lower EDE constraint
f1 = 2*omega_de_0/(omega_m_0*omega_EDE)
f2 = 2*omega_de_0/omega_m_0
z_c = ((f1-f2)**(1/3))

figure(figsize= ((15.0, 12.8)))
for i in range(3):

    #applies chosen parameters to the STSQ model
    etas = [0.1,0.5,0.9]
    stsq_model = model(z_i, z_c, N, etas[i], True)
    colors = ['g','b','c']

    #cosmo parameters calculator (w,H,a,...) for STSQ
    stsq = stsq_model.stsq()

    a_ax = stsq[2]
    z_ax = stsq[3]

    #Hubbles
    HSTSQ_ax = stsq[7]

    CDM = stsq_model.CDM(a_ax)
    HCDM_ax = CDM[3]

    DGP = stsq_model.DGP(a_ax)
    HDGP_ax = DGP[2]

    #growth factor STSQ
    gSTSQ_ax = stsq[4]


    #growth factor CDM
    gCDM_ax = CDM[0]


    #growth factor DGP
    gDGP_ax = DGP[0]


    #omegas
    omegamSTSQ_ax = stsq[6]
    omegamCDM_ax = CDM[2]
    omegadeSTSQ_ax = stsq[8]
    omegamDGP_ax = DGP[3]

    #gammas
    fgSTSQ_ax = stsq[5]

    gammaSTSQ = stsq_model.gammaSTSQ(a_ax, gSTSQ_ax, fgSTSQ_ax, HSTSQ_ax, omegamSTSQ_ax)
    gammaSTSQ_ax = gammaSTSQ[0]
    fgaSTSQ_ax = gammaSTSQ[1]


    fgCDM_ax = CDM[1]

    gammaCDM = stsq_model.gammaCDM(a_ax, gCDM_ax, fgCDM_ax, omegamCDM_ax)
    gammaCDM_ax = gammaCDM[0]

    fgaCDM_ax = gammaCDM[1]


    fgDGP_ax = DGP[1]

    gammaDGP = stsq_model.gammaDGP(a_ax,gDGP_ax,fgDGP_ax,omegamDGP_ax)
    gammaDGP_ax = gammaDGP[0]

    plot(z_ax, gammaSTSQ_ax,colors[i],linewidth=1)

    # this loop finds the z at which the aggregate precision is overtaken
    # if it never does, no value is given
    listone = list(gammaSTSQ_ax)
    for j in range(len(listone)):
        z = z_ax[j]
        if listone[j]>=0.543:
            print(z)
            break

plot(z_ax, gammaCDM_ax, 'r--',linewidth=1)
#plot(z_ax, gammaDGP_ax,'m--',linewidth=1)
tick_params(labelsize='20', direction='inout',length=5)
axhspan(0.543, 0.557, alpha=0.2, color='green')
left, right = xlim()
xlim(left, right) 
hlines(0.543, left, right,'g', 'dashed',label='', data=None, linewidth=1)
hlines(0.557, left, right,'g', 'dashed',label='', data=None, linewidth=1)
xlabel(r'$z$',fontsize=20)
ylabel(r'$\gamma(z)$',fontsize=20)
legend([r'$\eta(STSQ) = 0.1$',r'$\eta(STSQ) = 0.5$',r'$\eta(STSQ) = 0.9$',r'$\Lambda$''CDM','EUCLID aggregate error\'s span'],fontsize=20)
show()
#savefig('gammas(eta).pdf', format='pdf', dpi=1200)