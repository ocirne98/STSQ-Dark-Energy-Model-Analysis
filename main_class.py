from math import sin, cos, log, exp, sqrt, pi
from numpy import array,arange, linspace, subtract

class model:
    def __init__(self, z_i, z_c, N, eta, decision):

        #the class's variables: initial redshift z_i, reiterations N, eta parameter, critical redshift z_c, True/False for DGP model's boolean
        self.z_i = z_i 
        self.N = N
        self.eta = eta
        self.z_c = z_c
        self.decision = decision


        # SI to natural units
        self.metre = 5.07e15 #1/GeV
        self.second = 1.52e24 #1/GeV
        self.kg = 5.63e26 #GeV

        # Planck mass m, n=log(density)/log(scalefactor), c, current Hubble constant H_0 and hubble, omega density parameters
        self.m = 2.43e18 #GeV
        self.n = 3
        self.c = 3.00e8*self.metre/self.second
        self.H_0 = 1.44e-42 #GeV
        self.omega_de_0 = 0.68
        self.omega_m_0 = 0.32
        self.hubble = 0.674

        # critical and current densities of matter and dark energy
        self.rho_c_SI = (self.hubble**2)*(1.878e-29)*(10**3) #kg/m^3
        self.rho_c = self.rho_c_SI*self.kg/(self.metre**3) #GeV^4
        self.rho_m0 = self.omega_m_0*self.rho_c #GeV^4
        self.rho_de0 = self.omega_de_0*self.rho_c #GeV^4

        # current, initial and critical 'a' values
        self.a_0 = 1
        self.a_i = 1/(self.z_i+1) 
        self.a_c = 1/(self.z_c+1)

        # lambda parameter of STSQ
        self.lam = (self.n*(2*self.omega_de_0+(self.a_0/self.a_c)**self.n*self.omega_m_0)/(2*self.omega_de_0))**(1/2)

        # initial values of phi, x=dphi/dt, V, H, DE density, matter density
        self.phi_i = (self.m/self.lam)*log((2*(self.m**4)*((self.lam**2)/self.n)-1)/(self.rho_m0*((self.a_0/self.a_i)**self.n)))
        self.phi_c = (self.m/self.lam)*log(self.m**4/self.rho_de0)
        self.x_i = sqrt(2*(self.m**4)*(exp(-self.lam*(self.phi_i/self.m))))
        self.V_i = (self.m**4)*exp(-self.lam*(self.phi_i/self.m)) 
        self.rho_i = (self.x_i**2)/2 + self.V_i
        self.rho_mi = self.rho_m0*((self.a_0/self.a_i)**self.n)
        self.H_i = (1/(self.m*sqrt(3)))*sqrt(self.rho_i + self.rho_mi)


        # initial times for LCDM and STSQ models
        self.t_i = 2/(self.n*self.H_i)
        self.tCDM_i = (1/(self.H_0*self.n*(self.omega_de_0**0.5))) * log( ((self.omega_m_0*self.a_i+self.omega_de_0)**0.5 + self.omega_de_0**0.5)/((self.omega_m_0*self.a_i+self.omega_de_0)**0.5 - self.omega_de_0**0.5) )
        self.t_c = self.t_i*((self.a_c/self.a_i)**(self.n/2))
        self.t_f = 100*self.t_c  # choose t_f big enough, arbitrary

        # initial parameters for structure growth
        self.g_i = 1.0
        self.y_i = 0.0

        # time steps for LCDM and STSQ 
        self.h = (self.t_f-self.t_i)/self.N
        self.hCDM = (self.t_f-self.tCDM_i)/self.N

        # time ranges for LCDM and STSQ 
        self.tempo_ax = arange(self.t_i,self.t_f,self.h)
        self.tempoCDM_ax = arange(self.tCDM_i,self.t_f,self.hCDM)

    # STSQ's V for phi<phi_c
    def V(self,phi):
        return (self.m**4)*exp(-self.lam*(phi/self.m))
    
    # STSQ's V for phi>phi_c
    def V2(self,phi):
        return (self.m**4)*exp(-self.lam*(self.phi_c/self.m))*exp(self.eta*self.lam*(phi-self.phi_c)/self.m)

    # STSQ's dV/dphi for phi<phi_c
    def V_var(self,phi):
        return -self.lam*(self.V(phi)/self.m)

    # STSQ's dV/dphi for phi>phi_c
    def V_var2(self,phi):
        return self.eta*self.lam*(self.V2(phi)/self.m)

    # STSQ's DE density
    def rho_de(self,V, x):
        return ((x**2)/2)+V

    # STSQ's DE pressure
    def p_de(self,V, x):
        return ((x**2)/2)-V

    # matter density
    def rho_m(self,a):
        return self.rho_m0*((self.a_0/a)**self.n)

    # STSQ's H variation
    def H(self,V, x, a):
        return (1/(self.m*sqrt(3)))*sqrt(self.rho_m(a) + self.rho_de(V, x))

    # LCDM's H variation
    def H_CDM(self,a):
        return (1/(self.m*sqrt(3)))*sqrt(self.rho_m(a) + self.rho_de0)

    # DGP's H variation
    def H_DGP(self,a):
        H = 0.5*(self.H_0*(1-self.omega_m_0) + sqrt(((self.H_0*(1-self.omega_m_0))**2) + 4*self.rho_m(a)/(3*self.m**2)))
        return H

    # STSQ's barotrpic paramter
    def w(self,V,x):
        return self.p_de(V,x)/self.rho_de(V,x)

    # DE's density parameter for STSQ
    def omegade_STSQ(self,V, x, a):
        return self.rho_de(V, x)/(self.rho_de(V, x)+self.rho_m(a))

    # DE's density parameter for LCDM
    def omegade_CDM(self,a):
        return self.rho_de0/(self.rho_m(a)+self.rho_de0)

    # matter's density parameter for STSQ
    def omegam_STSQ(self,V,x,a):
        return self.rho_m(a)/(self.rho_m(a)+self.rho_de(V,x))

    # matter's density parameter for LCDM
    def omegam_CDM(self,a):
        return self.rho_m(a)/(self.rho_m(a)+self.rho_de0)

    # matter's density parameter for DGP
    def omegam_DGP(self,a):
        r_c = 1/(self.H_0*(1-self.omega_m_0))
        rho_acc = 3*(self.m**2)*self.H_DGP(a)/r_c
        return self.rho_m(a)/(self.rho_m(a)+rho_acc)

    # array updated at each reiteration by R-K algorithm for STSQ; f_i are derivatives in time
    def f(self,V, V_var, r):
        x = r[1]
        a = r[2]

        fphi = x
        fx = -(V_var + 3*self.H(V, x, a)*x)
        fa = a*self.H(V, x, a)

        g = r[3]
        y = r[4]

        fg = y
        fy = -(4*self.H(V,x,a)*y + self.rho_de(V,x)*(1-self.w(V,x))*g/(2*self.m**2))

        return array([fphi,fx,fa,fg,fy],float)
    

    # array updated at each reiteration by R-K algorithm for LCDM; f_i are derivatives in time
    def fCDM(self, a, rCDM):
        g = rCDM[0]
        y = rCDM[1]

        fg = y

        H_prime = -(self.n*self.rho_m0)/(6*(self.m**2)*(a**(self.n+1))*self.H_CDM(a))
        j = ((5/a) + (H_prime/self.H_CDM(a)))*y
        k = self.rho_de0*g/((self.m**2)*(a**2)*(self.H_CDM(a)**2))
        fy = -(j+k)

        return array([fg,fy],float)

    # array updated at each reiteration by R-K algorithm for DGP; f_i are derivatives in time
    def fDGP(self, a, rDGP):
        g = rDGP[0]
        y = rDGP[1]

        fg = y

        r_c = 1/(self.H_0*(1-self.omega_m_0))
        H_primed= -(self.n/(3*self.m**2)) * ( self.rho_m(a) / (a*sqrt((self.H_0/r_c) + (4*self.rho_m(a)/(3*self.m**2)))) )
        if self.decision == True:
            beta = 1 - 2*r_c*(self.H_DGP(a) + a*H_primed/3)
        else:
            beta = float("inf")
        a_dotdot = self.H_DGP(a)**2 + a*self.H_DGP(a)*H_primed

        j = ((5/a) + (H_primed/self.H_DGP(a)))*y
        k = (a_dotdot + 2*self.H_DGP(a)**2 - self.rho_m(a)*(1+1/(3*beta))/(2*self.m**2))*g/((a*self.H_DGP(a))**2)
        fy = -(j+k)

        return array([fg,fy],float) 

    # CDM's R-K algorithm
    def CDM(self, a_ax):
        rCDM = array([self.g_i, self.y_i])
        lista1 = a_ax[1:]
        lista2 = a_ax[0:-1]
        lista3 = subtract(lista1,lista2)
        gCDM_ax = []
        fg_ax = []
        omegam_ax = []
        HCDM_ax = []

        for k in range(len(a_ax[0:-1])): 
            gCDM_ax.append(rCDM[0])
            a = a_ax[k]
            da = lista3[k] 

            fg = self.fCDM(a,rCDM)
            fg_ax.append(fg[0])

            omegam = self.omegam_CDM(a)
            omegam_ax.append(omegam)

            HCDM_ax.append(self.H_CDM(a))

            k_1 = da*self.fCDM(a,rCDM)
            k_2 = da*self.fCDM(a,rCDM + 0.5*k_1)
            k_3 = da*self.fCDM(a,rCDM + 0.5*k_2)
            k_4 = da*self.fCDM(a,rCDM + k_3)
            rCDM += (1/6)*(k_1+2*k_2+2*k_3+k_4)
        #for the last values of g,y
        gCDM_ax.append(rCDM[0]) 
        HCDM_ax.append(self.H_CDM(a_ax[-1]))

        #for the last fg,omegam
        fg = self.fCDM(a,rCDM)
        fg_ax.append(fg[0])
        omegam = self.omegam_CDM(a)
        omegam_ax.append(omegam)
        return array([gCDM_ax,fg_ax,omegam_ax,HCDM_ax])


    # STSQ's R-K algorithm
    def stsq(self):
        r = array([self.phi_i, self.x_i, self.a_i, self.g_i, self.y_i])
        phi_ax =[]
        x_ax = [] 
        a_ax = []
        z_ax = []
        g_ax = []
        fg_ax = []
        omegam_ax = []
        H_ax = []
        omegade_ax = []
        V_ax = []
        w_ax = []
        dens_ax = []

        for k in range(self.N):

            if r[2]<=self.a_0:
                phi_ax.append(r[0])
                x_ax.append(r[1])
                a_ax.append(r[2])
                z_ax.append((1/r[2])-1)
                g_ax.append(r[3])



                if r[0] <= self.phi_c:
                    V_ax.append(self.V(r[0]))
                    w_ax.append(self.w(self.V(r[0]), r[1]))

                    omega = self.omegam_STSQ(self.V(r[0]), r[1], r[2]) 
                    omegam_ax.append(omega)

                    fg = self.f(self.V(r[0]), self.V_var(r[0]), r) 
                    fg_ax.append(fg[3])

                    H = self.H(self.V(r[0]), r[1], r[2])
                    H_ax.append(H)

                    omegade = self.omegade_STSQ(self.V(r[0]), r[1], r[2])
                    omegade_ax.append(omegade)



                    dens = self.rho_de(self.V(r[0]),r[1])
                    dens_ax.append(dens)
                    #hubble = self.H(pot, r[1], r[2])
                    #barpar = self.w(pot,r[1])
                    #dens_par = self.omega_de(pot, r[1], r[2])

                    k_1 = self.h*self.f(self.V(r[0]), self.V_var(r[0]), r)
                    r1 = r + 0.5*k_1
                    k_2 = self.h*self.f(self.V(r1[0]), self.V_var(r1[0]), r1)
                    r2 = r + 0.5*k_2
                    k_3 = self.h*self.f(self.V(r2[0]), self.V_var(r2[0]), r2)
                    r3 = r + k_3
                    k_4 = self.h*self.f(self.V(r3[0]), self.V_var(r3[0]), r3)
                    r += (1/6)*(k_1+2*k_2+2*k_3+k_4)

                else:
                    V_ax.append(self.V2(r[0]))
                    w_ax.append(self.w(self.V2(r[0]),r[1]))

                    omega = self.omegam_STSQ(self.V2(r[0]), r[1], r[2])
                    omegam_ax.append(omega)

                    fg = self.f(self.V2(r[0]), self.V_var2(r[0]), r)
                    fg_ax.append(fg[3])

                    H = self.H(self.V2(r[0]), r[1], r[2])
                    H_ax.append(H)

                    omegade = self.omegade_STSQ(self.V(r[0]), r[1], r[2])
                    omegade_ax.append(omegade)

                    dens = self.rho_de(self.V2(r[0]),r[1])   
                    dens_ax.append(dens)                 
                    #hubble = self.H(pot, r[1], r[2])
                    #barpar = self.w(pot,r[1])
                    #dens_par = self.omega_de(pot, r[1], r[2])

                    k_1 = self.h*self.f(self.V2(r[0]), self.V_var2(r[0]), r)
                    r1 = r + 0.5*k_1
                    k_2 = self.h*self.f(self.V2(r1[0]), self.V_var2(r1[0]), r1)
                    r2 = r + 0.5*k_2
                    k_3 = self.h*self.f(self.V2(r2[0]), self.V_var2(r2[0]), r2)
                    r3 = r + k_3
                    k_4 = self.h*self.f(self.V2(r3[0]), self.V_var2(r3[0]), r3)
                    r += (1/6)*(k_1+2*k_2+2*k_3+k_4)
            else:
                break

        return array([phi_ax, x_ax, a_ax, z_ax, g_ax, fg_ax, omegam_ax, H_ax,omegade_ax,V_ax,w_ax,dens_ax])

    # DGP's R-K algorithm
    def DGP(self,a_ax):
        rDGP = array([self.g_i, self.y_i])
        lista1 = a_ax[1:]
        lista2 = a_ax[0:-1]
        lista3 = subtract(lista1,lista2)
        gDGP_ax = []
        fgDGP_ax = [] 
        HDGP_ax = []
        omegamDGP_ax = []

        for k in range(len(a_ax[0:-1])):
            a = a_ax[k]
            da = lista3[k]

            gDGP_ax.append(rDGP[0])
            fgDGP_ax.append(rDGP[1])
            HDGP_ax.append(self.H_DGP(a))
            omegamDGP_ax.append(self.omegam_DGP(a))


            k_1 = da*self.fDGP(a,rDGP)
            k_2 = da*self.fDGP(a,rDGP + 0.5*k_1)
            k_3 = da*self.fDGP(a,rDGP + 0.5*k_2)
            k_4 = da*self.fDGP(a,rDGP + k_3)
            rDGP += (1/6)*(k_1+2*k_2+2*k_3+k_4)
        a = a_ax[-1]    
        gDGP_ax.append(rDGP[0])
        fgDGP_ax.append(rDGP[1])
        HDGP_ax.append(self.H_DGP(a))
        omegamDGP_ax.append(self.omegam_DGP(a))

        return array([gDGP_ax,fgDGP_ax,HDGP_ax,omegamDGP_ax])

    # luminosity distance function
    def distance(self, lista, lista_da, H_ax):
        s=0
        for j in range(len(lista_da)):
            aa = lista[j]
            da = lista_da[j]
            H = H_ax[j]
            s += da/((aa**2)*H)
        dist = self.c*self.a_0*s
        dista = dist/(self.metre*3.086e22)
        return dista
    
    # gamma factor for STSQ
    def gammaSTSQ(self,a_ax,g_ax,fg_ax,H_ax,omegam_ax):
        gamma_ax = []
        fga_ax = []
        for k in range(len(a_ax)):
            a = a_ax[k]
            g = g_ax[k]
            fg = fg_ax[k]
            H = H_ax[k]

            omegam = omegam_ax[k]
            gamma = log(1+(fg/(g*H)))/log(omegam)
            gamma_ax.append(gamma)

            fga = (1+(fg/(g*H)))*g*a
            fga_ax.append(fga)

        return array([gamma_ax,fga_ax])

    # gamma factor for CDM
    def gammaCDM(self,a_ax,g_ax,fg_ax,omegam_ax):
        gamma_ax = []
        fga_ax = []
        for k in range(len(a_ax)):
            a = a_ax[k]
            g = g_ax[k]
            fg = fg_ax[k]

            omegam = omegam_ax[k]
            gamma = log(1+(a*fg/g))/log(omegam)
            gamma_ax.append(gamma)

            fga = (1+(a*fg/g))*g*a
            fga_ax.append(fga)

        return array([gamma_ax,fga_ax])

    # gamma factor for DGP
    def gammaDGP(self,a_ax,g_ax,fg_ax,omegam_ax):
        gamma_ax = []
        fga_ax = []
        for k in range(len(a_ax)):
            a = a_ax[k]
            g = g_ax[k]
            fg = fg_ax[k]

            omegam = omegam_ax[k]
            gamma = log(1+(a*fg/g))/log(omegam)
            gamma_ax.append(gamma)

            fga = (1+(a*fg/g))*g*a
            fga_ax.append(fga)

        return array([gamma_ax,fga_ax])
