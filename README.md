# STSQ Dark Energy Model Analysis
This repository contains a class called 'main_class.py' and other files to solve the equation of motion for a Quintessence dark energy model called 'Sharp Transition Scaling Quintessence - STSQ' and to simulate the redshift evolution of these observables according to STSQ: luminosity distance 'd<sub>L</sub>', growth factor combination 'fσ<sub>m</sub>' and the factor 'γ'.
It can also be employed to constrain STSQ using the 'Early Dark Energy - EDE' constraints imposed by the angular position of the first acoustic peak of the spectrum of the 'Cosmic Microwave Background - CMB', as measured by the Planck probe.
By simulating the same observables for the ΛCDM model and only γ for the modified gravity 'DGP' model, the code also evaluates their fractional difference between different models. Some of the plots made with this program are attached in the file 'plots'.

### main_class.py
This file contains the Runge-Kutta algorithm(s) used to solve two differential equations: the single, real scalar field equation (which determines how the STSQ             scalar field evolves in time) for a FLRW universe and the matter density perturbation equation for a FLRW universe (determines how the growth factor evolves in time). The class parameters are 'z<sub>c</sub>', which is the critical redshift, 'η', which determines how steeply the potential decays, 'z<sub>i</sub>', which is the initial redshift chosen for the simulation and 'N', which is the number of points generated.
Solving the equation of motion gives us an observable like the luminosity distance while solving the perturbation equation gives us the growth factor combination fσ<sub>m</sub>; both observables can be computed for the STSQ and ΛCDM case. For the DGP model only γ is computed. This file also contains the functions that evaluate the fractional accuracies for the different variables.
This file does not need modification to run.
All the other files call the class contained in this file and the class parameters can be flexibly changed each time the class is called.

### lumdis(z).py
This file plots the fractional deviation in the luminosity distances of the STSQ and the ΛCDM models for different values of z<sub>c</sub>. One can also vary the initial redshift z<sub>i</sub> or N, the number of reiterations, for greater precision (and a greater number of points).

### lumdis(eta).py
This file plots the fractional deviation in the luminosity distances of the STSQ and the ΛCDM models for different values of η. One can also vary the initial redshift z<sub>i</sub> or N, the number of reiterations, for greater precision (and a greater number of points).

### gdeviation(eta).py
This file plots the fractional difference in the growth factor combination fσ<sub>m</sub> of the STSQ and the ΛCDM models, for different values of η. One can also vary the initial redshift z<sub>i</sub> and N, for greater precision (or a greater number of points).

### gdeviation(z).py
This file plots the fractional difference in the growth factor combination fσ<sub>m</sub> of the STSQ and the ΛCDM models, for different values of z<sub>c</sub>. One can also vary the initial redshift z<sub>i</sub> and N, for greater precision (or a greater number of points).

### gammas(eta).py
This file plots γ for the ΛCDM and the DGP model. γ is simulated also for the STSQ model, for different values of η. One can also vary the initial redshift z<sub>i</sub> or N, the number of reiterations, for greater precision (and a greater number of points).

### gammas(z).py
This file plots γ for the ΛCDM and the DGP model. γ is simulated also for the STSQ model, for different values of z<sub>c</sub>. One can also vary the initial redshift z<sub>i</sub> or N, the number of reiterations, for greater precision (and a greater number of points).

### EDE.py
This file calculates how the luminosity distance evolves from the redshift 'z<sub>dec</sub>' (z at decoupling) to z=0 (z now), to find whether the luminosity distance fractional deviation between ΛCDM and STSQ violates the constraint on the distance imposed by first acoustic peak of the CMB. 
