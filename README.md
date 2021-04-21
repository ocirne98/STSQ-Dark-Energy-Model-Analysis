# STSQ Dark Energy Model Analysis
This repository contains a class called 'main_class.py' and other files to solve the equation of motion for a Quintessence dark energy model called 'Sharp Transition Scaling Quintessence - STSQ' and to simulate the redshift evolution of these observables according to STSQ: luminosity distance 'd<sub>L</sub>', growth factor 'G' and the factor 'γ'. It can also be employed to constrain STSQ using the 'Early Dark Energy - EDE' constraints imposed by the angular position of the first acoustic peak of the spectrum of the 'Cosmic Microwave Background - CMB', as measured by the Planck probe. By simulating the same observables for the ΛCDM model and only γ for the modified gravity 'DGP' model, the code also evaluates their fractional difference between different models.

### main_class.py
This file contains the Runge-Kutta algorithm(s) used to solve two differential equations: the single, real scalar field equation (which determines how the STSQ             scalar field evolves in time) for a FLRW universe and the matter density perturbation equation for a FLRW universe (determines how the growth factor evolves in time). The class parameters are 'z<sub>c</sub>', which is the critical redshift, 'η', which determines how steeply the potential decays, 'z<sub>i</sub>', which is the initial redshift chosen for the simulation, 'N', which is the number of points generated and a boolean (True or False). Solving the equation of motion gives us an observable like the luminosity distance while solving the perturbation equation gives us the growth factor; both observables can be computed for the STSQ and ΛCDM case. For the DGP model only γ is computed. This file also contains the functions that evaluate the fractional accuracies for the different variables.
Functions of the file:
- 
This file does not need modification to run.

### gdeviation(z)
This file plots the fractional difference between the growth factor g in the STSQ model and in the Lambda-CDM model, for different values of z_c. The fractional deviation is compared to the WFIRST predicted aggregate error, to find the redshift ranges where the two models are distinguishable. This file can be just run to plot the evolution of the fractional deviation in time: one can vary the initial redshift 'z_i' or 'N', the number of reiterations, for greater precision.

### gdeviation(eta)
This file plots the fractional difference between the growth factor g in the STSQ model and in the Lambda-CDM model, for different values of eta. The fractional deviation is again compared to the WFIRST predicted aggregate error, to find the redshift ranges where the two models are distinguishable. This file can be just run to plot the evolution of the fractional deviation in time: one can vary the initial redshift 'z_i' or 'N', the number of reiterations, for greater precision.

### plotter.py
This file plots the fractional difference between the luminosity distance in the STSQ model and in the Lambda-CDM model. The fractional deviation is again compared to the WFIRST predicted aggregate error, to find the redshift ranges where the two models are distinguishable. This file can be just run to plot the evolution of the fractional deviation in time: one can vary the initial redshift 'z_i' or 'N', the number of reiterations, for greater precision.

### EDE.py
This file calculates how the luminosity distance evolves from the redshift 'z_dec' (z at decoupling) to z=0 (z now), to find whether the luminosity fractional deviation between Lambda-CDM and STSQ violates the constraint on the distance imposed by first acoustic peak of the CMB. 
