# STSQ Dark Energy Model Analysis
This repository contains a class called 'main_class.py' and other files to solve the equation of motion for a Quintessence dark energy model called 'Sharp Transition Scaling Quintessence - STSQ' and to simulate the redshift evolution of luminosity distances and growth factors for this type of dark energy. It can also be employed to constrain STSQ using the 'Early Dark Energy - EDE' constraints imposed by the angular position of the first acoustic peak of the spectrum of the 'Cosmic Microwave Background - CMB', as measured by the Planck probe.

### main_class.py
This file contains the Runge-Kutta algorithm(s) used to solve two differential equations: the single, real scalar field equation (which determines how the STSQ             scalar field evolves in time) for a FLRW universe and the matter density perturbation equation for a FLRW universe (determines how the growth factor evolves in time). What distinguishes the STSQ field from other scalar fields is its potential: it has an exponential functional form that changes at a critical redshift z<sub>c</sub> text. 'z_c'; it also has a paramter called 'eta', determining how steeply the potential decays. Solving the first equaiton gives us an observable like the luminosity distance, solving the second one gives us the growth factor; both observables can be computed for the STSQ and Lambda-CDM case. This file does not need modification to run.

### gdeviation(z)
This file plots the fractional difference between the growth factor g in the STSQ model and in the Lambda-CDM model, for different values of z_c. The fractional deviation is compared to the WFIRST predicted aggregate error, to find the redshift ranges where the two models are distinguishable. This file can be just run to plot the evolution of the fractional deviation in time: one can vary the initial redshift 'z_i' or 'N', the number of reiterations, for greater precision.

### gdeviation(eta)
This file plots the fractional difference between the growth factor g in the STSQ model and in the Lambda-CDM model, for different values of eta. The fractional deviation is again compared to the WFIRST predicted aggregate error, to find the redshift ranges where the two models are distinguishable. This file can be just run to plot the evolution of the fractional deviation in time: one can vary the initial redshift 'z_i' or 'N', the number of reiterations, for greater precision.

### plotter.py
This file plots the fractional difference between the luminosity distance in the STSQ model and in the Lambda-CDM model. The fractional deviation is again compared to the WFIRST predicted aggregate error, to find the redshift ranges where the two models are distinguishable. This file can be just run to plot the evolution of the fractional deviation in time: one can vary the initial redshift 'z_i' or 'N', the number of reiterations, for greater precision.

### EDE.py
This file calculates how the luminosity distance evolves from the redshift 'z_dec' (z at decoupling) to z=0 (z now), to find whether the luminosity fractional deviation between Lambda-CDM and STSQ violates the constraint on the distance imposed by first acoustic peak of the CMB. 
