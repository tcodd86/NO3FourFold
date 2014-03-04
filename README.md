NO3FourFold
===========

FourFold model for NO3 spectral Analysis

This model should be compiled to a .dll for use with the SpecView software package available at https://molspect.chemistry.ohio-state.edu/goes/GOESsoftware.html

Functions in this model are called by SpecView to build the rotational Hamiltonian and then diagonalize and plot a spectrum. SpecView then allows for transitions in the simulated spectrum to be assigned to experimental spectra and then fit rotational constants to reproduce experimental data.

This model in particular is for use in electronic states with Jahn-Teller coupling and a threefold axis of symmetry. These molecules will have their energy levels split into 3 symmetry components, a1, a2, and e for a total of four states (e = e_plus + e_minus). The assumption here is that all symmetries are either symmetric (') or antisymmetric (''). Interactions between ' and '' levels are only allowed at higher order interactions.

Coupling elements included are Jahn-Teller, Spin-Orbit, Spin-Rotation, and Coriolis. These are derived and listed in a to be published paper by Dmitry Melnik. The formula numbers listed next to the matrix elements in the model correspond to this document. A reference will be added once it is published. For explanations of the functions in this code download the SpecView software and see the chapter in the instruction manual about writing new models.

This model includes as a subset an oblate symmetric top Hamiltonian with spin rotation and centrifugal distortion. The oblate symmetric top model was written by Ming-Wei Chen as part of his doctoral research.
