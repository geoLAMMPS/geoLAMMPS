The Thornton_FCC directory contains an analytical verification in two steps.

Firstly, run the simulation contained in Step_1_Isotropic_Compression to create a small isotropically compressed face-centred cubic sample at 150 kPa and with mu = 0.25.

Secondly, copy the last restart file to the Step_2_Triaxial_Shearing directory and run that quick simulation. The major:minor stress ratio tends towards 3.333 as it ought to: analytically the ratio should become 2*(1+mu)/(1-mu) where mu is the interparticle friction coefficient.

The explanation for this is given by Thornton, C. (1979), The conditions for failure of a face-centered cubic array of uniform rigid spheres, Geotechnique, 29(4), 441–459, https://doi.org/10.1680/geot.1979.29.4.441