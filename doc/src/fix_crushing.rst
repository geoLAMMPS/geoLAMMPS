.. index:: fix crushing

fix crushing command
====================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID crushing outputflag seed m sigma0 d0 a b chi alpha redtype *reduction* constante *commlimit* *m2* *sigma02* *d02* *a2* *b2* *reallocate*

* ID, group-ID are documented in :doc:`fix <fix>` command
* crushing = style name of this fix command
* outputflag = a flag to indicate whether or not detailed crushing output is displayed: set this at 1 to display the output; 0 to suppress it
* seed = a positive integer which is used as the seed of the random number generator
* m = Weibull modulus for allocating strengths
* sigma0 = characteristic strength at which 37% of particles of diameter d0 survive
* d0 = nominal particle diameter in m
* a = slope of linear trendline of survival probability vs. normalised characteristic stress
* b = y-intercept of linear trendline of survival probability vs. normalised characteristic stress
* chi = parameter in brittle failure criterion of Christensen
* alpha = multiplicative factor by which to reduce compressive stress to give correspondence with experimental data
* redtype = an integer flag equal to 1 if particle diameters are reduced to just lose contact with neighbours; else 0
* *reduction* = if redtype = 0, the fractional reduction in radius upon failure between 0 and 1
* constante = an integer flag equal to 1 if void ratio is held constant; else 0
* *commlimit* = (optional): a comminution limit on radii in m
* *m2* = (optional): Weibull modulus for > 1 breakage
* *sigma02* = (optional): sigma0 for particles of diameter d02 after > 1 breakage
* *d02* = (optional): Nominal particle diameter in m for > 1 breakage
* *a2* = (optional): a for > 1 breakage
* *b2* = (optional): b for > 1 breakage
* *reallocate* = (optional): keyword meaning that strengths should be reallocated to particles


Examples
""""""""

.. parsed-literal::

   fix 1 all crushing 1 73547 1.0 2.5e7 1.29e-3 -0.5882 1.0 100 0.001 1 0 6e-5

Description
"""""""""""

This fix allows particle crushing during a simulation. The two-parameter failure 
criterion for brittle materials proposed by Christensen (2000) was chosen, which
is well described by Russell and Muir Wood (2009). In order to use this failure
criterion, two parameters are needed: chi and kappa. The *chi* parameter (to which
the model is fairly insensitive) is one of the input parameters to this fix. The
kappa parameter may be found using *m*\ , *sigma0*\ , *d0*\ , *a* and *b*\ . Experimental 
values of these parameters for a real sand are given by Nakata et al. (1999).

*alpha* is a fudge which allows the strengths to be reduced by a specified amount.
This was added because the particles tend to be far too strong (negligible failure
occurring at very high stresses) without this arbitrary reduction of strength. The
reasons for the excessively high strengths are discussed by Hanley et al. (2014). If
you do not wish to use this *alpha* reduction, set it to 1.0 to disable it.

There are two available options for the amount to reduce particle diameter
upon failure which can be selected using *redtype*\ . If *redtype* = 1, this means
that particle diameters are reduced upon failure so that they just lose contact 
with all surrounding particles. If *redtype* = 0, the following argument in
FixCrushing specifies the fraction by which particle diameters are reduced.

The specification of a comminution limit, *commlimit*\ , is optional but highly
recommended. Note that the timestep calculated using *'timestep auto'* is not
recalculated during a simulation. Therefore, if very small particles are added
to a simulation, stability is no longer assured.

It is mandatory to specify one set of *m*\ , *sigma0*\ , *d0*\ , *a* and *b*\ . Optionally,
a second set can be specified. If this is done, the first set of parameters is used
to allocate particle strengths before initial failure, while the second set of
parameters is used to allocated particle strengths after the first
failure. The logic is that the first failure might be attributed to asperity
breakage while the second and subsequent failures might be due to gross splitting.
If the second set of parameters is not specified, the first set of parameters is
used in all cases.

By default, the allocated particle strengths are stored in the restart file and
reallocated to the particles when the simulation is resumed from a restart file.
It may be desired to change the strengths part-way through a simulation, e.g., 
if the Weibull modulus is changed. The *reallocate* specifier is used to ignore
the strengths in the restart file and reallocate strengths based on the current
FixCrushing arguments.

If the multi neighbor style is used, it is important to know that this fix
inserts atoms with the same atom type and density as another atom on the
processor. Since you are inserting small atoms, you will probably want to
associate the newly created atoms with a specific atom type if using multi.
Therefore, you should consider modifying the relevant code in the
insert\_particles function to meet your needs.

Output info
"""""""""""

This fix writes detailed information about each particle failure to the screen 
if *outputflag* is set at 1; this detailed output is suppressed if 
*outputflag* is 0.

In addition, three per-particle quantities are available to write to dump
files using :doc:`dump\_custom <dump>`: (1) the uniaxial compressive
strength, (2) the uniaxial tensile strength and (3) the number of failures
experienced by the particle. These can be written to dump files using the 
following syntax:

.. parsed-literal::

   	fix 	5 all crushing 1 73547 1.0 2.5e7 1.29e-3 -0.5882 1.0 100 0.001 1 0 6e-5
   	dump	1 all custom 100000 dump.crushing_\* f_5[1] f_5[2] f_5[3]

Restrictions
""""""""""""

\* Fix crushing requires the use of the sphere atom\_style.

\* It is defined only for 3D simulations.

References
""""""""""

* Christensen, R. M. (2000). Yield functions, damage states, and intrinsic strength. Mathematics and Mechanics of Solids, 5, 285-300.
* Hanley, K.J., O'Sullivan, C. & Huang, X. (2014). Investigation of Christensen's two-parameter failure criterion for brittle materials. International Symposium on Geomechanics from Micro to Macro (TC105: IS-Cambridge 2014), Cambridge, UK.
* Nakata, Y., Hyde, A. F. L., Hyodo, M. & Murata, H. (1999). A probabilistic approach to sand particle crushing in the triaxial test. Geotechnique, 49(5), 567-583.
* Russell, A. R. & Muir Wood, D. (2009). Point load tests and strength measurements for brittle spheres. International Journal of Rock Mechanics & Mining Sciences, 46, 272-280.
