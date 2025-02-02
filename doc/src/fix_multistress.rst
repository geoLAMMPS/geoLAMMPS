.. index:: fix multistress

fix multistress command
=======================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID multistress N M R keyword

* ID, group-ID are documented in :doc:`fix <fix>` command
* multistress = style name of this fix command
* N = the percentage error in the target stress. See note 1.
* M = the number of steps to run within N% of the target stress before displaying that the target stress has been attained. See note 1.
* R = optional but highly recommended - this is the maximum engineering strain rate permissible for stress-controlled boundaries.
* one or more keyword/value pairs must be appended
  
  .. parsed-literal::
  
     keyword = *x* or *xx* or *y* or *yy* or *z* or *zz* or *xy* or *yx* or *xz* or *zx* or *yz* or *zy*
       *x\** or *y\** or *z\** = the boundary to be controlled. See note 2. Order doesn't matter. One of the following must be placed directly after each boundary specification:
     	*stress* = the target stress for stress control. Followed by ONE number, which is the target stress itself. See note 3.
     	*final* = the final position of the boundaries. Followed by TWO numbers for the orthogonal boundaries and ONE number for the triclinic. See note 4.
     	*delta* = the change in position of the boundaries relative to the starting position. Following numbers as for *final*\ . See note 4.
     	*scale* = a scaling factor by which to move the boundaries. Followed by ONE number for the orthogonal bounding box. See note 5.
     	*vel* = velocity at which the boundaries approach each other (-ve) or move apart (+ve). Followed by ONE number: the velocity. See note 6.
     	*erate* = engineering strain rate at which the boundary moves. Followed by ONE number: the rate. See note 7. 
     	*cyclicmean* = used for cyclic loading of a stress-controlled boundary. Followed by THREE numbers: mean stress, amplitude of sinusoid and period. See note 10. 
     	*cyclicdeviator* = similar to *cyclicmean* except variation is in q rather than mean stress. See note 10.
     	*constantb* = used to maintain a constant intermediate stress ratio. Followed by ONE number, b, and TWO characters: the boundaries. See note 12 for details.
  
  .. parsed-literal::
  
     keyword = *linkvolstress*
       *x y* or *x z* or *y z* or *y x* or *z x* or *z y* = used for volume control. See note 8.
  
  .. parsed-literal::
  
     keyword = *constantp*
       *x y* or *x z* or *y z* or *y x* or *z x* or *z y* = used for constant p' tests. See note 11.

  .. parsed-literal::
  
     keyword = *constantq*
       *x y* or *x z* or *y z* or *y x* or *z x* or *z y* = used for constant q tests. See note 11.

  .. parsed-literal::
  
     keyword = *units*
       *lattice* or *box* = the type of units to use. See note 9.

Examples
""""""""

.. parsed-literal::

   fix 1 all multistress 0.5 1000 0.05 x stress 100000 y stress 50000 z stress 30000 units box
   fix 2 all multistress 1 40000 1.5 x stress 100000 y stress 100000 z erate -0.01 units box
   fix 1 all multistress 1 40000 x erate 0.01 y erate 0.04 z erate -0.03 units box
   fix 3 all multistress 0.5 100000 10 z vel -0.01 linkvolstress x y units box
   fix 3 all multistress 0.5 100000 10 x stress 100000 y constantb 0.5 z x z erate -0.01 units box
   fix 3 all multistress 1 50000 1 z erate -0.01 constantp x y units box
   fix name all multistress 1 100000 1 x vel -0.5 constantq y z units box

Description
"""""""""""

Apply any reasonable permutation of strain and stress control to a
sample bounded by periodic boundaries, whether in an orthogonal or
triclinic bounding box. Can also be used to simulate an undrained
test by maintaining a constant volume. More on this in Note 7. This
fix is now also capable of doing constant p' and constant q tests, if this is useful to
you: details are given in Note 11. The capability to maintain a constant
intermediate stress ratio is also available.

If you want only simple strain control of a few periodic boundaries,
you would be better off choosing fix deform instead of fix multistress: it's more efficient
as it does not need to include lots of stress-control code and is not redefined
on every timestep.

Because fix multistress generally is used for controlling the stresses on at least
one boundary, it automatically sets up a stress compute (compute\_stress\_atom). The facility
has been included for the user to access this compute and write out mean stresses on the boundaries.
Furthermore, fix multistress allows the per-atom particle stresses to be dumped out. The name of the
stress compute is "{fix multistress ID}\_stress". As an example, the following lines define a
fix multistress, write out the mean stresses in a file, and dump out the per-atom stresses in another
file:

.. parsed-literal::

   fix		7 all multistress 0.5 100000 0.1 z vel -1e-3 linkvolstress x y units box

   variable 	step equal step
   variable 	xxstress equal f_7[1]
   variable 	yystress equal f_7[2]
   variable 	zzstress equal f_7[3]
   variable 	xystress equal f_7[4]
   variable 	xzstress equal f_7[5]
   variable 	yzstress equal f_7[6]
   fix 		8 all print 1000 "${step} ${xxstress} ${yystress} ${zzstress} ${xystress} ${xzstress} ${yzstress}" file mean_stresses.txt screen no

   dump		1 all custom 1000 dump.atom_stresses_$i id c_7_stress[1] c_7_stress[2] c_7_stress[3] c_7_stress[4] c_7_stress[5] c_7_stress[6]

Currently the maximum engineering strain rate, if defined, is linked to the proportional gain used in the stress control
algorithm. The gain is calculated as twice the maximum strain rate divided by the target stress. This has some logical
justification. Suppose that you are isotropically compressing a sample and the stresses are initially negligible (noting
that these stresses are all mean stresses, not the per-particle stresses).
The strain rate is calculated as the gain times the error (difference between target stress and actual calculated stress).
If the error is initially approximately equal to the target stress, the calculated strain rate would be twice the defined 
maximum strain rate. However the strain rate limitation is then imposed, so the actual strain rate would be restricted to the 
maximum strain rate. This remains true until the error becomes half the target stress. Thereafter the strain rate becomes less
than the maximum strain rate, and tends towards zero as the error does likewise. In effect, the boundary is moving under
strain control, at the defined maximum strain rate, until the mean stress on the boundary is halfway towards its target. Then
the strain rate drops gradually as the calculated stress approaches the target, i.e., the strain control mechanism 'takes over'.

I can't think of a better way to set the gain (and it avoids having a 'magic number' either in the code or in the input specification).
If you do, please let me know! If you wish to change this behaviour, the relevant code is given in the calc\_ctrl\_params function
in fix\_multistress.cpp, and looks like the following line:

.. parsed-literal::

   if (maxrateflag == 1) Kp[i] = 2.0\*maxrate[i]/starget[i];

**General Warnings:**

The following features are implemented but have not been thoroughly tested. Use at your own risk and recognise that
there may be bugs in these features:

  \*  Lattice units (which are the default for consistency with fix deform...) 

\*  All of the triclinic parts of the code ESPECIALLY the pre\_exchange flip flag bit

*\ **Note 1:**\ *

N and M are basically legacy remnants of a piece of code that I ended up deciding not
to use. Don't worry about their values. All they do is display an on-screen message to the
user as soon as M steps have been run while the percentage error on all stress-controlled boundaries
is less than N. Not important. (For interest, originally I had planned to automatically stop
the simulation when this criterion had been met, but decided against this as it could potentially be dangerous).
N is required to be positive, while M is required to be a positive integer.

*\ **Note 2:**\ *

Some of these are synonyms. The following pairs do the same thing: *x* and *xx*\ ; *y* and *yy*\ ; *z* and *zz*\ ;
*xy* and *yx*\ ; *xz* and *zx*\ ; *yz* and *zy*\ .

*\ **Note 3:**\ *

The stress specifier is used for stress control. If you are using standard SI units
(which I strongly recommend), the following number will be in Pa. Therefore

.. parsed-literal::

   x stress 100000

will control the stress on the x boundary at 100 kPa. The stress is a signed quantity: positive stresses are compressive.
An error will be returned if a negative stress is entered - there is no mechanism in the DEM code at present to generate
tensile stresses. Once a bonding model has been implemented, this might be reassessed.

*\ **Note 4:**\ *

*final* and *delta* allow the user to specify the exact positions of the boundaries, or alternatively the
change in positions of the boundaries from the initial position (initial == when fix multistress was
run in the script). Signed quantities, so make sure that the signs are correct.

For the keywords *x*\ , *y* and *z* (and synonyms), each direction has a pair of boundaries, e.g., the top and bottom for
the (vertical) z direction. Therefore you need to specify the final positions of both boundaries (for *final*\ ) or
the distance by which each boundary must be moved from the initial position (for *delta*\ ). It is important to note that
this specification is misleading (sorry). Historically, an early version of this fix allowed both boundaries to move
as you might expect. In the current version, the lower boundaries do not move regardless of the user input. The strain
rate will be correct though and because the sample is in a periodic space, it is equivalent to moving both boundaries.
For example, if the cell has a height of 1 m and you specify delta as (-0.1 0.1), the cell will have final dimensions
of 1.2 m as the upper boundary will move by 0.2 m.

For the keywords *xy*\ , *xz* and *yz* (and synonyms), only one number needs to be specified which refers to the final value
of or change in tilt factor (for *final* and for *delta*\ , respectively).

*\ **Note 5:**\ *

*scale* is defined only for the orthogonal (standard cuboid) bounding box. The behaviour is the same as for
fix deform, so look at the explanation given in :doc:`fix\_deform <fix_deform>`.

*\ **Note 6:**\ *

Straightforward enough. Simply sets the velocity of the boundaries. Signed, so be careful with signs. If the
velocity is negative, the upper boundary approaches the stationary lower boundary, while if the velocity is
positive, the boundaries move further apart. For triclinic keywords (\ *xy* etc.), *vel* sets the velocity at which
the tilt factor changes.

fix multistress allows you to specify a strain rate in five ways: *final*\ , *delta*\ , *vel*\ , *erate* and *scale*\ . 
The one that I recommend is *vel*\ , i.e., specify the velocity as this is fixed. All of the others require changing
if you declare fix multistress within a loop, or if you read in a restart file, or redefine the fix for some other reason.
To demonstrate this, suppose we have a periodic cell that is 1 m high (keeping the maths easy...) and you want a strain
rate of -0.2/s so you write 'z erate -0.2'. Run for 1 second of simulation time - the cell will now have a height of 0.8 m. 
Suppose now you redefine fix multistress and keep 'z erate -0.2' the same. Run for another second, and the cell height becomes
0.64 m. In effect, you have inadvertently introduced a discontinuity in velocity: the velocity for the first second was -0.2 m/s
and for the second second (don't often write that...) it was -0.16 m/s.

This is NOT a bug - from a technical standpoint, the strain rate **is** constant, but the dimension of the box has changed between
definitions of fix multistress. What you would need to do to keep the velocity the same is to change the fix multistress
specification to contain 'z erate -0.25'. This is annoying, and a similar problem occurs for every strain rate specification other than *vel*\ .

*\ **Note 7:**\ *

*erate*\ , the engineering strain rate is quite straightforward as well. Signed, and be careful when fix multistress is
redefined for the reasons discussed in Note 6. Note that this strain rate is defined based on the initial dimensions of
the bounding box.

*\ **Note 8:**\ *

*linkvolstress* is used for controlling the volume, i.e., run an 'undrained' simulation. In general, you will be using an
orthogonal bounding box, and will impose strain control on one boundary. Suppose you are compressing vertically (z direction).
You can use *linkvolstress* (short for 'link together for volume and stress') in the x and y directions as follows:

.. parsed-literal::

   z vel -0.01 linkvolstress x y

The *linkvolstress* option does 2 things. 1) It maintains a constant volume and 2) it keeps the stresses in the two specified
directions equal. In the example, the volume will be held constant and the stresses in the x and y directions will be the same.
This volume control should be very precise: the initial volume enclosed by the periodic boundaries is stored externally (in
domain->initialvolume) and if fix multistress is redefined, the volume is read back in from this variable to avoid 'drift'
of the value.

Obviously, you are not allowed to include a dimension in *linkvolstress* and use the same dimension with another keyword, e.g.,

.. parsed-literal::

   x vel -0.01 linkvolstress x y

makes no sense and you'll get a well-deserved error.

*\ **Note 9:**\ *

The default units are lattice units to maintain consistency with fix deform; however these have not been
thoroughly tested (see the warning above) and you will get a polite on-screen warning if you choose this option.
Therefore it is recommended that you use box units, i.e., put 'units box' at the end of your fix multistress
specification.

*\ **Note 10:**\ *

These two options both allow cyclic loading of the sample, but in different ways. For *cyclicmean*\ , the mean
stress on the specified boundary is cycled between the limits *mean\_stress* +/- *amplitude\_of\_sinusoid* in a
sinusoidal manner with a defined period. For example,

.. parsed-literal::

   x cyclicmean 150000 10000 100

would cause the mean stress on the x boundary to fluctuate between 140 and 160 kPa with a period of 100 s
(assuming SI units).

For *cyclicdeviator*\ , the deviatoric stress is cycled sinusoidally. This option only makes sense if
*linkvolstress* is active on the other two boundaries. For example,

.. parsed-literal::

   x cyclicdeviator 50000 10000 100 linkvolstress y z

would cause q to be maintained at a mean value of 50 kPa (the major principal stress is in the
x direction) with a sinusoidal variation of 10 kPa (period of 100 s) while the intermediate and
minor principal stresses are kept equal and the volume is maintained constant. This option is very
specific to soil mechanics!

*\ **Note 11:**\ *

*constantp* has the same syntax as *linkvolstress*\ , but instead of controlling the volume, this option allows you
to run a simulation at constant mean effective stress, p'. The *constantp* option does two things: 1) it maintains a
constant p' and 2) it keeps the stresses in the two specified directions equal. The initial mean effective stress
is stored externally (in domain->meaneffectivestress) and is also stored in the restart file so accuracy should be
quite good.

Note that it is quite easy to cause premature termination of the simulation if this option is used and the third
(usually strain-controlled) dimension is not controlled very carefully. This option aims to maintain a constant
p', and if this becomes impossible, an error is issued and the simulation exits. This might happen if, for example,
you begin shearing in the -z direction from an isotropic stress state of 100 kPa. If sigma\_3' exceeds 300 
kPa, the simulation will stop as the other two stresses cannot become negative to maintain p' constant.

**If** you want to modify this behaviour so that the simulation does not exit but continues anyway, search for
the following error message in fix\_multistress.cpp:

.. parsed-literal::

   error->all(FLERR,"The stress on the strain-controlled boundary is too high to allow p' to be maintained at its desired value");

...change the word 'all' to 'warning' and recompile. The *constantp* option will still be active though, so strange things
might happen with the strain rates of the two controlled boundaries.

*constantq* works in essentially an identical way, except the deviator stress is held constant rather than p'. As
you might expect, the *constantq* option does two things: 1) it maintains a constant q and 2) it keeps the stresses in the two 
specified directions equal. The initial deviator stress is stored externally (in domain->deviatorstress) and also in the restart file.

*\ **Note 12:**\ *

*constantb* has a similar syntax to *cyclicmean* or *cyclicdeviator* and allows the intermediate stress ratio to be maintained
at a constant value. It takes three mandatory arguments. The first is the b value which must be >= 0.0 and <= 1.0. The second and third
arguments are two of *x*\ , *y* and *z*\ . The idea is that the target stress for boundary i is calculated using the stresses for boundaries
j and k:

.. parsed-literal::

   b\*(stress for boundary j)+(1-b)\*(stress for boundary k)

The second and third arguments give the order of j and k, i.e., *x constantb 0.5 y z* assigns the mean stress in the y direction to
j in the equation above. Obviously having both of these arguments is unnecessary, but a total of three arguments is useful for 
consistency with *cyclicmean* and *cyclicdeviator*\ .

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix\_modify <fix_modify>` options
are relevant to this fix.

This fix produces lots of output information indirectly as it sets up an implicit
stress compute. This is described, and an example is provided, above.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

\* Fix multistress requires the use of the sphere atom\_style.

\* Boundaries that are controlled by fix multistress must be periodic.

\* The *linkvolstress*\ , *constantp*\ , *constantq* and *constantb* options are not defined for 2D simulations.

Defaults
""""""""

The option defaults are *units* = *lattice*\ , no maximum strain rate and no volume/p'/q control.
