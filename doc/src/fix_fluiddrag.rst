.. index:: fix fluiddrag

fix fluiddrag command
=====================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID fluiddrag i N style magnitude args

* ID, group-ID are documented in :doc:`fix <fix>` command
* fluiddrag = style name of this fix command
* i = hydraulic gradient
* N = *0* or *1*
  
  .. parsed-literal::
  
     *0* = subtract the velocity of the largest particle from
     the velocities of all particles in the system on each timestep
     *1* = do not perform this velocity adjustment

* magnitude = size of acceleration (force/mass units)
* style = *chute* or *spherical* or *gradient* or *vector*
  
  .. parsed-literal::
  
       *chute* args = angle
         angle = angle in +x away from -z or -y axis in 3d/2d (in degrees)
       *spherical* args = phi theta
         phi = azimuthal angle from +x axis (in degrees)
         theta = angle from +z or +y axis in 3d/2d (in degrees)
       *gradient* args = phi theta phi_grad theta_grad
         phi = azimuthal angle from +x axis (in degrees)
         theta = angle from +z or +y axis in 3d/2d (in degrees)
         phi_grad = rate of change of angle phi (full rotations per time unit)
         theta_grad = rate of change of angle theta (full rotations per time unit)
       *vector* args = x y z
         x y z = vector direction to apply the acceleration

Examples
""""""""

.. parsed-literal::

   fix 1 all fluiddrag 5.0 0 vector 1 0 0

Description
"""""""""""

A fix based on :doc:`fix gravity <fix_gravity>` which simulates the flow
of water through a packing of particles using a simple model. It is
important to note that SI units must be used as the specific weight
of water is hard-coded as 9810 kg/(m2.s2) for convenience.

A per-atom array is created by this fix which contains data about the
simulation which may be of interest. This can be written to a file using the
:doc:`dump custom <dump>` command, for example:

.. parsed-literal::

   fix		9 all fluiddrag 1.0 1 vector 1 0 0 
   dump		1 all custom 5000 dump.drag_forces f_9[1] f_9[2] f_9[3] f_9[4] f_9[5] f_9[6] f_9[7] f_9[8] f_9[9] f_9[10] f_9[11]

Of course, asterisks can be used to write separate files at each write interval
if so desired. There are 9 available columns in 2D and 11 in 3D which are ordered as follows:

*\ **2D:**\ *

* **1** The tag of particle *i*
* **2** The x-coordinate of particle *i*
* **3** The y-coordinate of particle *i*
* **4** The radius of particle *i*
* **5** The total volume enclosed by the walls or periodic boundaries
* **6** The x-component of the particle-fluid interaction force applied to particle *i*
* **7** The y-component of the particle-fluid interaction force
* **8** The porosity of the sample
* **9** The drag force (based on effective diameter)
*\ **3D:**\ *

* **1** The tag of particle *i*
* **2** The x-coordinate of particle *i*
* **3** The y-coordinate of particle *i*
* **4** The z-coordinate of particle *i*
* **5** The radius of particle *i*
* **6** The total volume enclosed by the walls or periodic boundaries
* **7** The x-component of the particle-fluid interaction force applied to particle *i*
* **8** The y-component of the particle-fluid interaction force
* **9** The z-component of the particle-fluid interaction force
* **10** The porosity of the sample
* **11** The drag force (based on effective diameter)
**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix\_modify <fix_modify>` options
are relevant to this fix. No parameter of this fix can
be used with the *start/stop* keywords of the :doc:`run <run>` command.
This fix is not invoked during :doc:`energy minimization <minimize>`.
The existence of per-atom quantities has already be described.

Restrictions
""""""""""""

\* Fix fluiddrag requires the use of the sphere atom\_style.

\* SI units must be used at present.

Related commands
""""""""""""""""

:doc:`fix gravity <fix_gravity>`

Default
"""""""
none
