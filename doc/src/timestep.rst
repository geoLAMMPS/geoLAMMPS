.. index:: timestep

timestep command
================

Syntax
""""""

.. parsed-literal::

   timestep dt

* dt = timestep size (time units) or auto

Examples
""""""""

.. code-block:: LAMMPS

   timestep 2.0
   timestep 0.003
   timestep auto

Description
"""""""""""

Set the timestep size for subsequent molecular dynamics simulations.
See the :doc:`units <units>` command for the time units associated with
each choice of units that LAMMPS supports.

The default value for the timestep size also depends on the choice of
units for the simulation; see the default values below.

An automatic timestep calculation capability has been added for
granular simulations where shear history is stored. For the Hertzian 
contact models, the stiffnesses are calculated based on the assumption
of a 5% overlap, so are somewhat conservative. Furthermore the
multiplicative factor is 0.1, so the timestep is 0.1*sqrt(m/k): again
this is quite conservative. These choices are hard-coded in the function
Input::auto_timestep, so the user can change these figures if so desired.

When the :doc:`run style <run_style>` is *respa*, dt is the timestep for
the outer loop (largest) timestep.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix dt/reset <fix_dt_reset>`, :doc:`run <run>`,
:doc:`run_style <run_style>` respa, :doc:`units <units>`

Default
"""""""

+--------------------------------+---------------+-----------------------+
| choice of :doc:`units <units>` | time units    | default timestep size |
+--------------------------------+---------------+-----------------------+
| lj                             | :math:`\tau`  | 0.005 :math:`\tau`    |
+--------------------------------+---------------+-----------------------+
| real                           | fs            | 1.0 fs                |
+--------------------------------+---------------+-----------------------+
| metal                          | ps            | 0.001 ps              |
+--------------------------------+---------------+-----------------------+
| si                             | s             | 1.0e-8 s (10 ns)      |
+--------------------------------+---------------+-----------------------+
| cgs                            | s             | 1.0e-8 s (10 ns)      |
+--------------------------------+---------------+-----------------------+
| electron                       | fs            | 0.001 fs              |
+--------------------------------+---------------+-----------------------+
| micro                          | :math:`\mu`\ s| 2.0 :math:`\mu`\ s    |
+--------------------------------+---------------+-----------------------+
| nano                           | ns            | 0.00045 ns            |
+--------------------------------+---------------+-----------------------+
