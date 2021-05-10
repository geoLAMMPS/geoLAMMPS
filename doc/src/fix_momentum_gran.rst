.. index:: fix momentum/gran

fix momentum/gran command
=========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID momentum/gran

* ID, group-ID are documented in :doc:`fix <fix>` command
* momentum/gran = style name of this fix command

Examples
""""""""

.. parsed-literal::

   fix 1 all momentum/gran

Description
"""""""""""

Zero the angular momentum of the group of atoms on every
timestep. This was originally developed to facilitate validation
simulations using the approach of Thornton (1979) (based on
face-centred cubic arrays). Usage could not be simpler as there are
no optional arguments.

The :doc:`fix momentum <fix_momentum>` option does weird things
with granular materials, which prompted development of this fix.

It is worth explaining briefly how this fix works, for which it
is necessary to understand how a timestep is ordered using the
velocity-Verlet scheme in LAMMPS.

Initially the velocities are updated using the forces from the previous
timestep and then the positions are updated based on these velocities. Both of
these operations are done within the initial\_integrate function of
fix\_nve\_sphere. The forces are updated by the compute function of the 
pairstyle. Finally the velocities are updated once more using these
updated forces. This second velocity updating step is done by the
final\_integrate function of fix\_nve\_sphere.

This fix zeroes two quantities at different points within the integration loop.
The angular velocities are set to zero immediately before the force computation
(as a pre\_force step). Therefore the angular velocities are zero in the
pairstyle when the forces are calculated.

This force computation updates the torque on the atoms using the shear
component of the forces. The torques are zeroed immediately after the force
computation (as a post\_force step). Therefore when the final\_integrate
function is run to update the velocities, the angular velocities (which are
calculated from the torques) are not changed. The same is true in the
initial\_integrate function on the next timestep as the torque will still
be zero.

The fact that the angular velocities are zeroed on every timestep despite
not being updated means that it would suffice to zero the angular velocities
only on the first timestep, and set the torques to zero on all timesteps.
This is more straightforward though, and the effect on efficiency is
completely negligible.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix\_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`. No parameter of this fix can
be used with the *start/stop* keywords of the :doc:`run <run>` command.
This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Default
"""""""
none
