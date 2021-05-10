.. index:: compute energy/gran

compute energy/gran command
===========================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID energy/gran term_1 term_2 ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* energy/gran = style name of this compute command
* input = one or more of the following energy terms:
  
  .. parsed-literal::
  
          friction = energy dissipated by friction
          rotational_kinetic = rotational kinetic energy
          translational_kinetic = translational kinetic energy
          kinetic = total kinetic energy
          boundary = total energy added/removed by boundary work
          volumetric = energy added/removed by volumetric work
          distortional = energy added/removed by distortional work
          normal_strain = normal spring contribution to strain energy
          shear_strain = shear spring contribution to strain energy
          strain = total strain energy
          local_damping = energy dissipated by local damping
          viscous_damping = energy dissipated by viscous damping
          damping = total energy dissipated by damping

Examples
""""""""

.. parsed-literal::

   compute 1 all energy/gran translational_kinetic friction kinetic rotational_kinetic normal_strain shear_strain boundary strain damping
   compute 1 all energy/gran strain friction kinetic

Description
"""""""""""

Define a computation that allows various specified components of
energy to be written out periodically during a simulation. This
compute requires a granular pairstyle to be defined. When *boundary*
and/or *volumetric* and/or *distortional* work terms are listed
among the energy terms, a helper fix, FixEnergyBoundary,
is implicitly created in the code.

WARNING: SINCE EXTRA PER-CONTACT QUANTITIES ARE STORED WHEN ENERGY
TRACING IS ACTIVE, ENERGY TRACING CANNOT BE ENABLED PART-WAY THROUGH
A SIMULATION; IT MUST BE ACTIVE FROM THE START OF THE SIMULATION
(TYPICALLY A CLOUD OF NON-CONTACTING PARTICLES).

Distortional and volumetric work are subdivisions of the total
boundary work. See, for example, Eq. 1.33 (p.21) of 'Soil Behaviour
and Critical State Soil Mechanics' by David Muir Wood (1990).

At present, this compute is applicable to all granular pairstyles and
boundary conditions with one exception. When rigid walls are present 
AND they are moving (an unusual combination), *boundary*\ , *volumetric*
and *distortional* are not currently 
calculated. If this capability is needed, the framework and all of the
hooks are already present in FixEnergyBoundary so it should be reasonably
easy to add. Strain energy for particle-wall overlaps is seamlessly
integrated with the strain energy for particle-particle overlaps.

The input energy terms can be listed in any order. A term cannot be
duplicated (but there is no reason to do this anyway).

Output info
"""""""""""

This compute provides three types of output. Firstly, the specified total
energy terms for the whole system are available to be written to a file
using :doc:`fix\_print <fix_print>`. For example:


.. parsed-literal::

   	compute		2 all energy/gran translational_kinetic friction kinetic normal_strain shear_strain boundary

   	variable 	step equal step
   	variable	tke equal c_2[1]
   	variable	friction equal c_2[2]
   	variable	ke equal c_2[3]
   	variable	nstr equal c_2[4]
   	variable	sstr equal c_2[5]
   	variable	boundary equal c_2[6]

   	fix 		5 all print 1000 "${step} ${friction} ${boundary} ${tke} ${ke} ${nstr} ${sstr}" file Energy_Terms.txt screen no

It is noted that the order of the variables in :doc:`fix\_print <fix_print>` does not
need to correspond to the order in the compute. The second type of
output is a per-atom vector/array of kinetic energies. If one or more
kinetic energy terms (*rotational\_kinetic*, *translational\_kinetic* or
*kinetic*\ ) are defined in the compute, these can be accessed using 
:doc:`dump\_custom <dump>`. The order of the per-particle kinetic energies is the 
same as their order among the compute arguments. In the example shown
above, two kinetic energies are defined: (1) translational\_kinetic; 
(2) kinetic. These can be written to dump files using the following
command:

.. parsed-literal::

   	dump		1 all custom 100000 dump.ke_\* c_2[1] c_2[2]

Finally, per-contact energy terms are available via :doc:`dump\_local <dump>` if the
*trace\_energy* option of :doc:`pair\_modify <pair_modify>` is enabled in the script.
As described
in `Interpreting\_ComputePairLocal\_Output.doc <USER/gran/Interpreting_ComputePairLocal_Output.doc>`_, 
the last four slots are 
available for energy components, though the last of these is currently 
unused. Continuing this example, write out the atom tags and three defined
energy components using the following syntax:

.. parsed-literal::

           pair_modify     trace_energy
   	compute 	1 all pair/local p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18
   	dump		2 all local 100000 dump.per_contact_energy_\* c_1[5] c_1[6] c_1[15] c_1[16] c_1[17]

A simple input script, in.compute\_energy\_gran, which demonstrates all of these output 
options and a MATLAB plotting script are provided with the other documentation.

Restrictions
""""""""""""
 A granular pairstyle must be defined.

Default
"""""""
none
