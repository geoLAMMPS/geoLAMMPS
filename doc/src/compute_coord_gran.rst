.. index:: compute coord/gran

compute coord/gran command
==========================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID coord/gran type1 type2 ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* coord/gran = style name of this compute command
* typeN = atom type for Nth coordination count (see asterisk form below)

Examples
""""""""

.. parsed-literal::

   compute 1 all coord/gran
   compute 1 all coord/gran 1 2
   compute 1 all coord/gran 2\*4 5\*8 \*

Description
"""""""""""

Define a computation that calculates one or more coordination numbers
for each atom in a group.

This compute is a modification of :doc:`compute coord/atom <compute_coord_atom>`.
This compute had defined the coordination number for a specific atom type as
"the number of neighbor atoms within the specified cutoff distance from the 
central atom". This differs from the standard geomechanics definition of
coordination number: the average number of contacts per particle. Therefore, 
this compute was created to comply with the geomechanics definition of
coordination number.

Strictly, the name is a misnomer; this compute returns the connectivity,
not the coordination number. However this is more useful, and it is
trivial to find the coordination number (or variants like the mechanical
coordination number) as the mean of all, or the appropriate subset, of the 
output connectivities.

The *typeN* keywords allow you to specify which atom types contribute
to each coordination number.  One coordination number is computed for
each of the *typeN* keywords listed.  If no *typeN* keywords are
listed, a single coordination number is calculated, which includes
atoms of all types (same as the "\*" format, see below).

The *typeN* keywords can be specified in one of two ways.  An explicit
numeric value can be used, as in the 2nd example above.  Or a
wild-card asterisk can be used to specify a range of atom types.  This
takes the form "\*" or "\*n" or "n\*" or "m\*n".  If N = the number of
atom types, then an asterisk with no numeric values means all types
from 1 to N.  A leading asterisk means all types from 1 to n
(inclusive).  A trailing asterisk means all types from n to N
(inclusive).  A middle asterisk means all types from m to n
(inclusive).

From a geomechanics perspective, this capability allows the user to
output connectivities only for specific ranges of particle size (for
example), so that separate connectivities could be output for fine and
coarse particles.

The value of all coordination numbers will be 0.0 for atoms not in the
specified compute group.

The neighbor list needed to compute this quantity is constructed each
time the calculation is performed (i.e. each time a snapshot of atoms
is dumped).  Thus it can be inefficient to compute/dump this quantity
too frequently.

Output info
"""""""""""

If single *type1* keyword is specified (or if none are specified),
this compute calculates a per-atom vector.  If multiple *typeN*
keywords are specified, this compute calculates a per-atom array, with
N columns.  These values can be accessed by any command that uses
per-atom values from a compute as input.  See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS output options.

The per-atom vector or array values will be a number >= 0.0, as
explained above.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute coord/atom <compute_coord_atom>`
:doc:`compute cluster/atom <compute_cluster_atom>`

Default
"""""""
none
