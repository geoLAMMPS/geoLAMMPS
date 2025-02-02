/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(gran/hooke/history/oldstyle,PairGranHookeHistoryOldstyle);
// clang-format on
#else

#ifndef LMP_PAIR_GRAN_HOOKE_HISTORY_OLDSTYLE_H
#define LMP_PAIR_GRAN_HOOKE_HISTORY_OLDSTYLE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairGranHookeHistoryOldstyle : public Pair {
 public:
  PairGranHookeHistoryOldstyle(class LAMMPS *);
  virtual ~PairGranHookeHistoryOldstyle();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style(); //~ Made virtual [KH - 23 November 2012]
  double init_one(int, int);
  virtual void write_restart(FILE *); //~ Made virtual
  virtual void read_restart(FILE *); //~ Made virtual
  virtual void write_restart_settings(FILE *); //~ Made virtual
  virtual void read_restart_settings(FILE *); //~ Made virtual
  void reset_dt();
  virtual double single(int, int, int, int, double, double, double, double &);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  void rolling_resistance(int, int, int, int, double, double, double, double, double, double, double, double, double **, double *, double *, double *, double *, double *); //~ Added these two functions [KH - 24 October 2013]
  void Deresiewicz1954_spin(int, int, int, int, double, double **, double *, double *, double &, double &, double *, double &, double &, double &, double &, double, double, double &, double, double); // Added D_spin model [MO - 30 November 2014]
  void add_old_omega_fix();
  double memory_usage();
  double atom2cut(int);
  double radii2cut(double, double);
  void *extract(const char *, int &);
  
 protected:
  double gamman, gammat;
  int dampflag;
  double dt;
  int freeze_group_bit;
  int history;

  int neighprev;
  double *onerad_dynamic, *onerad_frozen;
  double *maxrad_dynamic, *maxrad_frozen;

  double Geq, Poiseq, RMSf, Hp; // Added to extract for wall/gran/oldstyle.cpp [MO - 03 April 2015]
  int Model, THETA1;          // Added to extract for wall/gran/oldstyle.cpp [MO - 12 Sep 2015]

  class FixDummy *fix_dummy;
  class FixNeighHistory *fix_history;

  // storage of rigid body masses for use in granular interactions

  class Fix *fix_rigid;    // ptr to rigid body fix, nullptr if none
  double *mass_rigid;      // rigid mass for owned+ghost atoms
  int nmax;                // allocated size of mass_rigid

  void allocate();

  /*~ Used for adding fix old_omega when rolling resistance model
    is active [KH - 24 October 2013]*/
  class Fix *deffix;

  int lastwarning[2]; //~ Used to control frequencies at which warnings about failures to calculate contact stiffnesses are output in the rolling resistance model [KH - 6 November 2013]

  //~ Add quantities for tracing global energy [KH - 19 February 2014]	
  double dissipfriction, normalstrain, shearstrain, spinenergy;
  double gatheredf, gatheredss, gatheredse; //~ Two more added [KH - 17 October 2014]
  //~~ Two more added for D_spin [MO - 13 November 2014]

  //~ Added the function and ints below [KH - 16 July 2019]
  void transfer_history(double*, double*);
  int *history_transfer_factors;
  int numshearquants;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair granular requires atom attributes radius, rmass

The atom style defined does not have these attributes.

E: Pair granular requires ghost atoms store velocity

Use the comm_modify vel yes command to enable this.

E: Could not find pair fix neigh history ID

UNDOCUMENTED

U: Pair granular with shear history requires newton pair off

This is a current restriction of the implementation of pair
granular styles with history.

U: Could not find pair fix ID

A fix is created internally by the pair style to store shear
history information.  You cannot delete it.

*/
