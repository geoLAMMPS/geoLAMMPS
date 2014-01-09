/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(SHEAR_HISTORY,FixShearHistory)

#else

#ifndef LMP_FIX_SHEAR_HISTORY_H
#define LMP_FIX_SHEAR_HISTORY_H

#include "fix.h"
#include "my_page.h"

namespace LAMMPS_NS {

class FixShearHistory : public Fix {
  friend class Neighbor;
  friend class PairGranHookeHistory;
  friend class PairGranShmHistory; //~ Added this [KH - 23 November 2012]

 public:
  FixShearHistory(class LAMMPS *, int, char **);
  ~FixShearHistory();
  int setmask();
  void init();
  void setup_pre_exchange();
  virtual void pre_exchange();
  void min_setup_pre_exchange();
  void min_pre_exchange();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

 protected:
  int *npartner;                // # of touching partners of each atom
  int **partner;                // tags for the partners

  /*~ The number of shear quantities is not necessarily 3, but can
    be several different values. In here, define arrays with differing
    numbers of shear quantities. There are no space implications as
    memory will be allocated for only one of these later on (using new)
    [KH - 9 January 2014]*/
  double (**shearpartner3)[3]; //~ hooke/history or hertz/history
  double (**shearpartner4)[4]; //~ shm/history
  double (**shearpartner16)[16]; //~ hooke/history or hertz/history with rolling
  double (**shearpartner17)[17]; //~ shm/history with rolling

  int num_quants;               // the number of extra quantities for each partner (i.e. contact) modified GM
  int maxtouch;                 // max # of touching partners for my atoms

  class Pair *pair;
  int *computeflag;             // computeflag in PairGranHookeHistory

  int pgsize,oneatom;           // copy of settings in Neighbor
  MyPage<int> *ipage;           // pages of partner atom IDs

  /*~ As before, the number of shear quantities can vary [KH - 9 January 
    2014]*/
  MyPage<double[3]> *dpage3;     // pages of shear history with partners
  MyPage<double[4]> *dpage4;
  MyPage<double[16]> *dpage16;
  MyPage<double[17]> *dpage17;

  char *stringshearpartner, *stringdpage;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Pair style granular with history requires atoms have IDs

Atoms in the simulation do not have IDs, so history effects
cannot be tracked by the granular pair potential.

E: Too many touching neighbors - boost MAXTOUCH

A granular simulation has too many neighbors touching one atom.  The
MAXTOUCH parameter in fix_shear_history.cpp must be set larger and
LAMMPS must be re-built.

*/
