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

#ifdef FIX_CLASS
// clang-format off
FixStyle(crushing,FixCrushing);
// clang-format on
#else

#ifndef LMP_FIX_CRUSHING_H
#define LMP_FIX_CRUSHING_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCrushing : public Fix {
  friend class FixMultistress;

public:
  FixCrushing(class LAMMPS *, int, char **);
  ~FixCrushing();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void set_weibull_parameters(int);
  double strength_calculation(int,double,int);
  void setup(int);
  void pre_force(int);
  void end_of_step();
  double failure_occurs(int, double **, int);
  double reduce_radius(int, double **, int);
  double increase_rattler_diameter(int, double);
  void change_strengths(int,double,int=0); //~ Give default argument
  double insert_particles(int);
  void xyz_random(double *);
  double minimise_overlap(double *,double **,int,double);
  int adjust_position(double *,double **,int,double,double,double,double &);
  void print_optional_info(double **, int);

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int maxsize_restart();
  int size_restart(int);
  void write_restart(FILE *);
  void restart(char *);

 private:
  //~ Data specified by the user
  int seed,redtype,constante,reallocateflag,me,nprocs;
  double m,sigma0,d0,slopea,interceptb;
  double chiplusone,reduction,commlimit,alphafactor;

  int mstressid; //~ The ID of fix_multistress if active
  int displaymessages; //~ Whether or not to display user output
  double voidratio;
  double totalpvolume; //~ The total volume of particles
  double PI;
  double cumulredvolume; //~ The accumulated solid volume change
  double radiusparticletoinsert,volumeparticletoinsert;

  //~ m, sigma0, d0, a and b for first and (optionally) subsequent breakages
  double weibullparams[5][2];

  class NeighList *list; //~ Allow access to standard neighbor list
  class Compute *tstress;
  class Compute *tcompute;
  class RanPark *random;
 
  char *id_stress;
  double **cparams; //~ For storing strengths and numbers of breakages
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
