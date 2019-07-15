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

FixStyle(wall/gran/oldstyle,FixWallGranOldstyle)

#else

#ifndef LMP_FIX_WALL_GRAN_OLDSTYLE_H
#define LMP_FIX_WALL_GRAN_OLDSTYLE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallGranOldstyle : public Fix {
 public:
  FixWallGranOldstyle(class LAMMPS *, int, char **);
  virtual ~FixWallGranOldstyle();
  int setmask();
  virtual void init();
  void setup(int);
  virtual void post_force(int);
  virtual void post_force_respa(int, int, int);

  virtual double memory_usage();
  virtual void grow_arrays(int);
  virtual void copy_arrays(int, int, int);
  virtual void set_arrays(int);
  virtual int pack_exchange(int, double *);
  virtual int unpack_exchange(int, double *);
  virtual int pack_restart(int, double *);
  virtual void unpack_restart(int, int);
  virtual int size_restart(int);
  virtual int maxsize_restart();
  void reset_dt();
  int modify_param(int, char **);
  double compute_vector(int);
  void *extract(const char *, int &); //~ [KH - 20 February 2014]
  void write_restart(FILE *); //~ [KH - 20 February 2014]
  void restart(char *); //~ [KH - 20 February 2014]

  void hooke(double, double, double, double, double *,
             double *, double *, double *, double *, double, double, int);
  void hooke_history(double, double, double, double, double *,
                     double *, double *, double *, double *, double, double,
                     double *, int);
  void hertz_history(double, double, double, double, double *, double,
                     double *, double *, double *, double *, double, double,
                     double *, int);
  void bonded_history(double, double, double, double, double *, double,
                       double *, double *, double *, double *, double, double,
                       double *);
  void shm_history(double, double, double, double, double *, double,
		   double *, double *, double *, double *, double, double,
		   double *, int); //~ [KH - 30 October 2013]
  void CM_history(double, double, double, double, double *,
		    double *, double *, double *, double *, double, double,
		    double *, int); //~ [MO - 18 July 2014]
  void HMD_history(double, double, double, double, double *,
		     double *, double *, double *, double *, double, double,
		     double *, int); //~ [MO - 21 July 2014]
  void CMD_history(double, double, double, double, double *,
		     double *, double *, double *, double *, double, double,
		     double *, int); //~ [MO - 30 November 2014]
  
 protected:
  int sfound; // added to fetch meanstress in stresscontrol [MO - 13 August 2015]
  class Compute *stressatom; // added to fetch meanstress in stresscontrol [MO - 13 August 2015]
  int wallstyle,wiggle,wshear,axis,dampflag;
  double wcoordnos[1],wcoordnos_all[1]; // coordination number of wall [MO - 12 March 2014]
  int wtranslate,wscontrol; //flags for wall movement and wall stress control respectively
  double kn,kt,gamman,gammat,xmu,Geq,Poiseq,RMSf,Hp; // increased for CM & CMD models [MO - 18 July 2014]
  int Model,THETA1; // for CM & CMD models [MO - 12 June 2015] and also for HMD & CMD [MO - 12 Sep 2015]
  int wiggletype; // added to select sin or cos type of wiggle[MO - 09 May 2016]
  double *xmu_p,*Geq_p,*Poiseq_p,*RMSf_p,*Hp_p;    // parameters of the contacting particle [MO - 05 December 2014]

  int pairstyle,nlevels_respa;
  bigint time_origin;
  double E,G,SurfEnergy;
  double lo,hi,cylradius;
  double loINI,hiINI; // for wiggle only
  double vwall[3],fwall[3],fwall_all[3],w_ierates[3];
  double amplitude,period,omega,vshear;  
  double dt;
  double targetf,gain;
  char *fstr;
  int fvar;
  int ftvarying; // 1 if ftarget set through a variable
  int *rolling,*model_type; //~ Quantities for rolling resistance model [KH - 30 October 2013]
  int *D_spin,*D_switch; // Quantities for D_spin model [MO - 30 November 2014]
  double *rolling_delta,*kappa,*post_limit_index;
  int lastwarning[2]; //~ Used to control frequencies at which warnings about failures to calculate contact stiffnesses are output in the rolling resistance model [KH - 6 November 2013]

  //~ Add quantities for tracing global energy [KH - 20 February 2014] also for spinenergy [MO - 30 November 2014]
  int pairenergy, *trace_energy;
  double dissipfriction, normalstrain, shearstrain, spinenergy;

  /*~ Used for accessing fix old_omega when rolling resistance model
    is active [KH - 30 October 2013]*/
  class Fix *deffix;
  
  int *touch;

  char *idregion;

  int history;       // if particle/wall interaction stores history
  int shearupdate;   // flag for whether shear history is updated
  int sheardim;      // # of shear history values per contact

  // shear history for single contact per particle

  double **shearone;

  // rigid body masses for use in granular interactions

  class Fix *fix_rigid;    // ptr to rigid body fix, NULL if none
  double *mass_rigid;      // rigid mass for owned+ghost atoms
  int nmax;                // allocated size of mass_rigid

  void move_wall();
  void velscontrol();
  void ev_tally_wall(int, double, double, double, double, double, double, double);
  void rolling_resistance(int, int, double, double, double, double, double, double, double, double, double *, double *, double *, double *, double *); //~ Added this function [KH - 30 October 2013]
  void Deresiewicz1954_spin(int, int, double, double, double, double, double, double *, double *, double *, double &, double &, double *, double &, double &, double &, double &, double, double, double &, double, double); //~ Added this function [MO - 30 November 2014]
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix wall/gran/oldstyle requires atom style sphere

Self-explanatory.

E: Invalid fix wall/gran/oldstyle interaction style

UNDOCUMENTED

E: Cannot use wall in periodic dimension

Self-explanatory.

E: Cannot wiggle and shear fix wall/gran/oldstyle

Cannot specify both options at the same time.

E: Invalid wiggle direction for fix wall/gran/oldstyle

Self-explanatory.

E: Invalid shear direction for fix wall/gran/oldstyle

Self-explanatory.

E: Cannot wiggle or shear with fix wall/gran/region/oldstyle

UNDOCUMENTED

U: Fix wall/gran/oldstyle is incompatible with Pair style

Must use a granular pair style to define the parameters needed for
this fix.

*/
