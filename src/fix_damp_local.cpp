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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_damp_local.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDampLocal::FixDampLocal(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix damp local command");

  alpha = atof(arg[3]);
  if (alpha>1.0 || alpha<0.0) error->all(FLERR,"Please check local damping coefficient - should be between 0 and 1");

}

/* ---------------------------------------------------------------------- */

FixDampLocal::~FixDampLocal()
{
  //delete alpha;
}

/* ---------------------------------------------------------------------- */

int FixDampLocal::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDampLocal::init()
{
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixDampLocal::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixDampLocal::post_force(int vflag)
{
  // apply damping force to atoms in group
  // direction is opposed to velocity vector
  // magnitude is alpha times original force

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      f[i][0] -= alpha*fabs(f[i][0])*signofnum(v[i][0]);
      f[i][1] -= alpha*fabs(f[i][1])*signofnum(v[i][1]);
      f[i][2] -= alpha*fabs(f[i][2])*signofnum(v[i][2]);

      torque[i][0] -= alpha*fabs(torque[i][0])*signofnum(omega[i][0]);
      torque[i][1] -= alpha*fabs(torque[i][1])*signofnum(omega[i][1]);
      torque[i][2] -= alpha*fabs(torque[i][2])*signofnum(omega[i][2]);
    }
}

/* ---------------------------------------------------------------------- */

void FixDampLocal::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixDampLocal::signofnum(double numb)
{
  if (numb > 0.0) return 1.0;
  else if  (numb < 0.0) return -1.0;
  else return 0.0;
}

/* ---------------------------------------------------------------------- */

int FixDampLocal::modify_param(int narg, char **arg)
{
   fprintf(screen, "Local damping modified to %s \n",arg[0]);
   alpha=atof(arg[0]);
   int argsread=1;
   return argsread;

}
