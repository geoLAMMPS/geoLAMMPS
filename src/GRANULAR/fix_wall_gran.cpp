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

/* ----------------------------------------------------------------------
   Contributing authors: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_wall_gran.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "input.h"
#include "variable.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{XPLANE,YPLANE,ZPLANE,ZCYLINDER};    // XYZ PLANE need to be 0,1,2
enum{HOOKE,HOOKE_HISTORY,HERTZ_HISTORY};

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixWallGran::FixWallGran(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 10) error->all(FLERR,"Illegal fix wall/gran command");

  if (!atom->sphere_flag)
    error->all(FLERR,"Fix wall/gran requires atom style sphere");

  restart_peratom = 1;
  create_attribute = 1;

  vector_flag = 1;
  size_vector = 5;
  global_freq = 1;

  // wall/particle coefficients

  kn = atof(arg[3]);
  if (strcmp(arg[4],"NULL") == 0) kt = kn * 2.0/7.0;
  else kt = atof(arg[4]);

  gamman = atof(arg[5]);
  if (strcmp(arg[6],"NULL") == 0) gammat = 0.5 * gamman;
  else gammat = atof(arg[6]);

  xmu = atof(arg[7]);
  dampflag = atoi(arg[8]);
  if (dampflag == 0) gammat = 0.0;

  if (kn < 0.0 || kt < 0.0 || gamman < 0.0 || gammat < 0.0 ||
      xmu < 0.0 || dampflag < 0 || dampflag > 1)
    error->all(FLERR,"Illegal fix wall/gran command");

  // convert Kn and Kt from pressure units to force/distance^2 if Hertzian

  if (force->pair_match("gran/hertz/history",1)) {
    kn /= force->nktv2p;
    kt /= force->nktv2p;
  }

  // wallstyle args

  int iarg = 9;
  if (strcmp(arg[iarg],"xplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = XPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = atof(arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = atof(arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"yplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = YPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = atof(arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = atof(arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"zplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = ZPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = atof(arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = atof(arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"zcylinder") == 0) {
    if (narg < iarg+2) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = ZCYLINDER;
    lo = hi = 0.0;
    cylradius = atof(arg[iarg+1]);
    iarg += 2;
  }

  // check for trailing keyword/values

  wiggle = 0;
  wshear = 0;
  wtranslate = 0;
  wscontrol = 0;
  ftvarying = 0;
  fstr = NULL;
  velwall[0] = velwall[1] = velwall[2] = 0.0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"wiggle") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix wall/gran command");
      if (strcmp(arg[iarg+1],"x") == 0) axis = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) axis = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) axis = 2;
      else error->all(FLERR,"Illegal fix wall/gran command");
      amplitude = atof(arg[iarg+2]);
      period = atof(arg[iarg+3]);
      wiggle = 1;
      loINI = lo;
      hiINI = hi;
      iarg += 4;
    } else if (strcmp(arg[iarg],"shear") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix wall/gran command");
      if (strcmp(arg[iarg+1],"x") == 0) axis = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) axis = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) axis = 2;
      else error->all(FLERR,"Illegal fix wall/gran command");
      vshear = atof(arg[iarg+2]);
      wshear = 1;
      iarg += 3;
    } else if (strcmp(arg[iarg],"translate") == 0) {
      wtranslate = 1;
      velwall[0] = atof(arg[iarg+1]);
      velwall[1] = atof(arg[iarg+2]);
      velwall[2] = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"stresscontrol") == 0) {
      wscontrol = 1;
      wtranslate = 1;
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        ftvarying = 1;
        int nn = strlen(&arg[iarg+1][2]) + 1;
        fstr = new char[nn];
        strcpy(fstr,&arg[iarg+1][2]);
      } else targetf = atof(arg[iarg+1]);
      gain = atof(arg[iarg+2]);
      if (strcmp(arg[iarg+2],"auto") == 0) error->all(FLERR,"Illegal fix wall/gran command - more coding needed");
      iarg += 3;
    } else error->all(FLERR,"Illegal fix wall/gran command");
  }

  if (wscontrol == 1 && (lo != -BIG && hi != BIG)) error->all(FLERR,"Cannot have both lo and hi walls with stresscontrol"); // put warning message for fix output too?
  if (wallstyle == XPLANE && domain->xperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");
  if (wallstyle == YPLANE && domain->yperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");
  if (wallstyle == ZPLANE && domain->zperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");
  if (wallstyle == ZCYLINDER && (domain->xperiodic || domain->yperiodic))
    error->all(FLERR,"Cannot use wall in periodic dimension");

  if (wiggle && wshear)
    error->all(FLERR,"Cannot wiggle and shear fix wall/gran");
  if (wiggle && wallstyle == ZCYLINDER && axis != 2)
    error->all(FLERR,"Invalid wiggle direction for fix wall/gran");
  if (wshear && wallstyle == XPLANE && axis == 0)
    error->all(FLERR,"Invalid shear direction for fix wall/gran");
  if (wshear && wallstyle == YPLANE && axis == 1)
    error->all(FLERR,"Invalid shear direction for fix wall/gran");
  if (wshear && wallstyle == ZPLANE && axis == 2)
    error->all(FLERR,"Invalid shear direction for fix wall/gran");
  if (wtranslate && (lo != -BIG && hi != BIG))
    error->all(FLERR,"Cannot specify both top and bottom walls and translate for fix wall/gran");
  if (wtranslate && wallstyle == ZCYLINDER)
    error->all(FLERR,"Cannot use translate with cylinder fix wall/gran");
  if (wtranslate && (wiggle || wshear))
    error->all(FLERR,"Cannot translate and wiggle or shear fix wall/gran");

  // setup oscillations

  if (wiggle) omega = 2.0*MY_PI / period;

  // perform initial allocation of atom-based arrays
  // register with Atom class

  shear = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  // initialize as if particle is not touching wall

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    shear[i][0] = shear[i][1] = shear[i][2] = 0.0;

  time_origin = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

FixWallGran::~FixWallGran()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored arrays

  memory->destroy(shear);
  delete [] fstr;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::init()
{
  dt = update->dt;

  // check variables for Ftarget

  if (fstr) {
    fvar = input->variable->find(fstr);
    if (fvar < 0)
      error->all(FLERR,"Variable name for fix wall/gran does not exist");
    if (!input->variable->equalstyle(fvar)) error->all(FLERR,"Variable for fix wall/gran is invalid style");
  }

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // set pairstyle from granular pair style

  if (force->pair_match("gran/hooke",1))
    pairstyle = HOOKE;
  else if (force->pair_match("gran/hooke/history",1))
    pairstyle = HOOKE_HISTORY;
  else if (force->pair_match("gran/hooke/history/omp",1))
    pairstyle = HOOKE_HISTORY;
  else if (force->pair_match("gran/hertz/history",1))
    pairstyle = HERTZ_HISTORY;
  else if (force->pair_match("gran/hertz/history/omp",1))
    pairstyle = HERTZ_HISTORY;
  else error->all(FLERR,"Fix wall/gran is incompatible with Pair style");
}

/* ---------------------------------------------------------------------- */

void FixWallGran::setup(int vflag)
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

void FixWallGran::post_force(int vflag)
{
  double dx,dy,dz,del1,del2,delxy,delr,rsq;
  fwall[0] = fwall[1] = fwall[2] = fwall_all[0] = fwall_all[1] = fwall_all[2] = 0.0;

  // if wiggle or shear, set wall position and velocity accordingly
  // if wtranslate lo and hi track the wall position and velwall is set in the constructor
  if (wiggle) {
    double arg = omega * (update->ntimestep - time_origin) * dt;
    if (wallstyle == axis) {
      lo = loINI + amplitude - amplitude*cos(arg);
      hi = hiINI + amplitude - amplitude*cos(arg);
    }
    velwall[axis] = amplitude*omega*sin(arg);
  } else if (wtranslate || wscontrol) {
      if (wscontrol) velscontrol(); // velocty calculation for stress control
      move_wall(); // move_wall will update hi & lo
  } else if (wshear) velwall[axis] = vshear;

  // loop over all my atoms
  // rsq = distance from wall
  // dx,dy,dz = signed distance from wall
  // for rotating cylinder, reset vwall based on particle position
  // skip atom if not close enough to wall
  //   if wall was set to NULL, it's skipped since lo/hi are infinity
  // compute force and torque on atom if close enough to wall
  //   via wall potential matched to pair potential
  // set shear if pair potential stores history

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  shearupdate = 1;
  if (update->setupflag) shearupdate = 0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      dx = dy = dz = 0.0;

      if (wallstyle == XPLANE) {
        del1 = x[i][0] - lo;
        del2 = hi - x[i][0];
        if (del1 < del2) dx = del1;
        else dx = -del2;
      } else if (wallstyle == YPLANE) {
        del1 = x[i][1] - lo;
        del2 = hi - x[i][1];
        if (del1 < del2) dy = del1;
        else dy = -del2;
      } else if (wallstyle == ZPLANE) {
        del1 = x[i][2] - lo;
        del2 = hi - x[i][2];
        if (del1 < del2) dz = del1;
        else dz = -del2;
      } else if (wallstyle == ZCYLINDER) {
        delxy = sqrt(x[i][0]*x[i][0] + x[i][1]*x[i][1]);
        delr = cylradius - delxy;
        if (delr > radius[i]) dz = cylradius;
        else {
          dx = -delr/delxy * x[i][0];
          dy = -delr/delxy * x[i][1];
          if (wshear && axis != 2) {
            velwall[0] = vshear * x[i][1]/delxy;
            velwall[1] = -vshear * x[i][0]/delxy;
            velwall[2] = 0.0;
          }
        }
      }

      rsq = dx*dx + dy*dy + dz*dz;

      if (rsq > radius[i]*radius[i]) {
        if (pairstyle != HOOKE) {
          shear[i][0] = 0.0;
          shear[i][1] = 0.0;
          shear[i][2] = 0.0;
        }
      } else {
        if (pairstyle == HOOKE)
          hooke(rsq,dx,dy,dz,velwall,v[i],f[i],omega[i],torque[i],
                radius[i],rmass[i]);
        else if (pairstyle == HOOKE_HISTORY)
          hooke_history(rsq,dx,dy,dz,velwall,v[i],f[i],omega[i],torque[i],
                        radius[i],rmass[i],shear[i]);
        else if (pairstyle == HERTZ_HISTORY)
          hertz_history(rsq,dx,dy,dz,velwall,v[i],f[i],omega[i],torque[i],
                        radius[i],rmass[i],shear[i]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallGran::hooke(double rsq, double dx, double dy, double dz,
                        double *vwall, double *v,
                        double *f, double *omega, double *torque,
                        double radius, double mass)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double meff,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,ft,fs1,fs2,fs3,fx,fy,fz,rinv,rsqinv;

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

  // normal component

  vnnr = vr1*dx + vr2*dy + vr3*dz;
  vn1 = dx*vnnr * rsqinv;
  vn2 = dy*vnnr * rsqinv;
  vn3 = dz*vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity - removed
  //wr1 = radius*omega[0] * rinv; //radius should be substituted by r, so wr1==omega[0] etc

  // normal forces = Hookian contact + normal velocity damping

  meff = mass;
  damp = meff*gamman*vnnr*rsqinv;
  ccel = kn*(radius-r)*rinv - damp;

  // relative velocities

  vtr1 = vt1 - dz*omega[1]+dy*omega[2];
  vtr2 = vt2 - dx*omega[2]+dz*omega[0];
  vtr3 = vt3 - dy*omega[0]+dx*omega[1];
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // force normalization

  fn = xmu * fabs(ccel*r);
  fs = meff*gammat*vrel;
  if (vrel != 0.0) ft = MIN(fn,fs) / vrel;
  else ft = 0.0;

  // tangential force due to tangential velocity damping

  fs1 = -ft*vtr1;
  fs2 = -ft*vtr2;
  fs3 = -ft*vtr3;

  // forces & torques

  fx = dx*ccel + fs1;
  fy = dy*ccel + fs2;
  fz = dz*ccel + fs3;

  f[0] += fx;
  f[1] += fy;
  f[2] += fz;

  torque[0] -= dy*fs3 - dz*fs2;
  torque[1] -= dz*fs1 - dx*fs3;
  torque[2] -= dx*fs2 - dy*fs1;

  fwall[0] += fx;
  fwall[1] += fy;
  fwall[2] += fz;

}

/* ---------------------------------------------------------------------- */

void FixWallGran::hooke_history(double rsq, double dx, double dy, double dz,
                                double *vwall, double *v,
                                double *f, double *omega, double *torque,
                                double radius, double mass, double *shear)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double meff,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,fs1,fs2,fs3,fx,fy,fz;
  double shrmag,rsht,rinv,rsqinv;

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

  // normal component

  vnnr = vr1*dx + vr2*dy + vr3*dz;
  vn1 = dx*vnnr * rsqinv;
  vn2 = dy*vnnr * rsqinv;
  vn3 = dz*vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity - removed see above

  // normal forces = Hookian contact + normal velocity damping

  meff = mass;
  damp = meff*gamman*vnnr*rsqinv;
  ccel = kn*(radius-r)*rinv - damp;

  // relative velocities

  vtr1 = vt1 - dz*omega[1]+dy*omega[2];
  vtr2 = vt2 - dx*omega[2]+dz*omega[0];
  vtr3 = vt3 - dy*omega[0]+dx*omega[1];
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // shear history effects

  if (shearupdate) {
    shear[0] += vtr1*dt;
    shear[1] += vtr2*dt;
    shear[2] += vtr3*dt;
  }
  shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2]);

  // rotate shear displacements

  rsht = shear[0]*dx + shear[1]*dy + shear[2]*dz;
  rsht = rsht*rsqinv;
  if (shearupdate) {
    shear[0] -= rsht*dx;
    shear[1] -= rsht*dy;
    shear[2] -= rsht*dz;
  }

  // tangential forces = shear + tangential velocity damping

  fs1 = - (kt*shear[0] + meff*gammat*vtr1);
  fs2 = - (kt*shear[1] + meff*gammat*vtr2);
  fs3 = - (kt*shear[2] + meff*gammat*vtr3);

  // rescale frictional displacements and forces if needed

  fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  fn = xmu * fabs(ccel*r);

  if (fs > fn) {
    if (shrmag != 0.0) {
      shear[0] = (fn/fs) * (shear[0] + meff*gammat*vtr1/kt) -
        meff*gammat*vtr1/kt;
      shear[1] = (fn/fs) * (shear[1] + meff*gammat*vtr2/kt) -
        meff*gammat*vtr2/kt;
      shear[2] = (fn/fs) * (shear[2] + meff*gammat*vtr3/kt) -
        meff*gammat*vtr3/kt;
      fs1 *= fn/fs ;
      fs2 *= fn/fs;
      fs3 *= fn/fs;
    } else fs1 = fs2 = fs3 = 0.0;
  }

  // forces & torques

  fx = dx*ccel + fs1;
  fy = dy*ccel + fs2;
  fz = dz*ccel + fs3;

  f[0] += fx;
  f[1] += fy;
  f[2] += fz;

  torque[0] -= dy*fs3 - dz*fs2;
  torque[1] -= dz*fs1 - dx*fs3;
  torque[2] -= dx*fs2 - dy*fs1;

  fwall[0] += fx;
  fwall[1] += fy;
  fwall[2] += fz;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::hertz_history(double rsq, double dx, double dy, double dz,
                                double *vwall, double *v,
                                double *f, double *omega, double *torque,
                                double radius, double mass, double *shear)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double meff,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,fs1,fs2,fs3,fx,fy,fz;
  double shrmag,rsht,polyhertz,rinv,rsqinv;

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

  // normal component

  vnnr = vr1*dx + vr2*dy + vr3*dz;
  vn1 = dx*vnnr / rsq;
  vn2 = dy*vnnr / rsq;
  vn3 = dz*vnnr / rsq;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity - removed

  // normal forces = Hertzian contact + normal velocity damping

  meff = mass;
  damp = meff*gamman*vnnr*rsqinv;
  ccel = kn*(radius-r)*rinv - damp;
  polyhertz = sqrt((radius-r)*radius);
  ccel *= polyhertz;

  // relative velocities

  vtr1 = vt1 - dz*omega[1]+dy*omega[2];
  vtr2 = vt2 - dx*omega[2]+dz*omega[0];
  vtr3 = vt3 - dy*omega[0]+dx*omega[1];
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // shear history effects

  if (shearupdate) {
    shear[0] += vtr1*dt;
    shear[1] += vtr2*dt;
    shear[2] += vtr3*dt;
  }
  shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2]);

  // rotate shear displacements

  rsht = shear[0]*dx + shear[1]*dy + shear[2]*dz;
  rsht = rsht*rsqinv;
  if (shearupdate) {
    shear[0] -= rsht*dx;
    shear[1] -= rsht*dy;
    shear[2] -= rsht*dz;
  }

  // tangential forces = shear + tangential velocity damping

  fs1 = -polyhertz * (kt*shear[0] + meff*gammat*vtr1);
  fs2 = -polyhertz * (kt*shear[1] + meff*gammat*vtr2);
  fs3 = -polyhertz * (kt*shear[2] + meff*gammat*vtr3);

  // rescale frictional displacements and forces if needed

  fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  fn = xmu * fabs(ccel*r);

  if (fs > fn) {
    if (shrmag != 0.0) {
      shear[0] = (fn/fs) * (shear[0] + meff*gammat*vtr1/kt) -
        meff*gammat*vtr1/kt;
      shear[1] = (fn/fs) * (shear[1] + meff*gammat*vtr2/kt) -
        meff*gammat*vtr2/kt;
      shear[2] = (fn/fs) * (shear[2] + meff*gammat*vtr3/kt) -
        meff*gammat*vtr3/kt;
      fs1 *= fn/fs ;
      fs2 *= fn/fs;
      fs3 *= fn/fs;
    } else fs1 = fs2 = fs3 = 0.0;
  }

  // forces & torques

  fx = dx*ccel + fs1;
  fy = dy*ccel + fs2;
  fz = dz*ccel + fs3;

  f[0] += fx;
  f[1] += fy;
  f[2] += fz;

  torque[0] -= dy*fs3 - dz*fs2;
  torque[1] -= dz*fs1 - dx*fs3;
  torque[2] -= dx*fs2 - dy*fs1;

  fwall[0] += fx;
  fwall[1] += fy;
  fwall[2] += fz;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixWallGran::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += 3*nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixWallGran::grow_arrays(int nmax)
{
  memory->grow(shear,nmax,3,"fix_wall_gran:shear");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixWallGran::copy_arrays(int i, int j, int delflag)
{
  shear[j][0] = shear[i][0];
  shear[j][1] = shear[i][1];
  shear[j][2] = shear[i][2];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixWallGran::set_arrays(int i)
{
  shear[i][0] = shear[i][1] = shear[i][2] = 0.0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixWallGran::pack_exchange(int i, double *buf)
{
  buf[0] = shear[i][0];
  buf[1] = shear[i][1];
  buf[2] = shear[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   unpack values into local atom-based arrays after exchange
------------------------------------------------------------------------- */

int FixWallGran::unpack_exchange(int nlocal, double *buf)
{
  shear[nlocal][0] = buf[0];
  shear[nlocal][1] = buf[1];
  shear[nlocal][2] = buf[2];
  return 3;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixWallGran::pack_restart(int i, double *buf)
{
  int m = 0;
  buf[m++] = 4;
  buf[m++] = shear[i][0];
  buf[m++] = shear[i][1];
  buf[m++] = shear[i][2];
  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixWallGran::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  shear[nlocal][0] = extra[nlocal][m++];
  shear[nlocal][1] = extra[nlocal][m++];
  shear[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixWallGran::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixWallGran::size_restart(int nlocal)
{
  return 4;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- 
Allows the user to do a fix_modify at the input script and change the
parameters of the fix. Only allows a few things to be modified.
Returns the number of arguments read.
------------------------------------------------------------------------- */

int FixWallGran::modify_param(int narg, char **arg)
{
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    int argsread=0;

    if (strcmp(arg[argsread],"gamman") == 0) {// do while loop instead?
      fprintf(screen, "changed wall gamman from %f to ",gamman);
      gamman=atof(arg[argsread+1]);
      argsread+=2;
      fprintf(screen, "%f\n",gamman);
    }
    else if (strcmp(arg[argsread],"gammat") == 0) {
      fprintf(screen, "changed wall gammat from %f to ",gammat);
      gammat=atof(arg[argsread+1]);
      argsread+=2;
      fprintf(screen, "%f\n",gammat);
    }
    else if (strcmp(arg[argsread],"mu") == 0) {
      fprintf(screen, "changed wall friction coefficient from %f to ",xmu);
      xmu=atof(arg[argsread+1]);
      argsread+=2;
      fprintf(screen, "%f\n",xmu);
    }
    else if (strcmp(arg[argsread],"dampflag") == 0) {
      dampflag=atoi(arg[argsread+1]);
      argsread+=2;
      fprintf(screen, "changed wall dampflag to %d \n",dampflag);
    }
    else if (strcmp(arg[argsread],"translate") == 0) {
      if (strcmp(arg[argsread+1],"off") == 0) {
         wtranslate=0;
         velwall[0]=velwall[1]=velwall[2]=0.0;
         argsread+=2;
         fprintf(screen, "stopped wall translation\n");
      } else {
         wtranslate = 1;
         velwall[0] = atof(arg[argsread+1]);
         velwall[1] = atof(arg[argsread+2]);
         velwall[2] = atof(arg[argsread+3]);
         argsread+= 4;
         fprintf(screen, "changed wall velocity to [ %e %e %e ]\n",velwall[0],velwall[1],velwall[2]);
      }
    }
    else {
       fprintf(screen,"Argument %s not yet supported\n",arg[argsread]);
       error->all(FLERR,"Illegal fix modify wall/gran command");
    }
    if (argsread==narg) {
       if (gamman <0.0 || gammat <0.0 || xmu <0.0 ) error->all(FLERR,"Check the values for the modified granular wall parameters");
       if (wtranslate && (lo != -BIG && hi != BIG))
    error->all(FLERR,"Cannot specify both top and bottom walls and translate for fix wall/gran - check your fix_modify"); // this check will fail if the walls have moved..
       if (wtranslate && wallstyle == ZCYLINDER)
    error->all(FLERR,"Cannot use translate with cylinder fix wall/gran - check your fix_modify");
       if (wtranslate && (wiggle || wshear))
    error->all(FLERR,"Cannot translate and wiggle or shear fix wall/gran- check your fix_modify");
    }
    return argsread;
}

/* ---------------------------------------------------------------------- 
A function that implements wall movement
------------------------------------------------------------------------- */

void FixWallGran::move_wall() {

 lo+=velwall[wallstyle]*dt;
 hi+=velwall[wallstyle]*dt;

}

/* ---------------------------------------------------------------------- 
A function that calculates the wall velocity based on a target force and a 
specific controller
------------------------------------------------------------------------- */

void FixWallGran::velscontrol() {

 MPI_Allreduce(fwall,fwall_all,3,MPI_DOUBLE,MPI_SUM,world);
 if (ftvarying == 1) {
   //modify->clearstep_compute();needed???
   targetf = input->variable->compute_equal(fvar);
   //modify->addstep_compute(update->ntimestep + 1);needed???
 }
 velwall[wallstyle] = gain * (targetf - fwall_all[wallstyle]);

}

/* ---------------------------------------------------------------------- 
A function that calculates the outputs of the fix
fid[0]:low wall position
fid[1]:high wall position
fid[2-5]:forces on wall
------------------------------------------------------------------------- */

double FixWallGran::compute_vector(int n)
{

  if (n == 0) return lo;
  if (n == 1) return hi;
  // only sum across procs one time?? //if (eflag == 0) {??
  MPI_Allreduce(fwall,fwall_all,3,MPI_DOUBLE,MPI_SUM,world);
  if (n>4) error->all(FLERR,"Illegal fix_wall_gran output");
  return fwall_all[n-2];
}

