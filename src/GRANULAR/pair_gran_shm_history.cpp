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

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "pair_gran_shm_history.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "fix.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "comm.h"
#include "error.h"
#include "domain.h"
#include "modify.h" //~ This header file was added [KH - 9 November 2011]
#include "neigh_request.h" //~ Added four more header files [KH - 23 November 2012]
#include "fix.h"
#include "fix_neigh_history.h"
#include "memory.h"
#include "compute_energy_gran.h" //~ For energy tracing [KH - 20 February 2014]

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGranShmHistory::PairGranShmHistory(LAMMPS *lmp) :
  PairGranHookeHistory(lmp) {}

/* ---------------------------------------------------------------------- */

void PairGranShmHistory::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double ccel,tor1,tor2,tor3;
  double fslim,fs,fs1,fs2,fs3;
  double shsqmag,shsqnew,shratio,rsht,polyhertz;
  double wspinx,wspiny,wspinz,shint0,shint1,shint2,omdel;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int shearupdate = 1;
  if (update->setupflag) shearupdate = 0;

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body

  if (fix_rigid && neighbor->ago == 0) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"pair:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++)
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    comm->forward_comm_pair(this);
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double deltan,cri,crj;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = fix_history->firstflag;
  firstshear = fix_history->firstvalue;

  /*~ The following piece of code was added to determine whether or not
    any periodic boundaries, if present, are moving either via fix_
    multistress or fix_deform [KH - 9 November 2011]*/
  int velmapflag = 0;

  if (domain->xperiodic || domain->yperiodic || domain->zperiodic) {
    for (int q = 0; q < modify->nfix; q++)
      if (strcmp(modify->fix[q]->style,"multistress") == 0) {
	ierates = modify->fix[q]->param_export();
	velmapflag = 1;
	break;
      }
  
    if (velmapflag == 0)
      for (int q = 0; q < modify->nfix; q++) {
	if (strcmp(modify->fix[q]->style,"deform") == 0) {
	  ierates = modify->fix[q]->param_export();
	  velmapflag = 1;
	  break;
	}
      }
  }

  /*~ Ascertain whether or not energy tracing is active by checking
    for the presence of compute energy/gran and also the trace_energy
    flag which indicates per-contact tracing is active. If so, check if the
    tracked terms include those calculated in this pairstyle. The
    energy terms are updated only if necessary for efficiency.
    [KH - 6 March 2014]*/
  int pairenergy = trace_energy;
  if (!pairenergy)
    for (int q = 0; q < modify->ncompute; q++)
      if (strcmp(modify->compute[q]->style,"energy/gran") == 0) {
	pairenergy = ((ComputeEnergyGran *) modify->compute[q])->pairenergy;
	break;
      }

  //~ Initialise the non-accumulated strain energy terms to zero
  normalstrain = 0.0;

  /*~ The number of shear quantities is 19 if rolling resistance
    is active [KH - 29 July 2013]*/
  /*~ Another 4 shear quantities were added for per-contact energy
    tracing [KH - 6 March 2014]. Modified by MO [4 November 2014]*/
  int numshearquants = 4 + 15*rolling + 20*D_spin + 4*trace_energy;

  //~ Use tags to consider contacts only once [KH - 28 February 2014]
  tagint *tag = atom->tag; 

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    touch = firsttouch[i];
    allshear = firstshear[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      if (rsq >= radsum*radsum) {

        // unset non-touching neighbors

        touch[jj] = 0;
        shear = &allshear[numshearquants*jj]; // shear[] refers to a shear force here, not shear displacement
	for (int q = 0; q < numshearquants; q++)
	  shear[q] = 0.0; //~ Added the 'for' loop [KH - 23 October 2013]

      } else {
        r = sqrt(rsq);
        rinv = 1.0/r;
        rsqinv = 1.0/rsq;
        deltan = radsum-r;
	cri = radi-0.5*deltan*rinv*(2.0*radj-deltan);
        crj = radj-0.5*deltan*rinv*(2.0*radi-deltan);

        // relative translational velocity

        vr1 = v[i][0] - v[j][0];
        vr2 = v[i][1] - v[j][1];
        vr3 = v[i][2] - v[j][2];

	/*~ These relative velocity components have to be updated for the
	  periodic boundaries [KH - 14 November 2011]*/
	if (velmapflag == 1) {
	  vr1 += (ierates[0]*delx + ierates[3]*dely + ierates[4]*delz);
	  vr2 += (ierates[3]*delx + ierates[1]*dely + ierates[5]*delz);
	  vr3 += (ierates[4]*delx + ierates[5]*dely + ierates[2]*delz);
	}

        // normal component

        vnnr = vr1*delx + vr2*dely + vr3*delz;
        vn1 = delx*vnnr * rsqinv;
        vn2 = dely*vnnr * rsqinv;
        vn3 = delz*vnnr * rsqinv;

        // tangential component

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // relative rotational velocity

	wr1 = (cri*omega[i][0] + crj*omega[j][0]) * rinv;
	wr2 = (cri*omega[i][1] + crj*omega[j][1]) * rinv;
	wr3 = (cri*omega[i][2] + crj*omega[j][2]) * rinv;

        // relative velocities

        vtr1 = vt1 - (delz*wr2-dely*wr3);
        vtr2 = vt2 - (delx*wr3-delz*wr1);
        vtr3 = vt3 - (dely*wr1-delx*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);

        // shear history effects

        touch[jj] = 1;
        shear = &allshear[numshearquants*jj]; // shear[] refers to a shear force here, not shear displacement

        shsqmag = shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2];

        // rotate shear forces onto new contact plane conserving length

        rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
        rsht *= rsqinv;
        if (shearupdate) {
          shear[0] -= rsht*delx;
          shear[1] -= rsht*dely;
          shear[2] -= rsht*delz;
          shsqnew = shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2];
          if (shsqnew!=0.0) {
              shratio=sqrt(shsqmag/shsqnew);
              shear[0] *= shratio; // conserve shear force length
              shear[1] *= shratio;
              shear[2] *= shratio;
          }
        }

        // then perform rotation for rigid-body SPIN
        omdel=(omega[i][0]+omega[j][0])*delx+(omega[i][1]+omega[j][1])*dely+(omega[i][2]+omega[j][2])*delz;
        wspinx=0.5*rsqinv*delx*omdel;
        wspiny=0.5*rsqinv*dely*omdel;
        wspinz=0.5*rsqinv*delz*omdel;
        shint0 = shear[0];
        shint1 = shear[1];
        shint2 = shear[2];

        if (shearupdate) {
            shear[0]=shint0+shint1*(-wspinz*dt)+shint2*wspiny*dt;
            shear[1]=shint0*wspinz*dt+shint1+shint2*(-wspinx*dt);
            shear[2]=shint0*(-wspiny*dt)+shint1*wspinx*dt+shint2;
        }

        // normal force = Hertzian contact

        ccel = kn*(radsum-r)*rinv;
        polyhertz = sqrt((radsum-r)*radi*radj / radsum);
        ccel *= polyhertz;

        // tangential forces done incrementally

        if (shearupdate) {
	  shear[0] -= polyhertz*kt*vtr1*dt;//shear displacement =vtr*dt
	  shear[1] -= polyhertz*kt*vtr2*dt;
	  shear[2] -= polyhertz*kt*vtr3*dt;
	}

        // rescale frictional forces if needed

        fs = sqrt(shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2]);
        fslim = xmu * fabs(ccel*r); //replaced fn with fslim

        if (fs > fslim) {
          if (fs != 0.0) {
           if (shearupdate) {
             shear[0] *= fslim/fs;
             shear[1] *= fslim/fs;
             shear[2] *= fslim/fs;
           }
          } else shear[0] = shear[1] = shear[2] = 0.0;
        }

        // forces & torques

        fx = delx*ccel + shear[0];
        fy = dely*ccel + shear[1];
        fz = delz*ccel + shear[2];
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        tor1 = rinv * (dely*shear[2] - delz*shear[1]);
        tor2 = rinv * (delz*shear[0] - delx*shear[2]);
        tor3 = rinv * (delx*shear[1] - dely*shear[0]);
	torque[i][0] -= cri*tor1;
	torque[i][1] -= cri*tor2;
	torque[i][2] -= cri*tor3;

        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
	  torque[j][0] -= crj*tor1;
	  torque[j][1] -= crj*tor2;
	  torque[j][2] -= crj*tor3;
        }

	double dur[3], dus[3], localdM[3], globaldM[3]; //~ Pass by reference
	double effectivekt = kt*polyhertz;
	double aveshearforce, slipdisp, oldshearforce, newshearforce;
	double incdissipf, nstr, sstr, incrementaldisp, rkt;
	
	double dspin_i[3],dspin_stm,spin_stm,dM_i[3],dM,K_spin,theta_r,M_limit,Dspin_energy,a,N;
	a = polyhertz;
	N = ccel*r;

	int consideronce = 0; //~ Consider contacts only once
	if (j < nlocal || tag[j] < tag[i]) consideronce = 1;
	  
	//~ Call function for rolling resistance model [KH - 25 October 2013]
	if (rolling && shearupdate && consideronce) {
	  /*~ The first '0' indicates that the rolling_resistance function is
	    called by the compute rather than the single function*/
	  rolling_resistance(0,i,j,numshearquants,delx,dely,delz,r,rinv,ccel,
			     fslim,effectivekt,torque,shear,dur,dus,localdM,
			     globaldM);
	}

	//~~ Call function for twisting resistance model [MO - 04 November 2014]]
	if (D_spin && shearupdate && consideronce) {
	  Deresiewicz1954_spin(0,i,j,numshearquants,r,torque,shear,dspin_i,
			       dspin_stm,spin_stm,dM_i,dM,K_spin,theta_r,
			       M_limit,Geq,Poiseq,Dspin_energy,a,N);
	}

	//~ Add contributions to traced energy [KH - 20 February 2014]
	if (pairenergy) {
	  //~ Ensure rkt cannot become infinite [KH - 21 October 2014]
	  if (effectivekt > 1.0e-30) rkt = 1.0/effectivekt;
	  else rkt = 0.0;
	  /*~~ A tiny effective kt should reult in a negligile energy.rkt=0 is reasonable 
	    under this circumstance to avoid an enormous energy value [MO - 22 October 2014]~~*/

	  /*~ Increment the friction energy only if the slip condition
	    is invoked*/
	  oldshearforce = sqrt(shsqmag);
	  if (fs > fslim && fslim > 0.0) {
	    //~ current shear displacement = fslim/effectivekt;
	    slipdisp = rkt*(fs-fslim);
	    aveshearforce = 0.5*(oldshearforce + fslim);

	    //~ slipdisp and aveshearforce are both positive
	    incdissipf = aveshearforce*slipdisp;
	    if (consideronce) dissipfriction += incdissipf;
	    if (trace_energy) shear[4] += incdissipf;
	  }

	  /*~ Update the normal contribution to strain energy which 
	    doesn't need to be calculated incrementally*/
	  nstr = 0.4*kn*polyhertz*deltan*deltan;
	  if (consideronce) normalstrain += nstr;
	  if (trace_energy) shear[5] = nstr;

	  //~~ Update the spin contribution [MO - 13 November 2014]
	  if (D_spin && shearupdate) {
	    if (consideronce) spinenergy += Dspin_energy;
	    if (trace_energy) shear[7] += Dspin_energy;
	  }

	  //~ The shear component does require incremental calculation
	  if (shearupdate) {
	    newshearforce = sqrt(shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2]);
	    incrementaldisp = rkt*(newshearforce - oldshearforce);
	    sstr = 0.5*incrementaldisp*(newshearforce + oldshearforce);
	    if (consideronce) shearstrain += sstr;
	    if (trace_energy) shear[6] += sstr;
	  }
	}

	//~ Assign current polyhertz value to shear[3] [KH - 23 November 2012]
	shear[3] = polyhertz;

        if (evflag) ev_tally_gran(i,j,nlocal,newton_pair,fx,fy,fz,x[i][0],x[i][1],x[i][2],
				  radius[i],x[j][0],x[j][1],x[j][2],radius[j]);
      }
    }
  }

  //~ Accumulate per-processor energy terms [KH - 17 October 2014]
  //~~ Added for D_spin model [MO - 13 November 2014]
  gatheredf = gatheredss = gatheredse = 0.0;
  if (pairenergy) {
    MPI_Allreduce(&dissipfriction,&gatheredf,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&shearstrain,&gatheredss,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&spinenergy,&gatheredse,1,MPI_DOUBLE,MPI_SUM,world);
  }

  if (vflag_fdotr) virial_fdotr_compute();
}



/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranShmHistory::settings(int narg, char **arg)
{

  if (narg != 3) error->all(FLERR,"Illegal pair_style command");
  if (domain->dimension == 2) error->all(FLERR,"PairGranShmHistory formulae only valid for 3D simulations");

  Geq = force->numeric(FLERR,arg[0]);//0.5 * (Gi+Gj)
  Poiseq = force->numeric(FLERR,arg[1]);//0.5 * (Poisi+Poisj)
  xmu = force->numeric(FLERR,arg[2]);

  if (Geq < 0.0 || Poiseq < 0.0 || xmu < 0.0 || Poiseq > 0.5) error->all(FLERR,"Illegal Shm pair parameter values");

  kn = 4.0*Geq / (3.0*(1.0-Poiseq)); // calculate these here to save effort in the compute
  kt = 4.0*Geq / (2.0-Poiseq);

  // convert Kn and Kt from pressure units to force/distance^2

  kn /= force->nktv2p;
  kt /= force->nktv2p;

}

/* ---------------------------------------------------------------------- */

double PairGranShmHistory::single(int i, int j, int itype, int jtype,
                                    double rsq,
                                    double factor_coul, double factor_lj,
                                    double &fforce)
{
  /*~ Since the shear forces are stored in shear, they do
    not need to be recalculated [KH - 13 December 2012]*/

  double radi,radj,radsum;
  double r,rinv;
  double ccel,polyhertz;
  double fs1,fs2,fs3,fs,fn;

  double *radius = atom->radius;
  radi = radius[i];
  radj = radius[j];
  radsum = radi + radj;

  double **x = atom->x;
  tagint *tag = atom->tag; //~ Write out the atom tags

  /*~ The number of optional entries in svector for energy
    tracing and/or rolling resistance [KH - 6 March 2014]*/
  // D_spin was added [MO - 05 November 2014]
  int optionalq = 4*trace_energy + 27*rolling + 27*D_spin;

  if (rsq >= radsum*radsum) {
    fforce = 0.0;
    svector[0] = svector[1] = svector[2] = svector[3] = 0.0;
    //~ The tags, radii etc. will not be zero [KH - 10 January 2013]
    svector[4] = tag[i];
    svector[5] = tag[j];
    for (int q = 0; q < 3; q++)
      svector[q+6] = x[i][q];
    svector[9] = radi;
    for (int q = 0; q < 3; q++)
      svector[q+10] = x[j][q];
    svector[13] = radj;

    /*~ Four optional entries for energy tracing [KH - 6 March 2014]
      Order: energy dissipated by friction, normal component of
      strain energy, shear component of strain energy, other
      contributions (e.g., from rolling resistance model)

      27 optional entries for the rolling resistance model 
      [KH - 29 July 2014]. Order: 
      dUr[*], accumulated dUr[*], dUs[*], accumulated dUs[*], 
      localdM[*], accumulated localdM[*], globaldM[*], accumulated 
      globaldM[*], ksbar, (K+1)_rolling, (K+1)_twisting*/
    if (optionalq > 0)
      for (int q = 0; q < optionalq; q++) svector[q+14] = 0.0;
    
    return 0.0;
  }

  r = sqrt(rsq);
  rinv = 1.0/r;

  // normal force = Hertzian contact

  ccel = kn*(radsum-r)*rinv;
  polyhertz = sqrt((radsum-r)*radi*radj / radsum);
  ccel *= polyhertz;

  // shear history effects
  // neighprev = index of found neigh on previous call
  // search entire jnum list of neighbors of I for neighbor J
  // start from neighprev, since will typically be next neighbor
  // reset neighprev to 0 as necessary

  int *jlist = list->firstneigh[i];
  int jnum = list->numneigh[i];
  double *allshear = fix_history->firstvalue[i];

  for (int jj = 0; jj < jnum; jj++) {
    neighprev++;
    if (neighprev >= jnum) neighprev = 0;
    if (jlist[neighprev] == j) break;
  }

  /*~ The number of shear quantities is 19 if rolling resistance
    is active [KH - 29 July 2014]*/
  /*~ Another 4 shear quantities were added for per-contact energy
    tracing [KH - 6 March 2014]. Modified by MO [04 November 2014]*/
  int numshearquants = 4 + 15*rolling + 20*D_spin + 4*trace_energy;
  double *shear = &allshear[numshearquants*neighprev];

  //~ Call function for rolling resistance model [KH - 30 October 2013]
  double dur[3], dus[3], localdM[3], globaldM[3]; //~ Pass by reference
  double **torque = atom->torque;
  double delx, dely, delz, fslim; //~ Calculate these to pass to function
  fslim = xmu * fabs(ccel*r);
  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];

  double effectivekt = kt*polyhertz;
  if (rolling) {
    /*~ The first '1' indicates that the rolling_resistance function is
      called by the single function rather than the compute*/
    rolling_resistance(1,i,j,numshearquants,delx,dely,delz,r,rinv,ccel,
		       fslim,effectivekt,torque,shear,dur,dus,localdM,globaldM);
  }

  double dspin_i[3],dspin_stm,spin_stm,dM_i[3],dM,K_spin,theta_r,M_limit,Dspin_energy,a,N;
   
  //~~ Added for Deresiewicz1954_spin function [MO - 4 November 2014]
  if (D_spin) {
    a = polyhertz;
    N = ccel*r;
    Deresiewicz1954_spin(1,i,j,numshearquants,r,torque,shear,dspin_i,
			 dspin_stm,spin_stm,dM_i,dM,K_spin,theta_r,M_limit,Geq,Poiseq,Dspin_energy,a,N);
  }

  /*~ Some of the following are included only for convenience as
    the data could instead be obtained from a dump of the sphere
    coordinates [KH - 13 December 2011]*/
  fforce = ccel;
  svector[0] = shear[0];
  svector[1] = shear[1];
  svector[2] = shear[2];
  svector[3] = ccel;
  svector[4] = tag[i];
  svector[5] = tag[j];
  for (int q = 0; q < 3; q++)
    svector[q+6] = x[i][q];
  svector[9] = radi;
  for (int q = 0; q < 3; q++)
    svector[q+10] = x[j][q];
  svector[13] = radj;

  //~ Add for energy tracing [KH - 6 March 2014]
  int nq = 4*trace_energy;
  if (trace_energy)
    for (int q = 0; q < 4; q++)
      svector[q+14] = shear[q+4]; //~ 4 rather than 3...

  //~ Add for the rolling resistance model [KH - 30 October 2013]
  if (rolling) {
    for (int q = 0; q < 3; q++) {
      svector[q+nq+14] = dur[q];
      svector[q+nq+17] = shear[q+nq+4]; //~ 4 rather than 3...
      svector[q+nq+20] = dus[q];
      svector[q+nq+23] = shear[q+nq+7];
      svector[q+nq+26] = localdM[q];
      svector[q+nq+29] = shear[q+nq+10];
      svector[q+nq+32] = globaldM[q];
      svector[q+nq+35] = shear[q+nq+13];
    }
    svector[nq+38] = shear[numshearquants-3]; //~ ksbar
    svector[nq+39] = shear[numshearquants-2]; //~ stored rolling (kappa+1)
    svector[nq+40] = shear[numshearquants-1]; //~ stored twisting (kappa+1)
  }
  
  //~ Add for the Deresiewicz1954_spin model [MO - 19 November 2014]
  if (D_spin) {
    for (int q = 0; q < 3; q++) {
      svector[q+nq+14] = dspin_i[q]; // dspin_i
      svector[q+nq+17] = shear[numshearquants-3+q]; // spin_i
      svector[q+nq+22] = dM_i[q];
      svector[q+nq+25] = shear[numshearquants-16+q];// M_i;
    }
    svector[nq+20] = dspin_stm; // system value 
    svector[nq+21] = spin_stm;  // system value
    svector[nq+28] = dM; // system value 
    svector[nq+29] = shear[numshearquants-4]; // M, system value
    svector[nq+30] = K_spin; 
    svector[nq+31] = theta_r; 
    svector[nq+32] = shear[numshearquants-12]; // step
    svector[nq+33] = M_limit; 
    svector[nq+34] = shear[numshearquants-5]; // M_star1
    svector[nq+35] = shear[numshearquants-6]; // M_star2
    svector[nq+36] = 0; 
    svector[nq+37] = 0; 
    svector[nq+38] = 0;
    svector[nq+39] = 0;
    svector[nq+40] = 0;
  } 

  return 0.0;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGranShmHistory::init_style()
{
  int i;

  // error and warning checks

  if (!atom->sphere_flag)
    error->all(FLERR,"Pair granular requires atom style sphere");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair granular requires ghost atoms store velocity");

  /*~ Have 19 shear quantities if rolling resistance is included
    [KH - 29 July 2014]*/
  /*~ Another 4 shear quantities are needed for per-contact energy
    tracing [KH - 6 March 2014]. Modified by MO [4 November 2014]*/
  int numshearquants = 4 + 15*rolling + 20*D_spin + 4*trace_energy;

  // need a granular neigh list and optionally a granular history neigh list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->size = 1;
  if (history) neighbor->requests[irequest]->history = 1;

  dt = update->dt;

  // if first init, create Fix needed for storing shear history

  if (history && fix_history == NULL) {
    char dnumstr[16];
    sprintf(dnumstr,"%d",numshearquants); //~ Now variable [KH - 23 May 2017]
    char **fixarg = new char*[4];
    fixarg[0] = (char *) "NEIGH_HISTORY";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "NEIGH_HISTORY";
    fixarg[3] = dnumstr;
    modify->add_fix(4,fixarg,1);
    delete [] fixarg;
    fix_history = (FixNeighHistory *) modify->fix[modify->nfix-1];
    fix_history->pair = this;
  }

  /*~ If rolling resistance is active, implicitly set up a fix,
    fix_old_omega, to store values of omega from the preceding
    timestep [KH - 24 October 2013]*/
  if (rolling) add_old_omega_fix();

  // check for FixFreeze and set freeze_group_bit

  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"freeze") == 0) break;
  if (i < modify->nfix) freeze_group_bit = modify->fix[i]->groupbit;
  else freeze_group_bit = 0;

  // check for FixRigid so can extract rigid body masses

  fix_rigid = NULL;
  for (i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) break;
  if (i < modify->nfix) fix_rigid = modify->fix[i];

  // check for FixPour and FixDeposit so can extract particle radii

  int ipour;
  for (ipour = 0; ipour < modify->nfix; ipour++)
    if (strcmp(modify->fix[ipour]->style,"pour") == 0) break;
  if (ipour == modify->nfix) ipour = -1;

  int idep;
  for (idep = 0; idep < modify->nfix; idep++)
    if (strcmp(modify->fix[idep]->style,"deposit") == 0) break;
  if (idep == modify->nfix) idep = -1;

  // set maxrad_dynamic and maxrad_frozen for each type
  // include future FixPour and FixDeposit particles as dynamic

  int itype;
  for (i = 1; i <= atom->ntypes; i++) {
    onerad_dynamic[i] = onerad_frozen[i] = 0.0;
    if (ipour >= 0) {
      itype = i;
      onerad_dynamic[i] = 
        *((double *) modify->fix[ipour]->extract("radius",itype));
    }
    if (idep >= 0) {
      itype = i;
      onerad_dynamic[i] = 
        *((double *) modify->fix[idep]->extract("radius",itype));
    }
  }

  double *radius = atom->radius;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & freeze_group_bit)
      onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]],radius[i]);
    else
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]],radius[i]);

  MPI_Allreduce(&onerad_dynamic[1],&maxrad_dynamic[1],atom->ntypes,
                MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&onerad_frozen[1],&maxrad_frozen[1],atom->ntypes,
                MPI_DOUBLE,MPI_MAX,world);
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGranShmHistory::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++)
      fwrite(&setflag[i][j],sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGranShmHistory::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGranShmHistory::write_restart_settings(FILE *fp)
{
  fwrite(&kn,sizeof(double),1,fp);
  fwrite(&kt,sizeof(double),1,fp);
  fwrite(&Geq,sizeof(double),1,fp);
  fwrite(&Poiseq,sizeof(double),1,fp);
  fwrite(&xmu,sizeof(double),1,fp);

  //~ Added energy terms [KH - 28 February 2014]
  //~~ Added for D_spin model [MO - 13 November 2014]
  fwrite(&gatheredf,sizeof(double),1,fp);
  fwrite(&gatheredss,sizeof(double),1,fp);
  fwrite(&gatheredse,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGranShmHistory::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&kn,sizeof(double),1,fp);
    fread(&kt,sizeof(double),1,fp);
    fread(&Geq,sizeof(double),1,fp);
    fread(&Poiseq,sizeof(double),1,fp);
    fread(&xmu,sizeof(double),1,fp);

    /*~ Added energy terms. The total energy is read to the root
    proc and is NOT broadcast to all procs as only the total summed
    across all procs is of interest [KH - 28 February 2014]*/
    //~~ Added for D_spin model [MO - 13 November 2014]
    fread(&dissipfriction,sizeof(double),1,fp);
    fread(&shearstrain,sizeof(double),1,fp);
    fread(&spinenergy,sizeof(double),1,fp);
  }
  MPI_Bcast(&kn,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&kt,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&Geq,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&Poiseq,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&xmu,1,MPI_DOUBLE,0,world);
}
