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

#include "string.h"
#include "stdio.h"
#include "fix_shear_history.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "pair.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "stdlib.h"  //added in GM

using namespace LAMMPS_NS;
using namespace FixConst;

#define MAXTOUCH 15

/* ---------------------------------------------------------------------- */

FixShearHistory::FixShearHistory(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  restart_peratom = 1;
  create_attribute = 1;

  //~ Read in the number of shear quantities if available [KH - 21 November 2012]
  if (narg == 4) num_quants = force->inumeric(FLERR,arg[3]); // added in GM
  else num_quants = 3; //~ Assume the default value of 3 instead

  // perform initial allocation of atom-based arrays
  // register with atom class

  npartner = NULL;
  partner = NULL;
  shearpartner = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  // initialize npartner to 0 so neighbor list creation is OK the 1st time

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) npartner[i] = 0;
}

/* ---------------------------------------------------------------------- */

FixShearHistory::~FixShearHistory()
{
  // unregister this fix so atom class doesn't invoke it any more

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored arrays

  memory->destroy(npartner);
  memory->destroy(partner);
  memory->destroy(shearpartner);
}

/* ---------------------------------------------------------------------- */

int FixShearHistory::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= MIN_PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixShearHistory::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,
               "Pair style granular with history requires atoms have IDs");

  int dim;
  computeflag = (int *) pair->extract("computeflag",dim);
}

/* ----------------------------------------------------------------------
   called by setup of run or minimize
   called by write_restart as input script command
   only invoke pre_exchange() if neigh list stores more current history info
     than npartner/partner arrays in this fix
   that will only be case if pair->compute() has been invoked since
     update of npartner/npartner
   this logic avoids 2 problems:
     run 100; write_restart; run 100
       setup_pre_exchange is called twice (by write_restart and 2nd run setup)
       w/out a neighbor list being created in between
     read_restart; run 100
       setup_pre_exchange called by run setup whacks restart shear history info
------------------------------------------------------------------------- */

void FixShearHistory::setup_pre_exchange()
{
  if (*computeflag) pre_exchange();
  *computeflag = 0;
}

/* ----------------------------------------------------------------------
   copy shear partner info from neighbor lists to atom arrays
   so can be migrated or stored with atoms
------------------------------------------------------------------------- */

void FixShearHistory::pre_exchange()
{
  int i,j,ii,jj,m,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  // zero npartner for all current atoms

  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) npartner[i] = 0;

  // copy shear info from neighbor list atoms to atom arrays

  int *tag = atom->tag;
  NeighList *list = pair->list;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = list->listgranhistory->firstneigh;
  firstshear = list->listgranhistory->firstdouble;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    allshear = firstshear[i];
    jnum = numneigh[i];
    touch = firsttouch[i];

    for (jj = 0; jj < jnum; jj++) {
      if (touch[jj]) {
        shear = &allshear[num_quants*jj]; // was 3, modified GM
        j = jlist[jj];
        j &= NEIGHMASK;
        if (npartner[i] < MAXTOUCH) {
          m = npartner[i];
          partner[i][m] = tag[j];
	  for (int kk = 0; kk < num_quants; kk++)
            shearpartner[i][m][kk] = shear[kk]; // put the loop in, modified GM
        }
        npartner[i]++;
        if (j < nlocal) {
          if (npartner[j] < MAXTOUCH) {
            m = npartner[j];
            partner[j][m] = tag[i];
	    for (int kk = 0; kk < num_quants; kk++) {
	      /*~ Modified this when the rolling resistance shear
		parameters were stored in the shear array [KH - 5
		November 2013]*/
	      shearpartner[j][m][kk] = -shear[kk]; // put the loop in, modified GM
	    }
          }
          npartner[j]++;
        }
      }
    }
  }

  // test for too many touching neighbors

  int flag = 0;
  for (i = 0; i < nlocal; i++)
    if (npartner[i] >= MAXTOUCH) flag = 1;
  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all)
    error->all(FLERR,"Too many touching neighbors - boost MAXTOUCH");
}

/* ---------------------------------------------------------------------- */

void FixShearHistory::min_setup_pre_exchange()
{
  if (*computeflag) pre_exchange();
  *computeflag = 0;
}

/* ---------------------------------------------------------------------- */

void FixShearHistory::min_pre_exchange()
{
  pre_exchange();
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixShearHistory::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += nmax*MAXTOUCH * sizeof(int);
  bytes += nmax*MAXTOUCH*num_quants * sizeof(double); // was 3, modified GM
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixShearHistory::grow_arrays(int nmax)
{
  memory->grow(npartner,nmax,"shear_history:npartner");
  memory->grow(partner,nmax,MAXTOUCH,"shear_history:partner");
  memory->grow(shearpartner,nmax,MAXTOUCH,num_quants,"shear_history:shearpartner"); // changed from 3 to num_quants, the number of per-contact quantities required, modified GM
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixShearHistory::copy_arrays(int i, int j, int delflag)
{
  npartner[j] = npartner[i];
  for (int m = 0; m < npartner[j]; m++) {
    partner[j][m] = partner[i][m];
    for (int k = 0; k < num_quants; k++) // put the loop in, modified GM
      shearpartner[j][m][k] = shearpartner[i][m][k];
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixShearHistory::set_arrays(int i)
{
  npartner[i] = 0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixShearHistory::pack_exchange(int i, double *buf)
{
  int m = 0;
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
    for (int k = 0; k < num_quants; k++) // put the loop in, modified GM
      buf[m++] = shearpartner[i][n][k];
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixShearHistory::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  npartner[nlocal] = static_cast<int> (buf[m++]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<int> (buf[m++]);
    for (int k = 0; k < num_quants; k++) // put the loop in, modified GM
      shearpartner[nlocal][n][k] = buf[m++];
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixShearHistory::pack_restart(int i, double *buf)
{
  int m = 0;
  buf[m++] = (num_quants+1)*npartner[i] + 2; // changed from 4 to num_quants+1 , modified GM
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
    for (int k = 0; k < num_quants; k++) // put the loop in, modified GM
      buf[m++] = shearpartner[i][n][k];
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixShearHistory::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  npartner[nlocal] = static_cast<int> (extra[nlocal][m++]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<int> (extra[nlocal][m++]);
    for (int k = 0; k < num_quants; k++) // put the loop in, modified GM
      shearpartner[nlocal][n][k] = extra[nlocal][m++];
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixShearHistory::maxsize_restart()
{
  return (num_quants+1)*MAXTOUCH + 2; // changed from 4 to num_quants+1 , modified GM
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixShearHistory::size_restart(int nlocal)
{
  return (num_quants+1)*npartner[nlocal] + 2; // changed from 4 to num_quants+1 , modified GM
}
