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
   Contributing author: Kevin Hanley (Imperial)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_crushing.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "random_park.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "integrate.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "fix_multistress.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 1.0e-20
#define BIG   1.0e20

/* ----------------------------------------------------------------------
The syntax of the input command is as follows:

fix ID group crushing outputflag seed m sigma0 d0 chi redtype {reduction} constante {commlimit} {m2} {sigma02} {d02}

   outputflag Set this at 1 to display detailed crushing output; 0 to suppress this output
   seed       Seed of the random number generator
   m          Weibull modulus
   sigma0     Characteristic strength at which 37% of particles of diameter d0 survive
   d0         Nominal particle diameter in m
   chi        Parameter in brittle failure criterion of Christensen
   redtype    Flag = 1 if diameter reduced to just give touching contact; else 0
   reduction  If redtype = 0, the fractional reduction in radius between 0 and 1
   constante  Flag = 1 if void ratio is held constant; else 0
   commlimit  Optional: a comminution limit on radii in m
   m2         Optional: Weibull modulus for > 1 breakage
   sigma02    Optional: sigma0 for particles of diam. d02 after > 1 breakage
   d02        Optional: Nominal particle diameter in m for > 1 breakage
 ---------------------------------------------------------------------- */

FixCrushing::FixCrushing(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  restart_peratom = 1; //~ Global information is saved to the restart file
  nevery = 1; //~ Set how often to run the end_of_step function
  peratom_flag = 1;
  size_peratom_cols = 3; //~ m, sigma0, d0
  peratom_freq = 1;
  create_attribute = 1;

  if (narg < 11 || narg > 16)
    error->all(FLERR,"Illegal fix crushing command");

  if (strcmp(atom->atom_style,"sphere") != 0)
    error->all(FLERR,"Fix crushing is defined only for the sphere atom style");
  
  if (domain->dimension == 2)
    error->all(FLERR,"The 2D case has not been written for fix_crushing");

  //~ Issue a warning if SI units are not used
  if (strcmp(update->unit_style,"si") != 0)
    error->warning(FLERR,"SI units are best when using fix crushing");

  //~ Read in all user-specified data
  displaymessages = atoi(arg[3]); //~ Write information to the screen (1) or not (0)
  seed = atoi(arg[4]);
  weibullparams[0][0] = atof(arg[5]); //~ m for first breakage
  weibullparams[1][0] = atof(arg[6]); //~ sigma0 for first breakage
  weibullparams[2][0] = atof(arg[7]); //~ d0 for first breakage
  chiplusone = atof(arg[8])+1.0;
  redtype = atoi(arg[9]);

  int iarg = 10;
  if (redtype == 0) {
    reduction = atof(arg[iarg]);

    if (reduction <= 0 || reduction >= 1)
      error->all(FLERR,"Radius reduction in fix crushing must be between 0 and 1");

    iarg++;
  } else if (redtype != 1) error->all(FLERR,"Illegal redtype specification in fix crushing command");

  constante = atoi(arg[iarg]);

  if (narg == iarg+1 || narg == iarg+4) commlimit = -1.0;
  else if (narg == iarg+2 || narg == iarg+5) {
    commlimit = atof(arg[iarg+1]);
    iarg++;
  } else error->all(FLERR,"Illegal number of arguments in fix crushing command");
  
  if (narg == iarg+4) {
    weibullparams[0][1] = atof(arg[iarg+1]); //~ m for second and subsequent breakages
    weibullparams[1][1] = atof(arg[iarg+2]); //~ sigma0 for second and subsequent breakages
    weibullparams[2][1] = atof(arg[iarg+3]); //~ d0 for second and subsequent breakages
  } else { //~ Use the same values for the first, second... breakages
    weibullparams[0][1] = weibullparams[0][0];
    weibullparams[1][1] = weibullparams[1][0];
    weibullparams[2][1] = weibullparams[2][0];
  }

  //~ Now check that these inputs are sensible
  if (seed <= 0 || weibullparams[0][0] <= 0.0 || weibullparams[1][0] <= 0.0 || 
      weibullparams[2][0] <= 0.0 || chiplusone <= 1.0 || weibullparams[0][1] <= 0.0 ||
      weibullparams[1][1] <= 0.0 || weibullparams[2][1] <= 0.0)
    error->all(FLERR,"Fix crushing parameter must be positive");

  if (displaymessages != 0 && displaymessages != 1)
    error->all(FLERR,"outputflag in fix crushing must be 0 or 1");

  if (constante != 0 && constante != 1)
    error->all(FLERR,"constante flag in fix crushing must be 0 or 1");

  //~ Confirm that a granular pairstyle is used
  if (!force->pair_match("gran",0))
    error->all(FLERR,"A granular pairstyle must be used for fix crushing");

  random = new RanPark(lmp,seed); //~ Establish the random number generator

  // perform initial allocation of atom-based array
  // register with Atom class
  cparams = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  /*~ The cparams array stores three pieces of data for each
    atom: 1) the uniaxial compressive strength; 2) the uniaxial 
    tensile strength; 3) the number of particle breakage events.
    Initialise everything at zero regardless of masking.*/
  for (int i = 0; i < atom->nlocal; i++)
    cparams[i][0] = cparams[i][1] = cparams[i][2] = 0.0;

  PI = 4.0*atan(1.0); //~ Calculate pi once as it is used often
}

/* ---------------------------------------------------------------------- */

FixCrushing::~FixCrushing()
{
  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored array
  memory->destroy(cparams);

  //~ Delete the stress compute if fix_multistress is not active
  mstressid = -1;
  for (int q = 0; q < modify->nfix; q++)
    if (strcmp(modify->fix[q]->style,"multistress") == 0) mstressid = q;

  if (mstressid < 0) {
    modify->delete_compute(id_stress);
    delete [] id_stress;
  }

  delete random;
}

/* ---------------------------------------------------------------------- */

int FixCrushing::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCrushing::init()
{
  /*~ Check that all relevant (masked) particles have non-zero radii
    and if a comminution limit is specified, ensure that the initial
    radii are greater than or equal to this limit*/
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double minrad = BIG; //~ minrad contains the radius of the smallest atom

  for (int i = 0; i < nlocal; i++)
    if ((mask[i] & groupbit) && radius[i] < minrad)
      minrad = radius[i];

  if (minrad <= SMALL)
    error->all(FLERR,"Fix crushing requires extended particles");
  else if (commlimit > 0.0 && minrad < commlimit)
    error->all(FLERR,"Particles smaller than the comminution limit are in the system");

  /*~ Allocate initial strengths to all particles. This is run only
    if the data are not read in from a restart file, noting that the 
    compressive strengths are negative, if present.*/

  if (force->pair_match("gran/hertz",0) || force->pair_match("gran/shm",0)) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	if (cparams[i][0] > -1.0*SMALL) {
	  cparams[i][0] = strength_calculation(i,radius[i]);
	  cparams[i][1] = -1.0*cparams[i][0]/chiplusone;
	}
  } else
    error->all(FLERR,"Currently fix crushing requires a Hertzian pairstyle");
  
  /*~ Calculate the total volume of particles and the domain initially if the void
    ratio is to be held constant*/
  if (constante == 1) {
    double pvolume = 0.0; //~ Volume of particles on a processor
    double totalvolume = 1.0;
    double fourpioverthree = 4.0*PI/3.0;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	pvolume += fourpioverthree*radius[i]*radius[i]*radius[i];
  
    //~ Gather pvolume from all processors in totalpvolume
    MPI_Allreduce(&pvolume,&totalpvolume,1,MPI_DOUBLE,MPI_SUM,world);

    //~ Find the current size of the simulation domain
    for (int i = 0; i < domain->dimension; i++)
      totalvolume *= (domain->boxhi[i]-domain->boxlo[i]);
    
    voidratio = (totalvolume - totalpvolume)/totalpvolume;
  }

  //~ If redtype == 1, need to set up a mechanism for neighbour list updating
  if (redtype == 1) {
    // need a full neighbor list, built whenever re-neighboring occurs
    int irequest = neighbor->request((void *) this);
    neighbor->requests[irequest]->pair = 0;
    neighbor->requests[irequest]->fix = 1;
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }
}

/* ---------------------------------------------------------------------- */

void FixCrushing::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixCrushing::set_weibull_parameters(int numbreaks)
{
  m = weibullparams[0][numbreaks];
  sigma0 = weibullparams[0][numbreaks];
  d0 = weibullparams[2][numbreaks];
}

/* ---------------------------------------------------------------------- */

double FixCrushing::strength_calculation(int i, double radius)
{
  /*~ Assign a suitable uniaxial tensile and compressive strength to 
    particle i*/

  //~ Import data from pair and calculate G and Poisson's ratio
  double kn,kt,xmu,smod,pr;
  kn = force->pair->kn;
  kt = force->pair->kt;
  xmu = force->pair->xmu;
  smod = 0.75*kn*kt/(3.0*kn-kt);
  pr = (3.0*kn-2.0*kt)/(3.0*kn-kt);

  //~ Pick appropriate values for m, sigma0 and d0 from weibullparams
  if (static_cast<int> (cparams[i][2]) == 0) set_weibull_parameters(0);
  else set_weibull_parameters(1);

  double firstfixedterm = pow(sigma0,m);
  double secondfixedterm = 3.0*(3.0/32.0 + sqrt(2.0)/24.0 + xmu*(sqrt(2.0)/12.0 - 0.25) + xmu*xmu*(0.5 - sqrt(2.0)/3.0))/((2.0-sqrt(2.0))*(1.0+xmu)); //~ Russell & Muir-Wood, 2009, Eq. 19
  double sigmaf = fabs(pow(-1.0*firstfixedterm*log(random->uniform())*pow(0.5*d0/radius,3.0),(1.0/m)));
  double forcef = 4.0*sigmaf*radius*radius;

  double estar,contactarea,compressivestrength;

  /*~ Assume contact with a sphere of infinite radius and stiffness. Recall
    that Young's modulus is 2*smod*(1+pr).*/
  estar = 2.0*smod*(1+pr)/(1-pr*pr);
  contactarea = PI*pow(0.75*forcef*radius/estar,(2.0/3.0));
  
  //~ The -abs operation below is used to yield negative compressive strengths
  compressivestrength = -1.0*fabs(secondfixedterm*forcef/contactarea); //~ Russell & Muir-Wood, 2009, Eq. 19

  return compressivestrength;
}

/* ---------------------------------------------------------------------- */

void FixCrushing::setup(int vflag)
{
  //~ Set up or find an existing stress compute
  
  //~ Determine whether fix_multistress is active
  mstressid = -1;

  for (int q = 0; q < modify->nfix; q++)
    if (strcmp(modify->fix[q]->style,"multistress") == 0) {
      mstressid = q;
 
      if (constante == 1) {
	//~ Set constanteflag in fix multistress to 1 if lvstressflag == 1
	int lvstressflag = ((FixMultistress *) modify->fix[mstressid])->lvstressflag;
	
	if (lvstressflag == 0)
	  error->all(FLERR,"Constant e option of fix crushing requires the linkvolstress option of fix multistress to be active");
      }
      break;
    }

  /*~ If fix multistress is active, find the pointer to the 
    existing stress compute; otherwise set up a new stress compute*/
  if (mstressid >= 0) {
    id_stress = new char[20];
    strcpy(id_stress,modify->fix[mstressid]->id);
  } else {
    int n = strlen(id) + strlen("_stress") + 1;
    id_stress = new char[n];
    strcpy(id_stress,id);
  }

  strcat(id_stress,"_stress");

  //~ Check whether the compute already exists
  int icompute = modify->find_compute(id_stress);

  if (mstressid < 0 && icompute < 0) {
    char **snewarg = new char*[6];
    snewarg[0] = id_stress;
    snewarg[1] = (char *) "all";
    snewarg[2] = (char *) "stress/atom";
    snewarg[3] = (char *) "pair";
    snewarg[4] = (char *) "fix";
    snewarg[5] = (char *) "bond";

    modify->add_compute(6,snewarg);
    delete [] snewarg;
  }

  //~ Confirm that the compute exists now
  icompute = modify->find_compute(id_stress);
  if (icompute < 0) error->all(FLERR,"Stress ID does not exist in fix crushing");
  tstress = modify->compute[icompute];
}

/* ---------------------------------------------------------------------- */

void FixCrushing::pre_force(int vflag)
{
  if (update->integrate->vflag <= 2) update->integrate->vflag += 4;
}

/* ---------------------------------------------------------------------- */

void FixCrushing::end_of_step()
{
  if (tstress->invoked_peratom != update->ntimestep)
    tstress->compute_peratom(); //~ Update the atom stresses if necessary
  double **stresses = tstress->array_atom; //~ Fetch the stresses

  /*~ Note that tensile stresses are positive in Christensen (2000) (or
    Russell and Muir-Wood (2009)), but compressive stresses are positive
    in the stresses array. Swap the signs around in the calculations
    below to take this into account. All of the entries in the cstrength
    vector are negative (hence the need for the fabs operation below).*/

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double i1,i2,j2,kappa,lhs,rhs;
  double oneover3 = 1.0/3.0; //~ To minimise the number of divisions
  double oneoverroot3 = sqrt(oneover3);
  double chiplusonesqd = chiplusone*chiplusone;
  double oneoverchiplusone = 1.0/chiplusone;
  double tempredvolume = 0.0; //~ Reinitialise on each timestep

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      //~ Calculate the required invariants for particle i
      i1 = -1.0*(stresses[i][0] + stresses[i][1] + stresses[i][2]);
      i2 = stresses[i][0]*stresses[i][1] + stresses[i][0]*stresses[i][2]
	+ stresses[i][1]*stresses[i][2] - stresses[i][3]*stresses[i][3]
	- stresses[i][4]*stresses[i][4] - stresses[i][5]*stresses[i][5];
      j2 = i1*i1*oneover3 - i2;

      //~ Calculate kappa
      kappa = chiplusone*fabs(cparams[i][0])*oneoverroot3;

      //~ Now calculate the left and right sides of the failure condition
      lhs = (chiplusone-1.0)*kappa*i1*oneoverroot3 + chiplusonesqd*j2;
      rhs = kappa*kappa*oneoverchiplusone;

      //~ If failure has occurred, reduce the particle diameter
      if (lhs >= rhs) tempredvolume += failure_occurs(i);
    }
  }

  //~ Accumulate the volume reduction of all particles in redvolume
  double redvolume = 0.0;
  MPI_Allreduce(&tempredvolume,&redvolume,1,MPI_DOUBLE,MPI_SUM,world);

  /*~ Shrink the domain so that the void ratio is held constant if
    constante == 1*/
  if (redvolume > SMALL && constante == 1) {
    totalpvolume -= redvolume;
    domain->initialvolume = (voidratio+1)*totalpvolume;
  }
}

/* ---------------------------------------------------------------------- */

double FixCrushing::failure_occurs(int i)
{
  //~ Check whether the comminution limit has been reached, if present
  if (atom->radius[i] > commlimit) {
    //~ Increment the number of particle failures
    cparams[i][2]++;

    //~ Displaying an optional user message
    if (displaymessages)
      fprintf(screen,"\nFailure number %u of particle %u on timestep %u.\n",static_cast<int> (cparams[i][2]),atom->tag[i],update->ntimestep);

    //~ Reduce the radius of the particle
    double volred = reduce_radius(i);
    return volred;
  } else return 0.0;
}

/* ---------------------------------------------------------------------- */

double FixCrushing::reduce_radius(int i)
{
  double *radius = atom->radius;
  double oldradius = radius[i];
  double mindist = -1.0; //~ The amount by which to reduce the radius

  if (redtype == 1) {
    //~ Reduce the diameter to just lose contact with neighbouring particles 
    double **x = atom->x;
    int *mask = atom->mask;

    //~ Check for contacts using the neighbor list
    int i2,ii,j,jj,inum,jnum;
    int *ilist,*jlist,*numneigh,**firstneigh;
    double delx,dely,delz,r,rsq,radsum,overlap;

    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    for (ii = 0; ii < inum; ii++) {
      i2 = ilist[ii];
      jlist = firstneigh[i2];
      jnum = numneigh[i2];

      //~ Check whether i and i2 are the same particle
      //~ If so, loop through its neighbours and break
      if (radius[i] == radius[i2] && x[i][0] == x[i2][0]
	  && x[i][1] == x[i2][1] && x[i][2] == x[i2][2]) {
	for (jj = 0; jj < jnum; jj++) {
	  j = jlist[jj];
	  j &= NEIGHMASK;

	  delx = x[i][0] - x[j][0];
	  dely = x[i][1] - x[j][1];
	  delz = x[i][2] - x[j][2];
	  rsq = delx*delx + dely*dely + delz*delz;
	  r = sqrt(rsq);
	  radsum = radius[i] + radius[j];
	
	  if (rsq < radsum*radsum) {
	    overlap = radsum - r;
	    if (overlap > mindist) mindist = overlap;
	  }
	}
	/*~ Include a small tolerance to ensure that contact is lost
	  and the shear displacements are reset to zero*/
	mindist *= 1.0001;

	break;
      }
    }
  } else {
    //~ Reduce the radius by the fixed fraction 'reduction'
    mindist = reduction*radius[i];
  }
  
  /*~ Reduce the radius of particle i by mindist if the comminution
    limit has not been reached.*/
  if (radius[i]-mindist >= commlimit) {//~ Comminution limit inactive or not reached
    if (radius[i]-mindist <= 0.0)
      error->all(FLERR,"Check the simulation conditions: particles are inside of other particles");
    
    radius[i] -= mindist;
  } else {
    radius[i] = commlimit; //~ Set equal to the comminution limit
  
    if (displaymessages)
      fprintf(screen,"\tComminution limit reached for particle %u.\n",atom->tag[i]);
  }

  //~ Display optional messages for the user
  if (displaymessages)
    fprintf(screen,"\tRadius changed from %1.3e to %1.3e.\n",oldradius,radius[i]);

  /*~ Change the stored values of uniaxial tensile and compressive
    strength using the radius which has already been reduced*/
  change_strengths(i,radius[i]);

  //~ Return the change in solid volume
  double volchange = (4.0*PI/3.0)*(oldradius*oldradius*oldradius-radius[i]*radius[i]*radius[i]);

  return volchange;
}

/* ---------------------------------------------------------------------- */

void FixCrushing::change_strengths(int i, double radius)
{
  //~ Display first part of optional message for the user
  if (displaymessages)
    fprintf(screen,"\tStrengths reduced from \t%1.3e (compressive) and %1.3e (tensile)\n",cparams[i][0],cparams[i][1]);

  int counter = 0;
  int counterlimit = 100;
  double newcstrength = 0.0;

  //~ Force the new compressive strength to be larger in magnitude than the old value
  while (-1.0*newcstrength <= -1.0*cparams[i][0]) {//~ cparams[i][0] is negative
    newcstrength = strength_calculation(i,radius);
  
    counter++;
    if (counter > counterlimit) {
      newcstrength = cparams[i][0]; //~ Retain the old value
      break; //~ To prevent an infinite loop
    } 
  }

  cparams[i][0] = newcstrength;
  cparams[i][1] = -1.0*cparams[i][0]/chiplusone;

  //~ Complete the optional user message
  if (displaymessages)
    fprintf(screen,"\t\t\tto\t%1.3e (compressive) and %1.3e (tensile).\n",cparams[i][0],cparams[i][1]);

  if (counter > counterlimit)
    error->message(FLERR,"Loop exit condition invoked in FixCrushing to prevent an infinite loop");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixCrushing::memory_usage()
{
  double bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixCrushing::grow_arrays(int nmax)
{
  memory->grow(cparams,nmax,3,"fix_crushing:cparams");
  array_atom = cparams;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixCrushing::copy_arrays(int i, int j, int delflag)
{
  cparams[j][0] = cparams[i][0];
  cparams[j][1] = cparams[i][1];
  cparams[j][2] = cparams[i][2];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixCrushing::set_arrays(int i)
{
  cparams[i][0] = strength_calculation(i,atom->radius[i]);
  cparams[i][1] = -1.0*cparams[i][0]/chiplusone;
  cparams[i][2] = 0.0; //~ The number of breakage events
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixCrushing::pack_exchange(int i, double *buf)
{
  buf[0] = cparams[i][0];
  buf[1] = cparams[i][1];
  buf[2] = cparams[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixCrushing::unpack_exchange(int nlocal, double *buf)
{
  cparams[nlocal][0] = buf[0];
  cparams[nlocal][1] = buf[1];
  cparams[nlocal][2] = buf[2];
  return 3;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixCrushing::pack_restart(int i, double *buf)
{
  buf[0] = 4;
  buf[1] = cparams[i][0];
  buf[2] = cparams[i][1];
  buf[3] = cparams[i][2];
  return 4;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixCrushing::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  cparams[nlocal][0] = extra[nlocal][m++];
  cparams[nlocal][1] = extra[nlocal][m++];
  cparams[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixCrushing::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixCrushing::size_restart(int nlocal)
{
  return 4;
}
