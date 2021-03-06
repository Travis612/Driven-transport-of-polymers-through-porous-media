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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "fix_test_bc2.h"
#include "math_const.h"
#include "atom.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "input.h"
#include "force.h"
#include "comm.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace std;
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixTestBc2::FixTestBc2(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{


  if (narg != 4) error->all(FLERR,"Illegal fix test/bc2 command");


  rhow = force->numeric(FLERR,arg[3]);


  time_integrate = 1;
  maxatom=1;
  comm_reverse = 4;
  memory->create(nw,maxatom,3,"test/bc2:nw");
}

/* ---------------------------------------------------------------------- */

FixTestBc2::~FixTestBc2()
{
  memory->destroy(nw);
}



/* ---------------------------------------------------------------------- */

int FixTestBc2::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;

  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTestBc2::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;



  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->half = 1;
  neighbor->requests[irequest]->full = 0;


}

/* ---------------------------------------------------------------------- */

void FixTestBc2::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixTestBc2::initial_integrate(int vflag)
{
  double dtfm;

  // update v and x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;


    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];

      }

  // calculate boundary volume fraction

  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,r;
  int *ilist,*jlist,*numneigh,**firstneigh;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  Pair *pair = force->pair_match("dpd/bc2",1);
  if (pair == NULL) error->all(FLERR,"No pair dpd/bc2 for fix test/bc2");
  int jtmp;
  double *gcutoff = (double *) pair->extract("cut_global",jtmp);

  double cutg = *gcutoff;
  double cutg2 = pow(cutg,2);
  double absnw,dotp,enx,eny,enz;
  int n = atom->ntypes;
  int *ntw = (int *) pair->extract("nwt",jtmp);
  int nwt = *ntw;
  int nft = n - nwt;
  int newton_pair = force->newton_pair;


  // reallocate nw array if necessary

  if (atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(nw);
    memory->create(nw,maxatom,3,"test/bc2:nw");
  }



  int tmp;

  int index1 = atom->find_custom("phi", tmp);

  double *phi = atom->dvector[index1];


  for (i = 0; i < nall; i++){
    phi[i] = 0.0;
    nw[i][0] = 0.0;
    nw[i][1] = 0.0;
    nw[i][2] = 0.0;

  }

  for (int ii=0; ii < inum; ii++) {

    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutg2 && jtype > nft ) {
        r = sqrt(rsq);
        phi[i] += 105*(1+3*(r/cutg))*pow((1-r/cutg),3)/(16*MY_PI*rhow*pow(cutg,3));
        nw[i][0] += 315*delx*pow((1-(r/cutg)),2)/(4*MY_PI*rhow*pow(cutg,5));
        nw[i][1] += 315*dely*pow((1-(r/cutg)),2)/(4*MY_PI*rhow*pow(cutg,5));
        nw[i][2] += 315*delz*pow((1-(r/cutg)),2)/(4*MY_PI*rhow*pow(cutg,5));
      } else if (rsq < cutg2 && itype > nft  ) {
        r = sqrt(rsq);
        phi[j] += 105*(1+3*(r/cutg))*pow((1-r/cutg),3)/(16*MY_PI*rhow*pow(cutg,3));
        nw[j][0] -= 315*delx*pow((1-(r/cutg)),2)/(4*MY_PI*rhow*pow(cutg,5));
        nw[j][1] -= 315*dely*pow((1-(r/cutg)),2)/(4*MY_PI*rhow*pow(cutg,5));
        nw[j][2] -= 315*delz*pow((1-(r/cutg)),2)/(4*MY_PI*rhow*pow(cutg,5));
      }
    }
  }

   if (newton_pair) comm->reverse_comm_fix(this);

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (phi[i] > 0.5) {
          dtfm = dtf / mass[type[i]];
          x[i][0] -= dtv * v[i][0];
          x[i][1] -= dtv * v[i][1];
          x[i][2] -= dtv * v[i][2];
          absnw = sqrt(nw[i][0]*nw[i][0] + nw[i][1]*nw[i][1] + nw[i][2]*nw[i][2]);
          enx=nw[i][0]/absnw;
          eny=nw[i][1]/absnw;
          enz=nw[i][2]/absnw;
          dotp = v[i][0]*enx + v[i][1]*eny + v[i][2]*enz;
          v[i][0] = -v[i][0] + 2.0*max(0.0,dotp)*enx;
          v[i][1] = -v[i][1] + 2.0*max(0.0,dotp)*eny;
          v[i][2] = -v[i][2] + 2.0*max(0.0,dotp)*enz;
          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];

        }
      }


}

/* ---------------------------------------------------------------------- */

void FixTestBc2::final_integrate()
{
  double dtfm;

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;


    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }

}

/* ---------------------------------------------------------------------- */

int FixTestBc2::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  int tmp;
  int index1 = atom->find_custom("phi", tmp);
  double *phi = atom->dvector[index1];


  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = phi[i];
    buf[m++] = nw[i][0];
    buf[m++] = nw[i][1];
    buf[m++] = nw[i][2];

  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixTestBc2::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  int tmp;
  int index1 = atom->find_custom("phi", tmp);
  double *phi = atom->dvector[index1];


  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    phi[j] += buf[m++];
    nw[j][0] += buf[m++];
    nw[j][1] += buf[m++];
    nw[j][2] += buf[m++];

  }
}



/* ---------------------------------------------------------------------- */

void FixTestBc2::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixTestBc2::memory_usage()
{
  double bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}
