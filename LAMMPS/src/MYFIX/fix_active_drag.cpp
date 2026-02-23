#include "atom.h"
#include "atom_masks.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "random_mars.h"
#include "update.h"

#include "domain.h"
#include "force.h"
#include "lammps.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pointers.h"
#include "update.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>

#include "fix_active_drag.h"

using namespace std;
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace NeighConst;

static constexpr double SMALL = 0.001;

FixActiveDrag::FixActiveDrag(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  jgroup = group->find(arg[3]);
  jgroupbit = group->bitmask[jgroup];
  if (jgroup < 0) error->all(FLERR, "Invalid group ID");
  // prefactor of tengential active force (ct*vrelt)
  ct = utils::numeric(FLERR, arg[4], false, lmp);
  // prefactor of normal active force (cn*vreln)
  cn = utils::numeric(FLERR, arg[5], false, lmp);
  // cutoff distance for active drag
  Lcut = utils::numeric(FLERR, arg[6], false, lmp);
}

int FixActiveDrag::setmask()
{
  int mask = 0;
  // mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

void FixActiveDrag::init()
{
  NeighRequest *req;
  if(neighbor->cutneighmax < Lcut)
  {
   
    if (comm->get_comm_cutoff() < Lcut) {
      error->all(FLERR,
                 "Communication cutoff is smaller than the cutoff of coupling cutoff distance.");
    }
    req = neighbor->add_request(this, REQ_OCCASIONAL | REQ_FULL);
    req->set_cutoff(Lcut);
  }
  else
    req = neighbor->add_request(this, REQ_FULL);


}

void FixActiveDrag::init_list(int, NeighList *ptr)
{
  list = ptr;
}

void FixActiveDrag::setup(int)
{
  cut_comm = comm->get_comm_cutoff();
  if (Lcut > cut_comm) error->all(FLERR, "Lcut is larger than the communication cutoff.");
}

// void FixActiveDrag::pre_force(int)
// {
//   // build an occasional list if the cutoff distance is larger than the maximum cutoff distance
 
// }

void FixActiveDrag::post_force(int)
{
  // Compute relative velocity and apply active drag force
  if(neighbor->cutneighmax < Lcut)
  {
    neighbor->build_one(list);
  }
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;

  int inum = list->inum;
  int nlocal = atom->nlocal;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int *ilist = list->ilist;

  int i, j, jnum;
  int *jlist;
  double xtmp, ytmp, ztmp, vxtmp, vytmp, vztmp;
  double delx, dely, delz, delvx, delvy, delvz, rsq, r, rinv, r2inv, dot;
  double nx, ny, nz;

  for (int ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      vxtmp = v[i][0];
      vytmp = v[i][1];
      vztmp = v[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      for (int jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        if (mask[j] & jgroupbit) {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < Lcut * Lcut) {
            r = sqrt(rsq);
            rinv = 1.0 / r;
            delvx = vxtmp - v[j][0];
            delvy = vytmp - v[j][1];
            delvz = vztmp - v[j][2];

            nx = delx * rinv;
            ny = dely * rinv;
            nz = delz * rinv;

            double del_vn = nx * delvx + ny * delvy + nz * delvz;

            double vtx = delvx - del_vn * nx;
            double vty = delvy - del_vn * ny;
            double vtz = delvz - del_vn * nz;

            // normal force
            double fnx = cn * del_vn * nx;
            double fny = cn * del_vn * ny;
            double fnz = cn * del_vn * nz;

            // tangential force
            double ftx = -ct * vtx;
            double fty = -ct * vty;
            double ftz = -ct * vtz;

            double fxtmp = fnx + ftx;
            double fytmp = fny + fty;
            double fztmp = fnz + ftz;
            f[i][0] += fxtmp;
            f[i][1] += fytmp;
            f[i][2] += fztmp;
            if (force->newton_pair || j < nlocal) {
              f[j][0] -= fxtmp;
              f[j][1] -= fytmp;
              f[j][2] -= fztmp;
            }
          }
        }
      }
    }
  }
}