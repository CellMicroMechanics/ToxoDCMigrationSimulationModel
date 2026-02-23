#include "fix_tdpd_force_contractile.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
// #include "neigh_list.h"
// #include "neigh_request.h"
// #include "neighbor.h"
// #include "pair.h"
// #include "pair_tdpd.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
// using namespace NeighConst;

using namespace std;

/* ---------------------------------------------------------------------- */
#define EPSILON 1.0e-10
FixTDPDForceContractile::FixTDPDForceContractile(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{

  int iarg = 3;
  cc_index = utils::inumeric(FLERR, arg[iarg++], false, lmp);
  zeta = utils::numeric(FLERR, arg[iarg++], false, lmp);
}

int FixTDPDForceContractile::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

// void FixTDPDForceContractile::init()
// {
//   // Use half neighbor list; we apply equal/opposite forces with newton checks
//   neighbor->add_request(this);
//   // Check if the pair style is of type PairTDPD
//   pair = force->pair;
//   if (pair == nullptr) error->all(FLERR, "Pair style tdpd not found");
//   pair_tdpd = dynamic_cast<PairTDPD *>(pair);
//   if (pair_tdpd == nullptr) error->all(FLERR, "Pair style tdpd not found");
// }

// void FixTDPDForceContractile::init_list(int /*id*/, NeighList *ptr /*arg*/)
// {
//   list = ptr;
// }

/* ---------------------------------------------------------------------- */

void FixTDPDForceContractile::post_force(int /*vflag*/)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **f = atom->f;
  // double **cc = atom->cc;
  // double **x = atom->x;
  // int *type = atom->type;

  double **cc_grad = atom->cc_grad;

  int k = cc_index - 1;    // species index

  for(int i = 0; i < nlocal; i++)
  {
    if(mask[i] & groupbit)
    {
      f[i][0] -= zeta * cc_grad[i][3*k];
      f[i][1] -= zeta * cc_grad[i][3*k + 1];
      f[i][2] -= zeta * cc_grad[i][3*k + 2];
    }
  }

  // int inum = list->inum;
  // int *ilist = list->ilist;
  // int *numneigh = list->numneigh;
  // int **firstneigh = list->firstneigh;

  // const int newton_pair = force->newton_pair;

  // double **cutsq = pair_tdpd->cutsq;
  // double **cutcc = pair_tdpd->cutcc;
  // double ***powercc = pair_tdpd->powercc;

  // for (int ii = 0; ii < inum; ii++) {
  //   int i = ilist[ii];
  //   if (mask[i] & groupbit) {
  //     double xtmp = x[i][0];
  //     double ytmp = x[i][1];
  //     double ztmp = x[i][2];
  //     int itype = type[i];
  //     int *jlist = firstneigh[i];
  //     int jnum = numneigh[i];

  //     double Pi = -zeta * cc[i][k];

  //     for (int jj = 0; jj < jnum; jj++) {
  //       int j = jlist[jj];
  //       j &= NEIGHMASK;

  //       double fxtmp = 0.0;
  //       double fytmp = 0.0;
  //       double fztmp = 0.0;

  //       if (mask[j] & groupbit) {
  //         double Pj = -zeta * cc[j][k];
  //         double delx = xtmp - x[j][0];
  //         double dely = ytmp - x[j][1];
  //         double delz = ztmp - x[j][2];
  //         double rsq = delx * delx + dely * dely + delz * delz;
  //         int jtype = type[j];

  //         if (rsq < cutcc[itype][jtype] * cutcc[itype][jtype]) {
  //           double r = sqrt(rsq);
  //           if (r < EPSILON) continue;

  //           double wd_prime = -powercc[itype][jtype][k] / cutcc[itype][jtype] *
  //               pow(1.0 - r / cutcc[itype][jtype], powercc[itype][jtype][k] - 1.0);

  //           fxtmp += wd_prime * (delx / r) * (Pi + Pj);
  //           fytmp += wd_prime * (dely / r) * (Pi + Pj);
  //           fztmp += wd_prime * (delz / r) * (Pi + Pj);
  //         }

  //         f[i][0] += fxtmp;
  //         f[i][1] += fytmp;
  //         f[i][2] += fztmp;

  //         // For half neighbor lists only add to owned partner atoms
  //         if (newton_pair || j < nlocal) {
  //           f[j][0] -= fxtmp;
  //           f[j][1] -= fytmp;
  //           f[j][2] -= fztmp;
  //         }
  //       }
  //     }
  //   }
  // }
}
