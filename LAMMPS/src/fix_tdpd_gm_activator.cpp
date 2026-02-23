#include "fix_tdpd_gm_activator.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair_tdpd.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace NeighConst;

#define EPSILON 1.0e-10
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MAX(A, B) ((A) > (B) ? (A) : (B))

/* ---------------------------------------------------------------------- */

FixTDPDGMActivator::FixTDPDGMActivator(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix tdpd/source", error);

  int iarg = 3;
  cc_index1 = utils::inumeric(FLERR, arg[iarg++], false,
                             lmp);    // species of the polarity cue
  cc_index2 = utils::inumeric(
      FLERR, arg[iarg++], false, lmp);  // species of the inhibitor
  
  rho0 = utils::numeric(FLERR, arg[iarg++], false, lmp);    // constant production
  rhomax = utils::numeric(FLERR, arg[iarg++], false, lmp);    // the maximum production rate of the activator
  Ka = utils::numeric(FLERR, arg[iarg++], false, lmp);    // the half saturation constant for the activator
  na = utils::numeric(FLERR, arg[iarg++], false, lmp);    // the Hill coefficient for the activator

  Kb = utils::numeric(FLERR, arg[iarg++], false, lmp);    // the half saturation constant for the inhibitor
  nb = utils::numeric(FLERR, arg[iarg++], false, lmp);    // the Hill coefficient for the inhibitor
  mu0 = utils::numeric(FLERR, arg[iarg++], false, lmp);    // the base decay rate of the activator
  mumax = utils::numeric(FLERR, arg[iarg++], false, lmp);    // the decay rate of the activator
}

/* ---------------------------------------------------------------------- */

FixTDPDGMActivator::~FixTDPDGMActivator()
{
}

/* ---------------------------------------------------------------------- */

int FixTDPDGMActivator::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */




/* ---------------------------------------------------------------------- */

void FixTDPDGMActivator::post_force(int /*vflag*/)
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **cc_flux = atom->cc_flux;
  double **cc = atom->cc;

  int k1 = cc_index1 - 1; // activated
  int k2 = cc_index2 - 1; // inhibitor


  for(int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double R = rho0 + rhomax * pow(cc[i][k1], na) / (pow(Ka,na) + pow(cc[i][k1], na)); // the production rate of the activator

      double mu = mumax * pow(cc[i][k2],nb) / (pow(Kb,nb) + pow(cc[i][k2],nb)) + mu0;    // the inhibition effect of the inhibitor
      
      cc_flux[i][k1] += R - mu * cc[i][k1];
    }
  }


}