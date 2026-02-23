#include "fix_tdpd_fg_reaction.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "region.h"
#include "group.h"

#include <cmath>
#include <cstring>
#include <algorithm>
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;


/* ---------------------------------------------------------------------- */
constexpr double EPSILON = 1e-10;

FixTDPDFGReaction::FixTDPDFGReaction(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  
  int iarg = 3;
  cc_index = utils::inumeric(FLERR, arg[iarg++], false, lmp);
  rho = utils::numeric(FLERR, arg[iarg++], false, lmp);
  K = utils::numeric(FLERR, arg[iarg++], false, lmp);
  n = utils::numeric(FLERR, arg[iarg++], false, lmp);
  mu = utils::numeric(FLERR, arg[iarg++], false, lmp); 
  cc_index_pool = cc_index; cc_index_inhibitor = cc_index;
  npool = n; Kpool = K; Kinhibitor = K; ninhibitor = n;
  jgroupbit = groupbit;
  if(iarg < narg)
  {
    cc_index_pool = utils::inumeric(FLERR, arg[iarg++], false, lmp);
    Kpool = utils::numeric(FLERR, arg[iarg++], false, lmp);
    npool = utils::numeric(FLERR, arg[iarg++], false, lmp);
    if(iarg < narg)
    {
      int jgroup = group->find(arg[iarg++]); // local source group
      jgroupbit = group->bitmask[jgroup];
      if (jgroup < 0) error->all(FLERR, "Invalid group ID");
      if(iarg < narg) 
      {
        cc_index_inhibitor = utils::inumeric(FLERR, arg[iarg++], false, lmp);
        Kinhibitor = utils::numeric(FLERR, arg[iarg++], false, lmp);
        ninhibitor = utils::numeric(FLERR, arg[iarg++], false, lmp);
      }
    }
  }
    
}

/* ---------------------------------------------------------------------- */

FixTDPDFGReaction::~FixTDPDFGReaction()
{
}

/* ---------------------------------------------------------------------- */

int FixTDPDFGReaction::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void FixTDPDFGReaction::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **cc_flux = atom->cc_flux;
  double **cc = atom->cc;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;


  int k1 = cc_index - 1;
  int k2 = cc_index_pool - 1;
  int k3 = cc_index_inhibitor - 1;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      cc_flux[i][k1] -= mu * cc[i][k1];
      cc_flux[i][k2] += mu * cc[i][k1];
    }

    if (mask[i] & jgroupbit) {
      double R = rho * (pow(cc[i][k2], npool) / (pow(Kpool, npool) + pow(cc[i][k2], npool))) *
          (pow(Kinhibitor, ninhibitor) / (pow(Kinhibitor, ninhibitor) + pow(cc[i][k3], ninhibitor)))*cc[i][k2] * cc[i][k2];
      cc_flux[i][k1] += R;
      cc_flux[i][k2] -= R;
    }
  }
}
