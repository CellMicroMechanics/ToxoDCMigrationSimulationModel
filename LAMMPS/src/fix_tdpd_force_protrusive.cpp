#include "fix_tdpd_force_protrusive.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "pair_tdpd.h"
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace NeighConst;
using namespace std;

/* ---------------------------------------------------------------------- */
#define EPSILON 1.0e-10
FixTDPDForceProtrusive::FixTDPDForceProtrusive(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), list(nullptr), pair(nullptr), pair_tdpd(nullptr)
{

  int iarg = 3;
  cc_index = utils::inumeric(FLERR, arg[iarg++],false, lmp);
  zeta = utils::numeric(FLERR, arg[iarg++],false, lmp);

}



void FixTDPDForceProtrusive::init()
{
  neighbor->add_request(this,REQ_FULL);
   // Check if the pair style is of type PairTDPD
  pair = force->pair;
  if (pair == nullptr) error->all(FLERR, "Pair style tdpd not found");
  pair_tdpd = dynamic_cast<PairTDPD *>(pair);
  if (pair_tdpd == nullptr) error->all(FLERR, "Pair style tdpd not found");
}

int FixTDPDForceProtrusive::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

void FixTDPDForceProtrusive::init_list(int /*id*/, NeighList *ptr /*arg*/)
{
  list = ptr;
}



/* ---------------------------------------------------------------------- */

void FixTDPDForceProtrusive::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **cc = atom->cc;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **f = atom->f;
  int *type = atom->type;

  double **cc_grad = atom->cc_grad;

  int k = cc_index - 1; // species index
  // x component 3*k, y component 3*k+1, z component 3*k+2

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;


  double **cutsq = pair_tdpd->cutsq;
  double **cutcc = pair_tdpd->cutcc;
  double ***powercc = pair_tdpd->powercc;


  // compute force from the active stress
  for(int ii = 0; ii < inum; ii++)
  {
    int i = ilist[ii];
    if(mask[i] & groupbit) // apply to cytosol particles
    {
      double xtmp = x[i][0];
      double ytmp = x[i][1];
      double ztmp = x[i][2];
      int itype = type[i];
      int *jlist = firstneigh[i];
      int jnum = numneigh[i];
      
      double Ai = -zeta * cc[i][k] / max((cc_grad[i][3*k] * cc_grad[i][3*k] +
                                      cc_grad[i][3*k + 1] * cc_grad[i][3*k + 1] +
                                      cc_grad[i][3*k + 2] * cc_grad[i][3*k + 2]),EPSILON);

      double fxtmp = 0.0;
      double fytmp = 0.0;
      double fztmp = 0.0;

      for(int jj = 0; jj < jnum; jj++)
      {
        int j = jlist[jj];
        j &= NEIGHMASK;
    
        if(mask[j] & groupbit)
        {
          double delx = xtmp - x[j][0];
          double dely = ytmp - x[j][1];
          double delz = ztmp - x[j][2];
          double rsq = delx * delx + dely * dely + delz * delz;
          int jtype = type[j];

          if(rsq < cutcc[itype][jtype] * cutcc[itype][jtype])
          {
            double r = sqrt(rsq);
            if(r < EPSILON) continue;

              double wd_prime = -powercc[itype][jtype][k] / cutcc[itype][jtype] *
                   pow(1.0 - r / cutcc[itype][jtype],powercc[itype][jtype][k] - 1.0);

              double Aj = -zeta * cc[j][k] / max((cc_grad[j][3*k] * cc_grad[j][3*k] +
                                              cc_grad[j][3*k + 1] * cc_grad[j][3*k + 1] +
                                              cc_grad[j][3*k + 2] * cc_grad[j][3*k + 2]),EPSILON);

              
             fxtmp += wd_prime * (delx / r) * (Ai*cc_grad[i][3*k]*cc_grad[i][3*k] - Aj*cc_grad[j][3*k]*cc_grad[j][3*k]) 
             + wd_prime * (dely / r) * (Ai*cc_grad[i][3*k + 1]*cc_grad[i][3*k] - Aj*cc_grad[j][3*k + 1]*cc_grad[j][3*k])
             + wd_prime * (delz / r) * (Ai*cc_grad[i][3*k + 2]*cc_grad[i][3*k] - Aj*cc_grad[j][3*k + 2]*cc_grad[j][3*k]);
              
              fytmp += wd_prime * (delx / r) * (Ai*cc_grad[i][3*k]*cc_grad[i][3*k + 1] - Aj*cc_grad[j][3*k]*cc_grad[j][3*k + 1])
              + wd_prime * (dely / r) * (Ai*cc_grad[i][3*k + 1]*cc_grad[i][3*k + 1] - Aj*cc_grad[j][3*k + 1]*cc_grad[j][3*k + 1])
              + wd_prime * (delz / r) * (Ai*cc_grad[i][3*k + 2]*cc_grad[i][3*k + 1] - Aj*cc_grad[j][3*k + 2]*cc_grad[j][3*k + 1]);

              fztmp += wd_prime * (delx / r) * (Ai*cc_grad[i][3*k]*cc_grad[i][3*k + 2] - Aj*cc_grad[j][3*k]*cc_grad[j][3*k + 2])
              + wd_prime * (dely / r) * (Ai*cc_grad[i][3*k + 1]*cc_grad[i][3*k + 2] - Aj*cc_grad[j][3*k + 1]*cc_grad[j][3*k + 2])
              + wd_prime * (delz / r) * (Ai*cc_grad[i][3*k + 2]*cc_grad[i][3*k + 2] - Aj*cc_grad[j][3*k + 2]*cc_grad[j][3*k + 2]);
            }
        }
      }
    
      f[i][0] += fxtmp;
      f[i][1] += fytmp;
      f[i][2] += fztmp;
    }
  }




  






 
}
