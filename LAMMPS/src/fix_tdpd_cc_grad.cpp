#include "fix_tdpd_cc_grad.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair_tdpd.h"
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

FixTDPDCCGrad::FixTDPDCCGrad(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), list(nullptr), pair(nullptr), pair_tdpd(nullptr)
{
  if (narg < 4) {
    error->warning(FLERR,
                   "fix tdpd/cc/grad will compute the concentration gradients of all species");
    mode = 0;
  } else {
    mode = 1;
    cc_index_size = narg - 3;
    cc_index = new int[cc_index_size];
    for (int i = 0; i < narg - 3; i++) {
      cc_index[i] = utils::inumeric(FLERR, arg[i + 3], false, lmp);
    }
  }
}

FixTDPDCCGrad::~FixTDPDCCGrad() {}

int FixTDPDCCGrad::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

void FixTDPDCCGrad::init()
{
  // Check if the pair style is of type PairTDPD
  pair = force->pair;
  if (pair == nullptr) error->all(FLERR, "Pair style tdpd not found");
  pair_tdpd = dynamic_cast<PairTDPD *>(pair);
  if (pair_tdpd == nullptr) error->all(FLERR, "Pair style tdpd not found");


  

  NeighRequest *req = neighbor->add_request(this, REQ_FULL);

  int* mask = atom->mask;
  double **cc_grad = atom->cc_grad;

  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      for (int kk = 0; kk < atom->cc_species; kk++) {
        cc_grad[i][3 * kk] = 0.0;
        cc_grad[i][3 * kk + 1] = 0.0;
        cc_grad[i][3 * kk + 2] = 0.0;
      }
    }
  }
}

void FixTDPDCCGrad::init_list(int, NeighList *ptr)
{
  list = ptr;
}

void FixTDPDCCGrad::pre_force(int /*vflag*/)
{
  switch (mode) {
    case 0:
      for (int k = 0; k < atom->cc_species; k++) {
        compute_cc_grad(k);
      }
      break;
    case 1:
      for(int l = 0; l < cc_index_size; l++)
      {
        int k = cc_index[l] - 1;
        compute_cc_grad(k);
      }
      break;
  }
}

void FixTDPDCCGrad::compute_cc_grad(int k)
{
    int *mask = atom->mask;
    double **cc = atom->cc;
    double **cc_grad = atom->cc_grad;
    double **x = atom->x;
    int *type = atom->type;



    double **cutsq = pair_tdpd->cutsq;
    double **cutcc = pair_tdpd->cutcc;
    double ***powercc = pair_tdpd->powercc;

    int inum = list->inum;
    int *ilist = list->ilist;
    int *numneigh = list->numneigh;
    int **firstneigh = list->firstneigh;


    //FILE *f = fopen("cc_grad.txt", "w");

    for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii];
        if (mask[i] & groupbit) {
          double xtmp = x[i][0];
          double ytmp = x[i][1];
          double ztmp = x[i][2];
          int itype = type[i];
          int *jlist = firstneigh[i];
          int jnum = numneigh[i];
          double w = 0,gx = 0,gy = 0,gz = 0;
          for (int jj = 0; jj < jnum; jj++) {
            int j = jlist[jj];
            j &= NEIGHMASK;
            if (mask[j] & groupbit) {
              double delx = xtmp - x[j][0];
              double dely = ytmp - x[j][1];
              double delz = ztmp - x[j][2];
              double rsq = delx * delx + dely * dely + delz * delz;
              int jtype = type[j];
              if (rsq < cutcc[itype][jtype]*cutcc[itype][jtype]) {
                double r = sqrt(rsq);
                if (r < EPSILON) continue;    // r can be 0.0 in DPD systems
                // chemical concentration transport
                if (r < cutcc[itype][jtype]) { 
      

                  double wd_prime = -powercc[itype][jtype][k] / cutcc[itype][jtype] *
                   pow(1.0 - r / cutcc[itype][jtype],powercc[itype][jtype][k] - 1.0);
                  

                  double d=  wd_prime*(cc[j][k] -cc[i][k]);
                  gx += d * (delx / r);
                  gy += d *( dely / r) ;
                  gz += d * (delz / r) ;
                  // w += wd_prime;
                }
              }
            }
          }

          // gx = gx / MAX(w, EPSILON);
          // gy = gy / MAX(w, EPSILON);
          // gz = gz / MAX(w, EPSILON);

          cc_grad[i][3*k] = gx;
          cc_grad[i][3*k + 1] = gy;
          cc_grad[i][3*k + 2] = gz;
          //fprintf(f, "%d, %f,%f,%f,%f,%f\n", atom->tag[i], x[i][0], x[i][1], cc[i][k], cc_grad[i][3*k], cc_grad[i][3*k + 1]);

         // std::cout << x[i][0] << " " << x[i][1] << " " << cc[i][k] << " " << cc_grad[i][3*k] << " " << cc_grad[i][3*k + 1] << std::endl;
        }
      }
}