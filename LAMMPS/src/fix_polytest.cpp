#include "atom.h"
#include "atom_masks.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "random_mars.h"
#include "update.h"

#include "domain.h"
#include "fix_polytest.h"
#include "force.h"
#include "lammps.h"
#include "pointers.h"
#include "polygon2D.h"
#include "update.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <tuple>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

static constexpr double SMALL = 0.001;

FixPolyTest::FixPolyTest(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
 
    px = utils::numeric(FLERR, arg[3], false, lmp);
    py = utils::numeric(FLERR, arg[4], false, lmp);
    pz = utils::numeric(FLERR, arg[5], false, lmp);
    rate = utils::numeric(FLERR, arg[6], false, lmp);    // add-atom rate
    seed = utils::numeric(FLERR, arg[7], false, lmp);    // random seed
    nevery = utils::numeric(FLERR, arg[8], false, lmp);
    rng = new RanMars(lmp, seed);
    force_reneighbor = 1;
    next_reneighbor = update->ntimestep + 1;
    dim = lmp->domain->dimension;
}

void FixPolyTest::getCOM()
{
  double cx, cy, cz;
  double **x = atom->x;
  int num = 0;
  double unwrap[3];
  for (int i = 0; i < atom->nlocal; ++i)    // careful, this is only capable for serial version
  {
    if (atom->mask[i] & groupbit) {
      domain->unmap(x[i], atom->image[i], unwrap);
      cx += unwrap[0];
      cy += unwrap[1];
      cz += unwrap[2];
      num++;
    }
  }
  cx /= num;
  cy /= num;
  cz /= num;    // COM of the membrane

  xc[0] = cx;
  xc[1] = cy;
  xc[2] = cz;
}

void FixPolyTest::deleteOneAtomAtLargestNegativePolarity(vector<double> &proj, vector<int> &ind)
{
  // proj: atom projection length along the polarity;
  // ind: stored atom indices in the same order as the atom projection;
  double minProj = proj[0];
  int minI = ind[0];
  for (int i = 0; i < proj.size(); ++i) {
    if (proj[i] < minProj) {
      minProj = proj[i];
      minI = ind[i];
    }
  }

  tagint delTag = atom->tag[minI];
  int *dlist = new int[atom->nlocal];
  for (int i = 0; i < atom->nlocal; ++i) {
    if (i == minI)
      dlist[i] = 1;
    else
      dlist[i] = 0;
  }

  int nlocal = atom->nlocal;
  int j = 0;
  while (j < nlocal) {
    if (dlist[j]) {
      atom->avec->copy(nlocal - 1, j, 1);
      dlist[j] = dlist[nlocal - 1];
      nlocal--;
    } else
      j++;
  }

  atom->nlocal = nlocal;
  //printf("deleted one atom of tag %d \n", delTag);

  // Note that if multiple minima exist, only delete the first one.
}

void FixPolyTest::addOneAtomAtLargestPositivePolarity(vector<double> &proj,
                                                                  vector<int> &ind)
{
  double maxProj = proj[0];
  int maxI = ind[0];
  for (int i = 0; i < proj.size(); ++i) {
    if (proj[i] > maxProj) {
      maxProj = proj[i];
      maxI = ind[i];
    }
  }

  double xnew[3];    // change here for more stable atom addition
  int itype, imask;
  double irad;    // radius or interaction range;
  imageint img;

  tagint tagMax = atom->tag[maxI];



  // get the atom type to create
  itype = atom->type[maxI];
  imask = atom->mask[maxI];
  if (atom->radius != nullptr)
    irad = atom->radius[maxI];
    else
    irad = 0.2;

  img = atom->image[maxI];
  xnew[0] = atom->x[maxI][0] - px * irad;
  xnew[1] = atom->x[maxI][1] - py * irad;
  xnew[2] = atom->x[maxI][2] - pz * irad;

  if (!domain->inside(xnew)) { domain->remap(xnew, img); }

  atom->avec->create_atom(itype, xnew);
  int n = atom->nlocal - 1;
  atom->tag[n] =
      0;    // by default, in new versions of LAMMPS, newly created atom automatically has a tag of 0
  atom->mask[n] = imask;
  if (atom->radius != nullptr) atom->radius[n] = irad;
  atom->image[n] = img;

  atom->tag_extend();
  atom->tag_check();

}



void FixPolyTest::init()
{
  // normalize polarity vector;
  double lp = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2));
  px /= lp;
  py /= lp;
  pz /= lp;
}



int FixPolyTest::setmask()
{
  int mask;
  mask |= PRE_EXCHANGE;    // add and delete atoms here
  return mask;
}

void FixPolyTest::pre_exchange()
{
  if (next_reneighbor != update->ntimestep) return;    // not time to polymerize or reneighbor

  double prob = 1 - exp(-rate * (update->dt) * nevery);
  if (rng->uniform() <= prob) {

    getCOM();

    vector<int> ind;
    vector<double> proj;
    double unwrap[3];
    double **x = atom->x;
    for (int i = 0; i < atom->nlocal; ++i) {
      if (atom->mask[i] & groupbit) {
        domain->unmap(x[i], atom->image[i], unwrap);
        double projtmp =
            px * (unwrap[0] - xc[0]) + py * (unwrap[1] - xc[1]) + pz * (unwrap[2] - xc[2]);
        proj.push_back(projtmp);
        ind.push_back(i);
      }
    }    // projection along polarity direction. Positive -> cell-front, Negative->cell-back

    // depolymerization happens at the position with minimal projection along the palority
    // delete one particle at the negative polarity direction

    deleteOneAtomAtLargestNegativePolarity(proj, ind);

    // polymerization happens at the position with maximal projection along the polarity
    // add one particle at the center of particles having the largest proj and the second largest proj.

    addOneAtomAtLargestPositivePolarity(
        proj,
        ind);    // addTags[0]: tag of the newly created atom; addTags[1]: the tag of the atom where you want to add an atom


    atom->map_delete();
    atom->map_init(1);
    atom->map_set();
  }

  next_reneighbor += nevery;
}

