#include "atom.h"
#include "atom_masks.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "random_mars.h"
#include "update.h"

#include "fix_addatom2d.h"
#include "force.h"
#include "polygon2D.h"
#include "update.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

FixAddAtom2D::FixAddAtom2D(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  px = utils::numeric(FLERR, arg[3], false, lmp);
  py = utils::numeric(FLERR, arg[4], false, lmp);
  rate = utils::numeric(FLERR, arg[5], false, lmp);
  rng = new RanMars(lmp, 123546);
  nevery = utils::numeric(FLERR, arg[6], false, lmp);
  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
  cout << "This fix has initialized " <<endl;
}

int FixAddAtom2D::setmask()
{
  int mask;
  mask |= PRE_EXCHANGE;
  return mask;
}

void FixAddAtom2D::pre_exchange()
{
  if(next_reneighbor != update->ntimestep) return;

  // for deleting atoms at the back 
  int nlocal = atom->nlocal;
  int *dlist = new int[nlocal];
  for(int i = 0; i < nlocal; ++i)
  {
    dlist[i] = 0;
  }

  // cout << "position to generate atom " << sortedVertices[closestIndex].first<<", "<<sortedVertices[closestIndex].second<<endl;
  // create an atom
  double prob = 1 - exp(-rate * (update->dt) * nevery);
 // cout<< "probability is %f " << prob << endl;
  if (rng->uniform() <= prob) {
    cout << "An atom should be added" <<endl;
    

    double ang0 = atan2(py, px);

    double cx = 0.0, cy = 0.0;

    vector<pair<double, double>> vertices;

    for (int i = 0; i < nlocal; ++i) {
      if (atom->mask[i] & groupbit) {
        vertices.push_back({atom->x[i][0], atom->x[i][1]});
        cx += atom->x[i][0];
        cy += atom->x[i][1];
      }
    }

    cx /= vertices.size();
    cy /= vertices.size();

    vector<pair<double, double>> sortedVertices = polygon2D::sortVertices(vertices);
    vector<double> angs;
    for (int i = 0; i < sortedVertices.size(); ++i) {
      angs.push_back(atan2(sortedVertices[i].second - cy, sortedVertices[i].first - cx));
    }

    int closestIndex = 0;
    int furtherestIndex = 0;
    double smallestDiff = abs(angs[0] - ang0);
    double largestDiff = abs(angs[0] - ang0);
    for (int i = 0; i < angs.size(); ++i) {
      double currentDiff = abs(angs[i] - ang0);
      if (currentDiff < smallestDiff) {
        closestIndex = i;
        smallestDiff = currentDiff;
      }

      if(currentDiff > largestDiff) {
        furtherestIndex = i;
        largestDiff = currentDiff;
      }
    }

    // First delete an atom at the back of the membrane

    for(int i = 0; i < nlocal; ++i)
    {
      if(atom->mask[i]&groupbit)
      {
        pair<double, double> currentVertex = {atom->x[i][0], atom->x[i][1]};
        if (currentVertex == sortedVertices[furtherestIndex])
        {
          dlist[i] = 1;
        } 
      }
    }

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


    nlocal = atom->nlocal;
    
    for (int i = 0; i < nlocal; ++i) {
      if (atom->mask[i] & groupbit) {
        pair<double, double> currentVertex = {atom->x[i][0], atom->x[i][1]};
        if (currentVertex == sortedVertices[closestIndex]) {
          int followIndex = closestIndex + 1;
          if (followIndex >= sortedVertices.size()) { followIndex = 0; }
          double *xadd = new double[3];
          xadd[0] = 0.5 * (sortedVertices[closestIndex].first + sortedVertices[followIndex].first);
          xadd[1] =
              0.5 * (sortedVertices[closestIndex].second + sortedVertices[followIndex].second);
          xadd[2] = 0.0;
          atom->avec->create_atom(atom->type[i], xadd);
          int n = atom->nlocal - 1;
          atom->tag[n] = 0;
          atom->mask[n] = atom->mask[i];
          atom->radius[n] = atom->radius[i];
          cout << "an Atom is added" <<endl;

          delete[] xadd;

          break;
        }
      }
    }
  }
  atom->tag_extend();
  atom->tag_check();

  next_reneighbor += nevery;
}
