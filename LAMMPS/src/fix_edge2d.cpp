#include "atom.h"
#include "atom_masks.h"
#include "error.h"
#include "group.h"
#include "update.h"
#include "random_mars.h"
#include "comm.h"

#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <algorithm>
#include "polygon2D.h"
#include "fix_edge2d.h"
#include "force.h"
#include "update.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;



FixEdge2D::FixEdge2D(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  K = utils::numeric(FLERR, arg[3], false, lmp);
  L0 = utils::numeric(FLERR, arg[4], false, lmp);
}

int FixEdge2D::setmask()
{
  int mask;
  mask |= POST_FORCE;
  return mask;
}

void FixEdge2D::post_force(int)
{
  pos = atom->x;
  nlocal = atom->nlocal;
  f = atom->f;
  mask = atom->mask;

  // check update


  vector<pair<double, double>> vertices;
  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) { vertices.push_back({pos[i][0], pos[i][1]}); }
  }
  // cout << "Compute " << vertices.size() <<" particles in membrane group." << endl;
  vector<pair<double, double>> sortedVertices = polygon2D::sortVertices(vertices); // vertices are sorted in a counterclockwise order
  // apply spring force
  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      pair<double, double> currentVertex = {pos[i][0], pos[i][1]};
      for (int j = 0; j < sortedVertices.size(); j++) {
        if (currentVertex == sortedVertices[j]) {
          int prej = j - 1;
          int aftj = j + 1;

          if (prej < 0) prej = sortedVertices.size() - 1;
          if (aftj >= sortedVertices.size()) aftj = 0;

          double prevEdgeLength = sqrt(pow(sortedVertices[prej].first - currentVertex.first, 2) +
                                       pow(sortedVertices[prej].second - currentVertex.second, 2));
          double aftrEdgeLength = sqrt(pow(sortedVertices[aftj].first - currentVertex.first, 2) +
                                       pow(sortedVertices[aftj].second - currentVertex.second, 2));

          f[i][0] += K * (prevEdgeLength - L0) *
                  (sortedVertices[prej].first - currentVertex.first) / prevEdgeLength +
              K * (aftrEdgeLength - L0) * (sortedVertices[aftj].first - currentVertex.first) /
                  aftrEdgeLength;
          f[i][1] += K * (prevEdgeLength - L0) *
                  (sortedVertices[prej].second - currentVertex.second) / prevEdgeLength +
              K * (aftrEdgeLength - L0) * (sortedVertices[aftj].second - currentVertex.second) /
                  aftrEdgeLength;
        }
      }
    }
  }
}












