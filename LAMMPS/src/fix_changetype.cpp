#include "atom.h"
#include "atom_masks.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "random_mars.h"
#include "update.h"

#include "fix_changetype.h"
// #include "neighbor.h"
// #include "neigh_list.h"
// #include "neigh_request.h"
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

FixChangeType::FixChangeType(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  jgroup = group->find(arg[3]);
  jgroupbit = group->bitmask[jgroup];

  px = utils::numeric(FLERR, arg[4], false, lmp);
  py = utils::numeric(FLERR, arg[5], false, lmp);

  rate = utils::numeric(FLERR, arg[6], false, lmp);
  nevery = utils::numeric(FLERR, arg[7], false, lmp);

  rng = new RanMars(lmp, 123546);
  // force_reneighbor = 1;
  // next_reneighbor = update->ntimestep + 1;
  changeFlag = false;

  // neighbor->add_request(this,NeighConst::REQ_FULL);
}

// void FixChangeType::init_list(int, NeighList* ptr)
// {
//   list = ptr;
// }

int FixChangeType::setmask()
{
  int mask;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

void FixChangeType::deleteOutsideParticles()
{
  int nlocal = atom->nlocal;
  int *dlist = new int[nlocal];

  for (int i = 0; i < nlocal; ++i) { dlist[i] = 0; }
  vector<pair<double, double>> vertices;    // membrane vertices
  for (int i = 0; i < nlocal; ++i) {
    if (atom->mask[i] & jgroupbit) { vertices.push_back({atom->x[i][0], atom->x[i][1]}); }
  }

  // counterclockwise-sorted membrane particles
  vector<pair<double, double>> sortedVertices = polygon2D::sortVertices(vertices);
  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      pair<double, double> point = {atom->x[i][0], atom->x[i][1]};
      if (!polygon2D::isPointInsidePolygon(point, sortedVertices)) { dlist[i] = 1; }
    }
  }

  int i = 0;
  while (i < nlocal) {
    if (dlist[i]) {
      atom->avec->copy(nlocal - 1, i, 1);
      dlist[i] = dlist[nlocal - 1];
      nlocal--;
    } else
      i++;
  }

  atom->nlocal = nlocal;
}



void FixChangeType::pre_force(int)
{
  if(update->ntimestep % nevery) return;
  /*DELETE UNWANTED ATOMS*/
  // deleteOutsideParticles();
 
  // normalize polarity
  double px_norm = px / sqrt(pow(px, 2) + pow(py, 2));
  double py_norm = py / sqrt(pow(px, 2) + pow(py, 2));
  // switch an atom from type of jgroup to type of group
  double prob = 1 - exp(-rate * (update->dt) * nevery);
  if (rng->uniform() <= prob) {
    changeFlag = true;
    
    double cx,cy;
    vector<pair<double,double>> vertices; // cytoplasmic particles
    for(int i = 0; i < atom->nlocal; ++i)
    {
      if(atom->mask[i] & groupbit)
      {
        vertices.push_back({atom->x[i][0],atom->x[i][1]});
        cx += atom->x[i][0];
        cy += atom->x[i][1];
      }
    }

    cx /= vertices.size();
    cy /= vertices.size();


    
    for (int i = 0; i < atom->nlocal; ++i) {
      if (atom->mask[i] & jgroupbit) {
        toType = atom->type[i];
        toMask = atom->mask[i];
        toR = atom->radius[i];
        break;
      }
    }

    double maxProjectionAtPolarity = 0.0;
    int maxProjIndex;

    for (int i = 0; i < vertices.size(); ++i) {

      double projectionAtPolarity =
          (vertices[i].first - cx) * px_norm + (vertices[i].second - cy) * py_norm;
      if (projectionAtPolarity > maxProjectionAtPolarity) {
        maxProjectionAtPolarity = projectionAtPolarity;
        maxProjIndex = i;
      }
    }

    for (int i = 0; i < atom->nlocal; ++i) {
      if (atom->mask[i] & groupbit) {
        pair<double, double> currentVertex = {atom->x[i][0], atom->x[i][1]};
        if (currentVertex == vertices[maxProjIndex]) {
          atom->type[i] = toType;
          atom->mask[i] = toMask;
        }
      }
    }

    double cx2, cy2;
    // int nlocal = atom->nlocal;
    vector<pair<double, double>> vertices2tmp;
    for (int i = 0; i < atom->nlocal; ++i) {
      if (atom->mask[i] & jgroupbit) {
        vertices2tmp.push_back({atom->x[i][0], atom->x[i][1]});
        cx2 += atom->x[i][0];
        cy2 += atom->x[i][1];
      }
    }

    cx2 /= vertices2tmp.size();
    cy2 /= vertices2tmp.size();

    vector<pair<double, double>> vertices2 = polygon2D::sortVertices(vertices2tmp);

    double minProjectionAtPolarity = 0.0;
    int minProjIndex;

    for (int i = 0; i < vertices2.size(); ++i) {

      double projectionAtPolarity =
          (vertices2[i].first - cx2) * px_norm + (vertices2[i].second - cy2) * py_norm;
      if (projectionAtPolarity < minProjectionAtPolarity) {
        minProjectionAtPolarity = projectionAtPolarity;
        minProjIndex = i;
      }
    }

    for (int i = 0; i < atom->nlocal; ++i) {
      if (atom->mask[i] & groupbit) {
        toType2 = atom->type[i];
        toMask2 = atom->mask[i];
        toR2 = atom->radius[i];
        break;
      }
    }

    // change a membrane particle to cytoplasma particle

    for (int i = 0; i < atom->nlocal; ++i) {
      if (atom->mask[i] & jgroupbit) {
        pair<double, double> currentVertex = {atom->x[i][0], atom->x[i][1]};
        if (currentVertex == vertices2[minProjIndex]) {
          atom->type[i] = toType2;
          atom->mask[i] = toMask2;
          cout<< "Changed atom index is " << i << endl;
          // atom->f[i][0] += 1000;
          // // check if the changed particle is out side of the membrane
          // int preMinIndex = minProjIndex - 1;
          // int aftMinIndex = minProjIndex + 1;

          // if (preMinIndex < 0) preMinIndex = vertices2.size() - 1;
          // if (aftMinIndex > vertices2.size() - 1) aftMinIndex = 0;

          // double dx = vertices2[aftMinIndex].first - vertices2[preMinIndex].first;
          // double dy = vertices2[aftMinIndex].second - vertices2[preMinIndex].second;

          // pair<double, double> norm_in = {-dy / sqrt(pow(dx, 2) + pow(dy, 2)),
          //                                 dx / sqrt(pow(dx, 2) + pow(dy, 2))};

          // pair<double, double> aux = {vertices2[aftMinIndex].first - vertices2[minProjIndex].first,
          //                             vertices2[aftMinIndex].second -
          //                                 vertices2[minProjIndex].second};

          // double l = aux.first * norm_in.first + aux.second * norm_in.second;
          // double dd = 1e-3;
          // cout << l << endl;
          // if (l >= 0) 
          // {
          //   cout << "wanring: a particle goes out" << endl;
          // }
          // while (l >= 0) {
          //   vertices2[minProjIndex].first += dd * norm_in.first;
          //   vertices2[minProjIndex].second += dd * norm_in.second;
          //   aux = {vertices2[aftMinIndex].first - vertices2[minProjIndex].first,
          //          vertices2[aftMinIndex].second - vertices2[minProjIndex].second};
          //   l = aux.first * norm_in.first + aux.second * norm_in.second;
          // }
          // atom->x[i][0] = vertices2[minProjIndex].first;
          // atom->x[i][1] = vertices2[minProjIndex].second;
        }
      }
    }
  }
}


void FixChangeType::post_force(int)
{
  vector<pair<double,double>> vertices; // new membrane
  for(int i = 0; i < atom->nlocal; ++i)
  {
    if(atom->mask[i] & jgroupbit)
    {
      vertices.push_back({atom->x[i][0],atom->x[i][1]});
    }
  }
  vector<pair<double,double>> sortedVertices = polygon2D::sortVertices(vertices);

  for(int i = 0; i < atom->nlocal; ++i)
  {
    if(atom->mask[i] & groupbit)
    {
      pair<double,double> point = {atom->x[i][0],atom->x[i][1]};
      vector<double> dist = polygon2D::computeDistanceToEdgeCenter(point,sortedVertices);
      auto minIt = min_element(dist.begin(),dist.end());
      double minDist = *minIt;
      if( minDist <= toR + toR2)
      {
        int minIdx = distance(dist.begin(),minIt);
        int nextIdx = (minIdx + 1) % sortedVertices.size();
        double dx = sortedVertices[nextIdx].first - sortedVertices[minIdx].first;
        double dy = sortedVertices[nextIdx].second - sortedVertices[minIdx].second;
        pair<double,double> norm_in = {-dy / sqrt(pow(dx, 2) + pow(dy, 2)),dx / sqrt(pow(dx, 2) + pow(dy, 2))};
        if (polygon2D::isPointInsidePolygon(point, sortedVertices)) {
          atom->f[i][0] += 20 * abs(minDist - toR - toR2) * norm_in.first;
          atom->f[i][1] += 20 * abs(minDist - toR - toR2) * norm_in.second;
        } else {
          atom->f[i][0] += 20 * abs(minDist + toR + toR2) * norm_in.first;
          atom->f[i][1] += 20 * abs(minDist + toR + toR2) * norm_in.second;
        }
      }
    }
  }

}