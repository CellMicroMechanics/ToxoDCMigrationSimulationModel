#include "atom.h"
#include "atom_masks.h"
#include "atom_vec.h"
#include "error.h"
#include "group.h"

#include "comm.h"
#include "update.h"

#include "force.h"

#include "update.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>

#include "domain.h"
#include "lammps.h"
#include "pointers.h"

#include "fix_nucleus2d.h"
#include "polygon2D.h"
#include <tuple>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

static constexpr double SMALL = 0.001;

FixNucleus2D::FixNucleus2D(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  K = utils::numeric(FLERR, arg[3], false, lmp);
  L0 = utils::numeric(FLERR, arg[4], false, lmp);
  KC = utils::numeric(FLERR, arg[5], false, lmp);
  theta0 = utils::numeric(FLERR, arg[6], false, lmp);
  if(narg == 8)
    theta_thres = utils::numeric(FLERR, arg[7], false, lmp);
  else if(narg == 9)
  {
    theta_thres = 0.0;
    KA = utils::numeric(FLERR, arg[7], false, lmp);
    A0 = utils::numeric(FLERR, arg[8], false, lmp);
  }
  //nevery = utils::numeric(FLERR, arg[8], false, lmp);
}

void FixNucleus2D::initializeNuclearBonds()
{
  // initialize 2d bonds between nuclear particles
  vector<tuple<double, double, tagint>> vertices;
  double  unwrap[3];
  for (int i = 0; i < atom->nlocal; ++i)    // only works for serial program
  {
    if (atom->mask[i] & groupbit) {
      domain->unmap(atom->x[i], atom->image[i], unwrap);
      vertices.push_back({unwrap[0],unwrap[1],atom->tag[i]});
    }
  }
  vector<tuple<double,double,tagint>> sortedVertices = polygon2D::sortVertices(vertices);

  if (bondedTags.empty()) {
    for (int i = 0; i < sortedVertices.size(); ++i) {
      int j = (i + 1) % sortedVertices.size();
      int k = (i + sortedVertices.size() - 1) % sortedVertices.size();
      bondedTags.push_back({get<2>(sortedVertices[i]), get<2>(sortedVertices[k]), get<2>(sortedVertices[j])});
    }
  }

  for (int i = 0; i < bondedTags.size(); ++i) { isBreak.push_back(false); }
}


void FixNucleus2D::init()
{
  theta0 *= M_PI / 180;
  theta_thres *= M_PI / 180;
}

void FixNucleus2D::setup(int vflag)
{
  initializeNuclearBonds();
}

int FixNucleus2D::setmask()
{
  int mask;
  mask |= POST_FORCE;
  return mask;
}

void FixNucleus2D::post_force(int)
{
  double dx1, dx2, dy1, dy2, dz1, dz2, theta, l1, l2, l1sq, l2sq, stheta, ctheta;
  //if(update->ntimestep % nevery != 0) return;
  double **f = atom->f;

/* Prepare for Area elasticity */
  vector<pair<double, double>> vertices;
  double unwrap[3];
  for (int i = 0; i < atom->nlocal; ++i)    // only works for serial program
  {
    if (atom->mask[i] & groupbit) {
      domain->unmap(atom->x[i],atom->image[i],unwrap);
      vertices.push_back({unwrap[0], unwrap[1]});
    }
  }
  double currentArea = polygon2D::polyArea2d(vertices);
  vertices.clear();

  double dA = currentArea - A0;
  double fA = KA * dA;
/*****************************/

  double unwrap_next[3],unwrap_prev[3];
  for (int i = 0; i < bondedTags.size(); ++i) {
    int nextI = atom->map(bondedTags[i][2]);
    int thisI = atom->map(bondedTags[i][0]);
    int prevI = atom->map(bondedTags[i][1]);
    if(nextI == -1 || thisI == -1 || prevI == -1) 
    {
      printf("WARNING: Bonded atoms are not found, skipping the bond \n");
      continue;
    }
    //double dx1,dx2,dy1,dy2,dz1,dz2,theta,l1,l2,l1sq,l2sq,stheta,ctheta;

    dx1 = atom->x[prevI][0] - atom->x[thisI][0];
    dy1 = atom->x[prevI][1] - atom->x[thisI][1];
    dz1 = atom->x[prevI][2] - atom->x[thisI][2];
    domain->minimum_image(dx1, dy1, dz1);

    dx2 = atom->x[nextI][0] - atom->x[thisI][0];
    dy2 = atom->x[nextI][1] - atom->x[thisI][1];
    dz2 = atom->x[nextI][2] - atom->x[thisI][2];
    domain->minimum_image(dx2, dy2, dz2);

    l1sq = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
    l2sq = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;

    l1 = sqrt(l1sq);
    l2 = sqrt(l2sq);
    // printf("distance of the two edges l1,l2: %f, %f \n", l1,l2);
    ctheta = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (l1 * l2);    // cos
    // printf("distance of the two edges l1,l2 after cos calcuation: %f, %f \n", l1,l2);
    // printf("cosine: %f \n", ctheta);

    if (KC != 0) {
      if (ctheta > 1) ctheta = 1;
      if (ctheta < -1) ctheta = -1;
      // break the bending if the curvature is too large
      if (ctheta > cos(theta_thres) && update->ntimestep != 1)
        isBreak[i] = true;
      else
        isBreak[i] = false;    // instant bond recovery

      if (!isBreak[i]) {
        stheta = sqrt(1 - ctheta * ctheta);    // sin  // D[cos] = - 1/sin
        if (stheta < SMALL) stheta = SMALL;
        double dtheta = acos(ctheta) - theta0;

        double fmag = KC * dtheta;    //

        double a = -fmag / stheta;    // prefactor
        // printf("Bending resistance force magnitude %f \n", fmag);

        double a11, a12, a22;

        a11 = a * ctheta / l1sq;
        a12 = -a / (l1 * l2);
        a22 = a * ctheta / l2sq;    // matrix form of the force

        f[thisI][0] -= (a11 + a12) * dx1 + (a12 + a22) * dx2;
        f[thisI][1] -= (a11 + a12) * dy1 + (a12 + a22) * dy2;
        f[thisI][2] -= (a11 + a12) * dz1 + (a12 + a22) * dz2;

        f[prevI][0] += a11 * dx1 + a12 * dx2;
        f[prevI][1] += a11 * dy1 + a12 * dy2;
        f[prevI][2] += a11 * dz1 + a12 * dz2;

        f[nextI][0] += a12 * dx1 + a22 * dx2;
        f[nextI][1] += a12 * dy1 + a22 * dy2;
        f[nextI][2] += a12 * dz1 + a22 * dz2;
      }
    }

    // apply bonded forces
    double fspr1 = K * (l1 - L0);
    double fspr2 = K * (l2 - L0);

    //printf("Spring force magnitudes of both bonds %f, %f \n", fspr1,fspr2);

    f[thisI][0] += fspr1 * dx1 / l1 + fspr2 * dx2 / l2;
    f[thisI][1] += fspr1 * dy1 / l1 + fspr2 * dy2 / l2;
    f[thisI][2] += fspr1 * dz1 / l1 + fspr2 * dz2 / l2;

    f[prevI][0] -= fspr1 * dx1 / l1;
    f[prevI][1] -= fspr1 * dy1 / l1;
    f[prevI][2] -= fspr1 * dz1 / l1;

    f[nextI][0] -= fspr2 * dx2 / l2;
    f[nextI][1] -= fspr2 * dy2 / l2;
    f[nextI][2] -= fspr2 * dz2 / l2;

    // apply area elasticity （only works for 2d)
   
    if (KA != 0) {
      //printf("Area elasticity force magnitude %f \n", fA);   
      domain->unmap(atom->x[nextI],atom->image[nextI],unwrap_next);
      domain->unmap(atom->x[prevI],atom->image[prevI],unwrap_prev); 
      f[thisI][0] -= fA * (unwrap_next[1] - unwrap_prev[1]) * 0.5;
      f[thisI][1] -= fA * (unwrap_prev[0] - unwrap_next[0]) * 0.5;
      f[thisI][2] += 0;
    }
  }
}
