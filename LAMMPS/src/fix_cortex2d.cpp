#include "atom.h"
#include "atom_masks.h"
#include "atom_vec.h"
#include "error.h"
#include "group.h"
#include "comm.h"
#include "random_mars.h"
#include "update.h"

#include "fix_cortex2d.h"
#include "force.h"
#include "polygon2D.h"
#include "update.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <tuple>
#include <vector>
#include "pointers.h"
#include "domain.h"
#include "lammps.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;


static constexpr double SMALL = 0.001;

FixCortex2d::FixCortex2d(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  K = utils::numeric(FLERR, arg[3], false, lmp); // edge length elasticity
  L0 = utils::numeric(FLERR, arg[4], false, lmp);
  KC = utils::numeric(FLERR, arg[5], false, lmp);
  theta0 = utils::numeric(FLERR, arg[6], false, lmp); 
  px = utils::numeric(FLERR, arg[7], false, lmp);
  py = utils::numeric(FLERR, arg[8], false, lmp);
  pz = utils::numeric(FLERR, arg[9], false, lmp);
  rate = utils::numeric(FLERR, arg[10], false, lmp);    // add-atom rate
  nevery = utils::numeric(FLERR, arg[11], false, lmp);
  
  rng = new RanMars(lmp, 12345); // Initialize the seed with a valid value
  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  dim = lmp->domain->dimension;
}



void FixCortex2d::getCOM()
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
  cz /= num;    // COM of the cortex

  xc[0] = cx;
  xc[1] = cy;
  xc[2] = cz;
}

tagint FixCortex2d::deleteOneAtomAtLargestNegativePolarity(vector<double>& proj, vector<int>& ind)
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
  return delTag;

  // Note that if multiple minima exist, only delete the first one. 
}

vector<tagint> FixCortex2d::addOneAtomAtLargestPositivePolarity(vector<double> &proj, vector<int> &ind)
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
  double irad; // radius or interaction range;
  imageint img;

  tagint tagMax = atom->tag[maxI];

  tagint tagNextMax;
  for(int i = 0; i < bondedTags.size();++i)
  {
    if(tagMax == bondedTags[i][0])
    {
      tagNextMax = bondedTags[(i+1) % bondedTags.size()][0];
    }
  }
  int nextMaxI = atom->map(tagNextMax);

  // get the atom type to create
  itype = atom ->type[nextMaxI];
  imask = atom->mask[nextMaxI];
  if(atom->radius != nullptr)
     irad = atom->radius[nextMaxI];
  else
  {
    irad = 0.5 * L0;
  }

  img = atom->image[nextMaxI];
  xnew[0] = atom->x[nextMaxI][0] + px * 2.05 * irad;
  xnew[1] = atom->x[nextMaxI][1] + py * 2.05 * irad;
  xnew[2] = atom->x[nextMaxI][2] + pz * 2.05 * irad;
 

  if(!domain->inside(xnew)) 
  { 
    domain->remap(xnew,img); 
  } 

  


  atom->avec->create_atom(itype, xnew);
  int n = atom->nlocal - 1;
  atom->tag[n] =
      0;    // by default, in new versions of LAMMPS, newly created atom automatically has a tag of 0
  atom->mask[n] = imask;
  if(atom->radius != nullptr) atom->radius[n] = irad;
  atom->image[n] = img;

  atom->tag_extend();
  atom->tag_check();

  //printf("an atom with a tag of %d is added at position %f, %f,%f \n",atom->tag[n],atom->x[n][0],atom->x[n][1],atom->x[n][2]);
  return {atom->tag[n],tagMax};
}

void FixCortex2d::initializeMembraneBonds2D()
{
  // cx,cy COM position of membrane
    double unwrap[3];
    vector<tuple<double, double, tagint>> vertices;
    for (int i = 0; i < atom->nlocal; ++i)    // only works for serial program
    {
      if (atom->mask[i] & groupbit) {
        domain->unmap(atom->x[i],atom->image[i],unwrap);
        vertices.push_back({unwrap[0], unwrap[1], atom->tag[i]});
      }
    }
    vector<tuple<double, double, tagint>> sorted = polygon2D::sortVertices(vertices);

    if (bondedTags.empty()) {
      for (int i = 0; i < sorted.size(); ++i) {
        int j = (i + 1) % sorted.size();
        int k = (i + sorted.size() - 1) % sorted.size();
        bondedTags.push_back({get<2>(sorted[i]), get<2>(sorted[k]), get<2>(sorted[j])});
      }
    }
}

void FixCortex2d::updateMembraneBonds2D(vector<tagint> addTags, tagint delTag)
{
  tagint newTag = addTags[0];
  tagint posTag = addTags[1];
  vector<tagint> sortedTags; // this variable is originally sorted counterclockwisely -pi to pi
  for(int i = 0; i < bondedTags.size(); ++i)
  {
    if(bondedTags[i][0] != delTag)
    {
      sortedTags.push_back(bondedTags[i][0]);
      if(bondedTags[i][0] == posTag)
      {
        sortedTags.push_back(newTag);
      }
    }
  }

  bondedTags.clear();
  for(int i = 0; i < sortedTags.size(); ++i)
  {
    int j = (i + 1) % sortedTags.size();
    int k = (i + sortedTags.size() - 1) % sortedTags.size();
    bondedTags.push_back({sortedTags[i],sortedTags[k],sortedTags[j]});
  }
  



}




void FixCortex2d::init()
{
  // normalize polarity vector;
  double lp = sqrt(pow(px,2) + pow(py,2) + pow(pz,2));
  px /= lp;
  py /= lp;
  pz /= lp;

  theta0 *= M_PI / 180;

  xc = new double[3];

  seed = 3121545; // need to change 
  // initialize bonded particles
  if (dim == 2) {
    initializeMembraneBonds2D();
  } 
}

int FixCortex2d::setmask()
{
  int mask;
  mask |= PRE_EXCHANGE; // add and delete atoms here
  mask |= POST_FORCE; // apply force here
  mask |= END_OF_STEP;
  return mask;
}


void FixCortex2d::pre_exchange()
{
  isChanged = false;
  if(next_reneighbor != update->ntimestep) return;

  double prob = 1 - exp(-rate * (update->dt) * nevery);
  if(rng->uniform() <= prob)
  {
   
    getCOM();

    vector<int> ind; 
    vector<double> proj; 
    double unwrap[3];
    double **x = atom->x;
    for(int i = 0; i < atom->nlocal; ++i)
    {
      if(atom->mask[i] & groupbit)
      {
        domain ->unmap(x[i],atom->image[i],unwrap);
        double projtmp = px * (unwrap[0] - xc[0]) + py * (unwrap[1] - xc[1]) + pz * (unwrap[2] - xc[2]);
        proj.push_back(projtmp);
        ind.push_back(i);
      }
    } // projection along polarity direction. Positive -> cell-front, Negative->cell-back

    // depolymerization happens at the position with minimal projection along the palority
    // delete one particle at the negative polarity direction 

    tagint delTag = deleteOneAtomAtLargestNegativePolarity(proj,ind);


    // polymerization happens at the position with maximal projection along the polarity
    // add one particle at the center of particles having the largest proj and the second largest proj.
    
    vector<tagint> addTags = addOneAtomAtLargestPositivePolarity(proj,ind);  // addTags[0]: tag of the newly created atom; addTags[1]: the tag of the atom where you want to add an atom


    // update new bonds
    if(dim == 2)
    {
      updateMembraneBonds2D(addTags,delTag);
    }

    isChanged = true;
  }

  if(isChanged)
  {
    atom->map_delete();
    atom->map_init(1);
    atom->map_set();
  }


  next_reneighbor += nevery;
}



void FixCortex2d::post_force(int)
{
  double **f = atom->f;

  // apply bending rigidity
  double dx1,dx2,dy1,dy2,dz1,dz2,theta,l1,l2,l1sq,l2sq,stheta,ctheta;
  for(int i = 0; i < bondedTags.size(); ++i)
  {
        int nextI = atom->map(bondedTags[i][2]);
        int thisI = atom->map(bondedTags[i][0]);
        int prevI = atom->map(bondedTags[i][1]);


        dx1 = atom->x[prevI][0] - atom->x[thisI][0]; 
        dy1 = atom->x[prevI][1] - atom->x[thisI][1];
        dz1 = atom->x[prevI][2] - atom->x[thisI][2];
        domain->minimum_image(dx1,dy1,dz1);

        dx2 = atom->x[nextI][0] - atom->x[thisI][0];
        dy2 = atom->x[nextI][1] - atom->x[thisI][1];
        dz2 = atom->x[nextI][2] - atom->x[thisI][2];
        domain->minimum_image(dx2,dy2,dz2);

        l1sq = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
        l2sq = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;

        l1 = sqrt(l1sq);
        l2 = sqrt(l2sq);

        if (KC != 0) {
          ctheta = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (l1 * l2);    // cos

          if (ctheta > 1) ctheta = 1;
          if (ctheta < -1) ctheta = -1;

          stheta = sqrt(1 - ctheta * ctheta);    // sin  // D[cos] = - 1/sin
          if (stheta < SMALL) stheta = SMALL;
          double dtheta = acos(ctheta) - theta0;

          double fmag = KC * dtheta;    //

          double a = -fmag / stheta;    // prefactor

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
        // apply bonded forces
        double fspr1 = K * (l1 - L0);
        double fspr2 = K * (l2 - L0);
        

        f[thisI][0] += fspr1 * dx1 / l1 + fspr2 * dx2 / l2;
        f[thisI][1] += fspr1 * dy1 / l1 + fspr2 * dy2 / l2;
        f[thisI][2] += fspr1 * dz1 / l1 + fspr2 * dz2 / l2;

        f[prevI][0] -= fspr1 * dx1 / l1;
        f[prevI][1] -= fspr1 * dy1 / l1;
        f[prevI][2] -= fspr1 * dz1 / l1;

        f[nextI][0] -= fspr2 * dx2 / l2;
        f[nextI][1] -= fspr2 * dy2 / l2;
        f[nextI][2] -= fspr2 * dz2 / l2;

  }




}

void FixCortex2d::end_of_step()
{
  getCOM(); // update Center of mass
}