#include "atom.h"
#include "atom_masks.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "random_mars.h"
#include "update.h"

#include "memory.h"
#include "modify.h"

#include "domain.h"
#include "fix_membrane2d.h"
#include "force.h"
#include "lammps.h"
#include "math_extra.h"
#include "mpi.h"
#include "pointers.h"
#include "universe.h"
#include "update.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

static constexpr double SMALL = 0.001;

FixMembrane2d::FixMembrane2d(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  K = utils::numeric(FLERR, arg[3], false, lmp);    // edge length elasticity
  L0 = utils::numeric(FLERR, arg[4], false, lmp);
  KC = utils::numeric(FLERR, arg[5], false, lmp);
  theta0 = utils::numeric(FLERR, arg[6], false, lmp) * M_PI / 180;
  KA = utils::numeric(FLERR, arg[7], false, lmp);    // bending rigidity
  A0 = utils::numeric(FLERR, arg[8], false, lmp);    // prefered area
  if (narg == 10) nevery = utils::numeric(FLERR, arg[9], false, lmp);
  // line tension
  tension0 = 0;
  tension_grad0 = 0;
  tension_flag = 0;
  // adaptation parameter of depolymerization rate
  alpha = 20;
  if (narg > 10) {

    double px0 = utils::numeric(FLERR, arg[9], false, lmp);
    double py0 = utils::numeric(FLERR, arg[10], false, lmp);
    double pz0 = utils::numeric(FLERR, arg[11], false, lmp);
    // normalize polarity vector;
    double lp = sqrt(pow(px0, 2) + pow(py0, 2) + pow(pz0, 2));
    px = px0 / lp;
    py = py0 / lp;
    pz = pz0 / lp;
    rate = utils::numeric(FLERR, arg[12], false, lmp);    // add-atom rate
    seed = utils::numeric(FLERR, arg[13], false, lmp);    // random seed
    nevery = utils::numeric(FLERR, arg[14], false, lmp);
    if (rate != 0) isPolymerisationEnabled = true;
    rng = new RanMars(lmp, seed);
    force_reneighbor = 1;

    if (narg > 15) {
      tension0 = utils::numeric(FLERR, arg[15], false, lmp);
      tension_grad0 = utils::numeric(FLERR, arg[16], false, lmp);
      tension_flag = 1;
      if (narg > 17) {
         alpha = utils::numeric(FLERR, arg[17], false, lmp);
         tension_dynamic_flag = 1;
      }
    }
  }

  dim = lmp->domain->dimension;
  nprocs = lmp->universe->nprocs;

  xc = new double[3];
  next_reneighbor = update->ntimestep + 1;
  // to guaratee the depolymerization rate controlled by membrane tension is limited.
 
  rate_lo = 0.0;
  rate_hi = 2.0 * rate;    //  - log(0.01) / (nevery * update->dt);

  if (comm->me == 0 && isPolymerisationEnabled == true) {
    if (1 - exp(-rate * (update->dt) * nevery) > 0.99)
      printf("Warning: probability of polymerization is too high within given update frequency, it "
             "may cause numerical instability \n");
  }

  if (force->newton_bond != 1) force->newton_bond = 1;
  newton_bond = force->newton_bond;

  if (tension_flag) { Etension0 = 2 * sqrt(M_PI * A0) * tension0;}
  
  if(tension_dynamic_flag) {lref = 2 * sqrt(A0 / M_PI);}

  if (atom->membraneBond == nullptr)
    error->all(FLERR, "membraneBond is not defined, must use myatomic atom style");
}

void FixMembrane2d::getCOM()
{
  double xclocal[3] = {0, 0, 0};
  double **x = atom->x;
  int num = 0;    // total number of membrane atoms
  double unwrap[3];
  // MPI snippet
  int numlocal = 0;
  for (int i = 0; i < atom->nlocal; ++i) {
    if (atom->mask[i] & groupbit) {
      domain->unmap(x[i], atom->image[i], unwrap);
      xclocal[0] += unwrap[0];
      xclocal[1] += unwrap[1];
      xclocal[2] += unwrap[2];
      numlocal++;
    }
  }

  MPI_Allreduce(&numlocal, &num, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(xclocal, xc, 3, MPI_DOUBLE, MPI_SUM, world);

  xc[0] /= num;
  xc[1] /= num;
  xc[2] /= num;
  xc[2] = 0;    // force 2d
  domain->remap(xc);
}

void FixMembrane2d::findMaxEdgeLengthAndTag()
{
  double maxEdgeLength = -1;
  maxLengthTag1 = -1;
  maxLengthTag2 = -1;
  double **x = atom->x;
  double unwrap[3];
  // for (int i = 0; i < bondedTags.size(); i++) {
  //   int thisI = atom->map(bondedTags[i][0]);
  //   int nextI = atom->map(bondedTags[i][2]);
  //   if (thisI != -1 && nextI != -1 && thisI < atom->nlocal) {
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      int thisI = i;
      int nextI = atom->map(atom->membraneBond[i][1]);
      if (nextI != -1) {
        double dx = x[nextI][0] - x[thisI][0];
        double dy = x[nextI][1] - x[thisI][1];
        // double dx = x[thisI][0] - xc[0];
        // double dy = x[thisI][1] - xc[1];
        double dz = 0;    //x[nextI][2] - x[thisI][2];
        domain->minimum_image(dx, dy, dz);
        double edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
        // double edgeLength = dx * px + dy*py + dz*pz;
        if (edgeLength > maxEdgeLength) {
          maxEdgeLength = edgeLength;
          // maxLengthTag1 = bondedTags[i][0];
          // maxLengthTag2 = bondedTags[i][2];
          maxLengthTag1 = atom->tag[thisI];
          maxLengthTag2 = atom->membraneBond[thisI][1];
        }
      }
    }
  }
  // handle information from other processors
  struct {
    double val;
    int rank;
  } in, out;

  in.val = maxEdgeLength;
  in.rank = comm->me;
  // printf("I am proc %d, max edge length: %f, tags: %d, %d \n", comm->me, maxEdgeLength, maxLengthTag1, maxLengthTag2);

  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, world);
  MPI_Bcast(&maxLengthTag1, 1, MPI_LMP_TAGINT, out.rank, world);
  MPI_Bcast(&maxLengthTag2, 1, MPI_LMP_TAGINT, out.rank, world);
}

void FixMembrane2d::findMinEdgeLengthAndTag()
{
  double minEdgeLength = 9999;
  minLengthTag1 = -1;
  minLengthTag2 = -1;
  double **x = atom->x;
  double unwrap[3];
  // for (int i = 0; i < bondedTags.size(); i++) {
  //   int thisI = atom->map(bondedTags[i][0]);
  //   int nextI = atom->map(bondedTags[i][2]);
  //   if (thisI != -1 && nextI != -1 && thisI < atom->nlocal) {
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      int thisI = i;
      int nextI = atom->map(atom->membraneBond[i][1]);
      if (nextI != -1) {
        double dx = x[nextI][0] - x[thisI][0];
        double dy = x[nextI][1] - x[thisI][1];
        double dz = 0;    //x[nextI][2] - x[thisI][2];
        domain->minimum_image(dx, dy, dz);
        double edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
        if (edgeLength < minEdgeLength) {
          minEdgeLength = edgeLength;
          minLengthTag1 = atom->tag[thisI];
          minLengthTag2 = atom->membraneBond[thisI][1];
        }
      }
    }
  }
  // handle information from other processors
  struct {
    double val;
    int rank;
  } in, out;

  in.val = minEdgeLength;
  in.rank = comm->me;

  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, world);

  MPI_Bcast(&minLengthTag1, 1, MPI_LMP_TAGINT, out.rank, world);
  MPI_Bcast(&minLengthTag2, 1, MPI_LMP_TAGINT, out.rank, world);

  // if(comm->me == 0)
  //  cout<< "Min proj tag " << minLengthTag1 << "proj " << out.val <<endl;
}

void FixMembrane2d::applyMembraneDynamics(tagint tag1, tagint tag2, int mode_flag)
{
  // mode_flag == 1; add an atom;
  if (mode_flag == 1) {
    double xnew[3], vnew[3], fnew[3];
    if (tag1 != -1 || tag2 != -1) {

      int i1 = atom->map(tag1);
      int i2 = atom->map(tag2);
      int numlocal = atom->nlocal;
      

      tagint maxtag = 0;
      for (int i = 0; i < numlocal; i++) maxtag = MAX(maxtag, atom->tag[i]);
      tagint maxtag_all;
      MPI_Allreduce(&maxtag, &maxtag_all, 1, MPI_LMP_TAGINT, MPI_MAX, world);

      tagint newtag = maxtag_all + 1;
  
      if (i1 != -1 && i2 != -1 && i2 < numlocal) {
        double dx = atom->x[i1][0] - atom->x[i2][0];
        double dy = atom->x[i1][1] - atom->x[i2][1];
        double dz = 0;    //atom->x[i1][2] - atom->x[i2][2];
        domain->minimum_image(dx, dy, dz);
        xnew[0] = 0.5 * dx + atom->x[i2][0];
        xnew[1] = 0.5 * dy + atom->x[i2][1];
        xnew[2] = 0;
        vnew[0] = 0.5 * (atom->v[i1][0] + atom->v[i2][0]);
        vnew[1] = 0.5 * (atom->v[i1][1] + atom->v[i2][1]);
        vnew[2] = 0;
        fnew[0] = 0.5 * (atom->f[i1][0] + atom->f[i2][0]);
        fnew[1] = 0.5 * (atom->f[i1][1] + atom->f[i2][1]);
        fnew[2] = 0;

        imageint newImg = atom->image[i2];

        if (!domain->inside(xnew)) domain->remap(xnew, newImg);
        atom->avec->create_atom(atom->type[i2], xnew);
        atom->avec->copy_bonus(i2, atom->nlocal - 1, 1);
        atom->mask[atom->nlocal - 1] = atom->mask[i2];
        atom->image[atom->nlocal - 1] = newImg;
        if (atom->radius != nullptr) atom->radius[atom->nlocal - 1] = atom->radius[i2];
        atom->v[atom->nlocal - 1][0] = vnew[0];
        atom->v[atom->nlocal - 1][1] = vnew[1];
        atom->v[atom->nlocal - 1][2] = vnew[2];
        atom->f[atom->nlocal - 1][0] = fnew[0];
        atom->f[atom->nlocal - 1][1] = fnew[1];
        atom->f[atom->nlocal - 1][2] = fnew[2];
        atom->tag[atom->nlocal - 1] = newtag;
        atom->membraneBond[atom->nlocal - 1][0] = tag1;
        atom->membraneBond[atom->nlocal - 1][1] = tag2;
        atom->membraneBond[i2][0] = newtag;
      }

      if (i1 != -1 && i1 < numlocal) atom->membraneBond[i1][1] = newtag;
    }
  }
  // mode_flag = 0; delete an atom
  if (mode_flag == 0) {
    int *dlist = new int[atom->nlocal];
    tagint tag1_prev = -1;
    for (int i = 0; i < atom->nlocal; ++i) {
      if (atom->tag[i] == tag1) {
        dlist[i] = 1;
        tag1_prev = atom->membraneBond[i][0];
      } else {
        dlist[i] = 0;
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &tag1_prev, 1, MPI_LMP_TAGINT, MPI_MAX, world);

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

    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->tag[i] == tag2) { atom->membraneBond[i][0] = tag1_prev; }
      if (atom->tag[i] == tag1_prev) { atom->membraneBond[i][1] = tag2; }
    }
  }
}


void FixMembrane2d::initializeMembraneBonds2D()
{
  // initialize the bondedTags counter clockwise
  if(membraneBond_init_flag == 1) return;

  //
  vector<tagint> tags_local;
  vector<double> angs_local;
  int m = 0;

  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      double dx = atom->x[i][0] - xc[0];
      double dy = atom->x[i][1] - xc[1];
      double dz = 0;    //atom->x[i][2] - xc[2];
      domain->minimum_image(dx, dy, dz);

      double ang = atan2(dy, dx);
      tags_local.push_back(atom->tag[i]);
      angs_local.push_back(ang);
      m++;
    }
  }
  vector<int> data_sizes(nprocs, 0);
  vector<int> displacements(nprocs, 0);
  MPI_Allgather(&m, 1, MPI_INT, data_sizes.data(), 1, MPI_INT, world);
  int totalSize = accumulate(data_sizes.begin(), data_sizes.end(), 0);
  partial_sum(data_sizes.begin(), data_sizes.end() - 1, displacements.begin() + 1);

  vector<tagint> tags_total(totalSize);
  vector<double> angs_total(totalSize);
  MPI_Allgatherv(tags_local.data(), m, MPI_LMP_TAGINT, tags_total.data(), data_sizes.data(),
                 displacements.data(), MPI_LMP_TAGINT, world);
  MPI_Allgatherv(angs_local.data(), m, MPI_DOUBLE, angs_total.data(), data_sizes.data(),
                 displacements.data(), MPI_DOUBLE, world);

  vector<pair<tagint, double>> temp;
  for (int i = 0; i < totalSize; i++) { temp.push_back({tags_total[i], angs_total[i]}); }

  sort(temp.begin(), temp.end(), [](const pair<tagint, double> &a, const pair<tagint, double> &b) {
    return a.second < b.second;
  });

  for (int i = 0; i < temp.size(); i++) {
   // bondedTags.push_back({temp[i].first, temp[(i - 1 + temp.size()) % temp.size()].first,
   //                       temp[(i + 1) % temp.size()].first});

    tagint tagi = temp[i].first;
    tagint tagj = temp[(i - 1 + temp.size()) % temp.size()].first;
    tagint tagk = temp[(i + 1) % temp.size()].first;
    int thisI = atom->map(tagi);
    if(thisI != -1 && thisI < atom->nlocal)
    {
      atom->membraneBond[thisI][0] = tagj;
      atom->membraneBond[thisI][1] = tagk;
    }
  }
  membraneBond_init_flag = 1;

  // if (!bondedTags.empty()) {
  //   for (int i = 0; i < bondedTags.size(); i++) {
  //     int thisI = atom->map(bondedTags[i][0]);
  //     if (thisI != -1 && thisI < atom->nlocal) {
  //       atom->membraneBond[thisI][0] = bondedTags[i][1];
  //       atom->membraneBond[thisI][1] = bondedTags[i][2];
  //     }
  //   }
  // }
}

int FixMembrane2d::setmask()
{
  int mask;
  mask |= PRE_EXCHANGE;    // add and delete atoms here
  mask |= PRE_FORCE;       // calculate force here
  mask |= POST_FORCE;      // apply force here
  mask |= END_OF_STEP;
  return mask;
}

void FixMembrane2d::setup(int)
{
  atom->map_init(1);
  atom->map_set();
  getCOM();
  if (dim == 2 && membraneBond_init_flag == 0) initializeMembraneBonds2D();
  if (tension_flag != 0) updatePredefinedTension();
  init_update = update->ntimestep + 1;
}

void FixMembrane2d::updatePredefinedTension()
{
  int size = atom->nlocal + atom->nghost;
  tension_coeff.clear();
  tension_coeff.resize(size);

  double **x = atom->x;
  double max_proj_local = -1;
  int max_proj_local_index = -1;
  double min_proj_local = 9999;
  int min_proj_local_index = 9999;
  vector<pair<int, double>> proj;
  // printf("I am procs %d, xc: %f, %f, %f \n", comm->me, xc[0], xc[1], xc[2]);
  // for (int i = 0; i < bondedTags.size(); i++) {
  //     int thisI = atom->map(bondedTags[i][0]);
  for (int i = 0; i < atom->nlocal + atom->nghost; i++) {
    if (atom->mask[i] & groupbit) {
      int thisI = i;
      // if(thisI != -1)
      // {
      double dx = x[thisI][0] - xc[0];
      double dy = x[thisI][1] - xc[1];
      double dz = 0;    //x[thisI][2] - xc[2];
      domain->minimum_image(dx, dy, dz);
      double projtemp = dx * px + dy * py + dz * pz;
      proj.push_back({thisI, projtemp});
      if (projtemp > max_proj_local) {
        max_proj_local = projtemp;
        max_proj_local_index = thisI;
      }
      if (projtemp < min_proj_local) {
        min_proj_local = projtemp;
        min_proj_local_index = thisI;
      }
    }
  }

  struct {
    double val;
    int rank;
  } in_max, out_max, in_min, out_min;

  in_min.val = min_proj_local;
  in_min.rank = comm->me;
  in_max.val = max_proj_local;
  in_max.rank = comm->me;
  MPI_Allreduce(&in_min, &out_min, 1, MPI_DOUBLE_INT, MPI_MINLOC, world);
  MPI_Allreduce(&in_max, &out_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, world);
  double lproj = out_max.val - out_min.val;
  double scale = 1;   // scale the tension gradient dynamically based on the reference length of the membrane
  if(tension_dynamic_flag) scale = lproj / lref; // test with the ratio
  for (int i = 0; i < proj.size(); i++) {
    double tension =
        tension0 + tension_grad0 * scale * (1 + cos(M_PI * (proj[i].second - out_min.val) / lproj));
    tension_coeff[proj[i].first] = tension;
  }
}

void FixMembrane2d::computeTensionEnergy()
{
  double **x = atom->x;
  double Etension_local = 0;
  // for (int i = 0; i < bondedTags.size(); i++) {
  //   int thisI = atom->map(bondedTags[i][0]);
  //   int nextI = atom->map(bondedTags[i][2]);
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      int thisI = i;
      int nextI = atom->map(atom->membraneBond[i][1]);
      // if (thisI != -1 && nextI != -1 && thisI < atom->nlocal) {
      if (nextI != -1) {
        double dx = x[nextI][0] - x[thisI][0];
        double dy = x[nextI][1] - x[thisI][1];
        double dz = 0;    //x[nextI][2] - x[thisI][2];
        domain->minimum_image(dx, dy, dz);
        double edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
        Etension_local += tension_coeff[thisI] * edgeLength;
      } else
        error->one(FLERR, "Error: bonded atom not found");
    }
  }
  MPI_Allreduce(&Etension_local, &Etension, 1, MPI_DOUBLE, MPI_SUM, world);
}

void FixMembrane2d::pre_exchange()
{
  isChanged = 0;
  if (!isPolymerisationEnabled) return;                // no polymerization
  if (next_reneighbor != update->ntimestep) return;    // not time to polymerize or reneighbor

  double prob0 = 1 - exp(-rate * (update->dt) * nevery);    // probability of polymerization

  computeTensionEnergy();
  // if(next_reneighbor == init_update)  Etension0 = Etension;
  // if (comm->me == 0) { printf("fraction: %f \n", Etension / Etension0); }

  double rate_depoly = rate_lo + (rate_hi - rate_lo) / (1 + exp(-alpha * (Etension / Etension0 - 1)));
  // if(comm->me == 0)
  // {
  //   printf("de rate %f, fraction %f, \n", rate_depoly, Etension/Etension0);
  // }
  double prob1 =
      1 - exp(-rate_depoly * (update->dt) * nevery);    // probability of depolymerization

  bool isPolymerized = rng->uniform() <= prob0 ? true : false;
  bool isDepolymerized = rng->uniform() <= prob1 ? true : false;

  // new tags for new atoms;

  // normal polymerization

  if (isPolymerized) {
    // all polyerization mechanisms go here
    findMaxEdgeLengthAndTag();
    applyMembraneDynamics(maxLengthTag1, maxLengthTag2, 1);


    isChanged += 1;

    atom->map_init(1);
    atom->map_set();
  }

  if (isDepolymerized) {
    findMinEdgeLengthAndTag();
    applyMembraneDynamics(
        minLengthTag1, minLengthTag2,
        0);    // if deleting an atom, the minLengthTag1 is the tag of the atom to be deleted
    isChanged += 1;
    atom->map_init(1);
    atom->map_set();
  }

  if (isChanged != 0) {
    getCOM();
  }

  next_reneighbor += nevery;
  maxLengthTag1 = -1;
  maxLengthTag2 = -1;
  minLengthTag1 = -1;
  minLengthTag2 = -1;
}

void FixMembrane2d::pre_force(int)
{
  // incrementing age
  if (atom->age != nullptr) {
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->mask[i] & groupbit) { atom->age[i]++; }
    }
  }
  // for gather error info
  // int breakflag = 0;
  // tagint errtag1 = -1, errtag2 = -1;

  if (next_reneighbor != update->ntimestep && tension_flag != 0) updatePredefinedTension();

  double **f = atom->f;
  /* Prepare for Area elasticity */
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  if (fOnProc.size() != nlocal + nghost) {
    fOnProc.resize(nlocal + nghost);
    for (int i = 0; i < nlocal + nghost; ++i) { fOnProc[i].resize(3); }
  }
  for (int i = 0; i < nlocal + nghost; ++i) {
    fOnProc[i][0] = 0;
    fOnProc[i][1] = 0;
    fOnProc[i][2] = 0;
  }
  if (KA != 0) {

    //getCOM();
    double currentArea = 0.0;
    // for (int i = 0; i < bondedTags.size(); i++) {
    //   int thisI = atom->map(bondedTags[i][0]);
    //   int nextI = atom->map(bondedTags[i][2]);
    for (int i = 0; i < nlocal; i++) {
      if (atom->mask[i] & groupbit) {
        int thisI = i;
        int nextI = atom->map(atom->membraneBond[i][1]);
        // if (thisI != -1 && thisI < nlocal) {
        double x1, x2, x3, y1, y2, y3, z1, z2;
        x1 = atom->x[thisI][0] - xc[0];
        y1 = atom->x[thisI][1] - xc[1];
        z1 = 0;    //atom->x[thisI][2] - xc[2];

        domain->minimum_image(x1, y1, z1);
        if (nextI == -1) {

          printf("error in area computation due to nextI is not found in the atom list on current "
                 "processor, \n");
          printf("this happens when polymerization status is changed: %d \n", isChanged);

          if (minLengthTag1 != -1)
            printf("Potentially removed tags: %d, %d \n", minLengthTag1, minLengthTag2);

          if (maxLengthTag1 != -1)
            printf("tags related of adding particles: %d, %d \n", maxLengthTag1, maxLengthTag2);

          printf("current step %ld, this tag %d, next tag %d \n", update->ntimestep, atom->tag[i],
                 atom->membraneBond[i][1]);
          printf("processor %d, num local atoms %d, num ghost atoms %d \n", comm->me, nlocal,
                 nghost);
          printf("thisI: %d, nextI: %d \n", thisI, nextI);
          printf("xc %f, %f, %f \n", xc[0], xc[1], xc[2]);
          // breakflag = 1;
          // errtag1 = bondedTags[i][0];
          // errtag2 = bondedTags[i][2];
          // break;
        }

        x2 = atom->x[nextI][0] - xc[0];
        y2 = atom->x[nextI][1] - xc[1];
        z2 = 0;    //atom->x[nextI][2] - xc[2];

        domain->minimum_image(x2, y2, z2);

        x3 = 0;
        y3 = 0;
        currentArea += abs(0.5 * (x1 * y2 + x2 * y3 + x3 * y1 - x1 * y3 - x2 * y1 - x3 * y2));
      }
      // }
    }
    MPI_Allreduce(MPI_IN_PLACE, &currentArea, 1, MPI_DOUBLE, MPI_SUM, world);
    double dA = currentArea - A0;
    double fA = KA * dA;
    for (int i = 0; i < nlocal; i++) {
      if(atom->mask[i] & groupbit)
      {
        int thisI = i;
        int prevI = atom->map(atom->membraneBond[i][0]);
        int nextI = atom->map(atom->membraneBond[i][1]);
        if (prevI != -1 && nextI != -1) {
          double dx = atom->x[nextI][0] - atom->x[prevI][0];
          double dy = atom->x[nextI][1] - atom->x[prevI][1];
          double dz = 0;    //atom->x[nextI][2] - atom->x[prevI][2];
          domain->minimum_image(dx, dy, dz);
          double fx = fA * dy * 0.5;
          double fy = -fA * dx * 0.5;
          f[thisI][0] -= fx;
          f[thisI][1] -= fy;
          //printf("Atom tag: %d, area force: %f, %f \n", atom->tag[thisI], fx, fy);
        }
      }
    }
  }
  // MPI_Allreduce(MPI_IN_PLACE, &breakflag, 1, MPI_INT, MPI_MAX, world);
  // if(breakflag > 0 )
  // {
  //   MPI_Allreduce(MPI_IN_PLACE, &errtag1, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  //   MPI_Allreduce(MPI_IN_PLACE, &errtag2, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  //   printf("Procs %d owns %d atoms, map of invovled tags %d, %d and %d, %d \n", comm->me, atom->nlocal ,errtag1, atom->map(errtag1), errtag2,atom->map(errtag2));
  //   printf("Newly created atom after add on proc %d is at %f %f \n", comm->me, xnew_after_add[0],xnew_after_add[1]);
  //   printf("Newly created atom after rem on proc %d is at %f %f \n", comm->me, xnew_after_remove[0],xnew_after_remove[1]);
  //   error->all(FLERR, "Error in area computation due to nextI is not found in the atom list on current processor");
  // }

  // apply bending rigidity
  // for (int i = 0; i < bondedTags.size(); i++) {
  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      double dx1, dx2, dy1, dy2, dz1, dz2, theta, l1, l2, l1sq, l2sq, stheta, ctheta;
      //  int thisI = atom->map(bondedTags[i][0]);
      //  int prevI = atom->map(bondedTags[i][1]);
      //  int nextI = atom->map(bondedTags[i][2]);
      int thisI = i;
      int prevI = atom->map(atom->membraneBond[i][0]);
      int nextI = atom->map(atom->membraneBond[i][1]);
      //  if (thisI != -1 && prevI != -1 && nextI != -1 && thisI < nlocal) {
      if (prevI != -1 && nextI != -1) {
        dx1 = atom->x[prevI][0] - atom->x[thisI][0];
        dy1 = atom->x[prevI][1] - atom->x[thisI][1];
        dz1 = 0;    //atom->x[prevI][2] - atom->x[thisI][2];
        domain->minimum_image(dx1, dy1, dz1);
        dx2 = atom->x[nextI][0] - atom->x[thisI][0];
        dy2 = atom->x[nextI][1] - atom->x[thisI][1];
        dz2 = 0;    //atom->x[nextI][2] - atom->x[thisI][2];
        domain->minimum_image(dx2, dy2, dz2);
        l1sq = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
        l2sq = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;
        l1 = sqrt(l1sq);
        l2 = sqrt(l2sq);
        if (l1sq < SMALL) l1sq = SMALL;
        if (l2sq < SMALL) l2sq = SMALL;
        if (l1 < SMALL) l1 = SMALL;
        if (l2 < SMALL) l2 = SMALL;
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
          double f_this[3] = {0, 0, 0};
          double f_prev[3] = {0, 0, 0};
          double f_next[3] = {0, 0, 0};

          f_this[0] = (a11 + a12) * dx1 + (a12 + a22) * dx2;
          f_this[1] = (a11 + a12) * dy1 + (a12 + a22) * dy2;
          f_this[2] = (a11 + a12) * dz1 + (a12 + a22) * dz2;

          f_prev[0] = a11 * dx1 + a12 * dx2;
          f_prev[1] = a11 * dy1 + a12 * dy2;
          f_prev[2] = a11 * dz1 + a12 * dz2;

          f_next[0] = a12 * dx1 + a22 * dx2;
          f_next[1] = a12 * dy1 + a22 * dy2;
          f_next[2] = a12 * dz1 + a22 * dz2;
          //printf("%d, %f, %f \n", bondedTags[i][0], f_this[0], f_this[1]);
          if (newton_bond || thisI < nlocal) {
            fOnProc[thisI][0] -= f_this[0];
            fOnProc[thisI][1] -= f_this[1];
            fOnProc[thisI][2] -= f_this[2];
          }
          if (newton_bond || prevI < nlocal) {
            fOnProc[prevI][0] += f_prev[0];
            fOnProc[prevI][1] += f_prev[1];
            fOnProc[prevI][2] += f_prev[2];
          }
          if (newton_bond || nextI < nlocal) {
            fOnProc[nextI][0] += f_next[0];
            fOnProc[nextI][1] += f_next[1];
            fOnProc[nextI][2] += f_next[2];
          }
        }

        // calculate stretching force
        if (K != 0) {
          double fspr1 = K * (l1 - L0);
          double fspr2 = K * (l2 - L0);
          double f_this[3] = {0, 0, 0};
          double f_prev[3] = {0, 0, 0};
          double f_next[3] = {0, 0, 0};
          f_this[0] = fspr1 * dx1 / l1 + fspr2 * dx2 / l2;
          f_this[1] = fspr1 * dy1 / l1 + fspr2 * dy2 / l2;
          f_this[2] = fspr1 * dz1 / l1 + fspr2 * dz2 / l2;

          f_prev[0] = fspr1 * dx1 / l1;
          f_prev[1] = fspr1 * dy1 / l1;
          f_prev[2] = fspr1 * dz1 / l1;

          f_next[0] = fspr2 * dx2 / l2;
          f_next[1] = fspr2 * dy2 / l2;
          f_next[2] = fspr2 * dz2 / l2;

          if (newton_bond || thisI < nlocal) {
            fOnProc[thisI][0] += f_this[0];
            fOnProc[thisI][1] += f_this[1];
            fOnProc[thisI][2] += f_this[2];
          }
          if (newton_bond || prevI < nlocal) {
            fOnProc[prevI][0] -= f_prev[0];
            fOnProc[prevI][1] -= f_prev[1];
            fOnProc[prevI][2] -= f_prev[2];
          }
          if (newton_bond || nextI < nlocal) {
            fOnProc[nextI][0] -= f_next[0];
            fOnProc[nextI][1] -= f_next[1];
            fOnProc[nextI][2] -= f_next[2];
          }
        }

        // apply line tension
        if (tension_flag != 0) {
          double tension = tension_coeff[thisI];
          double f_this[3] = {0, 0, 0};
          double f_prev[3] = {0, 0, 0};
          double f_next[3] = {0, 0, 0};

          f_this[0] = tension_coeff[thisI] * dx2 / l2 + tension_coeff[prevI] * dx1 / l1;
          f_this[1] = tension_coeff[thisI] * dy2 / l2 + tension_coeff[prevI] * dy1 / l1;
          f_this[2] = tension_coeff[thisI] * dz2 / l2 + tension_coeff[prevI] * dz1 / l1;

          f_prev[0] = tension_coeff[prevI] * dx1 / l1;
          f_prev[1] = tension_coeff[prevI] * dy1 / l1;
          f_prev[2] = tension_coeff[prevI] * dz1 / l1;

          f_next[0] = tension_coeff[thisI] * dx2 / l2;
          f_next[1] = tension_coeff[thisI] * dy2 / l2;
          f_next[2] = tension_coeff[thisI] * dz2 / l2;

          if (newton_bond || thisI < nlocal) {
            fOnProc[thisI][0] += f_this[0];
            fOnProc[thisI][1] += f_this[1];
            fOnProc[thisI][2] += f_this[2];
          }

          if (newton_bond || prevI < nlocal) {
            fOnProc[prevI][0] -= f_prev[0];
            fOnProc[prevI][1] -= f_prev[1];
            fOnProc[prevI][2] -= f_prev[2];
          }

          if (newton_bond || nextI < nlocal) {
            fOnProc[nextI][0] -= f_next[0];
            fOnProc[nextI][1] -= f_next[1];
            fOnProc[nextI][2] -= f_next[2];
          }
        }
      }
    }
  }

  comm->reverse_comm(this, 3);
  // clear values
  // xnew_after_add[0] = 0; xnew_after_add[1] = 0; xnew_after_add[2] = 0;
  // xnew_after_remove[0] = 0; xnew_after_remove[1] = 0; xnew_after_remove[2] = 0;
}

void FixMembrane2d::post_force(int)
{
  double **f = atom->f;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      f[i][0] += fOnProc[i][0];
      f[i][1] += fOnProc[i][1];
      f[i][2] += fOnProc[i][2];
    }
  }
}

void FixMembrane2d::end_of_step()
{
  getCOM();    // this is called nevery steps for efficiency since COM does not change drastically,
}

int FixMembrane2d::pack_reverse_comm(int n, int first, double *buff)
{
  int i, m = 0;
  for (i = first; i < first + n; i++) {
    buff[m++] = fOnProc[i][0];
    buff[m++] = fOnProc[i][1];
    buff[m++] = fOnProc[i][2];
  }
  return m;
}

void FixMembrane2d::unpack_reverse_comm(int n, int *list, double *buff)
{
  int i, j, m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    fOnProc[j][0] += buff[m++];
    fOnProc[j][1] += buff[m++];
    fOnProc[j][2] += buff[m++];
  }
}