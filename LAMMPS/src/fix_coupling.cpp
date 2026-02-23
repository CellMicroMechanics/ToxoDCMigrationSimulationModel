#include "atom.h"
#include "atom_masks.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "random_mars.h"
#include "update.h"

#include "domain.h"
#include "fix_coupling.h"
#include "force.h"
#include "lammps.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pointers.h"
#include "update.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>

using namespace std;
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace NeighConst;

static constexpr double SMALL = 0.001;

FixCoupling::FixCoupling(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  jgroup = group->find(arg[3]);
  jgroupbit =
      group->bitmask
          [jgroup];    // by default jgroup is the group of atoms that are needed to check ages
  if (jgroup < 0) error->all(FLERR, "Invalid group ID");

  style = arg[4];
  if (strcmp(style, "harmonic") == 0) {

    K = utils::numeric(FLERR, arg[5], false, lmp);
    L0 = utils::numeric(FLERR, arg[6], false, lmp);
    Lcut = utils::numeric(FLERR, arg[7], false, lmp);
    if (narg >= 9) {
      nevery = utils::inumeric(FLERR, arg[8], false, lmp);
      check_age_flag = utils::inumeric(FLERR, arg[9], false, lmp);
      if (check_age_flag != 0) { 
        age_min = utils::inumeric(FLERR, arg[10], false, lmp); 
        age_max = (narg == 12) ? utils::inumeric(FLERR, arg[11], false, lmp) : numeric_limits<int>::max();
      }
    } else {
      nevery = 1;
    }

  } else {
    error->all(FLERR, "Invalid coupling style between desired groups");
  }

  if (K < 0) error->all(FLERR, "Invalid K value, spring constant must be positive");
  if (L0 > Lcut) error->all(FLERR, "Invalid L0 value, rest length must be less than cutoff length");
  if (check_age_flag != 0 && age_min <= nevery)
    error->all(FLERR,
               "Invalid age_min should be larger than update frequency of coupling bonds");
  L0sq = L0 * L0;
  Lcutsq = Lcut * Lcut;
  nProcs = comm->nprocs;
  init_update = update->ntimestep + 1;
  next_update = update->ntimestep + nevery;
}

int FixCoupling::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}
void FixCoupling::init()
{
  NeighRequest *req = neighbor->add_request(this, REQ_OCCASIONAL | REQ_FULL);
  double temp_cut = Lcut;

  if (comm->get_comm_cutoff() <= temp_cut) {
    error->all(FLERR,
               "Communication cutoff is smaller than the cutoff of coupling cutoff distance.");
  }

  req->set_cutoff(temp_cut);
}

void FixCoupling::init_list(int, NeighList *ptr)
{
  list = ptr;
}

void FixCoupling::setup(int)
{
  double cutoffmax = neighbor->cutneighmax;
  double skin = neighbor->skin;
  cutoff_plus_skin = cutoffmax + skin;
  if (cutoff_plus_skin < L0 && comm->me == 0)
    error->warning(FLERR, "L0 is larger than the distance of building neighbor list, \
   occasional neighbor lists will be built.");
  double cut_comm = comm->get_comm_cutoff();
  if (Lcut > cut_comm) error->all(FLERR, "Lcut is larger than the communication cutoff.");
}

void FixCoupling::pre_force(int)
{
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
  if (init_update == update->ntimestep) {
    if (check_age_flag == 0)
      initialize_bonds();
    else
      initialize_bonds_check_age();
    next_update += nevery;
  } else if (next_update == update->ntimestep) {
    if (check_age_flag == 0)
      update_bonds();
    else
      update_bonds_check_age();
    next_update += nevery;
    // if(!coupling_list.empty() && comm->me == 0)
    //   printf("Coupling list is activate.\n");
  }

  if (K == 0) return;
  for (int i = 0; i < coupling_list.size(); i++) {
    int i1 = atom->map(coupling_list[i][0]);
    if (i1 != -1 && i1 < nlocal) {
      for (int j = 1; j < coupling_list[i].size(); j++) {
        int i2 = atom->map(coupling_list[i][j]);
        if (i2 != -1) {
          double dx = atom->x[i2][0] - atom->x[i1][0];
          double dy = atom->x[i2][1] - atom->x[i1][1];
          double dz = 0;    //atom->x[i2][2] - atom->x[i1][2];
          domain->minimum_image(dx, dy, dz);
          double rsq = dx * dx + dy * dy + dz * dz;

          double r = sqrt(rsq);
          double dr = (r - L0) > 0 ? r - L0 : 0;    // cable bonds are only tensile
          // double dr = r - L0;
          double f = K * dr / r;
          if (force->newton_bond || i1 < nlocal) {
            fOnProc[i1][0] += f * dx;
            fOnProc[i1][1] += f * dy;
            fOnProc[i1][2] += f * dz;
          }

          if (force->newton_bond || i2 < nlocal) {
            fOnProc[i2][0] -= f * dx;
            fOnProc[i2][1] -= f * dy;
            fOnProc[i2][2] -= f * dz;
          }
        }
      }
    }
  }
}

void FixCoupling::post_force(int)
{
  int nlocal = atom->nlocal;
  double **f = atom->f;
  for (int i = 0; i < nlocal; ++i) {
    f[i][0] += fOnProc[i][0];
    f[i][1] += fOnProc[i][1];
    f[i][2] += fOnProc[i][2];
  }
}

int FixCoupling::pack_reverse_comm(int n, int first, double *buff)
{
  int i, m = 0;
  for (i = first; i < first + n; i++) {
    buff[m++] = fOnProc[i][0];
    buff[m++] = fOnProc[i][1];
    buff[m++] = fOnProc[i][2];
  }
  return m;
}

void FixCoupling::unpack_reverse_comm(int n, int *list, double *buff)
{
  int i, j, m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    fOnProc[j][0] += buff[m++];
    fOnProc[j][1] += buff[m++];
    fOnProc[j][2] += buff[m++];
  }
}

void FixCoupling::initialize_bonds()
{
  if (list == NULL) error->all(FLERR, "NeighList is not initialized in FixCoupling");

  neighbor->build_one(list);

  int *ilist = list->ilist;
  int nlocal = atom->nlocal;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int inum = list->inum;
  double **x = atom->x;
  int *mask = atom->mask;
  tagint *tag = atom->tag;

  coupling_list.clear();

  vector<tagint> coupling_me;    // use -1 as delimiter
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if (mask[i] & groupbit) {
      coupling_me.push_back(tag[i]);
      for (int jj = 0; jj < numneigh[i]; jj++) {
        int j = firstneigh[i][jj];
        j &= NEIGHMASK;
        if (mask[j] & jgroupbit) {
          double dx = x[j][0] - x[i][0];
          double dy = x[j][1] - x[i][1];
          double dz = 0;    //x[j][2] - x[i][2];
          // domain->minimum_image(dx, dy, dz);
          double rsq = dx * dx + dy * dy + dz * dz;
          if (rsq < L0sq) { coupling_me.push_back(tag[j]); }
        }
      }
      coupling_me.push_back(-1);
    }
  }

  // gather coupling_list_local to coupling_list
  vector<int> data_sizes(nProcs, 0);
  vector<int> displacements(nProcs, 0);
  int m = coupling_me.size();
  MPI_Allgather(&m, 1, MPI_INT, data_sizes.data(), 1, MPI_INT, world);
  int totalSize = accumulate(data_sizes.begin(), data_sizes.end(), 0);
  partial_sum(data_sizes.begin(), data_sizes.end() - 1, displacements.begin() + 1);

  vector<tagint> coupling_list_all(totalSize);
  MPI_Allgatherv(coupling_me.data(), m, MPI_LMP_TAGINT, coupling_list_all.data(), data_sizes.data(),
                 displacements.data(), MPI_LMP_TAGINT, world);

  int i = 0;
  vector<tagint> temp;
  while (i < totalSize) {
    if (coupling_list_all[i] == -1) {
      coupling_list.push_back(temp);
      temp.clear();
    } else
      temp.push_back(coupling_list_all[i]);

    i++;
  }
}

void FixCoupling::update_bonds()
{
  if (list == NULL) error->all(FLERR, "NeighList is not initialized in FixCoupling");
  neighbor->build_one(list);
  int inum = list->inum;
  int nlocal = atom->nlocal;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  double **x = atom->x;
  int *mask = atom->mask;
  // check if existing bonds are ruptured or not
  vector<tagint> coupling_me;
  for (int i = 0; i < coupling_list.size(); i++) {
    tagint tag1 = coupling_list[i][0];
    vector<tagint> tag_array = coupling_list[i];
    int id1 = atom->map(tag1);
    if (id1 != -1 && id1 < atom->nlocal) {
      coupling_me.push_back(tag1);
      // push back unruptured bonds
      for (int j = 1; j < tag_array.size(); j++) {
        int id2 = atom->map(tag_array[j]);
        if (id2 != -1) {
          double dx = x[id2][0] - x[id1][0];
          double dy = x[id2][1] - x[id1][1];
          double dz = 0;    //atom->x[id2][2] - atom->x[id1][2];
          // domain->minimum_image(dx, dy, dz);
          double rsq = dx * dx + dy * dy + dz * dz;
          if (rsq <= Lcutsq) { coupling_me.push_back(tag_array[j]); }
        }
      }
      // push back new bonds
      for (int j = 0; j < numneigh[id1]; j++) {
        int id2 = firstneigh[id1][j];
        id2 &= NEIGHMASK;
        if (mask[id2] & jgroupbit) {
          double dx = x[id2][0] - x[id1][0];
          double dy = x[id2][1] - x[id1][1];
          double dz = 0;    //atom->x[id2][2] - atom->x[id1][2];
          // domain->minimum_image(dx, dy, dz);
          double rsq = dx * dx + dy * dy + dz * dz;
          if (rsq < L0sq &&
              find(tag_array.begin(), tag_array.end(), atom->tag[id2]) == tag_array.end()) {
            coupling_me.push_back(atom->tag[id2]);
          }
        }
      }

      coupling_me.push_back(-1);
    }
  }

  // gather coupling_list_local to coupling_list
  vector<int> data_sizes(nProcs, 0);
  vector<int> displacements(nProcs, 0);
  int m = coupling_me.size();
  MPI_Allgather(&m, 1, MPI_INT, data_sizes.data(), 1, MPI_INT, world);
  int totalSize = accumulate(data_sizes.begin(), data_sizes.end(), 0);
  partial_sum(data_sizes.begin(), data_sizes.end() - 1, displacements.begin() + 1);

  vector<tagint> coupling_list_all(totalSize);
  MPI_Allgatherv(coupling_me.data(), m, MPI_LMP_TAGINT, coupling_list_all.data(), data_sizes.data(),
                 displacements.data(), MPI_LMP_TAGINT, world);

  coupling_me.clear();
  // cout<< "Coupling_me sizes: " << coupling_me.size() << endl;

  int i = 0;
  vector<tagint> temp;
  coupling_list.clear();
  while (i < totalSize) {
    if (coupling_list_all[i] == -1) {
      coupling_list.push_back(temp);
      temp.clear();
    } else
      temp.push_back(coupling_list_all[i]);

    i++;
  }

  // for (int i = 0; i < coupling_list.size(); i++) {
  //   for (int j = 0; j < coupling_list[i].size(); j++) { cout << coupling_list[i][j] << " "; }
  //   cout << endl;
  // }
  // cout << "Coupling_list sizes: " << coupling_list.size() << endl;
  // cout << "time stemp is: " << update->ntimestep << endl;
}

void FixCoupling::initialize_bonds_check_age()
{
  if (list == NULL) error->all(FLERR, "NeighList is not initialized in FixCoupling");

  neighbor->build_one(list);

  int *ilist = list->ilist;
  int nlocal = atom->nlocal;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int inum = list->inum;
  double **x = atom->x;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  int *age = atom->age;

  coupling_list.clear();

  vector<tagint> coupling_me;    // use -1 as delimiter
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if (mask[i] & groupbit) {
      coupling_me.push_back(tag[i]);
      for (int jj = 0; jj < numneigh[i]; jj++) {
        int j = firstneigh[i][jj];
        j &= NEIGHMASK;
        if (mask[j] & jgroupbit) {
          double dx = x[j][0] - x[i][0];
          double dy = x[j][1] - x[i][1];
          double dz = 0;    //x[j][2] - x[i][2];
          // domain->minimum_image(dx, dy, dz);
          double rsq = dx * dx + dy * dy + dz * dz;
          if (rsq < L0sq && age[j] >= age_min && age[j] <= age_max ) { coupling_me.push_back(tag[j]); }
        }
      }
      coupling_me.push_back(-1);
    }
  }

  // gather coupling_list_local to coupling_list
  vector<int> data_sizes(nProcs, 0);
  vector<int> displacements(nProcs, 0);
  int m = coupling_me.size();
  MPI_Allgather(&m, 1, MPI_INT, data_sizes.data(), 1, MPI_INT, world);
  int totalSize = accumulate(data_sizes.begin(), data_sizes.end(), 0);
  partial_sum(data_sizes.begin(), data_sizes.end() - 1, displacements.begin() + 1);

  vector<tagint> coupling_list_all(totalSize);
  MPI_Allgatherv(coupling_me.data(), m, MPI_LMP_TAGINT, coupling_list_all.data(), data_sizes.data(),
                 displacements.data(), MPI_LMP_TAGINT, world);

  int i = 0;
  vector<tagint> temp;
  while (i < totalSize) {
    if (coupling_list_all[i] == -1) {
      coupling_list.push_back(temp);
      temp.clear();
    } else
      temp.push_back(coupling_list_all[i]);

    i++;
  }
}

void FixCoupling::update_bonds_check_age()
{
  if (list == NULL) error->all(FLERR, "NeighList is not initialized in FixCoupling");
  neighbor->build_one(list);
  int inum = list->inum;
  int nlocal = atom->nlocal;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  double **x = atom->x;
  int *mask = atom->mask;
  int *age = atom->age;
  // check if existing bonds are ruptured or not
  vector<tagint> coupling_me;
  for (int i = 0; i < coupling_list.size(); i++) {
    tagint tag1 = coupling_list[i][0];
    vector<tagint> tag_array = coupling_list[i];
    int id1 = atom->map(tag1);
    if (id1 != -1 && id1 < atom->nlocal) {
      coupling_me.push_back(tag1);
      // push back unruptured bonds
      for (int j = 1; j < tag_array.size(); j++) {
        int id2 = atom->map(tag_array[j]);
        if (id2 != -1) {
          double dx = x[id2][0] - x[id1][0];
          double dy = x[id2][1] - x[id1][1];
          double dz = 0;    //atom->x[id2][2] - atom->x[id1][2];
          // domain->minimum_image(dx, dy, dz);
          double rsq = dx * dx + dy * dy + dz * dz;
          if (rsq <= Lcutsq && age[j] >= age_min && age[j] <= age_max ) { coupling_me.push_back(tag_array[j]); }
        }
      }
      // push back new bonds
      for (int j = 0; j < numneigh[id1]; j++) {
        int id2 = firstneigh[id1][j];
        id2 &= NEIGHMASK;
        if (mask[id2] & jgroupbit) {
          double dx = x[id2][0] - x[id1][0];
          double dy = x[id2][1] - x[id1][1];
          double dz = 0;    //atom->x[id2][2] - atom->x[id1][2];
          // domain->minimum_image(dx, dy, dz);
          double rsq = dx * dx + dy * dy + dz * dz;
          if (rsq < L0sq &&
              find(tag_array.begin(), tag_array.end(), atom->tag[id2]) == tag_array.end() && age[id2] >= age_min && age[id2]<= age_max) {
            coupling_me.push_back(atom->tag[id2]);
          }
        }
      }

      coupling_me.push_back(-1);
    }
  }

  // gather coupling_list_local to coupling_list
  vector<int> data_sizes(nProcs, 0);
  vector<int> displacements(nProcs, 0);
  int m = coupling_me.size();
  MPI_Allgather(&m, 1, MPI_INT, data_sizes.data(), 1, MPI_INT, world);
  int totalSize = accumulate(data_sizes.begin(), data_sizes.end(), 0);
  partial_sum(data_sizes.begin(), data_sizes.end() - 1, displacements.begin() + 1);

  vector<tagint> coupling_list_all(totalSize);
  MPI_Allgatherv(coupling_me.data(), m, MPI_LMP_TAGINT, coupling_list_all.data(), data_sizes.data(),
                 displacements.data(), MPI_LMP_TAGINT, world);

  coupling_me.clear();
  // cout<< "Coupling_me sizes: " << coupling_me.size() << endl;

  int i = 0;
  vector<tagint> temp;
  coupling_list.clear();
  while (i < totalSize) {
    if (coupling_list_all[i] == -1) {
      coupling_list.push_back(temp);
      temp.clear();
    } else
      temp.push_back(coupling_list_all[i]);

    i++;
  }
}