#include "atom.h"
#include "atom_masks.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "update.h"

#include "domain.h"
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

#include "fix_cyto_tether.h"
#include "fix_membrane2d.h"
#include "fix_membrane2dfull.h"
#include "modify.h"

using namespace std;
using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double SMALL = 0.001;

FixCytoTether::FixCytoTether(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg), fixmembrane(NULL), fixmembranefull(NULL), fixmembrane2(NULL), fixmembrane2full(NULL)
{
  fixID1 = string(arg[3]);
  fixID2 = string(arg[4]);
  k = utils::numeric(FLERR, arg[5], false, lmp);                   // spring constant
  if (narg > 6) L0 = utils::numeric(FLERR, arg[6], false, lmp);    // equilibrium length
}

int FixCytoTether::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

void FixCytoTether::init()
{
  fixIdx = modify->find_fix(
      fixID1);    // this fix find the membrane (fixID1) xc of the current group and
  fixIdx2 = modify->find_fix(fixID2);    // couples it to the xc of the target membrane (fixID2)
  if (fixIdx == -1) error->all(FLERR, "Invalid fix ID");
  Fix *fix = modify->fix[fixIdx];
  if (strcmp(fix->style, "membrane2d") == 0) {
    fixmembrane = dynamic_cast<FixMembrane2d *>(fix);
  } else if (strcmp(fix->style, "membrane2dfull") == 0) {
    fixmembranefull = dynamic_cast<FixMembrane2dFull *>(fix);
  } else {
    error->all(FLERR, "Found fix is not membrane2d style.");
  }

  if (fixmembrane == NULL && fixmembranefull == NULL)
    error->all(FLERR, "Failed to cast fix (Position 1) to membrane2d/membrane2dfull style.");

  if (fixIdx2 == -1) error->all(FLERR, "Invalid fix ID");
  Fix *fix2 = modify->fix[fixIdx2];
  if (strcmp(fix2->style, "membrane2d") == 0) {
    fixmembrane2 = dynamic_cast<FixMembrane2d *>(fix2);
  } else if (strcmp(fix2->style, "membrane2dfull") == 0) {
    fixmembrane2full = dynamic_cast<FixMembrane2dFull *>(fix2);
  } else {
    error->all(FLERR, "Found fix is not membrane2d style.");
  }

  if (fixmembrane2 == NULL && fixmembrane2full == NULL)
    error->all(FLERR, "Failed to cast fix (Position 2) to membrane2d/membrane2dfull style.");
}

void FixCytoTether::post_force(int)
{
  double *myxc, *hisxc;

  int jgroupbit; 
  if (fixmembrane != NULL)
  {
    myxc = fixmembrane->xc;
  }   
  else
  {
    myxc = fixmembranefull->xc;
  }
  if (fixmembrane2 != NULL)
  {
    hisxc = fixmembrane2->xc;
    jgroupbit = fixmembrane2->groupbit;
  }
  else
  {
    hisxc = fixmembrane2full->xc;
    jgroupbit = fixmembrane2full->groupbit;
  }

  double dx = hisxc[0] - myxc[0];
  double dy = hisxc[1] - myxc[1];
  double dz = hisxc[2] - myxc[2];

  domain->minimum_image(dx, dy, dz);
  double dr = sqrt(dx * dx + dy * dy + dz * dz);

  double fmag = k * (dr - L0);

  // get number of particles in the current group (membrane)

  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int nglocal = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) nglocal++;
  }

  int ngglobal;
  MPI_Allreduce(&nglocal, &ngglobal, 1, MPI_INT, MPI_SUM, world);

  double fmag1 = fmag / ngglobal;

  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      f[i][0] += fmag1 * dx / dr;
      f[i][1] += fmag1 * dy / dr;
      f[i][2] += fmag1 * dz / dr;
    }
  }

  

  int nglocal2 = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & jgroupbit) nglocal2++;
  }
  int ngglobal2;
  MPI_Allreduce(&nglocal2, &ngglobal2, 1, MPI_INT, MPI_SUM, world);
  double fmag2 = fmag / ngglobal2;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & jgroupbit) {
      f[i][0] -= fmag2 * dx / dr;
      f[i][1] -= fmag2 * dy / dr;
      f[i][2] -= fmag2 * dz / dr;
    }
  }
}