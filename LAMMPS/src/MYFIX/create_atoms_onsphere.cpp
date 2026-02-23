#include "create_atoms_onsphere.h"
#include "atom.h"
#include "atom_vec.h"
#include "create_atoms.h"
#include "domain.h"
#include "error.h"
#include "math_const.h"
#include <cmath>
#include <iostream>
#include "comm.h"

using namespace LAMMPS_NS;
using namespace std;
using MathConst::MY_2PI;
using MathConst::MY_PI;

CreateAtomsOnSphere::CreateAtomsOnSphere(LAMMPS *lmp) : CreateAtoms(lmp) {}

void CreateAtomsOnSphere::command(int narg, char **arg)
{
  // Check parameters
  atomType = utils::inumeric(FLERR, arg[0], false, lmp);
  nAtoms = utils::inumeric(FLERR, arg[1], false, lmp);
  x0 = utils::numeric(FLERR, arg[2], false, lmp);
  y0 = utils::numeric(FLERR, arg[3], false, lmp);
  z0 = utils::numeric(FLERR, arg[4], false, lmp);
  radius = utils::numeric(FLERR, arg[5], false, lmp);
  if(narg >= 7) radius_minor = utils::numeric(FLERR, arg[6], false, lmp);
  else radius_minor = radius;

  sublo[0] = domain->sublo[0];
  subhi[0] = domain->subhi[0];
  sublo[1] = domain->sublo[1];
  subhi[1] = domain->subhi[1];
  sublo[2] = domain->sublo[2];
  subhi[2] = domain->subhi[2];
  //cout<< "Simualtion dimension: " << domain->dimension << endl;
  if(domain->dimension == 3)
    create_atoms_on_sphere();
  else if(domain->dimension == 2)
    create_atoms_on_circle();
  else
    error->all(FLERR, "This command is only available in 2D or 3D simulations");
}

void CreateAtomsOnSphere::create_atoms_on_circle()
{
  // Implementation of create_atoms_on_circle goes here
  MPI_Barrier(world);
  double time1 = platform::walltime();

  bigint natoms_previous = atom->natoms;
  int nlocal_previous = atom->nlocal;

  if (atom->map_style != Atom::MAP_NONE) atom->map_clear();
  atom->nghost = 0;
  atom->avec->clear_bonus();

  double dtheta = MY_2PI / nAtoms;
  double xnew[3];

  for (int i = 0; i < nAtoms; i++) {
    double theta = i * dtheta;
    xnew[0] = x0 + radius * cos(theta);
    xnew[1] = y0 + radius_minor * sin(theta);
    xnew[2] = z0;
    add_one_atom(xnew, atomType);
  }

  atom->data_fix_compute_variable(nlocal_previous, atom->nlocal);
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT) error->all(FLERR, "Too many total atoms");
  
  atom->tag_extend();
  atom->tag_check();
  if (atom->map_style != Atom::MAP_NONE) {
    atom->map_init();
    atom->map_set();
  }
  if(comm->me == 0)
    cout<<"Total number of atoms created: "<<atom->natoms<<endl;
}

void CreateAtomsOnSphere::create_atoms_on_sphere()
{
  // Implementation of create_atoms_in_sphere goes here
  MPI_Barrier(world);
  double time1 = platform::walltime();

  bigint natoms_previous = atom->natoms;
  int nlocal_previous = atom->nlocal;

  if (atom->map_style != Atom::MAP_NONE) atom->map_clear();
  atom->nghost = 0;
  atom->avec->clear_bonus();

  // double ds = sqrt(4.0 * MY_PI * radius * radius / nAtoms);
  // double dtheta = ds / radius;
  // double dphi = ds / radius;
 
  double theta, phi;
  double xnew[3];

  double golden_angle = MY_PI * (3 - sqrt(5));  // Golden angle in radians

for (int i = 1; i <= nAtoms; i++) {
    double z = 1 - (2.0 * i / double(nAtoms));  // 1 - 2 * (i / nAtoms)
    double r = sqrt(1 - z*z);

    double phi = golden_angle * i;

    double x = cos(phi) * r;
    double y = sin(phi) * r;

    xnew[0] = x0 + x * radius;
    xnew[1] = y0 + y * radius;
    xnew[2] = z0 + z * radius;

      add_one_atom(xnew, atomType);
}
  // for (int i = 0; i < atom->nlocal; i++) {
  //   cout << "atom ID: " << atom->tag[i] << ", tpye: " << atom->type[i] << ", x: " << atom->x[i][0]
  //        << ", y: " << atom->x[i][1] << ", z: " << atom->x[i][2] << endl;
  // }
  atom->data_fix_compute_variable(nlocal_previous, atom->nlocal);
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT) error->all(FLERR, "Too many total atoms");



  atom->tag_extend();
  atom->tag_check();
  if (atom->map_style != Atom::MAP_NONE) {
    atom->map_init();
    atom->map_set();
  }
  if(comm->me == 0)
    cout<<"Total number of atoms created: "<<atom->natoms<<endl;
}

void CreateAtomsOnSphere::add_one_atom(double *xnew, int atomType)
{

  if (xnew[0] >= sublo[0] && xnew[0] < subhi[0] && xnew[1] >= sublo[1] && xnew[1] < subhi[1] &&
      xnew[2] >= sublo[2] && xnew[2] < subhi[2]) {
    imageint imagetmp = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
    domain->remap(xnew, imagetmp);
    atom->avec->create_atom(atomType, xnew);
  }
}