#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(create_atoms_onsphere,CreateAtomsOnSphere);
// clang-format on
#else

#ifndef LMP_CREATE_ATOMS_ONSPHERE_H
#define LMP_CREATE_ATOMS_ONSPHERE_H

#include "command.h"
#include "create_atoms.h"

namespace LAMMPS_NS {
class CreateAtomsOnSphere : public CreateAtoms {
 public:
  CreateAtomsOnSphere(class LAMMPS *);
  void command(int, char **) override;

 private:
  
  int nAtoms,atomType;
  double x0,y0,z0,radius,radius_minor;
  double sublo[3],subhi[3];
  void create_atoms_on_sphere();
  void create_atoms_on_circle();
  void add_one_atom(double*,int);
  
};
}
#endif
#endif