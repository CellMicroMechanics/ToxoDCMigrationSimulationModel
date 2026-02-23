#ifdef FIX_CLASS
// clang-format off
FixStyle(tdpd/gm/activator,FixTDPDGMActivator);
// clang-format on
#else

#ifndef LMP_FIX_TDPDGMACTIVATOR_H
#define LMP_FIX_TDPDGMACTIVATOR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTDPDGMActivator : public Fix {
 public:
  FixTDPDGMActivator(class LAMMPS *, int, char **);
  ~FixTDPDGMActivator() override;
  int setmask() override;
  void post_force(int) override;

 protected:
  int cc_index1; // activator

  int cc_index2; // inhibitor

  double rho0, rhomax, Ka, na, Kb,nb, mu0, mumax; // parameters for the activator

 



};
}    // namespace LAMMPS_NS
#endif
#endif
