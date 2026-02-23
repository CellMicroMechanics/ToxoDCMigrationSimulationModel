#ifdef FIX_CLASS
// clang-format off
FixStyle(tdpd/force/contractile,FixTDPDForceContractile);
// clang-format on
#else

#ifndef LMP_FIX_TDPD_FORCE_CONTRACTILE_H
#define LMP_FIX_TDPD_FORCE_CONTRACTILE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTDPDForceContractile : public Fix {
 public:
  FixTDPDForceContractile(class LAMMPS *, int, char **);
  int setmask() override;
  // void init() override;
  // void init_list(int,class NeighList*) override;
  void post_force(int) override;
 protected:
    int cc_index; // species of F actin
    double zeta; // active force coefficient
    // class NeighList *list;   
    // class Pair *pair;
    // class PairTDPD *pair_tdpd; // pointer to the pair style
};
}    // namespace LAMMPS_NS
#endif
#endif
