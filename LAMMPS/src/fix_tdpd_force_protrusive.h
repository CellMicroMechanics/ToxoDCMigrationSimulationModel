#ifdef FIX_CLASS
// clang-format off
FixStyle(tdpd/force/protrusive,FixTDPDForceProtrusive);
// clang-format on
#else

#ifndef LMP_FIX_TDPD_FORCE_PROTRUSIVE_H
#define LMP_FIX_TDPD_FORCE_PROTRUSIVE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTDPDForceProtrusive : public Fix {
 public:
  FixTDPDForceProtrusive(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;
  void init_list(int,class NeighList*) override;
  void post_force(int) override;
 protected:
    int cc_index; // species of F actin
    class NeighList *list;
    double zeta; // active force coefficient
    class Pair *pair;
    class PairTDPD *pair_tdpd; // pointer to the pair style
};
}    // namespace LAMMPS_NS
#endif
#endif
