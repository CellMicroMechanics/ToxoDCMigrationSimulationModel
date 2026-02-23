#ifdef FIX_CLASS
// clang-format off
FixStyle(tdpd/fg/reaction,FixTDPDFGReaction);
// clang-format on
#else

#ifndef LMP_FIX_TDPD_REACTION_H
#define LMP_FIX_TDPD_REACTION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTDPDFGReaction : public Fix {
 public:
  FixTDPDFGReaction(class LAMMPS *, int, char **);
  ~FixTDPDFGReaction() override;
  int setmask() override;
  void post_force(int) override;

 protected:
    int cc_index;
    double rho, n, K,
    mu; // rho: production rate, 
        // n: Hill coefficient
        // K: half saturation constant
        // mu: degradation rate
    int jgroupbit; // local source group
    int cc_index_pool; 
    double npool, Kpool;

    int cc_index_inhibitor; // inhibitor species (optional)
    double Kinhibitor, ninhibitor; // half saturation constant for inhibitor;

};
}    // namespace LAMMPS_NS
#endif
#endif
