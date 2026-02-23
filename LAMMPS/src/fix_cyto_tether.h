#ifdef FIX_CLASS
// clang-format off
FixStyle(cyto/tether, FixCytoTether);
// clang-format on
#else

#ifndef LMP_FIX_CYTO_TETHER_H
#define LMP_FIX_CYTO_TETHER_H

#include "fix.h"
#include "fix_membrane2d.h"
#include "fix_membrane2dfull.h"
#include <tuple>
#include <vector>

namespace LAMMPS_NS {
  class FixCytoTether : public Fix {
   public:
    FixCytoTether(class LAMMPS *, int, char **);
    int setmask() override;
    void init() override;
    void post_force(int) override;
  private:
    double k, xc;
    double L0=0;
    std::string fixID1,fixID2;
    int fixIdx,fixIdx2;
    FixMembrane2d *fixmembrane;
    FixMembrane2d *fixmembrane2;
    FixMembrane2dFull *fixmembranefull;  
    FixMembrane2dFull *fixmembrane2full;
  };

}    // namespace LAMMPS_NS

#endif
#endif