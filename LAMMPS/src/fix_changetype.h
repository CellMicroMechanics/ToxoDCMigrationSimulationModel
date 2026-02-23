#ifdef FIX_CLASS
// clang-format off
FixStyle(changetype, FixChangeType);
// clang-format on
#else

#ifndef LMP_FIX_CHANGETYPE_H
#define LMP_FIX_CHANGETYPE_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS {
class FixChangeType : public Fix {
 public:
  FixChangeType(class LAMMPS *, int, char **);
  int setmask() override;
  void pre_force(int) override;
  void post_force(int) override;
  // void init_list(int, NeighList*) override;
  // void post_force(int) override;
  // void initial_integrate(int) override;

 private:
  double px;      // polarity direction x
  double py;      // polarity direction y
  double rate;    // add atoms rate

  std::vector<int> inAtomIdx;
  bool changeFlag;
  class RanMars *rng;
  // methods
  int jgroup, jgroupbit;    // target group
  int toType;
  int toMask;
  double toR;

  int toType2;
  int toMask2;
  double toR2;

  //NeighList* list;

  void deleteOutsideParticles();
};

}    // namespace LAMMPS_NS

#endif
#endif
