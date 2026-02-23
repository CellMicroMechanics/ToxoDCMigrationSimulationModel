#ifdef FIX_CLASS
// clang-format off
FixStyle(active/drag, FixActiveDrag);
// clang-format on
#else

#ifndef LMP_FIX_ACTIVE_DRAG__H
#define LMP_FIX_ACTIVE_DRAG__H

#include "fix.h"

namespace LAMMPS_NS {
  class FixActiveDrag : public Fix {
   public:
    FixActiveDrag(class LAMMPS *, int, char **);
    int setmask() override;
    void init() override;
    void setup(int) override;
   // void pre_force(int) override;
    void post_force(int) override;
    void init_list(int, NeighList *) override;
   private:
    double Lcut; 
    NeighList *list;
    double ct, cn; // cofficients for tangential and normal components of the active drag (proportional to relative velocity) 
    int jgroup, jgroupbit;
    double cut_comm,cutmax;


  };

}    // namespace LAMMPS_NS

#endif
#endif
