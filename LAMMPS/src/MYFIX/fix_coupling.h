#ifdef FIX_CLASS
// clang-format off
FixStyle(coupling, FixCoupling);
// clang-format on
#else

#ifndef LMP_FIX_COUPLING_H
#define LMP_FIX_COUPLING_H

#include "fix.h"
#include <tuple>
#include <vector>

namespace LAMMPS_NS {
  class FixCoupling : public Fix {
   public:
    FixCoupling(class LAMMPS *, int, char **);
    int setmask() override;
    void init() override;
    void setup(int) override;
    void pre_force(int) override;
    void post_force(int) override;
    void init_list(int, NeighList *) override;
	int pack_reverse_comm(int, int, double *) override;
	void unpack_reverse_comm(int, int *, double *) override;
   private:
    char *style;
    double K, L0, Lcut;
	  double L0sq,Lcutsq;
    int jgroup, jgroupbit;
    class NeighList *list;
    double cutoff_plus_skin;
    int check_age_flag = 0, age_min = 0,age_max = 0;

	int nProcs;
	std::vector<std::vector<double> >fOnProc;
	int init_update,next_update;

    std::vector<std::vector<tagint>> coupling_list; // global variable: 0-th element is the tag of the central atom

   
	void initialize_bonds();
  void initialize_bonds_check_age();
	void update_bonds();
  void update_bonds_check_age();
  };

}    // namespace LAMMPS_NS

#endif
#endif
