#ifdef FIX_CLASS
// clang-format off
FixStyle(polytest, FixPolyTest);
// clang-format on
#else

#ifndef LMP_FIX_POLYTEST_H
#define LMP_FIX_POLYTEST_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS
{
	class FixPolyTest : public Fix
	{
	public:
		FixPolyTest(class LAMMPS *, int, char **);
		int setmask() override;
		void pre_exchange() override;
		void init() override;

	private:
        int jgroup,jgroupbit,seed,dim,next_reneighbor;
        class RanMars *rng;
        double px, py, pz, rate;

        double xc[3];

        void getCOM();
        void addOneAtomAtLargestPositivePolarity(std::vector<double> &, std::vector<int> &);
        void deleteOneAtomAtLargestNegativePolarity(std::vector<double> &, std::vector<int> &);
        
    };

} // namespace LAMMPS_NS

#endif
#endif
