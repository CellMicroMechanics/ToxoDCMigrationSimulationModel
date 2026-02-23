#ifdef FIX_CLASS
// clang-format off
FixStyle(cortex2d, FixCortex2d);
// clang-format on
#else

#ifndef LMP_FIX_CORTEX2D_H
#define LMP_FIX_CORTEX2D_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS
{
	class FixCortex2d : public Fix
	{
	public:
		FixCortex2d(class LAMMPS *, int, char **);
		int setmask() override;
		void post_force(int) override;
		void pre_exchange() override;
		void end_of_step() override;
		void init() override;

		std::vector<std::vector<tagint>> bondedTags; // store tags of bounded atoms for each atom
												     // sorted by the first element from -pi to pi regarding to COM  

		double *xc; // center of mass of the cortex;

		//void init_list(int, NeighList *) override;
		// void initial_integrate(int) override;

	private:
		// stretching
		double K; // spring constant
		double L0; // prefered spring length
		double KC, theta0;
		double px, py,pz; // predefined polarity

		int dim;
		bool isChanged;

		// for polymerization rate
		class RanMars *rng;
		int seed; // random seed
		double rate; // polymerization rate



		int maxnumBonds; 

		// methods
		tagint deleteOneAtomAtLargestNegativePolarity(std::vector<double>&, std::vector<int> &);
		std::vector<tagint> addOneAtomAtLargestPositivePolarity(std::vector<double>&, std::vector<int> &);
		void initializeMembraneBonds2D();
		void updateMembraneBonds2D(std::vector<tagint>,tagint);
		// void updateMembraneBonds3D();
		void getCOM();
		//int findAtomOfTag(tagint);
    };

} // namespace LAMMPS_NS

#endif
#endif
