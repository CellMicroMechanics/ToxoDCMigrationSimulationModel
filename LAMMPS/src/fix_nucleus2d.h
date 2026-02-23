#ifdef FIX_CLASS
// clang-format off
FixStyle(nucleus2d, FixNucleus2D);
// clang-format on
#else

#ifndef LMP_FIX_NUCLEUS2D_H
#define LMP_FIX_NUCLEUS2D_H

#include "fix.h"
#include <vector>
#include <tuple>

namespace LAMMPS_NS
{
    /* The difference between 2d nulceus and 2d cortex is the following.
    1. Nucleus does not have dynamically added or deleted particles
    2. Nucleus can determine whether bending is allowed or not
    3. Nuclear stiffness is achieved by large repulsive dpd particles
    */
	class FixNucleus2D : public Fix
	{
	public:
		FixNucleus2D(class LAMMPS *, int, char **);
        void init();
        void setup(int);
        int setmask();
        void post_force(int);

	private:

        double K; // elastic bonds
        double L0; // rest length of bonds
        double KC; // bending elasticity due to nuclear lamina
        double theta0; // rest angle
        double theta_thres; // threshold of bending rigidity. (theta > threshold -> bending, theta < threshold -> no bending, lamina rupture)
        double KA; //area elasticity 
        double A0; //rest area
        std::vector<std::vector<tagint>> bondedTags;
        std::vector<bool> isBreak; 
        void initializeNuclearBonds(); // Note that the topology of nuclear particles should remain the same during simulations
    };

} // namespace LAMMPS_NS

#endif
#endif