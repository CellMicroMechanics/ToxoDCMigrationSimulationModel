#ifdef FIX_CLASS
// clang-format off
FixStyle(addatom2d, FixAddAtom2D);
// clang-format on
#else

#ifndef LMP_FIX_ADDATOM2D_H
#define LMP_FIX_ADDATOM2D_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS
{
	class FixAddAtom2D : public Fix
	{
	public:
		FixAddAtom2D(class LAMMPS *, int, char **);
		int setmask() override;
		void pre_exchange() override;
		// void initial_integrate(int) override;

	private:
		double px; // polarity direction x
		double py; // polarity direction y
        double rate; // add atoms rate
        class RanMars *rng;
		// methods
		//void deleteAtomAtBack();

        
    };

} // namespace LAMMPS_NS

#endif
#endif
