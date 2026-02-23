#ifdef FIX_CLASS
// clang-format off
FixStyle(edge2d, FixEdge2D);
// clang-format on
#else

#ifndef LMP_FIX_EDGE2D_H
#define LMP_FIX_EDGE2D_H

#include "fix.h"

namespace LAMMPS_NS
{
	class FixEdge2D : public Fix
	{
	public:
		FixEdge2D(class LAMMPS *, int, char **);
		int setmask() override;
		void post_force(int) override;
		// void initial_integrate(int) override;

	private:
		double K; // spring constant
		double L0; // prefered spring length

		double **pos;
		int nlocal;
		double **f;
		int *mask;
		// methods
    };

} // namespace LAMMPS_NS

#endif
#endif
