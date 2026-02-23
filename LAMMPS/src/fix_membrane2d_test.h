#ifdef FIX_CLASS
// clang-format off
FixStyle(membrane2d/test, FixMembrane2dTest);
// clang-format on
#else

#ifndef LMP_FIX_MEMBRANE2DTEST_H
#define LMP_FIX_MEMBRANE2DTEST_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS
{
	class FixMembrane2dTest : public Fix
	{
	public:
		FixMembrane2dTest(class LAMMPS *, int, char **);
		int setmask() override;
		void setup(int) override;
		void pre_force(int) override;
		void post_force(int) override;
		void pre_exchange() override;
		void end_of_step() override;
		//void end_of_step() override;
		int pack_reverse_comm(int, int, double *) override;
		void unpack_reverse_comm(int, int *, double *) override;

		// this is a global variable
		// std::vector<std::vector<tagint>> bondedTags; // store tags of bounded atoms for each atom
												     // sorted by the first element from -pi to pi regarding to COM  
		double *xc; // center of mass of the membrane; (global variable)
		int nprocs; // number of processors

	private:
		// stretching
		double K; // spring constant
		double L0; // prefered spring length
		double KC, theta0;// bending constant and prefered bending angle
        double KA, A0; // area constant and prefered area
		double lref; // reference length of the membrane;

        bool isPolymerisationEnabled = false;
		double px, py,pz; // predefined polarity
		

		int membraneBond_init_flag = 0;

        int depolymerization_phase, depolymerization_phase_gap;

		int dim;
		int isChanged;

		std::vector<double> tension_coeff; // tension coefficient (size same as nghost+nlocal)

		// for polymerization rate
		class RanMars *rng;
		int seed; // random seed
		int seed_depoly; // random seed for depolymerization
		double rate; // polymerization rate
		int newton_bond; 
	
		double sublo[3], subhi[3];
		std::vector< std::vector<double>> fOnProc;

		

		void updatePredefinedTension(); // set predefined tension for each bond based on polarity
		double tension_grad0; // constant gradient of tension 		
		double tension0; // constant tension for testing
		int tension_flag;
		
		tagint maxLengthTag1, maxLengthTag2; // tags of the atoms that have the maximum edge length
		tagint minLengthTag1, minLengthTag2; // tags of the atoms that have the minimum edge length
		void findMaxEdgeLengthAndTag(); // find the maximum edge length of the membrane
		void findMinEdgeLengthAndTag();
		void applyMembraneDynamics(tagint,tagint,int); // apply polymerization and depolymerization
		// std::vector<tagint> extendAtomTags(std::vector<tagint>);

		void initializeMembraneBonds2D(); //initialize membrane bonds counterclockwise
		void getCOM();
		
		double Etension, Etension0; 
		int init_update;
		void computeTensionEnergy();
    };

} // namespace LAMMPS_NS

#endif
#endif
