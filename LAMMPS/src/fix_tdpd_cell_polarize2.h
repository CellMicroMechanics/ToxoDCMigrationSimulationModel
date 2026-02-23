#ifdef FIX_CLASS
// clang-format off
FixStyle(tdpd/cell/polarize2,FixTDPDCellPolarize2);
// clang-format on
#else

#ifndef LMP_FIX_TDPD_CELL_POLARIZE2_H
#define LMP_FIX_TDPD_CELL_POLARIZE2_H

#include "fix.h"
#include <array>

namespace LAMMPS_NS {

class FixTDPDCellPolarize2 : public Fix {
 public:
  FixTDPDCellPolarize2(class LAMMPS *, int, char **);
  ~FixTDPDCellPolarize2() override;
  int setmask() override;
  void init() override;
  void pre_force(int) override;

private:
  void ChooseOneReferenceParticle();
  void ChooseReferenceParticleCCThres();
  void UpdateReferenceInfo(std::vector<tagint> &, std::vector<double> &, std::vector<imageint> &);
  void SetReferenceParticleTags(std::vector<tagint> &, double);



  void InitializePoolTags();
  std::vector<tagint> ChooseRandomParticleFromPool(int);
  void SetSourceFromNReferenceParticles();

  //Dedicated stochastic process to select reference particles
  tagint MembraneRandomWalk(tagint, double);
  tagint MembraneOrnsteinUhlenbeck(tagint, tagint, double, double);
  void InitializeOrderedTags();


 protected:
    int seed;
    int j1groupbit;// cortex group, find a reference position
    double radius_ref = 0;
    double radius; // radius of the region containing the source particles
    int j2groupbit; // assign polarized particles to this group (in which the source term will be applied)
    class RanMars *rng;

    std::vector<tagint> pool_tags;
    std::vector<tagint> ref_tags;
    std::vector<double> xref_global;
    std::vector<imageint> imgref_global;

    int nprocs;
    int nsteps; // frequency of the reference position selection
    int next_selection; // next selection of the reference position

    int mode;
    imageint imgref;
    int cc_index;
    double cc_thres;

    // ref_tags by dedicated stochastic process (Ornstein-Uhlenbeck process)
    double p; // probability of membrane random walk
    double tau; // time constant for the Ornstein-Uhlenbeck process
    double noise; // noise amplitude for the Ornstein-Uhlenbeck process
    tagint polaritySite; // current position of the particle in the Ornstein-Uhlenbeck process
    tagint meanPolaritySite; // mean position of the particle in the Ornstein-Uhlenbeck process
    std::vector<tagint> orderedTags; // tags of the connected particles in the membrane ordered counterclockwisely
    int init_update;
};
}    // namespace LAMMPS_NS
#endif
#endif
