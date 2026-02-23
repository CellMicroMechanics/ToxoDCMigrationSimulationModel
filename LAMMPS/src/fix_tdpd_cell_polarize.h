#ifdef FIX_CLASS
// clang-format off
FixStyle(tdpd/cell/polarize,FixTDPDCellPolarize);
// clang-format on
#else

#ifndef LMP_FIX_TDPD_CELL_POLARIZE_H
#define LMP_FIX_TDPD_CELL_POLARIZE_H

#include "fix.h"
#include <array>

namespace LAMMPS_NS {

class FixTDPDCellPolarize : public Fix {
 public:
  FixTDPDCellPolarize(class LAMMPS *, int, char **);
  ~FixTDPDCellPolarize() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void pre_force(int) override;

private:
  void ChooseOneReferenceParticle();
  void ChooseReferenceParticleCCThres();
  void UpdateReferenceInfo(std::vector<tagint> &, std::vector<double> &, std::vector<imageint> &);
  void SetReferenceParticleTags(std::vector<tagint> &, double);



  void InitializePoolTags();
  std::vector<tagint> ChooseRandomParticleFromPool(int,std::vector<tagint>&);
  void SetSourceFromNReferenceParticles();


  //Dedicated stochastic process to select reference particles
  tagint MembraneSingleSite(tagint,tagint);
  void InitializeOrderedTags();
  std::vector<tagint> MembraneMultipleSites(std::vector<tagint>&, std::vector<tagint> &);

  double compute_pdeath(tagint,bool);

  std::vector<double> GetCOMOfGroup(int);
  tagint GetMaxProjectionTagInGroup(int, std::vector<double>&, double, double, double);

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
    std::vector<tagint> polaritySites;
    std::vector<tagint> meanPolaritySites; // mean polarity sites, used to calculate the mean tag of the reference particles
    std::vector<int> siteStates; // state: 1 -> active, 2-> cd
    std::vector<double> timer;
    std::vector<tagint> orderedTags; // tags of the connected particles in the membrane ordered counterclockwisely
    int init_update;
    // lifetime of a protrusion
    double tau_live0; // initial lifetime of a protrusion
    double tau_cd;
    double rate_birth,p_birth;
    int nMaxMeanPolaritySites; // maximum number of mean polarity sites
  
    // fluctuation of a local protrusion
    double tau; // lifetime of a reference
    double sigma; // noise strength

    bool bias_flag = false; // if true, the reference position is biased by the mean polarity site
    double *vest; // expoential moving average of the velocity of the com velocity.
    // double n; // hill coeff for tau_live saturation
    double K; // sigmoid response of the pseudopod mechanical bias

    // Check obstacles
    int check_obstacle_flag = 0; // 0: disabled, 1: check particle type, 2: check particle group

    std::vector<int> obstacle_types;
    std::vector<int> obstacle_groupbits;
    std::vector<tagint> excluded_tags;
    class NeighList *list;
    std::vector<tagint> FindExcludedTags(std::vector<tagint> &, int);

    // preferred direction
    double px0, py0, pz0;
    double px,py,pz;
    int direction_flag = 0;



};
}    // namespace LAMMPS_NS
#endif
#endif
