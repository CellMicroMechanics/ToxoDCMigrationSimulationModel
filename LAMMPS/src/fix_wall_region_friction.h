#ifdef FIX_CLASS
// clang-format off
FixStyle(wall/region/friction,FixWallRegionFriction);
// clang-format on
#else

#ifndef LMP_FIX_WALL_REGION_FRICTION_H
#define LMP_FIX_WALL_REGION_FRICTION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallRegionFriction : public Fix {
 public:
  FixWallRegionFriction(class LAMMPS *, int, char **);
  ~FixWallRegionFriction() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  void end_of_step() override;
  double compute_scalar() override;
  double compute_vector(int) override;

  class Region *region;


 private:
  int style;
  double epsilon, sigma,cutoff;
  double alpha;
  int eflag;
  double ewall[4], ewall_all[4];
  int ilevel_respa;
  char *idregion;
  

  double drag; // damp
  double mu; // static yield coeff

  double coeff1, coeff2, coeff3, coeff4, offset;
  double coeff5, coeff6, coeff7;
  double eng, fwall;

  int next_print, nevery_print;
  //char* filename;
  FILE *fp;


  void lj93(double);
  void lj126(double);
  void lj1043(double);
  void morse(double);
  void colloid(double, double);
  void harmonic(double);
};

}    // namespace LAMMPS_NS

#endif
#endif