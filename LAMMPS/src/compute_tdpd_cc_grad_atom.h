#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(tdpd/cc_grad/atom,ComputeTDPDCCGradAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_TDPD_CC_GRAD_ATOM_H
#define LMP_COMPUTE_TDPD_CC_GRAD_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTDPDCCGradAtom : public Compute {
 public:
  ComputeTDPDCCGradAtom(class LAMMPS *, int, char **);
  ~ComputeTDPDCCGradAtom() override;
  void init() override;
  void compute_peratom() override;
  double memory_usage() override;

 private:
  int nmax;
  int index;
  int icomponent;
  double *cc_grad_vector;
};

}    // namespace LAMMPS_NS

#endif
#endif
