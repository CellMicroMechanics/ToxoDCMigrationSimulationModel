#ifdef FIX_CLASS
FixStyle(tdpd/cc/grad,FixTDPDCCGrad);
#else
#ifndef LMP_FIX_TDPDCCGRAD_H
#define LMP_FIX_TDPDCCGRAD_H
#include "fix.h"

namespace LAMMPS_NS {
    class FixTDPDCCGrad : public Fix {
    public:
        FixTDPDCCGrad(class LAMMPS *, int, char **);
        ~FixTDPDCCGrad() override;
        int setmask() override;
        void init() override;
        void init_list(int, NeighList *) override;
        void pre_force(int) override;

    private:
        void compute_cc_grad(int);

    protected:
        int *cc_index; // of which concentration field to compute the gradient
        class NeighList *list; // neighbor list
        class Pair *pair;
        class PairTDPD *pair_tdpd; // pointer to the pair style
        int mode;
        int cc_index_size = 0;
    };
}
#endif
#endif