#include "CellMigration2d.h"
#include "error.h"


namespace LAMMPS_NS {
namespace CellMigration2d_NS {
    CellMigration2d& CellMigration2d::getInstance(LAMMPS *lmp, int narg, char **arg) {
        static CellMigration2d instance(lmp, narg, arg);
        return instance;
    }

    CellMigration2d::CellMigration2d(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp) {
    }
}
}