#ifndef LMP_CELLMIGRATION2D_H
#define LMP_CELLMIGRATION2D_H

#include "pointers.h"
#include <vector>

namespace LAMMPS_NS {
namespace CellMigration2d_NS {

    enum CortexTensionType{
        PREDEFINED = 0,
        DYNAMIC
    };



class CellMigration2d : protected Pointers {
 public: 
 static CellMigration2d& getInstance(class LAMMPS *,int, char**);
    ~CellMigration2d() override;

    CellMigration2d(const CellMigration2d &) = delete;
    void operator=(const CellMigration2d &) = delete;

   

protected:
    CellMigration2d(class LAMMPS *, int, char **);
    
private:
     int polymerization_type, cortex_tension_type;

};
}    // namespace cellmigration2d
}    // namespace lammps_ns
#endif