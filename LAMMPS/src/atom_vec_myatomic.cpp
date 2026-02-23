/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_vec_myatomic.h"

#include "atom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecMyAtomic::AtomVecMyAtomic(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = {"age","membraneBond","state"};
  fields_copy = {"age","membraneBond","state"};
  fields_comm = {"age","membraneBond","state"};
  fields_comm_vel = {"age","membraneBond","state"};
  fields_border = {"age","membraneBond","state"};
  fields_border_vel = {"age","membraneBond","state"};
  fields_exchange = {"age","membraneBond","state"};
  fields_restart = {"age","membraneBond","state"};
  fields_create = {"age","membraneBond","state"};

  fields_data_atom = {"id", "type", "x"};
  fields_data_vel = {"id", "v"};

  atom->age_flag = 1;
  atom->state_flag = 1;
  atom->membraneBond_flag = 1;
  

  setup_fields();
}
