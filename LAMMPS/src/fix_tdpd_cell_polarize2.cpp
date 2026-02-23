#include "fix_tdpd_cell_polarize2.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "region.h"
#include "group.h"
#include "random_mars.h"
#include "update.h"
#include "domain.h"
#include "modify.h"

#include <cmath>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <vector>
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;


/* ---------------------------------------------------------------------- */
constexpr double EPSILON = 1e-10;

FixTDPDCellPolarize2::FixTDPDCellPolarize2(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),rng(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix tdpd/source", error);

  int iarg = 3;
  if(strcmp(arg[iarg], "RandomInGroup")==0)
  {
    iarg++;
    int j1group = group->find(arg[iarg++]);
    j1groupbit = group->bitmask[j1group];
    if (j1group < 0)
      error->all(FLERR, "Cannot find the group to randomly select a polarized position.");
    radius_ref = utils::numeric(FLERR, arg[iarg++], false, lmp);
    radius = utils::numeric(FLERR, arg[iarg++], false, lmp);
    int j2group = group->find(arg[iarg++]);
    if (j2group < 0) error->all(FLERR, "Cannot find the group to assgin local source particles.");
    j2groupbit = group->bitmask[j2group];
    seed = utils::inumeric(FLERR, arg[iarg++], false, lmp);
    nsteps = update->nsteps;    // by default, nsteps is the total steps of a run (not changing the reference position after initialization)
    if (iarg < narg) nsteps = utils::numeric(FLERR, arg[iarg++], false, lmp);
    mode = 0;
  }
  else if(strcmp(arg[iarg], "CCThresInGroup") == 0)
  {
    iarg++;
    int j1group = group->find(arg[iarg++]);
    j1groupbit = group->bitmask[j1group];
    if (j1group < 0)
      error->all(FLERR, "Cannot find the group to randomly select a polarized position.");
    cc_index = utils::inumeric(FLERR, arg[iarg++], false, lmp);
    cc_thres = utils::numeric(FLERR, arg[iarg++], false, lmp);
    radius = utils::numeric(FLERR, arg[iarg++], false, lmp);
    // if(comm->me == 0) cout << iarg << endl;

    int j2group = group->find(arg[iarg++]);
    if (j2group < 0) error->all(FLERR, "Cannot find the group to assgin local source particles.");
    j2groupbit = group->bitmask[j2group];

    mode = 1;

    
  }
  else if(strcmp(arg[iarg], "PreDefinedStochastic")==0)
  {
    iarg++;
    int j1group = group->find(arg[iarg++]);
    j1groupbit = group->bitmask[j1group];
    
    p = utils::numeric(FLERR, arg[iarg++], false, lmp);
    tau = utils::numeric(FLERR, arg[iarg++], false, lmp);
    noise = utils::numeric(FLERR, arg[iarg++], false, lmp);

    radius_ref = utils::numeric(FLERR, arg[iarg++], false, lmp); // range of reference particles
    radius = utils::numeric(FLERR, arg[iarg++], false, lmp); // raidus of source 

    if(p > 0.1 * update->dt/tau) 
      error->warning(FLERR, "FixTDPDCellPolarize2: the OU process is too slow");


    int j2group = group->find(arg[iarg++]);
    if (j2group < 0) error->all(FLERR, "Cannot find the group to assgin local source particles.");
    j2groupbit = group->bitmask[j2group];
    seed = utils::inumeric(FLERR, arg[iarg++], false, lmp);

    mode = 2; 
  }
}

/* ---------------------------------------------------------------------- */

FixTDPDCellPolarize2::~FixTDPDCellPolarize2()
{
}

/* ---------------------------------------------------------------------- */

int FixTDPDCellPolarize2::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTDPDCellPolarize2::init()
{
  // initialization
  
  nprocs = comm->nprocs;
  InitializePoolTags();
  if(mode == 0)
  {
    next_selection = 1;
    rng = new RanMars(lmp, seed);
  }
  else if(mode == 2)
  {
    // mode = 2 : choose reference particles by dedicated stochastic process

    rng = new RanMars(lmp, seed);

    vector<tagint> select_one_tag = ChooseRandomParticleFromPool(1);
    
    polaritySite = select_one_tag[0];
    meanPolaritySite = polaritySite; // initialize mean position
    init_update = update->ntimestep+1; // initialize the update time for the stochastic process
    
  }
  
    
}

/* ---------------------------------------------------------------------- */

void FixTDPDCellPolarize2::pre_force(int /*vflag*/)
{
  if(mode == 0)
  {
    if(update->ntimestep == next_selection)
    {
      ChooseOneReferenceParticle();
      next_selection += nsteps;
    }

  }
  else if(mode == 1)
  {
    ChooseReferenceParticleCCThres();
  }
  else if (mode == 2)
  {
    // mode = 2 : choose reference particles by dedicated stochastic process
    if(update->ntimestep == init_update)
    {
      InitializeOrderedTags();
      // ning(FLERR, "FINISHED INIT");
    }

    meanPolaritySite = MembraneRandomWalk(meanPolaritySite,p);

    polaritySite = MembraneOrnsteinUhlenbeck(meanPolaritySite, polaritySite, tau, noise);


    // if(comm->me ==0) cout << polaritySite << " tag now; " << meanPolaritySite << " is tag ref" << endl;

    vector<tagint> selected_tags;
    selected_tags.push_back(polaritySite);
    if(radius_ref == 0)
     ref_tags = selected_tags;
    else 
      SetReferenceParticleTags(selected_tags, radius_ref);

    if(ref_tags.empty()) error->all(FLERR, "No reference particles found in the selected region.");
    
    
  }
  
  UpdateReferenceInfo(ref_tags, xref_global, imgref_global);

  SetSourceFromNReferenceParticles();

}

void FixTDPDCellPolarize2::ChooseOneReferenceParticle()
{
  if(radius_ref == 0)
    ref_tags = ChooseRandomParticleFromPool(1);
  else{
    vector<tagint> selected_tags = ChooseRandomParticleFromPool(1);
    SetReferenceParticleTags(selected_tags, radius_ref);
  }
}



void FixTDPDCellPolarize2::SetReferenceParticleTags(vector<tagint> &selected_tags, double r)
{
  if(selected_tags.empty()) return;

  vector<double> selected_xref;
  vector<imageint> selected_imgref;
  UpdateReferenceInfo(selected_tags, selected_xref, selected_imgref);
 
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  imageint *image = atom->image;

  double unwrap[3], unwrapref[3], xref[3];

  vector<tagint> ref_tags_me;
  for(int i = 0; i < selected_tags.size(); i++)
  {
    int ilocal = atom->map(selected_tags[i]);
    if (ilocal < 0 || ilocal >= nlocal + nghost) continue; // Ensure ilocal is valid
    xref[0] = selected_xref[i*3];
    xref[1] = selected_xref[i*3 + 1];
    xref[2] = selected_xref[i*3 + 2];
    imageint imgref = selected_imgref[i];
    domain->unmap(xref, imgref, unwrapref); // Unwrap the reference position
    for(int j = 0; j < nlocal+nghost; j++)
    {
      if(mask[j] & j1groupbit)
      {
        domain->unmap(x[j], image[j], unwrap); // Unwrap the current particle position
        double dx = unwrap[0] - unwrapref[0];
        double dy = unwrap[1] - unwrapref[1];
        double dz = unwrap[2] - unwrapref[2];
        if(dx*dx + dy*dy + dz*dz < r * r)
        {
          ref_tags_me.push_back(atom->tag[j]); // Assign the tag of the reference particle
        }
      }
    }
  }

  int mlocal = ref_tags_me.size();
  vector<int> data_sizes(nprocs, 0);
  vector<int> displacements(nprocs, 0);
  MPI_Allgather(&mlocal, 1, MPI_INT, data_sizes.data(), 1, MPI_INT, world);
  int totalSize = accumulate(data_sizes.begin(), data_sizes.end(), 0);
  partial_sum(data_sizes.begin(), data_sizes.end() - 1, displacements.begin() + 1);
  ref_tags.resize(totalSize);
  MPI_Allgatherv(ref_tags_me.data(), mlocal, MPI_LMP_TAGINT, ref_tags.data(), data_sizes.data(),
                 displacements.data(), MPI_LMP_TAGINT, world);

      // remove duplicate tags
  sort(ref_tags.begin(), ref_tags.end());
  ref_tags.erase(unique(ref_tags.begin(), ref_tags.end()), ref_tags.end());
}



void FixTDPDCellPolarize2::ChooseReferenceParticleCCThres()
{
  double **x = atom->x;
  int *mask = atom->mask;
  double **cc = atom->cc;
  int nlocal = atom->nlocal;
  int k = cc_index - 1; // cc_index is 1-based, so we need to subtract 1
  vector<tagint> ref_tags_me;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & j1groupbit) {    // j1group is the cortex group
        double cc_value = cc[i][k];
        if (cc_value >= cc_thres) {
          ref_tags_me.push_back(atom->tag[i]);
        }
      }
    }
  int mlocal = ref_tags_me.size();
  vector<int> data_sizes(nprocs, 0);
  vector<int> displacements(nprocs, 0);
  MPI_Allgather(&mlocal, 1, MPI_INT, data_sizes.data(), 1, MPI_INT, world);
  int totalSize = accumulate(data_sizes.begin(), data_sizes.end(), 0);
  partial_sum(data_sizes.begin(), data_sizes.end() - 1, displacements.begin() + 1); 
  ref_tags.resize(totalSize);
  MPI_Allgatherv(ref_tags_me.data(), mlocal, MPI_LMP_TAGINT, ref_tags.data(), data_sizes.data(),
                 displacements.data(), MPI_LMP_TAGINT, world);
}

void FixTDPDCellPolarize2::UpdateReferenceInfo(vector<tagint> &ref_tags, vector<double> &xref_global, vector<imageint> &imgref_global)
{
  if(ref_tags.empty()) return;
  xref_global.clear();
  imgref_global.clear();

  double **x = atom->x;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  vector<double> xref_me;
  vector<imageint> imgref_me;
  for (int i = 0; i < ref_tags.size(); i++) {
    int ilocal = atom->map(ref_tags[i]);
    if (ilocal < 0 || ilocal >= nlocal) continue;
    xref_me.push_back(x[ilocal][0]);
    xref_me.push_back(x[ilocal][1]);
    xref_me.push_back(x[ilocal][2]);
    imgref_me.push_back(image[ilocal]);
  }
  int mlocal = xref_me.size();
  vector<int> data_sizes(nprocs, 0);
  vector<int> displacements(nprocs, 0);
  MPI_Allgather(&mlocal, 1, MPI_INT, data_sizes.data(), 1, MPI_INT, world);
  int totalSize = accumulate(data_sizes.begin(), data_sizes.end(), 0);
  partial_sum(data_sizes.begin(), data_sizes.end() - 1, displacements.begin() + 1);
  xref_global.resize(totalSize);
  MPI_Allgatherv(xref_me.data(), mlocal, MPI_DOUBLE, xref_global.data(), data_sizes.data(),
                 displacements.data(), MPI_DOUBLE, world);

  mlocal = imgref_me.size();
  data_sizes.resize(nprocs, 0);
  displacements.resize(nprocs, 0);
  MPI_Allgather(&mlocal, 1, MPI_INT, data_sizes.data(), 1, MPI_INT, world);
  totalSize = accumulate(data_sizes.begin(), data_sizes.end(), 0);
  partial_sum(data_sizes.begin(), data_sizes.end() - 1, displacements.begin() + 1);
  imgref_global.resize(totalSize);
  MPI_Allgatherv(imgref_me.data(), mlocal, MPI_LMP_IMAGEINT, imgref_global.data(),
                 data_sizes.data(), displacements.data(), MPI_LMP_IMAGEINT, world);
  

}


/* ---------------------------------------------------------------------- */
void FixTDPDCellPolarize2::InitializePoolTags()
{
  pool_tags.clear();
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  vector<tagint> pool_tags_me;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & j1groupbit) { // j2group is the source group
      pool_tags_me.push_back(atom->tag[i]);
    }
  }
  int mlocal = pool_tags_me.size();
  vector<int> data_sizes(nprocs, 0);
  vector<int> displacements(nprocs, 0);

  MPI_Allgather(&mlocal, 1, MPI_INT, data_sizes.data(), 1, MPI_INT, world);
  int totalSize = accumulate(data_sizes.begin(), data_sizes.end(), 0);
  partial_sum(data_sizes.begin(), data_sizes.end() - 1, displacements.begin() + 1);
  pool_tags.resize(totalSize);
  MPI_Allgatherv(pool_tags_me.data(), mlocal, MPI_LMP_TAGINT, pool_tags.data(), data_sizes.data(),
                 displacements.data(), MPI_LMP_TAGINT, world);

}

vector<tagint> FixTDPDCellPolarize2::ChooseRandomParticleFromPool(int n)
{
  if (pool_tags.empty()) {
    error->all(FLERR, "FixTDPDCellPolarize2: No particles in the pool to select from.");
  }

  int totalSize = pool_tags.size();
  if (n < 1 || n > totalSize) {
    error->all(FLERR, "FixTDPDCellPolarize2: Invalid index for selecting a particle from the pool.");
  }

  vector<tagint> selected_tags; 
  selected_tags.reserve(n);
  for(int i = 0; i < n; i++)
  {
    selected_tags.push_back(pool_tags[static_cast<int>(rng->uniform() * totalSize)]);
  }

  return selected_tags;
}


void FixTDPDCellPolarize2::SetSourceFromNReferenceParticles()
{
  
  if(!ref_tags.empty())
  {
    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    imageint *image = atom->image;

    double unwrap[3],unwrapref[3];
    double xref_ptr[3];
    for(int i = 0; i < nlocal; i++)
    {
      if(mask[i] & groupbit)
      {
        int nsrcs = 0;
        domain->unmap(x[i], image[i], unwrap);
        for(int j = 0; j < ref_tags.size(); j++)
        {
          
          xref_ptr[0] = xref_global[3*j];
          xref_ptr[1] = xref_global[3*j + 1];
          xref_ptr[2] = xref_global[3*j + 2];
          domain->unmap(xref_ptr, imgref_global[j],unwrapref);

          double dx = unwrap[0] - unwrapref[0];
          double dy = unwrap[1] - unwrapref[1];
          double dz = unwrap[2] - unwrapref[2];
          double r2 = dx * dx + dy * dy + dz * dz;
          if(r2 < radius * radius) nsrcs++;
        }

        if(nsrcs > 0)
        {
          // Assign these particles to the source group
          if(!(mask[i] & j2groupbit)) 
          {
             mask[i] |= j2groupbit;
          }
        }
        else
        {
          // Remove these particles from the source group
          if(mask[i] & j2groupbit) 
          {
            mask[i] &= ~j2groupbit;
          }
        }
      }
    }
    
  }
}


void FixTDPDCellPolarize2::InitializeOrderedTags()
{
  orderedTags.clear();
  int nlocal = atom->nlocal;
  
  tagint current_tag = pool_tags[0]; 

  while(!pool_tags.empty())
  {
    orderedTags.push_back(current_tag);
    pool_tags.erase(remove(pool_tags.begin(), pool_tags.end(), current_tag), pool_tags.end());

    // Find the next tag in the counterclockwise direction
    int ilocal = atom->map(current_tag);
 
    // Get the connected tags
    tagint next_tag = -1;
    if(ilocal >= 0 && ilocal < nlocal)
    {
      next_tag = atom->membraneBond[ilocal][1]; // Assume the first connected tag is in the counterclockwise direction
    }
    MPI_Allreduce(MPI_IN_PLACE, &next_tag, 1, MPI_LMP_TAGINT, MPI_MAX, world);
    
    current_tag = next_tag;
  }
  // if (comm->me == 0) {
  //   for (int i = 0; i < orderedTags.size(); i++) 
  //   { 
  //     cout << "Ordered tag " << i << ": " << orderedTags[i] << endl;
  //   }
  // }
  // error->all(FLERR, "Force STOP");

}

tagint FixTDPDCellPolarize2::MembraneRandomWalk(tagint tag0, double p)
{
  int ilocal = atom->map(tag0);
  int nlocal = atom->nlocal;
  tagint tagout = -1;
  if (ilocal >= 0 && ilocal < nlocal)
  {
    // apply random walk
    double r1 = rng->uniform();
    if(r1 < p)
    {
      tagout = atom->membraneBond[ilocal][0];
    }
    else if( r1 >= p && r1 < 2*p)
    {
      tagout = atom->membraneBond[ilocal][1];
    }
    else // r1 > 2* p && r1 < 1 - 2*p
    {
      tagout = tag0;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &tagout, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  // if(tagout != tag0 && comm->me == 0)
  // cout<< "Polarity site changed from " << tag0 << " to " << tagout << endl;
  if(tagout == -1) error->all(FLERR, "FixTDPDCellPolarize2: Membrane random walk failed to find a valid tag.");
  
  return tagout;
}

tagint FixTDPDCellPolarize2::MembraneOrnsteinUhlenbeck(tagint tag0, tagint tagnow, double tau, double noise)
{
  return tag0;
}





