#include "fix_tdpd_cell_polarize.h"

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
#include <unordered_map>

#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;


/* ---------------------------------------------------------------------- */
constexpr double EPSILON = 1e-10;

FixTDPDCellPolarize::FixTDPDCellPolarize(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),rng(nullptr),vest(nullptr)
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
    tau_live0 = utils::numeric(FLERR, arg[iarg++], false, lmp); // lifetime of a protrusion
    tau_cd = utils::numeric(FLERR, arg[iarg++], false, lmp); // interval of the protrusion formation
    rate_birth = utils::numeric(FLERR, arg[iarg++], false, lmp); // rate of protrusion formation
    nMaxMeanPolaritySites = utils::inumeric(FLERR, arg[iarg++], false, lmp); // maximum number coexisting protrusions
    tau = utils::numeric(FLERR, arg[iarg++], false, lmp); // OU timescale
    sigma =  utils::numeric(FLERR, arg[iarg++], false, lmp); // noise strength
    // tau = utils::numeric(FLERR, arg[iarg++], false, lmp); // lifetime of a reference
    // tau1 = utils::numeric(FLERR, arg[iarg++], false, lmp); // reactivation cooldown time

    radius_ref = utils::numeric(FLERR, arg[iarg++], false, lmp); // range of reference particles
    radius = utils::numeric(FLERR, arg[iarg++], false, lmp); // raidus of source 

    nevery = utils::inumeric(FLERR, arg[iarg++], false, lmp); // Poisson process update frequency

    int j2group = group->find(arg[iarg++]);
    if (j2group < 0) error->all(FLERR, "Cannot find the group to assgin local source particles.");
    j2groupbit = group->bitmask[j2group];
    seed = utils::inumeric(FLERR, arg[iarg++], false, lmp);
    if (iarg < narg) {
      if (strcmp(arg[iarg], "bias") == 0) {
        iarg++;
        bias_flag = true; 
        K = utils::numeric(FLERR, arg[iarg++], false,
                           lmp);    // sigmoid response of the pseudopod mechanical bias
       // n = utils::numeric(FLERR, arg[iarg++], false, lmp);    // hill coeff for tau_live saturation
      } else if(strcmp(arg[iarg],"checkObstacles") == 0)
      {
        bias_flag = false;
        iarg++;
        if(strcmp(arg[iarg],"type") ==0)
        {
          // check type
          check_obstacle_flag = 1;
          iarg++;
          for(int i = iarg; i < narg; i++)
          {
            obstacle_types.push_back(utils::inumeric(FLERR, arg[i], false, lmp));
          }
          // if(comm->me == 0)
          // {
          //   cout << "Obstacle types to check: ";
          //   for(auto type : obstacle_types)
          //   {
          //     cout << type << " ";
          //   }
          //   cout << endl;
          // }

        }
        else if(strcmp(arg[iarg],"group") ==0)
        {
          // check group
          check_obstacle_flag = 2;
          iarg++;
          for(int i = iarg; i < narg; i++)
          {
            obstacle_groupbits.push_back(group->bitmask[group->find(arg[i])]);
          }
        }
      }
    }

    mode = 2; 

    p_birth = 1 - exp( - rate_birth * update->dt * nevery);
  }
  else if(strcmp(arg[iarg],"PreferredDirection") == 0)
  {
    iarg++;
    // cout<<comm->me << " : Using PreferredDirection mode to select reference particles." << endl;
    int j1group = group->find(arg[iarg++]); // cortex 
    j1groupbit = group->bitmask[j1group]; 

    px0 = utils::numeric(FLERR, arg[iarg++], false, lmp);
    py0 = utils::numeric(FLERR, arg[iarg++], false, lmp);
    pz0 = utils::numeric(FLERR, arg[iarg++], false, lmp);
    px = px0;
    py = py0;
    pz = pz0;
    direction_flag = 1;

    tau = utils::numeric(FLERR, arg[iarg++], false, lmp); // OU timescale
    sigma =  utils::numeric(FLERR, arg[iarg++], false, lmp); // noise strength

    radius_ref = utils::numeric(FLERR, arg[iarg++], false, lmp); // range of reference particles
    radius = utils::numeric(FLERR, arg[iarg++], false, lmp); // radius of source
    
    nevery = utils::inumeric(FLERR, arg[iarg++], false, lmp); 
    int j2group = group->find(arg[iarg++]);
    if (j2group < 0) error->all(FLERR, "Cannot find the group to assgin local source particles.");
    j2groupbit = group->bitmask[j2group];
    seed = utils::inumeric(FLERR, arg[iarg++], false, lmp);
    mode = 3;
  }
}

/* ---------------------------------------------------------------------- */

FixTDPDCellPolarize::~FixTDPDCellPolarize()
{
  //if (fp) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixTDPDCellPolarize::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTDPDCellPolarize::init()
{
  // initialization

  //fp = fopen("membrane_random_walk.log", "a");
  
  nprocs = comm->nprocs;
  InitializePoolTags();
  if(mode == 0)
  {
    init_update = update->ntimestep+1; // initialize the update time for the stochastic process
    next_selection = update->ntimestep+1;
    rng = new RanMars(lmp, seed);
  }
  else if(mode == 2)
  {
    // mode = 2 : choose reference particles by dedicated stochastic process
    if(check_obstacle_flag != 0)
    {
      neighbor->add_request(this, NeighConst::REQ_FULL);
    }
    rng = new RanMars(lmp, seed);

    vest = new double[3];
    vest[0] = 0.0;
    vest[1] = 0.0;
    vest[2] = 0.0;
    double **v = atom->v;
    int *mask = atom->mask;
    int ncount = 0;
    for(int i = 0; i < atom->nlocal; i++)
    {
      if(mask[i] & j1groupbit || mask[i] & groupbit)
      {
        vest[0] += v[i][0];
        vest[1] += v[i][1];
        vest[2] += v[i][2];
        ncount++;
      } 
    }
    MPI_Allreduce(MPI_IN_PLACE, vest, 3, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, &ncount, 1, MPI_INT, MPI_SUM, world);
    if(ncount > 0)
    {
      vest[0] /= ncount;
      vest[1] /= ncount;
      vest[2] /= ncount;
    }



    vector<tagint> select_one_tag = ChooseRandomParticleFromPool(1,pool_tags);
    
    polaritySites = select_one_tag;
    meanPolaritySites = polaritySites; // mean polarity sites, used to calculate the mean tag of the reference particles
    // if(comm->me == 0)
    // {
    //   cout << "Selected one tag for the polarity: " << meanPolaritySites[0] << endl;
    // }
    siteStates = vector<int>(meanPolaritySites.size(), 1); // state: 1 -> active, 2-> cd
    timer = vector<double>(meanPolaritySites.size(),tau_cd); // timer for the reprotrusion cd
    init_update = update->ntimestep+1; // initialize the update time for the stochastic process
    next_selection = init_update; // next update time for the stochastic process
    
  }
  else if(mode == 3)
  {
    rng = new RanMars(lmp, seed);
    init_update = update->ntimestep+1; // initialize the update time for the stochastic process
    next_selection = init_update; // next update time for the stochastic process
  }
  
    
}

/* ---------------------------------------------------------------------- */

void FixTDPDCellPolarize::init_list(int, NeighList *ptr)
{
  list = ptr;
}


/* ---------------------------------------------------------------------- */
void FixTDPDCellPolarize::pre_force(int /*vflag*/)
{
  if(mode == 0)
  {
    if(update->ntimestep == init_update)
    {
      InitializeOrderedTags();
    } 
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
    } 

    if(update->ntimestep == next_selection)
    {
      for(int i = 0; i < meanPolaritySites.size(); i++)
      {
  
        double p_death = compute_pdeath(meanPolaritySites[i], bias_flag);
        if(rng->uniform() < p_death && siteStates[i] == 1)
        {
          // deactivate the protrusion
          siteStates[i] = 2; // set the state to cd
        }
      }

      // set source sites near obstacles to death
      if(check_obstacle_flag != 0)
      {
        vector<tagint> tag_near_obstacles = FindExcludedTags(meanPolaritySites, check_obstacle_flag);
        for(int i = 0; i < tag_near_obstacles.size();i++)
        {
          auto it = find(meanPolaritySites.begin(), meanPolaritySites.end(), tag_near_obstacles[i]);
          if(it != meanPolaritySites.end())
          {
            int index = distance(meanPolaritySites.begin(), it);
            if(siteStates[index] == 1) 
              siteStates[index] = 2; // set the state to cd
          }
        }
      }








      
      int nActiveSites = 0;
      for(int i = 0; i < siteStates.size(); i++)
      {
        if(siteStates[i] == 1) nActiveSites++;
      }

      if(nActiveSites == 0 || (nActiveSites < nMaxMeanPolaritySites && rng->uniform() < p_birth))
      {
        // create a new protrusion
        vector<tagint> pool_tmp = orderedTags;
        for(int i = 0; i < meanPolaritySites.size(); i++)
        {
          auto it = find(pool_tmp.begin(), pool_tmp.end(), meanPolaritySites[i]);
          if(it != pool_tmp.end()) pool_tmp.erase(it); // remove the mean polarity site from the pool
        }
        
        // remove sites that are surrounding by obstacles
        if (check_obstacle_flag != 0) {
          // if(comm->me == 0)
          // {
          //   cout << "Checking obstacles for the pool of size " << pool_tmp.size() <<", "<< update->ntimestep << endl;
          // }
          excluded_tags = FindExcludedTags(pool_tmp, check_obstacle_flag);
          // if(comm->me == 0 && !excluded_tags.empty())
          // {
          //   cout << "Excluded tags: ";
          //   for(auto tag : excluded_tags)
          //   {
          //     cout << tag << " ";
          //   }
          //   cout << endl;
          // }
          // remove the excluded tags from the pool
          for (int i = 0; i < excluded_tags.size(); i++) {
            auto it = find(pool_tmp.begin(), pool_tmp.end(), excluded_tags[i]);
            if (it != pool_tmp.end()) pool_tmp.erase(it);
          }

          // 

        }

        // error->warning(FLERR,"Try to get a new protrusion");

        vector<tagint> new_mean_site = ChooseRandomParticleFromPool(1, pool_tmp);
        // error->warning(FLERR,"New protrusion selected: ", new_mean_site[0]);
        meanPolaritySites.push_back(new_mean_site[0]); // add the new mean polarity site
        polaritySites.push_back(new_mean_site[0]); // get the new polarity site
        siteStates.push_back(1); // set the state to active
        timer.push_back(tau_cd); // reset the timer for the new protrusion
      }
   


      next_selection += nevery; // next check time for the stochastic process
    }


    vector<tagint> selected_tags;

    polaritySites = MembraneMultipleSites(polaritySites,meanPolaritySites);
    
   
    for(int i = 0; i < meanPolaritySites.size();i++)
    {
      if(siteStates[i] == 1) // active site
      {
        selected_tags.push_back(polaritySites[i]);
      }
    }






    
    if(radius_ref == 0)
     ref_tags = selected_tags;
    else 
      SetReferenceParticleTags(selected_tags, radius_ref);
  }
  else if(mode == 3)
  {
    if(update->ntimestep == init_update)
    {
      InitializeOrderedTags();
    }
    // predefined direction
    if(update->ntimestep == next_selection)
    {
      meanPolaritySites.clear();
      polaritySites.clear();
      // find the furtherest cortex particle and set it as the reference particle
      // 1. compute the unwrapped COM position of the cortex particles
      // 2. project the cortex particles position in the COM frame to the predefined direction
      // 3. select the particle with the largest projection length
      vector<double> xc = GetCOMOfGroup(j1groupbit);
      tagint maxTag = GetMaxProjectionTagInGroup(j1groupbit, xc, px, py, pz);
      // update the direction by OU process
      if(tau != 0)
      {
        px = px0 + (px - px0)*exp(-nevery*update->dt/tau) + sigma*sqrt(1 - exp(-2*nevery*update->dt/tau))*rng->gaussian();
        py = py0 + (py - py0)*exp(-nevery*update->dt/tau) + sigma*sqrt(1 - exp(-2*nevery*update->dt/tau))*rng->gaussian();
        pz = pz0 + (pz - pz0)*exp(-nevery*update->dt/tau) + sigma*sqrt(1 - exp(-2*nevery*update->dt/tau))*rng->gaussian();
      }
      else 
      {
        px = px0;
        py = py0;
        pz = pz0;
      }
      // normalize the direction
      double pnorm = sqrt(px * px + py * py + pz * pz) + EPSILON;
      px /= pnorm;
      py /= pnorm;  
      pz /= pnorm;

      meanPolaritySites.push_back(maxTag);
      polaritySites.push_back(maxTag);
      next_selection += nevery;
    }

    // cout<< comm->me << " : Mean polarity site: " << meanPolaritySites[0] << endl;
    polaritySites = MembraneMultipleSites(polaritySites,meanPolaritySites);
    // cout << comm->me << " : Polarity site: " << polaritySites[0] << endl;
    if(radius_ref == 0)
     ref_tags = polaritySites;
    else 
      SetReferenceParticleTags(polaritySites, radius_ref);

    
  }

  UpdateReferenceInfo(ref_tags, xref_global, imgref_global);
  SetSourceFromNReferenceParticles();
  // cout<< comm->me << " : Number of source particles: " << ref_tags.size() << endl;
  if(mode != 2) return;
  // update timer for cd;
  for(int i = 0; i < meanPolaritySites.size(); i++)
  {
    if(siteStates[i] == 2) // cd site
    {
      timer[i] -= update->dt;
    }
  }
  
  // remove the mean polarity sites that have been cd
  // and update the polaritySites, siteStates, timer vectors accordingly
 int i = 0;
 while(true)
  {
    if(timer[i] <= 0)
    {
      meanPolaritySites.erase(meanPolaritySites.begin() + i);
      polaritySites.erase(polaritySites.begin() + i);
      siteStates.erase(siteStates.begin() + i);
      timer.erase(timer.begin() + i);
      i--; // adjust the index after erasing an element
    }
    i++;
    if(i >= meanPolaritySites.size()) break; // exit the loop when all elements
  }
}


vector<double> FixTDPDCellPolarize::GetCOMOfGroup(int targetGroupbit)
{
  vector<double> xc(3,0.0);
  double **x = atom->x;
  imageint *image = atom->image;
  int count = 0;
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & targetGroupbit) {
      double unwrap[3];
      domain->unmap(x[i], image[i], unwrap); // Unwrap the current particle position
      xc[0] += unwrap[0];
      xc[1] += unwrap[1];
      xc[2] += unwrap[2];
      count++;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, xc.data(), 3, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM, world);

  if (count > 0) {
    xc[0] /= count;
    xc[1] /= count;
    xc[2] /= count;
  } else {
    error->all(FLERR,"No particles found in the group for the COM calculation! Source: GetCOMOfGroup");
  }

  return xc;
}

tagint FixTDPDCellPolarize::GetMaxProjectionTagInGroup(int targetGroupbit,vector<double>& xc, double px, double py,double pz)
{
  tagint maxTag = -1;
  double maxProjection = -1.0;

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & targetGroupbit) {
      double unwrap[3];
      domain->unmap(atom->x[i], atom->image[i], unwrap);
      double delta_x = unwrap[0] - xc[0];
      double delta_y = unwrap[1] - xc[1];
      double delta_z = unwrap[2] - xc[2];
      double projection = delta_x * px + delta_y * py + delta_z * pz;
      if (projection > maxProjection) {
        maxProjection = projection;
        maxTag = atom->tag[i];
      }
    }
  }

  struct {
    double projection;
    int rank;
  } in,out;
  in.projection = maxProjection;
  in.rank = comm->me;
  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, world);
  MPI_Bcast(&maxTag, 1, MPI_LMP_TAGINT, out.rank, world);
  return maxTag;

}


vector<tagint> FixTDPDCellPolarize::FindExcludedTags(vector<tagint> &pool, int check_obstacles_flag)
{
  excluded_tags.clear();

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  int *mask = atom->mask;
  int *type = atom->type;
  tagint* tag = atom->tag;

  vector<tagint> excluded_tags_local; 
   for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if(mask[i] & j1groupbit)
    {
      auto it = find(pool.begin(), pool.end(), tag[i]);
      if(it == pool.end()) continue; // only check the particles in the pool
      int itype = type[i];
      int *jlist = firstneigh[i];
      int jnum = numneigh[i];
      for (int jj = 0; jj < jnum; jj++)
      {
        int j = jlist[jj];

        if(check_obstacle_flag == 1)
        {
          auto it = find(obstacle_types.begin(), obstacle_types.end(), type[j]);
          if(it != obstacle_types.end()) 
          {
            excluded_tags_local.push_back(tag[i]);
            break;
          }
        }
        else if(check_obstacle_flag == 2)
        {
          auto it = find_if(obstacle_groupbits.begin(), obstacle_groupbits.end(),
                                        [&](int groupbit) { return mask[j] & groupbit; });
          if(it != obstacle_groupbits.end()) 
          {
            excluded_tags_local.push_back(tag[i]);
            break;
          }
        }
      }
    }
   }


   int mlocal = excluded_tags_local.size();
   vector<int> data_sizes(nprocs, 0);
   vector<int> displacements(nprocs, 0);
   
   MPI_Allgather(&mlocal, 1, MPI_INT, data_sizes.data(), 1, MPI_INT, world);
   int totalSize = accumulate(data_sizes.begin(), data_sizes.end(), 0);
   partial_sum(data_sizes.begin(), data_sizes.end() - 1, displacements.begin() + 1);
   
   excluded_tags.resize(totalSize);
   MPI_Allgatherv(excluded_tags_local.data(), mlocal, MPI_LMP_TAGINT, excluded_tags.data(),
                  data_sizes.data(), displacements.data(), MPI_LMP_TAGINT, world);

  //  if(excluded_tags.size() >= pool.size())
  //  {
  //   if(comm->me == 0)
  //     cout << "All excluded tags are considered excluded" << endl;
  //   excluded_tags.clear();
  //  }
   return excluded_tags;

}

double FixTDPDCellPolarize::compute_pdeath(tagint tag, bool bias_flag)
{
  // tag: protrusive site tag
  // r: radius of the protrusive site
  double p_base = 1 - exp( - 1.0 / tau_live0 * update->dt * nevery);
  if(bias_flag == false) return p_base;


    double **x = atom->x;
    double **v = atom->v;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    imageint *image = atom->image;

    double unwrap[3],unwrapref[3];
    int i = atom->map(tag);
    int root = -1;
    if( i >= 0 && i < nlocal) {
      root = comm->me;
      domain->unmap(x[i], image[i], unwrapref); // Unwrap the position of the reference particle
    }
    MPI_Allreduce(MPI_IN_PLACE, &root, 1, MPI_INT, MPI_MAX, world); // Ensure all processes know the root
    MPI_Bcast(unwrapref, 3, MPI_DOUBLE, root, world); // Broadcast the reference position to all processes


      double vcx = 0.0,vcy = 0.0, vcz = 0.0, xc = 0.0, yc = 0.0, zc = 0.0;
      double *buff = new double[6] {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

      int ncounts = 0;
      for (int i = 0; i < nlocal; i++) {
        if ((mask[i] & groupbit) || (mask[i] & j1groupbit)) {
          domain->unmap(x[i], image[i], unwrap); // Unwrap the current particle position
          buff[0] += unwrap[0];
          buff[1] += unwrap[1];
          buff[2] += unwrap[2];
          buff[3] += v[i][0];
          buff[4] += v[i][1];
          buff[5] += v[i][2];
          ncounts++;
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, buff, 6, MPI_DOUBLE, MPI_SUM, world);
      MPI_Allreduce(MPI_IN_PLACE, &ncounts, 1, MPI_INT, MPI_SUM, world);

      if (ncounts > 0) {
        xc = buff[0] / ncounts;
        yc = buff[1] / ncounts;
        zc = buff[2] / ncounts;
        vcx = buff[3] / ncounts;
        vcy = buff[4] / ncounts;
        vcz = buff[5] / ncounts;
      } 
      delete [] buff;
      double *plr = new double[3];
      plr[0] = unwrapref[0] - xc;
      plr[1] = unwrapref[1] - yc;
      plr[2] = unwrapref[2] - zc;

      double plrnorm = MAX(sqrt(plr[0] * plr[0] + plr[1] * plr[1] + plr[2] * plr[2]),EPSILON);
      plr[0] /= plrnorm;
      plr[1] /= plrnorm;
      plr[2] /= plrnorm;

      double vcnorm = MAX(sqrt(vcx * vcx + vcy * vcy + vcz * vcz), EPSILON);

      double prj = plr[0] * vcx/vcnorm + plr[1] * vcy/vcnorm + plr[2] * vcz/vcnorm;

      delete [] plr;

      double tau_filter = 0.25 * tau_live0;
      vest[0] = (1 - update->dt * nevery / tau_filter) * vest[0] + vcx;
      vest[1] = (1 - update->dt * nevery / tau_filter) * vest[1] + vcy;
      vest[2] = (1 - update->dt * nevery / tau_filter) * vest[2] + vcz;
      double vest_norm = MAX(sqrt(vest[0] * vest[0] + vest[1] * vest[1] + vest[2] * vest[2]), EPSILON);
      double S = prj * vcnorm / vest_norm;
      double p_extra = 0.5*p_base * (1 + 1 /(1+exp(K*(S - 1))));
      return p_extra;
}


void FixTDPDCellPolarize::ChooseOneReferenceParticle()
{
  if(radius_ref == 0)
    ref_tags = ChooseRandomParticleFromPool(1,orderedTags);
  else{
    vector<tagint> selected_tags = ChooseRandomParticleFromPool(1,orderedTags);
    SetReferenceParticleTags(selected_tags, radius_ref);
  }
}



void FixTDPDCellPolarize::SetReferenceParticleTags(vector<tagint> &selected_tags, double r)
{
  if (selected_tags.empty()) return;
  if (orderedTags.empty())
    error->all(FLERR, "FixTDPDCellPolarize: orderedTags is empty; cannot assign reference neighbors.");

  int nn = static_cast<int>(r);
  if (nn < 0) nn = 0;
  const int nordered = static_cast<int>(orderedTags.size());
  nn = min(nn, nordered - 1); // avoid wrapping through the full ring

  vector<tagint> new_refs;
  new_refs.reserve(selected_tags.size() * (2 * nn + 1));

  for (auto tag : selected_tags) {
    auto it = find(orderedTags.begin(), orderedTags.end(), tag);
    if (it == orderedTags.end()) continue; // skip tags not on the membrane list
    int idx = static_cast<int>(distance(orderedTags.begin(), it));

    new_refs.push_back(tag); // include the selected site itself
    for (int k = 1; k <= nn; k++) {
      int prev_idx = (idx - k + nordered) % nordered;
      int next_idx = (idx + k) % nordered;
      new_refs.push_back(orderedTags[prev_idx]);
      new_refs.push_back(orderedTags[next_idx]);
    }
  }

  sort(new_refs.begin(), new_refs.end());
  new_refs.erase(unique(new_refs.begin(), new_refs.end()), new_refs.end());
  ref_tags.swap(new_refs);
  if(comm->me == 0)
  {
    cout << "Selected reference tags: ";
    for(auto tag : ref_tags)
    {
      cout << tag << " ";
    }
    cout << endl;
  }
}



void FixTDPDCellPolarize::ChooseReferenceParticleCCThres()
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

void FixTDPDCellPolarize::UpdateReferenceInfo(vector<tagint> &ref_tags, vector<double> &xref_global, vector<imageint> &imgref_global)
{
  if(ref_tags.empty()) return;
  xref_global.assign(3 * ref_tags.size(), 0.0);
  imgref_global.assign(ref_tags.size(), 0);

  double **x = atom->x;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  // Collect tag/pos/image for refs owned by this rank
  vector<tagint> tags_me;
  vector<double> xref_me;
  vector<imageint> imgref_me;
  tags_me.reserve(ref_tags.size());
  xref_me.reserve(3 * ref_tags.size());
  imgref_me.reserve(ref_tags.size());

  for (auto tag : ref_tags) {
    int ilocal = atom->map(tag);
    if (ilocal < 0 || ilocal >= nlocal) continue;
    tags_me.push_back(tag);
    xref_me.push_back(x[ilocal][0]);
    xref_me.push_back(x[ilocal][1]);
    xref_me.push_back(x[ilocal][2]);
    imgref_me.push_back(image[ilocal]);
  }

  int nlocal_tags = tags_me.size();
  vector<int> sizes_tags(nprocs, 0), disp_tags(nprocs, 0);
  MPI_Allgather(&nlocal_tags, 1, MPI_INT, sizes_tags.data(), 1, MPI_INT, world);
  int total_tags = accumulate(sizes_tags.begin(), sizes_tags.end(), 0);
  partial_sum(sizes_tags.begin(), sizes_tags.end() - 1, disp_tags.begin() + 1);

  vector<tagint> tags_all(total_tags);
  MPI_Allgatherv(tags_me.data(), nlocal_tags, MPI_LMP_TAGINT, tags_all.data(),
                 sizes_tags.data(), disp_tags.data(), MPI_LMP_TAGINT, world);

  vector<int> sizes_x(nprocs, 0), disp_x(nprocs, 0);
  for (int i = 0; i < nprocs; i++) sizes_x[i] = 3 * sizes_tags[i];
  partial_sum(sizes_x.begin(), sizes_x.end() - 1, disp_x.begin() + 1);
  int total_x = accumulate(sizes_x.begin(), sizes_x.end(), 0);

  vector<double> x_all(total_x);
  MPI_Allgatherv(xref_me.data(), 3 * nlocal_tags, MPI_DOUBLE, x_all.data(),
                 sizes_x.data(), disp_x.data(), MPI_DOUBLE, world);

  vector<imageint> img_all(total_tags);
  MPI_Allgatherv(imgref_me.data(), nlocal_tags, MPI_LMP_IMAGEINT, img_all.data(),
                 sizes_tags.data(), disp_tags.data(), MPI_LMP_IMAGEINT, world);

  unordered_map<tagint, int> tag_to_idx;
  tag_to_idx.reserve(total_tags * 2);
  for (int i = 0; i < total_tags; i++) tag_to_idx[tags_all[i]] = i;

  for (int j = 0; j < ref_tags.size(); j++) {
    auto it = tag_to_idx.find(ref_tags[j]);
    if (it == tag_to_idx.end()) continue;
    int idx = it->second;
    xref_global[3 * j    ] = x_all[3 * idx    ];
    xref_global[3 * j + 1] = x_all[3 * idx + 1];
    xref_global[3 * j + 2] = x_all[3 * idx + 2];
    imgref_global[j]      = img_all[idx];
  }
}


/* ---------------------------------------------------------------------- */
void FixTDPDCellPolarize::InitializePoolTags()
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

vector<tagint> FixTDPDCellPolarize::ChooseRandomParticleFromPool(int n, vector<tagint> &pool_tags)
{
  if (pool_tags.empty()) {
    error->all(FLERR, "FixTDPDCellPolarize: No particles in the pool to select from.");
  }

  int totalSize = pool_tags.size();
  if (n < 1 || n > totalSize) {
    error->all(FLERR, "FixTDPDCellPolarize: Invalid index for selecting a particle from the pool.");
  }

  vector<tagint> selected_tags; 
  selected_tags.reserve(n);
  for(int i = 0; i < n; i++)
  {
    selected_tags.push_back(pool_tags[static_cast<int>(rng->uniform() * totalSize)]);
  }

  return selected_tags;
}


void FixTDPDCellPolarize::SetSourceFromNReferenceParticles()
{
  
  if(!ref_tags.empty())
  {
    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    imageint *image = atom->image;

    double unwrap[3],unwrapref[3],unwraprefprev[3],unwraprefnext[3];
    double xref_ptr[3];

    vector<tagint> ref_tags_prev, ref_tags_next;
    vector<double> xref_tags_prev, xref_tags_next;
    vector<imageint> imgref_tags_prev, imgref_tags_next;

    for(int j = 0; j < ref_tags.size();j++)
    {
      auto it = find(orderedTags.begin(), orderedTags.end(), ref_tags[j]);
      if(it == orderedTags.end()) continue; // safety: skip if tag not found in ordered list
      int index = distance(orderedTags.begin(), it);
      int prevIndex = (index == 0) ? orderedTags.size() - 1 : index - 1;
      int nextIndex = (index == orderedTags.size() - 1) ? 0 : index + 1;
      ref_tags_prev.push_back(orderedTags[prevIndex]);
      ref_tags_next.push_back(orderedTags[nextIndex]);
    }

    UpdateReferenceInfo(ref_tags_prev, xref_tags_prev, imgref_tags_prev);
    UpdateReferenceInfo(ref_tags_next, xref_tags_next, imgref_tags_next);
    

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

          double xprev[3], xnext[3];
          xprev[0] = xref_tags_prev[3*j];
          xprev[1] = xref_tags_prev[3*j + 1];
          xprev[2] = xref_tags_prev[3*j + 2];

          xnext[0] = xref_tags_next[3*j];
          xnext[1] = xref_tags_next[3*j + 1];
          xnext[2] = xref_tags_next[3*j + 2];


          domain->unmap(xprev, imgref_tags_prev[j], unwraprefprev);
          domain->unmap(xnext, imgref_tags_next[j], unwraprefnext); 

          double tx = unwraprefnext[0] - unwraprefprev[0];
          double ty = unwraprefnext[1] - unwraprefprev[1];
          double tnorm = sqrt(tx*tx + ty*ty);
          double nx = 0.0, ny = 0.0;
          if (tnorm > EPSILON) {
            tx /= tnorm; ty /= tnorm;
            nx = -ty; ny = tx; // 90-deg CCW rotation gives inward normal for CCW membrane
          }
          // 4) side test: accept only if on inward side
          double sidedot = dx*nx + dy*ny;
          if(r2 < radius*radius && sidedot > 0.0) nsrcs++;
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


void FixTDPDCellPolarize::InitializeOrderedTags()
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
}

tagint FixTDPDCellPolarize::MembraneSingleSite(tagint tag0, tagint meantag)
{
  tagint tagout = -1;
  
  int istart = 0, imean = 0;
  for(int i = 0; i < orderedTags.size(); i++)
  {
    if(orderedTags[i] == tag0)
    {
      istart = i;
      break;
    }
  }

  for(int i = 0; i < orderedTags.size(); i++)
  {
    if(orderedTags[i] == meantag)
    {
      imean = i;
      break;
    }
  }

  double diff = (double) (istart - imean);
  if(diff > 0.5 * orderedTags.size())
    diff -= orderedTags.size();
  else if(diff < -0.5 * orderedTags.size())
    diff += orderedTags.size();


  double iend = (double) imean + exp(-update->dt/tau)* diff+ sqrt(sigma*tau * (1.0 - exp(-2.0*update->dt/tau))) * rng->gaussian();

  if(iend < 0) iend += orderedTags.size();
  if(iend >= orderedTags.size()) iend -= orderedTags.size();

  int iiend = static_cast<int>(floor(iend + 0.5)) % orderedTags.size(); // Round to nearest integer
  if(iiend < 0) iiend += orderedTags.size();
  tagout = orderedTags[iiend];












  // MPI_Allreduce(MPI_IN_PLACE, &tagout, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  // if(tagout != tag0 && comm->me == 0)
  //   cout<< "Polarity site changed from " << tag0 << " to " << tagout << endl;
  // if(tagout == -1) error->all(FLERR, "FixTDPDCellPolarize: Membrane random walk failed to find a valid tag.");
  
  return tagout;
}

vector<tagint> FixTDPDCellPolarize::MembraneMultipleSites(vector<tagint> &tags, vector<tagint> & meantags)
{
  vector<tagint> new_tags;
  new_tags.reserve(tags.size());
  for(int i = 0; i < tags.size(); i++)
  {
    tagint new_tag = MembraneSingleSite(tags[i],meantags[i]);
    new_tags.push_back(new_tag);
  }
  return new_tags;
}
