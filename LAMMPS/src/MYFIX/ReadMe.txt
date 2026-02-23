This is the simulation software for the manuscript.
To run the simulation, 
1. Open a terminal and go to the folder .../LAMMPS/src/
2. Compile the program:
    make yes-MYFIX  P.S.: MYFIX is the folder name of the customized LAMMPS package
    make mpi
3. Run the simulation with the example input file using the following command
        mpirun -np NUM_CORES PATH_TO_LAMMPS/LAMMPS/src/lmp_mpi -in ./in.test_3membranes_nuc_front \ 
         -var output_filename "dump.lammps"  -var log_filename "log.lammps" \
         -var seed1 123456 -var tgrad 0.8 -var nPV 300
4. variables: tgrad determines the cortical tension gradient, 
        nPV determines the number of particles of PV membrane.
        Other variables can be changed and tuned in the input script, in.test_3membranes_nuc_front, directly.