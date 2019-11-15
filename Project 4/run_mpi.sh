#! /bin/bash
make
NMC=1e3
nproc=4
OPT="-n ${nproc}"
prog=mpi_ising
for((L = 40; L <= 100; L+=20)); do
  string_run="outfile_MPI_L${L}_n1e3.txt"
  echo "Starting run for L = $L"
  mpirun -n 4 mpi_ising ${string_run} ${NMC} ${L}
done
