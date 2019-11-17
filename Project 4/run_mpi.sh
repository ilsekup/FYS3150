#! /bin/bash
make
NMC=1e7
nproc=4
OPT="-n ${nproc}"
prog=mpi_ising
for((L = 40; L <= 100; L+=20)); do
  string_run="outfile_MPI_L${L}.txt"
  echo "Starting run for L = $L"
  mpirun ${OPT} ${prog} ${string_run} ${NMC} ${L}
done
python plottingising.py
