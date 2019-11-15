#! /bin/bash
#make
NMC=1e3
for((L = 40; L <= 100; L+=20)); do
  string_run="outfile_MPI_L${L}_n1e3.txt"
  echo "Starting run for L = $L"
  ./isingmpi ${string_run} ${NMC} ${L}
done