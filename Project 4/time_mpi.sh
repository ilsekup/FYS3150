#! /bin/bash
make
NMC=1e5
nproc=1
prog=mpi_ising
n_test=100
for((nproc=1; nproc<=4; nproc++))do
    for((L = 40; L <= 100; L+=60));do
        echo "Starting run for L = $L"
        file="timeMPI_L${L}_nproc${nproc}.txt"
        for((i=1; i<=n_test; i++));do
            echo "Run number ${i} of ${n_test}"
            OPT="-n ${nproc}"
            mpirun ${OPT} ${prog} "test" ${NMC} ${L} 1 | tee --append ${file}
        done
    done
done
python plottingising.py
