#!/bin/bash
#SBATCH -J jacobiMPI
#SBATCH -N 2
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --exclusive 
#SBATCH --time=01:00:00
if command -v sinfo  2>/dev/null # if on cluster
then
    module load mpi/openmpi-x86_64
    module load pmi/pmix-x86_64
    mpiprocs=( 1 2 4 8 12 16 24 32 48 52 64 80 )
    resolutions=( 125 500 2000 4000 )
    folder="out/datacluster"
    mkdir -p $folder    
else  # if on local machine
    folder="out/datalocal"
    mkdir -p $folder    
    mpiprocs=( 1 2 3 4 )
    resolutions=( 200 )
fi

END=200
STEP=4
OUTFILE=../results/residuals_serial.txt
rm $OUTFILE
touch $OUTFILE
for ((iter = 6; iter <= END; iter += STEP)); do
    printf $iter >> $OUTFILE
    ../out/build/jacobiMPI 150 $iter | grep "|residual|" >> $OUTFILE
done

OUTFILE=../results/residuals_mpi.txt
rm $OUTFILE
touch $OUTFILE
for ((iter = 6; iter <= END; iter += STEP)); do
    printf $iter >> $OUTFILE
    mpirun -np 4 ../out/build/jacobiMPI 150 $iter | grep "|residual|" >> $OUTFILE
done


OUTFILE=../results/residuals_mpi2D.txt
rm $OUTFILE
touch $OUTFILE
for ((iter = 6; iter <= END; iter += STEP)); do
    printf $iter >> $OUTFILE
    mpirun -np 4 ../out/build/jacobiMPI 2D 150 $iter | grep "|residual|" >> $OUTFILE
done

