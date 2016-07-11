#!/bin/ksh
#BSUB -q "cccr-res"
#BSUB -x
#BSUB -o "a.out"

export I_MPI_FABRICS=shm:dapl
mpirun -n 4 ./greetings.x
