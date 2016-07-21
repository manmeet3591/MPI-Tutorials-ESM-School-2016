#!/bin/csh -f

#BSUB -q "cccr"
#BSUB -x


#export I_MPI_FABRICS=shm:dapl
#export OMP_NUM_THREADS=1
#export F_UFMTENDIAN=big
#export I_MPI_DEVICE=rdssm
#export I_MPI_RDMA_EAGER_THRESHOLD=65536
#export I_MPI_PIN_MODE=mpd
#export I_MPI_PIN_DOMAIN=auto

#ulimit -c unlimited
#set -xu

set OUTPUT_PATH = /iitm1/cccr/guest9/swathi/20072016 # change this to where you want the output to go
set MPPNCCOMBINE = "/iitm1/cccr/guest9/iitm-esm/src/postproc/mppncombine/mppnccombine"
set NPES = 16 # Try this on more number of PE's
set WORKDIR = $OUTPUT_PATH/tempdir
cd $WORKDIR
# Change what follows to your system template.
#
/bin/rm heat.nc*
# Change what follows to your system template.
#

cat > heat_eqn_nml <<EOF
&heat_eqn_nml
 ni = 12000,
 nj = 12000,
 del_t = 1.0e-9,
 nt_steps = 200,
 io_int=10,
 domain_layout=4,4
/
EOF
mpirun -np $NPES ./a.out >& results
exit



