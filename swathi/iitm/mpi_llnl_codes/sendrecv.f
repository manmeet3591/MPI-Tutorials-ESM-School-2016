
	IMPLICIT NONE
	INCLUDE "mpif.h"
	INTEGER a,b,myrank,nprocs,ierr
	integer istat(MPI_STATUS_SIZE)
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
	if (myrank.eq.0) then
	     a=1;b=3
  	else
     	     a=2;b=4
	endif

	if (myrank == 0) then
          call MPI_SENDRECV(b,1,MPI_REAL,1,0,
     .                     a,1,MPI_REAL,1,0,
     .                     MPI_COMM_WORLD,istat,ierr)
        elseif (myrank == 1) then
          call MPI_SENDRECV(b,1,MPI_REAL,0,0,
     .                     a,1,MPI_REAL,0,0,
     .                     MPI_COMM_WORLD,istat,ierr)
       end if
       if (myrank.eq.0) then
         write(*,*) b,a
       else
         write(*,*) a,b
       endif
       CALL MPI_FINALIZE(ierr)
       END
