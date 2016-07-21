!234567...............................................................
       program laplace_mpi
!234567...............................................................
       include 'mpif.h'
       parameter(id1=100)
       parameter(MAXPRC=1024)
       common /nali/itmax,iba(MAXPRC),ila(MAXPRC)
       integer*4 nprc, prc_id, ierr, rc

       write(*,*) " Starting program " 
       call MPI_INIT(ierr)
       if (ierr .ne. 0) then
          print *,'Error starting MPI program. Terminating.'
          call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
       end if

       call MPI_COMM_RANK(MPI_COMM_WORLD, prc_id, ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, nprc, ierr)
       print *, 'Number of tasks=',nprc,' My rank=',prc_id

!      if (prc_id.eq.0) then
!        write(*,*) "Reading ITMAX from file: inp"
!        open(unit=15,file='inp',status='old')
!          read(15,*) itmax
!        close(15)
!      endif
!      msg_count = 1
!      call MPI_BCAST(itmax,msg_count,MPI_INTEGER,0,
!    1                MPI_COMM_WORLD,ierr)

       itmax = 500

!      Define Input Load Distribution : 1-D decomposition
!      --------------------------------------------------
       inc = id1/nprc
       iba(1) = 1
       do i = 2,nprc+1
         iba(i) = iba(i-1)+inc
       enddo

       ila(1) = iba(1)+inc-1
       do i = 2,nprc - 1
         ila(i) = iba(i+1)-1
       enddo
       ila(nprc) = id1

       if (prc_id.eq.0) then 
          write(6,*) 'PE 0: iba and ila '
          write(6,*) (iba(i),ila(i),i=1,nprc)
       endif 
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!......Laplace equation solver ...............................
       call LAPLACE(nprc,prc_id)
!......Laplace equation solver ...............................

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       call Display_res(nprc,prc_id)

       call MPI_FINALIZE(ierr)
       stop
       end
!234567...............................................................
       SUBROUTINE LAPLACE(nprc,prc_id)

       include 'mpif.h'

       PARAMETER (ib=1,il=100)
       PARAMETER (jp=5,jl=2**jp+1)
       PARAMETER (xb=0.0,xl=0.25,yb=0.0,yl=1.0)
       DIMENSION exact(jl),comp(jl),err(jl)

       dimension s(jl)
       common /arrs1/ phi1(jl,il),x1(il),y(jl),err1(jl)

       COMMON /nalo/exact,comp,err

!    -------------------------------------------------------
!      Declarations for system configuration
       integer*4 prc_id
       parameter(MAXPRC=1024)
       common /nali/itmax,iba(MAXPRC),ila(MAXPRC)
!  ----------------------------------------------------------
       common /tarrays1/saf(jl),raf(jl)
       common /tarrays1/sab(jl),rab(jl)
!
! **** EXACT SOLUTION PHI(X,Y)=EXP(4X)*SIN(4Y)
       f(yy,xx) = EXP(4.0*xx)*SIN(4.0*yy)
!      f(YY,XX) = 1.0

!
       ib1= iba(prc_id+1)
       il1= ila(prc_id+1)

!      RELAXATION FACTOR :  RF
       rf = 1.5
       rf1 = 1.0-rf
!
!  ****** GRID DETAILS ******** 
       dx = (xl-xb)/(il-ib)
       dy = dx
       xb1 = (ib1-1)*dx
       yb1 = yb
!
!      SET UP X1 & Y GRID
       x1(ib1) = xb1
       ib1p1 = ib1+1
       do i=ib1p1,il1
         x1(i) = x1(i-1) +dx
       enddo 
       y(1) = yb1

       do j=2,jl
         y(j) = y(j-1) +dy
       enddo 
! **   END GRID DETAILS ***
!
!  ****INITIALISE ARRAYS ****
       do i=ib1,il1
         do j=1,jl
           phi1(j,i)=0
         enddo 
       enddo 
       do j=1,jl
         comp(j)=0.0
         err(j)=0.0
       enddo 
!
!  **** END INITIALISATION *****
!
! **** BOUNDARY CONDITIONS ***** 
!       ROWS
       do i=ib1,il1
          yy1 = y(1)
          xx1 = x1(i)
          yy2 = y(jl)
          xx2 = x1(i)
          phi1(1,i) = f(yy1,xx1)
          phi1(jl,i) = f(yy2,xx2)
       enddo 
!
       IF (ib1.EQ.1) THEN
!      LEFT COLUMN
         do j=1,jl
           yy1 = y(j)
           xx1 = x1(ib1)
           phi1(j,ib1) = f(yy1,xx1)
         enddo 
       END IF
!
       IF (il1.EQ.il) THEN
!      DEFINE RIGHT COLUMN
         DO j=1,jl
           yy1 = y(j)
           xx1 = x1(il1)
           phi1(j,il1) = f(yy1,xx1)
         enddo 
       END IF
!
       ie1 = (ib1+il1)/2
       il1m1 = il1-1
       jlm1 = jl-1

! ***** END BOUNDARY CONDITIONS ***
!
!
       iright =  1 
       ileft  = -1 
! ** ITERATION STARTS ****
       do ite=1,itmax
         if (mod(ite,10).eq.0) then
           WRITE (6,*) ' ITE =',ite
         endif
 
         do j=2,jlm1
           err1(j) = phi1(j,ie1)
         enddo 
!
! *****  MAIN COMPUTATIONAL LOOP DONE CONCURRENTLY ****
!
         do i=ib1p1,il1m1
           ip1 = i+1
           im1 = i-1
!
           s(1)  = phi1(1,i)
           s(jl) = phi1(jl,i)
           do j=2,jlm1
               s(j) = phi1(j,ip1) +phi1(j,im1)
           enddo 
           CALL PSOL1(s)
!
           DO  j=2,jlm1
             phi1(j,i) = rf*s(j) +rf1*phi1(j,i)
           enddo 
         enddo 
!
! **** END COMPUTATIONAL LOOP ***
!
! ***   FIND ERROR ON THE SELECTED LINE
        DO j=2,jlm1
          err1(j) = ABS(err1(j)-phi1(j,ie1))
        enddo 
        errmax = 0.0
        DO j=2,jlm1
          errmax = amax1(errmax,err1(j))
        enddo 
!       --------------------------------------
!       Interprocessor communications 

!---initialise the saf, raf, sab,rab arrays
        if (nprc.gt.1) then
          do j = 1,jl
            saf(j) = phi1(j,il1m1) 
            raf(j) = 0.0
            sab(j) =phi1(j,ib1p1) 
            rab(j) = 0.0
          enddo
          
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!          call IPCS_old (saf,raf, sab,rab, jl, jl,nprc,prc_id)
           call IPCS_new (saf,raf, sab,rab, jl, jl,nprc,prc_id)
         
          if(prc_id.ne.0) then
            do j = 1,jl
              phi1(j,ib1) = raf(j)
            end do
          endif

          if(prc_id.ne.nprc-1)then
            do j = 1,jl
              phi1(j,il1) = rab(j)
            end do
          endif
        endif

!      -----------------------------------------------------

      enddo 
! *** ITERATION END HERE
!
! *** COMPUTE THE DIFFERENCE BETWEEN EXACT & COMPUTED
        do j=1,jl
          yy1 = y(j)
          xx1 = x1(ie1)
          exact(j) = f(yy1,xx1)
          comp(j) = phi1(j,ie1)
          err(j) = ABS(exact(j)-comp(j))
        enddo 
        return
        end
!234567...............................................................
        SUBROUTINE PSOL1(s)
        PARAMETER (jp=5,jl1=2**jp,jl=2**jp+1)
        DIMENSION r(jl1),s(jl)
        a= 4.0 
        r(1) = 1.0/a
        m = 1
        DO i=2,jp
          n = m
          m = m+m
          do j=m+1,jl1,m
            s(j) = s(j-n) +s(j)*a +s(j+n)
          enddo 
          a = a*a -2.0
          r(i) = 1.0/a
        enddo 
        m = m+m
        do i=jp,1,-1
          m = m/2
          a = r(i)
          do j=m+1,jl1,m+m
            s(j) = (s(j-m) +s(j) +s(j+m))*a
          enddo 
        enddo 

        return
        END

!234567...............................................................
       subroutine Display_res(nprc,prc_id)
       include 'mpif.h'
       parameter (jp = 5, jl = 2**jp+1)
       integer*4 msg_id, prc_id
       integer istatus(MPI_STATUS_SIZE)

       parameter(MAXPRC=1024)
       common /nali/itmax,iba(MAXPRC),ila(MAXPRC)
       common /nalo/exact(jl),comp(jl),err(jl)
       common /nalo1/exact1(jl),comp1(jl),err1(jl)
       common /nalo2/exact2(jl),comp2(jl),err2(jl)

!      ---------------------------------------------------
       if (prc_id.eq.0) then 

        write(6,991) (prc_id, exact(j),comp(j),err(j),j=1,jl) 

        do i= 1, nprc-1
          msg_id = i
           msg_count = jl*3
           itag = 10
          call MPI_RECV(exact1,msg_count,MPI_REAL,msg_id,
     1       itag,MPI_COMM_WORLD,istatus,ierr)

          write(6,991) (i,exact1(j),comp1(j),err1(j),j=1,jl) 
        enddo
       else
         msg_id = 0 
         msg_count = jl*3
         itag = 10
         call MPI_SEND(exact,msg_count,MPI_REAL,msg_id,
     1                 itag,MPI_COMM_WORLD,ierr)

       endif 
!      ----------------------------------------------------------
!      Display the Results

 991   format(3x,'PE',I2,'  EXACT=',F10.4,5X,'COMP=',F10.4,5X,
     1       'ERR=', F10.4)

       return
       end
!234567...............................................................
       subroutine IPCS_old(sa,ra,sb,rb,msg_count,jdim,nprc,prc_id) 

       include 'mpif.h'
       integer*4 prc_id,nextp, ipc_start
       integer*4 istatus(MPI_STATUS_SIZE)
       dimension  sa(jdim),ra(jdim)
       dimension  sb(jdim),rb(jdim)

       ipc_start=0
       itag = 1

!--------------------------------------------------------------
         nextp=prc_id+1
         nextp=mod(nextp+nprc,nprc)
         if(prc_id .eq.ipc_start) then
           call MPI_SEND(sa,msg_count,MPI_REAL,nextp,
     1                   itag,MPI_COMM_WORLD,ierr)
           nextp = prc_id - 1
           nextp=mod(nextp+nprc,nprc)
           call MPI_RECV(ra,msg_count,MPI_REAL,nextp,
     1                   itag,MPI_COMM_WORLD,istatus,ierr)
         else 
           nextp = prc_id-1      
           nextp=mod(nextp+nprc,nprc)
           call MPI_RECV(ra,msg_count,MPI_REAL,nextp,
     1                   itag,MPI_COMM_WORLD,istatus,ierr)
           nextp = prc_id + 1
           nextp=mod(nextp+nprc,nprc)

           call MPI_SEND(sa,msg_count,MPI_REAL,nextp,
     1                   itag,MPI_COMM_WORLD,ierr)
         endif

!--------------------------------------------------------------

         if(prc_id.eq.ipc_start)then
           nextp = prc_id -1
           nextp=mod(nextp+nprc,nprc)
           call MPI_SEND(sb,msg_count,MPI_REAL,nextp, 
     1                   itag,MPI_COMM_WORLD,ierr)
           nextp = prc_id +1
           nextp=mod(nextp+nprc,nprc)
           call MPI_RECV (rb,msg_count,MPI_REAL,nextp,
     1                    itag,MPI_COMM_WORLD,istatus,ierr)
         else 
           nextp = prc_id+1
           nextp=mod(nextp+nprc,nprc)
           call MPI_RECV(rb,msg_count,MPI_REAL,nextp,
     1                   itag,MPI_COMM_WORLD,istatus,ierr)
           nextp = prc_id -1
           nextp=mod(nextp+nprc,nprc)

           call MPI_SEND(sb,msg_count,MPI_REAL,nextp,
     1        itag,MPI_COMM_WORLD,ierr)
         endif
!--------------------------------------------------------------

       return
       end
!234567...............................................................
       subroutine IPCS_new(sa,ra,sb,rb,msg_count,jdim,nprc,prc_id) 

       include 'mpif.h'
       integer*4 prc_id,nextp, ipc_start
       integer*4 istatus(MPI_STATUS_SIZE)
       dimension  sa(jdim),ra(jdim)
       dimension  sb(jdim),rb(jdim)

       itag = 1

       ileft  = prc_id-1
       if (ileft.lt.0) ileft = MPI_PROC_NULL
       iright = prc_id+1
       if (iright.eq.nprc) iright = MPI_PROC_NULL
!--------------------------------------------------------------
       call MPI_SENDRECV(sa,msg_count,MPI_REAL,iright,itag,
     1                   ra,msg_count,MPI_REAL,ileft,itag,
     2                   MPI_COMM_WORLD,istatus,ierr)

       call MPI_SENDRECV(sb,msg_count,MPI_REAL,ileft,itag,
     1                   rb,msg_count,MPI_REAL,iright,itag,
     2                   MPI_COMM_WORLD,istatus,ierr)


       return
       end
!234567...............................................................

