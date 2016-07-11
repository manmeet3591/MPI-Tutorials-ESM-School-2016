        program poisson_serial        
       
        
         !******************************************************************************/
         !*
         !   Purpose:
         !
         !   MAIN is the main program for POISSON_SERIAL.
         !
         !   Discussion:
         !
         !   POISSON_SERIAL is a program for solving the Poisson problem.
         !
         !   This program runs serially.  Its output is used as a benchmark for
         !   comparison with similar programs run in a parallel environment.
         !
         !   The Poisson equation
         !
         !   - DEL^2 U(X,Y) = F(X,Y)
         !
         !   is solved on the unit square [0,1] x [0,1] using a grid of NX by
         !   NX evenly spaced points.  The first and last points in each direction
         !   are boundary points.
         !
         !   The boundary conditions and F are set so that the exact solution is
         !
         !   U(x,y) = sin ( pi * x * y )
         !
         !   so that
         !
         !   - DEL^2 U(x,y) = pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )
         !
         !   The Jacobi iteration is repeatedly applied until convergence is detected.
         !
         !   For convenience in writing the discretized equations, we assume that NX = NY.
         !
         !   6 July 2016
         !
         !   Author:
         !
         !   Manmeet Singh
         !                                                        */
          integer,parameter :: NX = 11 
          integer,parameter :: NY = 11 
          integer ::  converged
          real :: diff
          real :: dx
          real :: dy
          real :: error
          real,dimension(NX,NY) :: f
          integer :: i
          integer :: it
          integer,parameter :: it_max = 1000
          integer :: j
          real :: tolerance = 0.000001
          real,dimension(NX,NY) :: u
          real :: u_norm
          real,dimension(NX,NY) :: udiff
          real,dimension(NX,NY) :: uexact
          real,dimension(NX,NY) :: unew
          real :: unew_norm
          real :: x
          real :: y
        
          dx = 1.0 / ( nx - 1 )
          dy = 1.0 / ( ny - 1 )
        !/*
        ! *   Print a message.
        ! *   */
          print *, " " 
          print *, "POISSON_SERIAL:" 
          print *, "  F90 version"
          print *, "  A program for solving the Poisson equation.\n" 
          print *, "" 
          print *, "  -DEL^2 U = F(X,Y)\n" 
          print *, " " 
          print *, "  on the rectangle 0 <= X <= 1, 0 <= Y <= 1." 
          print *, " " 
          print *, "  F(X,Y) = pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y)"

          print *, " " 
          print *, "  The number of interior X grid points is" , nx
          print *, "  The number of interior Y grid points is", ny
          print *, "  The X grid spacing is ", dx 
          print *, "  The Y grid spacing is ", dy 
        !/*
        ! *   Initialize the data.
        ! *   */
          call rhs ( nx, ny, f )
        print *,"check 1"
        !/*
        ! *   Set the initial solution estimate.
        ! *     We are "allowed" to pick up the boundary conditions exactly.
        ! *     */
          do j = 1,ny
          
            do i = 1,nx
            
              if ( (i .eq. 0) .or. (i .eq. nx - 1) .or. (j .eq. 0) .or. (j .eq. ny - 1) ) then
              
                unew(i,j) = f(i,j)
              
              else
              
                unew(i,j) = 0.0
              end if

            end do
          end do

          call r8mat_rms ( nx, ny, unew, unew_norm )
        !/*
        ! *   Set up the exact solution.
        ! *   */
          do j = 1,ny
          
            y =  float(j ) / float ( ny - 1 ) 
            do i = 1,nx
            
              x = float( i ) / float ( nx - 1 )
              call = u_exact ( x, y, uexact(i,j) )
            end do
          end do
          call r8mat_rms ( nx, ny, uexact, u_norm )
          print *, "  RMS of exact solution = ", u_norm
        !/*
        ! *   Do the iteration.
        ! *   */
          converged = 0
        
          print *," "
          print *, "  Step    ||Unew||     ||Unew-U||    ||Unew-Exact||"
          print *," " 
        
          do j = 1,ny
          
            do i = 1,nx
            
              udiff(i,j) = unew(i,j) - uexact(i,j)
            
            end do
          end do
          call r8mat_rms ( nx, ny, udiff, error )
          print *,  unew_norm, error 
        
          do it = 1,it_max
          
            do j = 1,ny
            
              do i = 1,nx
              
                u(i,j) = unew(i,j)
              
               end do

             end do

        !/*
        ! *   UNEW is derived from U by one Jacobi step.
        ! *   */
        !    sweep ( nx, ny, dx, dy, f, u, unew );
        !/*
        ! *   Check for convergence.
        ! *   */
            u_norm = unew_norm
            call r8mat_rms ( nx, ny, unew, unew_norm )
        
            do j = 1,ny
            
              do i = 1,nx
              
                udiff(i,j) = unew(i,j) - u(i,j)

              end do

            end do

            call r8mat_rms ( nx, ny, udiff, diff )
        
            do j = 1,ny
            
              do i = 1,nx
              
                udiff(i,j) = unew(i,j) - uexact(i,j)
              end do

            end do
            
            call r8mat_rms ( nx, ny, udiff, error )
        
            print *, it, unew_norm, diff, error 
        
            if ( diff .le. tolerance ) then
            
              converged = 1
              exit
            
            end if
        
        end do  
        
          if ( converged .eq. 1) then
          
            print *,"  The iteration has converged. with diff = ",diff,"&
         and tolerance = ",tolerance 
          
          else
          
            print *, "  The iteration has NOT converged."
          end if
        !/*
        ! *   Terminate.
        ! *   */
          print *, " " 
          print *, "POISSON_SERIAL:\n" 
          print *, "  Normal end of execution." 
          print *, " " 
        
        end program poisson_serial

!        List of subroutines        
!        double r8mat_rms ( int nx, int ny, double a[NX][NY] );
!        void rhs ( int nx, int ny, double f[NX][NY] );
!        void sweep ( int nx, int ny, double dx, double dy, double f[NX][NY], 
!        double u[NX][NY], double unew[NX][NY] );
!        double u_exact ( double x, double y );
!        double uxxyy_exact ( double x, double y );
 
!        /******************************************************************************/
        
        subroutine r8mat_rms ( NX, NY, a , v )
        
!        /******************************************************************************/
       ! 
       !     Purpose:
       !  
       !     R8MAT_RMS returns the RMS norm of a vector stored as a matrix.
       !  
       !     Modified:
       !  
       !     08 July 2016
       !  
       !     Author:
       !  
       !     Manmeet Singh
       !  
       !     Parameters:
       !  
       !     Input, int NX, NY, the number of rows and columns in A.
       !  
       !     Input, double A[NX][NY], the vector.
       !  
       !     Output, double R8MAT_RMS, the root mean square of the entries of A.
       !                                         */
          integer :: nx
          integer :: ny
          integer :: i
          integer :: j
          real :: v
          real,dimension(NX,NY) :: a       
 
          v = 0.0
        
          do j = 1,ny
          
            do i = 1,nx
            
              v = v + a(i,j) * a(i,j)

            end do

          end do 
         
          v = sqrt ( v / ( float  ( nx * ny )  ))
        
       end subroutine r8mat_rms 
!        /******************************************************************************/
        
        subroutine rhs ( nx, ny, f )
        
!        /******************************************************************************/
!        /*
!         *   Purpose:
!         *
!         *   RHS initializes the right hand side "vector".
!         *
!         *   Discussion:
!         *
!         *   It is convenient for us to set up RHS as a 2D array.  However, each
!         *   entry of RHS is really the right hand side of a linear system of the
!         *   form
!         *
!         *   A * U = F
!         *
!         *   In cases where U(I,J) is a boundary value, then the equation is simply
!         *
!         *   U(I,J) = F(i,j)
!         *
!         *   and F(I,J) holds the boundary data.
!         *
!         *   Otherwise, the equation has the form
!         *
!         *   (1/DX^2) * ( U(I+1,J)+U(I-1,J)+U(I,J-1)+U(I,J+1)-4*U(I,J) ) = F(I,J)
!         *
!         *   where DX is the spacing and F(I,J) is the value at X(I), Y(J) of
!         *
!         *   pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )
!         *   
!         *   Modified:
!         *
!         *   08 July 2016
!         *
!         *   Author:
!         *
!         *   Manmeet Singh
!         *
!         *   Parameters:
!         *
!         *   Input, int NX, NY, the X and Y grid dimensions.
!         *
!         *   Output, double F[NX][NY], the initialized right hand side data.
!         *                                                                                         */
       
          real,dimension(nx,ny) :: f
          integer :: nx
          integer :: ny 
          real :: fnorm
          integer :: i
          integer :: j
          real :: x
          real :: y
          print *,"Entered rhs"
!        /*
!         *   The "boundary" entries of F store the boundary values of the solution.
!         *     The "interior" entries of F store the right hand sides of the Poisson equation.
!         *     */
          do j = 1,ny
          
            y = ( ( float( j ) ) / ( float ( ny - 1 )) )
            do i = 1,nx
            
              x = ( float( i )) /(float ( nx - 1 ))
              if ( (i .eq. 0) .or. (i .eq. nx - 1) .or. (j .eq. 0) .or. (j .eq. ny - 1 ) ) then
              
                call u_exact ( x, y, f(i,j) )
              
              else
             
                call uxxyy_exact( x, y, f(i,j) ) 
                f(i,j) = - f(i,j)
              
              end if

            end do
           
          end do
            
          
          print *,"rhs check 2" 
          call r8mat_rms ( nx, ny, f , fnorm )
        
          print *, "  RMS of F =  ", fnorm
       
        end subroutine rhs
 
!        /******************************************************************************/
        
        subroutine sweep ( NX, NY, dx, dy, f, u, unew )
        
!        /******************************************************************************/
!        /*
!         *   Purpose:
!         *
!         *   SWEEP carries out one step of the Jacobi iteration.
!         *
!         *   Discussion:
!         *
!         *   Assuming DX = DY, we can approximate
!         *
!         *   - ( d/dx d/dx + d/dy d/dy ) U(X,Y) 
!         *
!         *   by
!         *
!         *   ( U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1) - 4*U(i,j) ) / dx / dy
!         *
!         *   The discretization employed below will not be correct in the general
!         *   case where DX and DY are not equal.  It's only a little more complicated
!         *   to allow DX and DY to be different, but we're not going to worry about 
!         *   that right now.
!         *
!         *   Modified:
!         *
!         *   08 July 2016
!         *
!         *   Author:
!         *
!         *   Manmeet Singh
!         *
!         *   Parameters:
!         *
!         *   Input, int NX, NY, the X and Y grid dimensions.
!         *
!         *   Input, double DX, DY, the spacing between grid points.
!         *
!         *   Input, double F[NX][NY], the right hand side data.
!         *
!         *   Input, double U[NX][NY], the previous solution estimate.
!         *
!         *   Output, double UNEW[NX][NY], the updated solution estimate.
!         *                                                                                    */
        
          integer :: i
          integer :: j
          integer :: NX
          integer :: NY
          real :: dx
          real :: dy
          real,dimension(NX,NY) :: f, u, unew
        
          do j = 1,ny
          
            do i = 1,nx
            
              if ( (i .eq. 1) .or. (j .eq. 1) .or. (i .eq. nx) .or. (j .eq. ny)) then
              
                unew(i,j) = f(i,j)
              
              else
               
                unew(i,j) = 0.25 * ( u(i-1,j) + u(i,j+1) + u(i,j-1) + u(i+1,j) + f(i,j)) * dx * dy 

              end if

            end do

          end do
        
        end subroutine sweep
!        /******************************************************************************/
       
        subroutine u_exact ( x, y, value )
        
!        /******************************************************************************/
!        /*
!         *   Purpose:
!         *
!         *   U_EXACT evaluates the exact solution.
!         *
!         *   Modified:
!         *
!         *   08 July 2016
!         *
!         *   Author:
!         *
!         *   Manmeet Singh
!         *
!         *   Parameters:
!         *
!         *   Input, double X, Y, the coordinates of a point.
!         *
!         *   Output, double U_EXACT, the value of the exact solution 
!         *   at (X,Y).
!         *                                       */
!        
          real,parameter :: pi = 3.141592653589793
          real :: value
          real :: x,y 
          value = sin ( pi * x * y )
               
        end subroutine u_exact

!        /******************************************************************************/
        
        subroutine uxxyy_exact ( x, y, value )
        
!        /******************************************************************************/
!        /*
!         *   Purpose:
!         *
!         *   UXXYY_EXACT evaluates ( d/dx d/dx + d/dy d/dy ) of the exact solution.
!         *
!         *   Modified:
!         *
!         *   08 July 2016
!         *
!         *   Author:
!         *
!         *   Manmeet Singh
!         *
!         *   Parameters:
!         *
!         *   Input, double X, Y, the coordinates of a point.
!         *
!         *   Output, double UXXYY_EXACT, the value of 
!         *   ( d/dx d/dx + d/dy d/dy ) of the exact solution at (X,Y).
!         *                                       */
        
          real,parameter :: pi = 3.141592653589793
          real :: value
          real ::x,y

          value = - pi * pi * ( x * x + y * y ) * sin ( pi * x * y )
        
        end subroutine uxxyy_exact
