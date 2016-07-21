c  
c     Numerical Analysis:
c     The Mathematics of Scientific Computing
c     D.R. Kincaid & E.W. Cheney
c     Brooks/Cole Publ., 1990
c
c     Section 4.7
c   
c     Example of preconditioned conjugate gradient method 
c     (conditioned with Jacobi method) 
c
c
c     file: pcg.f
c
      parameter (n=4)
      dimension a(n,n),r(n),b(n),v(n),z(n),x(n)
      data (a(1,i),i=1,n) /420.,210.,140.,105./
      data (a(2,i),i=1,n) /210.,140.,105.,84./
      data (a(3,i),i=1,n) /140.,105.,84.,70./
      data (a(4,i),i=1,n) /105.,84.,70.,60./
      data (b(i),i=1,n) /875.,539.,399.,319./
      data m/15/
      data (x(i),i=1,n) /0.,0.,0.,0./
      data delta /1.0e-20/
c
      print *
      print *,' Conjugate Gradient preconditioned with Jacobi matrix' 
      print *,' Section 4.5, Kincaid-Cheney'
      print *
c
      call residual(n,a,x,b,r)
      print *,' The residual vector corresponding to x(0) is:'
      print 7,r
      print *
c
      do 2 i=1,n
         z(i) = r(i)/a(i,i)
         v(i) = z(i)
 2    continue
c
      c = prod(n,z,r)
      do 6 k=1,M
         if (sqrt(prod(n,v,v)) .lt. delta) stop 
         call mult(n,a,v,z)
         t = c/prod(n,v,z)
         do 3 i=1,n
            x(i) = x(i) + t*v(i)
 3       continue
         do 4 i=1,n
            r(i) = r(i) - t*z(i)
            z(i) = r(i)/a(i,i)
 4       continue
c
         d = prod(n,z,r)
c
         do 5 i=1,n
            v(i) = z(i) + (d/c)*v(i)
 5       continue
c
         c = d
c
         print *,' At iteration',k
         print *,' Solution is:'
         print 7,x
         print *,' Residual vector is:'
         print 7,r
         print *
 6    continue
c
 7    format (3x,4(e13.6,2x))
      stop
      end
c
      function prod(n,x,y)
c
c     compute the vector product
c
      dimension x(n),y(n)
c
      sum = 0.0
      do 2 i=1,n
         sum = sum + x(i)*y(i)
 2    continue
      prod = sum
c
      return
      end
c
      subroutine residual(n,a,x,b,r)
c
c     compute the residual vector
c
      dimension a(n,n),x(n),b(n),r(n)
      call mult(n,a,x,r)
c
      do 2 i=1,n
         r(i) = b(i) - r(i)
 2    continue
c
      return
      end
c
      subroutine mult(n,a,x,y)
c
c     compute the matrix-vector product
c
      dimension a(n,n),x(n),y(n)
c
      do 3 i=1,n
         sum = 0.0
         do 2 j=1,n
            sum = sum + a(i,j)*x(j)
 2       continue
         y(i) = sum
 3    continue
c
      return
      end
