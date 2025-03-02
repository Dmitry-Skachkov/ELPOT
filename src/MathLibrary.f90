






    subroutine distance(X,Y,d)
     real(8)   :: X(3),Y(3),d
     d = dsqrt((X(1)-Y(1))**2+(X(2)-Y(2))**2+(X(3)-Y(3))**2)
    end subroutine distance

   




    real(8) function dlength(X)
     real(8)    :: X(3)
     dlength = dsqrt(X(1)**2+X(2)**2+X(3)**2)
    end function dlength







   subroutine create_list(N,Array,List,n_reg)         ! create list of numbers like 1,2,3,4,5
    integer        :: N
    integer        :: Array(N)
    character(*)   :: List
    integer        :: n_reg
    character(7)   :: fstr
    List = trim(adjustl(fstr(Array(1))))
    if(N>1) then
     do i=2,N
      List = trim(adjustl(List))//','//trim(adjustl(fstr(Array(i))))
     enddo
    endif
    n_reg = len(trim(adjustl(List)))
   end subroutine create_list









   subroutine sorting(NX,N)                          ! sort by bubble method
    integer NX(N)
    logical L_sort
    do 
     L_sort = .false.
     do j=1,N-1
      if(NX(j) > NX(j+1)) then
       L_sort = .true.
       N1 = NX(j)
       N2 = NX(j+1)
       NX(j) = N2
       NX(j+1) = N1
      endif
     enddo
     if(.not.L_sort) exit
    enddo
    print *
    print *,'after sorting'
    print *,(NX(j),j=1,N)
   end subroutine sorting
   
   








      subroutine Elements_Table_A(Aname,Zname)                                   ! convert the name of elements to the number
       character(*)   :: Aname
       integer        :: Zname
       if(trim(adjustl(Aname)).eq. 'Zn' ) then
         Zname = 30
       elseif(trim(adjustl(Aname)).eq. 'Ga'  ) then
         Zname = 31
       elseif(trim(adjustl(Aname)).eq. 'Ge' ) then
         Zname = 32
       elseif(trim(adjustl(Aname)).eq. 'N'  ) then
         Zname = 7
       elseif(trim(adjustl(Aname)).eq. 'O'  ) then
         Zname = 8
       elseif(trim(adjustl(Aname)).eq. 'P'  ) then
         Zname = 15
       elseif(trim(adjustl(Aname)).eq. 'E'  ) then
         Zname = 0
       elseif(trim(adjustl(Aname)).eq. 'Sn' ) then
         Zname = 50
       elseif(trim(adjustl(Aname)).eq. 'C' ) then
         Zname = 6
       elseif(trim(adjustl(Aname)).eq. 'Cl' ) then
         Zname = 17
       elseif(trim(adjustl(Aname)).eq. 'Al' ) then
         Zname = 13
       elseif(trim(adjustl(Aname)).eq. 'In' ) then
         Zname = 49
       elseif(trim(adjustl(Aname)).eq. 'Cu' ) then
         Zname = 29
       elseif(trim(adjustl(Aname)).eq. 'Mg' ) then
         Zname = 12
       elseif(trim(adjustl(Aname)).eq. 'Si' ) then
         Zname = 14
       elseif(trim(adjustl(Aname)).eq. 'As' ) then
         Zname = 33
       else
         Zname = 0
       endif
      end subroutine Elements_Table_A





      subroutine Elements_Table_Z(Zname,Aname)                                   ! convert the name of elements to the number
       character(2)   :: Aname
       integer        :: Zname
       if(Zname==30) then
         Aname = 'Zn'
       elseif(Zname==31) then
         Aname = 'Ga'
       elseif(Zname==32) then
         Aname = 'Ge'
       elseif(Zname==7) then
         Aname = 'N'
       elseif(Zname==8) then
         Aname = 'O'
       elseif(Zname==15) then
         Aname = 'P'
       elseif(Zname==0) then
         Aname = 'E'
       elseif(Zname==50) then
         Aname = 'Sn'
       elseif(Zname==6) then
         Aname = 'C'
       elseif(Zname==17) then
         Aname = 'Cl'
       elseif(Zname==13) then
         Aname = 'Al'
       elseif(Zname==49) then
         Aname = 'In'
       elseif(Zname==29) then
         Aname = 'Cu'
       elseif(Zname==12) then
         Aname = 'Mg'
       elseif(Zname==14) then
         Aname = 'Si'
       elseif(Zname==33) then
         Aname = 'As'
       elseif(Zname==0) then
         Aname = 'E '
       endif
      end subroutine Elements_Table_Z







      subroutine open_file(un,name,L_file)
       integer      :: un
       character(*) :: name
       logical      :: L_file
       if(L_check_file(trim(adjustl(name)))) then
        print *,' Open file ',name
        open(unit=un,file=trim(adjustl(name)))
        L_file = .true.
       else
        print *,'WARNING:  File ',trim(adjustl(name)),' does not exist'
        L_file = .false.
       endif
      end subroutine open_file







    logical function L_check_file(name)
     character(*)   name
     inquire(file=name,EXIST=L_check_file)
    end function L_check_file












       subroutine dist(a,b,c,xp,yp,zp,x1,y1,z1,ra,rmin)       ! calculate distance between two atoms with periodic condition
        real*8   :: a(3),b(3),c(3)
        real*8   :: xp,yp,zp
        real*8   :: x1,y1,z1
        real*8   :: ra(27)
        real*8   :: rmin
        real*8   :: P(3)
        integer  :: jjj
        integer  :: i,j,k
         rmin = 1.d100
         jjj=0
         do i=-1,1
          do j=-1,1
           do k=-1,1
            P = i*a + j*b + k*c
            jjj=jjj+1
            ra(jjj)=dsqrt((xp+P(1)-x1)**2+(yp+P(2)-y1)**2+(zp+P(3)-z1)**2)
            rmin=dmin1(rmin,ra(jjj))
           enddo
          enddo
         enddo
       end subroutine dist















     character(len=7) function fstr(k)                             !   Convert an integer to character*7
      integer, intent(in) :: k
      write (fstr,'(I7)') k
     end function fstr









     subroutine F_to_A(r,r_to_A,ndig,c_nreg)                    !   Convert an real to F10.ndig format
      character(20), intent(out) :: r_to_A
      real(8), intent(in)        :: r
      integer, intent(in)        :: ndig
      integer, intent(out)       :: c_nreg
      write (r_to_A,'(F20.<ndig>)') r
      c_nreg = len(trim(adjustl(r_to_A)))
     end subroutine F_to_A







     subroutine A_to_R(A,r_from_A)                              !   Convert an char variable (A) to real
      character(*), intent(in)   :: A
      real(8), intent(out)       :: r_from_A
      read (A,*) r_from_A
     end subroutine A_to_R







     subroutine A_to_I(A,i_from_A)                              !   Convert an char variable (A) to integer
      character(*), intent(in)   :: A
      integer, intent(out)       :: i_from_A
      read (A,*) i_from_A
     end subroutine A_to_I








       integer function c_nreg(i)
        integer    :: i
        integer    :: i2
        i2 = abs(i)
        if(        1<=i2.and.i2<=9) then
         c_nreg = 1
        elseif(   10<=i2.and.i2<=99) then
         c_nreg = 2
        elseif(  100<=i2.and.i2<=999) then
         c_nreg = 3
        elseif( 1000<=i2.and.i2<=9999) then
         c_nreg = 4
        elseif(10000<=i2.and.i2<=99999) then
         c_nreg = 5
        endif
        if(i<0) c_nreg = c_nreg+1
       end function c_nreg









      double precision function fmin(ax,bx,f,tol)
      implicit none
      double precision ax,bx,f,tol
!c
!c      an approximation  x  to the point where  f  attains a minimum  on
!c  the interval  (ax,bx)  is determined.
!c
!c
!c  input..
!c
!c  ax    left endpoint of initial interval
!c  bx    right endpoint of initial interval
!c  f     function subprogram which evaluates  f(x)  for any  x
!c        in the interval  (ax,bx)
!c  tol   desired length of the interval of uncertainty of the final
!c        result ( .ge. 0.0d0)
!c
!c
!c  output..
!c
!c  fmin  abcissa approximating the point where  f  attains a minimum
!c
!c
!c      the method used is a combination of  golden  section  search  and
!c  successive parabolic interpolation.  convergence is never much slower
!c  than  that  for  a  fibonacci search.  if  f  has a continuous second
!c  derivative which is positive at the minimum (which is not  at  ax  or
!c  bx),  then  convergence  is  superlinear, and usually of the order of
!  about  1.324....
!      the function  f  is never evaluated at two points closer together
!  than  eps*abs(fmin) + (tol/3), where eps is  approximately the square
!  root  of  the  relative  machine  precision.   if   f   is a unimodal
!  function and the computed values of   f   are  always  unimodal  when
!  separated by at least  eps*abs(x) + (tol/3), then  fmin  approximates
!  the abcissa of the global minimum of  f  on the interval  ax,bx  with
!  an error less than  3*eps*abs(fmin) + tol.  if   f   is not unimodal,
!  then fmin may approximate a local, but perhaps non-global, minimum to
!  the same accuracy.
!      this function subprogram is a slightly modified  version  of  the
!  algol  60 procedure  localmin  given in richard brent, algorithms for
!  minimization without derivatives, prentice - hall, inc. (1973).
!c
!c
      double precision  a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w
      double precision  fu,fv,fw,fx,x
      double precision  dabs,dsqrt,dsign
!c
!  c is the squared inverse of the golden ratio
!c
      c = 0.5d0*(3.d0 - dsqrt(5.0d0))
!c
!  eps is approximately the square root of the relative machine
!  precision.
!c
      eps = 1.0d00
   10 eps = eps/2.0d00
      tol1 = 1.0d0 + eps
      if (tol1 .gt. 1.0d00) go to 10
      eps = dsqrt(eps)
!c
!  initialization
!c
      a = ax
      b = bx
      v = a + c*(b - a)
      w = v
      x = v
      e = 0.0d0
      fx = f(x)
      fv = fx
      fw = fx
!c
!  main loop starts here
!c
   20 xm = 0.5d0*(a + b)
      tol1 = eps*dabs(x) + tol/3.0d0
      tol2 = 2.0d0*tol1
!c
!  check stopping criterion
!c
      if (dabs(x - xm) .le. (tol2 - 0.5d0*(b - a))) go to 90
!c
!c is golden-section necessary
!c
      if (dabs(e) .le. tol1) go to 40
!c
!  fit parabola
!c
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = 2.0d00*(q - r)
      if (q .gt. 0.0d0) p = -p
      q =  dabs(q)
      r = e
      e = d
!c
!  is parabola acceptable
!c
   30 if (dabs(p) .ge. dabs(0.5d0*q*r)) go to 40
      if (p .le. q*(a - x)) go to 40
      if (p .ge. q*(b - x)) go to 40

!  a parabolic interpolation step

      d = p/q
      u = x + d

!  f must not be evaluated too close to ax or bx

      if ((u - a) .lt. tol2) d = dsign(tol1, xm - x)
      if ((b - u) .lt. tol2) d = dsign(tol1, xm - x)
      go to 50

!  a golden-section step

   40 if (x .ge. xm) e = a - x
      if (x .lt. xm) e = b - x
      d = c*e

!  f must not be evaluated too close to x

   50 if (dabs(d) .ge. tol1) u = x + d
      if (dabs(d) .lt. tol1) u = x + dsign(tol1, d)
      fu = f(u)

!  update  a, b, v, w, and x

      if (fu .gt. fx) go to 60
      if (u .ge. x) a = x
      if (u .lt. x) b = x
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      go to 20
   60 if (u .lt. x) a = u
      if (u .ge. x) b = u
      if (fu .le. fw) go to 70
      if (w .eq. x) go to 70
      if (fu .le. fv) go to 80
      if (v .eq. x) go to 80
      if (v .eq. w) go to 80
      go to 20
   70 v = w
      fv = fw
      w = u
      fw = fu
      go to 20
   80 v = u
      fv = fu
      go to 20

!  end of main loop

   90 fmin = x
      return
      end function fmin

















FUNCTION zeroin( ax, bx, f, aerr) RESULT(fn_val)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-07-14  Time: 12:32:54
 
!-----------------------------------------------------------------------

!         FINDING A ZERO OF THE FUNCTION F(X) IN THE INTERVAL (AX,BX)

!                       ------------------------

!  INPUT...

!  F      FUNCTION SUBPROGRAM WHICH EVALUATES F(X) FOR ANY X IN THE
!         CLOSED INTERVAL (AX,BX).  IT IS ASSUMED THAT F IS CONTINUOUS,
!         AND THAT F(AX) AND F(BX) HAVE DIFFERENT SIGNS.
!  AX     LEFT ENDPOINT OF THE INTERVAL
!  BX     RIGHT ENDPOINT OF THE INTERVAL
!  AERR   THE ABSOLUTE ERROR TOLERANCE TO BE SATISFIED
!  RERR   THE RELATIVE ERROR TOLERANCE TO BE SATISFIED

!  OUTPUT...

!         ABCISSA APPROXIMATING A ZERO OF F IN THE INTERVAL (AX,BX)

!-----------------------------------------------------------------------
!  ZEROIN IS A SLIGHTLY MODIFIED TRANSLATION OF THE ALGOL PROCEDURE
!  ZERO GIVEN BY RICHARD BRENT IN ALGORITHMS FOR MINIMIZATION WITHOUT
!  DERIVATIVES, PRENTICE-HALL, INC. (1973).
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)

REAL (dp), INTENT(IN)  :: ax
REAL (dp), INTENT(IN)  :: bx
REAL (dp), INTENT(IN)  :: aerr
!REAL (dp), INTENT(IN)  :: rerr
REAL (dp)              :: rerr
REAL (dp)              :: fn_val

! EXTERNAL f
INTERFACE
  FUNCTION f(x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fn_val
  END FUNCTION f
END INTERFACE

REAL (dp)  :: a, b, c, d, e, eps, fa, fb, fc, tol, xm, p, q, r, s, atol, rtol

!  COMPUTE EPS, THE RELATIVE MACHINE PRECISION

eps = EPSILON(0.0_dp)

! INITIALIZATION

a = ax
b = bx
fa = f(a)
fb = f(b)
atol = 0.5 * aerr
rerr = aerr/(bx-ax)
rtol = MAX(0.5_dp*rerr, 2.0_dp*eps)

! BEGIN STEP

10 c = a
fc = fa
d = b - a
e = d
20 IF (ABS(fc) < ABS(fb)) THEN
  a = b
  b = c
  c = a
  fa = fb
  fb = fc
  fc = fa
END IF

! CONVERGENCE TEST

tol = rtol * MAX(ABS(b),ABS(c)) + atol
xm = 0.5 * (c-b)
IF (ABS(xm) > tol) THEN
  IF (fb /= 0.0) THEN
    
! IS BISECTION NECESSARY
    
    IF (ABS(e) >= tol) THEN
      IF (ABS(fa) > ABS(fb)) THEN
        
! IS QUADRATIC INTERPOLATION POSSIBLE
        
        IF (a == c) THEN
          
! LINEAR INTERPOLATION
          
          s = fb / fc
          p = (c-b) * s
          q = 1.0 - s
        ELSE
          
! INVERSE QUADRATIC INTERPOLATION
          
          q = fa / fc
          r = fb / fc
          s = fb / fa
          p = s * ((c-b)*q*(q-r)-(b-a)*(r-1.0))
          q = (q-1.0) * (r-1.0) * (s-1.0)
        END IF
        
! ADJUST SIGNS
        
        IF (p > 0.0) q = -q
        p = ABS(p)
        
! IS INTERPOLATION ACCEPTABLE
        
        IF (2.0*p < (3.0*xm*q-ABS(tol*q))) THEN
          IF (p < ABS(0.5*e*q)) THEN
            e = d
            d = p / q
            GO TO 30
          END IF
        END IF
      END IF
    END IF
    
! BISECTION
    
    d = xm
    e = d
    
! COMPLETE STEP
    
    30 a = b
    fa = fb
    IF (ABS(d) > tol) b = b + d
    IF (ABS(d) <= tol) b = b + SIGN(tol,xm)
    fb = f(b)
    IF (fb*(fc/ABS(fc)) > 0.0) GO TO 10
    GO TO 20
  END IF
END IF

! DONE

fn_val = b
RETURN
END FUNCTION zeroin







     real(8) function func(x)
      real(8) :: x
      func = dsin(x)
     end function func



     subroutine test_zeroin
      external func
      real(8)                          :: zeroin
      print *,'test zeroin: ',zeroin(1.5d0,4.5d0,func,1.d-9)
     end subroutine test_zeroin
