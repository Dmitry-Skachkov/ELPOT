




     Program main                                       ! convert potential (in QE, VASP, or DMol format) to TecPlot format
      use Parameters
      use TecPlot_module
      use SH_analysis
      call getarg(1,name)                               ! read name of the file
      call read_potential                               ! read electrostatic potential from file
      call calc_spline                                  ! calculate the potentials at necessary regions
      call write_TecPlot3D_box                          ! write potential in TecPlot format
      call write_TecPlot_1D                             ! write potential in TecPlot format
      call calc_stevens                                 ! calculate spherical harmonics and Stevens coefficients for potential
     end program main





     subroutine read_potential
      use Parameters
      use DMol_files
      use VASP_files
      use QE_files
      use TecPlot_module
      if(name=='LOCPOT') then                           ! read VASP potential
       PBC = .true.
       call read_vasp_potential
      elseif(name=='d.elpot0.xsf') then                 ! read QE potential
       PBC = .true.
       call read_qe_potential
       call reduce_qe_potential_to_box                  ! box 2x2x2 at (0,0,0)
      else                                              ! read DMol potential
       PBC = .false.
       call read_car(trim(adjustl(name))//'.car')                        ! read coordinates from .car file
       call read_DMol_potential(trim(adjustl(name))//'_potential.grd')   ! read potential from DMol file
      endif
     end subroutine read_potential



     subroutine calc_spline
      use Parameters
      use qshep                                         ! 3D splines
      use Potential
      implicit none
      real(8)              :: P(3,100)                  ! points for pathway
      integer, parameter   :: NR=129   !80              ! working parameters for qshep3
      integer, parameter   :: NW1=32 
      integer, parameter   :: NQ=17 
      integer              :: LCELL(NR,NR,NR)
      integer, allocatable :: LNEXT(:)
      real(8)              :: xyzmin(3),xyzdel(3)
      real(8)              :: rmax
      real(8), allocatable :: rsq(:)
      real(8), allocatable :: A1(:,:)
      integer              :: ier                       ! working parameters for qshep3
      real(8)              :: dx2,dy2,dz2
      integer              :: i,j,k
      real(8)              :: P1(3),P2(3)
      real(8)              :: Px,Py,Pz
      real(8)              :: theta,phi
      real(8)              :: df,rs
      integer              :: ir
      print *,'calc_spline:'
      print *,'Npot=',Npot
      allocate(LNEXT(Npot))
      allocate(rsq(Npot))
      allocate(A1(9,Npot))
      print *,'allocated LNEXT ',allocated(LNEXT)
      print *,'allocated rsq ',allocated(rsq)
      print *,'allocated A1 ',allocated(A1)
      print *,'allocated Pos ',allocated(Pos)
      print *,'allocated Pot ',allocated(Pot)
      call qshep3(Npot,Pos(1,1:Npot),Pos(2,1:Npot),Pos(3,1:Npot),Pot(1:Npot),NQ,NW1,NR,LCELL,LNEXT,xyzmin,xyzdel,rmax,rsq,A1,ier)
      if(ier/=0) then
       print *,'calc_barrier: error =',ier
       stop
      endif

      do j=1,3
       print *,'calculate path no.',j
       if(j==1) then                                ! path along X
        P1(1:3) = (/-1.d0, 0.d0, 0.d0 /)
        P2(1:3) = (/ 1.d0, 0.d0, 0.d0 /)
       elseif(j==2) then                            ! path along Y
        P1(1:3) = (/ 0.d0,-1.d0, 0.d0 /)
        P2(1:3) = (/ 0.d0, 1.d0, 0.d0 /)
       elseif(j==3) then                            ! path along Z
        P1(1:3) = (/ 0.d0, 0.d0,-1.d0 /)
        P2(1:3) = (/ 0.d0, 0.d0, 1.d0 /)
       endif
       call calc_pointsP1P2(P1,P2,P)                                         ! calculate 100 point on P1-P2 pathway
       do i=1,100                                                            ! calculate potential at all points
        Potw(i,j) = qs3val(P(1,i),P(2,i),P(3,i),Npot,Pos(1,1:Npot),Pos(2,1:Npot),Pos(3,1:Npot),Pot(1:Npot),NR,LCELL,LNEXT,xyzmin,xyzdel,rmax,rsq,A1)
        Pw(1:3,i,j) = P(1:3,i)
       enddo
       print *,'Potw(1:100,j)=',Potw(1:100,j)
      enddo
  
      do k=1,Nm
       do j=1,Nm
        do i=1,Nm                                 ! calculate potential in the box 2x2x2 at (0,0,0)
         Px = -2.d0 + (i-1)*0.2d0
         Py = -2.d0 + (j-1)*0.2d0
         Pz = -2.d0 + (k-1)*0.2d0
         Posm(1,i,j,k) =  Px
         Posm(2,i,j,k) =  Py
         Posm(3,i,j,k) =  Pz
         Potm(i,j,k) = qs3val(Px,Py,Pz,Npot,Pos(1,1:Npot),Pos(2,1:Npot),Pos(3,1:Npot),Pot(1:Npot),NR,LCELL,LNEXT,xyzmin,xyzdel,rmax,rsq,A1)
        enddo
       enddo
      enddo

     df = 180.d0/dfloat(ns)
     do ir=1,10                                    ! calculate potential on sphere rs on ns x ns grid
      print *,'ir=',ir
      rs = ir*0.1d0
      do j=1,ns
       do i=1,ns
        theta = (i-1)*df
        phi = (j-1)*(2.d0*df)
        Px = rs*dsin(theta)*dcos(phi)
        Py = rs*dsin(theta)*dsin(phi) 
        Pz = rs*dcos(theta)                                                       
        Pots(i,j,ir) = qs3val(Px,Py,Pz,Npot,Pos(1,1:Npot),Pos(2,1:Npot),Pos(3,1:Npot),Pot(1:Npot),NR,LCELL,LNEXT,xyzmin,xyzdel,rmax,rsq,A1)
       enddo
      enddo
     enddo 

      print *,'calc_spline: finish'
     end subroutine calc_spline





     subroutine calc_pointsP1P2(P1,P2,P)                                     ! calculate 100 points between P1 and P2
      real(8)              :: P(3,100)                                       ! points for pathway
      real(8)              :: P1(3),P2(3)
      real(8)              :: P12(3)
      real(8)              :: dv(3)
      print *,'call calc_pointsP1P2'
      P12(1:3) = P2(1:3) - P1(1:3)
      dv(1:3) = P12(1:3)/dfloat(101)
      do j=1,100
       P(1:3,j) = P1(1:3) + j*dv(1:3)
      enddo
     end subroutine calc_pointsP1P2






