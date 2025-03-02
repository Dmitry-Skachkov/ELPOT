




     Program pot_to_tecplot                             ! convert potential (in QE, VASP, or DMol format) to TecPlot format
      use Parameters
      use TecPlot_module
      call getarg(1,name)                               ! read name of the file
      call read_potential                               ! read electrostatic potential from file
      call calc_barrier                                 ! calculate the potentials at the pathway points
      call write_TecPlot_1D                             ! write potential in TecPlot format
      call write_TecPlot3D_box                             ! write potential in TecPlot format
     end program pot_to_tecplot





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
 !      call calc_natom                                  ! calculate atoms
 !      call set_path                                    ! set pathway for electron transfer
 !      call reduce_qe_potential
       call reduce_qe_potential_to_box                  ! box 1x1x1 at (0,0,0)
!      call write_tecplot3D                             ! TEST: write potential in TecPlot format
!      stop
      else                                              ! read DMol potential
       PBC = .false.
       call read_car(trim(adjustl(name))//'.car')                        ! read coordinates from .car file
       call read_DMol_potential(trim(adjustl(name))//'_potential.grd')   ! read potential from DMol file
      endif
     end subroutine read_potential




     subroutine calc_V_3D_box
      use Parameters
      use MyMathLibrary
      use Potential
      implicit none
      real(8)              :: Px,Py,Pz
      integer, parameter   :: NR=129   !80                     ! working parameters for qshep3
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
      print *,'calc_V_3D_box:'
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
      print *,'calc_barrier: qshep3'
      if(ier/=0) then
       print *,'calc_barrier: error =',ier
       stop
      endif
      do k=1,Nm
       do j=1,Nm                                                                 ! do all Nw pathways
        do i=1,Nm                                                            ! calculate potential at all points
         Px = (i-1)*0.1d0
         Py = (j-1)*0.1d0
         Pz = (k-1)*0.1d0
         Posm(1,i,j,k) =  Px
         Posm(2,i,j,k) =  Py
         Posm(3,i,j,k) =  Pz
         Potm(i,j,k) = qs3val(Px,Py,Pz,Npot,Pos(1,1:Npot),Pos(2,1:Npot),Pos(3,1:Npot),Pot(1:Npot),NR,LCELL,LNEXT,xyzmin,xyzdel,rmax,rsq,A1)
        enddo
       enddo
      enddo
      print *,'calc_V_3D_box: finish'
     end subroutine calc_V_3D_box




     subroutine calc_barrier
      use Parameters
      use MyMathLibrary
      use Potential
      real(8)              :: P(3,100)                  ! points for pathway
      integer, parameter   :: NR=129   !80                     ! working parameters for qshep3
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
      print *,'calc_barrier:'
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
      print *,'calc_barrier: qshep3'
      if(ier/=0) then
       print *,'calc_barrier: error =',ier
       stop
      endif
      do j=1,3   !Nw                                                              ! do all Nw pathways
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
      print *,'calc_barrier: finish'
     end subroutine calc_barrier








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
















