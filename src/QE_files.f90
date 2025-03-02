




  Module QE_files                                   ! Quantum Espresso files
   use Parameters
   use Potential
   implicit none
   logical              :: L_pathway
   real(8)              :: dx(3),dy(3),dz(3)        ! steps in x,y,z for potential
   real(8)              :: X0(3)                    ! the center of the molecule
   integer              :: Nxpot1,Nypot1,Nzpot1
   integer              :: Nx1,Nx2,Ny1,Ny2,Nz1,Nz2
  contains




  subroutine read_qe_potential
    real(8)            :: a,b,c
    integer            :: i,j,k
    integer            :: NX,NY,NZ
    character(2)       :: Atoms1(5)
    integer            :: jj,l
    real(8)            :: x,y,z
    real(8)            :: Point(3)
    real(8)            :: Center(3)
    logical            :: L_check_file
    print *,'read_qe_potential:',name
    open(unit=2,file=name)
    read(2,*)
    read(2,*) 
    read(2,*) Period(1:3,1)
    read(2,*) Period(1:3,2)
    read(2,*) Period(1:3,3)
    read(2,*) 
    read(2,*) Nat
    do i=1,Nat
     read(2,1) Atoms(i),Coord(1:3,i)
     print 11,i,Atoms(i),Coord(1:3,i)
    enddo
    Center(1:2) = Coord(1:2,1)               ! Cl 1st atom
    Center(3)   = Coord(3,9)                 ! Nb 1st atom
    read(2,*)
    read(2,*)
    read(2,*)
    read(2,*) Nxpot,Nypot,Nzpot
    print *,'Nxpot=',Nxpot
    print *,'Nypot=',Nypot
    print *,'Nzpot=',Nzpot
    read(2,*) Point
    print *,'Point=',Point
    read(2,*) Period(1:3,1)
    print *,'Period(1:3,1)=',Period(1:3,1)
    read(2,*) Period(1:3,2)
    print *,'Period(1:3,2)=',Period(1:3,2)
    read(2,*) Period(1:3,3)
    print *,'Period(1:3,3)=',Period(1:3,3)
    dx(1:3) = Period(1:3,1)/dfloat(Nxpot-1)
    dy(1:3) = Period(1:3,2)/dfloat(Nypot-1)
    dz(1:3) = Period(1:3,3)/dfloat(Nzpot-1)
    print *,'dx(1:3)=',dx(1:3)
    print *,'dy(1:3)=',dy(1:3)
    print *,'dz(1:3)=',dz(1:3)
    Npot = Nxpot*Nypot*Nzpot
    print *,'Npot=',Npot
    allocate(Pot(Npot))
    allocate(Pos(3,Npot))
    read(2,6) (((Pot(i+(j-1)*Nxpot+(k-1)*Nxpot*Nypot),i=1,Nxpot),j=1,Nypot),k=1,Nzpot)
    close(unit=2)
    Point = -Center                               ! make center = (0,0,0)
    print *,'Point=',Point
    jj = 0
    do k=1,Nzpot                                          ! calculate the real space mesh
     do j=1,Nypot
      do i=1,Nxpot
       jj = jj + 1
       if(jj > Npot) then
        print *,'read_potential: error jj > Npot'
        stop
       endif
       Pos(1:3,jj) = (i-1)*dx(1:3) + (j-1)*dy(1:3) + (k-1)*dz(1:3) + Point(1:3)
      enddo
     enddo
    enddo
    print *,'Nxpot=',Nxpot
    print *,'Nypot=',Nypot
    print *,'Nzpot=',Nzpot
    print *,'Npot=',Npot
!    if(L_check_file('d.elpot0.xyz')) then                ! for non-central molecule just need to save xyz file with coordinates
!     open(unit=2,file='d.elpot0.xyz')
!     read(2,*)
!     read(2,*)
!     do i=1,Nat
!      read(2,*) Atoms(i),Coord(1:3,i)
!     enddo
!     close(unit=2)
!    endif
    print *,'read_qe_potential:  finish'
1   format(A2,4x,3F15.9)
6   format(6E14.6)
11  format(I3,2x,A2,3F17.5) 
  end subroutine read_qe_potential



  subroutine reduce_qe_potential_to_box     ! box 1x1x1 A^3   at (0,0,0)
   integer              :: i,j,k
   integer              :: Npot1
   real(8), allocatable :: Pos1(:,:),Pot1(:) 
   integer              :: jj
   print *,'reduce_qe_potential to box 1x1x1 at (0,0,0):'
   allocate(Pos1(3,Npot))
   allocate(Pot1(Npot))
   Npot1 = 0
   do j=1,Npot                                      ! reduce points
    if(abs(Pos(1,j)) < 2.1d0 .and. & 
       abs(Pos(2,j)) < 2.1d0 .and. & 
       abs(Pos(3,j)) < 2.1d0) then                   ! points inside the box 1x1x1  
     Npot1 = Npot1 + 1
     Pos1(1:3,Npot1) = Pos(1:3,j)
     Pot1(Npot1)     = Pot(j)
    endif 
   enddo
   deallocate(Pos,Pot)
   print *,'reduce the size from',Npot,' to',Npot1
   Npot = Npot1
   allocate(Pos(3,Npot))
   allocate(Pot(Npot))
   do j=1,Npot
    Pos(1:3,j) = Pos1(1:3,j)
    Pot(j)     = Pot1(j)
   enddo
   deallocate(Pos1,Pot1)
  end subroutine reduce_qe_potential_to_box






  subroutine reduce_qe_potential
   integer              :: i,j,k
   integer              :: Npot1
   real(8), allocatable :: Pos1(:,:),Pot1(:) 
   integer              :: jj
   print *,'reduce_qe_potential:'
   print *,'Mn atoms:',Mnw
   print *,'O  atoms:',O_w
   do j=1,5
    print 1,j,Atoms(j),Coord(1:3,Mnw(j))
   enddo
!   call calc_box                                    ! calculate Xmin,Xmax
!   call calc_new_N                                  ! calculate Nx1,Nx2,Ny1,Ny2,Nz1,Nz2 starting and finishing points for the box
   allocate(Pos1(3,Npot))
   allocate(Pot1(Npot))
   Npot1 = 0
 if(.true.) then
   do j=1,Npot                                      ! reduce points
    call calc_points_near_Mn(Pos(1:3,j))            ! ccheck points around Mn atoms for pathway
!    call calc_points_in_box(Pos(1:3,j))             ! check points in box 
    if(L_pathway) then
 !    if(mod(j,13)==0) then
      Npot1 = Npot1 + 1
      Pos1(1:3,Npot1) = Pos(1:3,j)
      Pot1(Npot1)     = Pot(j)
 !    endif
    endif
   enddo
 else
   print *,'reduce Pot:'
   print *,'Nx1=',Nx1
   print *,'Nx2=',Nx2
   print *,'Ny1=',Ny1
   print *,'Ny2=',Ny2
   print *,'Nz1=',Nz1
   print *,'Nz2=',Nz2
   Npot1 = 0
   do k=Nz1,Nz2
    do j=Ny1,Ny2
     do i=Nx1,Nx2
      Npot1 = Npot1 + 1
      jj = i+(j-1)*Nxpot+(k-1)*Nxpot*Nypot
      Pos1(1:3,Npot1) = Pos(1:3,jj)
      Pot1(Npot1)     = Pot(jj)
     enddo
    enddo
   enddo
   Nxpot = Nx2-(Nx1-1)      
   Nypot = Ny2-(Ny1-1)      
   Nzpot = Nz2-(Nz1-1)      
   if(Nxpot*Nypot*Nzpot/=Npot1) then
    print *,'ERROR in Npot'
    stop
   endif
 endif
   deallocate(Pos,Pot)
   print *,'reduce the size from',Npot,' to',Npot1
   Npot = Npot1
   allocate(Pos(3,Npot))
   allocate(Pot(Npot))
   do j=1,Npot
    Pos(1:3,j) = Pos1(1:3,j)
    Pot(j)     = Pot1(j)
   enddo
   deallocate(Pos1,Pot1)
 1 format(I3,2x,A2,3F17.5)
  end subroutine reduce_qe_potential






   subroutine calc_points_near_Mn(Px)
    real(8)    :: Px(3)
    integer    :: i
    real(8)    :: d
    real(8)    :: Rcut
    Rcut = 3.5d0
    L_pathway = .false.
    do i=1,5                            ! cycle over all Mn atoms for pathway
     call distance(Px,Coord(1:3,Mnw(i)),d)
     if(d < Rcut) L_pathway = .true.
    enddo
   end subroutine calc_points_near_Mn











  subroutine calc_box
   integer      :: i,k
   print *
   print *,'calc_box:'
   do i=1,Nat
     print 11,i,Atoms(i),Coord(1:3,i)
11   format(I3,2x,A2,3F17.5) 
    enddo
   X0(1:3) = 0.d0
   do k=1,3                            ! calculate the center point of the molecule
    do i=1,Nat
     X0(k) = X0(k) + Coord(k,i) 
    enddo
    X0(k) = X0(k)/dfloat(Nat)
   enddo
   print *,'X0=',X0
   Xmin(1:3) = X0(1:3)
   Xmax(1:3) = X0(1:3)
   print *,'initial Xmin=',Xmin
   print *,'initial Xmax=',Xmax
   do k=1,3                            ! calculate min and max of the molecule
    do i=1,Nat
     Xmin(k) = dmin1(Xmin(k),Coord(k,i))
     Xmax(k) = dmax1(Xmax(k),Coord(k,i))
    enddo
   enddo
   print *,'Xmin(1)=',Xmin(1)
   print *,'Xmax(1)=',Xmax(1)
   print *,'Xmin(2)=',Xmin(2)
   print *,'Xmax(2)=',Xmax(2)
   print *,'Xmin(3)=',Xmin(3)
   print *,'Xmax(3)=',Xmax(3)
   print *,'calc_box: finish'
  end subroutine calc_box




  subroutine calc_new_N                                  ! calculate new Nxpot,Nypot,Pzpot
   real(8)    :: Xmin1
   print *,'calc_new_N:'
   Xmin1 = Xmin(1) + (Xmax(2)-Xmin(2))/dsqrt(3.d0)
   print *,'dx(1)=',dx(1)
   print *,'dy(2)=',dy(2)
   print *,'dz(3)=',dz(3)
   Nx1 = 180 !Xmin1/dx(1)               ! ne sovsem pravil'no po Nx1,Nx2
   Nx2 = 250 !Xmax(1)/dx(1)
   Ny1 = 180 !Xmin(2)/dy(2)
   Ny2 = 250 !Xmax(2)/dy(2)
   Nz1 = 180 !Xmin(3)/dz(3)
   Nz2 = 250 !Xmax(3)/dz(3)
   print *,'calc_new_N: finish'
  end subroutine calc_new_N





   subroutine calc_points_in_box(Px)             
    real(8)      :: Px(3)
    L_pathway = .false.
    if((Xmin(1)<Px(1).and.Px(1)<Xmax(1)) .and. &
       (Xmin(2)<Px(2).and.Px(2)<Xmax(2)) .and. &
       (Xmin(3)<Px(3).and.Px(3)<Xmax(3))) then
     L_pathway = .true.
    endif
   end subroutine calc_points_in_box












   subroutine make_central
    call calc_Pc                        ! calculate center of the molecule
    call move_parts                     ! combine all parts of the molecule around (0,0,0)
    call calc_center                    ! test of the new center
    call move_parts                     ! combine all parts of the molecule around (0,0,0)
    call make_shift                     ! shift to the middle
   end subroutine make_central




   subroutine calc_Pc                   ! calculate center of the molecule
    integer     :: i
    real(8)    :: Pc(3)
    Pc = 0.d0
    do i=1,Nat
     if(i==1 .or. &                     ! central Mn core of the molecule
        i==2 .or. &
        i==3 .or. &
        i==4) then
      Pc(1:3) = Pc(1:3) + Coord(1:3,i)
     endif
    enddo
    Pc = Pc/4.d0 - (Period(1:3,1) + Period(1:3,2) + Period(1:3,3))/3.d0
  !  print *,'Center of the molecule=',Pc
    do i=1,Nat
     Coord(1:3,i) = Coord(1:3,i) - Pc(1:3)
  !   print 1,i,Atoms(i),Coord(1:3,i)
    enddo
! 1  format(i3,2x,A2,3F19.5)
   end subroutine calc_Pc




   subroutine move_parts                     ! combine all parts of the molecule around (0,0,0)
    use Miller
    real(8)    :: u,v,w
    integer    :: i,j
!    print *
!    print *,'Miller indexes'
   do j=1,3
    do i=1,Nat
     call miller_indexes(Coord(1:3,i),Period(1:3,1),Period(1:3,2),Period(1:3,3),u,v,w)
     if(u>0.5d0) then
      Coord(1:3,i) = Coord(1:3,i) - Period(1:3,1)
     elseif(u<-0.5d0) then
      Coord(1:3,i) = Coord(1:3,i) + Period(1:3,1)
     endif
     if(v>0.5d0) then
      Coord(1:3,i) = Coord(1:3,i) - Period(1:3,2)
     elseif(v<-0.5d0) then
      Coord(1:3,i) = Coord(1:3,i) + Period(1:3,2)
     endif
     if(w>0.5d0) then
      Coord(1:3,i) = Coord(1:3,i) - Period(1:3,3)
     elseif(w<-0.5d0) then
      Coord(1:3,i) = Coord(1:3,i) + Period(1:3,3)
     endif
    enddo
   enddo
! 1  format(3F18.5)
   end subroutine move_parts



   subroutine calc_center                    ! test of the new center
    integer    :: i
    real(8)    :: Pc(3)
    Pc = 0.d0
    do i=1,4
     Pc(1:3) = Pc(1:3) + Coord(1:3,i)
    enddo
    Pc = Pc/4.d0
    do i=1,Nat
     Coord(1:3,i) = Coord(1:3,i) - Pc(1:3)
    enddo
!    print *,'new center=',Pc
   end subroutine calc_center   




   subroutine make_shift                     ! shift to the middle
    real(8)   :: X0(3)
    integer   :: i
    X0(1:3) = (Period(1:3,1) + Period(1:3,2) + Period(1:3,3))/3.d0
    do i=1,Nat
     Coord(1:3,i) = Coord(1:3,i) + X0(1:3)
    enddo
   end subroutine make_shift            








  end module QE_files











