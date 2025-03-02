



   Module Parameters
      integer, parameter   :: Nat_max = 1000
      integer, parameter   :: Nw = 6  !8                                     ! number of pathways
      integer              :: Kpw(2,Nw)
      real(8)              :: Pw(3,100,Nw)                                   ! 100 points in the pathway
      real(8)              :: Potw(100,Nw)
      real(8)              :: Period(3,3),Coord(3,Nat_max)
      logical              :: PBC
      character(20)        :: name
      character(2)         :: Atoms(Nat_max)
      integer              :: Mnw(4) ! 5)                                    ! Mn atoms for pathway
      integer              :: O_w(3) ! 4)                                    ! O  atoms for pathway
      real(8)              :: Xmin(3),Xmax(3)                                ! box size for molecule 
      integer              :: Nat
      integer              :: N_Mn,N_O,N_C,N_N,N_H,N_Cl,N_F                  ! number of atoms
      integer, parameter   :: Nm = 21                                        ! for Mathematica data file
      real(8)              :: Potm(Nm,Nm,Nm)
      real(8)              :: Posm(3,Nm,Nm,Nm)
      integer, parameter   :: lmax = 6                                       ! max number of shperical harmonics
      integer, parameter   :: ns = 2*(lmax+1)                                ! ns x ns grid on sphere
      real(8)              :: Pots(ns,ns,10)                                 ! potential on sphere for 10 radiuses 
   Contains






     subroutine set_path                          ! pathway on Mn12 molecule
! for -H, -CHCl2
!      Mnw = (/  8,  7, 12, 11, 6 /)              ! 5 Mn atoms
!      O_w = (/ 14, 15, 38, 39 /)                 ! 4 O atoms between Mn
! second way for -CHCl2
      Mnw = (/  9, 3, 4, 11 /)           
      O_w = (/  26, 13, 38 /)            
! second way for -C6H4F
      Mnw = (/  11, 4, 2, 7 /)           
      O_w = (/  12, 2, 7 /)            
! second way for -H
      Mnw = (/  9, 3, 4, 11 /)           
      O_w = (/  27, 13, 38 /)            
! for -C6H4F
!      Mnw = (/  8,  7,  6, 5, 12 /)              ! 5 Mn atoms
!      O_w = (/  8,  7,  6, 5 /)                  ! 4 O atoms between Mn
      O_w(1:4) = O_w(1:4) + N_Mn + N_H           ! calcuate serial number in xsf file (Mn,H,O,...)
      print *,'set_path: Mnw=',Mnw
      print *,'set_path: O_w=',O_w
      Kpw(1,1) = Mnw(1)                          ! pathway from Mn1 to O1
      Kpw(2,1) = O_w(1)
      Kpw(1,2) = O_w(1)                          ! pathway from O1 to Mn2
      Kpw(2,2) = Mnw(2)
      Kpw(1,3) = Mnw(2)                          ! pathway from Mn2 to O2
      Kpw(2,3) = O_w(2)
      Kpw(1,4) = O_w(2)                          ! pathway from O2 to Mn3
      Kpw(2,4) = Mnw(3)
      Kpw(1,5) = Mnw(3)                          ! pathway from Mn3 to O3
      Kpw(2,5) = O_w(3)
      Kpw(1,6) = O_w(3)                          ! pathway from O3 to Mn4
      Kpw(2,6) = Mnw(4)
!      Kpw(1,7) = Mnw(4)                          ! pathway from Mn4 to O5
!      Kpw(2,7) = O_w(4)
!      Kpw(1,8) = O_w(4)                          ! pathway from O5 to Mn6
!      Kpw(2,8) = Mnw(5)
     end subroutine set_path




  subroutine calc_natom
   print *,'calc_natom:'
   call calc_atom(N_Mn,'Mn')
   call calc_atom(N_O ,'O ')
   call calc_atom(N_C ,'C ')
   call calc_atom(N_Cl,'Cl')
   call calc_atom(N_N ,'N ')
   call calc_atom(N_H ,'H ')
   call calc_atom(N_F ,'F ')
   print *,'calc_natom: finish'
  end subroutine calc_natom




  subroutine calc_atom(N,At)
   integer        :: N
   character(2)   :: At
   integer        :: j
   N = 0
   do j=1,Nat
    if(Atoms(j)==At) then 
     N = N + 1
    endif
   enddo
   print *,'Number of ',At,' atoms is ',N
  end subroutine calc_atom




   end Module Parameters







