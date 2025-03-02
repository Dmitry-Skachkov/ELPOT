




  Module VASP_files
   use Parameters
   use Potential
   implicit none
  contains




  subroutine read_vasp_potential
    real(8)            :: a,b,c
    real(8)            :: dx(3),dy(3),dz(3)
    integer            :: i,j,k
    integer            :: NX,NY,NZ
    character(2)       :: Atoms1(5)
    integer            :: N1(5)
    integer            :: jj,l
    real(8)            :: Scale
    real(8)            :: x,y,z
    open(unit=2,file=name)
    read(2,*)
    read(2,*) Scale
    read(2,*) Period(1:3,1)
    read(2,*) Period(1:3,2)
    read(2,*) Period(1:3,3)
    read(2,*) Atoms1(1:5)
    read(2,*) N1(1:5)
    Nat = N1(1)+N1(2)+N1(3)+N1(4)+N1(5)
    read(2,*) 
    do i=1,Nat
     read(2,*) x,y,z
     Coord(1:3,i) = x*Period(1:3,1)+y*Period(1:3,2)+z*Period(1:3,3)
    enddo
    jj=0
    do l=1,5
     do i=1,N1(1)
      jj = jj+1
      Atoms(jj) = Atoms1(l)
     enddo
    enddo
    if(jj/=Nat) then 
     print *,'read_vasp_potential: jj/=Nat'
     stop
    endif
    read(2,*)
    read(2,*) Nxpot,Nypot,Nzpot
    dx(1:3) = Period(1:3,1)/dfloat(Nxpot)
    dy(1:3) = Period(1:3,2)/dfloat(Nypot)
    dz(1:3) = Period(1:3,3)/dfloat(Nzpot)
    Npot = Nxpot*Nypot*Nzpot
    allocate(Pot(Npot))
    allocate(Pos(3,Npot))
    do i=1,Npot-mod(Npot,5),5                             ! read potential by 5 points
     read(2,*) Pot(i),Pot(i+1),Pot(i+2),Pot(i+3),Pot(i+4)
    enddo
    if(mod(Npot,5)/=0) then                               ! read the remaining points 
     read(2,*) (Pot(i),i=Npot-mod(Npot,5)+1,Npot)
    endif
    jj = 0
    do k=1,Nzpot                                          ! calculate the real space mesh
     do j=1,Nypot
      do i=1,Nxpot
       jj = jj + 1
       if(jj > Npot) then
        print *,'read_potential: error jj > Npot'
        stop
       endif
       Pos(1:3,jj) = (i-1)*dx(1:3) + (j-1)*dy(1:3) + (k-1)*dz(1:3)
      enddo
     enddo
    enddo
  end subroutine read_vasp_potential



  end module VASP_files







