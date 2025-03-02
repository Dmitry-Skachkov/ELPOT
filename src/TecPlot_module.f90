

 Module TecPlot_module
        use Parameters
        use Potential
        implicit none
 Contains



  subroutine write_TecPlot_1D       ! write 3 pathways    X, Y, Z
   integer     :: i,j
   do j=1,3
    if(j==1) open(unit=2,file="Pot_X.dat")
    if(j==2) open(unit=2,file="Pot_Y.dat")
    if(j==3) open(unit=2,file="Pot_Z.dat")
    write(2,1) 100
    do i=1,100
     write(2,2) Pw(j,i,j),Potw(i,j)
    enddo
    close(unit=2)
   enddo
1  format('VARIABLES = "X", "Pot"'/'ZONE I=',I4,' F=POINT')
2  format(F14.7,E19.9)
  end subroutine write_TecPlot_1D



  subroutine write_TecPlot3D_box
   integer       :: i,j,k
   print *,'write_TecPlot3D_box:'
   open(unit=2,file="Pot.csv")                                    ! write data for Mathematica
   write(2,1)
   do k=1,Nm
    do j=1,Nm                                ! plot all pathways
     do i=1,Nm
      write(2,2) Posm(1:3,i,j,k),Potm(i,j,k) 
     enddo
    enddo
   enddo
   close(unit=2)
   print *,'write_TecPlot3D_box: finish'
1  format(' x, y, z, V(x,y,z)')
2  format(F5.1,',',F5.1,',',F5.1,',',F16.7)
  end subroutine write_TecPlot3D_box



  subroutine write_TecPlot3D
   integer       :: i,j,k
   print *,'write_TecPlot3D:'
   open(unit=2,file="Pot.dat")
!   write(2,1) Nxpot,Nypot,Nzpot
   write(2,2) Npot
   do i=1,Npot
    write(2,4) Pos(1:3,i),Pot(i)
   enddo
   call write_atoms(N_Mn,'Mn')
   call write_atoms(N_O, 'O ')
   call write_atoms(N_C, 'C ')
   call write_atoms(N_Cl,'Cl')
   call write_atoms(N_N, 'N ')
   call write_atoms(N_H, 'H ')
   call write_atoms(N_F, 'F ')
   do j=1,Nw                                ! plot all pathways
    write(2,*) 
    write(2,5) 'P',j,100 
    do k=1,100
     write(2,4) Pw(1:3,k,j),Potw(k,j) 
    enddo
   enddo
   close(unit=2)
   open(unit=2,file="Pathways.dat")
   do j=1,Nw                                ! plot all pathways
    write(2,*) 
    write(2,5) 'P',j,100 
    do k=1,100
     write(2,6) k+(j-1)*100,Potw(k,j)             !*(-27.211386d0) for DMol                     ! convert Ha/e to V 
    enddo
   enddo
   close(unit=2)
   print *,'write_TecPlot3D: finish'
1  format('VARIABLES = "X", "Y", "Z", "Pot"'/'ZONE I=',I4,' J=',I4,' K=',I4,' F=POINT')
2  format('VARIABLES = "X", "Y", "Z", "Pot"'/'ZONE I=',I8,' F=POINT')
4  format(3F14.7,E19.9)
5  format('VARIABLES = "I", "Pot"'/'ZONE T="',A1,I1,'" I=',I4,' F=POINT')
6  format(I5,E19.9)
  end subroutine write_TecPlot3D




  subroutine write_atoms(N,A1)
   character(*)  :: A1
   integer       :: N
   integer       :: j
   if(N > 0) then
    write(2,*)
    write(2,3) A1,N
    do j=1,Nat
     if(Atoms(j)==A1) then
      write(2,4) Coord(1:3,j),0.d0
     endif
    enddo
   endif
3  format('VARIABLES = "X", "Y", "Z", "Pot"'/'ZONE T="',A2,'" I=',I4,' F=POINT')
4  format(3F14.7,E19.9)
  end subroutine write_atoms








end Module TecPlot_module



