

 Module DMol_files
        use Parameters
        use Potential
        implicit none
        character(6)             :: frag(Nat_max)
        integer                  :: NZn(Nat_max)
        integer                  :: Nt
 Contains





   subroutine read_car(name)
	character(*), intent(in)    :: name
	character(1)                :: A1,A2
        character(2)                :: Ar     
	real(8)                     :: R(3)
        integer                     :: i
 	print *,'  read:  name=',name
	open(unit=2,file=name)
      read(2,*)
      read(2,*)
      read(2,*)
      read(2,*)
      if(PBC) read(2,*)
      ! calculate number of atoms
 	  Nat = 0
      do
       read(2,'(2A1)') A1,A2
       print *,'read_car: ',A1,A2,' Nat=',Nat
       if(A1=='e'.and.A2=='n') exit                ! we have found the end of .car file 
       Nat = Nat + 1
	  enddo
      rewind(2)
	  print *,' number of atoms =',Nat 
      read(2,*)
      read(2,*)
      read(2,*)
      read(2,*)
      if(PBC) read(2,920) Period(1,1),Period(2,2),Period(3,3)
      if(PBC) print 920,Period
      do i=1,Nat
       read(2,11) Ar,Coord(1:3,i)           ! read coordinates and name of atom
       if(Ar(2:2)=='1'.or.  &
          Ar(2:2)=='2'.or.  &
          Ar(2:2)=='3'.or.  &
          Ar(2:2)=='4'.or.  &
          Ar(2:2)=='5'.or.  &
          Ar(2:2)=='6'.or.  &
          Ar(2:2)=='7'.or.  &
          Ar(2:2)=='8'.or.  &
          Ar(2:2)=='9') then
        Atoms(i)=Ar(1:1)
       else
        Atoms(i) = Ar
       endif
	   print *,' Atoms=',Atoms(i)
	   print *,' Coord=',Coord(1:3,i)
      enddo
 11   format(A2,4x,3F15.9)
 920  format('PBC',3F10.4)
    end subroutine read_car








   subroutine write_car(name,comment)
      character(*),intent(in)  :: name
      character(*),intent(in)  :: comment
      character*3              :: A3            
      character*5              :: A5            
      character*7              :: fstr
      integer                  :: j
      if(Nat>999) then
       print *,'  error write_car: N>999'
       return
      endif
      open(64,file=name,status='unknown',form='formatted')
        write(64,905) 
        if(PBC) then
         write(64,904)
        else
         write(64,906)
        endif 
        write(64,915) comment,Nat
        write(64,912)
        if(PBC) write(64,920) Period(1,1),Period(2,2),Period(3,3) 
        do j=1,Nat
!          call chart_int(j,A3)
          A3 = fstr(j)
          A5 = trim(adjustl(Atoms(j)))//trim(adjustl(A3)) 
          write(64,901) A5,Coord(1:3,j),frag(j),Atoms(j)
        enddo
        write(64,911) 
        write(64,911) 
      close(64)
 905  format('!BIOSYM archive 3')
 904  format('PBC=ON ')
 906  format('PBC=OFF')
 915  Format(A60,6x,I5)
 912  format('!DATE     Jan 01 00:00:00 2020')
 901  format(A5,3F15.9,1x,'XXXX 'A6,' xx      ',A2,2x,'0.000')
 902  format(A5,1x,3F15.9,1x,A6,'      xx      ',A2,2x,'0.000')
 903  format(A5,1x,3F15.9,1x,A6,'      xx      ',A2,2x,'0.000')
 920  format('PBC',3F10.4,'   90.0      90.0      90.0')
 911  format('end')
   end subroutine write_car





   subroutine write_car_NZn(name,comment)
      character(*),intent(in)  :: name
      character(*),intent(in)  :: comment
      character*3              :: A3            
      character*5              :: A5            
      character(7)             :: fstr
      integer                  :: j
      if(Nat>999) then
       print *,'  error write_car: N>999'
       return
      endif
      open(64,file=name,status='unknown',form='formatted')
        write(64,905) 
        if(PBC) then
         write(64,904)
        else
         write(64,906)
        endif 
        write(64,915) comment,Nat
        write(64,912)
        if(PBC) write(64,920) Period(1,1),Period(2,2),Period(3,3) 
        do j=1,Nat
          if(Atoms(j)(1:1)=='N') then
!           call chart_int(NZn(j),A3)
           A3 = fstr(NZn(j))
          else
!           call chart_int(j,A3)
           A3 = fstr(j)
          endif
          A5 = trim(adjustl(Atoms(j)))//trim(adjustl(A3)) 
          write(64,901) A5,Coord(1:3,j),frag(j),Atoms(j)
        enddo
        write(64,911) 
        write(64,911) 
      close(64)
 905  format('!BIOSYM archive 3')
 904  format('PBC=ON ')
 906  format('PBC=OFF')
 915  Format(A60,6x,I5)
 912  format('!DATE     Jan 01 00:00:00 2014')
 901  format(A5,3F15.9,1x,'XXXX 'A6,' xx      ',A2,2x,'0.000')
 902  format(A5,1x,3F15.9,1x,A6,'      xx      ',A2,2x,'0.000')
 903  format(A5,1x,3F15.9,1x,A6,'      xx      ',A2,2x,'0.000')
 920  format('PBC',3F10.4,'   90.0      90.0      90.0')
 911  format('end')
   end subroutine write_car_NZn





   subroutine write_arc(name,comment)
      character(*),intent(in)  :: name
      character(*),intent(in)  :: comment
      character*3              :: A3            
      character*5              :: A5            
      character(7)             :: fstr
      integer                  :: kt,j
      if(Nat>999) then
       print *,'  error write_car: N>999'
       return
      endif
      print *,'you need to update Coord(:,:,:)'
      stop
      open(64,file=name,status='unknown',form='formatted')
        write(64,905) 
        if(PBC) then
         write(64,904)
        else
         write(64,906)
        endif 
        do kt=1,Nt
          write(64,915) comment,Nat
          write(64,912)
          if(PBC) write(64,920) Period 
          do j=1,Nat
!            call chart_int(j,A3)
            A3 = fstr(j)
            A5 = trim(adjustl(Atoms(j)))//trim(adjustl(A3)) 
            write(64,901) A5,Coord(1:3,j),frag(j),Atoms(j)    !Coord(1:3,j,kt),frag(j),Atoms(j)
          enddo
          write(64,911) 
          write(64,911) 
        enddo
      close(64)
 905  format('!BIOSYM archive 3')
 904  format('PBC=ON ')
 906  format('PBC=OFF')
 915  Format(A60,6x,I5)
 912  format('!DATE     Jan 01 00:00:00 2014')
 901  format(A5,3F15.9,1x,'XXXX 'A6,' xx      ',A2,2x,'0.000')
 902  format(A5,1x,3F15.9,1x,A6,'      xx      ',A2,2x,'0.000')
 903  format(A5,1x,3F15.9,1x,A6,'      xx      ',A2,2x,'0.000')
 920  format('PBC',3F10.4,'   90.0      90.0      90.0')
 911  format('end')
   end subroutine write_arc




  subroutine remove(C,A,N,Nremove)                 ! remove one atom
   real*8,      intent(inout)  :: C(3,N)     
   character*1, intent(inout)  :: A(N)     
   integer,     intent(inout)  :: N
   integer,     intent(in)     :: Nremove
   integer                     :: j
    if(Nremove > N) then
	  print *,'DMol remove: error Nremove > Natoms'
	else
      do j=Nremove+1,N
        C(1:3,j-1) = C(1:3,j)
        A(j-1)     = A(j)
	  enddo
	  N = N-1
	endif
  end subroutine remove



  subroutine read_DMol_potential(name)
    character(*)       :: name
    real(8)            :: a,b,c
    real(8)            :: dx,dy,dz
    integer            :: i,j,k
    integer            :: NA,N1,N2,N3,N4,N5,N6
    integer            :: jj
    open(unit=2,file=name)
    read(2,*)
    read(2,*)
    read(2,*) a,b,c
    read(2,*) Nxpot,Nypot,Nzpot
    dx = a/dfloat(Nxpot)
    dy = b/dfloat(Nypot)
    dz = c/dfloat(Nzpot)
    Nxpot = Nxpot + 1
    Nypot = Nypot + 1
    Nzpot = Nzpot + 1
    Npot = Nxpot*Nypot*Nzpot
    read(2,*) NA,N1,N2,N3,N4,N5,N6
    allocate(Pot(Npot))
    allocate(Pos(3,Npot))
    jj = 0
    do k=1,Nzpot
     do j=1,Nypot
      do i=1,Nxpot
       jj = jj + 1
       if(jj > Npot) then
        print *,'read_potential: error jj > Npot'
        stop
       endif
       Pos(1,jj) = (N1+i-1)*dx 
       Pos(2,jj) = (N3+j-1)*dy 
       Pos(3,jj) = (N5+k-1)*dz
       read(2,*) Pot(jj)
      enddo
     enddo
    enddo
  end subroutine read_DMol_potential




end Module DMol_files



