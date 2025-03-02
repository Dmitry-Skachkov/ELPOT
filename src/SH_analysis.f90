


  module SH_analysis                                              ! spherical harmonics
   use parameters
   use shtools
   use ftypes
   implicit none
   real(8)             :: cilm(2,lmax+1,lmax+1)                   ! spherical harmonics expansion coefficients
   real(8)             :: cilm2(2,lmax+1,lmax+1,10)
   integer             :: lmax2
   real(8)             :: Pots2(ns,ns)
  contains



  subroutine calc_stevens
   integer          :: ir,l,m,l0
   real(8)          :: rs
   do ir=1,10
    print *,'ir=',ir
    Pots2(:,:) = Pots(:,:,ir)                                     ! potential at theta_i, i=1,ns and phi_j, j=1,ns
    call SHExpandDH(Pots2, ns, cilm, lmax2, sampling = 1)         ! calculate spherical harmonics
    cilm2(:,:,:,ir) = cilm(:,:,:)
    if(lmax/=lmax2) print *,'WARNING:   lmax2=',lmax2
   enddo
   print *
   print *,'spherical harmonics'
   do l=1,lmax+1
    l0 = l-1                                                      ! index l of C_lm
    do m=-l0,l0                                                   ! index m of C_lm
     print *
     print *,'l=',l0,' m=',m
     do ir=1,10
      if(m>=0) print 1,cilm2(1,l,m+1,ir)
      if(m<0)  print 1,cilm2(2,l,-m+1,ir)
     enddo 
    enddo  
   enddo
   print *
   print *,'Stevens coefficients'
   do l=1,lmax+1
    l0 = l-1                                                     ! index l of C_lm
    do m=-l0,l0                                                  ! index m of C_lm
     print *
     print *,'l=',l0,' m=',m
     do ir=1,10
      rs = ir*0.1d0
      if(m>=0) print 1,cilm2(1,l,m+1,ir)/(rs**l0)
      if(m<0)  print 1,cilm2(2,l,-m+1,ir)/(rs**l0)
     enddo 
    enddo  
   enddo
1  format(E16.5)
  end subroutine calc_stevens



  end module SH_analysis


