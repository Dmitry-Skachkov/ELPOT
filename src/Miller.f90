




 Module Miller                                     ! calculate Miller indexes
  implicit none
 Contains




      subroutine miller_indexes(Qe,Aa,Ab,Ac,u,v,w)
       real(8)   :: Qe(3)
       real(8)   :: Aa(3),Ab(3),Ac(3)
       real(8)   :: u,v,w
       real(8)   :: M(3,3),Y(3)
       M(1,1:3) = Aa(1:3)
       M(2,1:3) = Ab(1:3)
       M(3,1:3) = Ac(1:3)
       call solve3x3(M,Qe,Y)                       ! solve 3x3 equation M*Y=Qe
       u = Y(1)
       v = Y(2)
       w = Y(3)
      end subroutine miller_indexes

 



     subroutine solve3x3(A,B,X)                    ! solve 3x3 linear equation AX=B https://en.wikipedia.org/wiki/Cramer%27s_rule
      real(8)  :: A(3,3)
      real(8)  :: B(3)
      real(8)  :: X(3)
      real(8)  :: D,D1,D2,D3
      real(8)  :: Aw(3,3)
      call det3x3(A,D)
      call make_dx(A,1,B,Aw)
      call det3x3(Aw,D1)
      call make_dx(A,2,B,Aw)
      call det3x3(Aw,D2)
      call make_dx(A,3,B,Aw)
      call det3x3(Aw,D3)
      X(1) = D1/D
      X(2) = D2/D
      X(3) = D3/D
     end subroutine solve3x3



     subroutine make_dx(A,j,B,Aw)                  ! substitute column j in A
      real(8)    :: A(3,3),B(3)
      integer    :: j
      real(8)    :: Aw(3,3)
      Aw(1:3,1:3) = A(1:3,1:3)
      Aw(j,1:3) = B(1:3)
     end subroutine make_dx



    
     subroutine det3x3(A,D)                        ! determinant of 3x3
      real(8)  :: A(3,3)
      real(8)  :: D
      D = A(1,1)*A(2,2)*A(3,3)  &
        + A(1,2)*A(2,3)*A(3,1)  &
        + A(1,3)*A(2,1)*A(3,2)  &
        - A(1,3)*A(2,2)*A(3,1)  &
        - A(1,1)*A(2,3)*A(3,2)  &
        - A(1,2)*A(2,1)*A(3,3)  
     end subroutine det3x3





 end module Miller










      
