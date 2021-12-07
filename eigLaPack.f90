module eigLaPack
    implicit none
    contains

!**************************************************************
!* To construct matrix for least squares plane fit and find
!* normal vector to plane
!--------------------------------------------------------------
    subroutine normalPlane(xs,normVector,meth)
            use parameters, only: dp,np,npd
            implicit none
            real(dp), intent(in) :: xs(3,np)
            real(dp), intent(out) :: normVector(3)
            integer, intent(out) :: meth
            real(dp) :: A(3,3),xsbar(3),lambda
            real(dp) :: u(3),v(3),det
            integer :: maxL(1),minL(1),xminL(1)

                ! Find center of mass:
                xsbar = sum(xs(:,:),DIM=2)/npd
                ! Entries of A take the following values
                ! See Least Squares Fitting of Data
                ! by David Eberly (Geomtric Tools
                A(1,1) = sum((xs(1,:)-xsbar(1))**2);
                A(2,2) = sum((xs(2,:)-xsbar(2))**2);
                A(3,3) = sum((xs(3,:)-xsbar(3))**2);
                A(1,2) = sum((xs(1,:)-xsbar(1))*(xs(2,:)-xsbar(2)));
                A(1,3) = sum((xs(1,:)-xsbar(1))*(xs(3,:)-xsbar(3)));
                A(2,3) = sum((xs(2,:)-xsbar(2))*(xs(3,:)-xsbar(3)));
                !symmetric matrix:
                A(2,1)=A(1,2); A(3,1)=A(1,3); A(3,2)=A(2,3)
                !Check determinant:
                det = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))-&
                      A(1,2)*(A(2,1)*A(3,3)-A(3,1)*A(2,3))+&
                      A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
                if (det==0.d0) then !co-planar
                   minL(1)=2; maxL(1)=3; xminL=1
                   u(1) = xs(1,maxL(1))-xs(1,xminL(1))
                   u(2) = xs(2,maxL(1))-xs(2,xminL(1))
                   u(3) = xs(3,maxL(1))-xs(3,xminL(1))
                   v(1) = xs(1,minL(1))-xs(1,xminL(1))
                   v(2) = xs(2,minL(1))-xs(2,xminL(1))
                   v(3) = xs(3,minL(1))-xs(3,xminL(1))
                   !cross product
                   normVector(1) = u(2)*v(3)-u(3)*v(2)
                   normVector(2) = u(3)*v(1)-u(1)*v(3)
                   normVector(3) = u(1)*v(2)-u(2)*v(1)
                   meth = 0
                else
                   call findEigVec(A,normVector,lambda)
                   if (lambda==0.d0) then ! det close to 0
                      call findUV(u,v,xs)
                      !cross product
                      normVector(1) = u(2)*v(3)-u(3)*v(2)
                      normVector(2) = u(3)*v(1)-u(1)*v(3)
                      normVector(3) = u(1)*v(2)-u(2)*v(1)
                      meth = 2
                   else
                      meth = 1
                   endif
                endif
                !Normalize:
                normVector = normVector/(SQRT(SUM(normVector**2)))
                return
    end subroutine normalPlane

!**************************************************************
!* Find indexes to calculate cross product more accurately
!**************************************************************
    subroutine findUV(u,v,xs)
            use parameters,only : dp,np
            implicit none
        
            real(dp), intent(in) :: xs(3,np)
            real(dp), intent(out) :: u(3),v(3)
            real(dp) :: m1(3),m2(3),m3(3)
            integer :: i1,i2,i3,i4
        
            !Define indexes: (want asymmetry in indexes) 
            i1 = np/5 
            i2 = np/3
            i3 = 3*np/5
            i4 = 4*np/5
            
            m1 = SUM(xs(:,1:i1),DIM=2)/DBLE(i1)
            m2 = SUM(xs(:,i2:i3),DIM=2)/DBLE(i3-i2+1)
            m3 = SUM(xs(:,i4:np),DIM=2)/DBLE(np-i4+1)
            

            u = m3-m1; 
            v = m2-m1; 
    end subroutine findUV

!**************************************************************
! To find smallest eigenvector for a 3x3 matrix A
!   Uses DSYEVD 
!**************************************************************
    subroutine findEigVec(A,v,lambda)
        use parameters, only: dp
        implicit none

        real(dp),intent(in) :: A(3,3)
        real(dp),intent(out) :: v(3),lambda
        integer :: N,LDA,LWMAX,INFO,LWORK,LIWORK
        integer :: IWORK(1000),posNum,ii
        real(dp) :: W(3),WORK(1000)

        N=3; LDA=3; LWMAX=1000;

        ! External Subroutines 
        !EXTERNAL         DSYEVD
        ! Intrinsic Functions 
        !INTRINSIC        INT, MIN

        ! Query the optimal workspace.
        LWORK = -1
        LIWORK = -1
        CALL DSYEVD('Vectors','Upper', N, A, LDA, W, WORK, LWORK,&
                  IWORK, LIWORK, INFO )
        LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
        LIWORK = MIN( LWMAX, IWORK( 1 ) )

        ! SOLVE EIGENVALUE PROBLEM:
        CALL DSYEVD( 'Vectors', 'Upper', N, A, LDA, W, WORK, LWORK,&
                 IWORK, LIWORK, INFO )
        !     Check for convergence.
        IF( INFO.GT.0 ) THEN
          !The algorithm failed to compute eigenvalues, return zeros
          lambda = 0.d0; v = 0.d0
          return
        ELSE
        ! Find min eigenvalue:
        lambda = W(1)
        posNum = 1
        do ii=1,2
           if (W(ii).lt.lambda) then
                lambda=W(ii)
                posNum=ii
           endif
        enddo
        ! Find associated evector:
        v = A(:,posNum)
        ENDIF
        return
      end subroutine findEigVec

end module eigLaPack
