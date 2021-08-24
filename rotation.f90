module rotation
!**************************************************************
! Rotation functions
!**************************************************************
    implicit none
     contains
    subroutine rotateFlag(xs,xbar,xsr,normVect,RMinv)
!==============================================================
! To rotate the flagellar planes to the x-y plane
! Idea is to translate data to be centered at (0,0,0)
! Planes Ax+By+Cz = 0 are known (stored in normVect)
! Rotate data about the axis which is cross product 
! of normVect and (0,0,1), by angle computed by dot
! product between normVect and (0,0,1). Also computes
! inverse rotation matrix and avg. xs values to be 
! used to translate data back to original frame.
!==============================================================
        use parameters, only: dp,np,npd,SpNum
        implicit none
        real(dp),intent(in) :: xs(3,np*SpNum),normVect(3,SpNum)
        real(dp),intent(out) :: xbar(3,SpNum),xsr(3,np*SpNum)
        real(dp),intent(out) :: RMinv(3,3,SpNum)
        real(dp) :: rV(3),cosang,sinang
        real(dp) :: xsc(3,np),RM(3,3)
        integer :: kk,ii,jj
        !-----------------------------------------
        !zero out quantities
        xsr = 0.d0; rV = 0.d0; RM = 0.d0; RMinv = 0.d0
        do kk=1,SpNum
           ! First make sure we need to rotate at all:
           if ((normVect(1,kk).ne.0.d0).or.(normVect(2,kk).ne.0.d0)) then
                ! Find average values
                xbar(:,kk) =SUM(xs(:,1+(kk-1)*np:kk*np),DIM=2)/npd
                ! Now center xs around (0,0,0)
                xsc(1,:) = xs(1,1+(kk-1)*np:kk*np)-xbar(1,kk)
                xsc(2,:) = xs(2,1+(kk-1)*np:kk*np)-xbar(2,kk)
                xsc(3,:) = xs(3,1+(kk-1)*np:kk*np)-xbar(3,kk)
                ! Desired rotated normal vector is (0,0,1)
                ! Cross product gives axis of rotation rV
                rV(1) = normVect(2,kk)
                rV(2) = -normVect(1,kk)
               ! Normalize
                rV = rV/(DSQRT(SUM(rV**2)))
                ! To rotate data:
                cosang = normVect(3,kk) ! from dot product
                sinang = DSQRT(1.d0-cosang**2)
                call rotateMat(rV,cosang,sinang,RM)
                !xsr(:,1+(kk-1)*np:kk*np) = MATMUL(RM,xsc)
                ! Some issue with Matmul for some reason?
                do jj =1,np
                   xsr(:,jj+(kk-1)*np) = MATMUL(RM,xsc(:,jj))
                enddo

                call rotateMat(rV,cosang,-sinang,RMinv(:,:,kk))
           else !No rotation
                RMinv(1,1,kk) =1.d0; 
                RMinv(2,2,kk) = 1.d0; 
                RMinv(3,3,kk) = 1.d0;
                xsr(:,1+(kk-1)*np:kk*np) = xs(:,1+(kk-1)*np:kk*np);
           endif
        end do !ending sperm Loop
        return
    end subroutine rotateFlag

!=======================================================================
    subroutine rotateInit(xs,theta)
!==============================================================
! To rotate the initial flagellar positions by the angle theta
! about (0,1,0)
!==============================================================
        use parameters, only: dp,np,npd,nptot,SpNum
        implicit none
        real(dp),intent(inout) :: xs(3,nptot)
        real(dp),intent(in) :: theta
        real(dp) :: xbar(3,SpNum),xsr(3,nptot)
        real(dp) :: rV(3),cosang,sinang
        real(dp) :: xsc(3,np),RM(3,3)
        integer :: kk,ii,jj
        !-----------------------------------------
        !initialize quantities
        xsr = xs; rV = 0.d0; RM = 0.d0; rV(2) = 1.d0

        do kk=1,SpNum
           ! Find average values
           xbar(1,kk) =SUM(xs(1,1+(kk-1)*np:kk*np))/npd
           xbar(2,kk) =SUM(xs(2,1+(kk-1)*np:kk*np))/npd
           xbar(3,kk) =SUM(xs(3,1+(kk-1)*np:kk*np))/npd
           ! Now center xs around (0,0,0)
           xsc(1,:) = xs(1,1+(kk-1)*np:kk*np)-xbar(1,kk)
           xsc(2,:) = xs(2,1+(kk-1)*np:kk*np)-xbar(2,kk)
           xsc(3,:) = xs(3,1+(kk-1)*np:kk*np)-xbar(3,kk)
           ! To rotate data:
           cosang = DCOS(theta)
           sinang = DSIN(theta)
           call rotateMat(rV,cosang,sinang,RM)
           do jj =1,np
              xs(:,jj+(kk-1)*np) = MATMUL(RM,xsc(:,jj))
           enddo
           !Translate back to xbar
           xs(1,1+(kk-1)*np:kk*np) = xs(1,1+(kk-1)*np:kk*np)+xbar(1,kk)
           xs(2,1+(kk-1)*np:kk*np) = xs(2,1+(kk-1)*np:kk*np)+xbar(2,kk)
           xs(3,1+(kk-1)*np:kk*np) = xs(3,1+(kk-1)*np:kk*np)+xbar(3,kk)
        end do !ending sperm Loop
        return
    end subroutine rotateInit


!=======================================================================
    subroutine rotateMat(rV,cosang,sinang,RM)
    !===================================================================
    ! This routine takes a unit vector rV (of length 3) the cosine of an 
    ! angle and the sine of an angle and returns RM (3x3), a rotation 
    ! matrix for rotating by this angle with an axis of rotation in the 
    ! direction of rV.
        use parameters, only: dp
        implicit none
        real(dp),intent(in) :: rV(3),cosang,sinang
        real(dp),intent(out) :: RM(3,3)
       
        RM = 0.d0
        ! Rotation matrix:
        RM(1,1) = rV(1)**2+(rV(2)**2+rV(3)**2)*cosang
        RM(1,2) = rV(1)*rV(2)*(1.d0-cosang)-rV(3)*sinang
        RM(1,3) = rV(1)*rV(3)*(1.d0-cosang)+rV(2)*sinang
        RM(2,1) = rV(1)*rV(2)*(1.d0-cosang)+rV(3)*sinang
        RM(2,2) = rV(2)**2+(rV(1)**2+rV(3)**2)*cosang
        RM(2,3) = rV(2)*rV(3)*(1.d0-cosang)-rV(1)*sinang
        RM(3,1) = rV(1)*rV(3)*(1.d0-cosang)-rV(2)*sinang
        RM(3,2) = rV(2)*rV(3)*(1.d0-cosang)+rV(1)*sinang
        RM(3,3) = rV(3)**2+(rV(1)**2+rV(2)**2)*cosang
        return
    end subroutine rotateMat


!=======================================================================
    subroutine rotateBack(xsr,xs,M)
    !===================================================================
    ! This routine simply takes rotation matrices M (3x3xSpNum) and an 
    ! array xsr (3xnpxSpNum) and multiplies. Defined this instead of
    ! MATMUL due to indexing issues.
    !===================================================================
        use parameters, only : dp,np,SpNum
        implicit none

        real(dp), intent(in) :: xsr(3,np*SpNum), M(3,3,SpNum)
        real(dp), intent(out) :: xs(3,np*SpNum)
        integer :: ii,jj,kk

        do kk=1,SpNum
           do ii=1,3
              do jj=1,np
                xs(ii,jj+(kk-1)*np)=SUM(M(ii,:,kk)*xsr(:,jj+(kk-1)*np))
              enddo
           enddo
        enddo
    end subroutine rotateBack        
!=======================================================================

end module rotation
