module changeInit
!**************************************************************
! To test different initial configurations
!**************************************************************
    implicit none
     contains

    subroutine moveXtoZ(xs)
!==============================================================
! To change flagellum data so that X init data is now in Z.
!==============================================================
        use parameters, only: dp,np,SpNum
        implicit none
        real(dp),intent(inout) :: xs(3,np*SpNum)
        real(dp) :: xsr(3,np*SpNum)
       
        xsr = xs; 
        xs(1,:) = xsr(3,:)
        xs(3,:) = xsr(1,:)

    end subroutine moveXtoZ

    subroutine nonStandPlane(xs)
!==============================================================
! To rotate flagellum data to a non-standard plane
!==============================================================
        use parameters, only: dp,np,SpNum
        implicit none
        real(dp),intent(inout) :: xs(3,np*SpNum)
        real(dp) :: RMinv(3,3,SpNum),xsr(3,np*SpNum)
        integer :: kk
        
        xsr = xs
        ! Rotation Matrix (example)
        RMinv(1,1,:)= 0.866025403784439d0
        RMinv(1,2,:)=  -0.5d0
        RMinv(1,3,:)= 0.d0
        RMinv(2,1,:)= 0.353553390593274d0
        RMinv(2,2,:)= 0.612372435695795d0  
        RMinv(2,3,:)= -0.707106781186547d0
        RMinv(3,1,:)= 0.353553390593274d0
        RMinv(3,2,:)= 0.612372435695794d0
        RMinv(3,3,:)= 0.707106781186548d0
        
        do kk=1,SpNum
           xs(:,1+(kk-1)*np:kk*np) = MATMUL(RMinv(:,:,kk),&
                xsr(:,1+(kk-1)*np:kk*np))
        enddo
   end subroutine nonStandPlane

end module changeInit
