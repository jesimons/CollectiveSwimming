module initialization

    implicit none
    contains
!=======================================================================
!=======================================================================
    subroutine sperminit(xb,dim1,dim2)
     use parameters, only:np,SpNum,SpDist,omega,cr,b,dp,ds,SpStart,Kp,pi
     implicit none
     real(dp), intent(out) :: xb(3,np*SpNum)
     real(dp) :: s
     integer :: ip,k,dim1,dim2

        loop_ip:  do ip=1,np
                xb(1,ip)=SpStart + dble(ip-1)*ds;
                !slightly non-linear start:
                xb(dim1,ip)=b*DSIN(Kp*(xb(1,ip)-SpStart));
                if(SpNum>1)then
                        do k=2,SpNum
                            xb(1,(k-1)*np+ip)=xb(1,ip)
                            xb(dim1,(k-1)*np+ip)=xb(dim1,ip)+SpDist*DCOS(2.d0*pi*DBLE(k-1)/DBLE(SpNum))
                            xb(dim2,(k-1)*np+ip)=xb(dim2,ip)+SpDist*DSIN(2.d0*pi*DBLE(k-1)/DBLE(SpNum))
                            !xb(dim1,(k-1)*np+ip)=xb(dim1,ip)
                            !xb(dim2,(k-1)*np+ip)=xb(dim2,ip)+SpDist*dble(k-1)                
                        end do
                end if
                xb(dim1,ip)=xb(dim1,ip)+SpDist*DCOS(2.d0*pi*DBLE(k-1)/DBLE(SpNum))
                xb(dim2,ip)=xb(dim2,ip)+SpDist*DSIN(2.d0*pi*DBLE(k-1)/DBLE(SpNum))
        end do loop_ip
        return
   end subroutine sperminit
!=======================================================================
!=======================================================================
    subroutine spermFromData1Sp(xb,nt,file1,dim1)
     use parameters, only:np,SpNum,dp,nptot,SpDist,SpStart
     implicit none
     character (LEN=*),intent(in) :: file1
     integer, intent(in) :: nt,dim1
     real(dp), intent(out) :: xb(3,nptot)
     integer :: ii,kk

        open(unit=4001,file=file1, status='old',action='read')
        ! read first nt lines and do nothing
        do ii=1,np*(nt-1)
           read(4001,*)
        enddo
        ! Retrieve time point desired:
        do ii=1,np
           read(4001,*) xb(:,ii)    !first sperm
        enddo
        !shift data so head is at 0,0,0:
        xb(1,1:np) = xb(1,1:np)-xb(1,1)+SpStart
        xb(2,1:np) = xb(2,1:np)-xb(2,1)
        xb(3,1:np) = xb(3,1:np)-xb(3,1)
        
        do kk=2,SpNum
           xb(:,(kk-1)*np+1:kk*np) = xb(:,1:np)
           xb(dim1,(kk-1)*np+1:kk*np) = xb(dim1,(kk-1)*np+1:kk*np)+SpDist*DBLE(kk-1)
        enddo

        
        return
   end subroutine spermFromData1Sp
!=======================================================================
!=======================================================================
    subroutine spermFromData2Sp(xb,nt,file1,file2)
     use parameters, only:np,SpNum,dp
     implicit none
     character (LEN=*),intent(in) :: file1,file2
     integer, intent(in) :: nt
     real(dp), intent(out) :: xb(3,np*SpNum)
     integer :: ii

        open(unit=4001,file=file1, status='old',action='read')
        open(unit=4002,file=file2, status='old',action='read')
        ! read first nt lines and do nothing
        do ii=1,np*(nt-1)
           read(4001,*)
           read(4002,*)
        enddo
        ! Retrieve time point desire:
        do ii=1,np
           read(4001,*) xb(:,ii)    !first sperm
           read(4002,*) xb(:,ii+np) !second sperm
        enddo
        return
   end subroutine spermFromData2Sp

!=======================================================================
!=======================================================================
    subroutine spermrestart(xb,normVect,numsteps)
     use parameters, only:pi,np,SpNum,SpDist,omega,cr,b,dp,ds
     implicit none
     real(dp), intent(out) :: xb(3,np*SpNum),normVect(3)
     real(dp) :: s
     integer :: kk,ks,vs,numsteps,numwrite,ii
     CHARACTER(LEN=20) :: str,strv


        open(unit=4002,file='progress.txt', status='old',action='read')
        read(4002,*) numsteps
        read(4002,*) numwrite
        close(4002)

        do kk = 1,SpNum
           ks = 40+kk-1
           vs = 90+kk-1
           !str = 'fort.40'
           write(str,"(a,i2)") "fort.",ks
           write(strv,"(a,i2)") "fort.",vs
           open(unit=4001,file=str, status='old',action='read')
           ! read first nt lines and do nothing
           do ii=1,np*(numwrite-1)
              read(4001,*)
           enddo
           ! Retrieve time point desired:
           do ii=1,np
              read(4001,*) xb(:,ii)    ! last sperm loc
              
           enddo
           close(4001)
           open(unit=4003,file=strv, status='old',action='read')
           ! read first nt lines and do nothing
           do ii=1,(numwrite-1)
              read(4003,*)
           enddo
           ! Retrieve normVect:
           read(4003,*) normVect    ! last sperm loc
           close(4003)
        enddo

     return
   end subroutine spermrestart
!=======================================================================
!=======================================================================

!=======================================================================
! ***** OTHER SAMPLE ROUTINES FOR CHANGING INITIAL CONFIGURATIONS *****
!=======================================================================

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


    subroutine switchDim(xs,dim1,dim2)
!================================================================
! To change flagellum data so that dim1 init data is now in dim2.
!================================================================
        use parameters, only: dp,np,SpNum
        implicit none
        real(dp),intent(inout) :: xs(3,np*SpNum)
        integer,intent(in) :: dim1,dim2
        real(dp) :: xsr(3,np*SpNum)
       
        xsr = xs; 
        xs(dim1,:) = xsr(dim2,:)
        xs(dim2,:) = xsr(dim1,:)

    end subroutine switchDim

    subroutine negate(xs,coord)
!==============================================================
! To negate flagellum data in one coordinate
!==============================================================
        use parameters, only: dp,nptot
        implicit none
        real(dp),intent(inout) :: xs(3,nptot)
        integer, intent(in) :: coord
        integer :: ii
       
        xs(coord,:) = -xs(coord,:)

    end subroutine negate

    subroutine nonStandPlane(xs)
!==============================================================
! To rotate flagellum data to a non-standard plane
!==============================================================
        use parameters, only: dp,np,nptot,SpNum
        implicit none
        real(dp),intent(inout) :: xs(3,nptot)
        real(dp) :: RMinv(3,3,SpNum),xsr(3,nptot)
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

    subroutine rotateData(xs,dim1,ang)
!==============================================================
! To rotate flagellum data to a non-standard plane
!==============================================================
        use parameters, only: dp,np,npd,SpNum
        implicit none
        real(dp),intent(inout) :: xs(3,np)
        integer, intent(in) :: dim1
        real(dp),intent(in) :: ang
        real(dp) :: RM(3,3),xsr(3,np),xbar(3)
        integer :: kk
        
        RM = 0.d0
        do kk = 1,3
           xbar(kk) = SUM(xs(kk,:))/npd
           xsr(kk,:) = xs(kk,:)-xbar(kk)
        enddo
        ! Rotation Matrix
        if (dim1==1) then ! in x
           RM(1,1)= 1
           RM(2,2)= DCOS(ang)
           RM(3,3)= DCOS(ang)
           RM(2,3)= -DSIN(ang)
           RM(3,2)= DSIN(ang)
        elseif (dim1==2) then ! in y
           RM(1,1)= DCOS(ang)
           RM(2,2)= 1
           RM(3,3)= DCOS(ang)
           RM(1,3)= DSIN(ang)
           RM(3,1)= -DSIN(ang)
        else  ! in z
           RM(1,1)= DCOS(ang)
           RM(2,2)= DCOS(ang)
           RM(3,3)= 1
           RM(1,2)= -DSIN(ang)
           RM(2,1)= DSIN(ang)
        endif
        xs = MATMUL(RM,xsr)
        do kk=1,3
           xs(kk,:) = xs(kk,:)+xbar(kk)
        enddo
   end subroutine rotateData


end module initialization
