module axonemeForces
!**************************************************************************
!**************************************************************************
    implicit none
    private :: normThree,Curve
    contains
!==========================================================================
!==========================================================================
    subroutine normThree(xk,xkk,norm3)
        use parameters, only: dp
        implicit none
        real(dp),intent(in) :: xk(3),xkk(3)
        real(dp),intent(out) :: norm3

        norm3=DSQRT(SUM((xk-xkk)**2));
        return
    end subroutine normThree
!==========================================================================
!==========================================================================
    subroutine CrossProd(xk,xkp,xkm,cross)
        use parameters, only: np,dp
        implicit none
        real(dp),intent(in) :: xk(2),xkp(2),xkm(2)
        real(dp),intent(out) :: cross

        cross=(xkp(2)-xk(2))*(xk(1)-xkm(1))-(xkp(1)-xk(1))*(xk(2)-xkm(2));
        return
    end subroutine CrossProd
!==========================================================================
!==========================================================================
    subroutine ArcLength(x,arcL,kk)
        use parameters, only: SpNum,np,dp
        implicit none
        real(dp),intent(in) :: x(2,np,SpNum)
        real(dp),intent(out) :: arcL
        integer :: i,kk

        arcL=0.d0;
        do i=1,(np-1)
                arcL=arcL+DSQRT(SUM((x(:,i+1,kk)-x(:,i,kk))**2.d0))
        end do
        return
    end subroutine ArcLength   
!==========================================================================
!==========================================================================
    subroutine Curve(ntime,ip,Ckt)
        use parameters, only: cr,np,b,Kp,omega,dp,motType,ds,dt
        implicit none
        real(dp) :: s,nt
        real(dp),intent(out) :: Ckt
        integer,intent(in) :: ip,ntime

        s=ds*DBLE(ip-1);nt=DBLE(ntime)
        Ckt=0.d0;
        Ckt=-Kp*Kp*dsin(Kp*s-omega*nt*dt)
        if (motType==0) then
        !for symmetric, CONSTANT amp, use next line of code
                Ckt=Ckt*b
        elseif (motType==1) then
        !for slight asymm amp
                if(Ckt>0)then
                   Ckt=Ckt*(b+0.05d0); !amp= 0.15
                else
                   Ckt=Ckt*b;  !amp=0.1
                end if
        elseif (motType==2) then
        ! Full asym
                if(Ckt>0)then
                   Ckt=Ckt*(b+0.15d0); !amp= 0.25
                else
                   Ckt=Ckt*b;  !amp=0.1
                end if
        endif
      return
    end subroutine Curve   
!=======================================================================
!=======================================================================
    subroutine spermF(x,xf,ntime)
        use parameters, only: SpNum,np,nptot,ds,S1,S2,S5,S6,S7,&
                              dt,how_often_wr_x,dp,penalty
        implicit none
        real(dp),intent(inout) :: x(3,nptot),xf(3,nptot)
        real(dp) :: crossXY,crossZX,crossYZ
        real(dp) :: Ckt,SC,SCz,SN,SNz,SQ
        real(dp) :: normL,normR,arcL,GM
        real(dp) :: xk(3),xkm(3),xkp(3),Deriv(2,np)
        real(dp) :: zx(2),zxp(2),zxm(2)
        real(dp) :: x2(2,np*SpNum)
        real(dp),dimension(np*SpNum) :: sk,skyz,skzx
        real(dp) :: dsk(np-1),dsk3(np-2),dsk3zx(np-2),dsk3yz(np-2),dsk5(np-2)
        real(dp) :: minDs,maxDs
        integer :: ip,ntime
        integer :: R1,R2,L1,L2,M,kk

        xf=0.d0
        !SC=2.d0*S2/(ds**5)
        SC=2.d0*S2 !unscaled here
        SCz=2.d0*S5/(ds**5)
        SN=2.d0*S1/ds 
        SNz=2.d0*S6/ds
        SQ=2.d0*S7/ds

        do kk=1,SpNum
           !finding delta s's: necessary for 3D curvature
           call getdsk(x(:,(kk-1)*np+1:kk*np),dsk)
           dsk3 = dsk(1:np-2)*dsk(2:np-1)*(dsk(1:np-2)+dsk(2:np-1))/2.d0
           dsk5 = dsk3*dsk(1:np-2)*dsk(2:np-1)
           do ip=1,np
              xk=x(:,ip+(kk-1)*np)
              if(ip==1)then
                !spring force
                xkp=x(:,ip+1+(kk-1)*np);
                call normThree(xk,xkp,normR);
                xf(:,ip+(kk-1)*np)=xf(:,ip+(kk-1)*np)+SN*(normR-ds)*((xkp-xk)/normR);
                ! out of plane penalty components:
                if (penalty==1) then
                   xf(3,ip+(kk-1)*np) = xf(3,ip+(kk-1)*np)-SNz*xk(3);
                elseif (penalty==3) then !penalty 2 taken care of below
                   xf(3,ip+(kk-1)*np)=xf(3,ip+(kk-1)*np)+SQ*(xkp(3)-xk(3));
                endif
              elseif(ip==np)then
                !spring force only
                xkm=x(:,ip-1+(kk-1)*np);
                call normThree(xk,xkm,normL)
                xf(:,ip+(kk-1)*np)=xf(:,ip+(kk-1)*np)+SN*(normL-ds)*((xkm-xk)/normL);
                ! out of plane penalty components:
                if (penalty==1) then
                   xf(3,ip+(kk-1)*np) = xf(3,ip+(kk-1)*np)-SNz*xk(3);
                elseif (penalty==3) then !penalty 2 taken care of below
                   xf(3,ip+(kk-1)*np)=xf(3,ip+(kk-1)*np)+SQ*(-xk(3)+xkm(3));
                endif
              else
                !spring from left, spring from right
                !curvature at interior points
                xkm=x(:,ip-1+(kk-1)*np);
                xkp=x(:,ip+1+(kk-1)*np);
                call normThree(xk,xkm,normL)
                call normThree(xk,xkp,normR);
                xf(:,ip+(kk-1)*np)=xf(:,ip+(kk-1)*np)+SN*(normL-ds)*((xkm-xk)/normL);
                xf(:,ip+(kk-1)*np)=xf(:,ip+(kk-1)*np)+SN*(normR-ds)*((xkp-xk)/normR);
                !Curve function calculates preferred curvature:
                call Curve(ntime,ip,Ckt)
                !Cross product comes into bending forces:
                call CrossProd(xk(1:2),xkp(1:2),xkm(1:2),crossXY)
                if (penalty==2) then
                   zx(1) = xk(1); zx(2) = xk(3);
                   zxp(1) = xkp(1); zxp(2) = xkp(3);
                   zxm(1) = xkm(1); zxm(2) = xkm(3);
                   call CrossProd(zx,zxp,zxm,crossZX)
                   call CrossProd(xk(2:3),xkp(2:3),xkm(2:3),crossYZ)
                   crossYZ=-crossYZ !to match notation in notes
                endif
                
                if (dsk3(ip-1).ne.0.d0) then !should rarely happen
                   xf(1,ip-1+(kk-1)*np)=xf(1,ip-1+(kk-1)*np)-&
                        SC*(crossXY-dsk3(ip-1)*Ckt)*(-(xkp(2)-xk(2)))/dsk5(ip-1)
                   xf(2,ip-1+(kk-1)*np)=xf(2,ip-1+(kk-1)*np)-&
                        SC*(crossXY-dsk3(ip-1)*Ckt)*((xkp(1)-xk(1)))/dsk5(ip-1)
                   xf(1,ip+(kk-1)*np)=xf(1,ip+(kk-1)*np)-&
                        SC*(crossXY-dsk3(ip-1)*Ckt)*((xkp(2)-xk(2))+(xk(2)-xkm(2)))/dsk5(ip-1)
                   xf(2,ip+(kk-1)*np)=xf(2,ip+(kk-1)*np)-&
                        SC*(crossXY-dsk3(ip-1)*Ckt)*(-(xkp(1)-xk(1))-(xk(1)-xkm(1)))/dsk5(ip-1)
                   xf(1,ip+1+(kk-1)*np)=xf(1,ip+1+(kk-1)*np)-&
                        SC*(crossXY-dsk3(ip-1)*Ckt)*(-(xk(2)-xkm(2)))/dsk5(ip-1)
                   xf(2,ip+1+(kk-1)*np)=xf(2,ip+1+(kk-1)*np)-&
                        SC*(crossXY-dsk3(ip-1)*Ckt)*((xk(1)-xkm(1)))/dsk5(ip-1)
                endif
                ! Additional forces for 3D:
                ! Simple penalty force in Z:
                if (penalty==1) then
                        xf(3,ip+(kk-1)*np) = xf(3,ip+(kk-1)*np)-SNz*xk(3);
                ! Curvature penalty force in Z:
                elseif (penalty==2) then
                    ! Now add to forces:
                    ! First ZX terms:
                    if (dsk3zx(ip-1).ne.0.d0) then
                       xf(1,ip-1+(kk-1)*np)=xf(1,ip-1+(kk-1)*np)+SCz*crossZX*(xkp(3)-xk(3))
                       xf(3,ip-1+(kk-1)*np)=xf(3,ip-1+(kk-1)*np)-SCz*crossZX*(xkp(1)-xk(1))
                       xf(1,ip+(kk-1)*np)=xf(1,ip+(kk-1)*np)-SCz*crossZX*(xkp(3)-xkm(3))
                       xf(3,ip+(kk-1)*np)=xf(3,ip+(kk-1)*np)+SCz*crossZX*(xkp(1)-xkm(1))
                       xf(1,ip+1+(kk-1)*np)=xf(1,ip+1+(kk-1)*np)+SCz*crossZX*(xk(3)-xkm(3))
                       xf(3,ip+1+(kk-1)*np)=xf(3,ip+1+(kk-1)*np)-SCz*crossZX*(xk(1)-xkm(1))
                    endif
                    !Now YZ terms:
                    if (dsk3yz(ip-1).ne.0.d0) then
                       xf(2,ip-1+(kk-1)*np)=xf(2,ip-1+(kk-1)*np)-SCz*crossYZ*(xkp(3)-xk(3))
                       xf(3,ip-1+(kk-1)*np)=xf(3,ip-1+(kk-1)*np)+SCz*crossYZ*(xkp(2)-xk(2))
                       xf(2,ip+(kk-1)*np)=xf(2,ip+(kk-1)*np)+SCz*crossYZ*(xkp(3)-xkm(3))
                       xf(3,ip+(kk-1)*np)=xf(3,ip+(kk-1)*np)-SCz*crossYZ*(xkp(2)-xkm(2))
                       xf(2,ip+1+(kk-1)*np)=xf(2,ip+1+(kk-1)*np)-SCz*crossYZ*(xk(3)-xkm(3))
                       xf(3,ip+1+(kk-1)*np)=xf(3,ip+1+(kk-1)*np)+SCz*crossYZ*(xk(2)-xkm(2))
                    endif
                ! Q penalty
                elseif (penalty==3) then
                   xf(3,ip+(kk-1)*np)=xf(3,ip+(kk-1)*np)+SQ*(xkp(3)-2*xk(3)+xkm(3));
                endif !penalties in z
              end if !ip cases
           end do !ip
        end do !kk
        return
    end subroutine spermF
!=======================================================================
!=======================================================================
    subroutine repulsion(xAll,fRep)
    ! Used to keep sperm a distance of ds away from each other
        use parameters, only: SpNum,dp,np,S3,repDist
        implicit none
        real(dp),intent(in) :: xAll(3,np*SpNum)
        real(dp),intent(out) :: fRep(3,np*SpNum)
        integer :: k,kk,i,ii
        real(dp) :: repForce(3),Pt2(3),Pt1(3),ptDist

        fRep = 0.d0
        do k=1,SpNum !outer sperm loop
           do i=1,np !outer point loop
              Pt1 = xAll(:,(k-1)*np+i)
              if (k<SpNum) then !interact with all other sperm
              do kk = k+1,SpNum !inner sperm loop
                 do ii = 1,np !inner point loop
                    Pt2 = xAll(:,(kk-1)*np+ii)
                    ptDist = DSQRT(SUM((Pt2-Pt1)**2))
                    
                    ! if distance within repDist then add repulsive force:
                    if (ptDist<repDist) then
                        repForce = S3*(Pt2-Pt1)*(1.d0/ptDist-1.d0/repDist);
                        ! Add this (with correct sign to xf):
                        fRep(:,(k-1)*np+i)= fRep(:,(k-1)*np+i)-repForce
                        fRep(:,(kk-1)*np+ii)= fRep(:,(kk-1)*np+ii)+repForce
                    endif
                 enddo !inner pt
              enddo    !inner sperm
              endif
           enddo       !outer pt
        enddo          !outer sperm
    end subroutine repulsion
!=======================================================================
!=======================================================================
    subroutine repulsionFlagMesh(xs,xm,fReps,fRepm,M)
    ! Used to keep sperm a distance of ds away from each other
        use parameters, only: SpNum,dp,np,S3,repDist
        implicit none
        integer :: M
        real(dp) :: xs(3,np*SpNum),xm(3,M)
        real(dp) :: fReps(3,np*SpNum),fRepm(3,M)
        integer :: k,i
        real(dp) :: repForce(3),Pt2(3),Pt1(3),ptDist

        fReps = 0.d0
        fRepm = 0.d0
        do k=1,np*SpNum !sperm loop
           do i=1,M !mesh loop
              Pt1 = xs(:,k)
              Pt2 = xm(:,i)
              ptDist = DSQRT(SUM((Pt2-Pt1)**2))
                    
              ! if distance within repDist then add repulsive force:
              if (ptDist<repDist) then
                 repForce = S3*(Pt2-Pt1)*(1.d0/ptDist-1.d0/repDist);
                 ! Add this (with correct sign):
                 fReps(:,k)= fReps(:,k)-repForce
                 fRepm(:,i)= fRepm(:,i)+repForce
               endif
           enddo       !mesh loop
        enddo          !sperm loop
    end subroutine repulsionFlagMesh
!=======================================================================
!=======================================================================
    subroutine attachFlag(xAll,fAtt)
    ! Used to keep sperm a distance of ds away from each other
        use parameters, only: SpNum,dp,np,S4,attDist,attThresh
        implicit none
        real(dp),intent(in) :: xAll(3,np*SpNum)
        real(dp),intent(out) :: fAtt(3,np*SpNum)
        integer :: k,kk,ii
        real(dp) :: attForce(3),Pt2(3),Pt1(3),ptDist

        fAtt = 0.d0
        do k=1,SpNum !outer sperm loop
           !head point
           Pt1 = xAll(:,(k-1)*np+1)
           do kk = 1,SpNum !inner sperm loop
              if (k/=kk) then !interact with all other sperm
                 do ii = 1,np !inner sperm  point loop
                    Pt2 = xAll(:,(kk-1)*np+ii)
                    ptDist = DSQRT(SUM((Pt2-Pt1)**2))
                    
                    ! if distance within attDist then add binding force:
                    if (ptDist<attThresh) then
                        attForce = S4*(Pt2-Pt1)*(ptDist-attDist)/ptDist;
                        ! Add this (with correct sign to xf):
                        fAtt(:,(k-1)*np+1)= fAtt(:,(k-1)*np+1)+attForce !head pt
                        fAtt(:,(kk-1)*np+ii)= fAtt(:,(kk-1)*np+ii)-attForce
                    endif
                 enddo !inner pt
              endif
           enddo    !inner sperm
        enddo       !outer sperm
    end subroutine attachFlag
!=======================================================================
!=======================================================================
    subroutine getdsk(x,dsk)
        use parameters, only: dp,np,ds
        implicit none
        real(dp),intent(in) :: x(3,np)
        real(dp),intent(out) :: dsk(np-1)
        real(dp) :: dsi(3)
        integer :: ii
        
        do ii=1,np-1
           dsi = x(:,ii+1)-x(:,ii) !vector direction we want
           dsi = ds*dsi/DSQRT(SUM(dsi**2)) !make length ds
           dsk(ii) = DSQRT(SUM(dsi(1:2)**2)) !length of projected dsi in xy
        enddo
    end subroutine getdsk


end module axonemeForces
