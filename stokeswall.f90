module stokeswall
!**************************************************************************
! Regularized Stokeslet Solver
!**************************************************************************
    implicit none
     contains
!==========================================================================
!=======================================================================
! Routine for stokes flow with images, for a wall at at x = wallLoc
    subroutine StokesletWall(xb,vb,xm,fm,nv,nf,delta)
        use parameters, only: mu,pi,dp,wallLoc
        implicit none
        include "omp_lib.h"

        integer,intent(in) :: nv,nf
        real(dp),intent(in) :: xb(3,nv),delta
        real(dp),intent(in) :: xm(3,nf),fm(3,nf)
        real(dp),intent(out) :: vb(3,nv)
        real(dp) :: xmim(3,nf)
        real(dp),dimension(nv) :: R,H1,H2,fac
        real(dp),dimension(nv) :: dx,dy,dz,fdotx
        real(dp),dimension(nv) :: pdotx,qdotx
        real(dp),dimension(nv) :: Hdp1,Hdp2,Hr1,Hdb2
        real(dp),dimension(nf) :: h0
        real(dp),dimension(3) :: p,q,L
        integer :: kk,ii
        !-----------------------------------------
        !zero out velocity vectors
        vb = 0.d0
        !vb(1,:)=0.0d0; vb(2,:)=0.0d0; vb(3,:)=0.0d0;
        dx=0.d0; dy=0.d0; dz=0.d0;
        fdotx=0.d0; H1=0.d0; H2=0.d0;
        !go through all points linearly
        !$OMP PARALLEL PRIVATE(dx,dy,dz,R,fac,fdotx,H1,H2)
        !$OMP DO REDUCTION (+:vb)
        do kk=1,nf
           !blob is: 15delta^4/(8pi(r^2+delta^2)^(7/2))
           dx=(xb(1,:)-xm(1,kk))
           dy=(xb(2,:)-xm(2,kk))
           dz=(xb(3,:)-xm(3,kk))
           R=SQRT(dx**2+dy**2+dz**2+delta**2)
           H1=(R**2+delta**2)/(2.d0*(R**3))
           H2=1.d0/(2.d0*(R**3))
           fdotx=(fm(1,kk)*dx+fm(2,kk)*dy+fm(3,kk)*dz)
           vb(1,:)=vb(1,:)+fm(1,kk)*H1+fdotx*dx*H2
           vb(2,:)=vb(2,:)+fm(2,kk)*H1+fdotx*dy*H2
           vb(3,:)=vb(3,:)+fm(3,kk)*H1+fdotx*dz*H2
        end do !ending nf force points interacting loop
        !$OMP END DO
        !$OMP END PARALLEL

        !Now the images

        h0 = xm(1,:) - wallLoc
        !image points:
        xmim = xm
        xmim(1,:) = 2.d0*wallLoc-xm(1,:)

        !$OMP PARALLEL PRIVATE(dx,dy,dz,R,fdotx,qdotx,pdotx,p,q,L,H1,H2,Hdp1,Hdp2,Hdb2,Hr1)
        !$OMP DO REDUCTION (+:vb)
        do kk=1,nf
           !these are the coefficients for the dipole and other elements
           p(1) = 2.d0*h0(kk)*fm(1,kk)
           p(2) = -2.d0*h0(kk)*fm(2,kk)
           p(3) = -2.d0*h0(kk)*fm(3,kk)
           q = -h0(kk)*p/2.d0
           L(1) = 0.d0
           L(2) = -2.d0*h0(kk)*fm(3,kk)
           L(3) = 2.d0*h0(kk)*fm(2,kk)

           dx = xb(1,:)-xmim(1,kk)
           dy = xb(2,:)-xmim(2,kk)
           dz = xb(3,:)-xmim(3,kk)

           !Needed functions
           R=SQRT(dx**2+dy**2+dz**2+delta**2)
           H1=(R**2+delta**2)/(2.d0*(R**3))
           H2=1.d0/(2.d0*(R**3))
           Hdp2 = -3.d0/(R**5)
           Hdp1 = 2.d0*H2 + (delta**2)*Hdp2
           Hr1 = -(delta**2)*Hdp2/2.d0
           Hdb2 = - (Hr1 + H2)
    
           !Stokeslet image
           fdotx = fm(1,kk)*dx + fm(2,kk)*dy + fm(3,kk)*dz
           vb(1,:)= vb(1,:) - fm(1,kk)*H1 - fdotx*dx*H2
           vb(2,:)= vb(2,:) - fm(2,kk)*H1 - fdotx*dy*H2
           vb(3,:)= vb(3,:) - fm(3,kk)*H1 - fdotx*dz*H2

           !Dipole image
           qdotx = q(1)*dx + q(2)*dy + q(3)*dz
           vb(1,:) = vb(1,:) + q(1)*Hdp1 + qdotx*dx*Hdp2
           vb(2,:) = vb(2,:) + q(2)*Hdp1 + qdotx*dy*Hdp2
           vb(3,:) = vb(3,:) + q(3)*Hdp1 + qdotx*dz*Hdp2
   
           !Doublet image
           pdotx = p(1)*dx + p(2)*dy + p(3)*dz
           vb(1,:)= vb(1,:) + p(1)*dx*H2 + p(1)*dx*H2 + pdotx*dx*dx*Hdp2/2.d0 + pdotx*Hdb2
           vb(2,:)= vb(2,:) + p(1)*dy*H2 + p(2)*dx*H2 + pdotx*dx*dy*Hdp2/2.d0
           vb(3,:)= vb(3,:) + p(1)*dz*H2 + p(3)*dx*H2 + pdotx*dx*dz*Hdp2/2.d0

           !Rotlets image
           vb(1,:) = vb(1,:) + (L(2)*dz-L(3)*dy)*Hr1
           vb(2,:) = vb(2,:) + (L(3)*dx-L(1)*dz)*Hr1
           vb(3,:) = vb(3,:) + (L(1)*dy-L(2)*dx)*Hr1
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        vb = vb/(4.d0*pi*mu)
        return
    end subroutine StokesletWall

end module stokeswall
