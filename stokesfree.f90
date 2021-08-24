module stokesfree
!**************************************************************************
! Regularized Stokeslet Solvers
!**************************************************************************
    implicit none
     contains
!==========================================================================
!=======================================================================
!==========================================================================
    subroutine StokesletFree(xb,vb,xm,fm,nv,nf,delta)
        use parameters, only: mu,pi,dp
        implicit none
        include "omp_lib.h"
        integer :: nv,nf
        real(dp),intent(in) :: xb(3,nv),xm(3,nf),fm(3,nf),delta
        real(dp),intent(out) :: vb(3,nv)
        real(dp),dimension(nv) :: r2,H1,H2,fac
        real(dp),dimension(nv) :: dx,dy,dz,fdotx
        integer :: kk,ii
        !-----------------------------------------
        !zero out vectors
        vb=0.0d0;
        dx=0.d0; dy=0.d0; dz=0.d0;
        fdotx=0.d0; H1=0.d0; H2=0.d0;
        !go through all points linearly
        !$OMP PARALLEL PRIVATE(dx,dy,dz,r2,fac,fdotx,H1,H2)
        !$OMP DO REDUCTION (+:vb)
        do kk=1,nf !go through all forces
           !blob is: 15delta^4/(8pi(r^2+delta^2)^(7/2))
           dx=(xb(1,:)-xm(1,kk))
           dy=(xb(2,:)-xm(2,kk))
           dz=(xb(3,:)-xm(3,kk))
           r2=(dx**2+dy**2+dz**2)
           fac=(r2+delta**2)
           H1=(r2+2.d0*delta**2)/(8.d0*pi*mu*((fac)**(1.5d0)))
           H2=1.d0/(8.d0*pi*mu*((fac)**(1.5d0)))
           fdotx=(fm(1,kk)*dx+fm(2,kk)*dy+fm(3,kk)*dz)
           vb(1,:)=vb(1,:)+fm(1,kk)*H1+fdotx*dx*H2
           vb(2,:)=vb(2,:)+fm(2,kk)*H1+fdotx*dy*H2
           vb(3,:)=vb(3,:)+fm(3,kk)*H1+fdotx*dz*H2
        end do !ending nf force points interacting loop
        !$OMP END DO
        !$OMP END PARALLEL
        return
    end subroutine StokesletFree

end module stokesfree
