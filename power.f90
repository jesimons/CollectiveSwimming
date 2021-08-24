module power
!**************************************************************************
!**************************************************************************
    implicit none
    contains
!==========================================================================
!==========================================================================
    subroutine findPower(us,vs,fs,ps)
        use parameters, only: dp,SpNum,np,nptot,ds,dt,omega,pi
        implicit none
        real(dp),intent(in),dimension(3,nptot) :: us,vs,fs
        real(dp),intent(inout) :: ps(SpNum)
        real(dp) :: psc(SpNum)
        real(dp),dimension(3,nptot) :: veltot
        integer :: kk,ii,ind,ind2

        veltot=us+vs
        psc=0.d0 !to find current power

        ! Riemann sum for computing power along entire flagellum:
        do kk=1,SpNum
           do ii = 2,np-1 !sum powers over length of flagellum, except endpts
              ind = (kk-1)*np+ii
              psc(kk) = psc(kk)+ds*DOT_PRODUCT(veltot(:,ind),fs(:,ind))
           enddo
           !End pts (trapezoidal rule):
           ind = (kk-1)*np+1; ind2 = kk*np;
           psc(kk) = psc(kk)+0.5d0*ds*(DOT_PRODUCT(veltot(:,ind),fs(:,ind))+&
                        DOT_PRODUCT(veltot(:,ind2),fs(:,ind2)))
        enddo
        ps = ps + (omega/(2.d0*pi))*psc*dt !Riemann integration over time
        return

    end subroutine findPower
!=======================================================================
!=======================================================================
!==========================================================================
!==========================================================================
    subroutine addPower(us,vs,fs,ps,ntime)
        use parameters, only: dp,SpNum,np,nptot,ds,dt,omega,pi,nbeat,&
                                how_often_wr_x
        use outputs;          !write outputs to files
        implicit none
        real(dp),intent(in),dimension(3,nptot) :: us,vs,fs
        real(dp),intent(inout) :: ps(SpNum)
        real(dp) :: psc(nptot), psum
        real(dp),dimension(3,nptot) :: veltot
        integer :: kk,ii,ind,ind2,ntime

        veltot=us+vs
        
        psc = sum(veltot*fs,1)
        do kk=1,SpNum
              ! first make trapezoidal rule for 1st and last point
              psc((kk-1)*np+1)=0.5d0*psc((kk-1)*np+1) !1st
              psc(kk*np)=0.5d0*psc(kk*np) !Last
        enddo
        ! Add to power but rescale to account for integration:
        do kk=1,SpNum
           psum = ds*sum(psc((kk-1)*np+1:kk*np)) !integrate ds
           ! integrate dt and average over beats (omega/2pi)
           ps(kk) = ps(kk)+(omega/(2.d0*pi))*dt*psum 
        enddo

        !ps = ps + psc
        !if (mod(ntime,how_oftien_wr_x).eq.0) then
        !   do kk=1,SpNum
        !      call oneouts(79+kk,ps((kk-1)*np+1:kk*np),np)
        !   enddo
        !endif
        ! Every beat, output power and reset:
        if (mod(ntime,nbeat).eq.0) then
           do kk=1,SpNum
              !call oneouts(82+kk,ps((kk-1)*np+1:kk*np),np)
              call oneout(79+kk,ps(kk))
           enddo
           ! reset power to zero:
           ps = 0.d0
        endif

        return

    end subroutine addPower
!=======================================================================
!=======================================================================


end module power
