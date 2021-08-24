!************************************************
!       3D Stokes Fluid, Method of Reg Stokes
!       Multiple Sperm with surfaces and mesh
!
!************************************************
program spermFlagRS
      use parameters;       !setup of most parameter values    
      use initialization;   !setup of sperm locations
      use mesh;             !setup of mesh if any
      use eigLaPack;        !to find normal plane
      use rotation;         !to do planar rotations (need for 3D)
      use axonemeForces;    !calculate forces
      use power;    !calculate power
      use stokesfree;       !to get fluid velocities in free space
      use stokeswall;       !to get fluid velocities in half space
      use stokessphere;     !to get fluid velocities with sphere
      use outputs;          !write outputs to files
      implicit none

!-----------------------------------------------------------
!Vectors for Fluid Equations
!Immersed Boundary points (moving & tethered), 
!and associated fluid velocity, IB velocity, pressure, forces
!-----------------------------------------------------------
      real(dp),dimension(3,nptot) :: xs,fs,vs,us,xsr,fsr,fRep,fAtt
      real(dp),dimension(3,SpNum) :: normVect,normVectOld,normVectR
      real(dp),dimension(SpNum) :: ps !power 
      real(dp) :: RMinv(3,3,SpNum),xbar(3,SpNum),xstemp(3,np)
      real(dp) :: dotCheck,randNum

!-----------------------------------------------------------
! Vectors for mesh
!-----------------------------------------------------------
    real(dp), allocatable, dimension ( : ) :: xm(:,:),vm(:,:)
    real(dp), allocatable, dimension ( : ) :: fm(:,:),fRepm(:,:)
    real(dp), allocatable, dimension ( : ) :: um(:,:),spg(:,:)
    integer :: M,Ns 
    double precision :: dm

!-----------------------------------------------------------
! Vectors for tracers
!-----------------------------------------------------------
    real(dp),dimension (3,nt) :: xt,vt,ut

!----------------------------------------
!Time step and cpu time, length step/discretization, etc.
!----------------------------------------
      real(dp) :: t1,t2,tot_t
      integer :: ntime,ii,jj,kk,meth,methR,ntimestart
      character(len=128) :: stand_t
      logical :: restart
      

!-----------------------------------------------------------------------
!     INITIALIZATION
!-----------------------------------------------------------------------
   INQUIRE(FILE="fort.40", EXIST=restart)
   if (restart) then  !only works for non-mesh case right now
          call spermrestart(xs,normVect,ntimestart)
          call rotateFlag(xs,xbar,xsr,normVect,RMinv)
          !call testout('test1.txt',DBLE(ntimestart))
        ! NEED TO DO:
        ! meshrestart(x,y,z,spg)
   else
      ntimestart = 0
      call cpu_time(t1)

!     ------------------------------------------------
!     Initialize mesh
!     ------------------------------------------------
      dm = ep ! if not sphere then use same blob size as for sperm
      M = 0;
      if (meshtype==3) then
         call meshsphere(xm,fm,fRepm,vm,um,spg,M,Ns,dm)
         dm = ep
      elseif (meshtype.gt.0) then
         call meshcube(xm,fm,fRepm,vm,um,spg,M,Ns)
      endif

      call paramsout('setup.txt','plt_data.m',restart,dm,M,Ns);
!     ------------------------------------------------
!     Initialize tracers
!     ------------------------------------------------
      !if (nt>0) then !!! NEED TO CREATE
      !  call tracerinit(xt) 
      !endif

!     ------------------------------------------------
!     Initialize sperm
!     ------------------------------------------------
      xs = 0.d0; fRep=0.d0; ps = 0.d0
      ! Initialize from sinusoidal curve:
      call sperminit(xs,2,3)
      ! Initialize from data for a single swimmer: choose 1-40 for second arg
      !call spermFromData1Sp(xs,1,'FreePos2pi',2)
      !xs(3,:) = -xs(3,:)     
      !call switchDim(xs,3,2)
      !call rotateInit(xs,theta) 
      do kk=1,SpNum
         xstemp = xs(:,1+(kk-1)*np:kk*np)
         call rotateData(xstemp,1,0.5d0*pi+sizeRot*DBLE(kk-1))
         call RANDOM_NUMBER(randNum)
         call rotateData(xstemp,1,sizeRand*(randNum-0.5d0))
         call RANDOM_NUMBER(randNum)
         call rotateData(xstemp,2,sizeRand*(randNum-0.5d0))
         call RANDOM_NUMBER(randNum)
         call rotateData(xstemp,3,sizeRand*(randNum-0.5d0))
         xs(:,1+(kk-1)*np:kk*np) = xstemp
      enddo 
!     ------------------------------------------------------------------
!     Rotate spatial data so that flagellum planes are x-y
!     ------------------------------------------------------------------
      do kk=1,SpNum
         call normalPlane(xs(:,1+(kk-1)*np:kk*np),normVect(:,kk),meth)
      enddo
      do kk=2,SpNum
         dotCheck = SUM(normVect(:,kk-1)*normVect(:,kk))
         if (dotCheck.gt.0.d0) then
            normVect(:,kk) = -normVect(:,kk)
         endif
      enddo  
      call rotateFlag(xs,xbar,xsr,normVect,RMinv)

!     ------------------------------------------------------------------
!     Output Initialized quantities
!     ------------------------------------------------------------------
       if (meshtype>0) then   ! Write out mesh positions
          call threeout(13,xm,M)
       endif
       if (nt>0) then ! write out tracer positions
          call threeout(14,xt,nt)  
       endif
      do kk=1,SpNum
          call threeout(39+kk,xs(:,1+(kk-1)*np:kk*np),np)  !positions
       !   call threeout(49+kk,xsr(:,1+(kk-1)*np:kk*np),np) !Lagrangian pos
          ! Write out normal plane vector: (if 3D)
          call threeout(89+kk,normVect(:,kk),1)
      end do
   endif


!=======================================================================
!     MAIN TIME INTEGRATION LOOP      
!=======================================================================
      main_loop: do ntime = ntimestart+1,ntimestart+tot_time_steps  

!        ---------------------------------------------------------------
!        Calculating the forces along the flagellum
!        ---------------------------------------------------------------
         call spermF(xsr,fsr,ntime)
!        ---------------------------------------------------------------
         !Rotate back:
         call rotateBack(fsr,fs,RMinv)
         
!        ---------------------------------------------------------------
!        Add Head attachment forces:
!        ---------------------------------------------------------------
         call attachFlag(xs,fAtt)
         fs = fs + fAtt
!        ---------------------------------------------------------------
!        Add repulsive forces:
!        ---------------------------------------------------------------
         call repulsion(xs,fRep)
         fs = fs + fRep
         call repulsionFlagMesh(xs,xm,fRep,fRepm,M)
         fs = fs + fRep
         fm = fm + fRepm         

!        ---------------------------------------------------------------
!        SOLVE Regularized Stokeslets EQNS
!        ---------------------------------------------------------------
         if (meshtype==0) then
            if (wallExist) then
               call StokesletWall(xs,us,xs,fs,nptot,nptot,ep)
            elseif (sphereExist) then
               call StokesletSphere(xs,us,xs,fs,nptot,nptot,ep)   
            else
               call StokesletFree(xs,us,xs,fs,nptot,nptot,ep)
            endif
            xs = xs + dt*us 
            if (nt>0) then! get sperm velocities resulting from all forces
               if (wallExist) then
                  call StokesletWall(xt,ut,xs,fs,nt,nptot,ep)
               elseif (sphereExist) then
                  call StokesletSphere(xt,ut,xs,fs,nt,nptot,ep)   
               else
                  call StokesletFree(xt,ut,xs,fs,nt,nptot,ep)
               endif
               xt = xt + dt*ut
            endif
            
        
!     ------------------------------------------------------------------
         else
!       Stokeslet evaluation split due to different blob size
!       on the mesh and the swimmer(s)
!       Also different time scale for mesh
            if (mod(ntime,slow_t).eq.0) then
                call meshForces(xm,fm,spg,Ns,M)
                ! get mesh velocities resulting from all forces
                if (wallExist) then
                   call StokesletWall(xm,um,xs,fs,M,nptot,ep)
                   call StokesletWall(xm,vm,xm,fm,M,M,dm)
                elseif (sphereExist) then
                   call StokesletSphere(xm,um,xs,fs,M,nptot,ep)   
                   call StokesletSphere(xm,vm,xm,fm,M,M,dm)   
                else
                   call StokesletFree(xm,um,xs,fs,M,nptot,ep)
                   call StokesletFree(xm,vm,xm,fm,M,M,dm)
                endif
                !-----------------------------------
                !        Update mesh Positions
                !-----------------------------------
                xm = xm + (vm+um)*slow_t*dt
            endif      
            ! get sperm velocities resulting from all forces
            if (wallExist) then
               call StokesletWall(xs,us,xs,fs,nptot,nptot,ep)
               call StokesletWall(xs,vs,xm,fm,nptot,M,dm)
            elseif (sphereExist) then
               call StokesletSphere(xs,us,xs,fs,nptot,nptot,ep)   
               call StokesletSphere(xs,vs,xm,fm,nptot,M,dm)   
               ! call testout('test3.txt',1.d0)
            else
               call StokesletFree(xs,us,xs,fs,nptot,nptot,ep)
               call StokesletFree(xs,vs,xm,fm,nptot,M,dm)
            endif
                !-----------------------------------
                !        Update Sperm Positions
                !-----------------------------------
            xs = xs + dt*(us+vs) 

         !---------------------------------------------------------------
         !        Update Positions of the tracers
         !---------------------------------------------------------------
            if (nt>0) then! get sperm velocities resulting from all forces
               if (wallExist) then
                  call StokesletWall(xt,ut,xs,fs,nt,nptot,ep)
                  call StokesletWall(xt,vt,xm,fm,nt,M,dm)
               elseif (sphereExist) then
                  call StokesletSphere(xt,ut,xs,fs,nt,nptot,ep)   
                  call StokesletSphere(xt,vt,xm,fm,nt,M,dm)   
               else
                  call StokesletFree(xt,ut,xs,fs,nt,nptot,ep)
                  call StokesletFree(xt,vt,xm,fm,nt,M,dm)
               endif
               xt = xt + dt*(ut+vt)
            endif
            
        endif !All possible stokes solving done

!     ------------------------------------------------------------------
!     Find power 
!     ------------------------------------------------------------------
       call addPower(us,vs,fs,ps,ntime) 

!     ------------------------------------------------------------------
!     Finding normal planes: 
!     ------------------------------------------------------------------
         normVectOld = normVect
         do kk=1,SpNum
            call normalPlane(xs(:,1+(kk-1)*np:kk*np),normVect(:,kk),meth)
         end do
         !CHECK TO SEE IF ORIENTATION IS CORRECT:
         do kk=1,SpNum
            dotCheck = SUM(normVectOld(:,kk)*normVect(:,kk))
            if (dotCheck.lt.0.d0) then
               normVect(:,kk) = -normVect(:,kk)
            endif
         enddo

         ! New rotated data
         call rotateFlag(xs,xbar,xsr,normVect,RMinv)

!        ---------------------------------------------------------------
!        OUTPUT SECTION
!        --------------------------------------------------------------- 
         if (mod(ntime,how_often_wr_x).eq.0) then
                if ((meshtype>0).AND.(mod(ntime,mwrite).eq.0)) then ! Write out mesh positions
                   call threeout(13,xm,M)
                endif
                if (nt>0) then ! Write out tracer positoons
                    call threeout(14,xt,nt)  
                endif
                do kk=1,SpNum
                    !write out the forces on IB
                    call threeout(29+kk,fs(:,1+(kk-1)*np:kk*np),np)
                    !call threeout(09+kk,fsr(:,1+(kk-1)*np:kk*np),np)
                    !write out the rep forces on IB
                    call threeout(19+kk,fRep(:,1+(kk-1)*np:kk*np),np)
                    call threeout(49+kk,fAtt(:,1+(kk-1)*np:kk*np),np)
                    !write out the Flagellum/IB/filament location
                    call threeout(39+kk,xs(:,1+(kk-1)*np:kk*np),np)
                    !call threeout(49+kk,xsr(:,1+(kk-1)*np:kk*np),np) !Lagrangian pos
                    call threeout(59+kk,us(:,(kk-1)*np+1:kk*np),np)
                    call threeout(69+kk,vs(:,(kk-1)*np+1:kk*np),np)
                    ! Write out normal plane vector: (if 3D)
                    call threeout(89+kk,normVect(:,kk),1)
               end do

               open(99,file='progress.txt',status='replace');
               write(99,*) ntime; 
               write(99,*) ntime/how_often_wr_x 
               close(99);
               !call cpu_time(t2);   tot_t = t2 - t1;
               !call std_time('run-time.txt',tot_t,stand_t) !record time
                
         end if
      end do main_loop !END OF MAIN TIME INTEGRATION LOOP
!!-----------------------------------------------------------------------
      call cpu_time(t2);   tot_t = t2 - t1;
      call std_time('run-time.txt',tot_t,stand_t) !record time
!-------------------------------------------------------------------------

end program spermFlagRS
