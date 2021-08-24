module outputs
!***********************************************************************
!***********************************************************************
   implicit none
   private :: i2ch
   contains
!=======================================================================
!Writing to file fort.* the value of a single number
!=======================================================================
   subroutine oneout(port,x)
!     this subroutine outputs x ==> fort. 
      use parameters, only: dp
      implicit none
      integer,intent(in) :: port
      real(dp),intent(in) :: x(1)
!     ------------------------------------------------------------------
1     format(2(e23.15))   
      write(port,1) x(1)
      return
   end subroutine oneout
!=======================================================================
!Writing to file fort.* the value of several single number
!=======================================================================
   subroutine oneouts(port,x,M)
!     this subroutine outputs x ==> fort. 
      use parameters, only: dp
      implicit none
      integer,intent(in) :: port,M
      real(dp),intent(in) :: x(M)
      integer :: ip
!     ------------------------------------------------------------------
1     format(2(e23.15))   
      do ip = 1,M
         write(port,1) x(ip)
      enddo
      return
   end subroutine oneouts
!=======================================================================
!Writing to file fort.* the values of x(1:2)
!=======================================================================
!Writing to file fort.* the values of x(1:2)
!=======================================================================
   subroutine twoout(port,x,M)
!     this subroutine outputs x ==> fort. 
      use parameters, only: dp
      implicit none
      integer,intent(in) :: port,M
      real(dp),intent(in) :: x(2)
      integer :: ip
!     ------------------------------------------------------------------
1     format(2(e23.15))   
      do ip = 1,M
         write(port,1) x(1), x(2)
      enddo
      return
   end subroutine twoout
!=======================================================================
!Writing to file fort.* the values of xm(1:3,:)
!=======================================================================
   subroutine threeout(port,xm,M)
!     this subroutine outputs x ==> fort. 
      use parameters, only: dp
      implicit none
      integer,intent(in) :: port,M
      real(dp),intent(in) :: xm(3,M)
      integer :: ip
      CHARACTER(LEN=20) :: str
!     ------------------------------------------------------------------

       write(str,"(a,i2)") "fort.",port
       open(unit=4001,file=str, access = 'append', action='write')
       do ip = 1,M
          write(4001,*) xm(1,ip), xm(2,ip), xm(3,ip)
       end do
       return
   end subroutine threeout
!=======================================================================
!=======================================================================
!Writing to file fort.* the values of xm(1:3,:)
!=======================================================================
   subroutine threeoutLong(port,xm,M)
!     this subroutine outputs x ==> fort. 
      use parameters, only: dp

      implicit none
      integer,intent(in) :: port,M
      real(dp),intent(in) :: xm(3,M)
      integer :: ip
      CHARACTER(LEN=20) :: str
!     ------------------------------------------------------------------

       write(str,"(a,i3)") "fort.",port
       open(unit=4001,file=str, access = 'append', action='write')
       do ip = 1,M
          write(4001,*) xm(1,ip), xm(2,ip), xm(3,ip)
       end do
       return
   end subroutine threeoutLong
!=======================================================================
!!=======================================================================
   subroutine paramsout(filename1,filename2,restart,dm,M,Ns)
      use parameters
      implicit none
      character(len=*) :: filename1,filename2
      logical :: restart
      integer :: M,Ns
      real(dp) :: dm

     open(10,file=filename1,status='replace');
     write(10,*) 'dt(time step) =  ',dt
     write(10,*) 'S1 = ',S1
     write(10,*) 'S2 = ',S2
     write(10,*) 'np = ',np
     write(10,*) 'nmax = ',tot_time_steps
     write(10,*) 'nwrite = ',how_often_wr_x
     write(10,*) 'mwrite = ',mwrite
     write(10,*) 'mfactor = ',mfactor
     write(10,*) 'restart =',restart
     write(10,*)
     write(10,*) '% network parameters - # layers, # nodes per layer,'
     write(10 ,*) '% link stiffnes and viscosity, sphere radius, mesh blob size'
     write(10 ,*) 'NL  =   ',NL
     write(10 ,*) 'N   =   ',N
     write(10 ,*) 'S   =  ',S
     write(10 ,*) 'eta =  ',eta
     write(10 ,*) 'srad   =  ', srad
     write(10 ,*) 'delta =',dm
     close(10)

     open(11,file=filename2,status='replace')
     write(11,*) 'np =',np,';'
     write(11,*) 'nwrite =',how_often_wr_x,';'
     write(11,*) 'mwrite =',mwrite,';'
     write(11,*) 'mfactor =',mfactor,';'
     write(11,*) 'SpNum=',SpNum,';'
     write(11,*) 'dt=',dt,';'
     if (sphereExist) then
         write(11,*) 'sphereExist = 1;'
     else
         write(11,*) 'sphereExist = 0;'
     endif
     if (wallExist) then
         write(11,*) 'wallExist = 1;'
     else
         write(11,*) 'wallExist = 0;'
     endif
     if (meshtype>0) then
         write(11,*) 'nm = ',M,';'
         write(11,*) 'ns = ',Ns,';'
     endif
     write(11,*) 'wallLoc = ',wallLoc,';'
     write(11,*) 'fid = fopen(''progress.txt'',''r'');'
     write(11,*) 'tot_timesteps = fscanf(fid,''%d'',1);'
     write(11,*) 'fclose(fid);'
     close(11)
   end subroutine paramsout
!=======================================================================
!=======================================================================
   subroutine std_time(fname,tot_time,stand_t,dd,hh,mm,ss)
   ! this subroutine convert time(in se) into standard form (ddhhmmss)
      use parameters, only: dp; 
      implicit none
      character(len=*),intent(in) :: fname
      real(dp),intent(in) :: tot_time
      character(len=11),intent(out) :: stand_t
      integer,intent(out),optional :: ss,mm,hh,dd
      real :: tm,th,td,ts;  integer :: s,m,h,d

      td = tot_time/3600./24.;               d = int(td);
      th = (tot_time - d*24.*3600)/3600.;    h = int(th);
      tm = (tot_time - (d*24.+h)*3600.)/60.; m = int(tm);
      ts = tot_time - ((d*24.+h)*60.+ m)*60; s = int(ts);
      stand_t = i2ch(d)//':'//i2ch(h)//':'//i2ch(m)//':'//i2ch(s)
      open(55,file=fname,status='replace')
      write(55,*) 'Time elapsed(dd:hh:mm:ss): ',stand_t; close(55);
      if(present(dd)) dd = d; if(present(hh)) hh = h;
      if(present(mm)) mm = m; if(present(ss)) ss = s;
   end subroutine std_time
!=======================================================================
   function i2ch(mm) result(ch)
      implicit none
      integer,intent(in) :: mm; character(len=2) :: ch
      ch = '';
      if (mm<0 .or. mm >100) then
         print*,'input number is out of range'; return
      end if
      ch = char(mm/10+48)//char(mm-(mm/10)*10+48)
   end function i2ch
!=======================================================================
   subroutine testout(filename1,msg1)
      use parameters, only: dp; 
      implicit none
      character(len=*) :: filename1
      real(dp) :: msg1

     open(11,file=filename1);
     write(11,*) 'got here, ', msg1
     close(11)
   end subroutine testout

end module outputs
