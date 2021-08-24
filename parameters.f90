module parameters

  implicit none
!----------------------------------------------------------------------
!general parameters: time steps, setting up, etc.
!----------------------------------------------------------------------
  integer,parameter :: sp=kind(1.0),dp=kind(1.d0) !selected_real_kind(30,32)
  real(dp),parameter :: pi=3.14159265358979323846264338327950288;
  real(dp),parameter :: mu=1.0d0
  real(dp),parameter :: dt=5.0d-7
  integer,parameter :: tot_time_steps = 2.0d+8
  integer,parameter :: how_often_wr_x = 5.0d+4
  !integer,parameter :: how_often_wr_x = 5.0d+4
  integer,parameter :: mfactor = 1.d+3
  integer,parameter :: mwrite = how_often_wr_x*mfactor
  integer,parameter :: nbeat = int(1.d0/dt)

!----------------------------------------------------------
!Sperm parameters:
!geometry, length, configuration, number, stiffness
!--------------------------------------------
   integer, parameter :: SpNum=2
! motType: 0 = sym, 1 = slight asym, 2 = asym
   integer, parameter :: motType=0
! penalty: 0 = none, 1 = z-only, 2 = curvature 3D, 3 = local z penalty
   integer, parameter :: penalty=3

! DOMAIN PARAMETERS, FLAGELLUM PARAMETERS:
   integer, parameter :: np=101
   integer, parameter :: nptot=np*SpNum
   real(dp),parameter :: npd = DBLE(np)
   real(dp),parameter :: cr=1.d0 !length of sperm
   real(dp),parameter :: omega=2.0d0*pi,Kp=2.0d0*pi/cr,b=0.1d0*cr
   real(dp),parameter :: ds=cr/DBLE(np-1) !spacing between points
   real(dp),parameter :: SpDist=4.d-2,repDist=3.d0*ds
   real(dp),parameter :: SpStart = 1.d0 !init location head (in x)
   !regularization parameter
   real(dp),parameter :: ep=cr/(3.d0*npd/4.d0)
   ! for rotation of initial positions:
   real(dp),parameter :: sizeRot = 2.d0*pi/SpNum 
   real(dp),parameter :: sizeRand = 0.d0
   real(dp),parameter :: rotBreak = -2.d0,theta=0.d0

! STIFFNESS PARAMETERS:
!s1=spring (inextensibility), s2=bending, s3=repulsive
   real(dp),parameter :: S1=1.0d+2,S2=5.0d-2,S3=5.0d-2
   !real(dp),parameter :: S1=1.0d+2,S2=0.05d0,S3=0.05d0
! bonding/attach parameters:
   real(dp),parameter :: S4=1.0d+2,attDist=6.d-2,attThresh=1.d-1
!s5=bending in z, s6=0 rest length spring in z, s7=local z penalty
   real(dp),parameter :: S6 = 1.0d0, S5=1.3d-1, S7 =5.0d+1
   !real(dp),parameter :: S6 = 1.d0, S5=1.3d-1, S7 = 50.d0

!--------------------------------------------
! Mesh parameters
!--------------------------------------------
    ! viscoelastic network parameters
    ! Meshtype: 0 (none), 1 (simple cube), 2 (diag cube), 3 (sphere)
    integer,parameter :: meshtype = 0
    integer, parameter :: slow_t = 1000
    ! Cube nodes per side:
    integer, parameter :: Ncube = 10
    ! Sphere number of layers, nodes per layer, slower time scale factor
    integer, parameter :: NL = 5, N = 1002
    ! stiffnes, viscosity, sphere radius
    !double precision, parameter :: S = 10.D0, eta = 0.D0, srad = 1.D0, El = 10.d0
    double precision, parameter :: S = 10.D0, eta = 0.D0, srad = 0.5D0, El = 10.d0
    double precision, parameter :: aa=0.5d0,pert=0.d0
    double precision, parameter :: hh = 2.d0*aa/dble(Ncube-1)
    double precision, parameter :: shift = 0.d0 !hh/2.d0
    integer, parameter :: rot = 0 !1

!--------------------------------------------
! Surface parameters       
!--------------------------------------------
    logical,parameter :: sphereExist = .false.
    logical,parameter :: wallExist = .false.
    double precision, parameter :: wallLoc = -0.5d0 
 
!--------------------------------------------
! Tracer points       
!--------------------------------------------
   integer, parameter :: nt=0
        
!==========================================================================
end module parameters
