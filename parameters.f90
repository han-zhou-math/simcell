!---------------------------------------------------------------------------  
! Set constants and parameters
! MODULE: parameters  
!> @author 
!> Lingxing Yao
!
!DESCRIPTION: 
!@brief
!> setting parameters, constants, and  data structure
!>  there are two groups: I) phsyical parameters II) computing parameters;
!>   in each of the categories there are two types too: 1) fluid/eulerian grid; 2) IB/lagrangian grid and IB objects
!> 
!> @param outcount the counter for saved data, outcount*nfreq = time steps the "outcount" set data stored
!> @param ment total numbers of IB objects (IB entities)
!> @param cent total numbers of IB objects (IB entities)
!> @param nring number of IB pts on each IB entity (same for different IB entity)
!> @param mpts total numbers of IB pts, mpts = ment*nring
!> @param xpt the x-coordinates of IB pts for ALL IB objects, in shape (mpts) 
!> @param ypt the y-coordinates of IB pts for ALL IB objects, in shape (mpts) 
!--------------------------------------------------------------------------- 
!
! all parameters set here except: kc, and active pump for chemical and J_actin
!
module parameters
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  integer, parameter :: dp = selected_real_kind(15,300)
  double precision, parameter :: cpi=3.141592653589793238D0
  double precision, parameter :: zero=0.d0
  double precision, parameter :: half=0.5d0
  double precision, parameter :: one=1.d0
  double precision, parameter :: two=2.d0
  double precision, parameter :: three=3.d0
  double precision, parameter :: four=4.d0
  double precision, parameter :: eight=8.d0
  double precision :: pi, tupi
!
  integer :: outcount
  integer,parameter:: ment=1                    !max number of IB entities
  integer,parameter:: cent=1
  integer, parameter :: nib=1 ! number of IB entities
  integer, parameter :: nent=3 ! legacy variable, nent = ment+2 = nib+2
  integer,parameter:: nfil=0                    !number of points/wall
#ifndef SIMCELL_NRING
#define SIMCELL_NRING 200
#endif
  integer,parameter:: nring=SIMCELL_NRING          !number of points/platelet
  ! Dualchem calls the same material-marker count npts.  Keep one storage
  ! owner and expose an alias instead of introducing a second parameter module.
  integer,parameter:: npts=nring
!  integer,parameter:: nring=80                    !number of points/platelet
!  integer,parameter:: nring=640                    !number of points/platelet
  integer,parameter:: mpts=ment*nring             !number of IB points
  integer,parameter:: mcoor=2*mpts                !number of IB coordinates
!=======================================================================
! for Cartesian grid!
#ifndef SIMCELL_L2NY
#define SIMCELL_L2NY 7
#endif
  integer,parameter:: l2ny=SIMCELL_L2NY                    !log2 of default y grid
  integer,parameter:: maspect=1                            !aspect ratio
#ifdef SIMCELL_NY
  integer,parameter:: ny=SIMCELL_NY                        !explicit y grid override
#else
  integer,parameter:: ny=2**l2ny                           !fluid grid sizes
#endif
#ifdef SIMCELL_NX
  integer,parameter:: nx=SIMCELL_NX                        !explicit x grid override
#else
  integer,parameter:: nx=maspect*ny                        !fluid grid sizes
#endif
  !!integer,parameter:: ny=32
  !!integer,parameter:: nx=128
  integer,parameter:: nxp2=nx+1                        !fluid grid sizes
  integer,parameter:: nyp2=ny+1                        !fluid grid sizes
!
!=======================================================================
  double precision :: CFL
  ! Dimensionless Stokes coefficients:
  !   nu=mu/(P0*T0),  bke/bki=b_e/b_i*L0^2/(P0*T0).
  ! The runtime scale metadata below do not perform this conversion.
  double precision :: nu
  double precision :: bke, bki
  double precision :: vco
  double precision, parameter :: rtc = 8.4d5
  double precision, parameter :: rtr = .0500d0, rho0=rtc*rtr ! R*T*rho0,scaled with RTC
  ! here rtr = {R*T*rho0}/{R*T*c(Na+)}
!=======================================================================
#ifndef SIMCELL_DOMAIN_HALF
#define SIMCELL_DOMAIN_HALF 1.0d0
#endif
#ifndef SIMCELL_X_HALF
#define SIMCELL_X_HALF SIMCELL_DOMAIN_HALF
#endif
#ifndef SIMCELL_Y_HALF
#define SIMCELL_Y_HALF SIMCELL_DOMAIN_HALF
#endif
  ! SIMCELL_X_HALF and SIMCELL_Y_HALF allow rectangular channels while the
  ! legacy SIMCELL_DOMAIN_HALF macro remains the square-domain default.
  double precision,parameter :: xmin = -SIMCELL_X_HALF ! shared FSI/chemical domain extent
  double precision,parameter :: xmax =  SIMCELL_X_HALF
  double precision,parameter :: ymin = -SIMCELL_Y_HALF
  double precision,parameter :: ymax =  SIMCELL_Y_HALF
  !!double precision,parameter :: ymax = 0.25
  double precision,parameter :: xlength = xmax-xmin
  double precision,parameter :: ylength = ymax-ymin
  double precision,parameter :: ytop = ymax
  double precision,parameter :: ybot = ymin
  ! Coordinates are absolute: half the domain length alone is not its center
  ! when xmin or ymin is nonzero.
  double precision,parameter :: xmid = xmin+half*(xmax-xmin)
  double precision,parameter :: ymid = ymin+half*(ymax-ymin)
  double precision,parameter :: hg=xlength/dble(nx)        ! spacestep
  !!double precision :: dlt = CFL*hg
  double precision :: dlt
  !double precision,parameter :: dlt = 0.02
  integer :: ntmax ! max time steps 
  integer :: nfreq ! data save frequency
!
!=======================================================================
! IB links parameters
  double precision,parameter :: cbstiff=.000000 ! bending stiffness
  double precision,parameter :: sb(1)=cbstiff   ! bending stiffness
  ! Runtime membrane stiffness.  Keep the legacy sw array as the value passed
  ! to the IB force routines, and synchronize it after reading input.par.
#ifndef SIMCELL_CLSTIFF
#define SIMCELL_CLSTIFF .0001d0
#endif
  double precision :: clstiff=SIMCELL_CLSTIFF    ! elastic modulus
  double precision :: sw(1)=[SIMCELL_CLSTIFF]    ! elastic modulus
  double precision,parameter :: dl =0.5d0*hg                 ! IB point separation
  !!double precision,parameter :: rsl=4.0*dble(nring/80)*dl
  double precision,parameter :: rsl=zero
  !!double precision,parameter :: rsl=.16
  double precision :: ctild                  !
  double precision :: cs(1)   ! bending angle?
  double precision :: xpt(mpts),ypt(mpts)  !IB coordinates (current coor)
  double precision :: xpk(mpts),ypk(mpts)  !IB coordinates (in iteration)
  double precision :: oxpt(mpts),oypt(mpts)!IB coordinates ("old" coor)
!  double precision :: xbp(mpts),ybp(mpts)! preferred shape of cell
!=======================================================================
  character(40)   :: runname  = './Data/frun'
!=======================================================================
  double precision, parameter :: xshift=0.5d0, yshift=0.5d0 ! cell center for chemicals
  double precision, parameter :: cmsi = 1.0d0, cpsi = 1.0d0 ! initial chemicals on two sides
  double precision, parameter :: nmsi = 1.0d0, npsi = 1.0d0 ! initial network volumes
  double precision :: kw(4), kc(4), dif(4)

  ! Minimal names required by the later imported dualchem modules.  Starter
  ! remains authoritative; scalar names are deliberately prefixed so they do
  ! not collide with the existing FSI kc(4)/kw(4) arrays above.
  real(dp), parameter :: mu = 1.0_dp
  real(dp) :: dualchem_kc = 0.01_dp
  real(dp) :: dualchem_pump_start_time = -1.0_dp
  real(dp) :: dualchem_kp = -1.0_dp
  ! Negative rear width selects the legacy symmetric profile exactly.
  real(dp) :: dualchem_rear_width = -1.0_dp
  real(dp) :: dualchem_rear_amplitude_ratio = 1.0_dp
  logical :: stage14_active_actin_feedback = .true.
  real(dp) :: dualchem_pump_width = 0.21_dp*cpi
  real(dp) :: stage12_initial_network_concentration = 1.0_dp
  real(dp) :: stage12_initial_free_concentration = 1.0_dp
  logical :: enforce_fixed_timestep = .false.
  ! Uniform concentration level C_0/C_S.  For dimensional initial data
  ! c_i=C_0+G_0*x and c_e=C_0-G_0*x, the polarization input is G_0*L_0/C_S;
  ! it is added with opposite signs on the two sides of the membrane.
  real(dp) :: dualchem_initial_concentration = 0.0_dp
  real(dp) :: dualchem_initial_polarization = 0.0_dp
  ! One physical chemical diffusion is configured for the whole run.  The
  ! invalid sentinel prevents an omitted namelist entry from silently becoming
  ! the historical D=1 sample value.
  real(dp) :: dualchem_diffusion = 0.0_dp
  ! Runtime reference scales.  Numerical unknowns and governing coefficients
  ! remain nondimensional; these values are metadata for validation, logging,
  ! and physical-unit postprocessing, not a physical-parameter converter.
  real(dp) :: scale_length_um = 1.0_dp
  ! This field is the neutral-solute scale C_S.  The current input contract has
  ! no separate metadata field for the actin scale C_A.
  real(dp) :: scale_concentration_millimolar = 1.0_dp
  real(dp) :: scale_velocity_um_s = 1.0_dp
  real(dp) :: scale_time_s = 1.0_dp
  real(dp) :: scale_stress_pa = 1.0_dp
  real(dp) :: dualchem_interface_gmres_rtol = 1.0e-10_dp
  real(dp), parameter :: dualchem_interface_gmres_atol_physical = 1.0e-14_dp
  real(dp), parameter :: dualchem_mg_rtol = 1.0e-12_dp
  real(dp), parameter :: dualchem_mg_atol = 1.0e-14_dp
  integer, parameter :: dualchem_mg_max_iterations = 1000
  ! Stage 12 one-way actin validation parameters.  They are declared once in
  ! the starter-owned parameter module; all reduced coefficients
  !
  !   a=eta/(eta+eta_s),  D_n=k_sigma/(eta+eta_s)
  !
  ! and both Robin coefficients are derived by actin_model_mod.  This set is
  ! chosen to exercise turnover and membrane exchange and is not presented as
  ! a calibrated biological parameter set.
  real(dp) :: stage12_actin_eta=2.0_dp
  real(dp) :: stage12_actin_eta_s=1.0_dp
  real(dp) :: stage12_actin_k_sigma=1.8_dp
  real(dp) :: stage12_actin_dc=1.1_dp
  real(dp) :: stage12_actin_gamma=0.5_dp
  real(dp) :: stage12_actin_jc=0.24_dp
  real(dp) :: stage12_actin_dw=1.0_dp
  real(dp) :: stage12_actin_theta0=1.0_dp
  logical :: stage12_localized_polymerization=.true.
  logical :: stage12_pnas_balanced_actin_profile=.false.
  logical :: stage14_remove_uniform_actin_stress=.false.
  logical :: stage14_recenter_interface=.false.
  integer :: fsi_max_outer_iterations=1
  real(dp) :: stage14_shape_force_scale=1.0_dp
  real(dp) :: stage14_adhesion=0.0_dp
  ! Signed net x-directed external force per unit out-of-plane depth,
  ! nondimensionalized by P0*L0.  Negative values oppose +x migration.
  real(dp) :: stage14_external_load_x=0.0_dp
  ! Direct dimensionless code inputs:
  !   water_stress_mobility  = M_mu = k_w*P0/U0
  !   water_osmotic_mobility = M_C  = k_w*R*T_abs*C_S/U0.
  ! A one-permeability physical model must satisfy
  ! M_C/M_mu=R*T_abs*C_S/P0; the code also permits a generalized independent
  ! pair for constitutive sensitivity studies.
  real(dp) :: water_stress_mobility=-1.0_dp
  real(dp) :: water_osmotic_mobility=-1.0_dp
  ! True only when an old input omits one of the explicit dimensionless
  ! mobilities and readpar supplies the historical kw(1) compatibility value.
  logical :: used_legacy_water_mobility_fallback=.false.
  logical :: enable_osmotic_feedback=.false.
  real(dp) :: initial_cell_radius=0.249_dp
  real(dp) :: initial_cell_axis_ratio=1.0_dp
  real(dp) :: initial_cell_shape_factor=0.0_dp
  real(dp) :: initial_cell_center_offset_x=0.0_dp
  real(dp) :: initial_cell_center_offset_y=0.0_dp
  real(dp) :: membrane_reference_metric=1.0_dp
  real(dp) :: scc(npts)
! index of Eulerian grid
  integer, dimension(1:nx,1:ny) :: IC_Tag
  integer, dimension(-1:nx+1,-1:ny+1) :: chker0, chker1, kdf, oid, ijd
  integer :: jgd(1:nx*ny)
  integer, dimension(-1:nx+1,-1:ny+1) :: idu,idv, jdu,jdv, idf,id, iid, idn, isf, ivf
  Data kc/1.00, 1.0, 0.0, 0.0/
  Data kw/.004, 0.0, 0.0, 0.0/ !here kw is scaled by RTC to become dimensionless
  Data dif/1.0, 1.0, 0.0, 0.0/
! Boundary smooth representation 
  double precision :: mko(7,nring,cent), mkn(7,nring,cent)
  double precision :: ndv(nring,2,cent) ! normal direction at each IB point
  double precision :: tdv(nring,2,cent) ! tangent direction at each IB point
  double precision :: mtd(nring,2,cent) ! normal direction at each IB point
  double precision :: ndx(mpts) ! x component of normal direction at (xpt,ypt)
  double precision :: ndy(mpts) ! x component of normal direction at (xpt,ypt)
  double precision :: bdPoly(7,mpts) ! boundary parametraization @(xpt,ypt)
  double precision :: vc_fe(nring,nent-2)
  logical :: iskip = .false.
!=======================================================================
  double precision :: mkp(7,nring,cent)
  double precision :: cc_pc(-1:nx+1,-1:ny+1)   ! chemical
  double precision :: cc_pn(-1:nx+1,-1:ny+1)   ! chemical
!=======================================================================
  integer :: olen(cent), qlen(cent)
  double precision :: xamin,xamax,yamin,yamax      ! domaina extent
  double precision :: unc(-1:nx+1,-1:ny+1)   ! x velocity for fluid & chemical
  double precision :: vnc(-1:nx+1,-1:ny+1)   ! x velocity for fluid & chemical
!========================================================================

contains

  pure elemental logical function is_finite_run_scalar(value) result(is_finite)
    real(dp),intent(in)::value
    integer(int64)::bits

    ! Binary64 exponent inspection cannot signal on a signaling NaN.  Keep
    ! this run-input predicate shared by readpar and direct initializers so
    ! neither gate depends on logical-expression evaluation order.
    if(storage_size(value)==64 .and. radix(value)==2 .and. digits(value)==53 .and. &
        minexponent(value)==-1021 .and. maxexponent(value)==1024)then
      bits=transfer(value,bits)
      is_finite=ibits(bits,52,11)/=int(z'7ff',int64)
    else
      is_finite=ieee_is_finite(value)
    endif
  end function is_finite_run_scalar

end module parameters
