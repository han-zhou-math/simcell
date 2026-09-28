! Stage 11 fully coupled actin boundary residual.
!
! This module is an orchestration layer, not a new Cartesian PDE solver.  Each
! scalar extension is still computed by the Stage 10 evaluator, which reaches
! the unchanged Stage 06 semi-periodic C++ multigrid backend.  The only new
! object here is the packed interface-density residual
!
!        z = [ psi_n ; psi_c ],
!
! where n denotes network/F-actin and c denotes free/G-actin.  For one
! residual callback we first solve the network extension, then use that SAME
! trial field in the free-actin volume source gamma*w_n.  With the physical
! interior identified as the plus side and [d_n w_j]=psi_j, the boundary
! equations are
!
!   F_n = psi_n/2 + {d_n w_n} + beta_n {w_n} - q {w_c},
!   F_c = psi_c/2 + {d_n w_c} + beta_c {w_c}.
!
! The marker source jump for gamma*w_n is zero because the auxiliary network
! value jump is [w_n]=0 on the current interface.  This is distinct from the
! old-network swept source used by the sequential algorithm.
module fsi_actin_block_harness_mod
  use parameters, only: dp
  use scalar_operator_mod, only: scalar_rhs_t,SCALAR_STATUS_OK, &
       validate_scalar_operator,is_finite_scalar
  use generic_gmres_mod, only: gmres_result_t,gmres_options_t, &
       GMRES_STATUS_SUCCESS,solve_affine_gmres
  use fsi_actin_sequential_harness_mod, only: ACTIN_SPECIES_NETWORK, &
       ACTIN_SPECIES_FREE,actin_extension_trial_t,actin_species_state_t, &
       actin_sequential_state_t,actin_sequential_problem_t, &
       actin_species_evaluator
  implicit none
  private

  integer,parameter,public::ACTIN_BLOCK_OK=0
  integer,parameter,public::ACTIN_BLOCK_INVALID=20
  integer,parameter,public::ACTIN_BLOCK_EVALUATOR_FAILED=21
  integer,parameter,public::ACTIN_BLOCK_SOLVE_FAILED=22
  integer,parameter,public::ACTIN_BLOCK_FINAL_FAILED=23

  type,public::actin_block_result_t
    integer::status=ACTIN_BLOCK_INVALID
    integer::evaluator_status=0
    type(gmres_result_t)::gmres
    real(dp),allocatable::final_residual(:)
  end type actin_block_result_t

  public::evaluate_actin_block_residual
  public::advance_actin_block

contains

  subroutine evaluate_actin_block_residual(problem,state,packed_density, &
       evaluate_species,evaluator_context,packed_residual,network_trial, &
       free_trial,status)
    type(actin_sequential_problem_t),intent(in)::problem
    type(actin_sequential_state_t),intent(in)::state
    real(dp),intent(in)::packed_density(:)
    procedure(actin_species_evaluator)::evaluate_species
    class(*),intent(inout)::evaluator_context
    real(dp),intent(out)::packed_residual(:)
    type(actin_extension_trial_t),intent(out)::network_trial,free_trial
    integer,intent(out)::status
    type(scalar_rhs_t)::network_rhs,free_rhs
    integer::nmarker,nx_local,ny_local,allocation_status,evaluator_status
    integer::marker,i,j
    real(dp)::source_value,term,value,next_value
    logical::ok

    ! Failure is atomic from the caller's perspective: the residual is zero,
    ! accepted state/problem are intent(in), and disposable trials are not
    ! published as valid unless status becomes ACTIN_BLOCK_OK.
    packed_residual=0.0_dp
    status=ACTIN_BLOCK_INVALID
    if(.not.valid_block_inputs(problem,state,packed_density,packed_residual, &
         nx_local,ny_local,nmarker))return

    allocate(network_rhs%volume(nx_local,ny_local), &
         network_rhs%source_jump(nmarker),free_rhs%volume(nx_local,ny_local), &
         free_rhs%source_jump(nmarker),stat=allocation_status)
    if(allocation_status/=0)return
    network_rhs%volume=0.0_dp
    network_rhs%source_jump=0.0_dp

    ! First member of the packed vector is the current network density jump.
    ! Its bulk operator already contains the implicit turnover reaction gamma.
    call evaluate_species(ACTIN_SPECIES_NETWORK, &
         packed_density(1:nmarker),problem%network_operator,network_rhs, &
         problem%dt,state%network%field, &
         state%network%correction_coefficients,network_trial, &
         evaluator_context,evaluator_status)
    if(evaluator_status/=0)then
      status=ACTIN_BLOCK_EVALUATOR_FAILED
      return
    endif
    if(.not.valid_trial(network_trial,nx_local,ny_local,nmarker))then
      status=ACTIN_BLOCK_EVALUATOR_FAILED
      return
    endif

    ! Fully implicit coupling: gamma*w_n^{m+1} is rebuilt from this callback's
    ! network trial, not from the accepted old field and not from a separately
    ! cached solve.  Checked multiplication makes extreme finite input fail
    ! before an IEEE overflow can escape the residual oracle.
    free_rhs%source_jump=0.0_dp
    do j=1,ny_local
      do i=1,nx_local
        call checked_multiply(problem%coefficients%turnover_rate, &
             network_trial%field(i,j),source_value,ok)
        if(.not.ok)return
        free_rhs%volume(i,j)=source_value
      enddo
    enddo
    call evaluate_species(ACTIN_SPECIES_FREE, &
         packed_density(nmarker+1:2*nmarker),problem%free_operator,free_rhs, &
         problem%dt,state%free%field,state%free%correction_coefficients, &
         free_trial,evaluator_context,evaluator_status)
    if(evaluator_status/=0)then
      status=ACTIN_BLOCK_EVALUATOR_FAILED
      return
    endif
    if(.not.valid_trial(free_trial,nx_local,ny_local,nmarker))then
      status=ACTIN_BLOCK_EVALUATOR_FAILED
      return
    endif

    ! Assemble one marker at a time so every product/sum has fail-closed
    ! arithmetic.  The -q*{w_c} term is the off-diagonal membrane coupling in
    ! the network Robin condition; the free condition has no network trace.
    do marker=1,nmarker
      call checked_multiply(0.5_dp,packed_density(marker),value,ok)
      if(.not.ok)return
      call checked_add(value,network_trial%traces%normal_average(marker), &
           next_value,ok)
      if(.not.ok)return
      value=next_value
      call checked_multiply(problem%beta_network(marker), &
           network_trial%traces%value_average(marker),term,ok)
      if(.not.ok)return
      call checked_add(value,term,next_value,ok)
      if(.not.ok)return
      value=next_value
      call checked_multiply(-problem%coupling_q(marker), &
           free_trial%traces%value_average(marker),term,ok)
      if(.not.ok)return
      call checked_add(value,term,packed_residual(marker),ok)
      if(.not.ok)then; packed_residual=0.0_dp; return; endif

      call checked_multiply(0.5_dp,packed_density(nmarker+marker),value,ok)
      if(.not.ok)then; packed_residual=0.0_dp; return; endif
      call checked_add(value,free_trial%traces%normal_average(marker), &
           next_value,ok)
      if(.not.ok)then; packed_residual=0.0_dp; return; endif
      value=next_value
      call checked_multiply(problem%beta_free(marker), &
           free_trial%traces%value_average(marker),term,ok)
      if(.not.ok)then; packed_residual=0.0_dp; return; endif
      call checked_add(value,term,packed_residual(nmarker+marker),ok)
      if(.not.ok)then; packed_residual=0.0_dp; return; endif
    enddo
    if(.not.all(is_finite_scalar(packed_residual)))then
      packed_residual=0.0_dp
      return
    endif
    status=ACTIN_BLOCK_OK
  end subroutine evaluate_actin_block_residual

  subroutine advance_actin_block(problem,state,evaluate_species, &
       evaluator_context,status,result)
    type(actin_sequential_problem_t),intent(in)::problem
    type(actin_sequential_state_t),intent(inout)::state
    procedure(actin_species_evaluator)::evaluate_species
    class(*),intent(inout)::evaluator_context
    integer,intent(out)::status
    type(actin_block_result_t),intent(out),optional::result
    type(actin_block_result_t)::work_result
    type(actin_sequential_state_t)::staged_state
    type(actin_extension_trial_t)::network_trial,free_trial
    real(dp),allocatable::packed_density(:),validation_residual(:)
    integer::nmarker,nx_local,ny_local,allocation_status,final_status

    status=ACTIN_BLOCK_INVALID
    work_result%status=ACTIN_BLOCK_INVALID
    work_result%evaluator_status=0
    if(.not.allocated(state%network%density))then
      if(present(result))result=work_result
      return
    endif
    nmarker=size(state%network%density)
    allocate(packed_density(2*nmarker),validation_residual(2*nmarker), &
         work_result%final_residual(2*nmarker),stat=allocation_status)
    if(allocation_status/=0)then
      if(present(result))result=work_result
      return
    endif
    packed_density=0.0_dp
    validation_residual=0.0_dp
    work_result%final_residual=0.0_dp
    if(.not.valid_block_inputs(problem,state,packed_density, &
         validation_residual,nx_local,ny_local,nmarker))then
      if(present(result))result=work_result
      return
    endif

    ! Generic GMRES only sees the affine map F(z).  Its matrix action is
    ! F(p)-F(0), so neither the species evaluator nor this harness assembles a
    ! dense 2N-by-2N boundary matrix.
    call solve_affine_gmres(block_callback,evaluator_context,packed_density, &
         problem%gmres_options,work_result%gmres)
    if(.not.work_result%gmres%converged .or. &
         work_result%gmres%status/=GMRES_STATUS_SUCCESS)then
      work_result%status=ACTIN_BLOCK_SOLVE_FAILED
      work_result%evaluator_status=work_result%gmres%callback_status
      status=work_result%status
      if(present(result))result=work_result
      return
    endif

    ! Re-evaluate the converged packed density once to obtain a matched pair
    ! of disposable trials.  We never commit whatever trial happened to be
    ! left by the last Arnoldi callback.
    call evaluate_actin_block_residual(problem,state,packed_density, &
         evaluate_species,evaluator_context,work_result%final_residual, &
         network_trial,free_trial,final_status)
    if(final_status/=ACTIN_BLOCK_OK .or. .not.final_residual_is_accepted( &
         work_result%final_residual,work_result%gmres, &
         problem%gmres_options))then
      work_result%status=ACTIN_BLOCK_FINAL_FAILED
      work_result%evaluator_status=final_status
      status=work_result%status
      if(present(result))result=work_result
      return
    endif

    staged_state=state
    call stage_species(staged_state%network,packed_density(1:nmarker), &
         network_trial)
    call stage_species(staged_state%free,packed_density(nmarker+1:2*nmarker), &
         free_trial)
    staged_state%geometry=problem%next_geometry

    ! This is the sole accepted-state write.  Every invalid, callback, Krylov,
    ! or final-pair failure above returns with both old species intact.
    state=staged_state
    work_result%status=ACTIN_BLOCK_OK
    work_result%evaluator_status=0
    status=ACTIN_BLOCK_OK
    if(present(result))result=work_result

  contains
    subroutine block_callback(density,residual,context,callback_status)
      real(dp),intent(in)::density(:)
      real(dp),intent(out)::residual(:)
      class(*),intent(inout)::context
      integer,intent(out)::callback_status
      type(actin_extension_trial_t)::disposable_network,disposable_free
      call evaluate_actin_block_residual(problem,state,density, &
           evaluate_species,context,residual,disposable_network, &
           disposable_free,callback_status)
    end subroutine block_callback
  end subroutine advance_actin_block

  subroutine stage_species(species,density,trial)
    type(actin_species_state_t),intent(inout)::species
    real(dp),intent(in)::density(:)
    type(actin_extension_trial_t),intent(in)::trial
    species%field=trial%field
    species%density=density
    species%correction_coefficients=trial%correction_coefficients
    species%correction_grid=trial%correction_grid
    species%traces=trial%traces
  end subroutine stage_species

  logical function final_residual_is_accepted(residual,gmres,options)
    real(dp),intent(in)::residual(:)
    type(gmres_result_t),intent(in)::gmres
    type(gmres_options_t),intent(in)::options
    real(dp)::residual_norm,relative_threshold,threshold
    logical::ok
    final_residual_is_accepted=.false.
    residual_norm=stable_norm(residual,ok)
    if(.not.ok)return
    call checked_multiply(options%relative_tolerance, &
         gmres%initial_residual_norm,relative_threshold,ok)
    if(.not.ok)relative_threshold=huge(1.0_dp)
    threshold=max(options%absolute_tolerance,relative_threshold, &
         gmres%absolute_residual)+1024.0_dp*epsilon(1.0_dp)
    final_residual_is_accepted=residual_norm<=threshold
  end function final_residual_is_accepted

  real(dp) function stable_norm(values,ok) result(norm_value)
    real(dp),intent(in)::values(:)
    logical,intent(out)::ok
    real(dp)::scale,ssq,magnitude,ratio,root_ssq
    integer::i
    scale=0.0_dp; ssq=1.0_dp; norm_value=0.0_dp; ok=.false.
    if(.not.all(is_finite_scalar(values)))return
    do i=1,size(values)
      magnitude=abs(values(i))
      if(magnitude/=0.0_dp)then
        if(scale<magnitude)then
          ratio=scale/magnitude
          ssq=1.0_dp+ssq*ratio*ratio
          scale=magnitude
        else
          ratio=magnitude/scale
          ssq=ssq+ratio*ratio
        endif
      endif
    enddo
    if(scale==0.0_dp)then; norm_value=0.0_dp; ok=.true.; return; endif
    root_ssq=sqrt(ssq)
    if(scale>huge(1.0_dp)/root_ssq)return
    norm_value=scale*root_ssq
    ok=is_finite_scalar(norm_value)
  end function stable_norm

  logical function valid_block_inputs(problem,state,packed_density, &
       packed_residual,nx_local,ny_local,nmarker)
    type(actin_sequential_problem_t),intent(in)::problem
    type(actin_sequential_state_t),intent(in)::state
    real(dp),intent(in)::packed_density(:),packed_residual(:)
    integer,intent(out)::nx_local,ny_local,nmarker
    valid_block_inputs=.false.; nx_local=0; ny_local=0; nmarker=0
    if(.not.allocated(state%network%field) .or. &
         .not.allocated(state%free%field))return
    nx_local=size(state%network%field,1); ny_local=size(state%network%field,2)
    if(nx_local<=0 .or. ny_local<=0)return
    if(any(shape(state%free%field)/=[nx_local,ny_local]))return
    if(.not.allocated(state%network%density) .or. &
         .not.allocated(state%free%density))return
    nmarker=size(state%network%density)
    if(nmarker<=0 .or. size(state%free%density)/=nmarker)return
    if(size(packed_density)/=2*nmarker .or. &
         size(packed_residual)/=2*nmarker)return
    if(.not.all(is_finite_scalar(packed_density)))return
    if(.not.valid_species(state%network,nx_local,ny_local,nmarker))return
    if(.not.valid_species(state%free,nx_local,ny_local,nmarker))return
    if(validate_scalar_operator(problem%network_operator,nx_local,ny_local, &
         nmarker)/=SCALAR_STATUS_OK)return
    if(validate_scalar_operator(problem%free_operator,nx_local,ny_local, &
         nmarker)/=SCALAR_STATUS_OK)return
    if(.not.is_finite_scalar(problem%dt))return
    if(problem%dt<=0.0_dp)return
    if(.not.is_finite_scalar(problem%coefficients%turnover_rate))return
    if(problem%coefficients%turnover_rate<0.0_dp)return
    if(.not.allocated(problem%coupling_q))return
    if(size(problem%coupling_q)/=nmarker)return
    if(.not.all(is_finite_scalar(problem%coupling_q)))return
    if(.not.allocated(problem%beta_network) .or. &
         .not.allocated(problem%beta_free))return
    if(size(problem%beta_network)/=nmarker .or. &
         size(problem%beta_free)/=nmarker)return
    if(.not.all(is_finite_scalar(problem%beta_network)) .or. &
         .not.all(is_finite_scalar(problem%beta_free)))return
    valid_block_inputs=.true.
  end function valid_block_inputs

  logical function valid_species(species,nx_local,ny_local,nmarker)
    type(actin_species_state_t),intent(in)::species
    integer,intent(in)::nx_local,ny_local,nmarker
    valid_species=.false.
    if(any(shape(species%field)/=[nx_local,ny_local]))return
    if(any(shape(species%density)/=[nmarker]))return
    if(.not.allocated(species%correction_coefficients))return
    if(any(shape(species%correction_coefficients)/=[nmarker,6]))return
    if(.not.all(is_finite_scalar(species%field)) .or. &
         .not.all(is_finite_scalar(species%density)) .or. &
         .not.all(is_finite_scalar(species%correction_coefficients)))return
    valid_species=.true.
  end function valid_species

  logical function valid_trial(trial,nx_local,ny_local,nmarker)
    type(actin_extension_trial_t),intent(in)::trial
    integer,intent(in)::nx_local,ny_local,nmarker
    valid_trial=.false.
    if(.not.allocated(trial%field))return
    if(any(shape(trial%field)/=[nx_local,ny_local]))return
    if(.not.allocated(trial%correction_coefficients))return
    if(any(shape(trial%correction_coefficients)/=[nmarker,6]))return
    if(.not.allocated(trial%correction_grid))return
    if(any(shape(trial%correction_grid)/=[nx_local,ny_local]))return
    if(.not.all(is_finite_scalar(trial%field)) .or. &
         .not.all(is_finite_scalar(trial%correction_coefficients)) .or. &
         .not.all(is_finite_scalar(trial%correction_grid)))return
    if(.not.valid_traces(trial,nmarker))return
    valid_trial=.true.
  end function valid_trial

  logical function valid_traces(trial,nmarker)
    type(actin_extension_trial_t),intent(in)::trial
    integer,intent(in)::nmarker
    valid_traces=.false.
    if(.not.allocated(trial%traces%value_average) .or. &
         .not.allocated(trial%traces%normal_average))return
    if(size(trial%traces%value_average)/=nmarker .or. &
         size(trial%traces%normal_average)/=nmarker)return
    if(.not.all(is_finite_scalar(trial%traces%value_average)) .or. &
         .not.all(is_finite_scalar(trial%traces%normal_average)))return
    valid_traces=.true.
  end function valid_traces

  subroutine checked_multiply(left,right,value,ok)
    real(dp),intent(in)::left,right
    real(dp),intent(out)::value
    logical,intent(out)::ok
    value=0.0_dp; ok=.false.
    if(.not.is_finite_scalar(left) .or. .not.is_finite_scalar(right))return
    if(left==0.0_dp .or. right==0.0_dp)then; ok=.true.; return; endif
    if(exponent(abs(left))+exponent(abs(right))>=maxexponent(value))return
    value=left*right
    ok=is_finite_scalar(value)
  end subroutine checked_multiply

  subroutine checked_add(left,right,value,ok)
    real(dp),intent(in)::left,right
    real(dp),intent(out)::value
    logical,intent(out)::ok
    value=0.0_dp; ok=.false.
    if(.not.is_finite_scalar(left) .or. .not.is_finite_scalar(right))return
    if(sign(1.0_dp,left)==sign(1.0_dp,right))then
      if(max(exponent(abs(left)),exponent(abs(right)))>= &
           maxexponent(value)-1)return
    endif
    value=left+right
    ok=is_finite_scalar(value)
  end subroutine checked_add

end module fsi_actin_block_harness_mod
