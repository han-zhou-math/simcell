! Stage 10 matrix-free affine GMRES utility.
!
! Provenance: adapted from the independently verified implementation at
! dualchem commit 29224de.  This module knows nothing about FSI, chemistry,
! actin species, or the C++ multigrid solver.  It solves A*x=-F(0) using
! callback actions A*p=F(p)-F(0), which is the boundary-density equation.
module generic_gmres_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use parameters, only: dp
  implicit none
  private

  integer, parameter, public :: GMRES_STATUS_SUCCESS = 0
  integer, parameter, public :: GMRES_STATUS_INVALID_ARGUMENT = 1
  integer, parameter, public :: GMRES_STATUS_CALLBACK_FAILED = 2
  integer, parameter, public :: GMRES_STATUS_BREAKDOWN = 3
  integer, parameter, public :: GMRES_STATUS_MAX_ITERATIONS = 4
  integer, parameter, public :: GMRES_STATUS_NONFINITE_CALLBACK = 5
  integer, parameter, public :: GMRES_STATUS_NUMERICAL_FAILURE = 6

  type, public :: gmres_options_t
    integer :: restart_length = 30
    integer :: max_iterations = 300
    real(dp) :: relative_tolerance = 1.0e-10_dp
    real(dp) :: absolute_tolerance = 0.0_dp
  end type gmres_options_t

  type, public :: gmres_result_t
    logical :: converged = .false.
    integer :: status = GMRES_STATUS_INVALID_ARGUMENT
    integer :: callback_status = 0
    integer :: iterations = 0
    integer :: restarts = 0
    real(dp) :: absolute_residual = 0.0_dp
    real(dp) :: relative_residual = 0.0_dp
    real(dp) :: initial_residual_norm = 0.0_dp
  end type gmres_result_t

  abstract interface
    subroutine affine_residual_callback(x, residual, context, status)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: residual(:)
      class(*), intent(inout) :: context
      integer, intent(out) :: status
    end subroutine affine_residual_callback
  end interface

  public :: affine_residual_callback
  public :: solve_affine_gmres

contains

  subroutine solve_affine_gmres(evaluate_residual, context, solution, options, &
       result)
    procedure(affine_residual_callback) :: evaluate_residual
    class(*), intent(inout) :: context
    real(dp), intent(out) :: solution(:)
    type(gmres_options_t), intent(in) :: options
    type(gmres_result_t), intent(out) :: result

    real(dp), allocatable :: affine_offset(:), current_residual(:)
    real(dp), allocatable :: callback_value(:), work(:), candidate(:)
    real(dp), allocatable :: basis(:,:), hessenberg(:,:), cosine(:), sine(:)
    real(dp), allocatable :: least_squares_rhs(:)
    real(dp) :: initial_norm, current_norm, column_norm, next_norm
    real(dp) :: residual_estimate, breakdown_scale
    integer :: dimension, restart_length, cycle, column, row, callback_status
    logical :: callback_ok, norm_ok, breakdown, update_ok

    solution = 0.0_dp
    call initialize_result(result)
    dimension = size(solution)
    if (.not. valid_options(options, dimension)) return

    restart_length = min(options%restart_length, dimension)
    allocate(affine_offset(dimension), current_residual(dimension))
    allocate(callback_value(dimension), work(dimension), candidate(dimension))
    allocate(basis(dimension, restart_length + 1))
    allocate(hessenberg(restart_length + 1, restart_length))
    allocate(cosine(restart_length), sine(restart_length))
    allocate(least_squares_rhs(restart_length + 1))

    call evaluate_residual(solution, affine_offset, context, callback_status)
    if (callback_status /= 0) then
      call set_callback_failure(result, callback_status)
      return
    endif
    if (.not. all(is_finite_gmres(affine_offset))) then
      result%status = GMRES_STATUS_NONFINITE_CALLBACK
      return
    endif

    current_residual = affine_offset
    initial_norm = stable_norm(current_residual, norm_ok)
    if (.not. norm_ok) then
      result%status = GMRES_STATUS_NUMERICAL_FAILURE
      return
    endif
    result%initial_residual_norm = initial_norm
    call set_residual_metrics(result, initial_norm, initial_norm)
    if (meets_tolerance(initial_norm, initial_norm, options)) then
      result%converged = .true.
      result%status = GMRES_STATUS_SUCCESS
      return
    endif

    cycle = 0
    do while (result%iterations < options%max_iterations)
      cycle = cycle + 1
      result%restarts = cycle - 1
      basis = 0.0_dp
      hessenberg = 0.0_dp
      cosine = 0.0_dp
      sine = 0.0_dp
      least_squares_rhs = 0.0_dp

      current_norm = stable_norm(current_residual, norm_ok)
      if (.not. norm_ok .or. current_norm <= 0.0_dp) then
        result%status = GMRES_STATUS_NUMERICAL_FAILURE
        return
      endif
      basis(:,1) = -current_residual/current_norm
      least_squares_rhs(1) = current_norm

      do column = 1, restart_length
        call evaluate_residual(basis(:,column), callback_value, context, &
             callback_status)
        if (callback_status /= 0) then
          call set_callback_failure(result, callback_status)
          return
        endif
        if (.not. all(is_finite_gmres(callback_value))) then
          result%status = GMRES_STATUS_NONFINITE_CALLBACK
          return
        endif
        if (.not. subtract_vectors(callback_value, affine_offset, work)) then
          result%status = GMRES_STATUS_NUMERICAL_FAILURE
          return
        endif
        column_norm = stable_norm(work, norm_ok)
        if (.not. norm_ok) then
          result%status = GMRES_STATUS_NUMERICAL_FAILURE
          return
        endif

        call orthogonalize_column(work, basis, column, hessenberg(:,column), &
             update_ok)
        if (.not. update_ok) then
          result%status = GMRES_STATUS_NUMERICAL_FAILURE
          return
        endif
        next_norm = stable_norm(work, norm_ok)
        if (.not. norm_ok) then
          result%status = GMRES_STATUS_NUMERICAL_FAILURE
          return
        endif
        hessenberg(column + 1,column) = next_norm
        breakdown_scale = 100.0_dp*epsilon(1.0_dp)*column_norm
        breakdown = next_norm <= 0.0_dp
        if (column_norm > 0.0_dp) then
          breakdown = breakdown .or. next_norm <= breakdown_scale
        endif
        if (.not. breakdown) basis(:,column + 1) = work/next_norm

        do row = 1, column - 1
          call apply_plane_rotation(hessenberg(row,column), &
               hessenberg(row + 1,column), cosine(row), sine(row))
        enddo
        call make_plane_rotation(hessenberg(column,column), &
             hessenberg(column + 1,column), cosine(column), sine(column), &
             update_ok)
        if (.not. update_ok) then
          result%status = GMRES_STATUS_NUMERICAL_FAILURE
          return
        endif
        call apply_plane_rotation(hessenberg(column,column), &
             hessenberg(column + 1,column), cosine(column), sine(column))
        call apply_plane_rotation(least_squares_rhs(column), &
             least_squares_rhs(column + 1), cosine(column), sine(column))

        result%iterations = result%iterations + 1
        residual_estimate = abs(least_squares_rhs(column + 1))
        if (meets_tolerance(residual_estimate, initial_norm, options) .or. &
            breakdown .or. column == restart_length .or. &
            result%iterations == options%max_iterations) then
          call update_candidate(solution, basis, hessenberg, &
               least_squares_rhs, column, candidate, update_ok)
          if (.not. update_ok) then
            call recompute_returned_residual(evaluate_residual, context, &
                 solution, initial_norm, result, callback_ok)
            if (callback_ok) result%status = GMRES_STATUS_BREAKDOWN
            return
          endif

          call recompute_returned_residual(evaluate_residual, context, &
               candidate, initial_norm, result, callback_ok, current_residual)
          if (.not. callback_ok) return
          solution = candidate
          current_norm = result%absolute_residual
          if (meets_tolerance(current_norm, initial_norm, options)) then
            result%converged = .true.
            result%status = GMRES_STATUS_SUCCESS
            return
          endif
          if (breakdown) then
            result%status = GMRES_STATUS_BREAKDOWN
            return
          endif
          if (result%iterations == options%max_iterations) then
            result%status = GMRES_STATUS_MAX_ITERATIONS
            return
          endif
          exit
        endif
      enddo
    enddo

    call recompute_returned_residual(evaluate_residual, context, solution, &
         initial_norm, result, callback_ok)
    if (callback_ok) result%status = GMRES_STATUS_MAX_ITERATIONS
  end subroutine solve_affine_gmres

  subroutine initialize_result(result)
    type(gmres_result_t), intent(out) :: result

    result%converged = .false.
    result%status = GMRES_STATUS_INVALID_ARGUMENT
    result%callback_status = 0
    result%iterations = 0
    result%restarts = 0
    result%initial_residual_norm = 0.0_dp
    result%absolute_residual = huge(1.0_dp)
    result%relative_residual = huge(1.0_dp)
  end subroutine initialize_result

  pure logical function valid_options(options, dimension) result(valid)
    type(gmres_options_t), intent(in) :: options
    integer, intent(in) :: dimension

    valid = .false.
    if (dimension <= 0) return
    if (options%restart_length <= 0) return
    if (options%max_iterations <= 0) return
    if (.not. is_finite_gmres(options%relative_tolerance)) return
    if (.not. is_finite_gmres(options%absolute_tolerance)) return
    if (options%relative_tolerance < 0.0_dp) return
    if (options%absolute_tolerance < 0.0_dp) return
    valid = .true.
  end function valid_options

  pure elemental logical function is_finite_gmres(value) result(is_finite)
    real(dp), intent(in) :: value
    integer(int64) :: bits

    if (storage_size(value) == 64 .and. radix(value) == 2 .and. &
        digits(value) == 53 .and. minexponent(value) == -1021 .and. &
        maxexponent(value) == 1024) then
      bits = transfer(value, bits)
      is_finite = ibits(bits, 52, 11) /= int(z'7ff', int64)
    else
      is_finite = ieee_is_finite(value)
    endif
  end function is_finite_gmres

  function stable_norm(vector, valid) result(norm_value)
    real(dp), intent(in) :: vector(:)
    logical, intent(out) :: valid
    real(dp) :: norm_value, scale, sum_squares, magnitude, ratio, root_sum
    integer :: index

    norm_value = 0.0_dp
    valid = .false.
    scale = 0.0_dp
    sum_squares = 1.0_dp
    do index = 1, size(vector)
      if (.not. is_finite_gmres(vector(index))) return
      magnitude = abs(vector(index))
      if (magnitude <= 0.0_dp) cycle
      if (scale < magnitude) then
        ratio = scale/magnitude
        sum_squares = 1.0_dp + sum_squares*ratio*ratio
        scale = magnitude
      else
        ratio = magnitude/scale
        sum_squares = sum_squares + ratio*ratio
      endif
    enddo

    if (scale <= 0.0_dp) then
      norm_value = 0.0_dp
    else
      root_sum = sqrt(sum_squares)
      if (exponent(abs(scale)) + exponent(abs(root_sum)) > &
          maxexponent(1.0_dp)) then
        norm_value = huge(1.0_dp)
      else
        norm_value = scale*root_sum
      endif
    endif
    valid = .true.
  end function stable_norm

  logical function subtract_vectors(left, right, difference) result(valid)
    real(dp), intent(in) :: left(:), right(:)
    real(dp), intent(out) :: difference(:)
    integer :: index

    valid = .false.
    difference = 0.0_dp
    if (size(left) /= size(right) .or. size(difference) /= size(left)) return
    do index = 1, size(left)
      if (.not. is_finite_gmres(left(index)) .or. &
          .not. is_finite_gmres(right(index))) return
      if (right(index) < 0.0_dp) then
        if (left(index) > huge(1.0_dp) + right(index)) return
      elseif (right(index) > 0.0_dp) then
        if (left(index) < -huge(1.0_dp) + right(index)) return
      endif
      difference(index) = left(index) - right(index)
    enddo
    valid = all(is_finite_gmres(difference))
  end function subtract_vectors

  subroutine orthogonalize_column(work, basis, column, hessenberg_column, valid)
    real(dp), intent(inout) :: work(:)
    real(dp), intent(in) :: basis(:,:)
    integer, intent(in) :: column
    real(dp), intent(inout) :: hessenberg_column(:)
    logical, intent(out) :: valid
    real(dp) :: projection
    integer :: pass, row

    valid = .false.
    do pass = 1, 2
      do row = 1, column
        projection = dot_product(basis(:,row), work)
        if (.not. is_finite_gmres(projection)) return
        hessenberg_column(row) = hessenberg_column(row) + projection
        if (.not. is_finite_gmres(hessenberg_column(row))) return
        work = work - projection*basis(:,row)
        if (.not. all(is_finite_gmres(work))) return
      enddo
    enddo
    valid = .true.
  end subroutine orthogonalize_column

  pure subroutine make_plane_rotation(first, second, cosine, sine, valid)
    real(dp), intent(in) :: first, second
    real(dp), intent(out) :: cosine, sine
    logical, intent(out) :: valid
    real(dp) :: scale, scaled_magnitude

    cosine = 1.0_dp
    sine = 0.0_dp
    valid = .false.
    if (.not. is_finite_gmres(first) .or. &
        .not. is_finite_gmres(second)) return
    scale = max(abs(first), abs(second))
    if (scale <= 0.0_dp) then
      valid = .true.
      return
    endif
    scaled_magnitude = sqrt((first/scale)**2 + (second/scale)**2)
    if (.not. is_finite_gmres(scaled_magnitude) .or. &
        scaled_magnitude <= 0.0_dp) return
    cosine = (first/scale)/scaled_magnitude
    sine = (second/scale)/scaled_magnitude
    valid = is_finite_gmres(cosine) .and. is_finite_gmres(sine)
  end subroutine make_plane_rotation

  pure subroutine apply_plane_rotation(first, second, cosine, sine)
    real(dp), intent(inout) :: first, second
    real(dp), intent(in) :: cosine, sine
    real(dp) :: temporary

    temporary = cosine*first + sine*second
    second = -sine*first + cosine*second
    first = temporary
  end subroutine apply_plane_rotation

  subroutine update_candidate(solution, basis, hessenberg, rhs, columns, &
       candidate, valid)
    real(dp), intent(in) :: solution(:), basis(:,:), hessenberg(:,:), rhs(:)
    integer, intent(in) :: columns
    real(dp), intent(out) :: candidate(:)
    logical, intent(out) :: valid
    real(dp) :: coefficients(columns), diagonal_scale
    integer :: row

    valid = .false.
    candidate = solution
    coefficients = rhs(1:columns)
    do row = columns, 1, -1
      if (row < columns) then
        coefficients(row) = coefficients(row) - &
             dot_product(hessenberg(row,row + 1:columns), &
             coefficients(row + 1:columns))
      endif
      diagonal_scale = maxval(abs(hessenberg(row,row:columns)))
      if (diagonal_scale <= 0.0_dp) return
      if (abs(hessenberg(row,row)) <= &
          100.0_dp*epsilon(1.0_dp)*diagonal_scale) return
      coefficients(row) = coefficients(row)/hessenberg(row,row)
      if (.not. is_finite_gmres(coefficients(row))) return
    enddo

    candidate = solution + matmul(basis(:,1:columns), coefficients)
    if (.not. all(is_finite_gmres(candidate))) return
    valid = .true.
  end subroutine update_candidate

  subroutine recompute_returned_residual(evaluate_residual, context, solution, &
       initial_norm, result, success, residual_copy)
    procedure(affine_residual_callback) :: evaluate_residual
    class(*), intent(inout) :: context
    real(dp), intent(in) :: solution(:), initial_norm
    type(gmres_result_t), intent(inout) :: result
    logical, intent(out) :: success
    real(dp), intent(out), optional :: residual_copy(:)
    real(dp) :: residual(size(solution)), residual_norm
    integer :: callback_status
    logical :: norm_ok

    success = .false.
    call evaluate_residual(solution, residual, context, callback_status)
    if (callback_status /= 0) then
      call set_callback_failure(result, callback_status)
      return
    endif
    if (.not. all(is_finite_gmres(residual))) then
      result%status = GMRES_STATUS_NONFINITE_CALLBACK
      return
    endif
    residual_norm = stable_norm(residual, norm_ok)
    if (.not. norm_ok) then
      result%status = GMRES_STATUS_NUMERICAL_FAILURE
      return
    endif
    if (present(residual_copy)) residual_copy = residual
    call set_residual_metrics(result, residual_norm, initial_norm)
    success = .true.
  end subroutine recompute_returned_residual

  pure subroutine set_residual_metrics(result, residual_norm, initial_norm)
    type(gmres_result_t), intent(inout) :: result
    real(dp), intent(in) :: residual_norm, initial_norm

    result%absolute_residual = residual_norm
    if (initial_norm <= 0.0_dp) then
      if (residual_norm <= 0.0_dp) then
        result%relative_residual = 0.0_dp
      else
        result%relative_residual = huge(1.0_dp)
      endif
    else
      if (initial_norm < 1.0_dp) then
        if (residual_norm > huge(1.0_dp)*initial_norm) then
          result%relative_residual = huge(1.0_dp)
          return
        endif
      endif
      result%relative_residual = residual_norm/initial_norm
    endif
  end subroutine set_residual_metrics

  pure logical function meets_tolerance(residual_norm, initial_norm, options) &
       result(converged)
    real(dp), intent(in) :: residual_norm, initial_norm
    type(gmres_options_t), intent(in) :: options
    real(dp) :: relative_threshold

    converged = residual_norm <= options%absolute_tolerance
    if (converged .or. options%relative_tolerance <= 0.0_dp) return
    if (options%relative_tolerance > 1.0_dp) then
      if (initial_norm > huge(1.0_dp)/options%relative_tolerance) then
        converged = .true.
        return
      endif
    endif
    relative_threshold = options%relative_tolerance*initial_norm
    converged = residual_norm <= relative_threshold
  end function meets_tolerance

  pure subroutine set_callback_failure(result, callback_status)
    type(gmres_result_t), intent(inout) :: result
    integer, intent(in) :: callback_status

    result%converged = .false.
    result%status = GMRES_STATUS_CALLBACK_FAILED
    result%callback_status = callback_status
  end subroutine set_callback_failure

end module generic_gmres_mod
