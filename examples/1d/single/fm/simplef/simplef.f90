program simplef

!
! Include the rjmcmc head file
!
include 'rjmcmc/rjmcmcf.h'

type(c_ptr) :: results

integer(kind = c_int), parameter ::burnin = 1000
integer(kind = c_int), parameter ::total = 10000

procedure(rjmcmc_uniform_rand), pointer :: random
procedure(rjmcmc_normal_rand), pointer :: normal

integer(kind = c_int), parameter :: nparameters = 1
type(forwardmodelparameter_t), dimension (nparameters) :: parameters
type(c_ptr) :: user_arg

procedure(single_fm_likelihood), pointer :: likelihood

real(kind = c_double) :: credible_interval
integer(kind = c_int) :: samples
integer(kind = c_int) :: requested_results

integer(kind = c_int), dimension(nparameters) :: propose, accept
integer(kind = c_int) :: t

real(kind = c_double), dimension(total) :: misfit
real(kind = c_double), dimension(total) :: history

integer :: i
integer :: ioerror

!
! Initialise the random number generator with a seed
!
call rjmcmc_seed(101)

!
! Set the simulation parameters
!
credible_interval = 0.95
samples = 100

!
! We want the mean and credible intervals calculated
! for us.
!
requested_results = RESULTSETFM_MEAN + RESULTSETFM_CREDIBLE

!
! Set the parameters
!
parameters(1)%fmin = 0.0
parameters(1)%fmax = 10.0
parameters(1)%fstd_value = 1.0
parameters(1)%fstd_bd = 0.0

!
! Use the rjmcmc random routines
!
random => rjmcmc_uniform
normal => rjmcmc_normal

!
! Point to our forward model function
!
likelihood => my_forwardmodel

!
! Run the regression
!
results = single_forwardmodel_f(burnin, total, random, normal, nparameters, &
     parameters, likelihood, user_arg, samples, credible_interval, &
     requested_results)

!
! Get the propose and acceptance counts
!
t = resultsetfm_get_propose_f(results, nparameters, propose)
t = resultsetfm_get_accept_f(results, nparameters, accept)

!
! Print them out
!
write (*,*) "Proposed", propose
write (*,*) "Accepted", accept

!
! Retrieve the log(likelihood)/misfit history and save it to a text
! file.
!
t = resultsetfm_get_misfit_f(results, total, misfit)

open(unit=8, file='misfit.txt', status='replace', action='write', &
     iostat=ioerror)

if (ioerror == 0) then

   do i = 1, total
      write(8, *) misfit(i)
   end do
   
   close(8)

end if

!
! Print the mean of the samples
!
write (*,*) "Mean", resultsetfm_get_parameter_mean(results, 0)

!
! Print the credible intervals
!
write (*,*), "Credible Min", resultsetfm_get_parameter_credible_min(results, 0)
write (*,*), "Credible Max", resultsetfm_get_parameter_credible_max(results, 0)

!
! Parameter histogram
!
t = resultsetfm_get_parameter_history_f(results, 0, total, history)

open(unit=8, file='history.txt', status='replace', action='write', &
     iostat=ioerror)

if (ioerror == 0) then

   do i = 1, total
      write (8, *) history(i)
   end do
   
   close(8)
   
end if


!
! Finally deallocate the results
!
call resultsetfm_destroy(results)

contains

!
! The forward model function
!
function my_forwardmodel(user_arg, n, values) bind(C)

  use, intrinsic :: iso_c_binding 

  !
  ! The return type is a real 
  !
  real (kind = c_double) :: my_forwardmodel

  !
  ! The user argument is a pointer that you can use to point to useful
  ! data or state information to use within this function (rather than
  ! relying on global variables. We don't need it in this case since
  ! the forward model is so simple
  !
  type(c_ptr), intent(in), value :: user_arg

  !
  ! The n parameter represents then dimension of the values array
  !
  integer (kind = c_int), intent(in), value :: n

  !
  ! The values array contains the current model that the log(likelihood)
  ! is required for.
  !
  real (kind = c_double), intent(in), dimension(n) :: values

  real (kind = c_double), parameter :: VALUE_TARGET = 5.0
  real (kind = c_double), parameter :: VALUE_SIGMA = 1.0

  !
  !
  my_forwardmodel = (values(1) - VALUE_TARGET)**2 / (2.0 * VALUE_SIGMA**2)
  
end function my_forwardmodel
 
end program simplef
