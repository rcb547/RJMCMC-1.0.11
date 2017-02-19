program functionfitf
!
! Include the rjmcmc head file
!
include 'rjmcmc/rjmcmcf.h'

type(c_ptr) :: results

!
! This data structure contains the synthetic data that is analysed.
! It consists of displacement (s) measurements versus time (t) and 
! the forward model attempts to reconstruct the true velocity (v).
! The true velocity as a function of time is defined in the 
! my_realfunction below.
!
integer, parameter :: datasize = 100
type :: mydata_t
   real(kind = c_double) :: minx
   real(kind = c_double) :: maxx
   real(kind = c_double) :: sigma
   real(kind = c_double), dimension(datasize) :: x
   real(kind = c_double), dimension(datasize) :: ry
   real(kind = c_double), dimension(datasize) :: sy
end type

type(mydata_t), target :: data

!
! There are always 4 processes for 1d forward models (birth, death,
! move, change local value, change global value).
!
integer(kind = c_int), parameter :: nprocesses = 5

integer(kind = c_int), parameter :: burnin = 1000
integer(kind = c_int), parameter :: total = 50000

procedure(rjmcmc_uniform_rand), pointer :: random
procedure(rjmcmc_normal_rand), pointer :: normal

integer(kind = c_int), parameter :: nparameters = 1
type(forwardmodelparameter_t), dimension (nparameters) :: parameters
type(c_ptr) :: user_arg

integer(kind = c_int), parameter :: nglobalparameters = 0
type(forwardmodelparameter_t), dimension (nglobalparameters) :: globalparameters

procedure(part1d_fm_likelihood), pointer :: likelihood

integer(kind = c_int) :: minpartitions
integer(kind = c_int) :: maxpartitions
real(kind = c_double) :: minx
real(kind = c_double) :: maxx
real(kind = c_double) :: credible_interval
integer(kind = c_int), parameter :: xsamples = 100
integer(kind = c_int), parameter :: ysamples = 100
real(kind = c_double) :: pd
integer(kind = c_int) :: requested_results

integer(kind = c_int), dimension(nprocesses) :: propose, accept
integer(kind = c_int) :: j

real(kind = c_double), dimension(total) :: misfit

real(kind = c_double), dimension(xsamples) :: sampled_x;
real(kind = c_double), dimension(xsamples) :: mean;
real(kind = c_double), dimension(xsamples) :: cred_min;
real(kind = c_double), dimension(xsamples) :: cred_max;

integer(kind = c_int), dimension(total) :: partitions
integer(kind = c_int), dimension(xsamples) :: partition_x_hist

integer :: ioerror
integer :: i

real(kind = c_double) :: x
real(kind = c_double) :: dx
integer(kind = c_int) :: t

!
! Initialise the random number generator with a seed
!
call rjmcmc_seed(101)

!
! Create some data from the real forward model and add some noise
! to the displacement measurements
!
data%minx = 0.0
data%maxx = 10.0
data%sigma = 5.0
x = data%minx
dx = (data%maxx - data%minx)/(datasize - 1)
do i = 1, datasize
   data%x(i) = x
   data%ry(i) = my_realfunction(x)
   data%sy(i) = data%ry(i) + rjmcmc_normal() * data%sigma 
   x = x + dx
enddo

!
! Set the simulation parameters
!
minpartitions = 2
maxpartitions = 10
credible_interval = 0.95
pd = 0.2

!
! We want the mean and credible intervals calculated
! for us.
!
requested_results = RESULTSET1DFM_MEAN + RESULTSET1DFM_CREDIBLE

!
! Set the parameters
!
parameters(1)%fmin = -50.0
parameters(1)%fmax = 50.0
parameters(1)%fstd_value = 1.0
parameters(1)%fstd_bd = 1.0

!
! Use the rjmcmc random routines
!
random => rjmcmc_uniform
normal => rjmcmc_normal

!
! Assign the user argument to point to our data
!
user_arg = c_loc(data)

!
! Point to our forward model function
!
likelihood => my_forwardmodel

!
! Run the regression
!
results = part1d_forwardmodel_f(burnin, total, &
     minpartitions, maxpartitions, &
     data%minx, data%maxx, xsamples, ysamples, credible_interval, pd, &
     random, normal, &
     nglobalparameters, globalparameters, &
     nparameters, parameters, &
     likelihood, user_arg, &
     requested_results)

!
! Get the propose and acceptance counts
!
t = resultset1dfm_get_propose_f(results, nprocesses, propose)
t = resultset1dfm_get_accept_f(results, nprocesses, accept)

!
! Print them out
!
write (*,*) "Proposed", propose
write (*,*) "Accepted", accept

!
! Retrieve the log(likelihood)/misfit history and save it to a text
! file.
!
t = resultset1dfm_get_misfit_f(results, total, misfit)

open(unit=8, file='misfit.txt', status='replace', action='write', &
     iostat=ioerror)

if (ioerror == 0) then

do i = 1, total
   write(8, *) misfit(i)
end do

close(8)

end if

!
! Retrieve the mean fit and the sampled x coordinates
!
t = resultset1dfm_get_xcoord_vector_f(results, xsamples, sampled_x)
t = resultset1dfm_get_local_parameter_mean_f(results, 0, xsamples, mean)

open(unit=8, file='mean.txt', status='replace', action='write', &
     iostat=ioerror)

if (ioerror == 0) then

do i = 1, xsamples
   write(8, *) sampled_x(i), mean(i)
end do

close(8)

end if

!
! Retrieve the credible intervals
!
t = resultset1dfm_get_local_parameter_credible_min_f(results, 0, xsamples, cred_min)
t = resultset1dfm_get_local_parameter_credible_max_f(results, 0, xsamples, cred_max)

open(unit=8, file='credible.txt', status='replace', action='write', &
     iostat=ioerror)

if (ioerror == 0) then

do i = 1, xsamples
   write(8, *) sampled_x(i), cred_min(i), cred_max(i)
end do

close(8)

end if

!
! Write the partition count history
!
t = resultset1dfm_get_partitions_f(results, total, partitions)

open(unit=8, file='partitions.txt', status='replace', action='write', &
     iostat=ioerror)

if (ioerror == 0) then

do i = 1, total
   write(8, *) partitions(i)
end do

close(8)

end if

!
! Write the partition location histogram
!
t = resultset1dfm_get_partition_x_histogram_f(results, xsamples, partition_x_hist)

open(unit=8, file='partition_x_hist.txt', status='replace', action='write', &
     iostat=ioerror)

if (ioerror == 0) then

do i = 1, xsamples
   write(8, *) sampled_x(i), partition_x_hist(i)
end do

close(8)

end if


!
! Write the actual curve out
!
open(unit=8, file='actual.txt', status='replace', action='write', &
     iostat=ioerror)

if (ioerror == 0) then

do i = 1, datasize
   write(8, *) data%x(i), data%ry(i), data%sy(i)
end do

close(8)

end if

!
! Finally deallocate the results
!
call resultset1dfm_destroy(results)

contains

!
! The forward model function
!
real (kind = c_double) function my_forwardmodel(user_arg, npartitions, &
     partitions, nglobalvalues, globalvalues, state, value_at_pointer, &
     gradient_at_pointer) bind(C)

  use, intrinsic :: iso_c_binding 

  !
  ! The user argument is a pointer that you can use to point to useful
  ! data or state information to use within this function (rather than
  ! relying on global variables. This points to our data structure 
  ! we created in the main program but before we can use it, we need to
  ! convert it to the write Fortran pointer (see below the call to 
  ! c_f_pointer).
  !
  type(c_ptr), intent(in), value :: user_arg

  !
  ! We are given an ordered list of the partition boundaries
  !
  integer (kind = c_int), intent(in), value :: npartitions
  real (kind = c_double), intent(in), dimension(npartitions) :: partitions

  !
  ! The list of current global values
  !
  integer (kind = c_int), intent(in), value :: nglobalvalues
  real (kind = c_double), intent(in), dimension(nglobalvalues) :: globalvalues
  
  !
  ! The opaque state information needed to pass to the callback functions
  !
  type(c_ptr), intent(in), value :: state
  
  !
  ! The two callback functions, we are passed C pointers to these functions
  ! however we need to convert them to Fortran procedure pointers to use
  ! them (see the calls to c_f_procpointer below).
  !
  type(c_funptr), intent(in), value :: value_at_pointer
  type(c_funptr), intent(in), value :: gradient_at_pointer

  procedure(part1d_fm_value_at), pointer :: value_at
  procedure(part1d_fm_value_at), pointer :: gradient_at

  !
  ! Our local variables we use for computing the forward model
  !
  real(kind = c_double), dimension(nparameters) :: values

  type(mydata_t), pointer :: data_pointer

  integer :: i
  integer :: j
  real(kind = c_double) :: py
  real(kind = c_double) :: sum

  !
  ! We need to convert the pointer values we receive from the C
  ! library to equivalent Fortran procedure pointers in order
  ! to use them. We can do this with the intrinsic c_f_procpointer
  ! subroutine.
  !
  call c_f_procpointer(value_at_pointer, value_at)
  call c_f_procpointer(gradient_at_pointer, gradient_at)

  !
  ! Similarly, when we passed in the user_arg pointer we used
  ! the c_loc function. To convert this value back into a
  ! Fortran pointer to our data, we undo this with the 
  ! intrinsic c_f_pointer subroutine.
  !
  call c_f_pointer(user_arg, data_pointer)

  !
  ! The actual forward model code
  !
  sum = 0.0
  do i = 1, datasize

     !
     ! We can use the value_at procedure pointer to query the proposed
     ! velocity profile versus time at different times. This method
     ! returns the number of local parameters (1 in this example) and
     ! fills the values array with the local parameter values at the
     ! point data_pointer%t(i). After this call is successfull, the
     ! proposed velocity at the given time will be the first element
     ! of the array values.
     !
     j = value_at(state, data_pointer%x(i), nparameters, values)
     
     py = values(1)

     ! 
     ! Accumulate the sum of the square errors between the displacement
     ! profile derived from the proposed velocity profile and the 
     ! measured displacement.
     !
     sum = sum + (py - data_pointer%sy(i))**2

  end do

  !
  ! Return the sum of the square errors.
  !

  !write (*,*) sum
  my_forwardmodel = sum/(2.0 * data_pointer%sigma**2)
  
end function my_forwardmodel

real(kind = c_double) function my_realfunction(t)
  
  real(kind = c_double) :: t

  !
  ! This is the true function we are trying to discern, three
  ! evenly spaced partitions of constant positive velocities
  ! with discontinuities.
  !
  if (t < 3.3) then
     my_realfunction = 20.0
  else if (t < 6.6) then
     my_realfunction = 30.0
  else
     my_realfunction = 10.0
  end if

end function my_realfunction
  
end program functionfitf
