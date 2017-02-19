program part1d_fortran

include 'rjmcmc/rjmcmcf.h'

integer, parameter :: MAX_DATA_SIZE = 100
integer, parameter :: MAX_PARTITIONS = 10

type :: data_t
   integer :: size
   real (kind = c_double), dimension(MAX_DATA_SIZE) :: x, y, n

   real (kind = c_double) :: xmin
   real (kind = c_double) :: xmax

end type

type(data_t), target :: data

integer :: i
integer :: status

type(c_ptr) :: user_arg

integer, parameter :: burnin = 10000
integer, parameter :: total = 100000
integer, parameter :: minpartitions = 2
integer, parameter :: maxpartitions = 10
integer, parameter :: xsamples = 100
integer, parameter :: ysamples = 100
integer, parameter :: nprocesses = 6

real (kind = c_double) :: pd
real (kind = c_double) :: credible_interval
integer :: requested_results

integer, parameter :: nparameters = 1
type(forwardmodelparameter_t), dimension(nparameters) :: parameters

integer, parameter :: nglobalparameters = 0
type(forwardmodelparameter_t), dimension (nglobalparameters) :: globalparameters

procedure(part1d_fm_likelihood), pointer :: likelihood
procedure(rjmcmc_uniform_rand), pointer :: random
procedure(rjmcmc_normal_rand), pointer :: normal

type(c_ptr) :: results

integer(kind = c_int), dimension(nprocesses) :: accept
integer(kind = c_int), dimension(nprocesses) :: propose
real(kind = c_double), dimension(total) :: misfit
real(kind = c_double), dimension(xsamples) :: sampled_x;
real(kind = c_double), dimension(xsamples) :: mean;

integer(kind = c_int) :: t
integer :: ioerror

!
! Initialise the random number generator
!
call rjmcmc_seed(101)

!
! Read the input file
!

data%size = 0
data%xmin = 1e37
data%xmax = -1e37
open(8, file = 'data.txt', status = 'old')

do i = 1, MAX_DATA_SIZE
   read (8, *, iostat = status) data%x(i), data%y(i), data%n(i)

   if (status .ne. 0) then
      if (status .eq. -1) then
         exit
      else
         write (*,*) "Unexpected error reading file"
         stop
      end if
   end if
   
   data%xmin = dmin1(data%xmin, data%x(i))
   data%xmax = dmax1(data%xmax, data%x(i))

   data%size = data%size + 1
end do
close(8)

write (*,*) data%xmin, data%xmax

!
! Set up all the simulation parameters
!
credible_interval = 0.95
pd = 1.0

!
! The local parameter is the VS value, ranging between 2.0 ... 5.0
!
parameters(1)%fmin = data%xmin
parameters(1)%fmax = data%xmax
parameters(1)%fstd_value = 0.1
parameters(1)%fstd_bd = 0.0

!
! Use the built in random number generators
!
random => rjmcmc_uniform
normal => rjmcmc_normal

!
! Point the user argument to our data
!
user_arg = c_loc(data)

!
! Point the likelihood function to our forward model
!
likelihood => my_forwardmodel

!
! Setup the requested results
!
requested_results = RESULTSET1DFM_MEAN

!
! Run the MCMC process
!
results = part1d_forwardmodel_f(burnin, &
     total, &
     minpartitions, &
     maxpartitions, &
     data%xmin, &
     data%xmax, &
     xsamples, &
     ysamples, &
     credible_interval, &
     pd, &
     random, &
     normal, &
     nglobalparameters, &
     globalparameters, &
     nparameters, &
     parameters, &
     likelihood, &
     user_arg, &
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

call resultset1dfm_destroy(results)

contains

!
! The PART1D_FORTRAN forward model function
!
function my_forwardmodel(user_arg, &
     npartitions, &
     partitions, &
     nglobalvalues, &
     globalvalues, &
     state, &
     value_at_pointer, &
     gradient_at_pointer) bind(C)

  use, intrinsic :: iso_c_binding

  ! Return value
  real (kind = c_double) :: my_forwardmodel

  ! Parameters
  type (c_ptr), intent(in), value :: user_arg
  integer (kind = c_int), intent(in), value :: npartitions
  real (kind = c_double), intent(in), dimension(npartitions) :: partitions
  integer (kind = c_int), intent(in), value :: nglobalvalues
  real (kind = c_double), intent(in), dimension(nglobalvalues) :: globalvalues
  type (c_ptr), intent(in), value :: state
  type(c_funptr), intent(in), value :: value_at_pointer
  type(c_funptr), intent(in), value :: gradient_at_pointer

  ! Local variables
  type(data_t), pointer :: data_pointer
  procedure(part1d_fm_value_at), pointer :: value_at
  procedure(part1d_fm_value_at), pointer :: gradient_at

  real (kind = c_double), dimension(nparameters) :: values

  integer :: i
  real (kind = c_double) :: phi, loglike

  ! Reference globalvalues parameter to remove warning
  i = size(globalvalues)

  call c_f_pointer(user_arg, data_pointer)
  call c_f_procpointer(value_at_pointer, value_at)
  call c_f_procpointer(gradient_at_pointer, gradient_at)

  loglike = 0.0
  do i = 1, data_pointer%size

     if (value_at(state, data_pointer%x(i), nparameters, values) .lt. 0) then
        write (*,*) "Failed to get value.", data_pointer%x(i)
        stop
     end if
     
     phi = values(1) - data_pointer%y(i)
     
     loglike = loglike + phi**2/(2.0 * data_pointer%n(i) ** 2)

  end do 

  my_forwardmodel = loglike

end function my_forwardmodel

end program part1d_fortran
