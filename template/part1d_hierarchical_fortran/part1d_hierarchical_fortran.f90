program part1d_hierarchical_fortran

include 'rjmcmc/rjmcmcf.h'

integer, parameter :: MAX_DATA_SIZE = 3000
integer, parameter :: MAX_PARTITIONS = 10

!
! We use this data structure for all out data and work space
! for computing the forward model results. The data arrays
! are fixed size (adjust MAX_DATA_SIZE if necessary) and
! the work arrays are allocated once the number of data points
! is known.
!
type :: data_t
   
   ! Actual size of the data
   integer :: size 

   ! Data
   real (kind = c_double), dimension(MAX_DATA_SIZE) :: x, y

   ! Data x-range
   real (kind = c_double) :: xmin
   real (kind = c_double) :: xmax

   real (kind = c_double) :: ymin
   real (kind = c_double) :: ymax

   ! Data Covariance Matrix
   real (kind = c_double), dimension(:,:), allocatable :: Cd

   ! SVD Decomposion of Data Covariance Matrix
   real (kind = c_double), dimension(:,:), allocatable :: U, V
   real (kind = c_double), dimension(:), allocatable :: S

   ! Work Arrays
   real (kind = c_double), dimension(:), allocatable :: B, A
   real (kind = c_double), dimension(:), allocatable :: WORK
   integer :: WORKSIZE

   integer, dimension(:), allocatable :: IWORK
   integer :: IWORKSIZE

   integer :: info

end type

type(data_t), target :: data

integer :: i
integer :: status

type(c_ptr) :: user_arg

!
! Parameters for the simulation
!
integer, parameter :: burnin = 10000
integer, parameter :: total = 100000

integer, parameter :: minpartitions = 2
integer, parameter :: maxpartitions = 10

integer, parameter :: xsamples = 100
integer, parameter :: ysamples = 100

integer, parameter :: nprocesses = 5

real (kind = c_double) :: pd
real (kind = c_double) :: credible_interval
integer :: requested_results

!
! Since this is a simple regression problem we have one parameter 
! (a single value at each partition boundary).
!
integer, parameter :: nparameters = 1
type(forwardmodelparameter_t), dimension(nparameters) :: parameters

integer, parameter :: nglobalparameters = 0
type(forwardmodelparameter_t), dimension (nglobalparameters) :: globalparameters

!
! For this example we are using 3 hierarchical parameters to construct
! a data covariance matrix, there parameters are sigma, p, q which
! are just the diagonal and off-diagonal values, ie the matrix will look
! like:
!
! | s p q 0 0 0 ... 0|
! | p s p q 0 0 ... 0|
! | q p s p q 0 ... 0|     
! | ................ |
! | 0 ... q p s p q 0|
! | 0 ... 0 q p s p q|
! | 0 ... 0 0 q p s p|
! | 0 ... 0 0 0 q p s|
!
integer, parameter :: nhierarchicalparameters = 3
type(forwardmodelparameter_t), dimension (nhierarchicalparameters) :: hierarchicalparameters

procedure(part1d_fm_hierarchical_likelihood), pointer :: likelihood
procedure(rjmcmc_uniform_rand), pointer :: random
procedure(rjmcmc_normal_rand), pointer :: normal

type(c_ptr) :: results

integer(kind = c_int), dimension(nprocesses) :: accept
integer(kind = c_int), dimension(nprocesses) :: propose
real(kind = c_double), dimension(total) :: misfit, hierarchical_sigma
real(kind = c_double), dimension(xsamples) :: sampled_x;
real(kind = c_double), dimension(xsamples) :: mean;

integer(kind = c_int) :: t

!
! Initialise the random number generator
!
call rjmcmc_seed(1012)

!
! Read the input file
!

data%size = 0
data%xmin = 1e37
data%xmax = -1e37
data%ymin = 1e37
data%ymax = -1e37
open(8, file = 'data.txt', status = 'old')

do i = 1, MAX_DATA_SIZE
   read (8, *, iostat = status) data%x(i), data%y(i)

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

   data%ymin = dmin1(data%ymin, data%y(i))
   data%ymax = dmax1(data%ymax, data%y(i))

   data%size = data%size + 1
end do
close(8)

!
! Allocate the matrices and work variables
!
allocate(data%Cd(data%size, data%size))
allocate(data%U(data%size, data%size))
allocate(data%V(data%size, data%size))

allocate(data%S(data%size))
allocate(data%B(data%size))
allocate(data%A(data%size))

data%IWORKSIZE = 8 * data%size
allocate(data%IWORK(data%IWORKSIZE))

!
! The SVD decomposition routine requires some workspace. To compute the
! size of this we set -1 as the size parameter and it returns the size
! required to the B vector here.
!
call dgesdd("A", data%size, data%size, data%Cd, data%size, data%S, data%U, &
     data%size, data%V, data%size, data%B, -1, data%IWORK, data%info)

data%WORKSIZE = data%B(1)
allocate(data%WORK(data%WORKSIZE))

!
! Set up all the simulation parameters
!
credible_interval = 0.95

!
! pd is the standard deviation of the normal distribution to sample from
! when moving a partition boundary
!
pd = 1.0

!
! The local parameter is the VS value, ranging between 2.0 ... 5.0
!
parameters(1)%fmin = data%ymin
parameters(1)%fmax = data%ymax
parameters(1)%fstd_value = 0.1
parameters(1)%fstd_bd = 0.1

!
! The hierarchical parameters
!
hierarchicalparameters(1)%fmin = 0.5
hierarchicalparameters(1)%fmax = 2.0
hierarchicalparameters(1)%fstd_value = 0.2
hierarchicalparameters(1)%fstd_bd = 0.0


hierarchicalparameters(2)%fmin = 0.0
hierarchicalparameters(2)%fmax = 0.1
hierarchicalparameters(2)%fstd_value = 0.01
hierarchicalparameters(2)%fstd_bd = 0.0

hierarchicalparameters(3)%fmin = 0.0
hierarchicalparameters(3)%fmax = 0.1
hierarchicalparameters(3)%fstd_value = 0.01
hierarchicalparameters(3)%fstd_bd = 0.0

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
results = part1d_forwardmodel_natural_hierarchical_f(burnin, &
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
     nhierarchicalparameters, &
     hierarchicalparameters, &
     likelihood, &
     user_arg, &
     requested_results)

if (.not. c_associated(results)) then
   write (*,*) "Failed to run regression"
   stop
end if

!
! Get the propose and acceptance counts
!
propose = 0
accept = 0
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
call save_vector("misfit.txt", misfit, total)

!
! Retrieve the mean fit and the sampled x coordinates
!
t = resultset1dfm_get_xcoord_vector_f(results, xsamples, sampled_x)
t = resultset1dfm_get_local_parameter_mean_f(results, 0, xsamples, mean)

call save_xy("mean.txt", sampled_x, mean, xsamples)

!
! Retrieve the hierarchical sigma history and save it to a text
! file (note we reuse the same array for the different hierarchical
! parameters).
!
t = resultset1dfm_get_hierarchical_f(results, 0, total, hierarchical_sigma)
call save_vector("hierarchical_sigma.txt", hierarchical_sigma, total)

t = resultset1dfm_get_hierarchical_f(results, 1, total, hierarchical_sigma)
call save_vector("hierarchical_p.txt", hierarchical_sigma, total)

t = resultset1dfm_get_hierarchical_f(results, 2, total, hierarchical_sigma)
call save_vector("hierarchical_q.txt", hierarchical_sigma, total)

!
! Destroy the results
!
call resultset1dfm_destroy(results)

!
! Deallocate the work arrays
!
deallocate(data%Cd)
deallocate(data%U)
deallocate(data%V)

deallocate(data%S)
deallocate(data%B)
deallocate(data%A)

deallocate(data%IWORK)

deallocate(data%WORK)

contains

!
! The PART1D_FORTRAN forward model function
!
function my_forwardmodel(user_arg, &
     npartitions, &
     partitions, &
     nglobalvalues, &
     globalvalues, &
     hierarchical, &
     nhierarchicalvalues, &
     hierarchicalvalues, &
     state, &
     value_at_pointer, &
     gradient_at_pointer, &
     logdetce) bind(C)

  use, intrinsic :: iso_c_binding

  ! Return value
  real (kind = c_double) :: my_forwardmodel

  ! Parameters
  type (c_ptr), intent(in), value :: user_arg
  integer (kind = c_int), intent(in), value :: npartitions
  real (kind = c_double), intent(in), dimension(npartitions) :: partitions
  integer (kind = c_int), intent(in), value :: nglobalvalues
  real (kind = c_double), intent(in), dimension(nglobalvalues) :: globalvalues
  integer(kind = c_int), intent(in), value :: hierarchical
  integer(kind = c_int), intent(in), value :: nhierarchicalvalues
  real(kind = c_double), intent(in), dimension(nhierarchicalvalues) :: hierarchicalvalues
  type (c_ptr), intent(in), value :: state
  type(c_funptr), intent(in), value :: value_at_pointer
  type(c_funptr), intent(in), value :: gradient_at_pointer

  real(kind = c_double), intent(out) :: logdetce

  ! Local variables
  type(data_t), pointer :: data_pointer
  procedure(part1d_fm_value_at), pointer :: value_at
  procedure(part1d_fm_value_at), pointer :: gradient_at

  real (kind = c_double), dimension(nparameters) :: values

  integer :: i
  real (kind = c_double) :: loglike, sigma, p, q
  real (kind = c_double) :: Strace

  ! Reference globalvalues parameter to remove warning
  i = size(globalvalues)

  call c_f_pointer(user_arg, data_pointer)
  call c_f_procpointer(value_at_pointer, value_at)
  call c_f_procpointer(gradient_at_pointer, gradient_at)

  sigma = hierarchicalvalues(1)
  p = hierarchicalvalues(2)
  q = hierarchicalvalues(3)

  !
  ! Create the error covariance matrix
  !
  call create_Cd(sigma, p, q, data_pointer%Cd, data_pointer%size)

  !
  ! Do the SVD decomposition of the matrix
  !
  call dgesdd("A", data_pointer%size, data_pointer%size, data_pointer%Cd, &
       data_pointer%size, data_pointer%S, data_pointer%U, data_pointer%size, &
       data_pointer%V, data_pointer%size, data_pointer%WORK, data_pointer%WORKSIZE, &
       data_pointer%IWORK, data_pointer%info)

  if (data_pointer%info .ne. 0) then
     !
     ! Numerical error computing SVD.
     !
     write (*,*) "Failed to compute SVD for Cd matrix."
     stop
  end if

  Strace = 1.0
  do i = 1, data_pointer%size
     Strace = Strace * data_pointer%S(i)
  end do

  !
  ! Compute the error vector between the model and the data
  !
  do i = 1, data_pointer%size

     if (value_at(state, data_pointer%x(i), nparameters, values) .lt. 0) then
        write (*,*) "Failed to get value.", data_pointer%x(i)
        stop
     end if

     data_pointer%B(i) = values(1) - data_pointer%y(i) 
     
  end do 

  !
  ! Use back substitution to compute Cd^-1 * G
  !
  call svdbacksub(data_pointer%U, data_pointer%S, data_pointer%V, &
       data_pointer%B, data_pointer%WORK, data_pointer%A, data_pointer%size)

  !
  ! Do the dot product to compute G^T * Cd^-1 * G
  !
  loglike = 0.0
  do i = 1, data_pointer%size
     
     loglike = loglike + data_pointer%A(i) * data_pointer%B(i)

  end do

  !
  ! If we have an hiearchical step we need the determinant and since
  ! we have the SVD decomposition, this is just the product of the singular
  ! values.
  !
  if (hierarchical .eq. 1) then
     logdetce = 1.0
     do i = 1, data%size
        logdetce = logdetce * data_pointer%S(i)
     end do

     logdetce = log(logdetce)

  end if

  !
  ! Finally divide by 2 for the log of the likelihood
  !
  if (loglike .lt. 0.0) then
     !
     ! It is possible (particularly in this example) for the Cd matrix
     ! to not be positive definite or due to numerical imprecision for 
     ! the loglikelihood to go below zero. In this case we return a very
     ! large log likelihood to ensure that this set of parameters is
     ! rejected.
     !
     ! Note if this warning is printing out a great deal then this may
     ! indicate a bug or some other error.
     !
     write (*,*) "Warning: Loglike < 0", sigma, p, q
     my_forwardmodel = 1.0d+30
  else
     my_forwardmodel = loglike/2.0
  end if
  

end function my_forwardmodel

subroutine create_Cd(sigma, p, q, A, N) 

  real (kind = c_double), intent(in) :: sigma, p, q
  real (kind = c_double), intent(out), dimension(N,N) :: A
  integer, intent(in) :: N

  integer :: i

  A = 0.0

  do i = 1, N
     
     A(i,i) = sigma*sigma

  end do

  if (N .gt. 1) then
     do i = 1, N - 1
        
        A(i + 1, i) = p
        A(i, i + 1) = p
        
     end do
  end if
     
  if (N .gt. 2) then
     do i = 1, N - 2
        
        A(i + 2, i) = q
        A(i, i + 2) = q
        
     end do
  end if

end subroutine

subroutine save_vector(filename, v, n)

  character(len = *), intent(in) :: filename
  real(kind = c_double), dimension(n), intent(in) :: v
  integer, intent(in) :: n

  integer :: ioerror
  integer :: i

  open(unit=8, file=filename, status='replace', action='write', &
     iostat=ioerror)
  
  if (ioerror == 0) then
     
     do i = 1, n
        write(8, *) v(i)
     end do
     
     close(8)
     
  end if

end subroutine

subroutine save_xy(filename, x, y, n)

  character(len = *), intent(in) :: filename
  real(kind = c_double), dimension(n), intent(in) :: x, y
  integer, intent(in) :: n

  integer :: ioerror
  integer :: i

  open(unit=8, file=filename, status='replace', action='write', &
     iostat=ioerror)
  
  if (ioerror == 0) then
     
     do i = 1, n
        write(8, *) x(i), y(i)
     end do
     
     close(8)
     
  end if

end subroutine

end program part1d_hierarchical_fortran
