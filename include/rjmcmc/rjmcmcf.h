! -*- fortran -*-

use, intrinsic :: iso_c_binding

implicit none

!
! forwardmodelparameter_t
!
type, bind(C) :: forwardmodelparameter_t
      real (kind = c_double) :: fmin
      real (kind = c_double) :: fmax
      real (kind = c_double) :: fstd_value
      real (kind = c_double) :: fstd_bd
end type

type, bind(C) :: bbox2d_t
      real (kind = c_double) :: xmin
      real (kind = c_double) :: xmax
      real (kind = c_double) :: ymin
      real (kind = c_double) :: ymax
end type

integer(kind = c_int), parameter :: RESULTSETFM_MEAN = 1
integer(kind = c_int), parameter :: RESULTSETFM_MEDIAN = 2
integer(kind = c_int), parameter :: RESULTSETFM_MODE = 4
integer(kind = c_int), parameter :: RESULTSETFM_CREDIBLE = 8

integer(kind = c_int), parameter :: RESULTSET1DFM_MEAN = 1
integer(kind = c_int), parameter :: RESULTSET1DFM_MEDIAN = 2
integer(kind = c_int), parameter :: RESULTSET1DFM_MODE = 4
integer(kind = c_int), parameter :: RESULTSET1DFM_CREDIBLE = 8

integer(kind = c_int), parameter :: RESULTSET2DFM_MEAN = 1
integer(kind = c_int), parameter :: RESULTSET2DFM_MEDIAN = 2
integer(kind = c_int), parameter :: RESULTSET2DFM_MODE = 4
integer(kind = c_int), parameter :: RESULTSET2DFM_CREDIBLE = 8

integer(kind = c_int), parameter :: POSITIONMAP2D_LINEAR = 0
integer(kind = c_int), parameter :: POSITIONMAP2D_DELAUNAY = 1
integer(kind = c_int), parameter :: POSITIONMAP2D_QUADTREE = 2

abstract interface

!
! Uniform random number generator interface
!
real (kind=c_double) function rjmcmc_uniform_rand() bind(C)
      
      use, intrinsic :: iso_c_binding

end function rjmcmc_uniform_rand

!
! Gaussian random number generator interface
!
real (kind=c_double) function rjmcmc_normal_rand() bind(C)
  use, intrinsic :: iso_c_binding
end function rjmcmc_normal_rand

!
! Single Forward Model likelihood function interface
!
real (kind=c_double) function single_fm_likelihood(user_arg, n, values) bind(C)

  use, intrinsic :: iso_c_binding 

  type(c_ptr), intent(in), value :: user_arg
  integer (kind = c_int), intent(in), value :: n
  real (kind = c_double), intent(in), dimension(n) :: values
  
end function single_fm_likelihood

!
! Single Forward Model Hierarchical likelihood function interface
!
real (kind=c_double) function single_fm_likelihood_hierarchical(&
      user_arg, &
      n, &
      values, &
      hierarchical, &
      nhierarchicalvalues, &
      hierarchicalvalues, &
      logdetce) bind(C)

      use, intrinsic :: iso_c_binding 
      
      type(c_ptr), intent(in), value :: user_arg
      integer (kind = c_int), intent(in), value :: n
      real (kind = c_double), intent(in), dimension(n) :: values
      integer (kind = c_int), intent(in), value :: hierarchical
      integer (kind = c_int), intent(in), value :: nhierarchicalvalues
      real (kind = c_double), intent(in), dimension(nhierarchicalvalues) &
      :: hierarchicalvalues

      real (kind = c_double), intent(out) :: logdetce

end function single_fm_likelihood_hierarchical

integer(kind=c_int) function part1d_fm_value_at(state, x, maxsize, values) bind(C)
      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: state
      real(kind = c_double), intent(in), value :: x
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), intent(out), dimension(maxsize) :: values

end function part1d_fm_value_at

!
! Part1d Forward Model likelihood function interface
!
function part1d_fm_likelihood(&
      user_arg, &
      npartitions, &
      partitions, &
      nglobalvalues, &
      globalvalues, &
      state, &
      value_at_pointer, &
      gradient_at_pointer) bind(C)

      use, intrinsic :: iso_c_binding 
      
      real (kind=c_double) part1d_fm_likelihood

      type(c_ptr), intent(in), value :: user_arg
      integer(kind = c_int), intent(in), value :: npartitions
      real(kind = c_double), intent(in), dimension(npartitions) :: partitions
      integer(kind = c_int), intent(in), value :: nglobalvalues
      real(kind = c_double), intent(in), dimension(nglobalvalues) :: globalvalues

      type(c_ptr), intent(in), value :: state

      type(c_funptr), intent(in), value :: value_at_pointer
      type(c_funptr), intent(in), value :: gradient_at_pointer

end function part1d_fm_likelihood

function part1d_fm_hierarchical_likelihood(&
      user_arg, &
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

      real(kind=c_double) :: part1d_fm_hierarchical_likelihood
      
      type(c_ptr), intent(in), value :: user_arg
      integer(kind = c_int), intent(in), value :: npartitions
      real(kind = c_double), intent(in), dimension(npartitions) :: partitions
      integer(kind = c_int), intent(in), value :: nglobalvalues
      real(kind = c_double), intent(in), dimension(nglobalvalues) :: globalvalues
      integer(kind = c_int), intent(in), value :: hierarchical
      integer(kind = c_int), intent(in), value :: nhierarchicalvalues
      real(kind = c_double), intent(in), dimension(nhierarchicalvalues) :: hierarchicalvalues

      type(c_ptr), intent(in), value :: state

      type(c_funptr), intent(in), value :: value_at_pointer
      type(c_funptr), intent(in), value :: gradient_at_pointer

      real(kind = c_double), intent(out) :: logdetce

end function part1d_fm_hierarchical_likelihood

!
! 2D Forward model callbacks
integer(kind=c_int) function part2d_fm_value_at(&
      state, &
      x, &
      y, &
      maxsize, &
      values) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: state
      real(kind = c_double), intent(in), value :: x
      real(kind = c_double), intent(in), value :: y
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), intent(out), dimension(maxsize) :: values

end function part2d_fm_value_at

!
! 2D Forward Model likelihood function interface
!
function part2d_fm_likelihood(user_arg, &
      nglobalvalues, &
      globalvalues, &
      state, &
      value_at_pointer, &
      gradient_at_pointer, &
      bound) bind(C)

      use, intrinsic :: iso_c_binding 
      import bbox2d_t
      
      real (kind=c_double) part2d_fm_likelihood

      type(c_ptr), intent(in), value :: user_arg
      integer(kind = c_int), intent(in), value :: nglobalvalues
      real(kind = c_double), intent(in), dimension(nglobalvalues) :: globalvalues
      type(c_ptr), intent(in), value :: state

      type(c_funptr), intent(in), value :: value_at_pointer
      type(c_funptr), intent(in), value :: gradient_at_pointer

      type(bbox2d_t), intent(in) :: bound

end function part2d_fm_likelihood

function part2d_fm_hierarchical_likelihood(user_arg, &
      nglobalvalues, &
      globalvalues, &
      hierarchical, &
      nhierarchicalvalues, &
      hierarchicalvalues, &
      state, &
      value_at_pointer, &
      gradient_at_pointer, &
      bound, &
      logdetce ) bind(C)

      use, intrinsic :: iso_c_binding 
      import bbox2d_t
      
      real (kind=c_double) part2d_fm_hierarchical_likelihood

      type(c_ptr), intent(in), value :: user_arg
      integer(kind = c_int), intent(in), value :: nglobalvalues
      
      real(kind = c_double), intent(in), dimension(nglobalvalues) :: globalvalues
      integer(kind = c_int), intent(in), value :: hierarchical
      integer(kind = c_int), intent(in), value :: nhierarchicalvalues
      real(kind = c_double), intent(in), dimension(nhierarchicalvalues) :: hierarchicalvalues

      type(c_ptr), intent(in), value :: state

      type(c_funptr), intent(in), value :: value_at_pointer
      type(c_funptr), intent(in), value :: gradient_at_pointer

      type(bbox2d_t), intent(in) :: bound
      real (kind = c_double), intent(out) :: logdetce

end function part2d_fm_hierarchical_likelihood

end interface

interface

!
! Random number generator seed
!
subroutine rjmcmc_seed(seed) bind(C)

      use, intrinsic :: iso_c_binding

      integer (kind = c_int), intent(in), value :: seed

end subroutine rjmcmc_seed

!
! Setting position map 2d method
!
subroutine position_map2d_set_type(t) bind (C)

      use, intrinsic :: iso_c_binding

      integer (kind = c_int), intent(in), value :: t

end subroutine position_map2d_set_type

!
! Uniform random number generator
!
real(kind = c_double) function rjmcmc_uniform() bind(C)

      use, intrinsic :: iso_c_binding

end function rjmcmc_uniform

!
! Gaussian random number generator
!
real(kind = c_double) function rjmcmc_normal() bind(C)

      use, intrinsic :: iso_c_binding

end function rjmcmc_normal

!
! Single partition forward model
!
type(c_ptr) function single_forwardmodel_f(&
      burnin, &
      total, &
      random, &
      normal, &
      nparameters, &
      parameters, &
      likelihood, &
      user_arg, &
      samples, &
      credible_interval, &
      requested_results) bind(C)

      use, intrinsic :: iso_c_binding
      
      import forwardmodelparameter_t
      
      integer(kind = c_int), intent(in), value :: burnin
      integer(kind = c_int), intent(in), value :: total
      
      procedure(rjmcmc_uniform_rand), pointer, intent(in), bind(C) :: random
      procedure(rjmcmc_normal_rand), pointer, intent(in), bind(C) :: normal
      
      integer(kind = c_int), intent(in), value :: nparameters
      type(forwardmodelparameter_t), dimension(nparameters), intent(in) :: parameters
      
      procedure(single_fm_likelihood), pointer, intent(in), bind(C) :: likelihood
      type(c_ptr), intent(in), value :: user_arg
      
      integer(kind = c_int), intent(in), value :: samples
      real(kind = c_double), intent(in), value :: credible_interval
      integer(kind = c_int), intent(in), value :: requested_results
      
end function single_forwardmodel_f

type(c_ptr) function single_forwardmodel_hierarchical_f(&
      burnin, &
      total, &
      random, &
      normal, &
      nparameters, &
      parameters, &
      nhierarchicalparameters, &
      hierarchicalparameters, &
      likelihood, &
      user_arg, &
      samples, &
      credible_interval, &
      requested_results) bind(C)

      use, intrinsic :: iso_c_binding
      
      import forwardmodelparameter_t
      
      integer(kind = c_int), intent(in), value :: burnin
      integer(kind = c_int), intent(in), value :: total
      
      procedure(rjmcmc_uniform_rand), pointer, intent(in), bind(C) :: random
      procedure(rjmcmc_normal_rand), pointer, intent(in), bind(C) :: normal
      
      integer(kind = c_int), intent(in), value :: nparameters
      type(forwardmodelparameter_t), dimension(nparameters), intent(in) :: parameters
      
      integer(kind = c_int), intent(in), value :: nhierarchicalparameters
      type(forwardmodelparameter_t), dimension(nparameters), intent(in) :: hierarchicalparameters

      procedure(single_fm_likelihood_hierarchical), pointer, intent(in), bind(C) :: likelihood
      type(c_ptr), intent(in), value :: user_arg
      
      integer(kind = c_int), intent(in), value :: samples
      real(kind = c_double), intent(in), value :: credible_interval
      integer(kind = c_int), intent(in), value :: requested_results
      
end function single_forwardmodel_hierarchical_f

!
! Single Partition Results Struture
!

subroutine resultsetfm_destroy(results) bind(C)

      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: results
      
end subroutine resultsetfm_destroy

integer(kind = c_int) function resultsetfm_get_propose_f(&
      results, &
      maxsize, &
      propose) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: maxsize
      integer(kind = c_int), dimension(maxsize), intent(out) :: propose

end function resultsetfm_get_propose_f

integer(kind = c_int) function resultsetfm_get_accept_f(&
      results, &
      maxsize, &
      accept) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: maxsize
      integer(kind = c_int), dimension(maxsize), intent(out) :: accept

end function resultsetfm_get_accept_f

integer(kind = c_int) function resultsetfm_get_misfit_f(&
      results, &
      maxsize, &
      misfit) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: misfit

end function resultsetfm_get_misfit_f

integer(kind = c_int) function resultsetfm_get_parameter_history_f(&
      results, &
      parameter_index, &
      maxsize, &
      phistory) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: parameter_index
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: phistory

end function resultsetfm_get_parameter_history_f

integer(kind = c_int) function resultsetfm_get_hierarchical_parameter_history_f(&
      results, &
      parameter_index, &
      maxsize, &
      phistory) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: parameter_index
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: phistory

end function resultsetfm_get_hierarchical_parameter_history_f

real(kind = c_double) function resultsetfm_get_parameter_mean(results, parameter_index) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: parameter_index
      
end function resultsetfm_get_parameter_mean

real(kind = c_double) function resultsetfm_get_parameter_mode(results, parameter_index) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: parameter_index
      
end function resultsetfm_get_parameter_mode

real(kind = c_double) function resultsetfm_get_parameter_median(results, parameter_index) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: parameter_index
      
end function resultsetfm_get_parameter_median

real(kind = c_double) function resultsetfm_get_parameter_credible_min(results, parameter_index) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: parameter_index
      
end function resultsetfm_get_parameter_credible_min

real(kind = c_double) function resultsetfm_get_parameter_credible_max(results, parameter_index) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: parameter_index
      
end function resultsetfm_get_parameter_credible_max

!
! 1D Partitioned Forward model
!
type(c_ptr) function part1d_forwardmodel_f(&
      burnin, &
      total, &
      minpart, &
      maxpart, &
      minx, &
      maxx, &
      xsamples, &
      ysamples, &
      credible_interval, &
      pd, &
      random, &
      normal, &
      nglobalparameters, &
      globalparameters, &
      nlocalparameters, &
      localparameters, &
      likelihood, &
      user_arg, &
      requested_results) bind(C)

      use, intrinsic :: iso_c_binding
      
      import forwardmodelparameter_t
      
      integer(kind = c_int), intent(in), value :: burnin
      integer(kind = c_int), intent(in), value :: total
      integer(kind = c_int), intent(in), value :: minpart
      integer(kind = c_int), intent(in), value :: maxpart
      real(kind = c_double), intent(in), value :: minx
      real(kind = c_double), intent(in), value :: maxx
      integer(kind = c_int), intent(in), value :: xsamples
      integer(kind = c_int), intent(in), value :: ysamples
      real(kind = c_double), intent(in), value :: credible_interval
      real(kind = c_double), intent(in), value :: pd

      procedure(rjmcmc_uniform_rand), pointer, intent(in), bind(C) :: random
      procedure(rjmcmc_normal_rand), pointer, intent(in), bind(C) :: normal
      
      integer(kind = c_int), intent(in), value :: nglobalparameters
      type(forwardmodelparameter_t), dimension(nglobalparameters), intent(in) :: globalparameters
      integer(kind = c_int), intent(in), value :: nlocalparameters
      type(forwardmodelparameter_t), dimension(nlocalparameters), intent(in) :: localparameters
      
      procedure(part1d_fm_likelihood), pointer, intent(in), bind(C) :: likelihood
      type(c_ptr), intent(in), value :: user_arg
      
      integer(kind = c_int), intent(in), value :: requested_results
      
end function part1d_forwardmodel_f

type(c_ptr) function part1d_forwardmodel_hierarchical_f(&
      burnin, &
      total, &
      minpart, &
      maxpart, &
      minx, &
      maxx, &
      xsamples, &
      ysamples, &
      credible_interval, &
      pd, &
      random, &
      normal, &
      nglobalparameters, &
      globalparameters, &
      nlocalparameters, &
      localparameters, &
      nhierarchicalparameters, &
      hierarchicalparameters, &
      likelihood, &
      user_arg, &
      requested_results) bind(C)

      use, intrinsic :: iso_c_binding
      
      import forwardmodelparameter_t
      
      integer(kind = c_int), intent(in), value :: burnin
      integer(kind = c_int), intent(in), value :: total
      integer(kind = c_int), intent(in), value :: minpart
      integer(kind = c_int), intent(in), value :: maxpart
      real(kind = c_double), intent(in), value :: minx
      real(kind = c_double), intent(in), value :: maxx
      integer(kind = c_int), intent(in), value :: xsamples
      integer(kind = c_int), intent(in), value :: ysamples
      real(kind = c_double), intent(in), value :: credible_interval
      real(kind = c_double), intent(in), value :: pd

      procedure(rjmcmc_uniform_rand), pointer, intent(in), bind(C) :: random
      procedure(rjmcmc_normal_rand), pointer, intent(in), bind(C) :: normal
      
      integer(kind = c_int), intent(in), value :: nglobalparameters
      type(forwardmodelparameter_t), dimension(nglobalparameters), intent(in) :: globalparameters
      integer(kind = c_int), intent(in), value :: nlocalparameters
      type(forwardmodelparameter_t), dimension(nlocalparameters), intent(in) :: localparameters
      
      integer(kind = c_int), intent(in), value :: nhierarchicalparameters
      type(forwardmodelparameter_t), dimension(nlocalparameters), intent(in) :: hierarchicalparameters

      procedure(part1d_fm_hierarchical_likelihood), pointer, intent(in), bind(C) :: likelihood
      type(c_ptr), intent(in), value :: user_arg
      
      integer(kind = c_int), intent(in), value :: requested_results
      
end function part1d_forwardmodel_hierarchical_f

type(c_ptr) function part1d_forwardmodel_natural_f(&
      burnin, &
      total, &
      minpart, &
      maxpart, &
      minx, &
      maxx, &
      xsamples, &
      ysamples, &
      credible_interval, &
      pd, &
      random, &
      normal, &
      nglobalparameters, &
      globalparameters, &
      nlocalparameters, &
      localparameters, &
      likelihood, &
      user_arg, &
      requested_results) bind(C)

      use, intrinsic :: iso_c_binding
      
      import forwardmodelparameter_t
      
      integer(kind = c_int), intent(in), value :: burnin
      integer(kind = c_int), intent(in), value :: total
      integer(kind = c_int), intent(in), value :: minpart
      integer(kind = c_int), intent(in), value :: maxpart
      real(kind = c_double), intent(in), value :: minx
      real(kind = c_double), intent(in), value :: maxx
      integer(kind = c_int), intent(in), value :: xsamples
      integer(kind = c_int), intent(in), value :: ysamples
      real(kind = c_double), intent(in), value :: credible_interval
      real(kind = c_double), intent(in), value :: pd

      procedure(rjmcmc_uniform_rand), pointer, intent(in), bind(C) :: random
      procedure(rjmcmc_normal_rand), pointer, intent(in), bind(C) :: normal
      
      integer(kind = c_int), intent(in), value :: nglobalparameters
      type(forwardmodelparameter_t), dimension(nglobalparameters), intent(in) :: globalparameters
      integer(kind = c_int), intent(in), value :: nlocalparameters
      type(forwardmodelparameter_t), dimension(nlocalparameters), intent(in) :: localparameters
      
      procedure(part1d_fm_likelihood), pointer, intent(in), bind(C) :: likelihood
      type(c_ptr), intent(in), value :: user_arg
      
      integer(kind = c_int), intent(in), value :: requested_results
      
end function part1d_forwardmodel_natural_f

type(c_ptr) function part1d_forwardmodel_natural_hierarchical_f(&
      burnin, &
      total, &
      minpart, &
      maxpart, &
      minx, &
      maxx, &
      xsamples, &
      ysamples, &
      credible_interval, &
      pd, &
      random, &
      normal, &
      nglobalparameters, &
      globalparameters, &
      nlocalparameters, &
      localparameters, &
      nhierarchicalparameters, &
      hierarchicalparameters, &
      likelihood, &
      user_arg, &
      requested_results) bind(C)

      use, intrinsic :: iso_c_binding
      
      import forwardmodelparameter_t
      
      integer(kind = c_int), intent(in), value :: burnin
      integer(kind = c_int), intent(in), value :: total
      integer(kind = c_int), intent(in), value :: minpart
      integer(kind = c_int), intent(in), value :: maxpart
      real(kind = c_double), intent(in), value :: minx
      real(kind = c_double), intent(in), value :: maxx
      integer(kind = c_int), intent(in), value :: xsamples
      integer(kind = c_int), intent(in), value :: ysamples
      real(kind = c_double), intent(in), value :: credible_interval
      real(kind = c_double), intent(in), value :: pd

      procedure(rjmcmc_uniform_rand), pointer, intent(in), bind(C) :: random
      procedure(rjmcmc_normal_rand), pointer, intent(in), bind(C) :: normal
      
      integer(kind = c_int), intent(in), value :: nglobalparameters
      type(forwardmodelparameter_t), dimension(nglobalparameters), intent(in) :: globalparameters
      integer(kind = c_int), intent(in), value :: nlocalparameters
      type(forwardmodelparameter_t), dimension(nlocalparameters), intent(in) :: localparameters
      
      integer(kind = c_int), intent(in), value :: nhierarchicalparameters
      type(forwardmodelparameter_t), dimension(nlocalparameters), intent(in) :: hierarchicalparameters

      procedure(part1d_fm_hierarchical_likelihood), pointer, intent(in), bind(C) :: likelihood
      type(c_ptr), intent(in), value :: user_arg
      
      integer(kind = c_int), intent(in), value :: requested_results
      
end function part1d_forwardmodel_natural_hierarchical_f

type(c_ptr) function part1d_forwardmodel_natural_hierarchical_init_f(&
      burnin, &
      total, &
      minpart, &
      maxpart, &
      minx, &
      maxx, &
      xsamples, &
      ysamples, &
      credible_interval, &
      pd, &
      random, &
      normal, &
      nglobalparameters, &
      globalparameters, &
      nlocalparameters, &
      localparameters, &
      nhierarchicalparameters, &
      hierarchicalparameters, &
      likelihood, &
      user_arg, &
      requested_results) bind(C)

      use, intrinsic :: iso_c_binding
      
      import forwardmodelparameter_t
      
      integer(kind = c_int), intent(in), value :: burnin
      integer(kind = c_int), intent(in), value :: total
      integer(kind = c_int), intent(in), value :: minpart
      integer(kind = c_int), intent(in), value :: maxpart
      real(kind = c_double), intent(in), value :: minx
      real(kind = c_double), intent(in), value :: maxx
      integer(kind = c_int), intent(in), value :: xsamples
      integer(kind = c_int), intent(in), value :: ysamples
      real(kind = c_double), intent(in), value :: credible_interval
      real(kind = c_double), intent(in), value :: pd

      procedure(rjmcmc_uniform_rand), pointer, intent(in), bind(C) :: random
      procedure(rjmcmc_normal_rand), pointer, intent(in), bind(C) :: normal
      
      integer(kind = c_int), intent(in), value :: nglobalparameters
      type(forwardmodelparameter_t), dimension(nglobalparameters), intent(in) :: globalparameters
      integer(kind = c_int), intent(in), value :: nlocalparameters
      type(forwardmodelparameter_t), dimension(nlocalparameters), intent(in) :: localparameters
      
      integer(kind = c_int), intent(in), value :: nhierarchicalparameters
      type(forwardmodelparameter_t), dimension(nlocalparameters), intent(in) :: hierarchicalparameters

      procedure(part1d_fm_hierarchical_likelihood), pointer, intent(in), bind(C) :: likelihood
      type(c_ptr), intent(in), value :: user_arg
      
      integer(kind = c_int), intent(in), value :: requested_results
      
end function part1d_forwardmodel_natural_hierarchical_init_f

integer(kind = c_int) function part1d_forwardmodel_natural_hierarchical_step_f(&
      state) bind(C)

      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: state

end function part1d_forwardmodel_natural_hierarchical_step_f

type(c_ptr) function part1d_forwardmodel_natural_hierarchical_finish_f(&
      state) bind(C)

      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: state

end function part1d_forwardmodel_natural_hierarchical_finish_f

subroutine resultset1dfm_destroy(results) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results

end subroutine resultset1dfm_destroy

integer(kind = c_int) function resultset1dfm_get_propose_f(&
      results, &
      maxsize, &
      propose) bind (C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: maxsize
      integer(kind = c_int), dimension(maxsize), intent(out) :: propose

end function resultset1dfm_get_propose_f

integer(kind = c_int) function resultset1dfm_get_accept_f(&
      results, &
      maxsize, &
      accept) bind (C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: maxsize
      integer(kind = c_int), dimension(maxsize), intent(out) :: accept

end function resultset1dfm_get_accept_f

integer(kind = c_int) function resultset1dfm_get_misfit_f(&
      results, &
      maxsize, &
      misfit) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: misfit

end function resultset1dfm_get_misfit_f

integer(kind = c_int) function resultset1dfm_get_global_parameter_f(&
      results, &
      pi, &
      maxsize, &
      global_parameter) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: pi
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: global_parameter

end function resultset1dfm_get_global_parameter_f

integer(kind = c_int) function resultset1dfm_get_local_parameter_mean_f(&
      results, &
      li, &
      maxsize, &
      mean) bind (C)
      
      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: li
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: mean

end function resultset1dfm_get_local_parameter_mean_f

integer(kind = c_int) function resultset1dfm_get_local_parameter_median_f(&
      results, &
      li, &
      maxsize, &
      median) bind (C)
      
      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: li
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: median

end function resultset1dfm_get_local_parameter_median_f


integer(kind = c_int) function resultset1dfm_get_local_parameter_mode_f(&
      results, &
      li, &
      maxsize, &
      mode) bind (C)
      
      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: li
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: mode

end function resultset1dfm_get_local_parameter_mode_f

integer(kind = c_int) function resultset1dfm_get_local_parameter_credible_min_f(&
      results, &
      li, &
      maxsize, &
      credible_min) bind (C)
      
      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: li
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: credible_min

end function resultset1dfm_get_local_parameter_credible_min_f


integer(kind = c_int) function resultset1dfm_get_local_parameter_credible_max_f(&
      results, &
      li, &
      maxsize, &
      credible_max) bind (C)
      
      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: li
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: credible_max

end function resultset1dfm_get_local_parameter_credible_max_f

integer(kind = c_int) function resultset1dfm_get_local_parameter_histogram_f(&
      results, &
      li, &
      xsamples, &
      ysamples, &
      histogram) bind (C)
      
      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: li
      integer(kind = c_int), intent(in), value :: xsamples
      integer(kind = c_int), intent(in), value :: ysamples
      integer(kind = c_int), dimension(xsamples, ysamples), intent(out) :: histogram

end function resultset1dfm_get_local_parameter_histogram_f

integer(kind = c_int) function resultset1dfm_get_xcoord_vector_f(&
      results, &
      maxsize, &
      xcoord) bind (C)
      
      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: xcoord

end function resultset1dfm_get_xcoord_vector_f

integer(kind = c_int) function resultset1dfm_get_ycoord_vector_f(&
      results, &
      li, &
      maxsize, &
      ycoord) bind (C)
      
      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: li
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: ycoord

end function resultset1dfm_get_ycoord_vector_f

integer(kind = c_int) function resultset1dfm_get_partitions_f(&
      results, &
      maxsize, &
      partitions) bind (C)
      
      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: maxsize
      integer(kind = c_int), dimension(maxsize), intent(out) :: partitions

end function resultset1dfm_get_partitions_f

integer(kind = c_int) function resultset1dfm_get_partition_x_histogram_f(&
      results, &
      maxsize, &
      histogram) bind (C)
      
      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: maxsize
      integer(kind = c_int), dimension(maxsize), intent(out) :: histogram

end function resultset1dfm_get_partition_x_histogram_f

integer(kind = c_int) function resultset1dfm_get_hierarchical_f(&
      results, &
      pi, &
      maxsize, &
      hierarchical) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: pi
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: hierarchical

end function resultset1dfm_get_hierarchical_f

!
! 2D Partitioned Forward Model
!

type(c_ptr) function part2d_forwardmodel_f(&
      burnin, &
      total, &
      thin, &
      minpart, &
      maxpart, &
      initpart, &
      minx, &
      maxx, &
      miny, &
      maxy, &
      xsamples, &
      ysamples, &
      zsamples, &
      credible_interval, &
      pdx, &
      pdy, &
      random, &
      normal, &
      nglobalparameters, &
      globalparameters, &
      nlocalparameters, &
      localparameters, &
      likelihood, &
      user_arg, &
      requested_results) bind(C)

      use, intrinsic :: iso_c_binding
      
      import forwardmodelparameter_t
      
      integer(kind = c_int), intent(in), value :: burnin
      integer(kind = c_int), intent(in), value :: total
      integer(kind = c_int), intent(in), value :: thin
      integer(kind = c_int), intent(in), value :: minpart
      integer(kind = c_int), intent(in), value :: maxpart
      integer(kind = c_int), intent(in), value :: initpart
      real(kind = c_double), intent(in), value :: minx
      real(kind = c_double), intent(in), value :: maxx
      real(kind = c_double), intent(in), value :: miny
      real(kind = c_double), intent(in), value :: maxy
      integer(kind = c_int), intent(in), value :: xsamples
      integer(kind = c_int), intent(in), value :: ysamples
      integer(kind = c_int), intent(in), value :: zsamples
      real(kind = c_double), intent(in), value :: credible_interval
      real(kind = c_double), intent(in), value :: pdx
      real(kind = c_double), intent(in), value :: pdy

      procedure(rjmcmc_uniform_rand), pointer, intent(in), bind(C) :: random
      procedure(rjmcmc_normal_rand), pointer, intent(in), bind(C) :: normal
      
      integer(kind = c_int), intent(in), value :: nglobalparameters
      type(forwardmodelparameter_t), dimension(nglobalparameters), intent(in) :: globalparameters
      integer(kind = c_int), intent(in), value :: nlocalparameters
      type(forwardmodelparameter_t), dimension(nlocalparameters), intent(in) :: localparameters
      
      procedure(part2d_fm_likelihood), pointer, intent(in), bind(C) :: likelihood
      type(c_ptr), intent(in), value :: user_arg
      
      integer(kind = c_int), intent(in), value :: requested_results
      
end function part2d_forwardmodel_f

type(c_ptr) function part2d_forwardmodel_hierarchical_f(&
      burnin, &
      total, &
      thin, &
      minpart, &
      maxpart, &
      initpart, &
      minx, &
      maxx, &
      miny, &
      maxy, &
      xsamples, &
      ysamples, &
      zsamples, &
      credible_interval, &
      pdx, &
      pdy, &
      random, &
      normal, &
      nglobalparameters, &
      globalparameters, &
      nlocalparameters, &
      localparameters, &
      nhierarchicalparameters, &
      hierarchicalparameters, &
      likelihood, &
      user_arg, &
      requested_results) bind(C)

      use, intrinsic :: iso_c_binding
      
      import forwardmodelparameter_t
      
      integer(kind = c_int), intent(in), value :: burnin
      integer(kind = c_int), intent(in), value :: total
      integer(kind = c_int), intent(in), value :: thin
      integer(kind = c_int), intent(in), value :: minpart
      integer(kind = c_int), intent(in), value :: maxpart
      integer(kind = c_int), intent(in), value :: initpart
      real(kind = c_double), intent(in), value :: minx
      real(kind = c_double), intent(in), value :: maxx
      real(kind = c_double), intent(in), value :: miny
      real(kind = c_double), intent(in), value :: maxy
      integer(kind = c_int), intent(in), value :: xsamples
      integer(kind = c_int), intent(in), value :: ysamples
      integer(kind = c_int), intent(in), value :: zsamples
      real(kind = c_double), intent(in), value :: credible_interval
      real(kind = c_double), intent(in), value :: pdx
      real(kind = c_double), intent(in), value :: pdy

      procedure(rjmcmc_uniform_rand), pointer, intent(in), bind(C) :: random
      procedure(rjmcmc_normal_rand), pointer, intent(in), bind(C) :: normal
      
      integer(kind = c_int), intent(in), value :: nglobalparameters
      type(forwardmodelparameter_t), dimension(nglobalparameters), intent(in) :: globalparameters
      integer(kind = c_int), intent(in), value :: nlocalparameters
      type(forwardmodelparameter_t), dimension(nlocalparameters), intent(in) :: localparameters
      
      integer(kind = c_int), intent(in), value :: nhierarchicalparameters
      type(forwardmodelparameter_t), dimension(nlocalparameters), intent(in) :: hierarchicalparameters

      procedure(part2d_fm_hierarchical_likelihood), pointer, intent(in), bind(C) :: likelihood
      type(c_ptr), intent(in), value :: user_arg
      
      integer(kind = c_int), intent(in), value :: requested_results
      
end function part2d_forwardmodel_hierarchical_f

!
! 2D Forward model results interface
!
subroutine resultset2dfm_destroy(results) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results

end subroutine resultset2dfm_destroy

integer(kind = c_int) function resultset2dfm_get_propose_f(&
      results, &
      maxsize, &
      propose) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: maxsize
      integer(kind = c_int), dimension(maxsize), intent(out) :: propose

end function resultset2dfm_get_propose_f

integer(kind = c_int) function resultset2dfm_get_accept_f(&
      results, &
      maxsize, &
      accept) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: maxsize
      integer(kind = c_int), dimension(maxsize), intent(out) :: accept

end function resultset2dfm_get_accept_f

integer(kind = c_int) function resultset2dfm_get_misfit_f(&
      results, &
      maxsize, &
      misfit) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: misfit

end function resultset2dfm_get_misfit_f

integer(kind = c_int) function resultset2dfm_get_global_parameter_f(&
      results, &
      maxsize, &
      global_parameter) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: global_parameter

end function resultset2dfm_get_global_parameter_f

integer(kind = c_int) function resultset2dfm_get_hierarchical_parameter_f(&
      results, &
      hi, &
      maxsize, &
      hierarchical_parameter) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: hi
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: hierarchical_parameter

end function resultset2dfm_get_hierarchical_parameter_f

integer(kind = c_int) function resultset2dfm_get_local_parameter_mean_f(&
      results, &
      li, &
      xsamples, &
      ysamples, &
      mean) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: li
      integer(kind = c_int), intent(in), value :: xsamples
      integer(kind = c_int), intent(in), value :: ysamples
      real(kind = c_double), dimension(xsamples, ysamples), intent(out) :: mean

end function resultset2dfm_get_local_parameter_mean_f

integer(kind = c_int) function resultset2dfm_get_local_parameter_variance_f(&
      results, &
      li, &
      xsamples, &
      ysamples, &
      variance) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: li
      integer(kind = c_int), intent(in), value :: xsamples
      integer(kind = c_int), intent(in), value :: ysamples
      real(kind = c_double), dimension(xsamples, ysamples), intent(out) :: variance

end function resultset2dfm_get_local_parameter_variance_f

integer(kind = c_int) function resultset2dfm_get_local_parameter_mode_f(&
      results, &
      li, &
      xsamples, &
      ysamples, &
      mode) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: li
      integer(kind = c_int), intent(in), value :: xsamples
      integer(kind = c_int), intent(in), value :: ysamples
      real(kind = c_double), dimension(xsamples, ysamples), intent(out) :: mode

end function resultset2dfm_get_local_parameter_mode_f

integer(kind = c_int) function resultset2dfm_get_local_parameter_median_f(&
      results, &
      li, &
      xsamples, &
      ysamples, &
      median) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: li
      integer(kind = c_int), intent(in), value :: xsamples
      integer(kind = c_int), intent(in), value :: ysamples
      real(kind = c_double), dimension(xsamples, ysamples), intent(out) :: median

end function resultset2dfm_get_local_parameter_median_f

integer(kind = c_int) function resultset2dfm_get_local_parameter_credible_min_f(&
      results, &
      li, &
      xsamples, &
      ysamples, &
      credible_min) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: li
      integer(kind = c_int), intent(in), value :: xsamples
      integer(kind = c_int), intent(in), value :: ysamples
      real(kind = c_double), dimension(xsamples, ysamples), intent(out) :: credible_min

end function resultset2dfm_get_local_parameter_credible_min_f

integer(kind = c_int) function resultset2dfm_get_local_parameter_credible_max_f(&
      results, &
      li, &
      xsamples, &
      ysamples, &
      credible_max) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: li
      integer(kind = c_int), intent(in), value :: xsamples
      integer(kind = c_int), intent(in), value :: ysamples
      real(kind = c_double), dimension(xsamples, ysamples), intent(out) :: credible_max

end function resultset2dfm_get_local_parameter_credible_max_f

integer(kind = c_int) function resultset2dfm_fill_xcoord_vector_f(&
      results, &
      maxsize, &
      xcoord) bind (C)
      
      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: xcoord

end function resultset2dfm_fill_xcoord_vector_f

integer(kind = c_int) function resultset2dfm_fill_ycoord_vector_f(&
      results, &
      maxsize, &
      ycoord) bind (C)
      
      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: maxsize
      real(kind = c_double), dimension(maxsize), intent(out) :: ycoord

end function resultset2dfm_fill_ycoord_vector_f

integer(kind = c_int) function resultset2dfm_get_partitions_f(&
      results, &
      maxsize, &
      partitions) bind(C)

      use, intrinsic :: iso_c_binding

      type(c_ptr), intent(in), value :: results
      integer(kind = c_int), intent(in), value :: maxsize
      integer(kind = c_int), dimension(maxsize), intent(out) :: partitions

end function resultset2dfm_get_partitions_f

integer(kind = c_int) function rjmcmc_vector_to_histogram(&
      s, &
      n, &
      v, &
      histsize, &
      vmin, &
      vmax, &
      hist) bind (C)

      use, intrinsic :: iso_c_binding

      integer(kind = c_int), intent(in), value :: s
      integer(kind = c_int), intent(in), value :: n
      real(kind = c_double), dimension(n), intent(in) :: v

      integer(kind = c_int), intent(in), value :: histsize
      real(kind = c_double), intent(in), value :: vmin
      real(kind = c_double), intent(in), value :: vmax
      integer(kind = c_int), intent(out), dimension(histsize) :: hist

end function rjmcmc_vector_to_histogram

end interface
