! -*- fortran -*-

interface 

!
! Single partition forward model
!
type(c_ptr) function MPI_single_forwardmodel_f(&
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
      requested_results, &
      mpisize, &
      mpirank, &
      root, &
      comm) bind(C)

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

      integer(kind = c_int), intent(in), value :: mpisize
      integer(kind = c_int), intent(in), value :: mpirank
      integer(kind = c_int), intent(in), value :: root
      integer(kind = c_int), intent(in), value :: comm
      
end function MPI_single_forwardmodel_f

type(c_ptr) function MPI_single_forwardmodel_hierarchical_f(&
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
      requested_results, &
      mpisize, &
      mpirank, &
      root, &
      comm) bind(C)

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
      
      integer(kind = c_int), intent(in), value :: mpisize
      integer(kind = c_int), intent(in), value :: mpirank
      integer(kind = c_int), intent(in), value :: root
      integer(kind = c_int), intent(in), value :: comm

end function MPI_single_forwardmodel_hierarchical_f

type(c_ptr) function MPI_part1d_forwardmodel_f(&
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
      requested_results, &
      mpisize, &
      mpirank, &
      root, &
      comm) bind(C)

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
      
      integer(kind = c_int), intent(in), value :: mpisize
      integer(kind = c_int), intent(in), value :: mpirank
      integer(kind = c_int), intent(in), value :: root
      integer(kind = c_int), intent(in), value :: comm

end function MPI_part1d_forwardmodel_f

type(c_ptr) function MPI_part1d_forwardmodel_hierarchical_f(&
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
      requested_results, &
      mpisize, &
      mpirank, &
      root, &
      comm) bind(C)

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
      
      integer(kind = c_int), intent(in), value :: mpisize
      integer(kind = c_int), intent(in), value :: mpirank
      integer(kind = c_int), intent(in), value :: root
      integer(kind = c_int), intent(in), value :: comm

end function MPI_part1d_forwardmodel_hierarchical_f

type(c_ptr) function MPI_part1d_forwardmodel_natural_f(&
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
      requested_results, &
      mpisize, &
      mpirank, &
      root, &
      comm) bind(C)

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
      
      integer(kind = c_int), intent(in), value :: mpisize
      integer(kind = c_int), intent(in), value :: mpirank
      integer(kind = c_int), intent(in), value :: root
      integer(kind = c_int), intent(in), value :: comm

end function MPI_part1d_forwardmodel_natural_f

type(c_ptr) function MPI_part1d_forwardmodel_natural_hierarchical_f(&
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
      requested_results, &
      mpisize, &
      mpirank, &
      root, &
      comm) bind(C)

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
      
      integer(kind = c_int), intent(in), value :: mpisize
      integer(kind = c_int), intent(in), value :: mpirank
      integer(kind = c_int), intent(in), value :: root
      integer(kind = c_int), intent(in), value :: comm

end function MPI_part1d_forwardmodel_natural_hierarchical_f

type(c_ptr) function MPI_part2d_forwardmodel_f(&
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
      requested_results, &
      mpisize, &
      mpirank, &
      root, &
      comm) bind(C)

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
      
      integer(kind = c_int), intent(in), value :: mpisize
      integer(kind = c_int), intent(in), value :: mpirank
      integer(kind = c_int), intent(in), value :: root
      integer(kind = c_int), intent(in), value :: comm

end function MPI_part2d_forwardmodel_f

type(c_ptr) function MPI_part2d_forwardmodel_hierarchical_f(&
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
      requested_results, &
      mpisize, &
      mpirank, &
      root, &
      comm) bind(C)

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
      
      integer(kind = c_int), intent(in), value :: mpisize
      integer(kind = c_int), intent(in), value :: mpirank
      integer(kind = c_int), intent(in), value :: root
      integer(kind = c_int), intent(in), value :: comm

end function MPI_part2d_forwardmodel_hierarchical_f

type(c_ptr) function MPI_part2d_forwardmodel_hierarchical_restartable_f(&
      infile_template, &
      outfile_template, &
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
      requested_results, &
      mpisize, &
      mpirank, &
      root, &
      comm) bind(C)

      use, intrinsic :: iso_c_binding
      
      import forwardmodelparameter_t
      
      character(kind=c_char) :: infile_template(*)
      character(kind=c_char) :: outfile_template(*)
      
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
      
      integer(kind = c_int), intent(in), value :: mpisize
      integer(kind = c_int), intent(in), value :: mpirank
      integer(kind = c_int), intent(in), value :: root
      integer(kind = c_int), intent(in), value :: comm

end function MPI_part2d_forwardmodel_hierarchical_restartable_f

end interface
