#ifndef rjmcmc_h
#define rjmcmc_h

#include <rjmcmc/rjmcmc_config.h>

/** \mainpage

\section Introduction

This library provides routines for running Reversible Jump Monte-Carlo
Markov chains for data regression and forward models. See the \subpage background section for a general overview and references for more details.

\section How How to use


In order to find the right function to use, you need to know the following
three things about your problem.

- What is the dimension of the problem? 
- Are the discontinuities in either value or gradient that you wish to
determine/model? 
- Is your problem a regression problem or a forward model problem?

The diagram below and the following sections help explain how to
classify your problem and the Directory below will help you identify
which method to use for your problem and have a link to an example
program using that method.

\image html main_tree.png
\image latex main_tree.pdf

\subsection Dimensionality

The rjmcmc library supports 1 and
2 dimensional regression/forward model problems. If your problem consists
of a regression/forward model of a single set of coordinates to a single 
set of values then it is a 1-dimensional problem. If your problem consists
of a regression/forward model of a single set of coordinate pairs (eg x, y)
to a single set of values then it is a 2-dimensional problem.

\subsection Partitioned Partitioned versus Single

If so you will want to use a partitioned or transdimensional
approach, otherwise you can use the single methods (ie single partition
methods).

\subsection Regression Regression versus Forward Model

A regression problem is one where there is a direct relationship
between your data values and the fit to be generated. An example might
be if your data consists of temperature measurements over time and you
wanted to generate the temperature over time curve, then this is a
regression problem. If instead you had a profile of the temperature
down a bore and wanted to reconstruct the temperature versus time at
the surface of the bore then this is a forward model problem. With
forward model problems you will need to provide your own code for the
forward model and calculate the log of the likelihood. In the case the
borehole example, the rjmcmc library will provide a trial surface
temperature record and the forward model will consist of a integration
of the heat-diffusion equation given this input to create a trial
borehole temperature profile. The likelihood is then calculated from
the sum of the squared errors between the measured and trial borehole
temperature profiles.

\section Directory

\subsection singleregress1d 1D Single Partition Regression

The following functions are available:

- ::single1d_regression uses an automatic prior to determine the 
best weights to be given to the order of the fitting polynomial.
- ::single1d_regression_with_prior uses a user supplied prior to
sample varying order polynomials to fit the data (deprecated/used
for testing).
- ::single1d_direct_regression use direct integration to determine
the best fit to the data (deprecated/used for testing).

The following example codes are available (available under the Examples
tab):

- 1d/single/regression/cubic/cubic.c performs a regression on a known
cubic function with added gaussian noise.

\subsection partregress1d 1D Partitioned Regression

The following functions are available:

- ::part1d_regression uses an automatic prior to determine the 
best weights for the polynomial order(s) within each partition.
- ::part1d_zero_regression use 0th order polynomials within
each partition.
- ::part1d_natural_regression use connected line segments between
each partition boundary to provide a continuous fit (with change
points becoming changes in gradients).

The following example codes are available (available under the Examples
tab):

- 1d/partitioned/regression/multiquad/multiquad.c performs a regression
on a function comprised of piece wise quadratic functions with discontinuties
with added gaussian noise.
- 1d/partitioned/regression/multistep/multistep.c
- 1d/partitioned/regression/sawtooth/sawtooth.c
- 1d/partitioned/regression/zeromultistep/zeromultistep.c

\subsection singlefm1d 1D Single Partition Forward Model

The following functions are available:

- ::single_forwardmodel performs a forward model analysis on an arbitrary
number of parameters.
- ::single_forwardmodel_f is a fortran 2003 interface to the ::single_forwardmodel function. 
- ::single_forwardmodel_hierarchical performs a forward model analysis on an
arbitrary number of parameters with a custom hierarchical peturbation of the 
covariance matrix.

The following example codes are available (available under the Examples
tab):

- 1d/single/fm/simplef/simplef.f90 
- 1d/single/fm/simpleimage/simpleimage.c
- 1d/single/fm/spherefit/spherefit.c


\subsection partfm1d 1D Partitioned Forward Model

The following functions are available:

- ::part1d_forwardmodel
- ::part1d_forwardmodel_f
- ::part1d_forwardmodel_hierarchical
- ::part1d_forwardmodel_natural
- ::part1d_forwardmodel_natural_hierarchical

The following example codes are available (available under the Examples
tab):

- 1d/partitioned/fm/functionfit/functionfit.c 
- 1d/partitioned/fm/functionfitf/functionfitf.f90
- 1d/partitioned/fm/regression/regression.c

\subsection partregress2d 2D Partitioned Regression

The following functions are available:

- ::part2d_regression

\subsection partfm2d 2D Partitioned Forward Model

The following functions are available:

- ::part2d_forwardmodel
- ::part2d_forwardmodel_hierarchical

*/

/** \page background Background

\section backgroundoverview Overview

This library provides routines for regression and forward model
problems using trans-dimensional Markov chain Monte Carlo methods.
Rather than producing a single fit for a dataset, these methods
produce an ensemble of results that give a sampled distribution of the
potential true fit. The original paper describing this approach is
by Green \cite green1995, and this software is primarily based on the work
of Bodin and Sambridge \cite sambridge2006A \cite bodinThesis \cite bodin2012A .

\section backgroundprior Prior

The general approach taken in all the routines within this library
is to use uniform priors on all values. 

For partition locations, we use a symmetric Dirichlet prior which
is described in Steininger \cite steininger2013 and references
therein.

\section backgroundmodel Model

For 1D applications we use zeroth order, natural, and polynomial models
within each partition.

Zeroth order 1D models are described in \cite bodinThesis \cite bodin2012A.

Natural uses jointed line segments between partition boundaries to 
create a C0 continuous curve over the domain of the entire model. It is
described in Hopcroft \cite hopcroft2007A .

For 1D regression problems only, we also use the data within each
partition to inform the selection of a suitable polynomial order (up
to a limit imposed by the user), ie trans-dimensional within each
partition. The mechanism for choosing this order is described in
Sambridge \cite sambridge2006A and the order prior is set to uniform
over a region determined from the data mean and standard deviation
within the partition.

For 2D applications at present only zeroth order partitions are used
as described in \cite bodin2012B and \cite bodin2012C.

\section backgroundproposal Proposal

All proposals used in this library use pertubations sampled from a 
Gaussian random variable. The Standard deviation is generally set
as a user parameter, the only exception to this is in the case of 
1D Regression routines where the standard deviation is obtained 
directly from the data.

\section backgroundhierarchical Hierarchical Parameter Estimation

Hierarchical parameters are parameters used to quantify the noise
in the data. In the simple regression case, we allow the use of 
a single hierarchical scaling factor called lambda that represents
a multiplier of the estimated data error.

In the general forward model case, we allow a general mechanism for
hierarchical parameter estimation requiring the forward model to
calculate the log of the determinant of the data covariance matrix.

For an example of more complex hierarchical parameter estimation, see
Bodin \cite bodin2012A .

 */

#endif /* rjmcmc_h */
