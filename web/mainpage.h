/**
 * \mainpage LIBEFP - The Effective Fragment Potential method implementation
 *
 * \section intro Introduction
 *
 * <i>LIBEFP</i> is a full implementation of the Effective Fragment Potential (EFP)
 * method (see <a href="http://dx.doi.org/10.1021/cr200093j"><i>Gordon et al.,
 * Chem. Rev., 2012</i></a>).
 * <i>LIBEFP</i> is designed to allow developers an easy way to add EFP support to
 * their favourite quantum chemistry software package.
 *
 * \section efpmd EFPMD package
 *
 * <a href="http://github.com/libefp/libefp/tree/master/efpmd"><i>EFPMD</i></a>
 * is a molecular simulation package based on <i>LIBEFP</i>. It allows you to run
 * EFP-only molecular simulations such as geometry optimization and
 * molecular dynamics. See
 * <a href="http://github.com/libefp/libefp/tree/master/efpmd/README.md"><i>README</i></a>
 * file for more information on how to use the package.
 *
 * \section uses Who uses LIBEFP
 *
 * <i>LIBEFP</i> is currently being used in
 * <a href="http://www.psicode.org/"><i>PSI</i></a>
 * (development version) and
 * <a href="http://github.com/libefp/libefp/tree/master/efpmd"><i>EFPMD</i></a>
 * packages.
 *
 * \section api Documentation
 *
 * Documentation for public API is available <a href="efp_8h.html"><i>here</i></a>.
 * For additional information see
 * <a href="http://github.com/libefp/libefp/blob/master/README.md"><i>README</i></a>
 * file.
 *
 * \section code Source code
 *
 * Latest development version of code can be found in git
 * <a href="http://github.com/libefp/libefp"><i>repository</i></a>.
 *
 * \section parallel Parallel performance
 *
 * The EFP code is parallelized using
 * <a href="http://www.openmp.org/"><i>OpenMP</i></a> to utilize the full power
 * of multi-core CPUs. Below is the plot of achieved parallel speedup in a
 * single point energy computation test with a system containing 3375 EFP water
 * molecules. The test was performed using
 * <a href="http://github.com/libefp/libefp/tree/master/efpmd"><i>EFPMD</i></a>
 * program on a machine with two Intel Xeon X5667 CPUs.
 *
 * \image html speedup-8-core.png
 *
 * \copyright Copyright (c) 2012 Ilya Kaliman.
 * Distributed under the terms of BSD 2-clause license. See
 * <a href="http://github.com/libefp/libefp/blob/master/LICENSE"><i>LICENSE</i></a>
 * file for details.
 */
