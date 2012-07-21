/**
 * \mainpage LIBEFP - The Effective Fragment Potential method implementation
 *
 * \section API Public API Documentation
 *
 * The library is still in active development and API is UNSTABLE.
 * Documentation for public API is available <a href="efp_8h.html">here</a>.
 * For additional information see
 * <a href="http://github.com/libefp/libefp/blob/master/README.md">README</a>
 * file.
 *
 * \section Repo Git Repository
 *
 * Latest development version of code can be found in git
 * <a href="http://github.com/libefp/libefp">repository</a>.
 *
 * \section EFPMD EFPMD package
 *
 * <a href="http://github.com/libefp/libefp/tree/master/efpmd">EFPMD</a>
 * is a molecular simulation package based on
 * <a href="http://libefp.github.com/">LIBEFP</a>. It allows to
 * perform EFP-only molecular simulation such as geometry optimization and
 * molecular dynamics. See
 * <a href="http://github.com/libefp/libefp/tree/master/efpmd/README.md">README</a>
 * file for more information on how to use the package.
 *
 * \section Parallel Parallel Performance
 *
 * The EFP code is parallelized using
 * <a href="http://www.openmp.org/">OpenMP</a> to utilize the full power of
 * multi-core CPUs. Below is the plot of achieved parallel speedup in a single
 * point energy computation test with a system containing 3375 EFP water
 * molecules. The test was performed using
 * <a href="http://github.com/libefp/libefp/tree/master/efpmd">EFPMD</a>
 * program on a machine with two Intel Xeon X5667 CPUs.
 *
 * \image html speedup-8-core.png
 *
 * \copyright Copyright (c) 2012 Ilya Kaliman.
 * Distributed under the terms of BSD 2-clause license. See
 * <a href="http://github.com/libefp/libefp/blob/master/LICENSE">LICENSE</a>
 * file for details.
 */
