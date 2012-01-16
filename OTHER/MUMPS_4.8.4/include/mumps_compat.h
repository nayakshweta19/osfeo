/*
 *
 *  This file is part of MUMPS 4.8.4, built on Mon Dec 15 15:31:38 UTC 2008
 *
 *
 *  This version of MUMPS is provided to you free of charge. It is public
 *  domain, based on public domain software developed during the Esprit IV
 *  European project PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL.
 *  Since this first public domain version in 1999, the developments are
 *  supported by the following institutions: CERFACS, ENSEEIHT-IRIT, and
 *  INRIA.
 *
 *  Main contributors are Patrick Amestoy, Iain Duff, Abdou Guermouche,
 *  Jacko Koster, Jean-Yves L'Excellent, and Stephane Pralet.
 *
 *  Up-to-date copies of the MUMPS package can be obtained
 *  from the Web pages:
 *  http://mumps.enseeiht.fr/  or  http://graal.ens-lyon.fr/MUMPS
 *
 *
 *   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
 *   EXPRESSED OR IMPLIED. ANY USE IS AT YOUR OWN RISK.
 *
 *
 *  User documentation of any code that uses this software can
 *  include this complete notice. You can acknowledge (using
 *  references [1], [2], and [3]) the contribution of this package
 *  in any scientific publication dependent upon the use of the
 *  package. You shall use reasonable endeavours to notify
 *  the authors of the package of this publication.
 *
 *   [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent,
 *   Multifrontal parallel distributed symmetric and unsymmetric solvers,
 *   in Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000).
 *
 *   [2] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
 *   A fully asynchronous multifrontal solver using distributed dynamic
 *   scheduling, SIAM Journal of Matrix Analysis and Applications,
 *   Vol 23, No 1, pp 15-41 (2001).
 *
 *   [3] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
 *   S. Pralet, Hybrid scheduling for the parallel solution of linear
 *   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).
 *
 */

/* Compatibility issues between various Windows versions */
#ifndef MUMPS_COMPAT_H
#define MUMPS_COMPAT_H


#if defined(_WIN32) && ! defined(__MINGW32__)
# define MUMPS_WIN32 1
#endif

#ifndef MUMPS_CALL
# ifdef MUMPS_WIN32
/* Modify/choose between next 2 lines depending
 * on your Windows calling conventions */
/* #  define MUMPS_CALL __stdcall */
#  define MUMPS_CALL
# else
#  define MUMPS_CALL
# endif
#endif

#if ( __STDC_VERSION__ >= 199901L )
# define MUMPS_INLINE inline
#else
# define MUMPS_INLINE
#endif


#endif /* MUMPS_COMPAT_H */
