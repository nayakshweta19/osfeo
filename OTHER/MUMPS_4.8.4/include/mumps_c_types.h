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


#ifndef MUMPS_C_TYPES_H
#define MUMPS_C_TYPES_H

#define MUMPS_INT int

#define SMUMPS_COMPLEX float
#define SMUMPS_REAL float

#define DMUMPS_COMPLEX double
#define DMUMPS_REAL double

/* Complex datatypes */
typedef struct {float r,i;} mumps_complex;
typedef struct {double r,i;} mumps_double_complex;

#define CMUMPS_COMPLEX mumps_complex
#define CMUMPS_REAL float

#define ZMUMPS_COMPLEX mumps_double_complex
#define ZMUMPS_REAL double


#ifndef mumps_ftnlen
/* When passing a string, what is the type of the extra argument
 * passed by value ? */
# define mumps_ftnlen int
#endif


#define MUMPS_ARITH_s 1
#define MUMPS_ARITH_d 2
#define MUMPS_ARITH_c 4
#define MUMPS_ARITH_z 8

#define MUMPS_ARITH_REAL   ( MUMPS_ARITH_s | MUMPS_ARITH_d )
#define MUMPS_ARITH_CMPLX  ( MUMPS_ARITH_c | MUMPS_ARITH_z )
#define MUMPS_ARITH_SINGLE ( MUMPS_ARITH_s | MUMPS_ARITH_c )
#define MUMPS_ARITH_DBL    ( MUMPS_ARITH_d | MUMPS_ARITH_z )


#endif /* MUMPS_C_TYPES_H */