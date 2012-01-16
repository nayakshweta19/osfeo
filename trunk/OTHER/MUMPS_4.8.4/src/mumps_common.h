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
#ifndef MUMPS_COMMON_H
#define MUMPS_COMMON_H
#include "mumps_compat.h"
#include "mumps_c_types.h"
/**
 * F_SYMBOL is a macro that converts a couple (lower case symbol, upper
 * case symbol) into the symbol defined by the compiler convention.
 * Example: For MUMPS_XXX, first define
 *   #define MUMPS_XXX F_SYMBOL(xxx,XXX) and then use
 *   MUMPS_XXX in the code to get rid of any symbol convention annoyance.
 *
 * NB: We need to provide both upper and lower case versions because to our
 *     knowledge, there is no way to perform the conversion with CPP
 *     directives only.
 */
#if defined(UPPER) || defined(MUMPS_WIN32)
# define F_SYMBOL(lower_case,upper_case) MUMPS_##upper_case
#elif defined(Add_)
# define F_SYMBOL(lower_case,upper_case) mumps_##lower_case##_
#elif defined(Add__)
# define F_SYMBOL(lower_case,upper_case) mumps_##lower_case##__
#else
# define F_SYMBOL(lower_case,upper_case) mumps_##lower_case
#endif
MUMPS_INT*
mumps_get_mapping();
#define MUMPS_AFFECT_MAPPING \
    F_SYMBOL(affect_mapping,AFFECT_MAPPING)
void MUMPS_CALL
MUMPS_AFFECT_MAPPING(MUMPS_INT *f77mapping);
#define MUMPS_NULLIFY_C_MAPPING F_SYMBOL(nullify_c_mapping,NULLIFY_C_MAPPING)
void MUMPS_CALL
MUMPS_NULLIFY_C_MAPPING();
MUMPS_INT*
mumps_get_pivnul_list();
#define MUMPS_AFFECT_PIVNUL_LIST \
    F_SYMBOL(affect_pivnul_list,AFFECT_PIVNUL_LIST)
void MUMPS_CALL
MUMPS_AFFECT_PIVNUL_LIST(MUMPS_INT *f77pivnul_list);
#define MUMPS_NULLIFY_C_PIVNUL_LIST \
    F_SYMBOL(nullify_c_pivnul_list,NULLIFY_C_PIVNUL_LIST)
void MUMPS_CALL
MUMPS_NULLIFY_C_PIVNUL_LIST();
MUMPS_INT*
mumps_get_uns_perm();
#define MUMPS_AFFECT_UNS_PERM \
    F_SYMBOL(affect_uns_perm,AFFECT_UNS_PERM)
void MUMPS_CALL
MUMPS_AFFECT_UNS_PERM(MUMPS_INT *f77sym_perm);
#define MUMPS_NULLIFY_C_UNS_PERM \
    F_SYMBOL(nullify_c_uns_perm,NULLIFY_C_UNS_PERM)
void MUMPS_CALL
MUMPS_NULLIFY_C_UNS_PERM();
MUMPS_INT*
mumps_get_sym_perm();
#define MUMPS_AFFECT_SYM_PERM \
    F_SYMBOL(affect_sym_perm,AFFECT_SYM_PERM)
void MUMPS_CALL
MUMPS_AFFECT_SYMBOL(MUMPS_INT *f77uns_perm);
#define MUMPS_NULLIFY_C_SYM_PERM \
    F_SYMBOL(nullify_c_sym_perm,NULLIFY_C_SYM_PERM)
void MUMPS_CALL
MUMPS_NULLIFY_C_SYM_PERM();
#endif /* MUMPS_COMMON_H */
