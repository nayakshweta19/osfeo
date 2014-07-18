#ifndef __EMP_CULA_SPARSE_COMMON_H__
#define __EMP_CULA_SPARSE_COMMON_H__

/*
 * Copyright (C) 2009-2012 EM Photonics, Inc.  All rights reserved.
 *
 * NOTICE TO USER:   
 *
 * This source code is subject to EM Photonics ownership rights under U.S. and
 * international Copyright laws.  Users and possessors of this source code may
 * not redistribute this code without the express written consent of EM
 * Photonics, Inc.
 *
 * EM PHOTONICS MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED
 * WARRANTY OF ANY KIND.  EM PHOTONICS DISCLAIMS ALL WARRANTIES WITH REGARD TO
 * THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY,
 * NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL EM
 * PHOTONICS BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL, OR CONSEQUENTIAL
 * DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
 * PROFITS,  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
 * ACTION,  ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
 * SOURCE CODE.  
 *
 * U.S. Government End Users.   This source code is a "commercial item" as that
 * term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of "commercial
 * computer  software"  and "commercial computer software documentation" as
 * such terms are  used in 48 C.F.R. 12.212 (SEPT 1995) and is provided to the
 * U.S. Government only as a commercial end item.  Consistent with 48
 * C.F.R.12.212 and 48 C.F.R. 227.7202-1 through 227.7202-4 (JUNE 1995), all
 * U.S. Government End Users acquire the source code with only those rights set
 * forth herein. 
 *
 * Any use of this source code in individual and commercial software must
 * include, in the user documentation and internal comments to the code, the
 * above Disclaimer and U.S. Government End Users Notice.
 *
 */

/// 
/// \file cula_sparse_common.h
/// 
/// Common API functions used by all versions of the CULA Sparse library. 
/// Includes all of the common initialize and error handling functions.
/// 

#include "cula_sparse_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/// Opaque structure holding the CULA Sparse library context
typedef int* culaSparseHandle;

/// Context Management
culaSparseStatus culaSparseCreate(culaSparseHandle* handle);
culaSparseStatus culaSparseDestroy(culaSparseHandle handle);

/// Version Information
culaSparseStatus culaSparseGetVersion(culaSparseHandle handle, int* version);
culaSparseStatus culaSparseGetCudaMinimumVersion(culaSparseHandle handle, int* version);
culaSparseStatus culaSparseGetCudaRuntimeVersion(culaSparseHandle handle, int* version);

/// Preinitialize Routines
culaSparseStatus culaSparsePreinitializeCuda(culaSparseHandle handle);

/// Configuration Initialization
culaSparseStatus culaSparseConfigInit(culaSparseHandle handle, culaSparseConfig* config);

/// Platform Initialization Routines
culaSparseStatus culaSparseHostOptionsInit(culaSparseHandle handle, culaSparseHostOptions* platformOpts);
culaSparseStatus culaSparseCudaOptionsInit(culaSparseHandle handle, culaSparseCudaOptions* platformOpts);
culaSparseStatus culaSparseCudaDeviceOptionsInit(culaSparseHandle handle, culaSparseCudaDeviceOptions* platformOpts);

/// Data Format Initialization Routines
culaSparseStatus culaSparseCooOptionsInit(culaSparseHandle handle, culaSparseCooOptions* formatOpts);
culaSparseStatus culaSparseCsrOptionsInit(culaSparseHandle handle, culaSparseCsrOptions* formatOpts);
culaSparseStatus culaSparseCscOptionsInit(culaSparseHandle handle, culaSparseCscOptions* formatOpts);
culaSparseStatus culaSparseMatrixFreeOptionsInit(culaSparseHandle handle, culaSparseMatrixFreeOptions* formatOpts);

/// Solver Initialization Routines
culaSparseStatus culaSparseBicgOptionsInit(culaSparseHandle handle, culaSparseBicgOptions* solverOpts);
culaSparseStatus culaSparseBicgstabOptionsInit(culaSparseHandle handle, culaSparseBicgstabOptions* solverOpts);
culaSparseStatus culaSparseBicgstablOptionsInit(culaSparseHandle handle, culaSparseBicgstablOptions* solverOpts);
culaSparseStatus culaSparseCgOptionsInit(culaSparseHandle handle, culaSparseCgOptions* solverOpts);
culaSparseStatus culaSparseGmresOptionsInit(culaSparseHandle handle, culaSparseGmresOptions* solverOpts);
culaSparseStatus culaSparseMinresOptionsInit(culaSparseHandle handle, culaSparseMinresOptions* solverOpts);

/// Preconditioner Initialization Routines
culaSparseStatus culaSparseEmptyOptionsInit(culaSparseHandle handle, culaSparseEmptyOptions* precondOpts);
culaSparseStatus culaSparseJacobiOptionsInit(culaSparseHandle handle, culaSparseJacobiOptions* precondOpts);
culaSparseStatus culaSparseBlockjacobiOptionsInit(culaSparseHandle handle, culaSparseBlockjacobiOptions* precondOpts);
culaSparseStatus culaSparseIlu0OptionsInit(culaSparseHandle handle, culaSparseIlu0Options* precondOpts);
culaSparseStatus culaSparseAinvOptionsInit(culaSparseHandle handle, culaSparseAinvOptions* precondOpts);
culaSparseStatus culaSparseFainvOptionsInit(culaSparseHandle handle, culaSparseFainvOptions* precondOpts);
culaSparseStatus culaSparseUserDefinedOptionsInit(culaSparseHandle handle, culaSparseUserDefinedOptions* precondOpts);

/// Print Helper Routines
culaSparseStatus culaSparseGetConfigString(culaSparseHandle handle, const culaSparseConfig* config, char* buf, int bufsize);
culaSparseStatus culaSparseGetLastStatusString(culaSparseHandle handle, char* buf, int bufsize);
culaSparseStatus culaSparseGetResultString(culaSparseHandle handle, const culaSparseResult* result, char* buf, int bufsize);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // __EMP_CULA_SPARSE_H__
