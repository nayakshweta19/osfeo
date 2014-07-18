#ifndef __EMP_CULA_SPARSE_H__
#define __EMP_CULA_SPARSE_H__

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

// 
// cula_sparse.h
// 

#include "cula_sparse_common.h"

#ifdef __cplusplus
extern "C" {
#endif

// Opaque structure holding the CULA Sparse execution plan
typedef int* culaSparsePlan;

//
// Plan Management
//

culaSparseStatus culaSparseCreatePlan(culaSparseHandle handle, culaSparsePlan* plan);
culaSparseStatus culaSparseDestroyPlan(culaSparsePlan plan);
culaSparseStatus culaSparseExecutePlan(culaSparseHandle handle, culaSparsePlan plan, const culaSparseConfig* config, culaSparseResult* result);

//
// Data & Platform Setting Methods
//

// CSR
culaSparseStatus culaSparseSetScsrData(culaSparseHandle handle, culaSparsePlan plan, const culaSparseCsrOptions* formatOpts, int n, int nnz, const float* a, const int* rowPtr, const int* colInd, float* x, const float* b);
culaSparseStatus culaSparseSetDcsrData(culaSparseHandle handle, culaSparsePlan plan, const culaSparseCsrOptions* formatOpts, int n, int nnz, const double* a, const int* rowPtr, const int* colInd, double* x, const double* b);
culaSparseStatus culaSparseSetCcsrData(culaSparseHandle handle, culaSparsePlan plan, const culaSparseCsrOptions* formatOpts, int n, int nnz, const culaSparseComplexFloat* a, const int* rowPtr, const int* colInd, culaSparseComplexFloat* x, const culaSparseComplexFloat* b);
culaSparseStatus culaSparseSetZcsrData(culaSparseHandle handle, culaSparsePlan plan, const culaSparseCsrOptions* formatOpts, int n, int nnz, const culaSparseComplexDouble* a, const int* rowPtr, const int* colInd, culaSparseComplexDouble* x, const culaSparseComplexDouble* b);

// CSC
culaSparseStatus culaSparseSetScscData(culaSparseHandle handle, culaSparsePlan plan, const culaSparseCscOptions* formatOpts, int n, int nnz, const float* a, const int* rowInd, const int* colPtr, float* x, const float* b);
culaSparseStatus culaSparseSetDcscData(culaSparseHandle handle, culaSparsePlan plan, const culaSparseCscOptions* formatOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b);
culaSparseStatus culaSparseSetCcscData(culaSparseHandle handle, culaSparsePlan plan, const culaSparseCscOptions* formatOpts, int n, int nnz, const culaSparseComplexFloat* a, const int* rowInd, const int* colPtr, culaSparseComplexFloat* x, const culaSparseComplexFloat* b);
culaSparseStatus culaSparseSetZcscData(culaSparseHandle handle, culaSparsePlan plan, const culaSparseCscOptions* formatOpts, int n, int nnz, const culaSparseComplexDouble* a, const int* rowInd, const int* colPtr, culaSparseComplexDouble* x, const culaSparseComplexDouble* b);

// COO
culaSparseStatus culaSparseSetScooData(culaSparseHandle handle, culaSparsePlan plan, const culaSparseCooOptions* formatOpts, int n, int nnz, const float* a, const int* rowInd, const int* colInd, float* x, const float* b);
culaSparseStatus culaSparseSetDcooData(culaSparseHandle handle, culaSparsePlan plan, const culaSparseCooOptions* formatOpts, int n, int nnz, const double* a, const int* rowInd, const int* colInd, double* x, const double* b);
culaSparseStatus culaSparseSetCcooData(culaSparseHandle handle, culaSparsePlan plan, const culaSparseCooOptions* formatOpts, int n, int nnz, const culaSparseComplexFloat* a, const int* rowInd, const int* colInd, culaSparseComplexFloat* x, const culaSparseComplexFloat* b);
culaSparseStatus culaSparseSetZcooData(culaSparseHandle handle, culaSparsePlan plan, const culaSparseCooOptions* formatOpts, int n, int nnz, const culaSparseComplexDouble* a, const int* rowInd, const int* colInd, culaSparseComplexDouble* x, const culaSparseComplexDouble* b);

// Matrix Free
culaSparseStatus culaSparseSetSmatrixfreeData(culaSparseHandle handle, culaSparsePlan plan, const culaSparseMatrixFreeOptions* formatOpts, int n, void (*ax)(char trans, const float* input, float* ouput, void* user), void* user, float* x, const float* b);
culaSparseStatus culaSparseSetDmatrixfreeData(culaSparseHandle handle, culaSparsePlan plan, const culaSparseMatrixFreeOptions* formatOpts, int n, void (*ax)(char trans, const double* input, double* ouput, void* user), void* user, double* x, const double* b);
culaSparseStatus culaSparseSetCmatrixfreeData(culaSparseHandle handle, culaSparsePlan plan, const culaSparseMatrixFreeOptions* formatOpts, int n, void (*ax)(char trans, const culaSparseComplexFloat* input, culaSparseComplexFloat* ouput, void* user), void* user, culaSparseComplexFloat* x, const culaSparseComplexFloat* b);
culaSparseStatus culaSparseSetZmatrixfreeData(culaSparseHandle handle, culaSparsePlan plan, const culaSparseMatrixFreeOptions* formatOpts, int n, void (*ax)(char trans, const culaSparseComplexDouble* input, culaSparseComplexDouble* ouput, void* user), void* user, culaSparseComplexDouble* x, const culaSparseComplexDouble* b);


//
// Platform Setting Routines
//

culaSparseStatus culaSparseSetHostPlatform(culaSparseHandle handle, culaSparsePlan plan, const culaSparseHostOptions* platformOpts);
culaSparseStatus culaSparseSetCudaPlatform(culaSparseHandle handle, culaSparsePlan plan, const culaSparseCudaOptions* platformOpts);
culaSparseStatus culaSparseSetCudaDevicePlatform(culaSparseHandle handle, culaSparsePlan plan, const culaSparseCudaDeviceOptions* platformOpts);

//
// Preconditioner Setting Methods
//

culaSparseStatus culaSparseSetNoPreconditioner(culaSparseHandle handle, culaSparsePlan plan, const culaSparseEmptyOptions* precondOpts);
culaSparseStatus culaSparseSetIlu0Preconditioner(culaSparseHandle handle, culaSparsePlan plan, const culaSparseIlu0Options* precondOpts);
culaSparseStatus culaSparseSetJacobiPreconditioner(culaSparseHandle handle, culaSparsePlan plan, const culaSparseJacobiOptions* precondOpts);
culaSparseStatus culaSparseSetBlockJacobiPreconditioner(culaSparseHandle handle, culaSparsePlan plan, const culaSparseBlockjacobiOptions* precondOpts);
culaSparseStatus culaSparseSetAinvPreconditioner(culaSparseHandle handle, culaSparsePlan plan, const culaSparseAinvOptions* precondOpts);
culaSparseStatus culaSparseSetFainvPreconditioner(culaSparseHandle handle, culaSparsePlan plan, const culaSparseFainvOptions* precondOpts);
culaSparseStatus culaSparseSetUserDefinedPreconditioner(culaSparseHandle handle, culaSparsePlan plan, const culaSparseUserDefinedOptions* precondOpts, void* user);

//
// Solver Setting Methods
//

culaSparseStatus culaSparseSetNoSolver(culaSparseHandle handle, culaSparsePlan plan, const culaSparseEmptyOptions* solverOpts);
culaSparseStatus culaSparseSetBicgSolver(culaSparseHandle handle, culaSparsePlan plan, const culaSparseBicgOptions* solverOpts);
culaSparseStatus culaSparseSetBicgstablSolver(culaSparseHandle handle, culaSparsePlan plan, const culaSparseBicgstablOptions* solverOpts);
culaSparseStatus culaSparseSetBicgstabSolver(culaSparseHandle handle, culaSparsePlan plan, const culaSparseBicgstabOptions* solverOpts);
culaSparseStatus culaSparseSetCgSolver(culaSparseHandle handle, culaSparsePlan plan, const culaSparseCgOptions* solverOpts);
culaSparseStatus culaSparseSetGmresSolver(culaSparseHandle handle, culaSparsePlan plan, const culaSparseGmresOptions* solverOpts);
culaSparseStatus culaSparseSetMinresSolver(culaSparseHandle handle, culaSparsePlan plan, const culaSparseMinresOptions* solverOpts);

//
// Plan Utilities
//

culaSparseStatus culaSparseGetPlanString(culaSparseHandle handle, culaSparsePlan plan, char* buf, int bufsize);
culaSparseStatus culaSparseGetPlanPreconditionerString(culaSparseHandle handle, culaSparsePlan plan, char* buf, int bufsize);
culaSparseStatus culaSparseGetPlanSolverString(culaSparseHandle handle, culaSparsePlan plan, char* buf, int bufsize);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // __EMP_CULA_SPARSE_H__
