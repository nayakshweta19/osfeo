#ifndef __EMP_CULA_SPARSE_CONVERT_H__
#define __EMP_CULA_SPARSE_CONVERT_H__

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

#include "cula_sparse_types.h"

#ifdef __cplusplus
extern "C" {
#endif

culaSparseStatus culaSparseCudaDeviceScoo2csr(int n, int nnz, const int* cooColInd, const int* cooRowInd, const float* cooA, int* csrColInd, int* csrRowPtr, float* csrA, int idxBaseIn, int idxBaseOut);
culaSparseStatus culaSparseCudaDeviceDcoo2csr(int n, int nnz, const int* cooColInd, const int* cooRowInd, const double* cooA, int* csrColInd, int* csrRowPtr, double* csrA, int idxBaseIn, int idxBaseOut);
culaSparseStatus culaSparseCudaDeviceCcoo2csr(int n, int nnz, const int* cooColInd, const int* cooRowInd, const culaSparseComplexFloat* cooA, int* csrColInd, int* csrRowPtr, culaSparseComplexFloat*  csrA, int idxBaseIn, int idxBaseOut);
culaSparseStatus culaSparseCudaDeviceZcoo2csr(int n, int nnz, const int* cooColInd, const int* cooRowInd, const culaSparseComplexDouble* cooA, int* csrColInd, int* csrRowPtr, culaSparseComplexDouble* csrA, int idxBaseIn, int idxBaseOut);

culaSparseStatus culaSparseHostScoo2csr(int n, int nnz, const int* cooColInd, const int* cooRowInd, const float* cooA, int* csrColInd, int* csrRowPtr, float* csrA, int idxBaseIn, int idxBaseOut);
culaSparseStatus culaSparseHostDcoo2csr(int n, int nnz, const int* cooColInd, const int* cooRowInd, const double* cooA, int* csrColInd, int* csrRowPtr, double* csrA, int idxBaseIn, int idxBaseOut);
culaSparseStatus culaSparseHostCcoo2csr(int n, int nnz, const int* cooColInd, const int* cooRowInd, const culaSparseComplexFloat* cooA, int* csrColInd, int* csrRowPtr, culaSparseComplexFloat*  csrA, int idxBaseIn, int idxBaseOut);
culaSparseStatus culaSparseHostZcoo2csr(int n, int nnz, const int* cooColInd, const int* cooRowInd, const culaSparseComplexDouble* cooA, int* csrColInd, int* csrRowPtr, culaSparseComplexDouble* csrA, int idxBaseIn, int idxBaseOut);

culaSparseStatus culaSparseCudaDeviceScsc2csr(int n, int nnz, const int* cscColPtr, const int* cscRowInd, const float* cscA, int* csrColInd, int* csrRowPtr, float* csrA, int idxBaseIn, int idxBaseOut);
culaSparseStatus culaSparseCudaDeviceDcsc2csr(int n, int nnz, const int* cscColPtr, const int* cscRowInd, const double* cscA, int* csrColInd, int* csrRowPtr, double* csrA, int idxBaseIn, int idxBaseOut);
culaSparseStatus culaSparseCudaDeviceCcsc2csr(int n, int nnz, const int* cscColPtr, const int* cscRowInd, const culaSparseComplexFloat* cscA, int* csrColInd, int* csrRowPtr, culaSparseComplexFloat*  csrA, int idxBaseIn, int idxBaseOut);
culaSparseStatus culaSparseCudaDeviceZcsc2csr(int n, int nnz, const int* cscColPtr, const int* cscRowInd, const culaSparseComplexDouble* cscA, int* csrColInd, int* csrRowPtr, culaSparseComplexDouble* csrA, int idxBaseIn, int idxBaseOut);

culaSparseStatus culaSparseHostScsc2csr(int n, int nnz, const int* cscColPtr, const int* cscRowInd, const float* cscA, int* csrColInd, int* csrRowPtr, float* csrA, int idxBaseIn, int idxBaseOut);
culaSparseStatus culaSparseHostDcsc2csr(int n, int nnz, const int* cscColPtr, const int* cscRowInd, const double* cscA, int* csrColInd, int* csrRowPtr, double* csrA, int idxBaseIn, int idxBaseOut);
culaSparseStatus culaSparseHostCcsc2csr(int n, int nnz, const int* cscColPtr, const int* cscRowInd, const culaSparseComplexFloat* cscA, int* csrColInd, int* csrRowPtr, culaSparseComplexFloat*  csrA, int idxBaseIn, int idxBaseOut);
culaSparseStatus culaSparseHostZcsc2csr(int n, int nnz, const int* cscColPtr, const int* cscRowInd, const culaSparseComplexDouble* cscA, int* csrColInd, int* csrRowPtr, culaSparseComplexDouble* csrA, int idxBaseIn, int idxBaseOut);

#ifdef __cplusplus
}
#endif

#endif // __EMP_CULA_SPARSE_CONVERT_H__
