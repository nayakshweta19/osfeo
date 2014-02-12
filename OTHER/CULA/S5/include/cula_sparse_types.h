#ifndef __EMP_CULA_SPARSE_TYPES_H__
#define __EMP_CULA_SPARSE_TYPES_H__

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
/// \file cula_sparse_types.h
/// 
/// Type definitions, structures, options, and macros for the CULA Sparse library
///

#define CULA_SPARSE_VERSION_NUMBER 5000
#define CULA_SPARSE_VERSION_STRING "S5 CUDA 5.0"

#ifdef __cplusplus
extern "C" {
#endif

/// 
/// Error Status Codes
/// 
/// All CULA Sparse routines return an error status indicating if is the routine
/// succeed or not.
/// 

typedef enum {
    culaSparseNoError,                      ///< The command completed successfully
    culaSparseNonConvergence,               ///< The iterative solver did not converge within the selected constraints
    culaSparsePreconditionerError,          ///< The selected preconditioner failed to generate
    culaSparseArgumentError,                ///< An invalid argument was passed to a function
    culaSparseDataFormatError,              ///< The matrix data format is malformed
    culaSparseInsufficientMemory,           ///< There is insufficient memory to continue
    culaSparseFeatureNotImplemented,        ///< The requested feature has not been implemented
    culaSparseRuntimeError,                 ///< A runtime error has occurred
    culaSparseAlignmentError,               ///< Unaligned data was encountered when required 
    culaSparseInteralError,                 ///< An internal error has occurred
    culaSparseHandleError,                  ///< An invalid reference handle was passed to a function
    culaSparsePlanError,                    ///< An invalid execution plan was passed to a function
    culaSparseUnspecifiedError              ///< An unspecified internal error has occurred
} culaSparseStatus;

/// 
/// Iterative Solver Flag Codes
/// 
/// CULA Sparse solver routines return additional information further describing
/// why the solver succeeded or failed.
/// 

typedef enum {
    culaSparseConvergedToRelative,          ///< The solver converged to the specified relative tolerance
    culaSparseConvergedToAbsolute,          ///< The solver converged to the specified absolute tolerance
    culaSparseDiverged,                     ///< The solver diverged by exceeding the specified divergence tolerance
    culaSparseMaxIterationsReached,         ///< Maximum iterations reached without convergence
    culaSparseMaxTimeReached,               ///< Maximum execution time reached without convergence
    culaSparseStagnation,                   ///< The solver stagnated
    culaSparseBreakdown,                    ///< An intermediate scalar value was out of range
    culaSparseInterrupted,                  ///< The routine was manually interrupted
    culaSparseSolveNotExecuted,             ///< The solver was not executed 
    culaSparseUnknownIterationError         ///< An unknown iteration error was encountered
} culaSparseFlag;

/// 
/// General Solver Configuration
///
/// Contains information that steers execution of an iterative solver routine.
/// These options are independent of the selected platform and solver.
/// 

typedef struct {
    /// The relative tolerance for which convergence is achieved.
    double relativeTolerance;

    /// The absolute tolerance for which convergence is achieved.
    double absoluteTolerance;

    /// The relative tolerance of the residual norm for which the iterative 
    /// solver will determine the server is diverging.
    double divergenceTolerance;

    /// The maximum number of iterations that the solver will attempt.
    int maxIterations;

    /// The maximum time, in seconds, at which the solver will not begin a new    
    /// iteration. If set to 0, the solver will run until convergence or the 
    /// maximum number of iterations has been reached.
    double maxRuntime;

    /// This parameter provides the means for a user to capture the relative 
    /// residual at each iteration. The specified array must be at least 
    /// maxIter+1 in length. This parameter may be NULL if these quantities 
    /// are not desired.
    double* residualVector;

    /// Indicates whether the 'x' vector in iterative solves should be used as
    /// given or ignored. When ignored, the 'x' vector is considered a zero.
    int useInitialResultVector;

    /// Indicates whether the 'x' vector in iterative solves should return the
    /// final answer or the best answer and its associated iteration number in 
    /// the case of non-convergence. 
    int useBestAnswer;

    /// Indicates whether to check whether the iterative solve stagnated. This
    /// option defaults to on; turning this option off will increase performance
    /// if a problem is certain not to stagnate.
    int useStagnationCheck;
} culaSparseConfig;

///
/// Residual Information
///
/// Contains information about the residual of an iterative function
///

typedef struct {
    double r0;                    ///< The norm of the initial preconditioned residual vector
    double relative;              ///< The relative residual obtained by the solver 
    double absolute;              ///< The absolute residual obtained by the solver
    double* byIteration;          ///< If requested, the residual at every step of iteration
} culaSparseResidualInfo;

/// 
/// Timing Information
/// 
/// Contains accurate timing information for the all stages of execution of an 
/// iterative function
///

typedef struct {
    double solve;                 ///< Time, in seconds, the solve portion of the iterative solver took to complete
    double preconditioner;        ///< Time, in seconds, the preconditioner generative portion of the iterative solver took to complete
    double overhead;              ///< Time, in seconds, of overhead needed by the iterative solver; includes any memory transfers
    double total;                 ///< Time, in seconds, the entire iterative solver took to complete
} culaSparseTiming;

///
/// Solver Results
///
/// Collects information regarding the execution of an iterative function. This
/// data is only valid if the solver converged or diverges.
///

typedef struct {
    culaSparseConfig config;          ///< Copy of the configuration a given solver was called with
    culaSparseFlag flag;              ///< Enumeration containing information about the success or failure of the iterative solver
    int iterations;                   ///< Number of iterations taken by the iterative solver
    culaSparseResidualInfo residual;  ///< Structure containing information about the residual
    culaSparseTiming timing;          ///< Structure containing timing information about the iterative solver and associated preconditioners
    unsigned int code;                ///< Internal information code used to retrieve additional information
} culaSparseResult;

///
/// Reordering Types
///
/// Enumerates the reordering strategies available
///

typedef enum {
    culaSparseNoReordering,         ///< Do not do any reordering
    culaSparseAmdReordering,        ///< Reorder using the column approximate minimum degree ordering method 
    culaSparseSymamdReordering      ///< Reorder using the symmetric minimum degree ordering method
} culaSparseReordering;

///
/// Pattern Types
///
/// Enumerates the preconditioner patterns used by the approximate preconditioners.
///

typedef enum {
    culaSparseAPattern,         ///< Uses the same sparsity pattern as A
    culaSparseA2Pattern,        ///< Uses a sparsity pattern of A^2
    culaSparseA3Pattern,        ///< Uses a sparsity pattern of A^3
    culaSparseA4Pattern         ///< Uses a sparsity pattern of A^4
} culaSparseSparsityPattern;

///
/// Complex Data Types
///
/// When using the CUDA device interface 16 bit alignment is required. On other interfaces
/// higher performance can be attained by using aligned types.
///

/// Single Precision (16 byte aligned)
typedef union {
    struct {
        float real;
        float imag;
    };
    long dummy;
} culaSparseComplexFloat;

/// Double Precision (16 byte aligned)
typedef struct {
    double real;
    double imag;
} culaSparseComplexDouble;


///
/// Storage Format Options
/// 

typedef struct {
    /// Indicates whether the sparse indexing arrays are represented using 0
    /// (C/C++) or 1 (FORTRAN) based indexing.
    int indexing;
} culaSparseCooOptions;

typedef struct {
    /// Indicates whether the sparse indexing arrays are represented using 0
    /// (C/C++) or 1 (FORTRAN) based indexing.
    int indexing;
} culaSparseCsrOptions;

typedef struct {
    /// Indicates whether the sparse indexing arrays are represented using 0
    /// (C/C++) or 1 (FORTRAN) based indexing.
    int indexing;
} culaSparseCscOptions;

typedef struct {
    int reserved;
} culaSparseMatrixFreeOptions;

///
/// Platform Options
///
 
/// Host
typedef struct {
    /// Selects a reordering method to apply to the input matrix prior to any
    /// processing. This operation is applied to a copy of the matrix rather 
    /// than actual input data. 
    culaSparseReordering reordering;

    /// Specifies whether to perform extra checks to aid in debugging.
    int debug;
} culaSparseHostOptions;

/// Cuda
typedef struct {
    /// Selects a reordering method to apply to the input matrix prior to any
    /// processing. This operation is applied to a copy of the matrix rather 
    /// than actual input data. 
    culaSparseReordering reordering;

    /// Select the CUDA device for which all operations will take place. The 
    /// device identifiers are enumerated by the CUDA API.
    int deviceId;

    /// Optionally use the hybrid matrix optimization which may result in higher
    /// performance for many matrix structures at the expense of higher memory
    /// usage.
    int useHybridFormat;

    /// Specifies whether to perform extra checks to aid in debugging.
    int debug;
} culaSparseCudaOptions;

/// Cuda Device
typedef struct {
    /// Select the CUDA device for which all operations will take place. The 
    /// device identifiers are enumerated by the CUDA API. The device id must
    /// be the same as where memory allocation took place using cudaMalloc.
    int deviceId;

    /// Optionally use the hybrid matrix optimization which may result in higher
    /// performance for many matrix structures at the expense of higher memory
    /// usage.
    int useHybridFormat; 

    /// Specifies whether to perform extra checks to aid in debugging.
    int debug;
} culaSparseCudaDeviceOptions;

///
/// Solver Options
///
/// Provides solver specific options
///

typedef struct {
     /// Duplicate and transpose input matrix to avoid internal transpose operations
    int avoidTranspose;    
} culaSparseBicgOptions;

typedef struct {
    int reserved;
} culaSparseBicgstabOptions;

typedef struct {
    /// Restart value
    int l;                 
} culaSparseBicgstablOptions;

typedef struct {
    int reserved;
} culaSparseCgOptions;

typedef struct {
    /// Restart value
    int restart;
} culaSparseGmresOptions;

typedef struct {
    int reserved;
} culaSparseMinresOptions;

///
/// Preconditioner Options
///
/// Provides solver specific options
///

typedef struct {
    int unused;
} culaSparseEmptyOptions;

typedef struct {
    int reserved;
} culaSparseJacobiOptions;

typedef struct {
    /// Size of the square diagonal blocks
    int blockSize;
} culaSparseBlockjacobiOptions;

typedef struct {
    int reserved;
} culaSparseIlu0Options;

typedef struct {
    culaSparseSparsityPattern pattern;
    double dropTolerance;
} culaSparseAinvOptions;

typedef struct {
    culaSparseSparsityPattern pattern;
    double dropTolerance;
} culaSparseFainvOptions;

typedef struct {
    /// function pointer to single precision real application function
    void (*sApply)(char trans, const float* input, float* output, void* user);     

    /// function pointer to double precision real application function
    void (*dApply)(char trans, const double* input, double* output, void* user);

    /// function pointer to single precision complex application function
    void (*cApply)(char trans, const culaSparseComplexFloat* input, culaSparseComplexFloat* output, void* user);

    /// function pointer to double precision complex application function
    void (*zApply)(char trans, const culaSparseComplexDouble* input, culaSparseComplexDouble* output, void* user);
} culaSparseUserDefinedOptions;

#ifdef __cplusplus
} // extern "C"
#endif

#endif // __EMP_CULA_SPARSE_TYPES_H__
