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

// This header includes all of CULA's sparse functionality

#include "cula_status.h"

#define CULA_SPARSE_VERSION_NUMBER 4000
#define CULA_SPARSE_VERSION_STRING "S4 CUDA 5.0"

#ifdef __cplusplus
extern "C" {
#endif

//==============================================================================
// Initialization and shutdown
//==============================================================================

/**
 * @brief Initializes CULA Sparse
 * Must be called before using any other function.  Some functions have an
 * exception to this rule:  culaGetDeviceCount, culaSelectDevice, and version
 * query functions
 *
 * @return culaNoError on a successful initialization or a culaStatus enum
 * that specifies an error
 */
culaStatus culaSparseInitialize();

/**
 * @brief Shuts down CULA Sparse
 *
 * Must be called to deallocate CULA Sparse internal data
 */
void culaSparseShutdown();

//==============================================================================
// Configuration
//==============================================================================

/**
 * @struct culaIterativeConfig
 * @brief Contains information that steers execution of iterative functions
 */
typedef struct
{
    // Indicates whether the sparse indexing arrays are represented using 0
    // (C/C++) or 1 (FORTRAN) based indexing.
    int indexing;

    // The tolerance is the point at which a lower residual will cause the
    // solver to determine that the solution has converged.
    double tolerance;

    // The maximum number of iterations that the solver will attempt.
    int maxIterations;

    // The maximum time, in seconds, at which the solver will not begin a new    
    // iteration. If set to 0, the solver will run untill convergence or the 
    // maximum number of iterations has been reached.
    double maxRuntime;

    // This parameter provides the means for a user to capture the residual at
    // each iteration. The specified array must be at least maxIter in length.
    // This parameter may be NULL if these quantities are not desired.
    double* residualVector;

    // Indicates whether the 'x' vector in iterative solves should be used as
    // given or ignored. When ignored, the 'x' vector is considered a zero.
    int useInitialResultVector;

    // Indicates whether the 'x' vector in iterative solves should return the
    // final answer or the best answer and its associated iteration number in 
    // the case of non-convergence. 
    int useBestAnswer;

    // Indicates whether to check whether the iterative solve stagnated. This
    // option defaults to on; turning this option off will increase performance
    // if a problem is certain not to stagnate.
    int useStagnationCheck;

    // Specifies whether to perform extra checks to aid in debugging
    int debug;
}culaIterativeConfig;

/**
 * @brief Constructs a config structure, initializing it to default values.
 *
 * @param config Pointer to the culaIterativeConfig struct that will be
 * initialized by this routine.
 *
 * @return culaNoError on success or culaArgumentError if a parameter is invalid
 */
culaStatus culaIterativeConfigInit(culaIterativeConfig* config);

/**
 * @brief Associates an iterative config structure with a readable config report
 *
 * @param config Pointer to the culaIterativeConfig struct that will be
 * analyzed by this routine
 * @param buf Pointer to a buffer into which information will be printed
 * (recommend size 256 or greater)
 * @param bufsize The size of buf, printed information will not exceed bufsize
 *
 * @return culaNoError on a successful config report or culaArgumentError on an
 * invalid argument to this function
 */
culaStatus culaIterativeConfigString(const culaIterativeConfig* config, char* buf, int bufsize);

//==============================================================================
// Results
//==============================================================================

/**
 * @struct culaIterativeFlag
 * @brief Enumerates error conditions for iterative functions
 */
typedef enum
{
    culaConverged,            // The solve converged successfully
    culaMaxItReached,         // Maximum iterations reached without convergence
    culaMaxTimeReached,       // Maximum execution time reached without convergance
    culaPreconditionerFailed, // The specified preconditioner failed
    culaStagnation,           // The iterative solve stagnated
    culaScalarOutOfRange,     // A scalar value was out of range
    culaReorderingError,      // Matrix reordering failed
    culaUnknownIterationError // An unknown iteration error was encountered
}culaIterativeFlag;

/**
 * @struct culaIterativeResidual
 * @brief Contains information about the residual of an iterative function
 */
typedef struct
{
    double relative;     // The relative residual obtained by the iterative solver when computation has completed or halted
    double* byIteration; // If requested, the residual at every step of iteration
}culaIterativeResidual;

/**
 * @struct culaIterativeTiming
 * @brief Contains timing for the execution of an iterative function
 */
typedef struct
{
    double solve;           // Time, in seconds, the solve portion of the iterative solver took to complete
    double preconditioner;  // Time, in seconds, the preconditioner generative portion of the iterative solver took to complete
    double overhead;        // Time, in seconds, of overhead needed by the iterative solver; includes memory transfers to-and-from the GPU
    double total;           // Time, in seconds, the entire iterative solver took to complete
}culaIterativeTiming;

/**
 * @struct culaIterativeResult
 * @brief Collects information about the execution of an iterative function
 *
 * This structure is only populated with valid information when culaNoError or
 * culaDataError is returned from a CULA Sparse function
 */
typedef struct
{
    culaIterativeConfig config;      // Matches the configuration that a given solver was called with
    unsigned long long code;         // Internal information code
    culaIterativeFlag flag;          // Enumeration containing information about the success or failure of the iterative solver
    int iterations;                  // Number of iterations taken by the iterative solver
    culaIterativeResidual residual;  // Structure containing information about the residual
    culaIterativeTiming timing;      // Structure containing timing information about the iterative solver and associated preconditioners
}culaIterativeResult;

/**
 * @brief Associates an iterative result structure with a readable result result
 *
 * @param e A culaStatus error code
 * @param i A pointer to a culaIterativeResult structure
 * @param buf Pointer to a buffer into which information will be printed
 * (recommend size 256 or greater)
 * @param bufsize The size of buf, printed information will not exceed bufsize
 *
 * @return culaNoError on a successful result report or culaArgumentError on an
 * invalid argument to this function
 */
culaStatus culaIterativeResultString(const culaIterativeResult* result, char* buf, int bufsize);

/**
 * @struct culaReordering
 * @brief Enumerates a reordering strategy
 */
typedef enum
{
    culaNoReordering,    // Do not do any reordering
    culaAmdReordering,   // Reorder using the column approximate minimum degree ordering method (AMD)
    culaSymamdReordering // Reorder using the symmetric minimum degree ordering method (SYMAMD)
}culaReordering;

/**
 * @struct culaEmptyOptions
 * @brief Used for parameters that need no configuration
 */
typedef struct
{
    int unused;
}culaEmptyOptions;
culaStatus culaEmptyOptionsInit(culaEmptyOptions* opts);

//==============================================================================
// Solver Options
//==============================================================================

typedef struct
{
    int avoidTranspose;     // Duplicate input matrix to avoid internal transpose operations
}culaBicgOptions;
culaStatus culaBicgOptionsInit(culaBicgOptions* opts);

typedef struct
{
    int reserved;
}culaBicgstabOptions;
culaStatus culaBicgstabOptionsInit(culaBicgstabOptions* opts);

typedef struct
{
    int l;
}culaBicgstablOptions;
culaStatus culaBicgstablOptionsInit(culaBicgstablOptions* opts);

typedef struct
{
    int reserved;
}culaCgOptions;
culaStatus culaCgOptionsInit(culaCgOptions* opts);

typedef struct
{
    int restart;
}culaGmresOptions;
culaStatus culaGmresOptionsInit(culaGmresOptions* opts);

typedef struct
{
    int reserved;
}culaMinresOptions;
culaStatus culaMinresOptionsInit(culaMinresOptions* opts);

//==============================================================================
// Preconditioner Options
//==============================================================================

typedef struct
{
    int reserved;
}culaJacobiOptions;
culaStatus culaJacobiOptionsInit(culaJacobiOptions* opts);

typedef struct
{
    int blockSize;
}culaBlockjacobiOptions;
culaStatus culaBlockjacobiOptionsInit(culaBlockjacobiOptions* opts);

typedef struct
{
    culaReordering reordering;
}culaIlu0Options;
culaStatus culaIlu0OptionsInit(culaIlu0Options* opts);

//==============================================================================
// Solvers
//==============================================================================

culaStatus culaDcscCg(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrCg(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooCg(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscCg(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrCg(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooCg(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscCgIlu0(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrCgIlu0(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooCgIlu0(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscCgIlu0(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrCgIlu0(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooCgIlu0(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscCgJacobi(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrCgJacobi(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooCgJacobi(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscCgJacobi(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrCgJacobi(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooCgJacobi(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscCgBlockjacobi(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrCgBlockjacobi(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooCgBlockjacobi(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscCgBlockjacobi(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrCgBlockjacobi(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooCgBlockjacobi(const culaIterativeConfig* config, const culaCgOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscBicg(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrBicg(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooBicg(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscBicg(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrBicg(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooBicg(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscBicgIlu0(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrBicgIlu0(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooBicgIlu0(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscBicgIlu0(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrBicgIlu0(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooBicgIlu0(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscBicgJacobi(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrBicgJacobi(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooBicgJacobi(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscBicgJacobi(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrBicgJacobi(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooBicgJacobi(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscBicgBlockjacobi(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrBicgBlockjacobi(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooBicgBlockjacobi(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscBicgBlockjacobi(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrBicgBlockjacobi(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooBicgBlockjacobi(const culaIterativeConfig* config, const culaBicgOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscBicgstab(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrBicgstab(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooBicgstab(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscBicgstab(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrBicgstab(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooBicgstab(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscBicgstabIlu0(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrBicgstabIlu0(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooBicgstabIlu0(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscBicgstabIlu0(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrBicgstabIlu0(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooBicgstabIlu0(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscBicgstabJacobi(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrBicgstabJacobi(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooBicgstabJacobi(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscBicgstabJacobi(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrBicgstabJacobi(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooBicgstabJacobi(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscBicgstabBlockjacobi(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrBicgstabBlockjacobi(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooBicgstabBlockjacobi(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscBicgstabBlockjacobi(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrBicgstabBlockjacobi(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooBicgstabBlockjacobi(const culaIterativeConfig* config, const culaBicgstabOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscBicgstabl(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrBicgstabl(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooBicgstabl(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscBicgstabl(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrBicgstabl(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooBicgstabl(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscBicgstablIlu0(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrBicgstablIlu0(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooBicgstablIlu0(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscBicgstablIlu0(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrBicgstablIlu0(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooBicgstablIlu0(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscBicgstablJacobi(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrBicgstablJacobi(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooBicgstablJacobi(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscBicgstablJacobi(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrBicgstablJacobi(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooBicgstablJacobi(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscBicgstablBlockjacobi(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrBicgstablBlockjacobi(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooBicgstablBlockjacobi(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscBicgstablBlockjacobi(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrBicgstablBlockjacobi(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooBicgstablBlockjacobi(const culaIterativeConfig* config, const culaBicgstablOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscGmres(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrGmres(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooGmres(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscGmres(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrGmres(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooGmres(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscGmresIlu0(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrGmresIlu0(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooGmresIlu0(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscGmresIlu0(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrGmresIlu0(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooGmresIlu0(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscGmresJacobi(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrGmresJacobi(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooGmresJacobi(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscGmresJacobi(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrGmresJacobi(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooGmresJacobi(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscGmresBlockjacobi(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrGmresBlockjacobi(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooGmresBlockjacobi(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscGmresBlockjacobi(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrGmresBlockjacobi(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooGmresBlockjacobi(const culaIterativeConfig* config, const culaGmresOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscMinres(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrMinres(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooMinres(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscMinres(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrMinres(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooMinres(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaEmptyOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscMinresIlu0(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrMinresIlu0(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooMinresIlu0(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscMinresIlu0(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrMinresIlu0(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooMinresIlu0(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaIlu0Options* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscMinresJacobi(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrMinresJacobi(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooMinresJacobi(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscMinresJacobi(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrMinresJacobi(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooMinresJacobi(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaJacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

culaStatus culaDcscMinresBlockjacobi(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* rowInd, const int* colPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcsrMinresBlockjacobi(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowPtr, double* x, const double* b, culaIterativeResult* result);
culaStatus culaDcooMinresBlockjacobi(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const double* a, const int* colInd, const int* rowInd, double* x, const double* b, culaIterativeResult* result);
culaStatus culaZcscMinresBlockjacobi(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* rowInd, const int* colPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcsrMinresBlockjacobi(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowPtr, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);
culaStatus culaZcooMinresBlockjacobi(const culaIterativeConfig* config, const culaMinresOptions* solverOpts, const culaBlockjacobiOptions* precondOpts, int n, int nnz, const culaDoubleComplex* a, const int* colInd, const int* rowInd, culaDoubleComplex* x, const culaDoubleComplex* b, culaIterativeResult* result);

//==============================================================================
// Auxiliary
//==============================================================================

/**
 * @brief Reports the CUSPARSE_VERSION that the running version of CULA was
 * compiled against, which indicates the minimum version of CUSPARSE that is
 * required to use this library
 *
 * @return An integer in the format XXXYY where XXX is the major version number
 * and YY is the minor version number of CUSPARSE that this version of CULA was
 * compiled against.  On error, a 0 is returned.
 */
culaVersion culaGetCusparseMinimumVersion();

/**
 * @brief Reports the version of the CUSPARSE runtime that operating system
 * linked against when the program was loaded
 *
 * This function is useful for discovering issues with an incorrectly set linker
 * path.
 *
 * @return An integer in the format XXXYY where XXX is the major version number
 * and YY is the minor version number of the CUSPARSE runtime.  On error, a 0 is
 * returned.
 */
culaVersion culaGetCusparseRuntimeVersion();

#ifdef __cplusplus
} // extern "C"
#endif

#endif // __EMP_CULA_SPARSE_H__

