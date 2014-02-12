#ifndef __EMP_CULA_STATUS_H__
#define __EMP_CULA_STATUS_H__

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

#include "cula_types.h"

#ifdef __cplusplus
extern "C" {
#endif


/**
 * @brief Enumerates error conditions
 */
typedef enum
{
    culaNoError,                       // No error
    culaNotInitialized,                // CULA has not been initialized
    culaNoHardware,                    // No hardware is available to run
    culaInsufficientRuntime,           // CUDA runtime or driver is not supported
    culaInsufficientComputeCapability, // Available GPUs do not support the requested operation
    culaInsufficientMemory,            // There is insufficient memory to continue
    culaFeatureNotImplemented,         // The requested feature has not been implemented
    culaArgumentError,                 // An invalid argument was passed to a function
    culaDataError,                     // An operation could not complete because of singular data
    culaBlasError,                     // A blas error was encountered
    culaRuntimeError,                  // A runtime error has occurred
    culaBadStorageFormat,              // An invalid storage format was used for a parameter
    culaInvalidReferenceHandle,        // An invalid reference handle was passed to a function
    culaUnspecifiedError               // An unspecified internal error has occurred
}culaStatus;


/**
 * @brief Initializes CULA
 * Must be called before using any other function.  Some functions have an
 * exception to this rule:  culaGetDeviceCount, culaSelectDevice, and version
 * query functions
 *
 * @return culaNoError on a successful initialization or a culaStatus enum
 * that specifies an error
 */
culaStatus culaInitialize();


/**
 * @brief Shuts down CULA
 *
 * Must be called to deallocate CULA internal data
 */
void culaShutdown();


/**
 * @brief Returns the last status code returned from a CULA function
 *
 * @return The last CULA status code
 */
culaStatus culaGetLastStatus();


/**
 * @brief Associates a culaStatus enum with a readable error string
 *
 * @param e A culaStatus error code
 *
 * @return A string that corresponds with the specified culaStatus enum
 */
const char* culaGetStatusString(culaStatus e);


/**
 * @brief Returns the culaStatus name as a string
 *
 * @param e A culaStatus error code
 *
 * @return A string that corresponds with the specified culaStatus enum
 */
const char* culaGetStatusAsString(culaStatus e);


/**
 * @brief This function is used to provide extended functionality that LAPACK's
 * info parameter typically provides
 *
 * @return Extended information about the last error or zero if it is
 * unavailable
 */
culaInfo culaGetErrorInfo();


/**
 * @brief Associates a culaStatus and culaInfo with a readable error string
 *
 * @param e A culaStatus error code
 * @param i An culaInfo error code
 * @param buf Pointer to a buffer into which information will be printed
 * @param bufsize The size of buf, printed information will not exceed bufsize
 *
 * @return culaNoError on a successful error report or culaArgumentError on an
 * invalid argument to this function
 */
culaStatus culaGetErrorInfoString(culaStatus e, culaInfo i, char* buf, int bufsize);


/**
 * @brief Releases any memory buffers stored internally by CULA
 */
void culaFreeBuffers();


/**
 * @brief Reports the version number of CULA
 *
 * This function provides a way to query the version of CULA in use at runtime.
 * Use the CULA_VERSION_NUMBER macro for a compile time query.
 *
 * @return An integer in the format XXXYY where XXX is the major version number
 * and YY is the minor version number of CULA.  On error, a 0 is returned.
 */
culaVersion culaGetVersion();


/**
 * @brief Reports the CUDA_VERSION that the running version of CULA was compiled
 * against, which indicates the minimum version of CUDA that is required to use
 * this library
 *
 * @return An integer in the format XXXYY where XXX is the major version number
 * and YY is the minor version number of CUDA that this version of CULA was
 * compiled against.  On error, a 0 is returned.
 */
culaVersion culaGetCudaMinimumVersion();


/**
 * @brief Reports the version of the CUDA runtime that the operating system
 * linked against when the program was loaded
 *
 * This function is useful for discovering issues with an incorrectly set linker
 * path.
 *
 * @return An integer in the format XXXYY where XXX is the major version number
 * and YY is the minor version number of the CUDA runtime.  On error, a 0 is
 * returned.
 */
culaVersion culaGetCudaRuntimeVersion();


/**
 * @brief Reports the version of the CUDA driver installed on the system
 *
 * @return An integer in the format XXXYY where XXX is the major version number
 * and YY is the minor version number of the CUDA driver currently installed on
 * the system. If no driver is installed, a 0 is returned.
 */
culaVersion culaGetCudaDriverVersion();


/**
 * @brief Reports the CUBLAS_VERSION that the running version of CULA was
 * compiled against, which indicates the minimum version of CUBLAS that is
 * required to use this library
 *
 * @return An integer in the format XXXYY where XXX is the major version number
 * and YY is the minor version number of CUBLAS that this version of CULA was
 * compiled against.  On error, a 0 is returned.
 */
culaVersion culaGetCublasMinimumVersion();


/**
 * @brief Reports the version of the CUBLAS runtime that operating system linked
 * against when the program was loaded
 *
 * This function is useful for discovering issues with an incorrectly set linker
 * path.
 *
 * @return An integer in the format XXXYY where XXX is the major version number
 * and YY is the minor version number of the CUBLAS runtime.  On error, a 0 is
 * returned.
 */
culaVersion culaGetCublasRuntimeVersion();


#ifdef __cplusplus
} // extern "C"
#endif

#endif  // __EMP_CULA_STATUS_H__

