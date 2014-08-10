/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.0 $
// $Date: 2014-07-14 20:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/CuSPSolver.cpp,v $

// Written: neallee@tju.edu.cn 
// Modified from XZ Lu's CuSPSolver for GPL
// Created: 14/05
//
// Description: This file contains the implementation for CuSPSolver

#include <CuSPSolver.h>
#include <cusp/blas.h>
#include <cusp/format.h>
#include <cusp/linear_operator.h>
#include <cusp/csr_matrix.h>
#include <cusp/monitor.h>
#include <cusp/exception.h>
#include <cusp/krylov/bicg.h>
#include <cusp/krylov/cg.h>
#include <cusp/krylov/gmres.h>
#include <cusp/krylov/bicgstab.h>
#include <cusp/precond/ainv.h>
#include <cusp/precond/aggregation/smoothed_aggregation.h>
#include <cusp/precond/diagonal.h>

#include <Windows.h>

CuSPSolver::CuSPSolver(void) :SparseGenRowLinSolver(SOLVER_TAGS_CuSP)
{
  single = 0;
  error = 0;

//  HINSTANCE hDLL = LoadLibrary("CuSPSolver.dll");
//  if (hDLL) {
//    SolveFunc = (CUSPSOLVE)GetProcAddress(hDLL, "CuSPSolve");
//
//    if (!SolveFunc) {
//      error = 1;
//      return;
//    }
//  }
//  else {
//    error = 2;
//    return;
//  }

  this->maxInteration = 100000;
  this->relTolerance = 1e-6;
  this->preCond = 0;  // 0 - none; 1 - diagonal; 2 - ainv
  this->solver = 0; // 0 - bicg; 1 - bicgstab; 2 - cg; 3 - gmres
}

CuSPSolver::CuSPSolver(int maxInt, double relTol, int pre, int solv) :SparseGenRowLinSolver(SOLVER_TAGS_CuSP)
{
  single = 0;
  error = 0;

//  HINSTANCE hDLL = LoadLibrary("CuSPSolver.dll");
//  if (hDLL) {
//    SolveFunc = (CUSPSOLVE)GetProcAddress(hDLL, "CuSPSolve");
//
//    if (!SolveFunc) {
//      error = 1;
//      return;
//    }
//  }
//  else {
//    error = 2;
//    return;
//  }

  this->maxInteration = maxInt;
  this->relTolerance = relTol;
  this->preCond = pre;  // 0 - none; 1 - diagonal; 2 - ainv
  this->solver = solv;  // 0 - bicg; 1 - bicgstab; 2 - cg; 3 - gmres
}

CuSPSolver::~CuSPSolver(void)
{
}

int
CuSPSolver::setSize()
{
  n = theSOE->size;
  nnz = theSOE->nnz;

  return 0;
}

int
CuSPSolver::setLinearSOE(SparseGenRowLinSOE &theLinearSOE)
{
  theSOE = &theLinearSOE;
  return 0;
}

int
CuSPSolver::sendSelf(int cTAg, Channel &theChannel)
{
  // doing nothing
  return 0;
}

int
CuSPSolver::recvSelf(int cTag,
Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  // nothing to do
  return 0;
}

int 
CuSPSolver::solve(void)
{
  if (theSOE == 0) {
    opserr << "WARNING SuperLU::solve(void)- ";
    opserr << " No LinearSOE object has been set\n";
    return -1;
  }

  n = theSOE->size;
  // check for quick return
  if (n == 0)
    return 0;

  nnz = theSOE->nnz;

  Bptr = theSOE->B;
  Aptr = theSOE->A;

  rowPtr = theSOE->rowStartA;
  colInd = theSOE->colA;

  // allocate storage for matrix A with NNZ non-zeros
  cusp::csr_matrix<int, double, cusp::host_memory> hostA(n, n, nnz);

  // initialize matrix entries on host
  for (int l = 0; l < nnz; l++) {
    hostA.values[l] = Aptr[l];
    hostA.column_indices[l] = colInd[l];
  }
  for (int k = 0; k < n+1; k++) {
    hostA.row_offsets[k] = rowPtr[k];
  }

  // copy to the device
  cusp::csr_matrix<int, double, cusp::device_memory> devA = hostA;

  // allocate storage for solution (x) and right hand side (b)
  cusp::array1d<double, cusp::device_memory> x(n, 0.0);
  cusp::array1d<double, cusp::device_memory> b(n);
  for (int k = 0; k < n; k++) {
//    x[k] = Xptr[k];
    b[k] = Bptr[k];
  }

  // set stopping criteria:
  //  iteration_limit    = 100
  //  relative_tolerance = 1e-6
  cusp::verbose_monitor<double> monitor(b, 1000, 1e-6);
  //cusp::default_monitor<double> monitor(b, 1000, 1e-6);
  //cusp::convergence_monitor<double> monitor(b, 1000, 1e-6);

  // set preconditioner (identity)
  cusp::identity_operator<double, cusp::device_memory> M(n, n);
  switch (preCond) // 0 - none; 1 - diagonal; 2 - ainv
  {
  case 0:
  {    
    break;
  }
  case 1:
  {
    cusp::precond::diagonal<double, cusp::device_memory> M(devA);
    break;
  }
  case 2:
  {
    //cusp::precond::scaled_bridson_ainv<double, cusp::device_memory> M(devA, 0.1);
    //cusp::precond:scaled_bridson_ainv<ValueType, MemorySpace> M(devA, 0, 10);
    cusp::precond::bridson_ainv<double, cusp::device_memory> M(devA, 0, -1, true, 2);
    break;
  }
  default:
  {
    opserr << "CuSPSolver::solve() - the wrong preCond type defined." << endln;
    break;
  }
  }
  
  // solve the linear system A x = b
  switch (solver) // 0 - bicg; 1 - bicgstab; 2 - cg; 3 - gmres
  {
  case 0:
  {
    // because both A and M are hermitian we can use 
    // them for their own conjugate transpose
    if (preCond != 0) cusp::krylov::bicgstab(devA, x, b, monitor, M);
    else cusp::krylov::bicgstab(devA, x, b);
    break;
  }
  case 1:
  {
    cusp::krylov::bicgstab(devA, x, b, monitor, M);
    break;
  }
  case 2:
  {
    cusp::krylov::cg(devA, x, b, monitor, M);
    break;
  }
  case 3:
  {
    int restart = 50;
    cusp::krylov::gmres(devA, x, b, restart, monitor, M);
    break;
  }
  default:
  {
    opserr << "CuSPSolver::solve() - the wrong solver type defined." << endln;
    break;
  }
  }

  for (int k = 0; k < n; k++) {
    theSOE->X[k] = x[k];
  }

  return 0;
}
