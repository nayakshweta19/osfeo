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

#include "CuSPSolver.h"
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
#include <cusp/precond/aggregation/smoothed_aggregation.h>
#include <cusp/precond/detail/ainv.inl>
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
  if (single == 1) {
    delete Xsingle;
    delete Bsingle;
    delete Asingle;
  }
}

int
CuSPSolver::setSize()
{
  int n = theSOE->size;
  int nnz = theSOE->nnz;

  if (single == 1) {
    Xsingle = new float[n];
    Bsingle = new float[n];
    Asingle = new float[nnz];
  }
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

  double *Xptr = theSOE->X;
  double *Bptr = theSOE->B;
  double *Aptr = theSOE->A;

  int *rowPtr = theSOE->rowStartA;
  int *colInd = theSOE->colA;

  // allocate storage for matrix A with NNZ non-zeros
  cusp::csr_matrix<int, float, cusp::host_memory> AA(n, n, nnz);

  // initialize matrix entries on host
  AA.row_offsets = rowPtr;
  AA.column_indices = colInd;
  AA.values = Aptr;

  // copy to the device
  cusp::csr_matrix<int, float, cusp::device_memory> AA = A;

  // allocate storage for solution (x) and right hand side (b)
  cusp::array1d<float, cusp::device_memory> x = Xptr;
  cusp::array1d<float, cusp::device_memory> b = Bptr;

  // set stopping criteria:
  //  iteration_limit    = 100
  //  relative_tolerance = 1e-6
  cusp::verbose_monitor<float> monitor(b, 1000, 1e-6);
  //cusp::default_monitor<float> monitor(b, 1000, 1e-6);
  //cusp::convergence_monitor<float> monitor(b, 1000, 1e-6);

  // set preconditioner (identity)
  switch (preCond) // 0 - none; 1 - diagonal; 2 - ainv
  {
//  case 0:
//  {
//    cusp::identity_operator<float, cusp::device_memory> M = AA;
//    break;
//  }
  case 1:
  {
    cusp::precond::diagonal<float, cusp::device_memory> M(AA);
    break;
  }
  case 2:
  {
    cusp::precond::aggregation:smoothed_aggregation<int, float, cusp::device_memory> M(AA);
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
    cusp::krylov::bicgstab(AA, AA, x, b, monitor, M, M);
    break;
  }
  case 1:
  {
    cusp::krylov::bicgstab(AA, x, b, monitor, M);
    break;
  }
  case 2:
  {
    cusp::krylov::cg(AA, x, b, monitor, M);
    break;
  }
  case 3:
  {
    cusp::krylov::gmres(AA, x, b, restart, monitor, M);
    break;
  }
  default:
  {
    opserr << "CuSPSolver::solve() - the wrong solver type defined." << endln;
    break;
  }
  }
//  int errorcode = SolveFunc(Aptr, Bptr, Xptr, n, nnz, rowPtr, colInd, maxInteration, relTolerance, preCond, solver);
//  if (errorcode != 0)
//  {
//    switch (errorcode) {
//    case -1:opserr << "Wrong Preconditioner! Please check the TCL file.\n"; break;
//    case -2:opserr << "Wrong Solver! Please check the TCL file.\n"; break;
//    }
//    return -1;
//  }

  return 0;
}
