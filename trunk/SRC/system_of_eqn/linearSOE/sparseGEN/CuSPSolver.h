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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/CuSPSolver.h,v $

// Written: neallee@tju.edu.cn 
// Modified from XZ Lu's CuSPSolver for GPL
// Created: 14/05
//
// Description: This file contains the definiation for CuSPSolver

#ifndef CuSPSolver_h
#define CuSPSolver_h

#include <SparseGenRowLinSolver.h>
#include <SparseGenRowLinSOE.h>
#include <SparseGenColLinSOE.h>
#include <SparseGenColLinSolver.h>

//typedef int(*CUSPSOLVE)(double* Aptr, double* Bptr, double* Xptr, int n, int nnz, int* rowPtr, int*  colInd, int maxInt, double relTol, int pre, int solv);

class CuSPSolver :
  public SparseGenRowLinSolver
{
public:
  CuSPSolver(void);
  CuSPSolver(int maxInt, double relTol, int pre, int solv);
  ~CuSPSolver(void);
  int solve(void);
  int setSize(void);

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  int setLinearSOE(SparseGenRowLinSOE &theLinearSOE);
  //int setLinearSOE(SparseGenColLinSOE &theLinearSOE);
  
private:
  int n;// order of matrix
  int nnz;// Number of no-zero ones in A;

  double *Aptr, *Bptr;//Matrix A, and B
  int *rowPtr, *colInd;

  double relTolerance;
  int maxInteration;
  int preCond;
  int solver;
  int single;
  int host;

  int error;

//  CUSPSOLVE SolveFunc;
};

#endif