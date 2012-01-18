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

// $Revision: 1.3 $
// $Date: 2003/04/02 22:02:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/fullGEN/FullGenRitzLapackSolver.cpp,v $


// File: ~/system_of_eqn/linearSOE/FullGEN/FullGenRitzLapackSolver.h
//
// Written: fmk 
// Created: Tue Sep 26 16:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// FullGenRitzLapackSolver. It solves the FullGenLinSOE object by calling
// Lapack routines.
//
// What: "@(#) FullGenRitzLapackSolver.h, revA"

#include <FullGenRitzLapackSolver.h>
#include <FullGenRitzSOE.h>
#include <f2c.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <AnalysisModel.h>
#include <DOF_GrpIter.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <FE_Element.h>
#include <Integrator.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#define SOLVER_TAGS_FullGenRitzLapackSolver     25

FullGenRitzLapackSolver::FullGenRitzLapackSolver()
	:FullGenRitzSOESolver(SOLVER_TAGS_FullGenRitzLapackSolver),iPiv(0),sizeIpiv(0),
	numRitz(0), eigenvalue(0), eigenvector(0),
	sortingID(0), eigenV(0)
{

}

FullGenRitzLapackSolver::~FullGenRitzLapackSolver()
{
	if (iPiv != 0)
		delete [] iPiv;
	if (eigenvalue != 0)
        delete [] eigenvalue;
    if (eigenvector != 0)
        delete [] eigenvector;
    if (sortingID != 0)
        delete [] sortingID;
    if (eigenV != 0)
        delete eigenV;
}


#ifdef _WIN32
extern "C" int  DGESV(int *N, int *NRHS, double *A, int *LDA, 
					  int *iPiv, double *B, int *LDB, int *INFO);

extern "C" int  DGETRS(char *TRANS,
					   int *N, int *NRHS, double *A, int *LDA, 
					   int *iPiv, double *B, int *LDB, int *INFO);
extern "C" int DGGEV(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                     double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                     double *BETA, double *VL, int *LDVL, double *VR,
                     int *LDVR, double *WORK, int *LWORK, int *INFO);

#else
extern "C" int dgesv_(int *N, int *NRHS, double *A, int *LDA, int *iPiv, 
					  double *B, int *LDB, int *INFO);

extern "C" int dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, 
					   int *iPiv, double *B, int *LDB, int *INFO);		       
extern "C" int dggev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                      double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                      double *BETA, double *VL, int *LDVL, double *VR,
                      int *LDVR, double *WORK, int *LWORK, int *INFO);

#endif

int
FullGenRitzLapackSolver::solve(int nRitz)
{
	if (theSOE == 0) {
		opserr << "WARNING FullGenRitzLapackSolver::solve(void)- ";
		opserr << " No LinearSOE object has been set\n";
		return -1;
	}

	int n = theSOE->size;

	// set the number of eigenvalues
    numRitz = nRitz;

	// check for quick return
	if (n == 0)
		return 0;

	// check iPiv is large enough
	if (sizeIpiv < n) {
		opserr << "WARNING FullGenRitzLapackSolver::solve(void)- ";
		opserr << " iPiv not large enough - has setSize() been called?\n";
		return -1;
	}	
	
    // output information
    int info = 0;
	double *Aptr = theSOE->A;
	double *Xptr = theSOE->X;
	double *Bptr = theSOE->B;
	double *Mptr = theSOE->M;
	/*int *iPIV = iPiv;
	
	//print K matrix
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++)
			opserr <<  *(Aptr+i*n+j) << " ";
		opserr << endln;
	}
	// print M matrix
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++)
			opserr <<  *(Mptr+i*n+j) << " ";
		opserr << endln;
	}
	// print B vector
	for (int i=0; i<n; i++) {
		opserr <<  *(Bptr+i) << " ";
		opserr << endln;
	}
	/*print X vector
	for (int i=0; i<n; i++) {
		opserr <<  *(Xptr+i) << " ";
	opserr << endln;
	}*/

    if (eigenvalue != 0)
        delete [] eigenvalue;

    eigenvalue = new double [n];

	// allocate memory for sorting index array
    if (sortingID != 0)
        delete [] sortingID;
    sortingID = new int [n];

    // allocate memory for right eigenvectors
    if (eigenvector != 0)
        delete [] eigenvector;
    eigenvector = new double [n*n];
	
	 // do not compute left eigenvalues and eigenvectors
    char *jobvl = "N";

    // compute right eigenvalues and eigenvectors
    char *jobvr = "V";

	// leading dimension of A
    int ldA = n;

    // leading dimension of M
    int ldM = n;

	// leading dimension of B
	int ldB = n;

    // allocate memory for eigenvalues
    double *alphaR = new double [n];
    double *alphaI = new double [n];
    double *beta   = new double [n];

    if (eigenvalue != 0)
        delete [] eigenvalue;

    eigenvalue = new double [n];

    // allocate memory for sorting index array
    if (sortingID != 0)
        delete [] sortingID;
    sortingID = new int [n];

    // dummy left eigenvectors
    double vl[1];

    // leading dimension of dummy left eigenvectors
    int ldvl = 1;

	// leading dimension of right eigenvectors
    int ldvr = n;

    // dimension of the workspace array
    int lwork = n*(8+64);

    // allocate memory for workspace array
    double *work = new double [lwork];

	// first copy B into X
	for (int i=0; i<n; i++)
	*(Xptr++) = *(Bptr++);
	Xptr = theSOE->X;

    // call the LAPACK eigenvalue subroutine
#ifdef _WIN32
    DGGEV(jobvl, jobvr, &n, Aptr, &ldA, Mptr, &ldM, alphaR, alphaI, beta,
          vl, &ldvl, eigenvector, &ldvr, work, &lwork, &info);
#else
    dggev_(jobvl, jobvr, &n, Aptr, &ldA, Mptr, &ldM, alphaR, alphaI, beta,
           vl, &ldvl, eigenvector, &ldvr, work, &lwork, &info);
#endif


    // check if successfull
	if (info != 0) {
		opserr << "WARNING FullGenRitzLapackSolver::solve()";
		opserr << " - lapack solver failed - " << info << " returned\n";
		return -info;
	}
	theSOE->factored = true;
	
	for (int i=0; i<n; i++) {
        double mag = sqrt(alphaR[i]*alphaR[i]+alphaI[i]*alphaI[i]);
        if (mag*DBL_EPSILON < fabs(beta[i])) {
            if (alphaI[i] == 0.0) {
                eigenvalue[i] = alphaR[i]/beta[i];
            }
            else {
                eigenvalue[i] = -mag/beta[i];
                opserr << "FullGenEigenSolver::solve() - the eigenvalue "
                    << i+1 << " is complex with magnitude "
                    << -eigenvalue[i] << endln;
            }
        }
        else {
	  eigenvalue[i] = DBL_MAX;
	}
        sortingID[i] = i;
    }


    // sort eigenvalues in ascending order and return sorting ID
    this->sort(n, eigenvalue, sortingID);

    for (int i=0; i<numRitz; i++) {
      if (eigenvalue[i] == DBL_MAX) {
	  opserr << "FullGenEigenSolver::solve() - the eigenvalue "
		 << i+1 << " is numerically undetermined or infinite\n";
      }
    }

    int lworkOpt = (int) work[0];
    if (lwork < lworkOpt) {
        opserr << "FullGenEigenSolver::solve() - optimal workspace size "
                << lworkOpt << " is larger than provided workspace size "
                << lwork << " consider increasing workspace\n";
    }

    // clean up the memory
    delete [] alphaR;
    delete [] alphaI;
    delete [] beta;
    delete [] work;

	return 0;
}



const Vector& FullGenRitzLapackSolver::getRitzvector(int mode)
{
	if (mode <= 0 || mode > numRitz) {
		opserr << "FullGenEigenSolver::getEigenVector() - mode "
			<< mode << " is out of range (1 - " << numRitz << ")\n";
		eigenV->Zero();
		return *eigenV;
	}

	int size = theSOE->size;
	int index = 0; //size*sortingID[mode-1];

	if (eigenvector != 0) {
		for (int i=0; i<size; i++) {
			(*eigenV)[i] = eigenvector[index++];
		}	
	}
	else {
		opserr << "FullGenEigenSolver::getEigenvector() - "
			<< "eigenvectors not computed yet\n";
		eigenV->Zero();
	}      

	return *eigenV;
}


double FullGenRitzLapackSolver::getRitzvalue(int mode)
{
	if (mode <= 0 || mode > numRitz) {
		opserr << "FullGenEigenSolver::getEigenvalue() - mode " 
			<< mode << " is out of range (1 - " << numRitz << ")\n";
		return 0.0;
	}

	if (eigenvalue != 0) {
		return eigenvalue[mode-1];
	}
	else {
		opserr << "FullGenEigenSolver::getEigenvalue() - "
			<< "eigenvalues not yet computed\n";
		return 0.0;
	}      
}


int
FullGenRitzLapackSolver::setSize()
{
	int n = theSOE->size;
	if (n > 0) {
		if (sizeIpiv < n) {
			if (iPiv != 0)
				delete [] iPiv;
			iPiv = new int[n];		
			if (iPiv == 0) {
				opserr << "WARNING FullGenRitzLapackSolver::setSize()";
				opserr << " - ran out of memory\n";
				return -1;
			}		
			sizeIpiv = n;
		}
	} else if (n == 0)
		return 0;
	else {
		opserr << "WARNING FullGenRitzLapackSolver::setSize()";
		opserr << " - ran out of memory\n";
		return -1;	
	}

	if (eigenV == 0 || eigenV->Size() != n) {
		if (eigenV != 0)
			delete eigenV;

		eigenV = new Vector(n);
		if (eigenV == 0 || eigenV->Size() != n) {
			opserr << "FullGenEigenSolver::setSize() ";
			opserr << " - ran out of memory for eigenVector of size ";
			opserr << theSOE->size << endln;
			return -2;	    
		}
	}

	return 0;
}

int
FullGenRitzLapackSolver::sendSelf(int commitTag,
								 Channel &theChannel)

{
	// nothing to do
	return 0;
}

int
FullGenRitzLapackSolver::recvSelf(int commitTag,
								 Channel &theChannel, 
								 FEM_ObjectBroker &theBroker)
{
	// nothing to do
	return 0;
}


void FullGenRitzLapackSolver::sort(int length, double *x, int *id)
{
    // this is an implementation of shell sort that
    // additionally keeps track of the sorting order
    int flag = 1;
    int d = length;
    int i, idTmp;
    double xTmp;
    
    while (flag || d>1) {
        flag = 0;
        d = (d+1)/2;
        for (i=0; i<(length-d); i++) {
            if (x[i+d] < x[i]) {
                // swap items at positions i+d and d
	      xTmp = x[i+d];  idTmp = id[i+d]; 
	      x[i+d] = x[i];  id[i+d] = id[i]; 
	      x[i] = xTmp;    id[i] = idTmp; 
	      // indicate that a swap has occurred
	      flag = 1;
            }
        }
    }

    return;
}