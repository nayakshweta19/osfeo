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
                                                                        
// $Revision: 1.4 $
// $Date: 2003/02/14 23:01:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlaneStressRCFiberMaterial.h,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition of PlaneStressRCFiberMaterial.
// The PlaneStressRCFiberMaterial class is a wrapper class that performs static
// condensation on a three-dimensional material model to give the 11, 12
// stress components which can then be integrated over an area to model a
// shear flexible 2D beam.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Vector.h>
#include <Matrix.h>
#include <ID.h> 
#include <NDMaterial.h>

class PlaneStressRCFiberMaterial: public NDMaterial {

  public:
    PlaneStressRCFiberMaterial(int tag, NDMaterial &theMat);
    PlaneStressRCFiberMaterial(void);
    virtual ~PlaneStressRCFiberMaterial(void);

    int setTrialStrain( const Vector &strainFromElement);
    const Vector& getStrain(void);
    const Vector& getStress(void);
    const Matrix& getTangent(void);
    const Matrix& getInitialTangent(void);

    double getRho(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    NDMaterial *getCopy(void);
    NDMaterial *getCopy(const char *type);
    const char *getType(void) const;
    int getOrder(void) const; 

    void Print(OPS_Stream &s, int flag);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  private:
    double Tstrain22;
    double Cstrain22;
	Matrix twoDtgLastCommit;

    NDMaterial *theMaterial;

    Vector strain;

    static Vector stress;
    static Matrix tangent;

	int    maxIters;               // maximum number of local iterations
	double tol;	                   // tolerance for relative energy norm for local

    //int indexMap(int i);
	static int iMap[3];
};





