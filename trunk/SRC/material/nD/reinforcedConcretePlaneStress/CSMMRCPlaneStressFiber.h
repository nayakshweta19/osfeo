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

#ifndef CSMMRCPlaneStressFiber_h
#define CSMMRCPlaneStressFiber_h

// $Revision: 1.2 $
// $Date: 2003/02/14 23:00:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/CSMMRCPlaneStressFiber.h,v $
// File: CSMMRCPlaneStressFiber.h
//
// Written: Li Ning (neallee@tju.edu.cn)
// Created: 2010.11
//
// Description: This file contains the class definition for 
// CSMMRCPlaneStressFiber material.

// For Detailed explanation of the model, please refer to the book
// entitled "Unified Theory of Concrete Structures,"
// by Thomas T.C. Hsu and Y.L. Mo, John Wiley & Sons, April 2010.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <UniaxialMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Tensor.h>

class CSMMRCPlaneStressFiber : public NDMaterial
{
  public:
    CSMMRCPlaneStressFiber ( int      tag, 
				      double   RHO,
				      UniaxialMaterial *s1,
				      UniaxialMaterial *s2,
				      UniaxialMaterial *c1,
				      UniaxialMaterial *c2,
				      double   ANGLE1,
				      double   ANGLE2,
				      double   ROU1,
				      double   ROU2,
				      double   FPC,
				      double   FY,
				      double   E,
				      double   EPSC0);					  
    CSMMRCPlaneStressFiber();
    ~CSMMRCPlaneStressFiber();				  
    
    double getRho(void);
    
    int setTrialStrain(const Vector &v); // really used one
    int setTrialStrain(const Vector &v, const Vector &r);
    int setTrialStrainIncr(const Vector &v);
    int setTrialStrainIncr(const Vector &v, const Vector &r);
    const Matrix &getTangent(void);
    const Matrix &getInitialTangent(void) {return this->getTangent();};

    const Vector &getStress(void);
    const Vector &getStrain(void);
    
    const Vector &getCommittedStress(void);
    const Vector &getCommittedStrain(void);    
    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    
    NDMaterial *getCopy(void);
    NDMaterial *getCopy(const char *type);

    void Print(OPS_Stream &s, int flag = 0);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    const char *getType(void) const { return "PlaneStress"; };
    int getOrder(void) const { return 3;};

	Response *setResponse (const char **argv, int argc, OPS_Stream &s);
	int getResponse (int responseID, Information &matInformation);

  protected:
    
  private:

    double   rho; 
    UniaxialMaterial **theMaterial; // pointer of the materials 
    Response **theResponses; // pointer to material responses needed for Concrete

	static double citaR;           // principal strain direction
	static double lastCitaR;       // last converged principle strain direction
	static bool   isSwapped;       // primary concrete direction has changed
	static int    lastDirStatus;
	static int    steelStatus;     // check if steel yield, 0 not yield, 1 yield
	static int    dirStatus;       // check if principle direction has exceed 90 degree, 1 yes, 0 no

	static double epslonOne;
	static double epslonTwo;
	static double halfGammaOneTwo;

	static double sigmaOneC;
	static double sigmaTwoC;

	static Vector strain_vec;
	static Vector stress_vec;
	static Matrix tangent_matrix;

	double   angle1;    // angel of the first steel layer to x coordinate 
    double   angle2;    // angel of the second steel layer to x coordinate
    double   rou1;      // steel ratio of the first steel layer
    double   rou2;      // steel ratio of the second steel layer
    double   fpc;       // compressive strength of the concrete
    double   fy;        // yield stress of the bare steel bar
    double   E0;        // young's modulus of the steel
    double   epsc0;     // compressive strain of the concrete
	Vector   Tstrain;   // Trial strains
    Vector   Tstress;   // Trial stresses
    Vector   lastStress;  // Last committed stresses, added for x, k
    
    double   citaStrain;      // principle strain direction
    double   citaStress;     // principle stress direction
    double   miu12;        // Hsu/Zhu ratio
    double   miu21;        // Hsu/Zhu ratio
    double   G12;          // Shear Modulus
    
    
    // for damgage factor D=1-0.4*epslonC'/epslon0; epslon0=0.002
    
    // Trial values
    int TOneReverseStatus;         // Trial reverse status for concrete One, 1 reverse, 0 no
    double TOneNowMaxComStrain;
    double TOneLastMaxComStrain;
    
    int TTwoReverseStatus;         // Trial reverse status for concrete Two, 1 reverse, 0 no
    double TTwoNowMaxComStrain;
    double TTwoLastMaxComStrain;
    
    // Converged values
    int COneReverseStatus;         // Converged reverse status for concrete One, 1 reverse, 0 no
    double COneNowMaxComStrain;
    double COneLastMaxComStrain;
    
    int CTwoReverseStatus;         // Converged reverse status for concrete Two, 1 reverse, 0 no
    double CTwoNowMaxComStrain;
    double CTwoLastMaxComStrain;
    
    double DDOne; // damage factor for concrete One
    double DDTwo; // damage factor for concrete Two

	Vector fiberStrain;
	Vector fiberStress;
	Matrix fiberTangent;

    int determineTrialStress(void);
    double getPrincipalStressAngle(double inputAngle);
    double getAngleError(double inputCita);
	void   determineConcreteStatus(int);
    
};

#endif
