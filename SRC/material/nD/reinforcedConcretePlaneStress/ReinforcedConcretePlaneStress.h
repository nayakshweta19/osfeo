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

#ifndef ReinforcedConcretePlaneStress_h
#define ReinforcedConcretePlaneStress_h

// File: ReinforcedConcretePlaneStress.h
//
// Written: JZhong
// Created: 2003.10
//
// Written: Lining
// Created: 2010.11
////
// Description: This file contains the class definition for 
// ReinforcedConcretePlaneStress material.
// Hsu's Model 2002
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

class ReinforcedConcretePlaneStress : public NDMaterial
{
  public:
    ReinforcedConcretePlaneStress (int      tag, 
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
    ReinforcedConcretePlaneStress();
    ~ReinforcedConcretePlaneStress();				  
    
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
    
    int RotatedMatrixAtCita(const double &cita, const Matrix &M);
    
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

	Response *setResponse(const char **argv, int argc, OPS_Stream &s);
	int getResponse(int responseID, Information &matInformation);

  protected:
    
  private:
    double   rho; 
    UniaxialMaterial **theMaterial; // pointer of the materials 
    Response **theResponses; // pointer to material responses needed for Concrete

	double citaR;           // principal strain direction
	double lastCitaR;       // last converged principle strain direction
	//bool   isSwapped;       // primary concrete direction has changed
	//int    lastDirStatus;
	int    steelStatus;     // check if steel yield, 0 not yield, 1 yield
	//int    dirStatus;       // check if principle direction has exceed 90 degree, 1 yes, 0 no
    double beta;            // citaOne - citaR

	double epslonOne;
	double epslonTwo;
	double halfGammaOneTwo;

	double sigmaOneC;
	double sigmaTwoC;

	Vector strain_vec;
	Vector stress_vec;
	Matrix tangent_matrix;
    
    Matrix DC;
    Matrix DC_bar;

    double   angle1;    // angel of the first steel layer to x coordinate 
    double   angle2;    // angel of the second steel layer to x coordinate
    double   rouL;      // steel ratio of the first steel layer
    double   rouT;      // steel ratio of the second steel layer
    double   fpc;       // compressive strength of the concrete
    double   fy;        // yield stress of the bare steel bar
    double   E0;        // young's modulus of the steel
    double   epsc0;     // compressive strain of the concrete
	Vector   Tstrain;     // Trial strains  //epslonx,epslony,0.5*gammaxy	
    Vector   Tstress;     // Trial stresses
    Vector   lastStress;  // Last committed stresses, added for x, k

    double   citaStrain;   // principle strain direction
    double   citaStress;   // principle stress direction
    double   miu12;        // Hsu/Zhu ratio
    double   miu21;        // Hsu/Zhu ratio
    double   G12;
    // for damgage factor D=1-0.4*epslonC'/epslon0; epslon0=0.002
    
    // Trial values
    int    TOneReverseStatus;         // Trial reverse status for concrete One, 1 reverse, 0 no
    double TOneNowMaxComStrain;
    double TOneLastMaxComStrain;
    
    int    TTwoReverseStatus;         // Trial reverse status for concrete Two, 1 reverse, 0 no
    double TTwoNowMaxComStrain;
    double TTwoLastMaxComStrain;
    
    // Converged values
    int    COneReverseStatus;         // Converged reverse status for concrete One, 1 reverse, 0 no
    double COneNowMaxComStrain;
    double COneLastMaxComStrain;
    
    int    CTwoReverseStatus;         // Converged reverse status for concrete Two, 1 reverse, 0 no
    double CTwoNowMaxComStrain;
    double CTwoLastMaxComStrain;
    
    double DDOne; // damage factor for concrete One
    double DDTwo; // damage factor for concrete Two
    
    double tt1, tt2;
    double xxx, kkk;
    
    int determineTrialStress(void);
//	void determineConcreteStatus(int);
};

#endif
