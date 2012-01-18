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

#ifndef MCFTRCPlaneStress_h
#define MCFTRCPlaneStress_h

// $Revision: 1.2 $
// $Date: 2011/08/31 23:00:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/MCFTRCPlaneStress.h,v $
// File: MCFTRCPlaneStress.h
//
// Written: Lining
// Created: 2011.8
//
// Description: This file contains the class definition for 
// MCFTRCPlaneStress material.
//
// For Detailed explanation of the model, please refer to the journal paper
// entitled "The modified compression-filed theory for reinforced concrete element
// subjected to shear, ACI Journal. s83-22. 1986. pp. 219-231"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <UniaxialMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Tensor.h>

class MCFTRCPlaneStress : public NDMaterial
{
  public:
    MCFTRCPlaneStress ( int      tag, 
				      double   RHO,
				      UniaxialMaterial *s1,
				      UniaxialMaterial *s2,
				      UniaxialMaterial *c1,
				      UniaxialMaterial *c2,
				      double   ANGLE1,
				      double   ANGLE2,
				      double   ROU1,
				      double   ROU2,
					  double   DB1,
					  double   DB2,
				      double   FPC,
				      double   FYX,
					  double   FYY,
				      double   E,
				      double   EPSC0,
					  double   AGGR,
					  double   XD,
					  double   YD);
    MCFTRCPlaneStress();
    ~MCFTRCPlaneStress();				  
    
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

    double   angle1;    // angel of the first steel layer to x coordinate 
    double   angle2;    // angel of the second steel layer to x coordinate
    double   roux;      // steel ratio of the first steel layer
    double   rouy;      // steel ratio of the second steel layer
	double   db1;       // steel diameter of the first steel layer
	double   db2;       // steel diameter of the second steel layer
    double   fpc;       // compressive strength of the concrete
    double   fyx;       // yield stress of the bare steel bar x
	double   fyy;       // yield stress of the bare steel bar y
    double   E0;        // young's modulus of the steel
    double   epsc0;     // compressive strain of the concrete
	double   aggr;      // aggregate size
	double   xd;        // x- crack spacing
    double   yd;        // y- crack spacing

    Vector   Tstress;  // Trial stresses
    Vector   lastStress;  // Last committed stresses, added for x, k
    
    int      steelStatus;  // check if steel yield, 0 not yield, 1 yield
    int      dirStatus;    // check if principle direction has exceed 90 degree, 1 yes, 0 no
    
    // for damgage factor D=1-0.4*epslonC'/epslon0; epslon0=0.002
    
    // Trial values
	Vector strain_vec;
	Vector epsC_vec;
	Vector epsCp_vec;
	Vector epsC12p_prevec; // pre time step eps
	Vector epsC12p_nowvec; // this time step eps
	
	Vector epsCcm_vec;
	Vector epsCtm_vec;
	Vector epsC12cm_vec;
	Vector epsC12tm_vec;
	
    // Converged values
	Vector CepsC_vec;
	Vector CepsCp_vec;

	Vector CepsCcm_vec;
	Vector CepsCtm_vec;
    Vector CepsC12cm_vec;
	Vector CepsC12tm_vec;

    double DDOne; // damage factor for concrete One
    double DDTwo; // damage factor for concrete Two

	Vector stress0_vec;
	Vector stress_vec;
	Matrix tangent_matrix;

    Vector determineTrialStress(Vector strain);
    double getPrincipalStressAngle(double inputAngle);
    double getAngleError(double inputCita);
 
	double kupferenvelop(double Tstrain, double sig_p, double eps_p);
	int determineTangent(void);
};

#endif
