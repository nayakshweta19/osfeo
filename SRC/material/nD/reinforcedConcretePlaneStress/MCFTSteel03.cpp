// MCFTSteel03.cpp
// Written: neallee@tju.edu.cn
// Created: 2011.7
//
// Description: This file contains the class implementation for
// MCFTSteel03


#include "MCFTSteel03.h"
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>

#include <MaterialResponse.h>
#include <elementAPI.h>
#define OPS_Export 

OPS_Export void *
OPS_NewMCFTSteel03Material(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs < 4 || numRemainingArgs > 12) {
    opserr << "Invalid Args want: uniaxialMaterial MCFTSteel03 tag? fy? E? b? "
		   << "<R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
    return 0;	
  }

  int    iData[1];
  double dData[12];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial MCFTSteel03 tag" << endln;
    return 0;
  }

  numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 3) {
    if (OPS_GetDoubleInput(&numRemainingArgs, dData) != 0) {
      opserr << "Invalid Args want: uniaxialMaterial MCFTSteel03 tag? fy? E? b? " << endln;
      return 0;	
    } else
      theMaterial = new MCFTSteel03(iData[0], dData[0], dData[1], dData[2]);
  } else if (numRemainingArgs == 6) {
    if (OPS_GetDoubleInput(&numRemainingArgs, dData) != 0) {
      opserr << "Invalid Args want: uniaxialMaterial MCFTSteel03 tag? fy? E? b? "
		     << " R0 cR1 cR2" << endln;
      return 0;	
    } else
      theMaterial = new MCFTSteel03(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);
  } else if (numRemainingArgs == 10) {
	if (OPS_GetDoubleInput(&numRemainingArgs, dData) != 0) {
      opserr << "Invalid Args want: uniaxialMaterial MCFTSteel03 tag? fy? E? b? "
		     << " R0 cR1 cR2 a1 a2 a3 a4" << endln;
      return 0;	
    } else
      theMaterial = new MCFTSteel03(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
	                                dData[6], dData[7], dData[8], dData[9]);
  } else if (numRemainingArgs == 11) {
	if (OPS_GetDoubleInput(&numRemainingArgs, dData) != 0) {
      opserr << "Invalid Args want: uniaxialMaterial MCFTSteel03 tag? fy? E? b? "
		     << " R0 cR1 cR2 a1 a2 a3 a4 sigini" << endln;
      return 0;	
    } else
      theMaterial = new MCFTSteel03(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
	                                dData[6], dData[7], dData[8], dData[9], dData[10]);
  } else {
	return 0;
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type MCFTSteel03\n";
    return 0;
  }

  return theMaterial;
}

MCFTSteel03::MCFTSteel03(int tag,
	double _Fy, double _E0, double _b,
	double _R0, double _cR1, double _cR2,
	double _a1, double _a2, double _a3, double _a4, double sigInit):
UniaxialMaterial(tag, MAT_TAG_MCFTSteel03),
	Fy(_Fy), E0(_E0), b(_b), R0(_R0), cR1(_cR1), cR2(_cR2), a1(_a1), a2(_a2), a3(_a3), a4(_a4), 
	sigini(sigInit)
{
  konP = 0;
  kon = 0;
  eP = E0;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;
  e = E0;

  epsmaxP = Fy/E0;
  epsminP = -epsmaxP;
  epsplP = 0.0;
  epss0P = 0.0;
  sigs0P = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;

  if (sigini != 0.0) {
    epsP = sigini/E0;
    sigP = sigini;
  }
}

MCFTSteel03::MCFTSteel03(int tag,
	double _Fy, double _E0, double _b,
	double _R0, double _cR1, double _cR2):
UniaxialMaterial(tag, MAT_TAG_MCFTSteel03),
	Fy(_Fy), E0(_E0), b(_b), R0(_R0), cR1(_cR1), cR2(_cR2), sigini(0.0)
{
  konP = 0;

  // Default values for no isotropic hardening
  a1 = 0.0;
  a2 = 1.0;
  a3 = 0.0;
  a4 = 1.0;

  eP = E0;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;
  e = E0;

  epsmaxP = Fy/E0;
  epsminP = -epsmaxP;
  epsplP = 0.0;
  epss0P = 0.0;
  sigs0P = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;
}

MCFTSteel03::MCFTSteel03(int tag,
	double _Fy, double _E0, double _b):
UniaxialMaterial(tag, MAT_TAG_MCFTSteel03),
	Fy(_Fy), E0(_E0), b(_b), sigini(0.0)
{
  konP = 0;

  // Default values for elastic to hardening transitions
  R0 = 15.0;
  cR1 = 0.925;
  cR2 = 0.15;

  // Default values for no isotropic hardening
  a1 = 0.0;
  a2 = 1.0;
  a3 = 0.0;
  a4 = 1.0;

  eP = E0;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;
  e = E0;

  epsmaxP = Fy/E0;
  epsminP = -epsmaxP;
  epsplP = 0.0;
  epss0P = 0.0;
  sigs0P = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;
}

MCFTSteel03::MCFTSteel03(void):
UniaxialMaterial(0, MAT_TAG_MCFTSteel03)
{
  konP = 0;
}

MCFTSteel03::~MCFTSteel03()
{
	// does nothing;
}

UniaxialMaterial*
MCFTSteel03::getCopy(void)
{
  MCFTSteel03 *theCopy = new MCFTSteel03(this->getTag(), Fy, E0, b, R0, cR1, cR2, a1, a2, a3, a4, sigini);
  
  return theCopy;
}

double
MCFTSteel03::getInitialTangent(void)
{
  return E0;
}

int
MCFTSteel03::setTrialStrain(double trialStrain, double strainRate)
{
  double Esh = b * E0;
  double epsy = Fy / E0;

  // modified C-P. Lamarche 2006
  if (sigini != 0.0) {
    double epsini = sigini/E0;
    eps = trialStrain+epsini;
  } else
    eps = trialStrain;
  // modified C-P. Lamarche 2006

  double deps = eps - epsP;
  
  epsmax = epsmaxP;
  epsmin = epsminP;
  epspl  = epsplP;
  epss0  = epss0P;  
  sigs0  = sigs0P; 
  epsr   = epssrP;  
  sigr   = sigsrP;  
  kon = konP;

  if (kon == 0 || kon == 3) { // modified C-P. Lamarche 2006


    if (fabs(deps) < 10.0*DBL_EPSILON) {

      e = E0;
      sig = sigini;                // modified C-P. Lamarche 2006
      kon = 3;                     // modified C-P. Lamarche 2006 flag to impose initial stess/strain
      return 0;

    } else {

      epsmax = epsy;
      epsmin = -epsy;
      if (deps < 0.0) {
	kon = 2;
	epss0 = epsmin;
	sigs0 = -Fy;
	epspl = epsmin;
      } else {
	kon = 1;
	epss0 = epsmax;
	sigs0 = Fy;
	epspl = epsmax;
      }
    }
  }
  
  // in case of load reversal from negative to positive strain increment, 
  // update the minimum previous strain, store the last load reversal 
  // point and calculate the stress and strain (sigs0 and epss0) at the 
  // new intersection between elastic and strain hardening asymptote 
  // To include isotropic strain hardening shift the strain hardening 
  // asymptote by sigsft before calculating the intersection point 
  // Constants a3 and a4 control this stress shift on the tension side 
  
  if (kon == 2 && deps > 0.0) {


    kon = 1;
    epsr = epsP;
    sigr = sigP;
    //epsmin = min(epsP, epsmin);
    if (epsP < epsmin)
      epsmin = epsP;
    double d1 = (epsmax - epsmin) / (2.0*(a4 * epsy));
    double shft = 1.0 + a3 * pow(d1, 0.8);
    epss0 = (Fy * shft - Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
    sigs0 = Fy * shft + Esh * (epss0 - epsy * shft);
    epspl = epsmax;

  } else if (kon == 1 && deps < 0.0) {
    
    // update the maximum previous strain, store the last load reversal 
    // point and calculate the stress and strain (sigs0 and epss0) at the 
    // new intersection between elastic and strain hardening asymptote 
    // To include isotropic strain hardening shift the strain hardening 
    // asymptote by sigsft before calculating the intersection point 
    // Constants a1 and a2 control this stress shift on compression side 

    kon = 2;
    epsr = epsP;
    sigr = sigP;
    //      epsmax = max(epsP, epsmax);
    if (epsP > epsmax)
      epsmax = epsP;
    
    double d1 = (epsmax - epsmin) / (2.0*(a2 * epsy));
    double shft = 1.0 + a1 * pow(d1, 0.8);
    epss0 = (-Fy * shft + Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
    sigs0 = -Fy * shft + Esh * (epss0 + epsy * shft);
    epspl = epsmin;
  }

  
  // calculate current stress sig and tangent modulus E 

  double xi     = fabs((epspl-epss0)/epsy);
  double R      = R0*(1.0 - (cR1*xi)/(cR2+xi));
  double epsrat = (eps-epsr)/(epss0-epsr);
  double dum1  = 1.0 + pow(fabs(epsrat),R);
  double dum2  = pow(dum1,(1/R));

  sig   = b*epsrat +(1.0-b)*epsrat/dum2;
  sig   = sig*(sigs0-sigr)+sigr;

  e = b + (1.0-b)/(dum1*dum2);
  e = e*(sigs0-sigr)/(epss0-epsr);

  return 0;
}

double 
MCFTSteel03::getStrain(void)
{
  return eps;
}

double 
MCFTSteel03::getStress(void)
{
  return sig;
}

double 
MCFTSteel03::getTangent(void)
{
  return e;
}

double
MCFTSteel03::getSecant ()
{
	if ( eps == 0.0 ) {
		return E0;
	} else {
		return sig/eps;
	}
}

int 
MCFTSteel03::commitState(void)
{
  epsminP = epsmin;
  epsmaxP = epsmax;
  epsplP = epspl;
  epss0P = epss0;
  sigs0P = sigs0;
  epssrP = epsr;
  sigsrP = sigr;
  konP = kon;
  
  eP = e;
  sigP = sig;
  epsP = eps;

  return 0;
}

int 
MCFTSteel03::revertToLastCommit(void)
{
  epsmin = epsminP;
  epsmax = epsmaxP;
  epspl = epsplP;
  epss0 = epss0P;
  sigs0 = sigs0P;
  epsr = epssrP;
  sigr = sigsrP;
  kon = konP;
  
  e = eP;
  sig = sigP;
  eps = epsP;
  return 0;
}

int 
MCFTSteel03::revertToStart(void)
{
  eP = E0;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;
  e = E0;  

  konP = 0;
  epsmaxP = Fy/E0;
  epsminP = -epsmaxP;
  epsplP = 0.0;
  epss0P = 0.0;
  sigs0P = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;

  if (sigini != 0.0) {
	  epsP = sigini/E0;
	  sigP = sigini;
   } 

  return 0;
}

int 
MCFTSteel03::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(23);
  data(0) = Fy;
  data(1) = E0;
  data(2) = b;
  data(3) = R0;
  data(4) = cR1;
  data(5) = cR2;
  data(6) = a1;
  data(7) = a2;
  data(8) = a3;
  data(9) = a4;
  data(10) = epsminP;
  data(11) = epsmaxP;
  data(12) = epsplP;
  data(13) = epss0P;
  data(14) = sigs0P;
  data(15) = epssrP;
  data(16) = sigsrP;
  data(17) = konP;  
  data(18) = epsP;  
  data(19) = sigP;  
  data(20) = eP;    
  data(21) = this->getTag();
  data(22) = sigini;

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "MCFTSteel03::sendSelf() - failed to sendSelf\n";
    return -1;
  }
  return 0;
}

int 
MCFTSteel03::recvSelf(int commitTag, Channel &theChannel, 
	     FEM_ObjectBroker &theBroker)
{
  static Vector data(23);

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "MCFTSteel03::recvSelf() - failed to recvSelf\n";
    return -1;
  }

  Fy = data(0);
  E0 = data(1);
  b = data(2); 
  R0 = data(3);
  cR1 = data(4);
  cR2 = data(5);
  a1 = data(6); 
  a2 = data(7); 
  a3 = data(8); 
  a4 = data(9); 
  epsminP = data(10);
  epsmaxP = data(11);
  epsplP = data(12); 
  epss0P = data(13); 
  sigs0P = data(14); 
  epssrP = data(15); 
  sigsrP = data(16); 
  konP = int(data(17));   
  epsP = data(18);   
  sigP = data(19);   
  eP   = data(20);   
  this->setTag(int(data(21)));
  sigini = data(22);

  e = eP;
  sig = sigP;
  eps = epsP;
  
  return 0;
}

Response* 
MCFTSteel03::setResponse(const char **argv, int argc,
			 OPS_Stream &theOutput)
{
  Response *theResponse = 0;

  if (strcmp(argv[0],"setVar") == 0) {
    theResponse = new MaterialResponse(this, 100, Vector(6));
  } else if (strcmp(argv[0],"getVar") == 0) {
    theResponse = new MaterialResponse(this, 101, Vector(6));
  } else
    return this->UniaxialMaterial::setResponse(argv, argc, theOutput);

  return theResponse;
}

int 
MCFTSteel03::getResponse(int responseID, Information &matInfo)
{
  if (responseID == 100) { // setVar
    //matInfo.theDouble = 0.0;
	Vector *theVector = matInfo.theVector;
	theVector->Zero();
	//ecminP = (*theVector)(0); 
	//deptP  = (*theVector)(1); 
	epsP   = (*theVector)(2); 
	sigP   = (*theVector)(3); 
	eP     = (*theVector)(4); 
	//epscp  = (*theVector)(5); 
  } else if (responseID == 101){ // get var
    Vector *theVector = matInfo.theVector;
    //(*theVector)(0) = ecminP;
    //(*theVector)(1) = deptP;
    (*theVector)(2) = epsP;
    (*theVector)(3) = sigP;
    (*theVector)(4) = eP;
	//(*theVector)(5) = epscp;
  } else
    return this->UniaxialMaterial::getResponse(responseID, matInfo);

  return 0;
}

void 
MCFTSteel03::Print(OPS_Stream &s, int flag)
{
  s << "MCFTSteel03:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
}
