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
// $Date: 2007/06/08 00:38:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nd/reinforcedConcretePlaneStress/MCFTConcrete03.cpp,v $

// Written: neal(neallee@tju.edu.cn)
// Created: 07/11
//
// Description: This file contains the class definition for 
// MCFTConcrete03. MCFTConcrete03 is based on an f2c of the FEDEAS material


#include "MCFTConcrete03.h"
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include <MaterialResponse.h>
#include <elementAPI.h>

using namespace std;

#define OPS_Export 

OPS_Export void *
OPS_NewMCFTConcrete03Material(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[7];
  int numData = 1;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs < 8) {
    opserr << "Want: uniaxialMaterial MCFTConcrete03 tag? fpc? epsc0? fpcu? epscu? rat? ft? Ets?" << endln;
    return 0;	
  }

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial MCFTConcrete03 tag?  fpc? epsc0? fpcu? epscu? rat? ft? Ets?" << endln;
    return 0;
  }

  numData = 7;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial MCFTConcrete03 tag?  fpc? epsc0? fpcu? epscu? rat? ft? Ets?" << endln;
    return 0;	
  }

  // Parsing was successful, allocate the material
  theMaterial = new MCFTConcrete03(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type MCFTConcrete03\n";
    return 0;
  }

  return theMaterial;
}

MCFTConcrete03::MCFTConcrete03(int tag, double _fc, double _epsc0, double _fcu,
		       double _epscu, double _rat, double _ft, double _Ets):
  UniaxialMaterial(tag, MAT_TAG_MCFTConcrete03),
  fc(_fc), epsc0(_epsc0), fcu(_fcu), epscu(_epscu), rat(_rat), ft(_ft), Ets(_Ets)
{
  ecminP = 0.0;
  ecmaxP = 0.0;
  sigminP = 0.0;
  sigmaxP = 0.0;

  deptP = 0.0;
  eP = 2.0*fc/epsc0;
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
  e = 2.0*fc/epsc0;

  epscp = 0.0;
  betaD = 1.0;
  K = 1.0;

  ept = 0.0;
  epsro = 0.0;
  sigro =0.0;
}

MCFTConcrete03::MCFTConcrete03(void):
  UniaxialMaterial(0, MAT_TAG_MCFTConcrete03),
  fc(0), epsc0(0), fcu(0), epscu(0), rat(0), ft(0), Ets(0)
{

}

MCFTConcrete03::~MCFTConcrete03(void)
{
  // Does nothing
}

UniaxialMaterial*
MCFTConcrete03::getCopy(void)
{
  MCFTConcrete03 *theCopy = new MCFTConcrete03(this->getTag(), fc, epsc0, fcu, epscu, rat, ft, Ets);
  
  return theCopy;
}

double
MCFTConcrete03::getInitialTangent(void)
{
  return 2.0*fc/epsc0;
}

int
MCFTConcrete03::setTrialStrain(double trialStrain, double strainRate)
{
  double  ec0 = fc * 2. / epsc0;

  // retrieve concrete history variables

  ecmin = ecminP;
  dept = deptP;

  // calculate current strain

  eps = trialStrain;
  double deps = eps - epsP;

  if (fabs(deps) < DBL_EPSILON)
    return 0;

  // if the current strain is less than the smallest previous strain 
  // call the monotonic envelope in compression and reset minimum strain 

  if (eps < ecmin) {
    this->Compr_Envlp(eps, sig, e);
    ecmin = eps;
	this->determineCompEpscp(ecmin);

  } else {

    // else, if the current strain is between the minimum strain and ept 
    // (which corresponds to zero stress) the material is in the unloading- 
    // reloading branch and the stress remains between sigmin and sigmax 
    
    // calculate strain-stress coordinates of point R that determines 
    // the reloading slope according to Fig.2.11 in EERC Report 
    // (corresponding equations are 2.31 and 2.32 
    // the strain of point R is epsR and the stress is sigmR 
    
    //double epsr = (fcu - rat * ec0 * epscu) / (ec0 * (1.0 - rat));
    //double sigmr = ec0 * epsr;

    // calculate the previous minimum stress sigmm from the minimum 
    // previous strain ecmin and the monotonic envelope in compression 
    
    double sigmm;
    double dumy;
    this->Compr_Envlp(ecmin, sigmm, dumy);
    
    // calculate current reloading slope Er (Eq. 2.35 in EERC Report) 
    // calculate the intersection of the current reloading slope Er 
    // with the zero stress axis (variable ept) (Eq. 2.36 in EERC Report) 
    
    //double er = (sigmm - sigmr) / (ecmin - epsr);
    //ept = ecmin - sigmm / er;

	/*this->determineCompEpscp(ecmin);*/
	ept = epscp;

    if (eps <= ept) {
	  
	  // if the current strain is between the minimum strain and ept 
	  // (which corresponds to zero stress) the material is in the unloading- 
	  // reloading branch and the stress remains between sigmin and sigmax 

    //  double sigmin = sigmm + er * (eps - ecmin);
    //  double sigmax = er * .5f * (eps - ept);
    //  sig = sigP + ec0 * deps;
    //  e = ec0;
    //  if (sig <= sigmin) {
	//sig = sigmin;
	//e = er;
    //  }
    //  if (sig >= sigmax) {
	//sig = sigmax;
	//e = 0.5 * er;
    //  }
	  if (deps >= 0) { //unloading
    
	double slop2 = ec0*K*betaD;
	double slop3 = 0.071 * ec0*K*betaD;
	double N = (slop2-slop3)*(epscp-ecmin)/(sigmm+slop2*(epscp-ecmin));
	//N = (N < 0.0 ? 0 : N);
	double dStrain = eps - ecmin;
	if (N>0.0 && pow(epscp-ecmin,N-1.0) != 0)
	  sig = sigmm+slop2*dStrain+(slop3-slop2)*pow(dStrain,N)/N/pow(epscp-ecmin,N-1.0);
	else 
	  sig = sigmm+slop2*dStrain;
	e = slop2 + (slop3-slop2)*pow(dStrain/(epscp-ecmin), N-1.0);
	sigro = sig;
	epsro = eps;

      } else { //reloading

    double beta;  
	if (abs(eps) < abs(epsc0)) beta = 1./(1.+0.1*pow((ecmin-epsro)/epsc0,0.5));
    else                       beta = 1./(1.+0.175*pow((ecmin-epsro)/epsc0,0.6));

	e = (beta*sigmm-sigro)/(ecmin-epsro);
	sig = sigro+e*(eps-epsro);

	  }

	} else {  // eps > ept   eps > epscp
      
      // else, if the current strain is between ept and epn 
      // (which corresponds to maximum remaining tensile strength) 
      // the response corresponds to the reloading branch in tension 
      // Since it is not saved, calculate the maximum remaining tensile 
      // strength sicn (Eq. 2.43 in EERC Report) 
      
      // calculate first the strain at the peak of the tensile stress-strain 
      // relation epn (Eq. 2.42 in EERC Report) 
      
      double epn = ept + dept;
      double sicn;
      if (eps <= epn) {
	this->Tens_Envlp(dept, sicn, e);
	/*this->determineTensEpscp(eps);*/

	if (dept != 0.0) {
	  e = sicn / dept;
	} else {
	  e = ec0;
	}
	sig = e * (eps - ept);
      } else {
	
	// else, if the current strain is larger than epn the response 
	// corresponds to the tensile envelope curve shifted by ept 
	
	double epstmp = eps - ept;
	this->Tens_Envlp(epstmp, sig, e);
	dept = eps - ept;
	/*this->determineTensEpscp(eps);*/
      }
    }
  }

  //opserr << this->getTag() << "Plastic offset strain: " << epscp << "\t" << ept << endln;
  if (e> 0.0 && e <DBL_EPSILON) 
	e = 1.e-10;
  
  return 0;
}

void
MCFTConcrete03::determineCompEpscp(double eps)
{
  //epscp = eps - epsc0*K*betaD*(0.868*(eps/epsc0/K/betaD)-0.166*pow(eps/epsc0/K/betaD, 2.0));
  epscp = eps - epsc0*(0.868*(eps/epsc0)-0.166*pow(eps/epsc0, 2.0));
}

void
MCFTConcrete03::determineTensEpscp(double eps)
{
  epscp = 146.0 * pow(eps, 2.0) + 0.523 * eps;
}

double 
MCFTConcrete03::getStrain(void)
{
  return eps;
}

double 
MCFTConcrete03::getStress(void)
{
  return sig;
}

double 
MCFTConcrete03::getTangent(void)
{
  return e;
}

/*double
MCFTConcrete03::getSecant(void)
{
  if ( abs(eps) <= DBL_EPSILON ) {
    return e;
  } else {
    return sig/(eps-ept);
  }
}*/

int 
MCFTConcrete03::commitState(void)
{
  ecminP = ecmin;
  ecmaxP = ecmax;
  sigminP = sigmin;
  sigmaxP = sigmax;
  
  deptP = dept;
  eptP = ept;

  eP = e;
  sigP = sig;
  epsP = eps;

  return 0;
}

int 
MCFTConcrete03::revertToLastCommit(void)
{
  ecmin = ecminP;
  ecmax = ecmaxP;
  sigmin = sigminP;
  sigmax = sigmaxP;
  
  dept = deptP;
  ept = eptP;

  e = eP;
  sig = sigP;
  eps = epsP;

  return 0;
}

int 
MCFTConcrete03::revertToStart(void)
{
  ecminP = 0.0;
  ecmaxP = 0.0;
  sigminP = 0.0;
  sigmaxP = 0.0;
  ecmin = 0.0;
  ecmax = 0.0;
  sigmin = 0.0;
  sigmax = 0.0;
  dept = 0.0;
  ept = 0.0;

  betaD = 1.0;
  K = 1.0;

  eP = 2.0*fc/epsc0;
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
  e = 2.0*fc/epsc0;

  return 0;
}

int 
MCFTConcrete03::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(20);
  data(0)  = fc;    
  data(1)  = epsc0; 
  data(2)  = fcu;   
  data(3)  = epscu; 
  data(4)  = rat;   
  data(5)  = ft;    
  data(6)  = Ets;   
  data(7)  = ecminP;
  data(8)  = ecmaxP;
  data(9)  = sigminP;
  data(10) = sigmaxP;
  data(11) = eptP;
  data(12) = deptP; 
  data(13) = epsP;  
  data(14) = sigP; 
  data(15) = eP;
  data(16) = epscp;
  data(17) = betaD;
  data(18) = K;
  data(19) = this->getTag();

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "MCFTConcrete03::sendSelf() - failed to sendSelf\n";
    return -1;
  }

  return 0;
}

int 
MCFTConcrete03::recvSelf(int commitTag, Channel &theChannel, 
	     FEM_ObjectBroker &theBroker)
{

  static Vector data(20);

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "MCFTConcrete03::recvSelf() - failed to recvSelf\n";
    return -1;
  }

  fc      = data(0);
  epsc0   = data(1);
  fcu     = data(2);
  epscu   = data(3);
  rat     = data(4);
  ft      = data(5);
  Ets     = data(6);
  ecminP  = data(7);
  ecmaxP  = data(8);
  sigminP = data(9);
  sigmaxP = data(10);
  eptP    = data(11);
  deptP   = data(12);
  epsP    = data(13);
  sigP    = data(14);
  eP      = data(15);
  epscp   = data(16);
  betaD   = data(17);
  K       = data(18);
  this->setTag(int(data(19)));

  e = eP;
  sig = sigP;
  eps = epsP;

  return 0;
}

void 
MCFTConcrete03::Print(OPS_Stream &s, int flag)
{
  s << "\t\t MCFTConcrete03:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
  s << "\t\t  epscp:  " << epscp << endln;
}

void
MCFTConcrete03::Tens_Envlp (double epsc, double &sigc, double &Ect)
{
/*-----------------------------------------------------------------------
! monotonic envelope of concrete in tension (positive envelope)
!
!   ft    = concrete tensile strength
!   Ec0   = initial tangent modulus of concrete 
!   Ets   = tension softening modulus
!   eps   = strain
!
!   returned variables
!    sigc  = stress corresponding to eps
!    Ect  = tangent concrete modulus
!-----------------------------------------------------------------------*/
  
  double Ec0  = 2.0*fc/epsc0;
  double eps0 = ft/Ec0;
  double epsu = ft*(1.0/Ets+1.0/Ec0);

  if (epsc<=eps0) {
    sigc = epsc*Ec0;
    Ect  = Ec0;
  } else {
    if (epsc<=epsu) {
      Ect  = -Ets;
      sigc = ft-Ets*(epsc-eps0);
    } else {
      Ect  = 1.0e-10; //   Ect  = 0.0
      sigc = 0.0;
    }
  }
  return;
}
  
void
MCFTConcrete03::Compr_Envlp (double epsc, double &sigc, double &Ect) 
{
/*-----------------------------------------------------------------------
! monotonic envelope of concrete in compression (negative envelope)
!
!   fc    = concrete compressive strength
!   epsc0 = strain at concrete compressive strength
!   fcu   = stress at ultimate (crushing) strain 
!   epscu = ultimate (crushing) strain
!   Ec0   = initial concrete tangent modulus
!   epsc  = strain
!
!   returned variables
!   sigc  = current stress
!   Ect   = tangent concrete modulus
-----------------------------------------------------------------------*/

  double Ec0  = 2.0*fc/epsc0;
  double ratLocal = epsc/epsc0/K/betaD;
  
  if (epsc>=epsc0*K*betaD) {
    sigc = fc*K*betaD*ratLocal*(2.0-ratLocal);
    Ect  = Ec0*(1.0-ratLocal);
  } else {
    //   linear descending branch between epsc0 and epscu
    if (epsc>epscu) {
      sigc = (fcu-fc*K*betaD)*(epsc-epsc0*K*betaD)/(epscu-epsc0*K*betaD)+fc*K*betaD;
      Ect  = (fcu-fc*K*betaD)/(epscu-epsc0*K*betaD);
    } else {
      // flat friction branch for strains larger than epscu
      sigc = fcu;
      Ect  = 1.0e-10; //Ect  = 0.0
    }
  }
  return;
}

Response* 
MCFTConcrete03::setResponse(const char **argv, int argc,
			 OPS_Stream &theOutput)
{
  Response *theResponse = 0;

  if (strcmp(argv[0],"setVar") == 0) {
    theResponse = new MaterialResponse(this, 100, Vector(8));
  } else if (strcmp(argv[0],"getVar") == 0) {
    theResponse = new MaterialResponse(this, 101, Vector(8));
  } else
    return this->UniaxialMaterial::setResponse(argv, argc, theOutput);

  return theResponse;
}
 
int 
MCFTConcrete03::getResponse(int responseID, Information &matInfo)
{
  if (responseID == 100) { // setVar
    //matInfo.theDouble = 0.0;
	Vector *theVector = matInfo.theVector;
	ecminP  = (*theVector)(0);
	ecmaxP  = (*theVector)(1);
	sigminP = (*theVector)(2);
	sigmaxP = (*theVector)(3);
	eptP    = (*theVector)(4); 
	epscp   = (*theVector)(5); 
	betaD   = (*theVector)(6); 
	K       = (*theVector)(7); 
  } else if (responseID == 101){ // get var
    Vector *theVector = matInfo.theVector;
	(*theVector)(0) = ecminP;
	(*theVector)(1) = ecmaxP;
	(*theVector)(2) = sigminP;
	(*theVector)(3) = sigmaxP;
	(*theVector)(4) = eptP;
	//this->determineCompEpscp(eps);
	(*theVector)(5) = epscp;
	(*theVector)(6) = betaD;
	(*theVector)(7) = K;
  } else
    return this->UniaxialMaterial::getResponse(responseID, matInfo);

  return 0;
}