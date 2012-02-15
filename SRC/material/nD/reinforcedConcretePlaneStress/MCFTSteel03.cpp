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
// $Date: 2006/05/24 21:44:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel03.cpp,v $
// MCFTSteel03.cpp

// Written: neallee@tju.edu.cn
// Created: 2011.7

// Description: This file contains the class implementation for
// MCFTSteel03, Based on Steel03


#include "MCFTSteel03.h"
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
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
  if (numRemainingArgs < 4 || numRemainingArgs > 11) {
    opserr << "Invalid Args want: uniaxialMaterial MCFTSteel03 tag? fy? E0? b? "
		   << "<r0£¿ coeffR1? coeffR2? a1? a2? a3? a4?>" << endln;
    return 0;	
  }

  int    iData[1];
  double dData[10];
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
    } else {
	  double r, coeffR1, coeffR2;
      r=18.0;
      coeffR1 = 0.925;
      coeffR2 =.15;
      theMaterial = new MCFTSteel03(iData[0], dData[0], dData[1], dData[2], r, coeffR1, coeffR2);
	}
  } else if (numRemainingArgs == 10) {
    if (OPS_GetDoubleInput(&numRemainingArgs, dData) != 0) {
      opserr << "Invalid Args want: uniaxialMaterial MCFTSteel03 tag? fy? E? b? "
		     << "r? coeffR1? coeffR2? a1? a2? a3? a4?" << endln;
      return 0;	
    } else
      theMaterial = new MCFTSteel03(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7],
	                                dData[8], dData[9]);
  } else {
    opserr << "Invalid Args want: uniaxialMaterial MCFTSteel03 tag? fy? E? b? "
		   << "r? coeffR1? coeffR2? a1? a2? a3? a4?" << endln;
	return 0;
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type MCFTSteel03\n";
    return 0;
  }

  return theMaterial;
}

MCFTSteel03::MCFTSteel03(
  int tag, double FY, double E, double B, double R,  double CR1, double CR2,
  double A1, double A2, double A3, double A4):
UniaxialMaterial(tag, MAT_TAG_MCFTSteel03),
	fy(FY), E0(E), b(B), r0(R),  coeffR1(CR1), coeffR2(CR2), a1(A1), a2(A2), a3(A3), a4(A4)
{
   // Sets all history and state variables to initial values
   // History variables
   CminStrain = 0;
   CmaxStrain = 0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;
   CYieldStrain = 0.0;
   CYieldStress = 0.0;
   CReverStrain = 0.0;
   CReverStress = 0.0;
   CPlasticExcursion = 0.0;

   TminStrain = 0;
   TmaxStrain = 0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;
   TYieldStrain = 0.0;
   TYieldStress = 0.0;
   TReverStrain = 0.0;
   TReverStress = 0.0;
   TPlasticExcursion = 0.0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = E0;
   CcurR = getR(0);

   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = E0;
   TcurR = getR(0);
// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////
}

MCFTSteel03::MCFTSteel03(void):
UniaxialMaterial(0, MAT_TAG_MCFTSteel03),
 fy(0.0), E0(0.0), b(0.0), r0(0.0), coeffR1(0.0), coeffR2(0.0), a1(0.0), a2(0.0), a3(0.0), a4(0.0)
{

}

MCFTSteel03::~MCFTSteel03()
{
  if (SHVs != 0)  delete SHVs;
}

UniaxialMaterial*
MCFTSteel03::getCopy(void)
{
   MCFTSteel03* theCopy = new MCFTSteel03(this->getTag(), fy, E0, b, r0, coeffR1, coeffR2, a1, a2);

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CmaxStrain = CmaxStrain;
   theCopy->Cloading = Cloading;
   theCopy->CYieldStrain = CYieldStrain;
   theCopy->CYieldStress = CYieldStress;
   theCopy->CReverStrain = CReverStrain;
   theCopy->CReverStress = CReverStress;
   theCopy->CPlasticExcursion = CPlasticExcursion;

   // Trial history variables
   theCopy->TminStrain = TminStrain;
   theCopy->TmaxStrain = TmaxStrain;
   theCopy->Tloading = Tloading;
   theCopy->TYieldStrain = TYieldStrain;
   theCopy->TYieldStress = TYieldStress;
   theCopy->TReverStrain = TReverStrain;
   theCopy->TReverStress = TReverStress;
   theCopy->TPlasticExcursion = TPlasticExcursion;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   // Trial state variables
   theCopy->Tstrain = Tstrain;
   theCopy->Tstress = Tstress;
   theCopy->Ttangent = Ttangent;

   return theCopy;
}

int
MCFTSteel03::setTrialStrain(double strain, double strainRate)
{
   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   TYieldStrain = CYieldStrain;
   TYieldStress = CYieldStress;
   TReverStrain = CReverStrain;
   TReverStress = CReverStress;
   TPlasticExcursion = CPlasticExcursion;
   TcurR = CcurR;

   // Set trial strain
   Tstrain = strain;

   // Determine change in strain from last converged state
   double dStrain = Tstrain - Cstrain;

   // Calculate the trial state given the trial strain
   if (fabs(dStrain) > DBL_EPSILON)   
     determineTrialState (dStrain);
   
   return 0;
}

int
MCFTSteel03::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   TYieldStrain = CYieldStrain;
   TYieldStress = CYieldStress;
   TReverStrain = CReverStrain;
   TReverStress = CReverStress;
   TPlasticExcursion = CPlasticExcursion;
   TcurR = CcurR;

   // Set trial strain
   Tstrain = strain;

   // Determine change in strain from last converged state
   double dStrain;
   dStrain = Tstrain - Cstrain;

   // Calculate the trial state given the trial strain
   if (fabs(dStrain) > DBL_EPSILON) 
     determineTrialState (dStrain);

   stress = Tstress;
   tangent = Ttangent;

   return 0;
}

double 
MCFTSteel03::getR (double x_in)
{
    // maybe modify this later, but it gives us better degradation at smaller strains
    x_in = fabs(x_in);
    double temp_r = r0;
    
    // new input parameters are supposed to match Steel02 which look like 
    // Dhakal and Maekawa values: cr1 = 0.925, cr2 = 0.15
    // where 0.925 comes from 18.5/20.0 where R0=20 and 18.5 is the old coefficient provided
    //
    // so for old Dhakal and Maekawa R0 = 20, cr1 = 18.5/R0 = 0.925, cr2 = 0.15
    // so for old Gomes and Appleton R0 = 20, cr1 = 19.0/R0 = 0.95, cR2 = 0.3
    // so for my old model, now just use R0 = 20, cR1 = 0, cR2 = 0
    if (coeffR1 < 0.1 && coeffR2 < 0.1) {
        // Mackie, rough trilinear fit to the tangent to the x_in-r first 
        // quadrant circle.  Try using with values of R0 like 20 to 30
        temp_r = r0*2.0/20.0;
        double t1 = -x_in/7+15/7*temp_r;
        double t2 = -4*x_in+6*temp_r;
        if (t1 > temp_r)
            temp_r = t1;
        if (t2 > temp_r)
            temp_r = t2;
        //opserr << "xin = " << x_in << " rout = " << temp_r << endln;
    } else {
    	temp_r = r0 * (1.0 - coeffR1*x_in/(coeffR2+x_in));
        if (temp_r < 0)
            temp_r = 1.0e-8;
    }
    
    return temp_r;
}

void 
MCFTSteel03::determineTrialState (double dStrain)
{
      double fyOneMinusB = fy * (1.0 - b);
      double Esh = b*E0;
	  double epsy = fy/E0;

      double c1 = Esh*Tstrain;
      double c2 = TshiftN*fyOneMinusB;
      double c3 = TshiftP*fyOneMinusB;
      double c = Cstress + E0*dStrain;

      // Determine initial loading condition
      if (Tloading == 0 && dStrain != 0.0) {
		TmaxStrain = epsy;
        TminStrain = -epsy;
        if (dStrain > 0.0) {
	  Tloading = 1;
	  TYieldStrain = TmaxStrain;
	  TYieldStress = fy;
	  TPlasticExcursion = TmaxStrain;
		}
		else {
	  Tloading = -1;
	  TYieldStrain = TminStrain;
	  TYieldStress = -fy;
	  TPlasticExcursion = TminStrain;
		}

        double intval = 1+pow(fabs(Tstrain/epsy),TcurR);
        Tstress = c1+(1-b)*E0*Tstrain/pow(intval,1/TcurR);
        Ttangent = Esh+E0*(1-b)/pow(intval,1+1/TcurR);
	 }

      // Transition from loading to unloading, i.e. positive strain increment
      // to negative strain increment
      if (Tloading == 1 && dStrain < 0.0) {
	    Tloading = -1;
	    if (Cstrain > TmaxStrain)
          TmaxStrain = Cstrain;

		TPlasticExcursion=TminStrain;
	    TshiftN = 1 + a1*pow((TmaxStrain-TminStrain)/(2.0*a2*epsy),0.8);
	    TReverStrain = Cstrain;
	    TReverStress = Cstress;
		TYieldStrain = (c2+c)/E0/(b-1)+Tstrain/(1-b);
		TYieldStress = 1/(b-1)*(b*c2+b*c-c1)-c2;
        TcurR = getR((TYieldStrain-TminStrain)/epsy);
      }

      // Transition from unloading to loading, i.e. negative strain increment
      // to positive strain increment
      if (Tloading == -1 && dStrain > 0.0) {
	    Tloading = 1;
        if (Cstrain < CminStrain)
		  TminStrain = Cstrain;

		TPlasticExcursion=TmaxStrain;
        TshiftP = 1 + a3*pow((TmaxStrain-TminStrain)/(2.0*a4*epsy),0.8);
		TReverStrain = Cstrain;
		TReverStress = Cstress;
		TYieldStrain = (c3-c)/E0/(1-b)+Tstrain/(1-b);
		TYieldStress = 1/(1-b)*(b*c3-b*c+c1)+c3;
	  }
      
      if (Cloading != 0) {
        double c4 = TYieldStrain - TReverStrain;
        double c5 = TYieldStress - TReverStress;
        double c6 = Tstrain - TReverStrain;
        double c4c5 = c5/c4;
        double intval = 1+pow(fabs(c6/c4),TcurR);
        
        Tstress = TReverStress+b*c4c5*c6+(1-b)*c4c5*c6/pow(intval,1/TcurR);
        Ttangent = c4c5*b+c4c5*(1-b)/pow(intval,1+1/TcurR);
      }
}

double 
MCFTSteel03::getStrain(void)
{
  return Tstrain;
}

double 
MCFTSteel03::getStress(void)
{
  return Tstress;
}

double 
MCFTSteel03::getTangent(void)
{
  return Ttangent;
}

double
MCFTSteel03::getSecant(void)
{
	if ( abs(Tstrain) < DBL_EPSILON ) {
		return E0;
	} else {
		return Tstress/Tstrain;
	}
}

int 
MCFTSteel03::commitState(void)
{
   // History variables
   CminStrain = TminStrain;
   CmaxStrain = TmaxStrain;
   CshiftP = TshiftP;
   CshiftN = TshiftN;
   Cloading = Tloading;
   CYieldStrain = TYieldStrain;
   CYieldStress = TYieldStress;
   CReverStrain = TReverStrain;
   CReverStress = TReverStress;
   CPlasticExcursion = TPlasticExcursion;
   CcurR = TcurR;

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int 
MCFTSteel03::revertToLastCommit(void)
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   TYieldStrain = CYieldStrain;
   TYieldStress = CYieldStress;
   TReverStrain = CReverStrain;
   TReverStress = CReverStress;
   TPlasticExcursion = CPlasticExcursion;

   // Reset trial state variables to last committed state
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;
   TcurR = CcurR;

   return 0;
}

int 
MCFTSteel03::revertToStart(void)
{
   // History variables
   CminStrain = 0;
   CmaxStrain = 0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;
   CYieldStrain = 0.0;
   CYieldStress = 0.0;
   CReverStrain = 0.0;
   CReverStress = 0.0;
   CPlasticExcursion = 0.0;

   TminStrain = 0;
   TmaxStrain = 0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;
   TYieldStrain = 0.0;
   TYieldStress = 0.0;
   TReverStrain = 0.0;
   TReverStress = 0.0;
   TPlasticExcursion = 0.0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = E0;
   CcurR = getR(0);

   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = E0;
   TcurR = getR(0);

   //////////////////// sensitivity ///////////////////
    parameterID=0;
	if (SHVs !=0) (*SHVs).Zero();


   return 0;
}

int 
MCFTSteel03::sendSelf(int commitTag, Channel &theChannel)
{
   int res = 0;
   static Vector data(25);
   data(0) = this->getTag();

   // Material properties
   data(1) = fy;
   data(2) = E0;
   data(3) = b;
   data(4) = r0;
   data(5) = coeffR1;
   data(6) = coeffR2;
   data(7) = a1;
   data(8) = a2;
   data(9) = a3;
   data(10) = a4;

   // History variables from last converged state
   data(11) = CminStrain;
   data(12) = CmaxStrain;
   data(13) = CshiftP;
   data(14) = CshiftN;
   data(15) = Cloading;
   data(16) = CYieldStrain;
   data(17) = CYieldStress;
   data(18) = CReverStrain;
   data(19) = CReverStress;
   data(20) = CPlasticExcursion;

   // State variables from last converged state
   data(21) = Cstrain;
   data(22) = Cstress;
   data(23) = Ctangent;
   data(24) = CcurR;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "MCFTSteel03::sendSelf() - failed to send data\n";

   return res;
}

int 
MCFTSteel03::recvSelf(int commitTag, Channel &theChannel, 
	     FEM_ObjectBroker &theBroker)
{
   int res = 0;
   static Vector data(25);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
   if (res < 0) {
      opserr << "MCFTSteel03::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties
      fy = data(1);
      E0 = data(2);
      b = data(3);
      r0 = data(4);
	  coeffR1=data(5);
	  coeffR2=data(6);
      a1 = data(7);
      a2 = data(8);
      a3 = data(9);
      a4 = data(10);

      // History variables from last converged state
      CminStrain = data(11);
      CmaxStrain = data(12);
      CshiftP = data(13);
      CshiftN = data(14);
      Cloading = int(data(15));
      CYieldStrain = data(16);
      CYieldStress = data(17);
      CReverStrain = data(18);
      CReverStress = data(19);
      CPlasticExcursion = data(20);

      // Copy converged history values into trial values since data is only
      // sent (received) after convergence
      TminStrain = CminStrain;
      TmaxStrain = CmaxStrain;
      TshiftP = CshiftP;
      TshiftN = CshiftN;
      Tloading = Cloading;
      TYieldStrain = CYieldStrain;
      TYieldStress = CYieldStress;
      TReverStrain = CReverStrain;
      TReverStress = CReverStress;
      TPlasticExcursion = CPlasticExcursion;

      // State variables from last converged state
      Cstrain = data(21);
      Cstress = data(22);
      Ctangent = data(23);
      CcurR = data(24);

      // Copy converged state values into trial values
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
      TcurR = CcurR;
   }
    
   return res;
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
	Tstrain   = (*theVector)(2); 
	Tstress   = (*theVector)(3); 
	Ttangent     = (*theVector)(4); 
	//epscp  = (*theVector)(5); 
  } else if (responseID == 101){ // get var
    Vector *theVector = matInfo.theVector;
    //(*theVector)(0) = ecminP;
    //(*theVector)(1) = deptP;
    (*theVector)(2) = Tstrain;
    (*theVector)(3) = Tstress;
    (*theVector)(4) = Ttangent;
	//(*theVector)(5) = epscp;
  } else
    return this->UniaxialMaterial::getResponse(responseID, matInfo);

  return 0;
}

void 
MCFTSteel03::Print(OPS_Stream &s, int flag)
{
   s << "\t\tMCFTSteel03 tag: " << this->getTag() << endln;
   s << "\t\t  fy: " << fy << " ";
   s << "\t\t  E0: " << E0 << " ";
   s << "\t\t  b: " << b << " ";
   s << "\t\t  r0:  " << r0 << endln;
   s << "\t\t  CR1: " << coeffR1 << " ";
   s << "\t\t  CR2: " << coeffR2 << " ";   
   s << "\t\t  a1: " << a1 << " ";
   s << "\t\t  a2: " << a2 << " ";
   s << "\t\t  a3: " << a3 << " ";
   s << "\t\t  a4: " << a4 << " ";
   s << endln;
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
MCFTSteel03::setParameter(const char **argv, int argc, Parameter &param)
{
	if (argc < 1)
		return -1;

	if ((strcmp(argv[0],"sigmaY") == 0)||(strcmp(argv[0],"fy") == 0)) {
		return param.addObject(1, this);
		
	}
	if (strcmp(argv[0],"E") == 0) {
		return param.addObject(2, this);
	}
	if (strcmp(argv[0],"b") == 0) {
		return param.addObject(3, this);
	}
	
	else
		opserr << "WARNING: Could not set parameter in MCFTSteel03. " << endln;
               
	return -1;
}

int
MCFTSteel03::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		this->fy = info.theDouble;
		break;
	case 2:
		this->E0 = info.theDouble;
		break;
	case 3:
		this->b = info.theDouble;
		break;
	default:
		return -1;
	}

//	Ttangent = E0;          // Initial stiffness
	revertToStart();

	return 0;
}

int
MCFTSteel03::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}

double MCFTSteel03::getStrainSensitivity(int gradNumber){
	if (SHVs ==0) {
			opserr<<"warning:MCFTSteel03::getStrainsSensitivity, SHVs =0 "<<endln;	
			return 0.0;
		} 
	double  Tsensitivity =(*SHVs)(0,(gradNumber));
	return Tsensitivity;
}

double
MCFTSteel03::getStressSensitivity(int gradNumber, bool conditional){

/*	if (conditional==false)     {
		if (SHVs ==0) {
			opserr<<"warning:MCFTSteel03::getStressSensitivity, SHVs =0 "<<endln;	
			return 0.0;
		} 
		double  Tsensitivity =(*SHVs)(1,(gradNumber-1));
		return Tsensitivity;
	}
*/
	// Initialize return value

	// Pick up sensitivity history variables
	double strainSensitivity = 0.0;
	double stressSensitivity = 0.0;
	double minStrainSensitivity = 0.0;
	double maxStrainSensitivity = 0.0;
	double plasticExcursionSensitivity = 0.0;
	double yieldStrainSensitivity = 0.0;
	double yieldStressSensitivity = 0.0;
	double reverStrainSensitivity = 0.0;
	double reverStressSensitivity = 0.0;

	if (SHVs != 0) {
		strainSensitivity = (*SHVs)(0,(gradNumber));
		stressSensitivity = (*SHVs)(1,(gradNumber));
	
		minStrainSensitivity = (*SHVs)(2,(gradNumber));
		maxStrainSensitivity = (*SHVs)(3,(gradNumber));
		plasticExcursionSensitivity = (*SHVs)(4,(gradNumber));
		yieldStrainSensitivity = (*SHVs)(5,(gradNumber));
		yieldStressSensitivity = (*SHVs)(6,(gradNumber));
		reverStrainSensitivity = (*SHVs)(7,(gradNumber));
		reverStressSensitivity = (*SHVs)(8,(gradNumber));

	}

	// Assign values to parameter derivatives (depending on what's random)
	double fySensitivity = 0.0;
	double E0Sensitivity = 0.0;
	double bSensitivity = 0.0;	  

	if (parameterID == 1) {
		fySensitivity = 1.0;
	}
	else if (parameterID == 2) {
		E0Sensitivity = 1.0;
	}
	else if (parameterID == 3) {
		bSensitivity = 1.0;
	}

//--------------------------------------------------------------------------
	  double strain =Tstrain;
	  double dStrain = strain - Cstrain;
//	  if ((Cloading == 0) && (dStrain ==0)) {/*opserr<<"warning: MCFTSteel03::getStressSensitivity, dStrain is zero! "<<endln; return 0.0;*/} 
      if (fabs(dStrain) < DBL_EPSILON) {  return stressSensitivity-Ctangent*strainSensitivity;}

	  double Esh = b*E0;
	  double epsy = fy/E0;
      double EshSensitivity = bSensitivity*E0+b*E0Sensitivity;
	  double epsySensitivity = fySensitivity/E0-E0Sensitivity*fy/pow(E0,2);

	  double epMax;
	  double epMaxSensitivity;

	  if (fabs(CminStrain)>fabs(CmaxStrain)){
	      epMax = fabs(CminStrain);

		  double sign;
		  if (CminStrain>=0) sign =1;
		  else sign=-1;
		  epMaxSensitivity = sign*minStrainSensitivity;
	  }
	  else {
		  epMax =fabs(CmaxStrain);
		  double sign;
		  if (CmaxStrain>=0) sign =1;
		  else sign=-1;
		  epMaxSensitivity = sign*maxStrainSensitivity;
	  }

	  double reverStrain=CReverStrain ;
	  double reverStress=CReverStress ;
	  double yieldStrain=CYieldStrain;
	  double yieldStress=CYieldStress;
	  double plasticExcursion=CPlasticExcursion;

      // Determine initial loading condition
      if (Cloading == 0 && dStrain != 0.0) {
		
		minStrainSensitivity = -epsySensitivity;
		maxStrainSensitivity =  epsySensitivity;
		epMaxSensitivity = epsySensitivity;

		if (dStrain > 0.0) {
		  //Tloading = 1;
		  yieldStrain = epMax;
		  yieldStrainSensitivity = epMaxSensitivity;
		  yieldStress = fy;
		  yieldStressSensitivity = fySensitivity;
		  plasticExcursion = epMax;
		  plasticExcursionSensitivity = epMaxSensitivity;
		}
		else {
		  yieldStrain = -epMax;
		  yieldStrainSensitivity = -epMaxSensitivity;
		  yieldStress = -fy;
		  yieldStressSensitivity = -fySensitivity;
		  plasticExcursion = -epMax;
		  plasticExcursionSensitivity = -epMaxSensitivity;
		}

		double epStar = strain/yieldStrain;
		double epStarSensitivity = -strain *yieldStrainSensitivity/pow(yieldStrain,2);

		double sigmaStar = b * epStar+(1-b)*epStar/pow((1+pow(fabs(epStar),r0)),1/r0);
		
		double sign;
		if (epStar>=0) sign =1;
		else sign=-1;

		double sigmaStarSensitivity = bSensitivity*epStar+b*epStarSensitivity-bSensitivity*epStar/pow((1+pow(fabs(epStar),r0)),1/r0)
		  +(1-b)*epStarSensitivity/pow((1+pow(fabs(epStar),r0)),1/r0)
		  -(1-b)*epStar/pow((1+pow(fabs(epStar),r0)),1/r0)*
		  pow(fabs(epStar),r0)/r0/(1+pow(fabs(epStar),r0)*sign*epStarSensitivity*r0/fabs(epStar));

		stressSensitivity= sigmaStarSensitivity * yieldStress+sigmaStar * yieldStressSensitivity;
	
	 }
  
      // Transition from loading to unloading, i.e. positive strain increment
      // to negative strain increment

      if (Cloading == 1 && dStrain < 0.0) {

		  reverStrain = Cstrain;
		  reverStress = Cstress;
		  reverStrainSensitivity = strainSensitivity;
		  reverStressSensitivity = stressSensitivity;
		  
		  double maxStrain=CmaxStrain;
		  if (reverStrain > maxStrain) {
			  maxStrain = reverStrain;
			  maxStrainSensitivity = reverStrainSensitivity;
		  }

	  	  if (fabs(CminStrain)>fabs(maxStrain)){
			double sign;
			if (CminStrain>=0) sign =1;
			else sign=-1;

			epMax = fabs(CminStrain);
			epMaxSensitivity = sign*minStrainSensitivity;
		  }
		  else {
			  double sign;
			  if (maxStrain>=0) sign =1;
			  else sign=-1;

			  epMax =fabs(maxStrain);
			  epMaxSensitivity = sign*maxStrainSensitivity;
		  }

		  double sigmaShift = fy * a1 *(epMax/epsy-a2);
		  double sigmaShiftSensitivity;
		  if (sigmaShift <=0) {
			  sigmaShift = 0;
			  sigmaShiftSensitivity =0.0;
		  }

		  else{
			sigmaShiftSensitivity = fySensitivity * a1 *(epMax/epsy-a2) +fy * a1 *epMaxSensitivity/epsy - fy * a1 *epMax*epsySensitivity/pow(epsy,2);
		  }
		  
		  yieldStrain = (reverStress+fy+sigmaShift-(E0*reverStrain+Esh*epsy))/(Esh-E0);
		  yieldStress = Esh * (yieldStrain+epsy)-fy-sigmaShift;
		  plasticExcursion=CminStrain;

		  yieldStrainSensitivity = (reverStressSensitivity+fySensitivity+sigmaShiftSensitivity-(E0Sensitivity*reverStrain
			  +E0*reverStrainSensitivity +EshSensitivity*epsy+Esh*epsySensitivity))
			  /(Esh-E0)-(EshSensitivity-E0Sensitivity)*yieldStrain/(Esh-E0);

		  yieldStressSensitivity=EshSensitivity*(yieldStrain+epsy)+Esh*(yieldStrainSensitivity+epsySensitivity)-fySensitivity-sigmaShiftSensitivity;

		  plasticExcursionSensitivity=minStrainSensitivity;
		  
      }

      // Transition from unloading to loading, i.e. negative strain increment
      // to positive strain increment
      if (Cloading == -1 && dStrain > 0.0) {

		  reverStrain = Cstrain;
		  reverStress = Cstress;
		  reverStrainSensitivity = strainSensitivity;
		  reverStressSensitivity = stressSensitivity;

		  double minStrain = CminStrain;
		  if (reverStrain < minStrain) {
			  minStrain = reverStrain;
			  minStrainSensitivity = reverStrainSensitivity;
		  }

	  	  if (fabs(minStrain)>fabs(CmaxStrain)){
			  double sign;
			  if (minStrain>=0) sign =1;
			  else sign=-1;

			  epMax = fabs(minStrain);
			  epMaxSensitivity = sign*minStrainSensitivity;
		  }
		  else {
  			  double sign;
			  if (CmaxStrain>=0) sign =1;
			  else sign=-1;

			  epMax =fabs(CmaxStrain);
			  epMaxSensitivity = sign*maxStrainSensitivity;
		  }

		  double sigmaShift = fy * a1 *(epMax/epsy-a2);
		  double sigmaShiftSensitivity;
		  if (sigmaShift <=0) {
			  sigmaShift = 0;
			  sigmaShiftSensitivity =0.0;
		  }

		  else{
			sigmaShiftSensitivity = fySensitivity * a1 *(epMax/epsy-a2) +fy * a1 *epMaxSensitivity/epsy - fy * a1 *epMax*epsySensitivity/pow(epsy,2);
		  }
		  
		  yieldStrain = (reverStress+Esh*epsy-(E0*reverStrain+fy+sigmaShift))/(Esh-E0);
		  yieldStress = Esh * (yieldStrain-epsy)+fy+sigmaShift;
		  plasticExcursion=CmaxStrain;

		  yieldStrainSensitivity = (reverStressSensitivity+EshSensitivity*epsy+Esh*epsySensitivity-
			  (E0Sensitivity*reverStrain +E0*reverStrainSensitivity + fySensitivity+sigmaShiftSensitivity))/ (Esh-E0)
			  -(EshSensitivity-E0Sensitivity)*yieldStrain/(Esh-E0);

		  yieldStressSensitivity=EshSensitivity*(yieldStrain-epsy)+Esh*(yieldStrainSensitivity-epsySensitivity)+fySensitivity+sigmaShiftSensitivity;

		  plasticExcursionSensitivity=maxStrainSensitivity;

	  }
      
      if (Cloading != 0) {

		  double Xi = fabs((plasticExcursion-yieldStrain)/epsy);
		  double R  = r0-(coeffR1*Xi)/(coeffR2+Xi);
		  
		  double sign;
		  if (((plasticExcursion-yieldStrain)/epsy)>=0) sign =1;
		  else sign=-1;

		  double XiSensitivity = sign*((plasticExcursionSensitivity-yieldStrainSensitivity)/epsy-(plasticExcursion-yieldStrain)*epsySensitivity/pow(epsy,2));
		  double RSensitivity  = -(coeffR1*XiSensitivity)/(coeffR2+Xi)+(coeffR1*Xi)*XiSensitivity/pow((coeffR2+Xi),2);

		  double epStar = (strain-reverStrain)/(yieldStrain-reverStrain);
		  double epStarSensitivity = (-reverStrainSensitivity)/(yieldStrain-reverStrain)-(strain-reverStrain)*(yieldStrainSensitivity-reverStrainSensitivity)/
			  pow((yieldStrain-reverStrain),2);

		  double sigmaStar = b * epStar+(1-b)*epStar/pow((1+pow(fabs(epStar),R)),1/R);

  		  
		  if (epStar>=0) sign =1;
		  else sign=-1;

		  double sigmaStarSensitivity = bSensitivity*epStar+b*epStarSensitivity-bSensitivity*epStar/pow((1+pow(fabs(epStar),R)),1/R)+(1-b)*epStarSensitivity/pow((1+pow(fabs(epStar),R)),1/R)
			  -(1-b)*epStar/pow((1+pow(fabs(epStar),R)),1/R)*
			  (-RSensitivity*log(1+pow(fabs(epStar),R))/pow(R,2)+pow(fabs(epStar),R)/R/(1+pow(fabs(epStar),R))*(log(fabs(epStar))*RSensitivity
			  +sign*epStarSensitivity*R/fabs(epStar)));

		  stressSensitivity= sigmaStarSensitivity * (yieldStress-reverStress)+sigmaStar * (yieldStressSensitivity-reverStressSensitivity)+reverStressSensitivity;

      }

	return stressSensitivity;
}

/*
double
MCFTSteel03::getInitialTangentSensitivity(int gradNumber)
{
	// For now, assume that this is only called for initial stiffness 
	if (parameterID == 2) {
		return 1.0; 
	}
	else {
		return 0.0;
	}
}*/

int
MCFTSteel03::commitSensitivity(double TstrainSensitivity, int gradNumber, int numGrads)
{
	if (SHVs == 0) {
		SHVs = new Matrix(9,numGrads);
		SHVs->Zero();
	}

	// Pick up sensitivity history variables
	double strainSensitivity = 0.0;
	double stressSensitivity = 0.0;
	double minStrainSensitivity = 0.0;
	double maxStrainSensitivity = 0.0;
	double plasticExcursionSensitivity = 0.0;
	double yieldStrainSensitivity = 0.0;
	double yieldStressSensitivity = 0.0;
	double reverStrainSensitivity = 0.0;
	double reverStressSensitivity = 0.0;

	if (SHVs != 0) {
	  strainSensitivity = (*SHVs)(0,(gradNumber));
	  stressSensitivity = (*SHVs)(1,(gradNumber));
	  
	  minStrainSensitivity = (*SHVs)(2,(gradNumber));
	  maxStrainSensitivity = (*SHVs)(3,(gradNumber));
	  plasticExcursionSensitivity = (*SHVs)(4,(gradNumber));
	  yieldStrainSensitivity = (*SHVs)(5,(gradNumber));
	  yieldStressSensitivity = (*SHVs)(6,(gradNumber));
	  reverStrainSensitivity = (*SHVs)(7,(gradNumber));
	  reverStressSensitivity = (*SHVs)(8,(gradNumber));
	}

	// Assign values to parameter derivatives (depending on what's random)
	double fySensitivity = 0.0;
	double E0Sensitivity = 0.0;
	double bSensitivity = 0.0;	  

	if (parameterID == 1) {
		fySensitivity = 1.0;
	}
	else if (parameterID == 2) {
		E0Sensitivity = 1.0;
	}
	else if (parameterID == 3) {
		bSensitivity = 1.0;
	}

//--------------------------------------------------------------------------
	  double strain =Tstrain;
	  double dStrain = strain - Cstrain;
//	  if (dStrain ==0){opserr<<"warning: MCFTSteel03::getStressSensitivity, dStrain is zero! "<<endln;} 
//	  if ((Cloading == 0) && (dStrain ==0)) {/*opserr<<"warning: MCFTSteel03::getStressSensitivity, dStrain is zero! "<<endln; return 0.0;*/} 
 //     if (fabs(dStrain) > DBL_EPSILON) return 0;
	  if (fabs(dStrain) < DBL_EPSILON) {  return 0;}
 
	  double Esh = b*E0;
	  double epsy = fy/E0;
      double EshSensitivity = bSensitivity*E0+b*E0Sensitivity;
	  double epsySensitivity = fySensitivity/E0-E0Sensitivity*fy/pow(E0,2);

	  double epMax;
	  double epMaxSensitivity;

	  if (fabs(CminStrain)>fabs(CmaxStrain)){
	    epMax = fabs(CminStrain);
	    
	    double sign;
	    if (CminStrain>=0) sign =1;
	    else sign=-1;
	    epMaxSensitivity = sign*minStrainSensitivity;
	  }
	  else {
	    epMax =fabs(CmaxStrain);
	    double sign;
	    if (CmaxStrain>=0) sign =1;
	    else sign=-1;
	    epMaxSensitivity = sign*maxStrainSensitivity;
	  }
	  double reverStrain=CReverStrain ;
	  double reverStress=CReverStress ;
	  double yieldStrain=CYieldStrain;
	  double yieldStress=CYieldStress;
	  double plasticExcursion=CPlasticExcursion;

      // Determine initial loading condition
      if (Cloading == 0 && dStrain != 0.0) {
		
		minStrainSensitivity = -epsySensitivity;
		maxStrainSensitivity =  epsySensitivity;
		epMaxSensitivity = epsySensitivity;

		if (dStrain > 0.0) {
		  //Tloading = 1;
		  yieldStrain = epMax;
		  yieldStrainSensitivity = epMaxSensitivity;
		  yieldStress = fy;
		  yieldStressSensitivity = fySensitivity;
		  plasticExcursion = epMax;
		  plasticExcursionSensitivity = epMaxSensitivity;
		}
		else {
		  yieldStrain = -epMax;
		  yieldStrainSensitivity = -epMaxSensitivity;
		  yieldStress = -fy;
		  yieldStressSensitivity = -fySensitivity;
		  plasticExcursion = -epMax;
		  plasticExcursionSensitivity = -epMaxSensitivity;
		}

		double epStar = strain/yieldStrain;

		double epStarSensitivity = TstrainSensitivity/yieldStrain-strain*yieldStrainSensitivity/pow(yieldStrain,2);
		
		double sigmaStar = b * epStar+(1-b)*epStar/pow((1+pow(fabs(epStar),r0)),1/r0);
		
		double sign;
		if (epStar>=0) sign =1;
		else sign=-1;

/*		double sigmaStarSensitivity = bSensitivity*epStar + b*epStarSensitivity - bSensitivity*epStar/pow((1+pow(fabs(epStar),r0)),1/r0) 
		  + (1-b)*epStarSensitivity/pow((1+pow(fabs(epStar),r0)),1/r0)
		  -(1-b)*epStar/pow((1+pow(fabs(epStar),r0)),1/r0)*sign*epStarSensitivity*r0/fabs(epStar);

*/	
		double sigmaStarSensitivity = bSensitivity*epStar+b*epStarSensitivity-bSensitivity*epStar/pow((1+pow(fabs(epStar),r0)),1/r0)
		  +(1-b)*epStarSensitivity/pow((1+pow(fabs(epStar),r0)),1/r0)
		  -(1-b)*epStar/pow((1+pow(fabs(epStar),r0)),1/r0)*
		  pow(fabs(epStar),r0)/r0/(1+pow(fabs(epStar),r0)*sign*epStarSensitivity*r0/fabs(epStar));

		stressSensitivity= sigmaStarSensitivity * yieldStress+sigmaStar * yieldStressSensitivity;
	
	 }
  
      // Transition from loading to unloading, i.e. positive strain increment
      // to negative strain increment

      if (Cloading == 1 && dStrain < 0.0) {
	
		  
		  reverStrain = Cstrain;
		  reverStress = Cstress;
		  reverStrainSensitivity = strainSensitivity;
		  reverStressSensitivity = stressSensitivity;
		  
		  double maxStrain=CmaxStrain;
		  if (reverStrain > maxStrain) {
			  maxStrain = reverStrain;
			  maxStrainSensitivity = reverStrainSensitivity;
		  }

		
	  	  if (fabs(CminStrain)>fabs(maxStrain)){
			double sign;
			if (CminStrain>=0) sign =1;
			else sign=-1;

			epMax = fabs(CminStrain);
			epMaxSensitivity = sign*minStrainSensitivity;
		  }
		  else {
			  double sign;
			  if (maxStrain>=0) sign =1;
			  else sign=-1;

			  epMax =fabs(maxStrain);
			  epMaxSensitivity = sign*maxStrainSensitivity;
		  }


		  double sigmaShift = fy * a1 *(epMax/epsy-a2);
		  double sigmaShiftSensitivity;
		  if (sigmaShift <=0) {
			  sigmaShift = 0;
			  sigmaShiftSensitivity =0.0;
		  }

		  else{
			sigmaShiftSensitivity = fySensitivity * a1 *(epMax/epsy-a2) +fy * a1 *epMaxSensitivity/epsy - fy * a1 *epMax*epsySensitivity/pow(epsy,2);
		  }
		  
		  yieldStrain = (reverStress+fy+sigmaShift-(E0*reverStrain+Esh*epsy))/(Esh-E0);
		  yieldStress = Esh * (yieldStrain+epsy)-fy-sigmaShift;
		  plasticExcursion=CminStrain;

		  yieldStrainSensitivity = (reverStressSensitivity+fySensitivity+sigmaShiftSensitivity-(E0Sensitivity*reverStrain
			  +E0*reverStrainSensitivity +EshSensitivity*epsy+Esh*epsySensitivity))
			  /(Esh-E0)-(EshSensitivity-E0Sensitivity)*yieldStrain/(Esh-E0);

		  yieldStressSensitivity=EshSensitivity*(yieldStrain+epsy)+Esh*(yieldStrainSensitivity+epsySensitivity)-fySensitivity-sigmaShiftSensitivity;

		  plasticExcursionSensitivity=minStrainSensitivity;
		  
      }

      // Transition from unloading to loading, i.e. negative strain increment
      // to positive strain increment
      if (Cloading == -1 && dStrain > 0.0) {

		  reverStrain = Cstrain;
		  reverStress = Cstress;
		  reverStrainSensitivity = strainSensitivity;
		  reverStressSensitivity = stressSensitivity;

		  double minStrain = CminStrain;
		  if (reverStrain < minStrain) {
			  minStrain = reverStrain;
			  minStrainSensitivity = reverStrainSensitivity;
		  }

		
	  	  if (fabs(minStrain)>fabs(CmaxStrain)){
			  double sign;
			  if (minStrain>=0) sign =1;
			  else sign=-1;

			  epMax = fabs(minStrain);
			  epMaxSensitivity = sign*minStrainSensitivity;
		  }
		  else {
  			  double sign;
			  if (CmaxStrain>=0) sign =1;
			  else sign=-1;

			  epMax =fabs(CmaxStrain);
			  epMaxSensitivity = sign*maxStrainSensitivity;
		  }

		  double sigmaShift = fy * a1 *(epMax/epsy-a2);
		  double sigmaShiftSensitivity;
		  if (sigmaShift <=0) {
			  sigmaShift = 0;
			  sigmaShiftSensitivity =0.0;
		  }

		  else{
			sigmaShiftSensitivity = fySensitivity * a1 *(epMax/epsy-a2) +fy * a1 *epMaxSensitivity/epsy - fy * a1 *epMax*epsySensitivity/pow(epsy,2);
		  }
		  
		  yieldStrain = (reverStress+Esh*epsy-(E0*reverStrain+fy+sigmaShift))/(Esh-E0);
		  yieldStress = Esh * (yieldStrain-epsy)+fy+sigmaShift;
		  plasticExcursion=CmaxStrain;

		  yieldStrainSensitivity = (reverStressSensitivity+EshSensitivity*epsy+Esh*epsySensitivity-
			  (E0Sensitivity*reverStrain +E0*reverStrainSensitivity + fySensitivity+sigmaShiftSensitivity))/ (Esh-E0)
			  -(EshSensitivity-E0Sensitivity)*yieldStrain/(Esh-E0);

		  yieldStressSensitivity=EshSensitivity*(yieldStrain-epsy)+Esh*(yieldStrainSensitivity-epsySensitivity)+fySensitivity+sigmaShiftSensitivity;

		  plasticExcursionSensitivity=maxStrainSensitivity;

	  }
      
      if (Cloading != 0) {

		  double Xi = fabs((plasticExcursion-yieldStrain)/epsy);
		  double R  = r0-(coeffR1*Xi)/(coeffR2+Xi);
		  
		  double sign;
		  if (((plasticExcursion-yieldStrain)/epsy)>=0) sign =1;
		  else sign=-1;

		  double XiSensitivity = sign*((plasticExcursionSensitivity-yieldStrainSensitivity)/epsy-(plasticExcursion-yieldStrain)*epsySensitivity/pow(epsy,2));
		  double RSensitivity  = -(coeffR1*XiSensitivity)/(coeffR2+Xi)+(coeffR1*Xi)*XiSensitivity/pow((coeffR2+Xi),2);

		  double epStar = (strain-reverStrain)/(yieldStrain-reverStrain);
		  double epStarSensitivity = (TstrainSensitivity-reverStrainSensitivity)/(yieldStrain-reverStrain)-(strain-reverStrain)*(yieldStrainSensitivity-reverStrainSensitivity)/
			  pow((yieldStrain-reverStrain),2);

		  double sigmaStar = b * epStar+(1-b)*epStar/pow((1+pow(fabs(epStar),R)),1/R);

  		  
		  if (epStar>=0) sign =1;
		  else sign=-1;

		  double sigmaStarSensitivity = bSensitivity*epStar+b*epStarSensitivity-bSensitivity*epStar/pow((1+pow(fabs(epStar),R)),1/R)+(1-b)*epStarSensitivity/pow((1+pow(fabs(epStar),R)),1/R)
			  -(1-b)*epStar/pow((1+pow(fabs(epStar),R)),1/R)*
			  (-RSensitivity*log(1+pow(fabs(epStar),R))/pow(R,2)+pow(fabs(epStar),R)/R/(1+pow(fabs(epStar),R))*(log(fabs(epStar))*RSensitivity
			  +sign*epStarSensitivity*R/fabs(epStar)));
		  stressSensitivity= sigmaStarSensitivity * (yieldStress-reverStress)+sigmaStar * (yieldStressSensitivity-reverStressSensitivity)+reverStressSensitivity;


      }






	// Commit history variables
	(*SHVs)(0,(gradNumber))=TstrainSensitivity ;
	(*SHVs)(1,(gradNumber))=stressSensitivity ;
	
	(*SHVs)(2,(gradNumber))=minStrainSensitivity;
	(*SHVs)(3,(gradNumber))=maxStrainSensitivity;
	(*SHVs)(4,(gradNumber))=plasticExcursionSensitivity;
	(*SHVs)(5,(gradNumber))=yieldStrainSensitivity;
	(*SHVs)(6,(gradNumber))=yieldStressSensitivity;
	(*SHVs)(7,(gradNumber))=reverStrainSensitivity;
	(*SHVs)(8,(gradNumber))=reverStressSensitivity;
	return 0;
}

double
MCFTSteel03::getInitialTangentSensitivity(int gradNumber)
{
	if (parameterID == 2) {
		return 1.0;
	}
	else {
		return 0.0;
	}
}

// AddingSensitivity:END /////////////////////////////////////////////