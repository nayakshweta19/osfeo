/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.6 $
// $Date: 2008-08-26 16:22:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/BoucWen3DMaterial.cpp,v $

#include <BoucWen3DMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>

#include <MaterialResponse.h>
#include <DummyStream.h>
#include <elementAPI.h>

#define OPS_Export 

static int numBoucWen3DMaterials = 0;

OPS_Export void *
OPS_NewBoucWen3DMaterial()
{
  if (numBoucWen3DMaterials == 0) {
    numBoucWen3DMaterials++;
  }

  // Pointer to a uniaxial material that will be returned
  NDMaterial *theMaterial = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 12) {
    opserr << "Want: NDMaterial BoucWen3D matTag? rho? tag? alpha? ko? n? gamma? epsc0? beta? Ao? deltaA? deltaNu? deltaEta?>\n";
	return 0;
  }

  int tag;
  double rho;
  int    iData[4];
  double dData[8];
  int numData = 0;

  numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid NDMaterial BoucWen3D tag" << endln;
    return 0;
  }

  numData = 9;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "Invalid Arg rho: NDMaterial BoucWen3D tag? alpha? ko? n? gamma? beta? Ao? deltaA? deltaNu? deltaEta?" << endln;
    return 0;	
  }

  // Check if the user has given a tolerance for the Newton scheme
  double tolerance = 1.0e-8;
  numData = 1;
  if (numRemainingArgs > 12) {
    if (OPS_GetDouble(&numData, &tolerance) != 0) {
      opserr << "Invalid Arg rho: NDMaterial BoucWen3D tag? alpha? ko? n? gamma? beta? Ao? deltaA? deltaNu? deltaEta?" << endln;
      return 0;	
    }
  }
    
  // Check if the user has given a maxNumIter for the Newton scheme
  int maxNumIter = 20;
  if (numRemainingArgs > 13) {
    if (OPS_GetInt(&numData, &maxNumIter) != 0) {
      opserr << "Invalid Arg rho: NDMaterial BoucWen3D tag? alpha? ko? n? gamma? beta? Ao? deltaA? deltaNu? deltaEta?" << endln;
      return 0;	
    }
  }

  theMaterial = new BoucWen3DMaterial (tag, 
						     dData[0],
						     dData[1],
						     dData[2],
						     dData[3],
						     dData[4],
						     dData[5],
						     dData[6],
						     dData[7],
							 dData[8],
							 tolerance, maxNumIter);
       
  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory creating material\n";
    opserr << "BoucWen3D: " << tag << endln;
    return 0;
  }

  return theMaterial;
}

BoucWen3DMaterial::BoucWen3DMaterial(int tag, 
					double p_alpha,
					double p_ko,
					double p_n,
					double p_gamma,
					double p_beta,
					double p_Ao,
					double p_deltaA,
					double p_deltaNu,
					double p_deltaEta,
					double ptolerance,
					int pMaxNumIter)
:NDMaterial(tag, MAT_TAG_BoucWen3D),
alpha(p_alpha), ko(p_ko), n(p_n), gamma(p_gamma), beta(p_beta), Ao(p_Ao), 
deltaA(p_deltaA), deltaNu(p_deltaNu), deltaEta(p_deltaEta), tolerance(ptolerance),
maxNumIter(pMaxNumIter)
{
	parameterID = 0;
	SHVs = 0;

	// Initialize variables
    this->revertToStart();
}

BoucWen3DMaterial::~BoucWen3DMaterial()
{
	if (SHVs != 0) 
		delete SHVs;
}

double 
BoucWen3DMaterial::signum(double value)
{
	if (value > 0.0) {
		return 1.0;
	}
	else {
		return -1.0;
	}
}

int 
BoucWen3DMaterial::setTrialStrain (const Vector &strain_from_element)
{
	// Set trial strain and compute strain increment
	/*Tstrain = strain_from_element;
	Vector dStrain = Tstrain - Cstrain;

	// Initial declarations (make sure not to declare class variables here!)
	double TA, Tnu, Teta, Psi, Phi, f, Te_, TA_, Tnu_;
	double Teta_, Phi_, f_, Tznew, Tzold, sign;

	// Newton-Raphson scheme to solve for z_{i+1} := z1
	int count = 0;
	double startPoint = 0.01;
	Tz = startPoint;
	Tzold = startPoint;
	Tznew = 1.0;  
	while ( ( fabs(Tzold-Tznew) > tolerance ) && count<maxNumIter) {

		Te = Ce + (1-alpha)*ko*dStrain*Tz;
		TA = Ao - deltaA*Te;
		Tnu = 1.0 + deltaNu*Te;
		Teta = 1.0 + deltaEta*Te;
		sign = signum(dStrain*Tz);
		Psi = gamma + beta*sign;
		Phi = TA - pow(fabs(Tz),n)*Psi*Tnu;
		f = Tz - Cz - Phi/Teta*dStrain;

		// Evaluate function derivative f' (underscore:=prime)
		Te_ = (1.0-alpha)*ko*dStrain;
		TA_ = -deltaA*Te_;
		Tnu_ = deltaNu*Te_;
		Teta_ = deltaEta*Te_;
		sign = signum(Tz);
		double pow1;
		double pow2;
		if (Tz == 0.0) {
			pow1 = 0.0;
			pow2 = 0.0;
		}
		else {
			pow1 = pow(fabs(Tz),(n-1));
			pow2 = pow(fabs(Tz),n);
		}
		Phi_ = TA_ - n*pow1*sign*Psi*Tnu - pow2*Psi*Tnu_;
		f_ = 1.0 - (Phi_*Teta-Phi*Teta_)/pow(Teta,2.0)*dStrain;

		// Issue warning if derivative is zero
		if ( fabs(f_)<1.0e-10 ) {
			opserr << "WARNING: BoucWen3DMaterial::setTrialStrain() -- zero derivative " << endln
				<< " in Newton-Raphson scheme" << endln;
		}

		// Take a Newton step
		Tznew = Tz - f/f_;

		// Update the root (but the keep the old for convergence check)
		Tzold = Tz; 
		Tz = Tznew;

		// Update counter
		count++;

		// Issue warning if we didn't converge
		if (count == maxNumIter) {
			opserr << "WARNING: BoucWen3DMaterial::setTrialStrain() -- did not" << endln
				<< " find the root z_{i+1}, after " << maxNumIter << " iterations" << endln
				<< " and norm: " << fabs(Tzold-Tznew) << endln;
		}

		// Compute stress
		Tstress = alpha*ko*Tstrain + (1-alpha)*ko*Tz;

		// Compute deterioration parameters
		Te = Ce + (1-alpha)*ko*dStrain*Tz;
		TA = Ao - deltaA*Te;
		Tnu = 1.0 + deltaNu*Te;
		Teta = 1.0 + deltaEta*Te;

		// Compute tangent
		if (Tz != 0.0) {
			Psi = gamma + beta*signum(dStrain*Tz);
			Phi = TA - pow(fabs(Tz),n)*Psi*Tnu;
			double b1, b2, b3, b4, b5;
			b1  = (1-alpha)*ko*Tz;
			b2  = (1-alpha)*ko*dStrain;
			b3  = dStrain/Teta;
			b4  = -b3*deltaA*b1 - b3*pow(fabs(Tz),n)*Psi*deltaNu*b1 
				- Phi/(Teta*Teta)*dStrain*deltaEta*b1 + Phi/Teta;
			b5  = 1.0 + b3*deltaA*b2 + b3*n*pow(fabs(Tz),(n-1))*signum(Tz)*Psi*Tnu
				+ b3*pow(fabs(Tz),n)*Psi*deltaNu*b2
				+ Phi/(Teta*Teta)*dStrain*deltaEta*b2;
			double DzDeps = b4/b5;
			Ttangent = alpha*ko + (1-alpha)*ko*DzDeps;
		}
		else {
			Ttangent = alpha*ko + (1-alpha)*ko;
		}
	}
	*/
    return 0;
}

//unused trial strain functions
int
BoucWen3DMaterial::setTrialStrain( const Vector &v, const Vector &r )
{ 
   return this->setTrialStrain( v ) ;
} 

//this is mike's problem
int
BoucWen3DMaterial::setTrialStrain(const Tensor &v) 
{
  return -1 ;
}

int
BoucWen3DMaterial::setTrialStrain(const Tensor &v, const Tensor &r)     
{
  return -1 ;
}

int
BoucWen3DMaterial::setTrialStrainIncr(const Vector &v) 
{/*
  static Vector newStrain(6);
  newStrain(0) = Cstrain(0,0) + v(0);
  newStrain(1) = strain(1,1) + v(1);
  newStrain(2) = strain(2,2) + v(2);
  newStrain(3) = 2.0*strain(0,1) + v(3);
  newStrain(4) = 2.0*strain(1,2) + v(4);
  newStrain(5) = 2.0*strain(2,0) + v(5);
  
  return this->setTrialStrain(newStrain);*/
  return 0;
}

int
BoucWen3DMaterial::setTrialStrainIncr(const Vector &v, const Vector &r) 
{
  return this->setTrialStrainIncr(v);
}

const Vector& 
BoucWen3DMaterial::getStress(void)
{
    return Tstress;
}

const Matrix& 
BoucWen3DMaterial::getInitialTangent(void)
{
    //return ( alpha*ko + (1-alpha)*ko*Ao );
	return Ttangent;
}

const Matrix& 
BoucWen3DMaterial::getTangent(void)
{
    return Ttangent;
}

const Vector& 
BoucWen3DMaterial::getStrain(void)
{
    return Tstrain;
}

int 
BoucWen3DMaterial::commitState(void)
{
    // Commit trial history variables
    Cstrain = Tstrain;
	Cz = Tz;
	Ce = Te;

    return 0;
}

int 
BoucWen3DMaterial::revertToLastCommit(void)
{
	// Nothing to do here
    return 0;
}

int 
BoucWen3DMaterial::revertToStart(void)
{
	Tstrain.Zero();
	Cstrain.Zero();
	Tz.Zero();
	Cz.Zero();
	Te.Zero();
	Ce.Zero();
	Tstress.Zero();
	Ttangent = alpha*ko + (1-alpha)*ko*Ao;

	if (SHVs != 0) 
		SHVs->Zero();

    return 0;
}

NDMaterial *
BoucWen3DMaterial::getCopy()
{
    BoucWen3DMaterial *theCopy =
	new BoucWen3DMaterial(this->getTag(), alpha, ko, n, gamma,
						beta, Ao, deltaA, deltaNu, deltaEta,tolerance,maxNumIter);
    	
    theCopy->Tstrain = Tstrain;
    theCopy->Cstrain = Cstrain;
    theCopy->Tz = Tz;
    theCopy->Cz = Cz;
    theCopy->Te = Te;
    theCopy->Ce = Ce;
    theCopy->Tstress = Tstress;
    theCopy->Ttangent = Ttangent;

    return theCopy;
}

//send back type of material
const char*
BoucWen3DMaterial::getType( ) const 
{
  return "ThreeDimensional" ;
}

int 
BoucWen3DMaterial::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int 
BoucWen3DMaterial::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}

void 
BoucWen3DMaterial::Print(OPS_Stream &s, int flag)
{
    s << "BoucWen3DMaterial, tag: " << this->getTag() << endln;
    s << "  alpha: " << alpha << endln;
    s << "  ko: " << ko << endln;
    s << "  n: " << n << endln;
    s << "  gamma: " << gamma << endln;
    s << "  beta: " << beta << endln;
    s << "  Ao: " << Ao << endln;
    s << "  deltaA: " << deltaA << endln;
    s << "  deltaNu: " << deltaNu << endln;
    s << "  deltaEta: " << deltaEta << endln;
}

int
BoucWen3DMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"alpha") == 0)
    return param.addObject(1, this);
  
  if (strcmp(argv[0],"ko") == 0)
    return param.addObject(2, this);
  
  if (strcmp(argv[0],"n") == 0)
    return param.addObject(3, this);
  
  if (strcmp(argv[0],"gamma") == 0)
    return param.addObject(4, this);
    
  if (strcmp(argv[0],"beta") == 0)
    return param.addObject(5, this);
  
  if (strcmp(argv[0],"Ao") == 0)
    return param.addObject(6, this);
  
  if (strcmp(argv[0],"deltaA") == 0)
    return param.addObject(7, this);
    
  if (strcmp(argv[0],"deltaNu") == 0)
    return param.addObject(8, this);
  
  if (strcmp(argv[0],"deltaEta") == 0)
    return param.addObject(9, this);

  return -1;
}

int
BoucWen3DMaterial::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case 1:
		this->alpha = info.theDouble;
		return 0;
	case 2:
		this->ko = info.theDouble;
		return 0;
	case 3:
		this->n = info.theDouble;
		return 0;
	case 4:
		this->gamma = info.theDouble;
		return 0;
	case 5:
		this->beta = info.theDouble;
		return 0;
	case 6:
		this->Ao = info.theDouble;
		return 0;
	case 7:
		this->deltaA = info.theDouble;
		return 0;
	case 8:
		this->deltaNu = info.theDouble;
		return 0;
	case 9:
		this->deltaEta = info.theDouble;
		return 0;
	default:
		return -1;
	}
}

int
BoucWen3DMaterial::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}

const Vector & 
BoucWen3DMaterial::getStrainSensitivity(int gradIndex)
{
	if (SHVs ==0) return 0.0;
	else{
		double sensitivity =(*SHVs)(4,gradIndex-1); // unconditional stress sensitivity
		return sensitivity;
	}
}

const Vector & 
BoucWen3DMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
	if (conditional == false) {  // return stress sensitivity for recorder purpose
		if (SHVs ==0) return 0.0;
		else {
			double sensitivity =(*SHVs)(3,gradIndex-1); // unconditional stress sensitivity
			return sensitivity;
		}
	}

	// First set values depending on what is random
	double SigmaYSensitivity = 0.0;
	double ESensitivity = 0.0;
	double HkinSensitivity = 0.0;
	double HisoSensitivity = 0.0;

	if (parameterID == 1) {  // sigmaY
		SigmaYSensitivity = 1.0;
	}
	else if (parameterID == 2) {  // E
		ESensitivity = 1.0;
	}
	else if (parameterID == 3) {  // Hkin
		HkinSensitivity = 1.0;
	}
	else if (parameterID == 4) {  // Hiso
		HisoSensitivity = 1.0;
	}
	else {
		// Nothing random here, but may have to return something in any case
	}

	double TStrainSensitivity = 0.0; 

	// Then pick up history variables for this gradient number
	double CPlasticStrainSensitivity = 0.0;
	double CBackStressSensitivity	 = 0.0;
	double CAccumulatedPlasticStrainSensitivity	 = 0.0;
	if (SHVs != 0) {
		CPlasticStrainSensitivity = (*SHVs)(0,gradIndex);
		CBackStressSensitivity	 = (*SHVs)(1,gradIndex);
		CAccumulatedPlasticStrainSensitivity = (*SHVs)(2,gradIndex);
	}

    // ------ Elastic trial -------
/*
	double deltaLambda = 0.0;
	TPlasticStrain = CPlasticStrain;
	TBackStress = CBackStress;
	TAccumulatedPlasticStrain = CAccumulatedPlasticStrain;
    TStress = E * (TStrain-CPlasticStrain);
	double TStressSensitivity = E*(TStrainSensitivity-CPlasticStrainSensitivity)+ESensitivity*(TStrain-CPlasticStrain);
	double CSigmaY = sigmaY+Hiso*CAccumulatedPlasticStrain;

    // Compute trial stress relative to committed back stress
    double xsi = TStress - TBackStress;

    // Compute yield criterion
    double f = fabs(xsi) - (sigmaY + Hiso*TAccumulatedPlasticStrain);

	double sensitivity;


    // Elastic step ... no updates required
    if (f <= -DBL_EPSILON * E) {
		//if (f <= 1.0e-8) {
		// Set trial tangent
		TTangent = E;
		sensitivity = TStressSensitivity;
	}

    // ------- Plastic corrector ... perform return mapping algorithm ---
    else {
      deltaLambda = (fabs(xsi)-CSigmaY)/(E+Hkin+Hiso);
      
	  // Find sign of xsi
      int sign = (xsi < 0) ? -1 : 1;
	  TPlasticStrain = CPlasticStrain +deltaLambda*sign;

	  TBackStress = CBackStress + Hkin*deltaLambda*sign;

	  TAccumulatedPlasticStrain = CAccumulatedPlasticStrain + deltaLambda;

	  TStress = E * (TStrain - TPlasticStrain);
	
      TTangent = E*(Hkin+Hiso) / (E+Hkin+Hiso);

	  double CSigmaYSensitivity = SigmaYSensitivity + HisoSensitivity*CAccumulatedPlasticStrain
		   + Hiso*CAccumulatedPlasticStrainSensitivity;

	  double deltaLambdaSensitivity = ((TStressSensitivity-CBackStressSensitivity)*sign-CSigmaYSensitivity)/(E+Hkin+Hiso)
		  - (ESensitivity+HkinSensitivity+HisoSensitivity)*((E * (TStrain-CPlasticStrain)-CBackStress)*sign-CSigmaY)/pow((E+Hkin+Hiso),2.0);

	  double TPlasticStrainSensitivity = CPlasticStrainSensitivity + deltaLambdaSensitivity*sign;

	  sensitivity = E * (TStrainSensitivity - TPlasticStrainSensitivity)+ESensitivity * (TStrain - TPlasticStrain);

   }*/
	return 0; //sensitivity;
}

const Matrix &
BoucWen3DMaterial::getInitialTangentSensitivity(int gradIndex)
{
	// For now, assume that this is only called for initial stiffness 
	if (parameterID == 2) {
		//return 1.0; 
	}
	else {
		//return 0.0;
	}
	return Ttangent;
}

int
BoucWen3DMaterial::commitSensitivity(Vector & strainGradient, int gradIndex, int numGrads)
{
	if (SHVs == 0) {
		SHVs = new Matrix(5,numGrads);
		SHVs->Zero();
	}

	// First set values depending on what is random
	double SigmaYSensitivity = 0.0;
	double ESensitivity = 0.0;
	double HkinSensitivity = 0.0;
	double HisoSensitivity = 0.0;

	if (parameterID == 1) {  // sigmaY
		SigmaYSensitivity = 1.0;
	}
	else if (parameterID == 2) {  // E
		ESensitivity = 1.0;
	}
	else if (parameterID == 3) {  // Hkin
		HkinSensitivity = 1.0;
	}
	else if (parameterID == 4) {  // Hiso
		HisoSensitivity = 1.0;
	}
	else {
		// Nothing random here, but may have to return something in any case
	}

	// Then pick up history variables for this gradient number

	double CPlasticStrainSensitivity = (*SHVs)(0,gradIndex);
	double CBackStressSensitivity	 = (*SHVs)(1,gradIndex);
	double CAccumulatedPlasticStrainSensitivity = (*SHVs)(2,gradIndex);


    // ------ Elastic trial -------
/*
	double deltaLambda = 0.0;
	TPlasticStrain = CPlasticStrain;
	TBackStress = CBackStress;
	TAccumulatedPlasticStrain = CAccumulatedPlasticStrain;
    TStress = E * (TStrain-CPlasticStrain);
	double TStressSensitivity = E*(TStrainSensitivity-CPlasticStrainSensitivity)+ESensitivity*(TStrain-CPlasticStrain);
	double CSigmaY = sigmaY+Hiso*CAccumulatedPlasticStrain;

    // Compute trial stress relative to committed back stress
    double xsi = TStress - TBackStress;

    // Compute yield criterion
    double f = fabs(xsi) - (sigmaY + Hiso*TAccumulatedPlasticStrain);

	double sensitivity;

    // Elastic step ... no updates required
    if (f <= -DBL_EPSILON * E) {
		//if (f <= 1.0e-8) {
		// Set trial tangent
		TTangent = E;
		sensitivity = TStressSensitivity;
	}

    // ------- Plastic step ... perform return mapping algorithm ---
    else {
      deltaLambda = (fabs(xsi)-CSigmaY)/(E+Hkin+Hiso);
      
	  // Find sign of xsi
      int sign = (xsi < 0) ? -1 : 1;
	  TPlasticStrain = CPlasticStrain +deltaLambda*sign;

	  TBackStress = CBackStress + Hkin*deltaLambda*sign;

	  TAccumulatedPlasticStrain = CAccumulatedPlasticStrain + deltaLambda;

	  TStress = E * (TStrain - TPlasticStrain);
	
      TTangent = E*(Hkin+Hiso) / (E+Hkin+Hiso);

	  double CSigmaYSensitivity = SigmaYSensitivity + HisoSensitivity*CAccumulatedPlasticStrain
		   + Hiso*CAccumulatedPlasticStrainSensitivity;

	  double deltaLambdaSensitivity = ((TStressSensitivity-CBackStressSensitivity)*sign-CSigmaYSensitivity)/(E+Hkin+Hiso)
		  - (ESensitivity+HkinSensitivity+HisoSensitivity)*((E * (TStrain-CPlasticStrain)-CBackStress)*sign-CSigmaY)/pow((E+Hkin+Hiso),2.0);

	  double TPlasticStrainSensitivity = CPlasticStrainSensitivity + deltaLambdaSensitivity*sign;

	  sensitivity = E * (TStrainSensitivity - TPlasticStrainSensitivity)+ESensitivity * (TStrain - TPlasticStrain);

	  double TAccumulatedPlasticStrainSensitivity = CAccumulatedPlasticStrainSensitivity + deltaLambdaSensitivity;

	  double TBackStressSensitivity = CBackStressSensitivity + HkinSensitivity*deltaLambda*sign + Hkin*deltaLambdaSensitivity*sign;

	  (*SHVs)(0,gradIndex) = TPlasticStrainSensitivity;
	  (*SHVs)(1,gradIndex) = TBackStressSensitivity;
	  (*SHVs)(2,gradIndex) = TAccumulatedPlasticStrainSensitivity;
	  (*SHVs)(3,gradIndex) = sensitivity;      // for recorder purpose
	  (*SHVs)(4,gradIndex) = TStrainSensitivity;  // for recorder purpose
    }
	*/
 	return 0;
}
// AddingSensitivity:END /////////////////////////////////////////////

