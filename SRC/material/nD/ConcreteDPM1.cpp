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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ConcreteDPM1.cpp,v $

#include <ConcreteDPM1.h>
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
#define matMode 6
#define _3dMat 6

static int numConcreteDPM1 = 0;

OPS_Export void *
OPS_NewConcreteDPM1()
{
  if (numConcreteDPM1 == 0) {
    numConcreteDPM1++;
    //opserr << "ConcreteDPM1 nDmaterial - Written: Ning Li, Tianjin University, China\n";
  }

  // Pointer to a nD material that will be returned
  NDMaterial *theMaterial = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 4) {
    opserr << "Want: NDMaterial ConcreteDPM1 tag? fc? ft? epsc0? <ecc? nu? AHard? BHard? CHard? DHard? yieldHardInitial? ASoft? helem? href?>\n";
	return 0;
  }

  int tag;
  //double rho;
  double dData[13];
  int numData = 0;

  numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid NDMaterial ConcreteDPM1 tag" << endln;
    return 0;
  }

  numData = numRemainingArgs - 1;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "Invalid Arg: NDMaterial ConcreteDPM1 tag? fc? ft? epsc0? <ecc? nu? AHard? BHard? CHard? DHard? yieldHardInitial? ASoft? helem? href?>" << endln;
    return 0;	
  }
  
  // default parameters
  double ecc = 0.525;
  double nu = 0.2;
  double yieldHardInitial = 0.1;
  double AHard = 8.e-2;
  double BHard = 3.e-3;
  double CHard = 2.;
  double DHard = 1.e-6;
  double ASoft = 15.;
  double helem = 0.;

#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
  double href = 0.;
#endif

  if (numRemainingArgs == 4) {
    theMaterial = new ConcreteDPM1(tag,
                                    dData[0],
                                    dData[1],
                                    dData[2],
                                    ecc,
                                    nu,
                                    yieldHardInitial,
                                    AHard,
                                    BHard,
                                    CHard,
                                    DHard,
                                    ASoft,
                                    helem,
                                    href);
  }
  else if (numRemainingArgs == 14) {
    theMaterial = new ConcreteDPM1(tag,
                                    dData[0],
                                    dData[1],
                                    dData[2],
                                    dData[3],
                                    dData[4],
                                    dData[5],
                                    dData[6],
                                    dData[7],
                                    dData[8],
                                    dData[9],
                                    dData[10],
                                    dData[11],
                                    dData[12]);
  }
  else {
    opserr << "OPS_NewConcreteDPM1, wrong material parameters!"<< endln;
  }


  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory creating material\n";
    opserr << "ConcreteDPM1: " << tag << endln;
    return 0;
  }

  return theMaterial;
}

ConcreteDPM1::ConcreteDPM1(int tag, 
					double fc,
					double ft,
                    double epsc0,
					double ecc,
					double nu,
                    double AHard,
                    double BHard,
                    double CHard,
                    double DHard,
                    double yieldHardInitial,
                    double ASoft,
                    double helem,
                    double href)
:NDMaterial(tag, MAT_TAG_ConcreteDPM1),
fc(fc), ft(ft), epsc0(epsc0), ecc(ecc), nu(nu), AHard(AHard), BHard(BHard),
CHard(CHard), DHard(DHard), yieldHardInitial(yieldHardInitial), ASoft(ASoft),
helem(helem), href(href)
{
	parameterID = 0;
	SHVs = 0;

	// Initialize variables
    this->revertToStart();
}

ConcreteDPM1::~ConcreteDPM1()
{
	if (SHVs != 0) 
		delete SHVs;
}

double 
ConcreteDPM1::signum(double value)
{
	if (value > 0.0) {
		return 1.0;
	}
	else {
		return -1.0;
	}
}

int 
ConcreteDPM1::setTrialStrain (const Vector &strain_from_element)
{
  // Set trial strain and compute strain increment
  Tstrain = strain_from_element;
  return 0;
}

//unused trial strain functions
int
ConcreteDPM1::setTrialStrain( const Vector &v, const Vector &r )
{ 
   return this->setTrialStrain( v ) ;
}

int
ConcreteDPM1::setTrialStrainIncr(const Vector &v) 
{
  return -1;
}

int
ConcreteDPM1::setTrialStrainIncr(const Vector &v, const Vector &r) 
{
  return -1;
}

const Vector& 
ConcreteDPM1::getStress(void)
{
  giveRealStressVector(Tstress, Tstrain);
  return Tstress;
}

const Matrix& 
ConcreteDPM1::getInitialTangent(void)
{
  give3dMaterialStiffnessMatrix(Ttangent);
  return Ttangent;
}

const Matrix& 
ConcreteDPM1::getTangent(void)
{
    return Ttangent;
}

const Vector& 
ConcreteDPM1::getStrain(void)
{
    return Tstrain;
}

int 
ConcreteDPM1::commitState(void)
{
    // Commit trial history variables
    Cstrain = Tstrain;
	Cz = Tz;
	Ce = Te;

    return 0;
}

int 
ConcreteDPM1::revertToLastCommit(void)
{
	// Nothing to do here
    return 0;
}

int 
ConcreteDPM1::revertToStart(void)
{
	Tstrain.Zero();
	Cstrain.Zero();
	Tz.Zero();
	Cz.Zero();
	Te.Zero();
	Ce.Zero();
	Tstress.Zero();
	//Ttangent *= (1-omega);
    /*Move from initializeFrom()*/
    E = eM = 1.4*fc / epsc0;
    gM = eM / (2. * (1. + nu));
    kM = eM / (3. * (1. - 2. * nu));

    double Gf;  // need to be revised
    ef = Gf / ft; // convert fracture energy to crack opening

    // default parameters
    ecc = 0.525;
    yieldHardInitial = 0.1;
    AHard = 8.e-2;
    BHard = 3.e-3;
    CHard = 2.;
    DHard = 1.e-6;
    ASoft = 15.;
    helem = 0.;
#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    href = 0.;
#endif

    //Compute m
    m = 3. * (pow(fc, 2.) - pow(ft, 2.)) / (fc * ft) * ecc / (ecc + 1.);
    //Compute default value of dilationConst
    dilationConst = -0.85;
    yieldTol = 1.e-10;
    newtonIter = 100;

	if (SHVs != 0) 
		SHVs->Zero();

    return 0;
}

NDMaterial *
ConcreteDPM1::getCopy()
{
    ConcreteDPM1 *theCopy =
      new ConcreteDPM1(this->getTag(), fc, ft, epsc0, ecc, nu, AHard, BHard, CHard, DHard, yieldHardInitial, ASoft, helem, href);
    	
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

NDMaterial *
ConcreteDPM1::getCopy(const char *type)
{
  if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) {
    ConcreteDPM1 *clone;
    clone = new ConcreteDPM1(this->getTag(), fc, ft, epsc0, ecc, nu, AHard, BHard, CHard, DHard, yieldHardInitial, ASoft, helem, href);
    return clone;
  }
  else if (strcmp(type, "PlaneStrain2D") == 0 || strcmp(type, "PlaneStrain") == 0) {
    ConcreteDPM1 *clone;
    clone = new ConcreteDPM1(this->getTag(), fc, ft, epsc0, ecc, nu, AHard, BHard, CHard, DHard, yieldHardInitial, ASoft, helem, href);
    return clone;
  }
  else {
    opserr << "ConcreteDPM1::getCopy failed to get copy: " << type << endln;
    return 0;
  }
}

//send back type of material
const char*
ConcreteDPM1::getType( ) const 
{
  return "ThreeDimensional" ;
}

int 
ConcreteDPM1::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int 
ConcreteDPM1::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}

void 
ConcreteDPM1::Print(OPS_Stream &s, int flag)
{
    s << "ConcreteDPM1, tag: " << this->getTag() << endln;
    s << "  fc: " << fc << endln;
    s << "  ft: " << ft << endln;
    s << "  epsc0: " << epsc0 << endln;
    s << "  nu: " << nu << endln;
    s << "  AHard: " << AHard << endln;
    s << "  BHard: " << BHard << endln;
    s << "  CHard: " << CHard << endln;
    s << "  DHard: " << DHard << endln;
    s << "  yieldHardInitial: " << yieldHardInitial << endln;
    s << "  ASoft: " << ASoft << endln;
    s << "  helem: " << helem << endln;
    s << "  href: " << href << endln;
}

int
ConcreteDPM1::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"fc") == 0)
    return param.addObject(1, this);
  
  if (strcmp(argv[0],"ft") == 0)
    return param.addObject(2, this);
  
  if (strcmp(argv[0],"epsc0") == 0)
    return param.addObject(3, this);
  
  if (strcmp(argv[0],"ecc") == 0)
    return param.addObject(4, this);
    
  if (strcmp(argv[0],"nu") == 0)
    return param.addObject(5, this);
  
  if (strcmp(argv[0],"AHard") == 0)
    return param.addObject(6, this);
  
  if (strcmp(argv[0],"BHard") == 0)
    return param.addObject(7, this);
    
  if (strcmp(argv[0],"CHard") == 0)
    return param.addObject(8, this);
  
  if (strcmp(argv[0],"DHard") == 0)
    return param.addObject(9, this);

  return -1;
}

int
ConcreteDPM1::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case 1:
		this->fc = info.theDouble;
		return 0;
	case 2:
		this->ft = info.theDouble;
		return 0;
	case 3:
		this->epsc0 = info.theDouble;
		return 0;
	case 4:
		this->ecc = info.theDouble;
		return 0;
	case 5:
		this->nu = info.theDouble;
		return 0;
	case 6:
		this->AHard = info.theDouble;
		return 0;
	case 7:
		this->BHard = info.theDouble;
		return 0;
	case 8:
		this->CHard = info.theDouble;
		return 0;
	case 9:
		this->DHard = info.theDouble;
		return 0;
	default:
		return -1;
	}
}

int
ConcreteDPM1::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;
	return 0;
}

const Vector & 
ConcreteDPM1::getStrainSensitivity(int gradIndex)
{
	if (SHVs ==0) return 0.0;
	else{
		double sensitivity =(*SHVs)(4,gradIndex-1); // unconditional stress sensitivity
		return sensitivity;
	}
}

const Vector & 
ConcreteDPM1::getStressSensitivity(int gradIndex, bool conditional)
{
  	return 0; //sensitivity;
}

const Matrix &
ConcreteDPM1::getInitialTangentSensitivity(int gradIndex)
{
	// For now, assume that this is only called for initial stiffness 
	return Ttangent;
}

int
ConcreteDPM1::commitSensitivity(Vector & strainGradient, int gradIndex, int numGrads)
{
 	return 0;
}

//   ******************************************************
//   *** CLASS CONCRETE DAMAGE-PLASTIC MATERIAL MODEL   ***
//   ******************************************************

//#define DPM_ITERATION_LIMIT 1.e-8

ConcreteDPM1::ConcreteDPM1() :
effectiveStress(0)
{
  tempKappaP = kappaP = 0.;
  tempKappaD = kappaD = 0.;
  tempDamage = damage = 0.;
  yieldTol = 0.;
  newtonIter = 0;
  helem = 0.;
#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
  href = 0.;
#endif
}

void
ConcreteDPM1::giveRealStressVector(Vector &answer, const Vector &strainVector)
{
  Vector reducedTotalVector = strainVector;

  //if (effectiveStress.giveStressStrainMode() == _Unknown) {
  //  effectiveStress.letStressStrainModeBe(matMode);
  //}

  //ConcreteDPMStatus *status = giveStatus(gp);

  // Initialize temp variables for this gauss point
  //initTempStatus();

  // subtract stress-independent part of strain
  // (due to temperature changes, shrinkage, etc.)
  //giveStressDependentPartOfStrainVector(reducedTotalVector, strainVector);
  Vector strain(reducedTotalVector);

  // perform plasticity return
  performPlasticityReturn(strain);

  // compute damage
  tempDamage = computeDamage(strain);

  // compute elastic strains and effective stress
  Vector elasticStrain = strain;
  //Vector tempPlasticStrain = tempPlasticStrain;
  elasticStrain.addVector(1.0, tempPlasticStrain, -1.0);
  //elasticStrain.applyElasticStiffness(effectiveStress, eM, nu);
  if (nu >= 0.5) {
    opserr <<"StrainVector::applyElasticStiffness: nu must be less than 0.5"<<endln;
  }

  double factor = eM / ((1. + nu) * (1. - 2. * nu));
  effectiveStress(0) = factor * ((1. - nu) * elasticStrain(0) + nu * elasticStrain(1) + nu * elasticStrain(2));
  effectiveStress(1) = factor * (nu * elasticStrain(0) + (1. - nu) * elasticStrain(1) + nu * elasticStrain(2));
  effectiveStress(2) = factor * (nu * elasticStrain(0) + nu * elasticStrain(1) + (1. - nu) * elasticStrain(2));
  effectiveStress(3) = factor * (((1. - 2. * nu) / 2.) * elasticStrain(3));
  effectiveStress(4) = factor * (((1. - 2. * nu) / 2.) * elasticStrain(4));
  effectiveStress(5) = factor * (((1. - 2. * nu) / 2.) * elasticStrain(5));
  // compute the nominal stress
  Vector stress(6);
  stress = effectiveStress;

  stress *= (1. - tempDamage);

  //letTempKappaDBe(tempKappaD);
  //letTempDamageBe(tempDamage);

  //letTempVectorBe(strainVector);
  //letTempVectorBe(stress);

  //assignStateFlag(gp);

  answer = stress;

#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
  // check whether the second-order work is negative
  // (must be done !!!before!!! the update of strain and stress)
  if (epsloc < 0.) {
    Vector strainIncrement, stressIncrement;
    strainIncrement = Tstrain - strainVector;
    stressIncrement = Tstress - stress;
    //int n = strainIncrement.Size();
    double work = strainIncrement^stressIncrement;
    //printf(" work : %g\n", work);
    if (work < 0.) {
      //double E = gp->giveMaterial()->give('E', gp);
      //double ft = gp->giveMaterial()->give(ft_strength, gp);
      double tmpEpsloc = kappaD + damage * ft / E;
      //letTempEpslocBe(tmpEpsloc);
    }
  }
#endif
}

double 
ConcreteDPM1::computeDamage(const Vector &strain)
{
  double equivStrain;
  computeEquivalentStrain(equivStrain, strain);
  double f = equivStrain - kappaD;
  if (f <= 0.0) {
    // damage does not grow
    tempKappaD = kappaD;
    tempDamage = damage;
  }
  else {
    // damage grows
    tempKappaD = equivStrain;
#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    if ((href <= 0.) || (epsloc >= 0.)) {
      // evaluate and store the effective element size (if not known yet)
      this->initDamaged(tempKappaD, strain);
    }

#else
    this->initDamaged(tempKappaD, strain, gp);
#endif
    tempDamage = computeDamageParam(tempKappaD);
  }

  return tempDamage;
}

void
ConcreteDPM1::computeEquivalentStrain(double &tempEquivStrain, const Vector &strain)
{
  //The equivalent  strain is based on the volumetric plastic strain
  double deltaEquivStrain = 0.;
  if (tempKappaP <= 1.0 || tempKappaP == kappaP) {
    tempEquivStrain = equivStrain;
    return;
  }
  else if (tempKappaP > 1.0 && tempKappaP != kappaP) {
    //Vector plasticStrain = givePlasticStrain();
    //Vector tempPlasticStrain = giveTempPlasticStrain();

    double volumetricPlasticStrain = plasticStrain(0) + plasticStrain(1) + plasticStrain(2);
    double tempVolumetricPlasticStrain = tempPlasticStrain(0) + tempPlasticStrain(1) + tempPlasticStrain(2);
    if (kappaP < 1.0) {
      //compute volumetric plastic strain at peak
      double peakVolumetricPlasticStrain = (1. - kappaP) / (tempKappaP - kappaP) *
        (tempVolumetricPlasticStrain - volumetricPlasticStrain) + volumetricPlasticStrain;
      if (peakVolumetricPlasticStrain < 0.) {
        peakVolumetricPlasticStrain = 0.;
      }

      deltaEquivStrain = tempVolumetricPlasticStrain - peakVolumetricPlasticStrain;
      tempEquivStrain = deltaEquivStrain / computeDuctilityMeasureDamage(strain);
      if (tempEquivStrain < 0.) {
        tempEquivStrain = 0.;
      }
    }
    else {
      deltaEquivStrain = (tempVolumetricPlasticStrain - volumetricPlasticStrain);
      if (deltaEquivStrain < 0.) {
        deltaEquivStrain = 0.;
      }

      tempEquivStrain = equivStrain +
        deltaEquivStrain / computeDuctilityMeasureDamage(strain);
    }
  }

  //letTempEquivStrainBe(tempEquivStrain);
  //letDeltaEquivStrainBe(deltaEquivStrain);
}

double
ConcreteDPM1::computeInverseDamage(double dam)
{
  if (le == 0.) {
    if (helem > 0.) {
      le = helem;
    }
    else {
      le = computeMeanSize();
    }
  }

  double answer = -(ef / le) * log(1. - dam) - dam * ft / eM;
  return answer;
}

#define DPM_DAMAGE_TOLERANCE 1.e-8

double
ConcreteDPM1::computeDamageParam(double kappa)
{
  double omega = 0.;
  if (kappa > 0.) {
    // iteration to achieve mesh objectivity
    // this is only valid for tension
    int nite = 0;
    double R, Lhs, aux, aux1;

    double h = le; // effective element size

#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    // the procedure is different before and after localization
    //double epsloc = giveEpsLoc();
    if ((href <= 0.) || (epsloc < 0.)) { // before localization, the reference element size href is used (but only if it is really specified by the user)
      if (href > 0.) {
        h = href; // reference element size
      }

#endif
      // standard damage evaluation procedure
      aux1 = (ft / eM) * h / ef;
      if (aux1 > 1) {
        opserr << "computeDamageParam: ft=" << ft << ", E = " << eM << ", wf = " << ef << ", hmax = E*wf / ft = " << eM * ef / ft << ", h = " << h << endln;
        opserr << "element too large" << endln;
      }

      do {
        nite++;
        aux = exp(-h * (omega * ft / eM + kappa) / ef);
        R = 1. - omega - aux;
        Lhs = -1. + aux1 * aux;
        omega -= R / Lhs;
        if (nite > 40) {
          opserr << "algorithm not converging" << endln;
        }
      } while (fabs(R) >= DPM_DAMAGE_TOLERANCE);

#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    }
    else {       // after localization, more complicated formula
      if (helem > 0.) {
        h = helem; // predefined element size
      }
      else {
        h = le; // effective element size
      }

      aux1 = (ft / eM) * h / ef;
      do {
        nite++;
        aux = exp(-((href - h) * epsloc + h * (omega * ft / eM + kappa)) / ef);
        R = 1. - omega - aux;
        Lhs = -1. + aux1 * aux;
        omega -= R / Lhs;
        if (nite > 40) {
          opserr << "computeDamageParam: algorithm not converging (part 2)" <<endln;
        }
      } while (fabs(R) >= DPM_DAMAGE_TOLERANCE);
    }

#endif
    if ((omega > 1.0) || (omega < 0.0)) {
      opserr << "internal error, omega = " << omega <<"." <<endln;
    }
  }

  return omega;
}

void
ConcreteDPM1::initDamaged(double kappaD, const Vector &strain)
{
  if (kappaD <= 0.) {
    return;
  }

  int i, indx = 1;
  double le;
  Vector principalStrains, crackPlaneNormal(3);
  Matrix principalDir(3, 3);

  if (helem > 0.) {
    le = helem;
  }
  else if (damage == 0.) {
    strain.computePrincipalValDir(principalStrains, principalDir);
    // find index of max positive principal strain
    for (i = 2; i <= 3; i++) {
      if (principalStrains(i) > principalStrains(indx)) {
        indx = i;
      }
    }

    for (i = 1; i <= 3; i++) {
      crackPlaneNormal(i) = principalDir(i, indx);
    }

    // Warning: This would give the element size divided by the number of Gauss points
    // le = gp->giveElement()->giveCharacteristicLenght (gp, crackPlaneNormal);

    // this gives the projected element size
    //le = giveLenghtInDir(crackPlaneNormal);
    //if (le == 0.) {
    //  le = computeMeanSize();
    //}

  }
  else if (le == 0.) {
    // this happens if the status is initialized from a file
    // with nonzero damage
    // le determined as square root of element area or cube root of el. volume
    le = computeMeanSize();
  }
}

double
ConcreteDPM1::computeDuctilityMeasureDamage(const Vector &strain)
{
  //plasticStrain = givePlasticStrain();
  // tempPlasticStrain = giveTempPlasticStrain();
  tempPlasticStrain -= plasticStrain;
  Vector principalStrain(6);
  double ductilityMeasure;
  double volStrain = tempPlasticStrain(0) + tempPlasticStrain(1) + tempPlasticStrain(2);
  //compute sum of negative principal strains
  tempPlasticStrain.computePrincipalValues(principalStrain);
  double negativeVolStrain = 0.;
  for (int i = 0; i < principalStrain.Size(); i++) {
    if (principalStrain(i) < 0.) {
      negativeVolStrain += principalStrain(i);
    }
  }

  //compute ductility measure
  double Rs = -negativeVolStrain / volStrain;
  if (Rs < 1.0) {
    ductilityMeasure = 1. + ASoft *pow(Rs, 2.);
  }
  else {
    ductilityMeasure = 1. - 3. * ASoft + 4. *ASoft *sqrt(Rs);
  }

  return ductilityMeasure;
}

void
ConcreteDPM1::performPlasticityReturn(Vector &strain)
{
  //initTempStatus();

  //get temp plastic strain and tempKappa
  //Vector tempPlasticStrain = giveTempPlasticStrain();
  //tempKappaP = giveTempKappaP();

  // compute elastic strains and trial stress
  Vector elasticStrain = strain;
  elasticStrain -= tempPlasticStrain;

  Matrix D(6, 6);
  double mu2 = eM / (1.0 + nu);
  double lam = nu*mu2 / (1.0 - 2.0*nu);
  double mu = 0.50*mu2;
  mu2 += lam;

  D(0, 0) = D(1, 1) = D(2, 2) = mu2;
  D(0, 1) = D(1, 0) = D(0, 2) = D(2, 0) = D(1, 2) = D(2, 1) = lam;
  D(3, 3) = mu;
  D(4, 4) = mu;
  D(5, 5) = mu;

  //elasticStrain..applyElasticStiffness(effectiveStress, eM, nu);
  effectiveStress.addMatrixVector(0.0, D, elasticStrain, 1.0);
  
  //Compute trial coordinates
  computeTrialCoordinates(effectiveStress);

  double yieldValue = computeYieldValue(sig, rho, thetaTrial, tempKappaP);
  // choose correct stress return and update state flag
  if (yieldValue  > yieldTol) {
    double apexStress = 0.;
    checkForVertexCase(apexStress, sig, tempKappaP);

    //Make the appropriate return
    if (vertexType == VT_Tension || vertexType == VT_Compression) {
      performVertexReturn(effectiveStress, apexStress);
      if (vertexType == VT_Regular) {
        //This was no real vertex case
        //get the original tempKappaP and stress
        //tempKappaP = giveTempKappaP();
        //elasticStrain.applyElasticStiffness(effectiveStress, eM, nu);
        effectiveStress.addMatrixVector(0.0, D, elasticStrain, 1.0);
      }
    }

    if (vertexType == VT_Regular) {
      performRegularReturn(effectiveStress);
    }
  }

  // compute the plastic strains
  //effectiveStress.applyElasticCompliance(elasticStrain, eM, nu);
  effectiveStress.addMatrixVector(0.0, D, elasticStrain, 1.0);
  tempPlasticStrain = strain;
  tempPlasticStrain -= elasticStrain;
}

bool
ConcreteDPM1::checkForVertexCase(double &answer, const double sig, const double tempKappa)
{
  //Compute sigZero and compare with actual sig
  int l = 0;
  double FZero = 1.;
  double dFZeroDSigZero = 0.;
  double sigZero;
  const double yieldHardOne = computeHardeningOne(tempKappa);
  double apexCompression = 0.;

  //compressive apex
  if ((tempKappa < 1.) && (sig < 0.)) {
    sigZero = -15. * fc;
    while (fabs(FZero) > yieldTol && l <= newtonIter) {
      l++;
      FZero = pow((1. - yieldHardOne), 2.) * pow((sigZero / fc), 4.) +
        pow(yieldHardOne, 2.) * m * (sigZero / fc) - pow(yieldHardOne, 2.);

      dFZeroDSigZero = pow((1. - yieldHardOne), 2.) * 4. * pow((sigZero / fc), 3.) / fc +
        pow(yieldHardOne, 2.) * m / fc;

      sigZero = sigZero - FZero / dFZeroDSigZero;
    }

    if (l < 15 && sigZero < 0.) {
      apexCompression = sigZero;
    }
    else {
      apexCompression = -15. * fc;
    }
  }

  if ((sig > 0. && tempKappa < 1.) || (sig > fc / m && tempKappa >= 1.)) {
    vertexType = VT_Tension;
    answer = 0.;
  }
  else if (sig < apexCompression) {
    vertexType = VT_Compression;
    answer = apexCompression;
  }
  else {
    vertexType = VT_Regular;
  }

  return false;
}

void
ConcreteDPM1::performRegularReturn(Vector &effectiveStress)
{
  //Define status
  //ConcreteDPMStatus *status = static_cast< ConcreteDPMStatus * >(this->giveStatus(gp));
  //MaterialMode matMode = gp->giveMaterialMode();
  //Variables
  deltaLambda = 0.;
  double deltaLambdaIncrement = 0.;
  double yieldValue;
  double residualNorm;
  double rhoTest = 0.;
  double deltaLambdaIncrementNew = 0.;
  double tempKappaPTest = 0.;
  int iterationCount = 0;
  int negativeRhoFlag = 0;
  Vector elasticStrain(matMode);
  Vector dGDInv(2);
  Vector dFDInv(2);
  Vector residual(3);
  Vector derivativesOfYieldSurface(3);
  Vector flowRules(3);
  Vector answerIncrement(3);
  //  double residual2;
  double dKappaDDeltaLambda;
  double dFDKappa = 0.;
  Matrix aMatrix;
  Vector helpVectorA(3);
  Vector helpVectorB(3);
  Vector helpIncrement(matMode);
  Vector plasticStrainIncrement(matMode);
  Vector dGDStressDeviatoric(matMode);
  Vector deviatoricStress(matMode);

  //compute the principal directions of the stress
  Vector help;
  Matrix stressPrincipalDir;
  effectiveStress.computePrincipalValDir(help, stressPrincipalDir);

  //compute invariants from stress state
  effectiveStress.computeDeviatoricVolumetricSplit(deviatoricStress, sig);
  rho = deviatoricStress.computeSecondCoordinate();

  const Vector initialStress = effectiveStress;

  const double volumetricPlasticStrain = 3. * giveVolumetricPlasticStrain();

  const double deviatoricPlasticStrainNorm = giveDeviatoricPlasticStrainNorm();

  double tempVolumetricPlasticStrain = volumetricPlasticStrain;
  double tempDeviatoricPlasticStrainNorm = deviatoricPlasticStrainNorm;

  //const double kappaP = giveKappaP();
  tempKappaP = kappaP;

  yieldValue = computeYieldValue(sig, rho, thetaTrial, tempKappaP);
  residualNorm = fabs(yieldValue);

  while (residualNorm > yieldTol) {
    if (++iterationCount == newtonIter) {
      opserr << "Closest point projection did not converge.\n";
    }

    //compute the stress, yield value and residuals
    computeDGDInv(dGDInv, sig, rho, tempKappaP);
    computeDFDInv(dFDInv, sig, rho, tempKappaP);
    dFDKappa = computeDFDKappa(sig, rho, tempKappaP);
    dKappaDDeltaLambda = computeDKappaDDeltaLambda(sig, rho, tempKappaP);
    yieldValue = computeYieldValue(sig, rho, thetaTrial, tempKappaP);
    //The residual volumetric plastic strain is defined as eps1+eps2+eps3
    residual(0) = volumetricPlasticStrain - tempVolumetricPlasticStrain + deltaLambda *dGDInv(0);
    residual(1) = deviatoricPlasticStrainNorm - tempDeviatoricPlasticStrainNorm + deltaLambda *dGDInv(1);
    residual(2) = kappaP - tempKappaP + deltaLambda * dKappaDDeltaLambda;

    // weighted norm
    residualNorm = pow(residual(0), 2.) + pow(residual(1), 2.) + pow(residual(2), 2.) + pow(yieldValue, 2.);
    residualNorm = sqrt(residualNorm);
    //    printf("\n residualNorm= %e\n", residualNorm);
    if (residualNorm > yieldTol) {
      computeAMatrix(aMatrix, sig, rho, tempKappaP);

      // assemble the derivatives of the yield surface
      for (int i = 0; i < 2; i++) {
        derivativesOfYieldSurface(i) = dFDInv(i);
      }

      derivativesOfYieldSurface(2) = dFDKappa;
      //assemble flow rules
      for (int i = 0; i < 2; i++) {
        flowRules(i) = dGDInv(i);
      }

      flowRules(2) = dKappaDDeltaLambda;

      //compute the deltaLambdaIncrement
      deltaLambdaIncrement = yieldValue;
      helpVectorA.addMatrixVector(0.0, aMatrix, residual, 1.0);
      for (int i = 0; i < 3; i++) {
        deltaLambdaIncrement -= derivativesOfYieldSurface(i) * helpVectorA(i);
      }

      helpVectorB.addMatrixVector(0.0, aMatrix, flowRules, 1.0);
      double denominator = 0.;
      for (int i = 0; i < 3; i++) {
        denominator += derivativesOfYieldSurface(i) * helpVectorB(i);
      }

      deltaLambdaIncrement /= denominator;

      //compute increment of answer
      Vector helpVectorC;
      helpVectorC = residual;
      helpVectorC.addVector(1.0, flowRules, deltaLambdaIncrement);
      answerIncrement.addMatrixVector(1.0, aMatrix, helpVectorC, 1.0);
      answerIncrement *= (-1.0);
      rhoTest = rho + answerIncrement(1);

      //Special case, if rho changes sign
      if (rhoTest < 0. && negativeRhoFlag == 0) {
        //Determine deltaLambdaIncrement, so that rho is equal to zero
        answerIncrement(1) = -rho;

        deltaLambdaIncrementNew =
          (-aMatrix(1, 0) * residual(0) - aMatrix(1, 1) * residual(1) -
          aMatrix(1, 2) * residual(2) - answerIncrement(1)) /
          (flowRules(0) * aMatrix(1, 0) + flowRules(1) * aMatrix(1, 1) +
          flowRules(2) * aMatrix(1, 2));

        // Special case, if deltaLambdaIncrement is equal to zero.
        if (fabs(deltaLambdaIncrementNew) < yieldTol * 1.e3) {
          negativeRhoFlag = 1;
          deltaLambdaIncrementNew = deltaLambdaIncrement;
        }

        helpVectorC = residual;
        helpVectorC.addVector(1.0, flowRules, deltaLambdaIncrementNew);
        answerIncrement.addMatrixVector(1.0, aMatrix, helpVectorC, 1.0);
        answerIncrement *= (-1.);

        sig += answerIncrement(0);
        rho += answerIncrement(1);

        tempKappaPTest = tempKappaP;
        tempKappaP += answerIncrement(2);
        deltaLambda += deltaLambdaIncrementNew;
      }
      else {
        tempKappaPTest = tempKappaP;
        tempKappaP += answerIncrement(2);
        sig += answerIncrement(0);
        rho += answerIncrement(1);

        deltaLambda += deltaLambdaIncrement;
      }

      //Special case, if deltaKappaP is negative
      if ((tempKappaP - kappaP) < 0.) {
        tempKappaP = tempKappaPTest;
      }

      tempVolumetricPlasticStrain -= answerIncrement(0) / (kM);
      tempDeviatoricPlasticStrainNorm -= answerIncrement(1) / (2. * gM);
    }

    //if (gp->giveElement()->giveNumber() == 1301){
    //  printf("%g %g %g\n",sig,rho,tempKappaP);
    //}
  }

  tempVolumetricPlasticStrain *= (1./3.);
  if (deltaLambda < 0) {
    opserr << "deltaLambda = " << deltaLambda << "." << endln;
    opserr << "plastic multiplier less than zero" << endln;
  }

  // printf("\nnumber of iterations = %i\n", iterationCount);
  Vector stressPrincipal(3);
  Vector stressTemp(6);
  stressPrincipal.Zero();

  stressPrincipal(0) = sig + sqrt(2. / 3.) * rho * cos(thetaTrial);
  stressPrincipal(1) = sig + sqrt(2. / 3.) * rho * cos(thetaTrial - 2. * PI / 3.);
  stressPrincipal(2) = sig + sqrt(2. / 3.) * rho * cos(thetaTrial + 2. * PI / 3.);

  //transformVectorTo(stressTemp, stressPrincipalDir, stressPrincipal, 1);
  stressTemp.addMatrixTransposeVector(0.0, stressPrincipalDir, stressPrincipal, 1.0);

  if (matMode == 4) { //_PlaneStrain
    effectiveStress(0) = stressTemp(0);
    effectiveStress(1) = stressTemp(1);
    effectiveStress(2) = stressTemp(2);
    effectiveStress(3) = stressTemp(5);
  }
  else {
    effectiveStress = stressTemp;
  }
}

void
ConcreteDPM1::performVertexReturn(Vector &effectiveStress, double apexStress)
{
  Vector deviatoricStressTrial(matMode);
  double sigTrial;
  Vector stressTemp(matMode);
  Vector deviatoricStress(matMode);
  double yieldValue = 0.;
  double yieldValueMid = 0.;
  double sig2 = 0.;
  double dSig;
  double sigMid;
  double sigAnswer;
  double ratioPotential;

  effectiveStress.computeDeviatoricVolumetricSplit(deviatoricStressTrial, sigTrial);
  const double rhoTrial = deviatoricStressTrial.computeSecondCoordinate();

  Vector deltaPlasticStrain(matMode);

  //tempKappaP = giveTempKappaP();
  const double kappaInitial = tempKappaP;

  if (vertexType == VT_Tension) {
    sig2 = -0.1 * ft;
  }
  else if (vertexType == VT_Compression) {
    sig2 = apexStress;
  }

  tempKappaP = computeTempKappa(kappaInitial, sigTrial, rhoTrial, sigTrial);

  yieldValue = computeYieldValue(sigTrial, 0., 0., tempKappaP);

  tempKappaP = computeTempKappa(kappaInitial, sigTrial, rhoTrial, sig2);

  yieldValueMid = computeYieldValue(sig2, 0., 0., tempKappaP);

  if (yieldValue * yieldValueMid >= 0.) {
    //vertexType = VT_Regular;
    return;
  }

  if (yieldValue < 0.0) {
    dSig = sig2 - sigTrial;
    sigAnswer = sig2;
  }
  else {
    dSig = sigTrial - sig2;
    sigAnswer = sig2;
  }

  for (int j = 0; j < 250; j++) {
    dSig = 0.5 * dSig;

    sigMid = sigAnswer + dSig;
    
    tempKappaP = computeTempKappa(kappaInitial, sigTrial, rhoTrial, sigMid);

    yieldValueMid = computeYieldValue(sigMid, 0., 0., tempKappaP);

    if (yieldValueMid <= 0.) {
      sigAnswer = sigMid;
    }

    if (fabs(yieldValueMid) < yieldTol && yieldValueMid <= 0.) {
      for (int i = 0; i < 3; i++) {
        effectiveStress(i) = sigAnswer;
      }

      for (int i = 3; i < effectiveStress.Size(); i++) {
        effectiveStress(i) = 0.;
      }

      ratioPotential = computeRatioPotential(sigAnswer, tempKappaP);

      double ratioTrial = rhoTrial / (sigTrial - sigAnswer);

      if ((((ratioPotential >= ratioTrial) && vertexType == VT_Tension)) ||
        ((ratioPotential <= ratioTrial) && vertexType == VT_Compression)) {
        return;
      }
      else {
        //vertexType = VT_Regular;
        return;
      }
    }
  }

  vertexType = VT_Regular;
}

double
ConcreteDPM1::computeTempKappa(const double kappaInitial, const double sigTrial, const double rhoTrial, const double sig)
{
  //This function is called, if stress state is in vertex case
  double equivalentDeltaPlasticStrain;
  Vector deltaPlasticStrainPrincipal(3);
  rho = 0.;
  equivalentDeltaPlasticStrain = sqrt(1. / 9. * pow((sigTrial - sig) / (kM), 2.) +
    pow(rhoTrial / (2. * gM), 2.));

  double thetaVertex = 0.;
  double ductilityMeasure = computeDuctilityMeasure(sig, rho, thetaVertex);

  return kappaInitial + equivalentDeltaPlasticStrain / ductilityMeasure;
}

double
ConcreteDPM1::computeYieldValue(const double sig, const double rho, const double theta, const double tempKappa) const
{
  //compute yieldHard
  const double yieldHardOne = computeHardeningOne(tempKappa);
  const double yieldHardTwo = computeHardeningOne(tempKappa);

  //  compute elliptic function r
  const double rFunction = (4. * (1. - pow(ecc, 2.)) * pow(cos(theta), 2.) +
    pow((2. * ecc - 1.), 2.)) /
    (2. * (1. - pow(ecc, 2.)) * cos(theta) +
    (2. * ecc - 1.) * sqrt(4. * (1. - pow(ecc, 2.)) * pow(cos(theta), 2.)
    + 5. * pow(ecc, 2.) - 4. * ecc));

  //compute help function Al
  const double Al = (1. - yieldHardOne) * pow((sig / fc + rho / (sqrt(6.) * fc)), 2.) +
    sqrt(3. / 2.) * rho / fc;

  //Compute yield equation
  return pow(Al, 2.) +
    pow(yieldHardOne, 2.) * m * (sig / fc + rho * rFunction / (sqrt(6.) * fc)) -
    pow(yieldHardTwo, 2.);
}

double
ConcreteDPM1::computeDFDKappa(const double sig, const double rho, const double tempKappa)
{
  const double theta = thetaTrial;
  //compute yieldHard and yieldSoft
  const double yieldHardOne = computeHardeningOne(tempKappa);
  const double yieldHardTwo = computeHardeningOne(tempKappa);
  // compute the derivative of the hardening and softening laws
  const double dYieldHardOneDKappa = computeHardeningOnePrime(tempKappa);
  const double dYieldHardTwoDKappa = computeHardeningOnePrime(tempKappa);

  //compute elliptic function r
  const double rFunction =
    (4. * (1. - ecc * ecc) * cos(theta) * cos(theta) + (2. * ecc - 1.) * (2. * ecc - 1.)) /
    (2 * (1. - ecc * ecc) * cos(theta) + (2. * ecc - 1.) *
    sqrt(4. * (1. - ecc * ecc) * cos(theta) * cos(theta) + 5. * ecc * ecc - 4. * ecc));

  //compute help functions Al, Bl
  const double Al = (1. - yieldHardOne) * pow((sig / fc + rho / (sqrt(6.) * fc)), 2.) + sqrt(3. / 2.) * rho / fc;

  const double Bl = sig / fc + rho / (fc * sqrt(6.));

  const double dFDYieldHardOne = -2. *Al *pow(Bl, 2.)
    + 2. * yieldHardOne * m * (sig / fc + rho * rFunction / (sqrt(6.) * fc));

  const double dFDYieldHardTwo = -2. * yieldHardTwo;

  // compute dFDKappa
  double dFDKappa = dFDYieldHardOne * dYieldHardOneDKappa +
    dFDYieldHardTwo * dYieldHardTwoDKappa;
  /*
  * set dFDKappa to zero, if it becomes greater than zero.
  * dFDKappa can only be negative or zero in the converged state for
  * the case of hardenig and perfect plasticity. For trial stresses, however,
  * it might be negative, which may lead to convergency problems. Therefore,
  * it is set to zero in this cases.
  */
  if (dFDKappa > 0.) {
    dFDKappa = 0.;
  }

  return dFDKappa;
}

void
ConcreteDPM1::computeDFDInv(Vector &answer, const double sig, const double rho, const double tempKappa) const
{
  const double theta = thetaTrial;

  //compute yieldHard
  const double yieldHardOne = computeHardeningOne(tempKappa);

  //compute elliptic function r
  const double rFunction = (4. * (1. - ecc * ecc) * cos(theta) * cos(theta) + (2. * ecc - 1.) * (2. * ecc - 1.)) /
    (2. * (1. - ecc * ecc) * cos(theta) + (2. * ecc - 1.) * sqrt(4. * (1. - ecc * ecc) * cos(theta) * cos(theta)
    + 5. * ecc * ecc - 4. * ecc));

  //compute help functions AL, BL
  const double AL = (1. - yieldHardOne) * pow((sig / fc + rho / (sqrt(6.) * fc)), 2.) + sqrt(3. / 2.) * rho / fc;
  const double BL = sig / fc + rho / (fc * sqrt(6.));

  //compute dfdsig
  const double dfdsig = 4. * (1. - yieldHardOne) / fc * AL * BL + yieldHardOne * yieldHardOne * m / fc;
  //compute dfdrho
  const double dfdrho = AL / (sqrt(6.) * fc) * (4. * (1. - yieldHardOne) * BL + 6.) + rFunction * m * yieldHardOne * yieldHardOne / (sqrt(6.) * fc);

  answer(0) = dfdsig;
  answer(1) = dfdrho;
}

double
ConcreteDPM1::computeDKappaDDeltaLambda(const double sig, const double rho, const double tempKappa)
{
  //Variables
  double equivalentDGDStress;
  Vector dGDInv(2);
  computeDGDInv(dGDInv, sig, rho, tempKappa);
  Vector dGDStressPrincipal(3);

  equivalentDGDStress = sqrt(1. / 3. * pow(dGDInv(0), 2.) + pow(dGDInv(1), 2.));

  double ductilityMeasure = computeDuctilityMeasure(sig, rho, this->thetaTrial);
  double dKappaDDeltaLambda = equivalentDGDStress / ductilityMeasure;
  return dKappaDDeltaLambda;
}

void
ConcreteDPM1::computeDDKappaDDeltaLambdaDInv(Vector &answer, const double sig, const double rho, const double tempKappa)
{
  //Variables
  double equivalentDGDStress;
  Vector dGDInv(2);
  Matrix dDGDDInv(2, 2);
  Vector dGDStressPrincipal(3);
  Vector dEquivalentDGDStressDInv(2);
  Vector helpA(3);

  //Compute first and second derivative of plastic potential
  computeDGDInv(dGDInv, sig, rho, tempKappa);
  computeDDGDDInv(dDGDDInv, sig, rho, tempKappa);

  //Compute equivalentDGDStress
  equivalentDGDStress = sqrt(1. / 3. * pow(dGDInv(0), 2.) + pow(dGDInv(1), 2.));

  //computeDuctilityMeasure
  double ductilityMeasure = computeDuctilityMeasure(sig, rho, this->thetaTrial);

  //Compute dEquivalentDGDStressDInv
  dEquivalentDGDStressDInv(0) =
    (2. / 3. * dGDInv(0) * dDGDDInv(0, 0) + 2. * dGDInv(1) * dDGDDInv(1, 0)) / (2. * equivalentDGDStress);
  dEquivalentDGDStressDInv(1) =
    (2. / 3. * dGDInv(0) * dDGDDInv(0, 1) + 2. * dGDInv(1) * dDGDDInv(1, 1)) / (2. * equivalentDGDStress);

  answer.Zero();

  // compute the derivative of
  Vector dDuctilityMeasureDInv(2);
  computeDDuctilityMeasureDInv(dDuctilityMeasureDInv, sig, rho, tempKappa);

  answer(0) = (dEquivalentDGDStressDInv(0) * ductilityMeasure - equivalentDGDStress * dDuctilityMeasureDInv(0)) / pow(ductilityMeasure, 2.);

  answer(1) = (dEquivalentDGDStressDInv(1) * ductilityMeasure - equivalentDGDStress * dDuctilityMeasureDInv(1)) / pow(ductilityMeasure, 2.);
}

double
ConcreteDPM1::computeDDKappaDDeltaLambdaDKappa(const double sig, const double rho, const double tempKappa)
{
  //Variables
  double equivalentDGDStress;
  Vector dGDInv(2);
  Vector dDGDInvDKappa(2);
  Vector dGDStressPrincipal(3);
  Vector helpA(3);
  double dEquivalentDGDStressDKappa;

  //Compute first and second derivative of plastic potential
  computeDGDInv(dGDInv, sig, rho, tempKappa);
  computeDDGDInvDKappa(dDGDInvDKappa, sig, rho, tempKappa);

  equivalentDGDStress = sqrt(1. / 3. * pow(dGDInv(0), 2.) + pow(dGDInv(1), 2.));

  //computeDuctilityMeasure
  double ductilityMeasure = computeDuctilityMeasure(sig, rho, this->thetaTrial);

  //Compute dEquivalentDGDStressDKappa
  dEquivalentDGDStressDKappa = (2. / 3. * dGDInv(0) * dDGDInvDKappa(0) + 2. * dGDInv(1) * dDGDInvDKappa(1)) / (2. * equivalentDGDStress);

  // compute the derivative of
  double dDuctilityMeasureDKappa = 0.;

  double dDKappaDDeltaLambdaDKappa = (dEquivalentDGDStressDKappa * ductilityMeasure -
    equivalentDGDStress * dDuctilityMeasureDKappa) / pow(ductilityMeasure, 2.);

  return dDKappaDDeltaLambdaDKappa;
}

double
ConcreteDPM1::computeDuctilityMeasure(const double sig, const double rho, const double theta)
{
  double thetaConst = pow(2. * cos(theta), 2.);
  double ductilityMeasure;
  double x = -(sig + fc / 3) / fc;
  if (x < 0.) {
    /*Introduce exponential help function which results in a smooth
    * transition. */
    double fZero = BHard;
    double fPrimeZero = -(BHard - AHard) / (CHard);
    double CHelp = DHard;
    double AHelp = fZero - CHelp;
    double BHelp = (fZero - CHelp) / fPrimeZero;
    ductilityMeasure = (AHelp * exp(x / BHelp) + CHelp) / thetaConst;
  }
  else {
    ductilityMeasure = (AHard + (BHard - AHard) * exp(-x / (CHard))) / thetaConst;
  }

  if (ductilityMeasure <= 0.) {
    opserr << "ductilityMeasure is zero or negative" << endln;
  }

  return ductilityMeasure;
}

void
ConcreteDPM1::computeDDuctilityMeasureDInv(Vector &answer, const double sig, const double rho, const double tempKappa)
{
  double thetaConst = pow(2. * cos(thetaTrial), 2.);
  double x = (-(sig + fc / 3)) / fc;
  if (x < 0.) {
    double dXDSig = -1. / fc;
    /* Introduce exponential help function which results in a
    * smooth transition. */
    double fZero = BHard;
    double fPrimeZero = -(BHard - AHard) / (CHard);
    double CHelp = DHard;
    double AHelp = fZero - CHelp;
    double BHelp = (fZero - CHelp) / fPrimeZero;
    double dDuctilityMeasureDX = AHelp / BHelp *exp(x / BHelp) / thetaConst;
    answer(0) = dDuctilityMeasureDX * dXDSig;
    answer(1) = 0.;
  }
  else {
    double dXDSig = -1. / fc;
    double dDuctilityMeasureDX = -(BHard - AHard) / (CHard) / thetaConst *exp(-x / (CHard));
    answer(0) = dDuctilityMeasureDX * dXDSig;
    answer(1) = 0.;
  }
}

void
ConcreteDPM1::computeDGDInv(Vector &answer, const double sig, const double rho, const double tempKappa)
{
  //Compute dilation parameter
  const double R = (sig - ft / 3.) / fc;
  const double mQTension = (3. * ft / fc + m / 2.);
  const double mQCompression = (1. + 2. * dilationConst) / (dilationConst - 1.) * (3. + m / 2.);
  const double AConst = -(ft + fc) / (3. * fc) / log(mQCompression / mQTension);
  const double mQ = mQTension * exp(R / AConst);

  //compute yieldHard and yieldSoft
  const double yieldHardOne = computeHardeningOne(tempKappa);

  const double Bl = sig / fc + rho / (fc * sqrt(6.));

  const double Al = (1. - yieldHardOne) * pow(Bl, 2.) + sqrt(3. / 2.) * rho / fc;

  const double dgdsig = 4. * (1. - yieldHardOne) / fc * Al * Bl + yieldHardOne * yieldHardOne * mQ / fc;

  const double dgdrho = Al / (sqrt(6.) * fc) * (4. * (1. - yieldHardOne) * Bl + 6.) +
    m *pow(yieldHardOne, 2.) / (sqrt(6.) * fc);

  answer(0) = dgdsig;
  answer(1) = dgdrho;
}

double
ConcreteDPM1::computeRatioPotential(const double sig, const double tempKappa)
{
  //Compute dilation parameter
  const double R = (sig - ft / 3.) / fc;
  const double mQTension = (3. * ft / fc + m / 2.);
  const double mQCompression =
    (1. + 2. * dilationConst) / (dilationConst - 1.) * (3. + m / 2.);
  const double AConst = -(ft + fc) / (3. * fc) / log(mQCompression / mQTension);
  const double mQ = mQTension * exp(R / AConst);

  //compute yieldHard and yieldSoft
  const double yieldHardOne = computeHardeningOne(tempKappa);

  const double Bl = sig / fc;

  const double Al = (1. - yieldHardOne) * pow(Bl, 2.);

  const double dgdsig = 4. * (1. - yieldHardOne) / fc * Al * Bl + yieldHardOne * yieldHardOne * mQ / fc;

  const double dgdrho = Al / (sqrt(6.) * fc) * (4. * (1. - yieldHardOne) * Bl + 6.) +
    m *pow(yieldHardOne, 2.) / (sqrt(6.) * fc);

  /* Debug: Introduce the factor 2G/K. I am not sure if this is correct.
  * There might be too many stress states in the vertex case. */
  return dgdrho / dgdsig * 3. * (1. - 2. * nu) / (1. + nu);
}

void
ConcreteDPM1::computeDDGDInvDKappa(Vector &answer, const double sig, const double rho, const double tempKappa)
{
  //Compute dilation parameter
  const double R = (sig - ft / 3.) / fc;
  const double mQTension = (3. * ft / fc + m / 2.);
  const double mQCompression = (1. + 2. * dilationConst) / (dilationConst - 1.) * (3. + m / 2.);
  const double AConst = -(ft + fc) / (3. * fc) / log(mQCompression / mQTension);
  const double mQ = mQTension * exp(R / AConst);

  //compute yieldHard and yieldSoft
  const double yieldHardOne = computeHardeningOne(tempKappa);
  const double dYieldHardOneDKappa = computeHardeningOnePrime(tempKappa);

  const double Bl = sig / fc + rho / (fc * sqrt(6.));

  const double Al = (1. - yieldHardOne) * pow(Bl, 2.) + sqrt(3. / 2.) * rho / fc;

  const double dAlDYieldHard = -pow(Bl, 2.);

  const double dDGDSigDKappa =
    (-4. * Al * Bl / fc + 4. * (1 - yieldHardOne) / fc * dAlDYieldHard * Bl +
    2. * yieldHardOne * mQ / fc) * dYieldHardOneDKappa;

  const double dDGDRhoDKappa =
    (dAlDYieldHard / (sqrt(6.) * fc) * (4. * (1. - yieldHardOne) * Bl + 6.) -
    4. * Al / (sqrt(6.) * fc) * Bl + 2. * m * yieldHardOne / (sqrt(6.) * fc)) * dYieldHardOneDKappa;

  answer(0) = dDGDSigDKappa;
  answer(1) = dDGDRhoDKappa;
}

void
ConcreteDPM1::computeDDGDDInv(Matrix &answer, const double sig, const double rho, const double tempKappa)
{
  //Compute dilation parameter
  const double R = (sig - ft / 3.) / fc;
  const double mQTension = (3. * ft / fc + m / 2.);
  const double mQCompression =
    (1. + 2. * dilationConst) / (dilationConst - 1.) * (3. + m / 2.);
  const double AConst = -(ft + fc) / (3. * fc) / log(mQCompression / mQTension);
  //  const double mQ = mQTension*exp(R/AConst);
  const double dMQDSig = mQTension / (AConst * fc) * exp(R / AConst);

  //compute yieldHardOne and yieldSoft
  const double yieldHardOne = computeHardeningOne(tempKappa);

  //compute help parameter Al and Bl and the corresponding derivatives
  const double Bl = sig / fc + rho / (fc * sqrt(6.));

  const double Al = (1. - yieldHardOne) * pow(Bl, 2.) + sqrt(3. / 2.) * rho / fc;

  const double dAlDSig = 2. * (1. - yieldHardOne) * Bl / fc;
  const double dBlDSig = 1. / fc;

  const double dAlDRho = 2. * (1. - yieldHardOne) * Bl / (fc * sqrt(6.)) + sqrt(3. / 2.) / fc;
  const double dBlDRho = 1. / (fc * sqrt(6.));

  //compute second derivatives of g
  const double ddgddSig = 4. * (1. - yieldHardOne) / fc * (dAlDSig * Bl + Al * dBlDSig) +
    yieldHardOne * yieldHardOne * dMQDSig / fc;

  const double ddgddRho = dAlDRho / (sqrt(6.) * fc) * (4. * (1. - yieldHardOne) * Bl + 6.) +
    Al * dBlDRho * 4. * (1. - yieldHardOne) / (sqrt(6.) * fc);

  const double ddgdSigdRho = 4. * (1. - yieldHardOne) / fc * (dAlDRho * Bl + Al * dBlDRho);

  const double ddgdRhodSig = dAlDSig / (sqrt(6.) * fc) * (4. * (1. - yieldHardOne) * Bl + 6.)
    + Al / (sqrt(6.) * fc) * (4. * (1. - yieldHardOne) * dBlDSig);

  answer(0, 0) = ddgddSig;
  answer(0, 1) = ddgdSigdRho;
  answer(1, 0) = ddgdRhodSig;
  answer(1, 1) = ddgddRho;
}

void
ConcreteDPM1::computeAMatrix(Matrix &answer, const double sig, const double rho, const double tempKappa)
{
  Matrix aMatrixInverse(3, 3);
  Vector dDKappaDDeltaLambdaDInv(2);
  Matrix dDGDDInv(2, 2);
  Vector dDGDInvDKappa(2);
  double dDKappaDDeltaLambdaDKappa;
  answer.Zero();
  computeDDGDDInv(dDGDDInv, sig, rho, tempKappa);
  computeDDKappaDDeltaLambdaDInv(dDKappaDDeltaLambdaDInv, sig, rho,
    tempKappa);
  dDKappaDDeltaLambdaDKappa =
    computeDDKappaDDeltaLambdaDKappa(sig, rho, tempKappa);
  computeDDGDInvDKappa(dDGDInvDKappa, sig, rho, tempKappa);


  aMatrixInverse(0, 0) = 1. / (kM)+deltaLambda *dDGDDInv(0, 0);
  aMatrixInverse(0, 1) = deltaLambda * dDGDDInv(0, 1);
  aMatrixInverse(0, 2) = deltaLambda * dDGDInvDKappa(0);

  aMatrixInverse(1, 0) = deltaLambda * dDGDDInv(1, 0);
  aMatrixInverse(1, 1) = 1. / (2. * gM) + deltaLambda *dDGDDInv(1, 1);
  aMatrixInverse(1, 2) = deltaLambda * dDGDInvDKappa(1);

  aMatrixInverse(2, 0) = deltaLambda * dDKappaDDeltaLambdaDInv(0);
  aMatrixInverse(2, 1) = deltaLambda * dDKappaDDeltaLambdaDInv(1);
  aMatrixInverse(2, 2) = -1. + deltaLambda * dDKappaDDeltaLambdaDKappa;

  aMatrixInverse.Invert(answer);
}

double
ConcreteDPM1::computeHardeningOne(const double kappa) const
{
  if (kappa <= 0.) {
    return yieldHardInitial;
  }
  else if (kappa > 0. && kappa < 1.) {
    return yieldHardInitial + (1. - yieldHardInitial) *
      kappa * (pow(kappa, 2.) - 3. * kappa + 3.);
  }
  else {
    return 1.;
  }
}

double
ConcreteDPM1::computeHardeningOnePrime(const double kappa) const
{
  if (kappa < 0.) {
    return 3. * (1. - yieldHardInitial);
  }
  else if (kappa >= 0. && kappa < 1.) {
    return (1. - yieldHardInitial) * (3. * pow(kappa, 2.) - 6. * kappa + 3.);
  }
  else {
    return 0.;
  }
}

void
ConcreteDPM1::give3dMaterialStiffnessMatrix(Matrix &answer)
{
  double omega = tempDamage;
  if (omega > 0.9999) {
    omega = 0.9999;
  }
  answer.resize(6, 6);
  double mu2 = eM / (1.0 + nu);
  double lam = nu*mu2 / (1.0 - 2.0*nu);
  double mu = 0.50*mu2;
  mu2 += lam;

  answer(0, 0) = answer(1, 1) = answer(2, 2) = mu2;
  answer(0, 1) = answer(1, 0) = answer(0, 2) = answer(2, 0) = answer(1, 2) = answer(2, 1) = lam;
  answer(3, 3) = mu;
  answer(4, 4) = mu;
  answer(5, 5) = mu;

  //this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, mode);
  answer *= (1. - omega);
  
}

void
ConcreteDPM1::computeTrialCoordinates(const Vector &stress)
{
  Vector deviatoricStress(effectiveStress.Size());
  effectiveStress.computeDeviatoricVolumetricSplit(deviatoricStress,
    sig);
  rho = deviatoricStress.computeSecondCoordinate();
  thetaTrial = deviatoricStress.computeThirdCoordinate();
}

void
ConcreteDPM1::assignStateFlag()
{
  //Get kappaD from status to define state later on
  //damage = giveDamage();
  //tempDamage = giveTempDamage();
  //kappaP = giveKappaP();
  //tempKappaP = giveTempKappaP();
  /*
  if (tempKappaP > kappaP) {
    if (tempDamage > damage) {
      //letTempStateFlagBe(ConcreteDPMStatus::ConcreteDPM_PlasticDamage);
    }
    else {
      //letTempStateFlagBe(ConcreteDPMStatus::ConcreteDPM_Plastic);
    }
  }
  else {
    const int state_flag = giveStateFlag();
    if (state_flag == ConcreteDPMStatus::ConcreteDPM_Elastic) {
      if (tempDamage > damage) {
        //letTempStateFlagBe(ConcreteDPMStatus::ConcreteDPM_Damage);
      }
      else {
        //letTempStateFlagBe(ConcreteDPMStatus::ConcreteDPM_Elastic);
      }
    }
    else {
      if (tempDamage > damage) {
        //letTempStateFlagBe(ConcreteDPMStatus::ConcreteDPM_Damage);
      }
      else {
        //letTempStateFlagBe(ConcreteDPMStatus::ConcreteDPM_Unloading);
      }
    }
  }*/
}

void
ConcreteDPM1::computeDRhoDStress(Vector &answer, const Vector &stress) const
{
  int size = 6;
  //compute volumetric deviatoric split
  Vector deviatoricStress(6);
  double volumetricStress;
  stress.computeDeviatoricVolumetricSplit(deviatoricStress, volumetricStress);
  double rho = deviatoricStress.computeSecondCoordinate();

  //compute the derivative of J2 with respect to the stress
  Vector dJ2DStress;
  dJ2DStress = deviatoricStress;
  for (int i = 3; i < size; i++) {
    dJ2DStress(i) = deviatoricStress(i) * 2.0;
  }

  //compute the derivative of rho with respect to stress
  Vector dRhoDStress;
  dRhoDStress = dJ2DStress;
  dRhoDStress *= (1. / rho);

  answer = dRhoDStress;
}

void
ConcreteDPM1::computeDSigDStress(Vector &answer) const
{
  int size = 6;
  for (int i = 0; i < 3; i++) {
    answer(i) = 1. / 3.;
  }

  for (int i = 3; i < size; i++) {
    answer(i) = 0.;
  }
}

void
ConcreteDPM1::computeDDRhoDDStress(Matrix &answer, const Vector &stress) const
{
  int size = 6;

  //compute volumetric deviatoric split
  Vector deviatoricStress(_3dMat);
  double volumetricStress;
  stress.computeDeviatoricVolumetricSplit(deviatoricStress, volumetricStress);
  double rho = deviatoricStress.computeSecondCoordinate();


  //compute first derivative of J2
  Vector dJ2dstress;
  dJ2dstress = deviatoricStress;
  for (int i = 3; i < deviatoricStress.Size(); i++) {
    dJ2dstress(i) = deviatoricStress(i) * 2.;
  }

  //compute second derivative of J2
  Matrix ddJ2ddstress(size, size);
  ddJ2ddstress.Zero();
  for (int i = 0; i < size; i++) {
    if (i < 3) {
      ddJ2ddstress(i, i) = 2. / 3.;
    }

    if (i > 2) {
      ddJ2ddstress(i, i) = 2.;
    }
  }

  ddJ2ddstress(0, 1) = -1. / 3.;
  ddJ2ddstress(0, 2) = -1. / 3.;
  ddJ2ddstress(1, 0) = -1. / 3.;
  ddJ2ddstress(1, 2) = -1. / 3.;
  ddJ2ddstress(2, 0) = -1. / 3.;
  ddJ2ddstress(2, 1) = -1. / 3.;

  //compute square of the first derivative of J2
  Matrix dJ2DJ2(size, size);
  for (int v = 0; v < size; v++) {
    for (int w = 0; w < size; ++w) {
      dJ2DJ2(v, w) = dJ2dstress(v) * dJ2dstress(w);
    }
  }

  //compute the second derivative of rho
  Matrix ddRhoddStress;
  ddRhoddStress = ddJ2ddstress;
  ddRhoddStress *= (1. / rho);
  Matrix help1;
  help1 = dJ2DJ2;
  help1 *= (-1. / (rho * rho * rho));
  ddRhoddStress += help1;
  answer = ddRhoddStress;
}

void
ConcreteDPM1::computeDCosThetaDStress(Vector &answer, const Vector &stress) const
{
  int size = stress.Size();

  //compute volumetric deviatoric split
  Vector deviatoricStress(6);
  double volumetricStress;
  stress.computeDeviatoricVolumetricSplit(deviatoricStress, volumetricStress);

  //compute the coordinates
  double rho = deviatoricStress.computeSecondCoordinate();

  //compute principal stresses and directions
  Vector principalDeviatoricStress;
  Matrix principalDir;
  deviatoricStress.computePrincipalValDir(principalDeviatoricStress, principalDir);

  //compute the derivative of s1 with respect to the cartesian stress
  Vector ds1DStress(size);
  ds1DStress(0) = principalDir(0, 0) * principalDir(0, 0) - 1. / 3.;
  ds1DStress(1) = principalDir(1, 0) * principalDir(1, 0) - 1. / 3.;
  ds1DStress(2) = principalDir(2, 0) * principalDir(2, 0) - 1. / 3.;
  ds1DStress(3) = 2. * principalDir(1, 0) * principalDir(2, 0);
  ds1DStress(4) = 2. * principalDir(2, 0) * principalDir(0, 0);
  ds1DStress(5) = 2. * principalDir(0, 0) * principalDir(1, 0);

  //compute dCosThetaDStress
  Vector dCosThetaDStress;
  dCosThetaDStress = ds1DStress;
  dCosThetaDStress *= (sqrt(3. / 2.) * rho / pow(rho, 2.));
  Vector help(size);
  computeDRhoDStress(help, stress);
  help *= (-sqrt(3. / 2.) * principalDeviatoricStress(0) / pow(rho, 2.));
  dCosThetaDStress += help;
  answer = dCosThetaDStress;
}

double
ConcreteDPM1::computeDRDCosTheta(const double theta, const double ecc) const
{
  double ACostheta = 4. * (1. - ecc * ecc) * cos(theta) * cos(theta) +
    (2. * ecc - 1.) * (2. * ecc - 1.);
  double BCostheta = 2. * (1. - ecc * ecc) * cos(theta) +
    (2. * ecc - 1.) * sqrt(4. * (1. - ecc * ecc) * cos(theta) * cos(theta)
    + 5. * ecc * ecc - 4. * ecc);
  double A1Costheta = 8. * (1. - pow(ecc, 2.)) * cos(theta);
  double B1Costheta = 2. * (1. - pow(ecc, 2.)) +
    4. * (2. * ecc - 1.) * (1. - pow(ecc, 2.)) * cos(theta) /
    sqrt(4. * (1. - pow(ecc, 2.)) * pow(cos(theta), 2.) +
    5. * pow(ecc, 2.) - 4. * ecc);
  double dRDCostheta = A1Costheta / BCostheta - ACostheta / pow(BCostheta, 2.) * B1Costheta;
  return dRDCostheta;
}

double
ConcreteDPM1::computeMeanSize()
{
  return 0.0;
}

