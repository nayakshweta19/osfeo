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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete10.cpp,v $
                                                                      
// Written: fmk
// Created: 03/06
//
// Description: This file contains the class implementation of Concrete10. 
// This Concrete10 is based on an f2c of the FEDEAS material
// Concr2.f which is:
//-----------------------------------------------------------------------
// concrete model with damage modulus    
//       by MOHD YASSIN (1993)
// adapted to FEDEAS material library
// by D. Sze and Filip C. Filippou in 1994
//-----------------------------------------------------------------------


#include <stdlib.h>
#include <Concrete10.h>
#include <OPS_Globals.h>
#include <math.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>

#define MAT_TAG_Concrete10 112244

#include <MaterialResponse.h>
#include <elementAPI.h>

#define OPS_Export 

OPS_Export void *
OPS_NewConcrete10Material(void)
{
	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial *theMaterial = 0;
	int numRemainingArgs = OPS_GetNumRemainingInputArgs();

    if (numRemainingArgs < 8) {
	  opserr << "WARNING invalid number of arguments\n";	
	  opserr << "Want: uniaxialMaterial Concrete10 tag? fpc? epsc0? fpcu? epscu? rat? ft? Ets?" << endln;
	  return 0;
	}

    double fpc, epsc0, fpcu, epscu;
    double rat, ft, Ets;
    int tag;
	int numData = 1;

    if (OPS_GetIntInput(&numData, &tag) != 0) {
	  opserr << "WARNING invalid uniaxialMaterial Concrete10 tag" << endln;
	  return 0;
	}      

    if (OPS_GetDoubleInput(&numData, &fpc) != 0) {
	  opserr << "WARNING invalid fpc\n";
	  return 0;	
	}

    if (OPS_GetDoubleInput(&numData, &epsc0) != 0) {
	  opserr << "WARNING invalid epsc0\n";
	  return 0;	
	}

    if (OPS_GetDoubleInput(&numData, &fpcu) != 0) {
	  opserr << "WARNING invalid fpcu\n";
	  return 0;	
	}

    if (OPS_GetDoubleInput(&numData, &epscu) != 0) {
	  opserr << "WARNING invalid epscu\n";
	  return 0;	
	}

    if (OPS_GetDoubleInput(&numData, &rat) != 0) {
	  opserr << "WARNING invalid rat\n";
	  return 0;	
	}

    if (OPS_GetDoubleInput(&numData, &ft) != 0) {
	  opserr << "WARNING invalid ft\n";
	  return 0;	
	}

    if (OPS_GetDoubleInput(&numData, &Ets) != 0) {
	  opserr << "WARNING invalid Ets\n";
	  return 0;	
	}
	// Parsing was successful, allocate the material
	theMaterial = new Concrete10(tag, fpc, epsc0, fpcu, epscu, rat, ft, Ets);

	if (theMaterial == 0) {
      opserr << "WARNING could not create uniaxialMaterial of type ConcreteL01\n";
      return 0;
    }

	return theMaterial;
}

Concrete10::Concrete10(int tag, double _fc, double _epsc0, double _fcu,
		       double _epscu, double _rat, double _ft, double _Ets):
  UniaxialMaterial(tag, MAT_TAG_Concrete10),
  fc(_fc), epsc0(_epsc0), fcu(_fcu), epscu(_epscu), rat(_rat), ft(_ft), Ets(_Ets)
{
  // Make all concrete parameters negative
  if (fc > 0.0)
    fc = -fc;

  if (epsc0 > 0.0)
    epsc0 = -epsc0;

  ecminP = 0.0;
  deptP = 0.0;

  eP = 2.0*fc/epsc0;
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
  e = 2.0*fc/epsc0;

  slip = 0.82;
  M = 0.4;
  Yint = 2.2;
  epslonTP = 0.001;

  rot1p = 0.00115;

  double Ct = 70.0/2.2337/100;
  double ratio = 0.025;
  double dim = 10;

  Alpha = Ct*(ratio*100)/(dim*25.4);
  friction = ft * Alpha;

}

Concrete10::Concrete10(void):
  UniaxialMaterial(0, MAT_TAG_Concrete10)
{
 
}

Concrete10::~Concrete10(void)
{
  // Does nothing
}

UniaxialMaterial*
Concrete10::getCopy(void)
{
  Concrete10 *theCopy = new Concrete10(this->getTag(), fc, epsc0, fcu, epscu, rat, ft, Ets);
  // History variables
  theCopy->TloadingState = TloadingState;

  // State variables
  theCopy->zeta = zeta;
  theCopy->epslonTP = epslonTP;
  theCopy->D = D;
  return theCopy;
}

double
Concrete10::getInitialTangent(void)
{
  return 2.0*fc/epsc0;
}

int
Concrete10::setTrialStrain(double trialStrain, double strainRate)
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

	/* <Region #1>- 1 - ecmin - 2 - ept - 3 - epn - 4  */

    this->Compr_Envlp(eps, epslonTP, sig, e);
    ecmin = eps;

  } else { /* <Region #2,3,4>- 1 - ecmin - 2 - ept - 3 - epn - 4  */

  // else, if the current strain is between the minimum strain and ept 
  // (which corresponds to zero stress) the material is in the unloading- 
  // reloading branch and the stress remains between sigmin and sigmax 
    
    // calculate strain-stress coordinates of point R that determines 
    // the reloading slope according to Fig.2.11 in EERC Report (corresponding eqn 2.31, 2.32.)
    // the strain of point R is epsR and the stress is sigmR 

    double epsr = (fcu - rat * ec0 * epscu) / (ec0 * (1.0 - rat));
    double sigmr = ec0 * epsr;
    
    // calculate the previous minimum stress sigmm from the minimum 
    // previous strain ecmin and the monotonic envelope in compression

    double sigmm;
    double dumy;
    this->Compr_Envlp(ecmin, epslonTP, sigmm, dumy);
    
    // calculate current reloading slope Er (Eq. 2.35 in EERC Report) 
    // calculate the intersection of the current reloading slope Er 
    // with the zero stress axis (variable ept) (Eq. 2.36 in EERC Report) 
    double er = (sigmm - sigmr) / (ecmin - epsr);
    double ept = ecmin - sigmm / er;
    
    if (eps <= ept) { /* <Region #2>- 1 - ecmin - 2 - ept - 3 - epn - 4  */

	  double sigmin = sigmm + er * (eps - ecmin);
	  double sigmax = er * .5f * (eps - ept);
	  
	  // compression unloading
	  if( deps > 0 ){
	  	
	  	double sigTemp1 = sigP + ec0 * deps;
	  	double eTemp1 = ec0;
	  	
	  	if (sigTemp1 <= sigmin) {
	  sigTemp1 = sigmin;
	  eTemp1 = er;
	  	}

	  	if (sigTemp1 >= sigmax) {
	  sigTemp1 = sigmax;
	  eTemp1 = 0.5 * er;
	  	}
	  	
	  	// Partial loading

	  	double eTemp2 = (sigP - 0.0)/(epsP - ept);
	  	double sigTemp2 = sigP + eTemp2 * deps;
	  	
	  	sig = __min(sigTemp1, sigTemp2);

	  	if (sig==sigTemp1) {
	  e = eTemp1;
	  	}

	  	if (sig==sigTemp2) {
	  e = eTemp2;
	  	}
	  
	  // compression reloading
	  } else { // if (deps <= 0)
	  	double sigTemp1;
	  	double eTemp1;
	  	this->Compr_Reloading(eps, epslonTP, epsP, sigP, ecmin, dept, friction,
	  					      Alpha, sigTemp1, eTemp1);
	  	
	  	double eTemp2 = (sigP - sigmm)/(epsP - ecmin);
	  	double sigTemp2 = sigP + eTemp2 * deps;
	  	
	  	sig = __max(sigTemp1, sigTemp2);

	  	if (sig==sigTemp1) {
	  e = eTemp1;
	  	} 

	  	if (sig==sigTemp2) {
	  e = eTemp2;
	  	}
	  }
	} else { 
	/* <Region #3,4>- 1 - ecmin - 2 - ept - 3 - epn - 4  */

	  // else, if the current strain is between ept and epn
	  // (which corresponds to maximum remaining tensile strength)
	  // the response corresponds to the reloading branch in tension
	  // Since it is not saved, calculate the maximum remaining tensile
	  // strength sicn (Eq. 2.43 in EERC Report)
	  
	  // calculate first the strain at the peak of the tensile stress-strain
	  // relation epn (Eq. 2.42 in EERC Report)

	  double epn = ept + dept;
	  double sicn;
	  if(eps<=epn){ 
	/* <Region #3>- 1 - ecmin - 2 - ept - 3 - epn - 4  */
	    
	    this->Tens_Envlp(dept, epsP, Alpha, sicn, e);
	    
	    if (dept != 0.0) {
	    // Tension reloading // Negative side
	      
	      if(eps < DBL_EPSILON){
	        if (deps > 0) { 
	      e = sicn / dept;
	      sig = e * (eps - ept);
	      
		//Partial loading
	      double sigTemp;
	      double eTemp;
	      this->Compr_Envlp(epsP, epslonTP, sigTemp, eTemp);
	      
	      if(sigP == sigTemp){
	        this->Compr_Envlp(eps, epslonTP, sig, e);
	      }

	      //if(sigP<-0.2){
	      //  this->Compr_Envlp(eps, epslonTP, sig, e);
	      //}

	        } else {
	      this->Tens_UnloadingN(eps, epslonTP, epsP, sigP, ecmin,
	          					dept, friction, Alpha, sig, e);
	        } //Partial loading
	        
	        //double eps0 = ft/ec0;
	        //if ( epn < eps0 ) {
	        //  e = 1.0 * ec0;
	        //  sig = sigP + 1.0*ec0 * deps;
	        //}
	      }

		  // Tension reloading // Positive side
		  else {
	        if ( deps > 0 ) {
	          double eTemp1 = sicn / dept;
	          double sigTemp1 = eTemp1 * (eps - ept);
	          double sigTemp2 = sigP + ec0*deps;
	          double eTemp2 = ec0;
	          
	          //Partial loading
	          sig = __min(sigTemp1, sigTemp2);

	          if (sig==sigTemp1) {
	        e = eTemp1;
	          } 

			  if (sig==sigTemp2) {
	        e = eTemp2;
	          }
	          
	      //Tension unloading - Positive side
	        } else {
	        //Stability - Tension Unloading slope
	          e = 1.0 * ec0;
	          sig = sigP + 1.0*ec0 * deps;
	          this->Tens_UnloadingP (eps, epslonTP, epsP, sigP, ecmin, dept, friction, Alpha, sig, e);
	        }
	        
	        double eps0 = ft/ec0; //Stability

	        if ( epn < eps0 ){
	          e = 1.0 * ec0;
	          sig = sigP + 1.0*ec0 * deps;
	        }
	      }

		//if (ept==0), linear
	    } else { 
	      e = ec0;
	      sig = e * (eps - ept);
	    }
	    
	    if ( ecmin == 0 ) {
	    // 1st iteration
	      e = sicn / dept;
	      sig = e * (eps - ept);
	    }
	  	
	  }else {
	  /* <Region #4>- 1 - ecmin - 2 - ept - 3 - epn - 4  */
	    
	    // else, if the current strain is larger than epn the response
	    // corresponds to the tensile envelope curve shifted by ept

	    double epstmp = eps - ept;
	    this->Tens_Envlp(epstmp, epsP, Alpha, sig, e);
	    dept = eps - ept;

	  }
	}
  }

  return 0;
}

double 
Concrete10::getStrain(void)
{
  return eps;
}

double 
Concrete10::getStress(void)
{
  return sig;
}

double 
Concrete10::getTangent(void)
{
  return e;
}

int 
Concrete10::commitState(void)
{
  ecminP = ecmin;
  deptP = dept;
  
  eP = e;
  sigP = sig;
  epsP = eps;
  return 0;
}

int 
Concrete10::revertToLastCommit(void)
{
  ecmin = ecminP;;
  dept = deptP;
  
  e = eP;
  sig = sigP;
  eps = epsP;
  return 0;
}

int 
Concrete10::revertToStart(void)
{
  ecminP = 0.0;
  deptP = 0.0;

  eP = 2.0*fc/epsc0;
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
  e = 2.0*fc/epsc0;
  
  // History variables

  TloadingState = 0;
  CloadingState = 0;
  reloadPath = 0;

  zeta     = 1.0;
  epslonTP = 0.0;
  beta     = 1.0;
  Wp       = 1.0;
  D        = 1.0;

  return 0;
}

int 
Concrete10::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(21);
  data(0) =fc;    
  data(1) =epsc0; 
  data(2) =fcu;   
  data(3) =epscu; 
  data(4) =rat;   
  data(5) =ft;    
  data(6) =Ets;   
  data(7) =ecminP;
  data(8) =deptP; 
  data(9) =epsP;  
  data(10) =sigP; 
  data(11) =eP;   
  data(12) = this->getTag();

  data(13) = zeta;     // Softening effect
  data(14) = epslonTP; // Strain in the perpendicular direction, needed to get the zeta
  data(15) = D;        // Damage factor for strength, get from parameter
  data(16) = X;        // for normal stresses 
  data(17) = K;        // for normal stresses
  data(18) = beta;     // Parameter needed for calculating zeta
  data(19) = fbeta;    // function of beta
  data(20) = Wp;       // pre-stressing factor

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "Concrete10::sendSelf() - failed to sendSelf\n";
    return -1;
  }
  return 0;
}

int 
Concrete10::recvSelf(int commitTag, Channel &theChannel, 
	     FEM_ObjectBroker &theBroker)
{
  static Vector data(21);

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "Concrete10::recvSelf() - failed to recvSelf\n";
    return -1;
  }

  fc = data(0);
  epsc0 = data(1);
  fcu = data(2);
  epscu = data(3);
  rat = data(4);
  ft = data(5);
  Ets = data(6);
  ecminP = data(7);
  deptP = data(8);
  epsP = data(9);
  sigP = data(10);
  eP = data(11);
  this->setTag(data(12));

  zeta    = data(13); // Softening effect
  epslonTP= data(14); // Strain in the perpendicular direction, needed to get the zeta
  D       = data(15); // Damage factor for strength, get from parameter
  X       = data(16); // for normal stresses 
  K       = data(17); // for normal stresses
  beta    = data(18); // Parameter needed for calculating zeta
  fbeta   = data(19); // function of beta
  Wp      = data(20); // pre-stressing factor

  e = eP;
  sig = sigP;
  eps = epsP;
  
  return 0;
}

void 
Concrete10::Print(OPS_Stream &s, int flag)
{
  s << "Concrete10:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
}

void
Concrete10::Tens_Envlp (double epsc, double epscP, double Alpha, double &sigc, double &Ect)
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

  if (epsc<=eps0) {
    sigc = epsc*Ec0;
    Ect  = Ec0;
  } else {
    double Lamda = __min(270./sqrt(Alpha),1000.);
	double Pwr = -1.*Lamda * (epsc - eps0);
	sigc = ft*(1-Alpha)*exp(Pwr)+Alpha;
	double PwrPrime = -1.*Lamda;
	Ect = ft*(1-Alpha)*PwrPrime*exp(Pwr);
  }
  return;
}
  
void
Concrete10::Compr_Envlp (double epsc, double epsc2, double &sigc, double &Ect)
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
  double zeta = 1./(0.8-0.34*epsc2/epsc0);
  
  /*softening effect- off*/
  if(epsc2<0){
	zeta=1.;
  }
  if(zeta>1.){
	zeta=1.;
  }
  /*end*/
  
  double ratLocal = epsc/epsc0;
  if (epsc>=epsc0) {
  	double Re = epsc/epsc0/zeta;
    sigc = zeta*fc*(2*Re-Re*Re);
    Ect  = zeta*fc*(2*(1./zeta/epsc0)-2.*(epsc/zeta/zeta/epsc0/epsc0));
  } else {
    
    //   linear descending branch between epsc0 and epscu
    if (epsc>epscu) {
      double fcp = zeta*fc;
      double epsc0p=zeta*epsc0;
      sigc = (fcu-fcp)*(epsc-epsc0p)/(epscu-epsc0p)+fc;
      Ect  = (fcu-fcp)/(epscu-epsc0p);
    } else {
	   
      // flat friction branch for strains larger than epscu
      sigc = fcu;
      Ect  = 1.0e-10; //       Ect  = 0.0
    }
  }
  return;
}

int
Concrete10::getVariable(const char *varName, Information &theInfo)
{
  if (strcmp(varName,"ec") == 0) {
    theInfo.theDouble = epsc0;
    return 0;
  } else
    return -1;
}

Response* 
Concrete10::setResponse(const char **argv, int argc,
OPS_Stream &theOutput)
{
  Response *theResponse = 0;
  if (strcmp(argv[0],"getPD") == 0) {
    double data = 0.0;
    theResponse = new MaterialResponse(this, 100, data);
  } else if (strcmp(argv[0],"setWallVar") == 0) {
    theResponse = new MaterialResponse(this, 101, Vector(5));
  } else
    return this->UniaxialMaterial::setResponse(argv, argc, theOutput);
  return theResponse;
}

int 
Concrete10::getResponse(int responseID, Information &matInfo)
{
  if (responseID == 100) {
    matInfo.theDouble = this->getPD();
  } else if (responseID == 101){
    Vector *theVector = matInfo.theVector;
    X = (*theVector)(0);
    K = (*theVector)(1);
    D = (*theVector)(2);
    beta = (*theVector)(3);
    epslonTP = (*theVector)(4);
  } else
    return this->UniaxialMaterial::getResponse(responseID, matInfo);

  return 0;
}

double
Concrete10::getPD ()
{
  double PD = 0.0;

  double tempRatio;	
  if ( epslonTP <= 0.0 )  //compression
  {
    PD = 0.0;
  }
  else  //tension for epslonTP, and will degrade compression strength
  {
    D = 1.0-0.4*(abs(epslonTP/epsc0));
	//opserr << "D=" << D << endln;
    if ( TloadingState == 1 )   //ascending
    {
	  tempRatio = eps/(zeta*epsc0);
      PD = - D * fbeta * Wp * 1160.0 * sqrt(-fc)* pow((1+400.0*epslonTP), -1.5)
         * pow(tempRatio,2.0);
	  //-\frac{1160. D \sqrt{\text{fc}} \text{$\epsilon $2}^2}{\text{$\epsilon $0}^2 \zeta ^2 \sqrt{1+\frac{400 \text{$\epsilon $1}}{\text{$\eta $p}}} (400. \text{$\epsilon $1}+\text{$\eta $p})}
	  //PD = -(1160. * D *sqrt(-fc)*pow(eps,2))/(pow(epsc0,2) * pow(zeta,2) * sqrt(1+(400*epslonTP)/beta) * (400.*epslonTP+beta));
	  //opserr << "PD=" << PD << endln;
    }
    else if ( TloadingState == 2 )   //descending
    {
      if ( abs(e) < 1e-8 ) // at the end platum part of descending branch FMK CHANGED FROM = 0.0
      {
        PD = 0.0;
      }
      else
      {
        tempRatio = eps/(zeta*epsc0);
        PD = - D * fbeta * Wp *1160.0 * sqrt(-fc)* pow((1+400.0*epslonTP), -1.5) 
           * (1.0 -(tempRatio-1)/pow(4.0/zeta-1.0,3.0)*(1.0-12.0/zeta+(4.0/zeta+1.0)*tempRatio)); 
        // 1160, 400
	    //-\frac{1160. D \sqrt{\text{fc}} \left(-16. \text{$\epsilon $0} \text{$\epsilon $2} \zeta +\text{$\epsilon $2}^2 (4.@ +\zeta )+\text{$\epsilon $0}^2 (-64.+48. \zeta )\right)}{\text{$\epsilon $0}^2 (-4.+\zeta )^3 \sqrt{1+\frac{400 \text{$\epsilon $1}}{\text{$\eta $p}}} (400. \text{$\epsilon $1}+\text{$\eta $p})}
	    //PD = -((1160.* D*sqrt(-fc)*(-16.*epsc0*eps*zeta+pow(eps,2)*(4.+zeta)+pow(epsc0,2)*(-64.+48.*zeta))))/(pow(epsc0,2)*pow((-4.+zeta),3)*sqrt(1+(400*epslonTP)/beta)*(400.*epslonTP+beta));
	    //opserr << "PD=" << PD << endln;
      }
    }		
    else
    {
      PD = 0.0;
    }
    if ( abs(zeta - 0.9) < 1e-8  || abs(zeta - 0.25) < 1e-8 ) // zeta = max or min value
    {
      PD = 0.0;
    }
  }
  
  return PD;
}

void
Concrete10::Tens_UnloadingP (double epsc, double epsc2, double epscP, double sigcP, double ecmin, double dept, double friction,
							double Alpha, double &sigc, double &Ect)
{
	//Stability - Tension Unloading slope
	double depsc = epsc - epscP;
	double ec0 = fc * 2. / epsc0;
	
	double epsr = (fcu - rat * ec0 * epscu) / (ec0 * (1.0 - rat));
	double sigmr = ec0 * epsr;
	double sigmm;
	double dumy;
	
	this->Compr_Envlp(ecmin, epsc2, sigmm, dumy);
	double er = (sigmm - sigmr) / (ecmin - epsr);
	double ept = ecmin - sigmm / er;
	double sicn;
	
	this->Tens_Envlp(dept, epscP, Alpha, sicn, Ect);
	double epn = ept + dept;
	
	if(epn<rot1p){
	  double eZEROsig = epn - sicn/ec0;
	  double ecmin4 = ecmin/8.;
	  double sigmm4;
	  
	  double Edummy;
	  this->Compr_Envlp(ecmin4, epsc2, sigmm4, Edummy);
	  double ex = (sigmm4 - sigcP) / (ecmin4 - epscP);
	  double sigtemp = sigmm4 + ex * (epsc - ecmin4);
	  double etemp = ex;
	  
	  if (sigc <=0) {
	    if (epsc>ecmin4){
	  sigc = sigtemp;
	  Ect = etemp;
	    } else {
	  this->Compr_Envlp(epsc, epsc2, sigc, Ect);
	    }
	  }
	} else {
	  double epstmp = epsc - ept;
	  double sigtmin = friction-0.000001*ec0*depsc;
	  double etmin = -0.000001*ec0;
	  
	  double epn = ept + dept;
	  double offset = - epsc0 - epn + slip * epn;
	  double ratLocal = -( 0.00001 + offset ) / epsc0;
	  double sigtcr = - ( fc * ratLocal * ( 2.0 - ratLocal )) - friction + fc;
	  
	  if ( -sigtmin > sigc) {
	    sigc = - sigtmin;
	    Ect = - etmin;
	    
	    double Ec0 = 2.0*fc/epsc0;
	    double eSlip0 = (1-slip)*epn; // slip strain
	    double offset = - epsc0 - epn + slip * epn;
	    
	    if ( eSlip0 >= epsc ) {
	      double eFF = -M*eSlip0;
	      double sigFF;
	      double Edummy;
	      this->Compr_Envlp(eFF, epsc2, sigFF, Edummy);
	      double Sige0 = sigFF/Yint;
	      double AA, BB, CC;
	      if ( abs(Sige0) < DBL_EPSILON ) {
	    AA = 0.0;
	    BB = 0.0;
	    CC = 0.0;
	      } else {
	    CC = Sige0;
	    BB = ((-friction) + (-sigFF/M)/M - Sige0 + Sige0/M/M)/(eSlip0+eSlip0/M);
	    AA = ((-friction) - BB*eSlip0 - Sige0)/eSlip0/eSlip0;
	      }
	      
	      double sigtmin2 = AA*epsc*epsc + BB*epsc + CC;
	      double etmin2 = 2*AA*epsc + BB;
	      
	      if (sigtmin2 < sigc) {
	    sigc = sigtmin2;
	    Ect = etmin2;
	      }
	    
	    }
	  }
	}
}

void
Concrete10::Tens_UnloadingN (double epsc, double epsc2, double epscP, double sigcP, double ecmin, double dept, double friction,
							double Alpha, double &sigc, double &Ect)
{
	//Stability - Tension Unloading slope
	double depsc = epsc - epscP;
	double ec0 = fc * 2. / epsc0;
	double epsr = (fcu - rat * ec0 * epscu) / (ec0 * (1.0 - rat));
	double sigmr = ec0 * epsr;
	double sigmm;
	double dumy;
	
	this->Compr_Envlp(ecmin, epsc2, sigmm, dumy);
	double er = (sigmm - sigmr) / (ecmin - epsr);
	double ept = ecmin - sigmm / er;
	double sicn;
	this->Tens_Envlp(dept, epscP, Alpha, sicn, Ect);
	double epn = ept + dept;
	
	double eZEROsig = epn - sicn/ec0;
	double ecmin4 = ecmin/8.;
	double sigmm4;
	double Edummy;
	this->Compr_Envlp(ecmin4, epsc2, sigmm4, Edummy);
	double ex = (sigmm4 - sigcP) / (ecmin4 - epscP);
	double sigtemp = sigmm4 + ex * (epsc - ecmin4);
	double etemp = ex;
	
	if (epn<rot1p) {
	  if (sigc <= 0) {
	if (epsc>ecmin4){
	  sigc = sigtemp;
	  Ect = etemp;
	} else {
	  this->Compr_Envlp(epsc, epsc2, sigc, Ect);
	}
	  } else {
	Ect = ec0;
	sigc = sigcP+Ect*depsc;
	  }
	} else {
	  double eSlip0 = (1-slip)*epn; // slip strain
	  double eFF = -M*eSlip0;
	  double sigFF;
	  double Edummy;
	  this->Compr_Envlp(eFF, epsc2, sigFF, Edummy);
	  double Sige0 = sigFF/Yint;
	  double AA, BB, CC;
	  //if ( abs(Sige0) < DBL_EPSILON) {
	  //	AA = 0.0;
	  //	BB = 0.0;
	  //	CC = 0.0;
	  //} else {
	    CC = Sige0;
	    BB = ((-friction) + (-sigFF/M)/M - Sige0 + Sige0/M/M)/(eSlip0+eSlip0/M);
	    AA = ((-friction) - BB*eSlip0 - Sige0)/eSlip0/eSlip0;
	  //}
	  
	  if (epsc >= eFF ){
	sigc = AA*epsc*epsc + BB*epsc +CC;
	Ect= 2*AA*epsc + BB;
	  }else{
	this->Compr_Envlp(epsc, epsc2, sigc, Ect);
	  }
	}
}

void
Concrete10::Compr_Reloading (double epsc, double epsc2, double epscP, double sigcP, double ecmin, double dept, double friction,
							double Alpha, double &sigc, double &Ect)
{
	//Stability - Tension Unloading slope
	double depsc = epsc - epscP;
	double ec0 = fc * 2. / epsc0;
	double epsr = (fcu - rat * ec0 * epscu) / (ec0 * (1.0 - rat));
	double sigmr = ec0 * epsr;
	double sigmm;
	double dumy;
	
	this->Compr_Envlp(ecmin, epsc2, sigmm, dumy);
	double er = (sigmm - sigmr) / (ecmin - epsr);
	double ept = ecmin - sigmm / er;
	
	double sicn;
	this->Tens_Envlp(dept, epscP, Alpha, sicn, Ect);
	double epn = ept + dept;
	double eSlip0 = (1-slip)*epn; // slip strain
	
	double eFF = -M * eSlip0;
	double sigFF;
	double Edummy;
	this->Compr_Envlp(eFF, epsc2, sigFF, Edummy);
	double Sige0 = sigFF/Yint;
	double AA, BB, CC;
	//if ( abs(Sige0) < DBL_EPSILON) {
	//	AA = 0.0;
	//	BB = 0.0;
	//	CC = 0.0;
	//} else {
	  CC = Sige0;
	  BB = ((-friction) + (-sigFF/M)/M - Sige0 + Sige0/M/M)/(eSlip0+eSlip0/M);
	  AA = ((-friction) - BB*eSlip0 - Sige0)/eSlip0/eSlip0;
	//}
	
	if( epsc >= eFF ){
	  sigc = AA*epsc*epsc + BB*epsc +CC;
	  Ect= 2*AA*epsc + BB;
	}else{
	  this->Compr_Envlp(epsc, epsc2, sigc, Ect);
	}
	
	if(epn<rot1p){
	  //double epn = ept + dept;
	  double sicn;
	  this->Tens_Envlp(dept, epscP, Alpha, sicn, Ect);
	  double eZEROsig = epn - sicn/ec0;
	  double DeltaEps = (ecmin - epscP);
	  if(DeltaEps==0){
	    DeltaEps =1.0e-10;
	  }
	  Ect = (sigmm - sigcP) / DeltaEps;
	  sigc = sigmm + Ect * (epsc - ecmin);
	}
}
