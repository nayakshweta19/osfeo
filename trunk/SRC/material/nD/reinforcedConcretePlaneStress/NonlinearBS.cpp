//# REF: (C) Copyright 1999, The Regents of the University of California
#include <stdlib.h>
#include <NonlinearBS.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <math.h>

#include <Information.h>
#include <ID.h>
#include <MaterialResponse.h>

#include <iostream>
#include <fstream>

#include <DummyStream.h>
#include <elementAPI.h>

using namespace std;

Vector NonlinearBS :: sig_p(3);
Vector NonlinearBS :: eps_p(3);
Vector NonlinearBS :: strain_vec(3);
Vector NonlinearBS :: eSlip(3);

Vector NonlinearBS :: sig_XY(3);
Vector NonlinearBS :: stress_vec(3);
Vector NonlinearBS :: sigSlip(3);
Matrix NonlinearBS :: tangent_matrix(3,3);

//material stiffness
Matrix NonlinearBS :: Ds(3,3);
Matrix NonlinearBS :: Dc(3,3);

//transformation
Matrix NonlinearBS :: T(3,3);
Matrix NonlinearBS :: Ts1(3,3);
Matrix NonlinearBS :: Ts2(3,3);

//tangent
Matrix NonlinearBS :: Dxy(3,3);


#define OPS_Export 

static int numNonlinearBSMaterials = 0;

OPS_Export void *
OPS_NewNonlinearBSMaterial()
{
  if (numNonlinearBSMaterials == 0) {
    numNonlinearBSMaterials++;
  }

  // Pointer to a uniaxial material that will be returned
  NDMaterial *theMaterial = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs != 27 && numRemainingArgs != 31) {
    opserr << "Invalid Args want: NDMaterial NonlinearBS matTag? fc? epsc0? fcu? epscu? rat? ft? Ets? reinfD1? reinfR1? reinf_dmm1? reinfD2? reinfR2? reinf_dmm2? m1p? r1p? m2p? r2p? <m3p? r3p?> m1n? r1n? m2n? r2n? <m3n? r3n?> px? py? d1? d2? b? \n";
    return 0;	
  }

  int tag;
  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid nDMaterial NonlinearBS tag" << endln;
    return 0;
  }

  if ( numRemainingArgs == 27) {
	double dData[26];
	int numData = 26;
	if (OPS_GetDouble(&numData, dData) != 0) {
	  opserr << "WARNING invalid data NonlinearBS tag: " << tag << endln;
	  return 0;
	}

	//now create the NonlinearBS
	theMaterial = new NonlinearBS (tag,
		dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
		dData[6], dData[7], dData[8], dData[9],dData[10],dData[11],
		dData[12],dData[13],dData[14],dData[15],dData[16],dData[17],
		dData[18],dData[19],dData[20],dData[21],dData[22],dData[23],
		dData[24],dData[25]);
  }
  if (numRemainingArgs == 31) {
	
	double dData[30];
	int numData = 30;
	if (OPS_GetDouble(&numData, dData) != 0) {
	  opserr << "WARNING invalid data NonlinearBS tag: " << tag << endln;
	  return 0;
	}

	//now create the NonlinearBS
	theMaterial = new NonlinearBS (tag,
		dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
		dData[6], dData[7], dData[8], dData[9],dData[10],dData[11],
		dData[12],dData[13],dData[14],dData[15],dData[16],dData[17],
		dData[18],dData[19],dData[20],dData[21],dData[22],dData[23],
		dData[24],dData[25],dData[26],dData[27],dData[28],dData[29]);
  }

  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory creating material\n";
    opserr << "NonlinearBS: " << tag << endln;
    return 0;
  }

  return theMaterial;
}

NonlinearBS::NonlinearBS(int tag, double fc, double epsc0, double fcu,
					double epscu, double rat, double ft, double Ets,
					double reinfD1, double reinfR1, double reinf_dmm1, double reinfD2,double reinfR2, double reinf_dmm2,
					double m1p, double r1p, double m2p, double r2p, double m3p, double r3p,
					double m1n, double r1n, double m2n, double r2n, double m3n, double r3n,
					double px, double py, double d1, double d2, double b):
NDMaterial(tag,ND_TAG_NonlinearBS),
  fc(fc), epsc0(epsc0), fcu(fcu), epscu(epscu), rat(rat), ft(ft), Ets(Ets),
  reinfD1(reinfD1), reinfR1(reinfR1), reinf_dmm1(reinf_dmm1), reinfD2(reinfD2), reinfR2(reinfR2), reinf_dmm2(reinf_dmm2),
  pinchX(px), pinchY(py), damfc1(d1), damfc2(d2), beta(b),
  mom1p(m1p), rot1p(r1p), mom2p(m2p), rot2p(r2p), mom3p(m3p), rot3p(r3p),
  mom1n(m1n), rot1n(r1n), mom2n(m2n), rot2n(r2n), mom3n(m3n), rot3n(r3n)
{
	ecminPX = 0.0;
	deptPX = 0.0;
	
	ecminPY = 0.0;
	deptPY = 0.0;
	
	ePX = 2.0*fc/epsc0;
	epsPX = 0.0;
	sigPX = 0.0;
	epsX = 0.0;
	sigX = 0.0;
	eX = ePX;
	
	ePY = 2.0*fc/epsc0;
	epsPY = 0.0;
	sigPY = 0.0;
	epsY = 0.0;
	sigY = 0.0;
	eY = ePY;
	
	double Ct = 70.0/2.2337/100; //10;
	Tens_Alpha1 = Ct*(reinfR1*100)/(reinf_dmm1);
	Tens_Alpha2 = Ct*(reinfR2*100)/(reinf_dmm2);
	frictionX = ft*Tens_Alpha1;
	frictionY = ft*Tens_Alpha2;
	slip = 0.82;
	
	//reduced steel stress envelop.
	ReducedL = 0.96 - frictionX/mom1p;
	ReducedT = 0.96 - frictionY/mom1p;
	
	//Shear retention
	ShearBeta = 0.08; //SW4
	
	//Crack closing curve:
	M = 0.4; //Multiplier
	Yint = 2.2;//Yint.
	thetaSP = 0;
	beta = 0;
	
	bool error = false;
	// Positive backbone parameters
	if(rot1p<=0.0)
		error = true;
	
	if(rot2p<=rot1p)
		error = true;
	
	if(rot3p<=rot2p)
		error = true;
	
	// Negative backbone parameters
	if(rot1n>=0.0)
		error = true;
	if(rot2n>=rot1n)
		error = true;
	if (rot3n >= rot2n)
		error = true;
	
	if (error) {
		opserr << "NonlinearBS::NonlinearBS material - input backbone is not unique (one-to-one)\n";
		exit(-1);
	}
	
	//opserr << mom1p << "\t" << rot1p << endln;
	//opserr << mom2p << "\t" << rot2p << endln;
	//opserr << mom3p << "\t" << rot3p << endln;
	//opserr << mom1n << "\t" << rot1n << endln;
	//opserr << mom2n << "\t" << rot2n << endln;
	//opserr << mom3n << "\t" << rot3n << endln;
	//opserr << pinchX << "\t" << pinchY << endln;
	
	energyA = 0.5 * (rot1p*mom1p + (rot2p-rot1p)*(mom2p+mom1p) + (rot3p-rot2p)*(mom3p+mom2p) +
			  rot1n*mom1n + (rot2n-rot1n)*(mom2n+mom1n) + (rot3n-rot2n)*(mom3n+mom2n));
	
	// Set envelope slopes
	this->setEnv();
	
	// Initialize history variables
	this->revertToStart();
	this->revertToLastCommit();
}

NonlinearBS::NonlinearBS(int tag, double fc, double epsc0, double fcu,
						double epscu, double rat, double ft, double Ets,
						double reinfDl, double reinfRl, double reinf_dmm1, double reinfD2, double reinfR2, double reinf_dmm2,
						double m1p, double r1p, double m2p, double r2p,
						double m1n, double r1n, double m2n, double r2n,
						double px, double py, double d1, double d2, double b):
NDMaterial(tag, ND_TAG_NonlinearBS),
fc(fc), epsc0(epsc0), fcu(fcu), epscu(epscu), rat(rat), ft(ft), Ets(Ets),
reinfD1(reinfD1), reinfR1(reinfR1), reinfD2(reinfD2), reinfR2(reinfR2),
pinchX(px), pinchY(py), damfc1(d1), damfc2(d2), beta(b),
mom1p(m1p), rot1p(r1p), mom3p(m2p), rot3p(r2p),
mom1n(m1n), rot1n(r1n), mom3n(m2n), rot3n(r2n)
{
	bool error = false;
	// Positive backbone parameters
	if(rot1p<=0.0)
		error = true;
	
	if(rot3p<=rot1p)
		error = true;
	
	// Negative backbone parameters
	if(rot1n>=0.0)
		error = true;
	
	if (rot3n >= rot1n)
		error = true;
	
	if (error) {
		opserr << "NonlinearBS::NonlinearBS Material - input backbone is not unique (one-to-one)\n";
		exit(-1);
	}
	
	energyA = 0.5 * (rot1p*mom1p + (rot3p-rot1p)*(mom3p+mom1p) +
			  rot1n*mom1n + (rot3n-rot1n)*(mom3n+mom1n) );
	
	mom2p = 0.5*(mom1p+mom3p);
	mom2n = 0.5*(mom1n+mom3n);
	
	rot2p = 0.5*(rot1p+rot3p);
	rot2n = 0.5*(rot1n+rot3n);
	
	//opserr << mom1p << "\t" << rot1p << endln;
	//opserr << mom2p << "\t" << rot2p << endln;
	//opserr << mom3p << "\t" << rot3p << endln;
	//opserr << mom1n << "\t" << rot1n << endln;
	//opserr << mom2n << "\t" << rot2n << endln;
	//opserr << mom3n << "\t" << rot3n << endln;

	// Set envelope slopes
	this->setEnv();
	
	// Initialize history variables
	this->revertToStart();
	this->revertToLastCommit();
}

NonlinearBS::NonlinearBS(void):
NDMaterial(0, ND_TAG_NonlinearBS),
fc(0.0), epsc0(0.0), fcu(fcu), epscu(0.0), rat(0.0), ft(0.0), Ets(0.0),
reinfD1(0.0), reinfR1(0.0), reinf_dmm1(0.0), reinfD2(0.0), reinfR2(0.0), reinf_dmm2(0.0),
pinchX(0.0), pinchY(0.0), damfc1(0.0), damfc2(0.0), beta(0.0),
mom1p(0.0), rot1p(0.0), mom2p(0.0), rot2p(0.0), mom3p(0.0), rot3p(0.0),
mom1n(0.0), rot1n(0.0), mom2n(0.0), rot2n(0.0), mom3n(0.0), rot3n(0.0)
{
	
}

NonlinearBS::~NonlinearBS(void)
{
// Does nothing
}

NDMaterial*
NonlinearBS::getCopy(void)
{
	NonlinearBS *theCopy = new NonlinearBS(this->getTag(), fc, epsc0, fcu, epscu, rat, ft, Ets,
	reinfD1, reinfR1, reinf_dmm1, reinfD2, reinfR2, reinf_dmm2,
	mom1p, rot1p, mom2p, rot2p, mom3p, rot3p,
	mom1n, rot1n, mom2n, rot2n, mom3n, rot3n,
	pinchX, pinchY, damfc1, damfc2, beta);
	
	return theCopy;
}

NDMaterial*
NonlinearBS :: getCopy (const char *type)
{
	if (strcmp(type,"PlaneStress2D") == 0
		|| strcmp(type,"PlaneStress") == 0){
		NonlinearBS *clone;
		clone = new NonlinearBS(this->getTag(), fc, epsc0, fcu, epscu, rat, ft, Ets,
								reinfD1, reinfR1, reinf_dmm1, reinfD2, reinfR2, reinf_dmm2,
								mom1p, rot1p, mom2p, rot2p, mom3p, rot3p,
								mom1n, rot1n, mom2n, rot2n, mom3n, rot3n,
								pinchX, pinchY, damfc1, damfc2, beta);
		
		clone ->CrotMaxL = CrotMaxL;
		clone ->CrotMinL = CrotMinL;
		clone ->CrotPuL = CrotPuL;
		clone ->CrotNuL = CrotNuL;
		clone ->CenergyDL = CenergyDL;
		clone ->CloadIndicatorL = CloadIndicatorL;
		clone ->CstressL = CstressL;
		clone ->CstrainL = CstrainL;
		clone ->TtangentL = TtangentL;
		
		clone ->CrotMaxT = CrotMaxT;
		clone ->CrotMinT = CrotMinT;
		clone ->CrotPuT = CrotPuT;
		clone ->CrotNuT = CrotNuT;
		clone ->CenergyDT = CenergyDT;
		clone ->CloadIndicatorT = CloadIndicatorT;
		clone ->CstressT = CstressT;
		clone ->CstrainT = CstrainT;
		clone ->TtangentT = TtangentT;
		
		return clone;
	}
	
	// Handle other cases
	else {
		opserr << "NonlinearBS::getModel failed to get model: " << type << endln;
		return 0;
	}
}

const char*
NonlinearBS::getType() const
{
	return "PlaneStress";
}

int
NonlinearBS::setTrialStrain(const Vector &strain_from_element)
{
	double ec0 = fc * 2. / epsc0;
	
	// retrieve concrete history variables
	ecminX = ecminPX;
	deptX = deptPX;
	
	ecminY = ecminPY;
	deptY = deptPY;
	
	// calculate current strain
	eOXX = strain_from_element(0);
	eOYY = strain_from_element(1);
	gOXY = strain_from_element(2);
	
	//double thetaEP=0.;
	
	double A = (eOXX-eOYY)/2.;
	double B = gOXY/2.;
	double C = pow(A,2.)+pow(B,2.);
	
	//principal strain dir
	double ratio, thetaEP;
	if ( abs(eOXX-eOYY) < DBL_EPSILON && abs(gOXY) > DBL_EPSILON) {
	  ratio = 1e12; // for the pure shear state
	} else {
	  ratio = gOXY/(eOXX-eOYY);
	}
	thetaEP = atan(ratio)/2.;

	//<Deviation angle between principal stresses and strains>
	betaR = thetaSP - thetaEP;
	
	//Concrete transformation angle - 1st crack direction
	thetaC = thetaSP;
	
	//Transformation Matrix//
	T(0,0) = pow(cos(thetaC),2.);
	T(0,1) = pow(sin(thetaC),2.);
	T(0,2) = 1.0*sin(thetaC)*cos(thetaC);
	T(1,0) = pow(sin(thetaC),2.);
	T(1,1) = pow(cos(thetaC),2.);
	T(1,2) = -1.0*sin(thetaC)*cos(thetaC);
	T(2,0) = -2.0*sin(thetaC)*cos(thetaC);
	T(2,1) = 2.0*sin(thetaC)*cos(thetaC);
	T(2,2) = pow(cos(thetaC),2.)-pow(sin(thetaC),2.);
	
	//T1(0,0) = pow(cos(-thetaC),2.);
	//T1(0,1) = pow(sin(-thetaC),2.);
	//T1(0,2) = 1.0*sin(-thetaC)*cos(-thetaC);
	//T1(1,0) = pow(sin(-thetaC),2.);
	//T1(1,1) = pow(cos(-thetaC),2.);
	//T1(1,2) = -1.0*sin(-thetaC)*cos(-thetaC);
	//T1(2,0) = -2.0*sin(-thetaC)*cos(-thetaC);
	//T1(2,1) = 2.0*sin(-thetaC)*cos(-thetaC);
	//T1(2,2) = pow(cos(-thetaC),2.)-pow(sin(-thetaC),2.);
	
	//Strains in the first crack direction // eps_p =  = T*strain_from_element;

	eps_p.addMatrixVector(0.0, T, strain_from_element, 1.0);
	
	//Concrete 02 model (principal strain-stress relationship)
	epsX = eps_p(0);
	double depsX = epsX - epsPX;
	
	epsY = eps_p(1);
	double depsY = epsY - epsPY;
	
	gammaXY = eps_p(2);
	
	//<CONCRETE>-------------------------------------------------------------------
	
	/* --- 1 ---- ecmin ---- 2 ---- ept ---- 3 ---- epn ---- 4 */
	
	// ### <principal d-crack dir> ############################
	
	// if the current strain is less than the smallest previous strain
	// call the monotonic envelope in compression and reset minimum strain
	
	if( epsX < ecminX ) {
		
	/* <Region # 1 > --- 1 ---- ecmin ---- 2 ---- ept ---- 3 ---- epn ---- 4  */
		this->Compr_Envlp(epsX, epsY, sigX, eX);
		ecminX = epsX;
	} else {
	/* <Region #2,3, and 4>--- 1 ---- ecmin ---- 2 ---- ept ---- 3 ---- epn ---- 4  */
			
		// else, if the current strain is between the minimum strain and ept
		// (which corresponds to zero stress) the material is in the unloading-
		// reloading branch and the stress remains between sigmin and sigmax
		
		// calculate strain-stress coordinates of point R that determines
		// the reloading slope according to Fig.2.11 in EERC Report
		// (corresponding equations are 2.31 and 2.32
		// the strain of point R is epsR and the stress is sigmR
		
		double epsrX = (fcu - rat * ec0 * epscu) / (ec0 * (1.0 - rat));
		double sigmrX = ec0 * epsrX;
		
		// calculate the previous minimum stress sigmm from the minimum
		// previous strain ecmin and the monotonic envelope in compression
		
		double sigmmX;
		double dumyX;
		
		this->Compr_Envlp(ecminX, epsY, sigmmX, dumyX);
		
		// calculate current reloading slope Er (Eq. 2.35 in EERC Report)
		// calculate the intersection of the current reloading slope Er
		// with the zero stress axis (variable ept) (Eq. 2.36 in EERC Report)
		
		double erX = (sigmmX - sigmrX) / (ecminX - epsrX);
		double eptX = ecminX - sigmmX / erX;
		
		if (epsX <= eptX) {
			
	/* <Region #2>--- 1 ---- ecmin ---- 2 ---- ept ---- 3 ---- epn ---- 4  */
		
			double sigminX = sigmmX + erX * (epsX - ecminX);
			double sigmaxX = erX * .5f * (epsX - eptX);
			
			// compression unloading
			if( depsX > 0 ){
				
				double sigTemp1X = sigPX + ec0 * depsX;
				double eTemp1X = ec0;
				
				if (sigTemp1X <= sigminX) {
			sigTemp1X = sigminX;
			eTemp1X = erX;
				}
				
				if (sigTemp1X >= sigmaxX) {
			sigTemp1X = sigmaxX;
			eTemp1X = 0.5 * erX;
				}
				
				//Partial loading

				double eTemp2X = (sigPX - 0.0)/(epsPX - eptX);
				double sigTemp2X = sigPX + eTemp2X * depsX;
				
				sigX = min(sigTemp1X, sigTemp2X);
				
				if (sigX==sigTemp1X) {
			eX = eTemp1X;
				}
				
				if (sigX==sigTemp2X) {
			eX = eTemp2X;
				}
				
				
			// compression reloading
			} else {
				double sigTemp1X;
				double eTemp1X;
				this->Compr_Reloading(epsX, epsY, epsPX, sigPX, ecminX, deptX, frictionX,
										Tens_Alpha1, sigTemp1X, eTemp1X);
				
				double eTemp2X = (sigPX - sigmmX)/(epsPX - ecminX);
				double sigTemp2X = sigPX + eTemp2X * depsX;
				
				sigX = max(sigTemp1X, sigTemp2X);
				
				if (sigX==sigTemp1X) {
			eX = eTemp1X;
				}
				
				if (sigX==sigTemp2X) {
			eX = eTemp2X;
				}
			}
		} else {
	/* <Region #3 and #4 >--- 1 ---- ecmin ---- 2 ---- ept ---- 3 ---- epn ---- 4  */
			
			// else, if the current strain is between ept and epn
			// (which corresponds to maximum remaining tensile strength)
			// the response corresponds to the reloading branch in tension
			// Since it is not saved, calculate the maximum remaining tensile
			// strength sicn (Eq. 2.43 in EERC Report)
			
			// calculate first the strain at the peak of the tensile stress-strain
			// relation epn (Eq. 2.42 in EERC Report)
			
			double epnX = eptX + deptX;
			double sicnX;
			if(epsX<=epnX){
			/* <Region #3 >--- 1 ---- ecmin ---- 2 ---- ept ---- 3 ---- epn ---- 4  */
				
				this->Tens_Envlp(deptX, epsPX, Tens_Alpha1, sicnX, eX);
				
				if (deptX != 0.0) {
				// Tension reloading
				
					// Negative side
					if(epsX<1e-10){
						
						if (depsX>0) {
							eX = sicnX / deptX;
							sigX = eX * (epsX - eptX);
							
							//Partial loading
							double sigTempX;
							double eTempX;
							this->Compr_Envlp(epsPX, epsPY, sigTempX, eTempX);
							
							if(sigPX == sigTempX){
								this->Compr_Envlp(epsX, epsY, sigX, eX);
							}
							
							//if(sigPX<-0.2){
							//	this->Compr_Envlp(epsX, epsY, sigX, eX);
							//}
							
						} else {
							this->Tens_UnloadingN (epsX, epsY, epsPX, sigPX, ecminX,
												deptX, frictionX, Tens_Alpha1, sigX, eX);
							//Partial loading
						}
						
						//double eps0 = ft/ec0;
						//
						//if ( epnX < eps0 ) {
						//	eX=1.0*ec0;
						//	sigX = sigPX + 1.0*ec0 * depsX;
						//}
						
					} else {
					// Positive side
						if ( depsX > 0 ) {
							double eTemp1X = sicnX / deptX;
							double sigTemp1X = eTemp1X * (epsX - eptX);
							double sigTemp2X = sigPX + ec0*depsX;
							double eTemp2X = ec0;
							
							//Partial loading
							
							sigX = min(sigTemp1X, sigTemp2X);
							
							if (sigX==sigTemp1X) {
								eX = eTemp1X;
							}
							
							if (sigX==sigTemp2X) {
								eX = eTemp2X;
							}
						
					//Tension unloading - Positive side
						} else {
						//Stability - Tension Unloading slope
						
							eX=1.0*ec0;
							sigX = sigPX + 1.0*ec0 * depsX;
							this->Tens_UnloadingP (epsX, epsY, epsPX, sigPX, ecminX, deptX, frictionX, Tens_Alpha1, sigX, eX);
						}
						
						double eps0 = ft/ec0; //Stability
						
						if ( epnX < eps0 ){
							eX=1.0*ec0;
							sigX = sigPX + 1.0*ec0 * depsX;
						}
					}
				//if (eptY==0), linear
				} else {
					eX = ec0;
					sigX = eX * (epsX - eptX);
				}
				
				if ( ecminX == 0 ) {
					// 1st iteration
					eX = sicnX / deptX;
					sigX = eX * (epsX - eptX);
				}
				
			}else {
			/* <Region #4 >--- 1 ---- ecmin ---- 2 ---- ept ---- 3 ---- epn ---- 4  */
				
				// else, if the current strain is larger than epn the response
				// corresponds to the tensile envelope curve shifted by ept
				
				double epstmpX = epsX - eptX;
				this->Tens_Envlp(epstmpX, epsPX, Tens_Alpha1, sigX, eX);
				deptX = epsX - eptX;
			}
		}
	}
	
	//<< end of d dir >
	
	/* --- 1 ---- ecmin ---- 2 ---- ept ---- 3 ---- epn ---- 4 */
	
	//### <principal r-crack dir> #############################
	
	// if the current strain is less than the smallest previous strain
	// call the monotonic envelope in compression and reset minimum strain
	
	if( epsY < ecminY ) {
		
	/* <Region # 1 > --- 1 ---- ecmin ---- 2 ---- ept ---- 3 ---- epn ---- 4  */
		this->Compr_Envlp(epsY, epsX, sigY, eY);
		ecminY = epsY;
	} else {
	/* <Region #2,3, and 4>--- 1 ---- ecmin ---- 2 ---- ept ---- 3 ---- epn ---- 4  */
		double epsrY = (fcu - rat * ec0 * epscu) / (ec0 * (1.0 - rat));
		double sigmrY = ec0 * epsrY;
		double sigmmY;
		double dumyY;
		
		this->Compr_Envlp(ecminY, epsX, sigmmY, dumyY);
		
		double erY = (sigmmY - sigmrY) / (ecminY - epsrY);
		double eptY = ecminY - sigmmY / erY;
		
		if (epsY <= eptY) {
			
	/* <Region #2>--- 1 ---- ecmin ---- 2 ---- ept ---- 3 ---- epn ---- 4  */
		
			double sigminY = sigmmY + erY * (epsY - ecminY);
			double sigmaxY = erY * .5f * (epsY - eptY);
			
			// compression unloading
			if( depsY > 0 ){
				
				double sigTemp1Y = sigPY + ec0 * depsY;
				double eTemp1Y = ec0;
				
				if (sigTemp1Y <= sigminY) {
					sigTemp1Y = sigminY;
					eTemp1Y = erY;
				}
				
				if (sigTemp1Y >= sigmaxY) {
					sigTemp1Y = sigmaxY;
					eTemp1Y = 0.5 * erY;
				}
				
				//Partial loading
				double eTemp2Y = (sigPY - 0.0)/(epsPY - eptY);
				double sigTemp2Y = sigPY + eTemp2Y * depsY;
				
				sigY = min(sigTemp1Y, sigTemp2Y);
				
				if (sigY==sigTemp1Y) {
					eY = eTemp1Y;
				}
				
				if (sigY==sigTemp2Y) {
					eY = eTemp2Y;
				}

			// compression reloading
			} else {
				double sigTemp1Y;
				double eTemp1Y;
				this->Compr_Reloading(epsY, epsX, epsPY, sigPY, ecminY, deptY, frictionY,
									  Tens_Alpha2, sigTemp1Y, eTemp1Y);
				
				double eTemp2Y = (sigPY - sigmmY)/(epsPY - ecminY);
				double sigTemp2Y = sigPY + eTemp2Y * depsY;
				
				sigY = max(sigTemp1Y, sigTemp2Y);
				
				if (sigY==sigTemp1Y) {
					eY = eTemp1Y;
				}
				
				if (sigY==sigTemp2Y) {
					eY = eTemp2Y;
				}
			}
		} else {
	/* <Region #3 and #4 >--- 1 ---- ecmin ---- 2 ---- ept ---- 3 ---- epn ---- 4  */
			
			double epnY = eptY + deptY;
			double sicnY;
			if(epsY<=epnY){
			/* <Region #3 >--- 1 ---- ecmin ---- 2 ---- ept ---- 3 ---- epn ---- 4  */
				
				this->Tens_Envlp(deptY, epsPY, Tens_Alpha2, sicnY, eY);
				
				if (deptY != 0.0) {
				// Tension reloading
				
					// Negative side
					if(epsY<1e-10){
						
						if (depsY>0) {
							eY = sicnY / deptY;
							sigY = eY * (epsY - eptY);
							
							//Partial loading
							double sigTempY;
							double eTempY;
							this->Compr_Envlp(epsPY, epsPX, sigTempY, eTempY);
							
							if(sigPY == sigTempY){
								this->Compr_Envlp(epsY, epsX, sigY, eY);
							}
							
							//if(sigPY<-0.2){
							//	this->Compr_Envlp(epsY, epsX, sigY, eY);
							//}
							
						} else {
							this->Tens_UnloadingN (epsY, epsX, epsPY, sigPY, ecminY,
												deptY, frictionY, Tens_Alpha2, sigY, eY);
							//Partial loading
						}
						
						//double eps0 = ft/ec0;
						//
						//if ( epnY < eps0 ) {
						//	eY=1.0*ec0;
						//	sigY = sigPY + 1.0*ec0 * depsY;
						//}
						
					} else {
					// Positive side
						if ( depsY > 0 ) {
							double eTemp1Y = sicnY / deptY;
							double sigTemp1Y = eTemp1Y * (epsY - eptY);
							double sigTemp2Y = sigPY + ec0*depsY;
							double eTemp2Y = ec0;
							
							//Partial loading
							
							sigY = min(sigTemp1Y, sigTemp2Y);
							
							if (sigY==sigTemp1Y) {
								eY = eTemp1Y;
							}
							
							if (sigY==sigTemp2Y) {
								eY = eTemp2Y;
							}
						
					//Tension unloading - Positive side
						} else {
						//Stability - Tension Unloading slope
						
							eY=1.0 * ec0;
							sigY = sigPY + 1.0 * ec0 * depsY;
							this->Tens_UnloadingP (epsY, epsX, epsPY, sigPY, ecminY, deptY, frictionY,
								                   Tens_Alpha2, sigY, eY);
						}
						
						double eps0 = ft/ec0; //Stability
						
						if ( epnY < eps0 ){
							eY=1.0*ec0;
							sigY = sigPY + 1.0*ec0 * depsY;
						}
					}
				//if (eptY==0), linear
				} else {
					eY = ec0;
					sigY = eY * (epsY - eptY);
				}
				
				if ( ecminY == 0 ) {
					// 1st iteration
					eY = sicnY / deptY;
					sigY = eY * (epsY - eptY);
				}
				
			}else {
			/* <Region #4 >--- 1 ---- ecmin ---- 2 ---- ept ---- 3 ---- epn ---- 4  */
				
				// else, if the current strain is larger than epn the response
				// corresponds to the tensile envelope curve shifted by ept
				
				double epstmpY = epsY - eptY;
				this->Tens_Envlp(epstmpY, epsPY, Tens_Alpha2, sigY, eY);
				deptY = epsY - eptY;
			}
		}
	}
	
	// END OF r-dir
	
	
	// <CONCRETE SHEAR in d-r coord> ---------------------------------------------------
	
	/* SW4 */
	tauXY = 2.0*fc/epsc0*ShearBeta*gammaXY/1.2/2.0;
	
	// <Cracking Criteria>  ------------------------------------------------------------
	double Ec0 = 2.0*fc/epsc0;
	double eps0 = ft/Ec0;
	/* 1. principal stress dir*/
	double thetaP = (atan(2.*tauXY/(sigX-sigY)))/2.;
	beta = thetaP - thetaSP;
	
	/*2. 1st crack dir*/
	if(eps_p(0)>eps0 && thetaSP==0){
		thetaSP = (atan(2.*tauXY/(sigX-sigY)))/2.;
	}
	/*3. initial shear stiffness*/
	if(thetaSP==0){
		tauXY = 2.0*fc/epsc0*gammaXY/1.00/2.0;
		//tauXY = 2.0*fc/epsc0*gammaXY/1.20/2.0;
	}
	
	//<STEEL>----------------------------------------------------------------------------
	epsL = strain_from_element(0);
	epsT = strain_from_element(1);
	gammaLT = strain_from_element(2);
	
	//L dir//
	if (TloadIndicatorL == 0 && epsL == 0.0)
		return 0;
	
	TrotMaxL = CrotMaxL;
	TrotMinL = CrotMinL;
	TenergyDL = CenergyDL;
	TrotPuL = CrotPuL;
	TrotNuL = CrotNuL;
	
	TstrainL = epsL;
	double dStrainL = TstrainL - CstrainL;
	TloadIndicatorL = CloadIndicatorL;
	
	if (TloadIndicatorL == 0)
		TloadIndicatorL = (dStrainL < 0.0) ? 2:1;
	
	if (TstrainL >= CrotMaxL) {
		TrotMaxL = TstrainL;
		TtangentL = posEnvlpTangent(TstrainL, ReducedL);
		TstressL = posEnvlpStress(TstrainL, ReducedL);
	}
	else if (TstrainL <= CrotMinL) {
		TrotMinL = TstrainL;
		TtangentL = negEnvlpTangent(TstrainL, ReducedL);
		TstressL = negEnvlpStress(TstrainL, ReducedL);
	}
	else {
		if (dStrainL < 0.0)
			negativeIncrementL(dStrainL, ReducedL);
		else if (dStrainL > 0.0)
			positiveIncrementL(dStrainL, ReducedL);
	}
	
	TenergyDL = CenergyDL + 0.5*(CstressL+TstressL)*dStrainL;
	//end of L dir//
	
	
	//T dir//
	if (TloadIndicatorT == 0 && epsT == 0.0)
		return 0;
	
	TrotMaxT = CrotMaxT;
	TrotMinT = CrotMinT;
	TenergyDT = CenergyDT;
	TrotPuT = CrotPuT;
	TrotNuT = CrotNuT;
	TstrainT = epsT;
	double dStrainT = TstrainT - CstrainT;
	TloadIndicatorT = CloadIndicatorT;
	
	if (TloadIndicatorT == 0)
		TloadIndicatorT = (dStrainT < 0.0) ? 2:1;
	
	if (TstrainT >= CrotMaxT) {
		TrotMaxT = TstrainT;
		TtangentT = posEnvlpTangent(TstrainT, ReducedT);
		TstressT = posEnvlpStress(TstrainT, ReducedT);
	}
	else if (TstrainT <= CrotMinT) {
		TrotMinT = TstrainT;
		TtangentT = negEnvlpTangent(TstrainT, ReducedT);
		TstressT = negEnvlpStress(TstrainT, ReducedT);
	}
	else {
		if (dStrainT < 0.0)
			negativeIncrementT(dStrainT, ReducedT);
		else if (dStrainT > 0.0)
			positiveIncrementT(dStrainT, ReducedT);
	}
	
	TenergyDT = CenergyDT + 0.5*(CstressT+TstressT)*dStrainT;
	//end of T dir//
	
	return 0;
}

const Matrix&
NonlinearBS::getTangent(void)
{
	dsigXdeX = eX;
	dsigXdeY = 0.;
	dsigYdeX = 0.;
	dsigYdeY = eY;
	/* SW4 */
	dtauXYdgXY = 2.0*fc/epsc0*ShearBeta/1.2/2.0;
	/* Initial shear stiffness */
	if(thetaSP==0){
		dtauXYdgXY = 2.0*fc/epsc0/1.00/2.0;
		//dtauXYdgXY = 2.0*fc/epsc0/1.20/2.0;
	}
	
	Dc(0,0) = dsigXdeX;
	Dc(0,1) = dsigXdeY;
	Dc(0,2) = 0.;
	Dc(1,0) = dsigYdeX;
	Dc(1,1) = dsigYdeY;
	Dc(1,2) = 0.;
	Dc(2,0) = 0.;
	Dc(2,1) = 0.;
	Dc(2,2) = dtauXYdgXY;
	
	Ds(0,0) = reinfR1*TtangentL;
	Ds(0,1) = 0.;
	Ds(0,2) = 0.;
	Ds(1,0) = 0.;
	Ds(1,1) = reinfR2*TtangentT;
	Ds(1,2) = 0.;
	Ds(2,0) = 0.;
	Ds(2,1) = 0.;
	Ds(2,2) = 0.; //reinfR1*0.0001 *mom1p/rot1p;
	
	//Dxy = T1*Dc*T + Ds;
	Dxy.addMatrixTripleProduct(0.0, T, Dc, 1.0);
	tangent_matrix = Dxy + Ds;
	
	return tangent_matrix;
}

const Vector&
NonlinearBS::getStrain(void)
{
	strain_vec(0) = eOXX;
	strain_vec(1) = eOYY;
	strain_vec(2) = gOXY;
	
	return strain_vec;
}

const Vector&
NonlinearBS::getStress(void)
{
	sigL = reinfR1 * TstressL;
	sigT = reinfR2 * TstressT;
	tauLT = 0.;
	
	/*conc+steel*/
	sigOXX = (sigX*pow(cos(thetaC),2.)+sigY*pow(sin(thetaC),2.)+tauXY*sin(thetaC)*cos(thetaC))+sigL;
	sigOYY = sigY*pow(cos(thetaC),2.)-sigX*pow(sin(thetaC),2.)+tauXY*sin(thetaC)*cos(thetaC)+sigT;
	tauOXY = (-sigY+sigX)*sin(thetaC)*cos(thetaC)+tauXY*(pow(cos(thetaC),2.)-pow(sin(thetaC),2.));
	
	stress_vec(0) = sigOXX;
	stress_vec(1) = sigOYY;
	stress_vec(2) = tauOXY;
	
	return stress_vec;
}

int
NonlinearBS::commitState(void)
{
	ecminPX = ecminX;
	deptPX = deptX;
	
	ecminPY = ecminY;
	deptPY = deptY;
	
	ePX = eX;
	sigPX = sigX;
	epsPX = epsX;
	
	ePY = eY;
	sigPY = sigY;
	epsPY = epsY;
	
	CrotMaxL = TrotMaxL;
	CrotMinL = TrotMinL;
	CrotPuL = TrotPuL;
	CrotNuL = TrotNuL;
	CenergyDL = TenergyDL;
	CloadIndicatorL = TloadIndicatorL;
	
	CstressL = TstressL;
	CstrainL = TstrainL;
	
	CrotMaxT = TrotMaxT;
	CrotMinT = TrotMinT;
	CrotPuT = TrotPuT;
	CrotNuT = TrotNuT;
	CenergyDT = TenergyDT;
	CloadIndicatorT = TloadIndicatorT;
	
	CstressT = TstressT;
	CstrainT = TstrainT;
	return 0;
}

int
NonlinearBS::revertToLastCommit(void)
{
	ecminX = ecminPX;
	deptX = deptPX;
	
	ecminY = ecminPY;
	deptY = deptPY;
	
	eX = ePX;
	sigX = sigPX;
	epsX = epsPX;
	
	eY = ePY;
	sigY = sigPY;
	epsY = epsPY;
	
	TrotMaxL = CrotMaxL;
	TrotMinL = CrotMinL;
	TrotPuL = CrotPuL;
	TrotNuL = CrotNuL;
	TenergyDL = CenergyDL;
	TloadIndicatorL = CloadIndicatorL;
	
	TstressL = CstressL;
	TstrainL = CstrainL;
	
	TrotMaxT = CrotMaxT;
	TrotMinT = CrotMinT;
	TrotPuT = CrotPuT;
	TrotNuT = CrotNuT;
	TenergyDT = CenergyDT;
	TloadIndicatorT = CloadIndicatorT;
	
	TstressT = CstressT;
	TstrainT = CstrainT;
	
	return 0;
}

int
NonlinearBS::revertToStart(void)
{
	ecminPX = 0.0;
	deptPX = 0.0;
	ecminPY = 0.0;
	deptPY = 0.0;
	
	ePX = 2.0*fc/epsc0;
	epsPX = 0.0;
	sigPX = 0.0;
	epsX = 0.0;
	sigX = 0.0;
	eX = 2.0*fc/epsc0;

	ePY = 2.0*fc/epsc0;
	epsPY = 0.0;
	sigPY = 0.0;
	epsY = 0.0;
	sigY = 0.0;
	eY = 2.0*fc/epsc0;
	
	CrotMaxL = 0.0;
	CrotMinL = 0.0;
	CrotPuL = 0.0;
	CrotNuL = 0.0;
	CenergyDL = 0.0;
	CloadIndicatorL = 0;
	
	CstressL = 0.0;
	CstrainL = 0.0;
	
	TstrainL = 0;
	TstressL = 0;
	TtangentL = E1p;
	
	CrotMaxT = 0.0;
	CrotMinT = 0.0;
	CrotPuT = 0.0;
	CrotNuT = 0.0;
	CenergyDT = 0.0;
	CloadIndicatorT = 0;
	
	CstressT = 0.0;
	CstrainT = 0.0;
	
	TstrainT = 0;
	TstressT = 0;
	TtangentT = E1p;
	
	return 0;
}

int
NonlinearBS::sendSelf(int commitTag, Channel &theChannel)
{
	static Vector data(60);
	data(0) = this->getTag();
	data(1)	= fc;
	data(2) = epsc0;
	data(3) = fcu;
	data(4) = epscu;
	data(5) = rat;
	data(6) = ft;
	data(7) = Ets;
	data(8) = reinfD1;
	data(9) = reinfR1;
	data(10) = reinf_dmm1;
	data(11) = reinfD2;
	data(12) = reinfR2;
	data(13) = reinf_dmm2;
	data(14) = ecminPX;
	data(15) = deptPX;
	data(16) = epsPX;
	data(17) = sigPX;
	data(18) = ePX;
	data(19) = ecminPY;
	data(20) = deptPY;
	data(21) = epsPY;
	data(22) = sigPY;
	data(23) = ePY;
	data(24) = thetaSP;
	data(25) = mom1p;
	data(26) = rot1p;
	data(27) = mom2p;
	data(28) = rot2p;
	data(29) = mom3p;
	data(30) = rot3p;
	data(31) = mom1n;
	data(32) = rot1n;
	data(33) = mom2n;
	data(34) = rot2n;
	data(35) = mom3n;
	data(36) = rot3n;
	data(37) = pinchX;
	data(38) = pinchY;
	data(39) = damfc1;
	data(40) = damfc2;
	data(41) = beta;
	data(42) = CrotMaxL;
	data(43) = CrotMinL;
	data(44) = CrotPuL;
	data(45) = CrotNuL;
	data(46) = CenergyDL;
	data(47) = CloadIndicatorL;
	data(48) = CstressL;
	data(49) = CstrainL;
	data(50) = TtangentL;
	data(51) = CrotMaxT;
	data(52) = CrotMinT;
	data(53) = CrotPuT;
	data(54) = CrotNuT;
	data(55) = CenergyDT;
	data(56) = CloadIndicatorT;
	data(57) = CstressT;
	data(58) = CstrainT;
	data(59) = TtangentT;

	if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "NonlinearBS::sendSelf() - failed to sendSelf\n";
		return -1;
	}
	return 0;
}

int
NonlinearBS::recvSelf(int commitTag, Channel &theChannel,
					FEM_ObjectBroker &theBroker)
{
	static Vector data(60);
	if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "NonlinearBS::recvSelf() - failed to recvSelf\n";
		return -1;
	}
	
	this->setTag(int(data(0)));
	fc              = data(1);
	epsc0           = data(2);
	fcu             = data(3);
	epscu           = data(4);
	rat             = data(5);
	ft              = data(6);
	Ets             = data(7);
	reinfD1         = data(8);
	reinfR1         = data(9);
	reinf_dmm1      = data(10);
	reinfD2         = data(11);
	reinfR2         = data(12);
	reinf_dmm2      = data(13);
	ecminPX         = data(14);
	deptPX          = data(15);
	epsPX           = data(16);
	sigPX           = data(17);
	ePX             = data(18);
	ecminPY         = data(19);
	deptPY          = data(20);
	epsPY           = data(21);
	sigPY           = data(22);
	ePY             = data(23);
	thetaSP         = data(24);
	mom1p           = data(25);
	rot1p           = data(26);
	mom2p           = data(27);
	rot2p           = data(28);
	mom3p           = data(29);
	rot3p           = data(30);
	mom1n           = data(31);
	rot1n           = data(32);
	mom2n           = data(33);
	rot2n           = data(34);
	mom3n           = data(35);
	rot3n           = data(36);
	pinchX          = data(37);
	pinchY          = data(38);
	damfc1          = data(39);
	damfc2          = data(40);
	beta            = data(41);
	CrotMaxL        = data(42);
	CrotMinL        = data(43);
	CrotPuL         = data(44);
	CrotNuL         = data(45);
	CenergyDL       = data(46);
	CloadIndicatorL = int(data(47));
	CstressL        = data(48);
	CstrainL        = data(49);
	TtangentL       = data(50);
	CrotMaxT        = data(51);
	CrotMinT        = data(52);
	CrotPuT         = data(53);
	CrotNuT         = data(54);
	CenergyDT       = data(55);
	CloadIndicatorT = int(data(56));
	CstressT        = data(57);
	CstrainT        = data(58);
	TtangentT       = data(59);
	
	eX = ePX;
	sigX = sigPX;
	epsX = epsPX;
	
	eY = ePY;
	sigY = sigPY;
	epsY = epsPY;
	
	// set the trial values
	TrotMaxL = CrotMaxL;
	TrotMinL = CrotMinL;
	TrotPuL = CrotPuL;
	TrotNuL = CrotNuL;
	TenergyDL = CenergyDL;
	TloadIndicatorL = CloadIndicatorL;
	TstressL = CstressL;
	TstrainL = CstrainL;
	
	TrotMaxT = CrotMaxT;
	TrotMinT = CrotMinT;
	TrotPuT = CrotPuT;
	TrotNuT = CrotNuT;
	TenergyDT = CenergyDT;
	TloadIndicatorT = CloadIndicatorT;
	TstressT = CstressT;
	TstrainT = CstrainT;
	
	return 0;
}

void
NonlinearBS::Print(OPS_Stream &s, int flag)
{
    s << "NonlinearBS: strain:\t" << strain_vec << endln;
    s << "             stress:\t" << stress_vec << endln;
    s << "            tangent:\t" << tangent_matrix << endln;
}

void
NonlinearBS::Tens_Envlp (double epsc, double epscP, double Alpha, double &sigc, double &Ect)
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
    double Lamda = min(270./sqrt(Alpha),1000.);
	double Pwr = -1.*Lamda * (epsc - eps0);
	sigc = ft*(1-Alpha)*exp(Pwr)+Alpha;
	double PwrPrime = -1.*Lamda;
	Ect = ft*(1-Alpha)*PwrPrime*exp(Pwr);
  }
  return;
}


void
NonlinearBS::Compr_Envlp (double epsc, double epsc2, double &sigc, double &Ect)
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
      Ect  = 1.0e-10;
      //       Ect  = 0.0
    }
  }
  return;
}

void
NonlinearBS::setEnv(void)
{
	E1p = mom1p/rot1p;
	E2p = (mom2p-mom1p)/(rot2p-rot1p);
	E3p = (mom3p-mom2p)/(rot3p-rot2p);
	E1n = mom1n/rot1n;
	E2n = (mom2n-mom1n)/(rot2n-rot1n);
	E3n = (Rmom3n-mom2n)/(rot3n-rot2n);
}

double
NonlinearBS::setEnvelope(double Reduced)
{
	Rmom1p = Reduced*mom1p; Rrot1p = Reduced*rot1p;
	Rmom2p = Reduced*mom2p; Rrot2p = rot2p;
	Rmom3p = Reduced*mom3p; Rrot3p = rot3p;
	Rmom1n = Reduced*mom1n; Rrot1n = Reduced*rot1n;
	Rmom2n = Reduced*mom2n; Rrot2n = rot2n;
	Rmom3n = Reduced*mom3n; Rrot3n = rot3n;
	rE1p = Rmom1p/Rrot1p;
	rE2p = (Rmom2p-Rmom1p)/(Rrot2p-Rrot1p);
	rE3p = (Rmom3p-Rmom2p)/(Rrot3p-Rrot2p);
	rE1n = Rmom1n/Rrot1n;
	rE2n = (Rmom2n-Rmom1n)/(Rrot2n-Rrot1n);
	rE3n = (Rmom3n-Rmom2n)/(Rrot3n-Rrot2n);
	return Rmom1p, Rmom2p, Rmom3p, Rmom1n, Rmom2n, Rmom3n, Rrot1p,
			Rrot2p, Rrot3p, Rrot1n, Rrot2n, Rrot3n, rE1p, rE2p, rE2p, rE1n, rE2n, rE2n;
}

double
NonlinearBS::posEnvlpStress(double strain, double Reduced)
{
	this->setEnvelope(Reduced);
	if(strain<=0.0)
		return 0.0;
	else if (strain <= Rrot1p)
		return rE1p*strain;
	else if (strain <= Rrot2p)
		return Rmom1p + rE2p*(strain-Rrot1p);
	else if (strain <= Rrot3p || rE3p > 0.0)
		return Rmom2p + rE3p*(strain-Rrot2p);
	else
		return Rmom3p;
}

double
NonlinearBS::negEnvlpStress(double strain, double Reduced)
{
	this->setEnvelope(Reduced);
	if (strain >= 0.0)
	return 0.0;
	else if (strain >= Rrot1n)
		return rE1n*strain;
	else if (strain >= Rrot2n)
		return Rmom1n + rE2n*(strain-Rrot1n);
	else if (strain >= Rrot3n || rE3n > 0.0)
		return Rmom2n + rE3n*(strain-Rrot2n);
	else
		return Rmom3n;
}

double
NonlinearBS::posEnvlpTangent(double strain, double Reduced)
{
	this->setEnvelope(Reduced);
	if(strain<0.0)
		return rE1p*1.0e-9;
	else if (strain <= Rrot1p)
		return rE1p;
	else if (strain <= Rrot2p)
		return rE2p;
	else if (strain <= Rrot3p || rE3p > 0.0)
		return rE3p;
	else
		return rE1p* 1.0e-9;
}

double
NonlinearBS::negEnvlpTangent(double strain, double Reduced)
{
	this->setEnvelope(Reduced);
	if (strain > 0.0)
		return rE1n*1.0e-9;
	else if (strain >= Rrot1n)
		return rE1n;
	else if (strain >= Rrot2n)
		return rE2n;
	else if (strain >= Rrot3n || rE3n > 0.0)
		return rE3n;
	else
		return rE1n*1.0e-9;
}

double
NonlinearBS::posEnvlpRotlim(double strain, double Reduced)
{
	this->setEnvelope(Reduced);
	double strainLimit = POS_INF_STRAIN;
	
	if(strain<=Rrot1p)
		return POS_INF_STRAIN;
	if (strain > Rrot1p && strain <= Rrot2p && rE2p < 0.0)
		strainLimit = Rrot1p - Rmom1p/rE2p;
	if (strain > Rrot2p && rE3p < 0.0)
		strainLimit = Rrot2p - Rmom2p/rE3p;

	if (strainLimit == POS_INF_STRAIN)
		return POS_INF_STRAIN;
	else if (posEnvlpStress(strainLimit, Reduced) > 0)
		return POS_INF_STRAIN;
	else
		return strainLimit;
}

double
NonlinearBS::negEnvlpRotlim(double strain, double Reduced)
{
	this->setEnvelope(Reduced);
	
	double strainLimit = NEG_INF_STRAIN;
	
	if (strain >=Rrot1n)
		return NEG_INF_STRAIN;
	if (strain < Rrot1n && strain >= Rrot2n && rE2n < 0.0)
		strainLimit = Rrot1n - Rmom1n/rE2n;
	if (strain < Rrot2n && rE3n < 0.0)
		strainLimit = Rrot2n - Rmom2n/rE3n;

	if (strainLimit == NEG_INF_STRAIN)
		return NEG_INF_STRAIN;
	else if (negEnvlpStress(strainLimit, Reduced) < 0)
		return NEG_INF_STRAIN;
	else
		return strainLimit;
}

void
NonlinearBS::positiveIncrementL(double dStrain, double Reduced)
{
	this->setEnvelope(Reduced);
	
	double RenergyA = 0.5 * (Rrot1p*Rmom1p + (Rrot2p-Rrot1p)*(Rmom2p+Rmom1p) + (Rrot3p-Rrot2p)*(Rmom3p+Rmom2p) +
		Rrot1n*Rmom1n + (Rrot2n-Rrot1n)*(Rmom2n+Rmom1n) + (Rrot3n-Rrot2n)*(Rmom3n+Rmom2n));
	
	double kn = pow(CrotMinL/Rrot1n, beta);
	kn = (kn<1.0)?1.0:1.0/kn;
	double kp = pow(CrotMaxL/Rrot1p, beta);
	kp = (kp<1.0)?1.0:1.0/kp;
	
	if (TloadIndicatorL == 2) {
		TloadIndicatorL = 1;
		if (CstressL <= 0.0) {
			TrotNuL = CstrainL - CstressL/(rE1n*kn);
			double energy = CenergyDL - 0.5*CstressL/(rE1n*kn)*CstressL;
			double damfc = 0.0;
			if(CrotMinL<Rrot1n){
				damfc = damfc2*energy / RenergyA;
				damfc += damfc1*(CrotMinL-Rrot1n)/Rrot1n;
			}
			TrotMaxL = CrotMaxL*(1.0+damfc);
		}
	}
	
	TloadIndicatorL = 1;
	TrotMaxL = (TrotMaxL > Rrot1p) ? TrotMaxL : Rrot1p;
	
	double maxmom = posEnvlpStress(TrotMaxL, Reduced);
	double rotlim = negEnvlpRotlim(CrotMinL, Reduced);
	double rotrel = (rotlim > TrotNuL) ? rotlim : TrotNuL;
	
	double rotmp2 = TrotMaxL - (1.0-pinchY)*maxmom/(rE1p*kp);
	double rotch = rotrel + (rotmp2-rotrel)*pinchX;
	double tmpmo1;
	double tmpmo2;
	
	if (TstrainL < TrotNuL) {
		TtangentL = rE1n*kn;
		TstressL = CstressL + TtangentL*dStrain;
		if(TstressL>=0.0){
			TstressL = 0.0;
			TtangentL = rE1n* 1.0e-9;
		}
	}
	else if (TstrainL >= TrotNuL && TstrainL < rotch) {
		if (TstrainL <= rotrel) {
			TstressL = 0.0;
			TtangentL = rE1p*1.0e-9;
		}
		else {
			TtangentL = maxmom*pinchY/(rotch-rotrel);
			tmpmo1 = CstressL + rE1p*kp*dStrain;
			tmpmo2 = (TstrainL-rotrel)*TtangentL;
			if (tmpmo1 < tmpmo2) {
				TstressL = tmpmo1;
				TtangentL = rE1p*kp;
			}
			else
				TstressL = tmpmo2;
		}
	}
	else {
		TtangentL = (1.0-pinchY)*maxmom/(TrotMaxL-rotch);
		tmpmo1 = CstressL + rE1p*kp*dStrain;
		tmpmo2 = pinchY*maxmom + (TstrainL-rotch)*TtangentL;
		if (tmpmo1 < tmpmo2) {
			TstressL = tmpmo1;
			TtangentL = rE1p*kp;
		}
		else
			TstressL = tmpmo2;
	}
}

void
NonlinearBS::positiveIncrementT(double dStrain, double Reduced)
{
	this->setEnvelope(Reduced);
	double RenergyA = 0.5 * (Rrot1p*Rmom1p + (Rrot2p-Rrot1p)*(Rmom2p+Rmom1p) + (Rrot3p-Rrot2p)*(Rmom3p+Rmom2p) +
		Rrot1n*Rmom1n + (Rrot2n-Rrot1n)*(Rmom2n+Rmom1n) + (Rrot3n-Rrot2n)*(Rmom3n+Rmom2n));
	
	double kn = pow(CrotMinT/Rrot1n,beta);
	kn = (kn < 1.0) ? 1.0: 1.0/kn;
	double kp = pow(CrotMaxT/Rrot1p,beta);
	kp = (kp < 1.0) ? 1.0: 1.0/kp;
	
	if (TloadIndicatorT == 2) {
		TloadIndicatorT = 1;
		if(CstressT<=0.0){
			TrotNuT = CstrainT- CstressT/(rE1n*kn);
			double energy = CenergyDT - 0.5*CstressT/(rE1n*kn)*CstressT;
			double damfc = 0.0;
			if(CrotMinT<Rrot1n) {
				damfc = damfc2 * energy / RenergyA;
				damfc += damfc1*(CrotMinT-Rrot1n)/Rrot1n;
			}
			TrotMaxT = CrotMaxT*( 1.0+damfc);
		}
	}
	
	TloadIndicatorT = 1;
	
	TrotMaxT = (TrotMaxT > Rrot1p) ? TrotMaxT : Rrot1p;

	double maxmom = posEnvlpStress(TrotMaxT, Reduced);
	double rotlim = negEnvlpRotlim(CrotMinT, Reduced);
	double rotrel = (rotlim > TrotNuT) ? rotlim : TrotNuT;
	double rotmp2 = TrotMaxT - (1.0-pinchY)*maxmom/(rE1p*kp);
	double rotch = rotrel + (rotmp2-rotrel)*pinchX;

	double tmpmo1;
	double tmpmo2;
	
	if (TstrainT < TrotNuT) {
		TtangentT = rE1n*kn;
		TstressT = CstressT + TtangentT*dStrain;
		if(TstressT>=0.0){
			TstressT = 0.0;
			TtangentT = rE1n*1.0e-9;
		}
	}
	else if (TstrainT >= TrotNuT && TstrainT < rotch) {
		if (TstrainT <= rotrel) {
			TstressT = 0.0;
			TtangentT =rE1p*1.0e-9;
		} else {
			TtangentT = maxmom*pinchY/(rotch-rotrel);
			tmpmo1 = CstressT + rE1p*kp*dStrain;
			tmpmo2 = (TstrainT-rotrel)*TtangentT;
			if (tmpmo1 <tmpmo2) {
				TstressT = tmpmo1;
				TtangentT = rE1p*kp;
			}
			else
				TstressT = tmpmo2;
		}
	}
	else {
		TtangentT = (1.0-pinchY)*maxmom/(TrotMaxT-rotch);
		tmpmo1 = CstressT + rE1p*kp*dStrain;
		tmpmo2 = pinchY*maxmom + (TstrainT-rotch)*TtangentT;
		if (tmpmo1 < tmpmo2) {
			TstressT = tmpmo1;
			TtangentT = rE1p*kp;
		}
		else
			TstressT = tmpmo2;
	}
}

void
NonlinearBS::negativeIncrementL(double dStrain, double Reduced)
{
	this->setEnvelope(Reduced);
	double RenergyA = 0.5 * (Rrot1p*Rmom1p + (Rrot2p-Rrot1p)*(Rmom2p+Rmom1p) + (Rrot3p-Rrot2p)*(Rmom3p+Rmom2p) +
		Rrot1n*Rmom1n + (Rrot2n-Rrot1n)*(Rmom2n+Rmom1n) + (Rrot3n-Rrot2n)*(Rmom3n+Rmom2n));
	
	double kn = pow(CrotMinL/Rrot1n,beta);
	kn = (kn<1.0) ? 1.0 : 1.0/kn;
	double kp = pow(CrotMaxL/Rrot1p,beta);
	kp = (kp<1.0) ? 1.0 : 1.0/kp;

	if(TloadIndicatorL == 1) {
		TloadIndicatorL = 2;
		if(CstressL>=0.0){
			TrotPuL = CstrainL - CstressL/(rE1p*kp);
			double energy = CenergyDL - 0.5*CstressL/(rE1p*kp)*CstressL;
			double damfc = 0.0;
			if(CrotMaxL>Rrot1p){
				damfc = damfc2*energy/RenergyA;
				damfc += damfc1*(CrotMaxL-Rrot1p)/Rrot1p;
			}
			TrotMinL = CrotMinL*(1.0+damfc);
		}
	}
	TloadIndicatorL = 2;
	TrotMinL = (TrotMinL < Rrot1n) ? TrotMinL : Rrot1n;
	
	double minmom = negEnvlpStress(TrotMinL, Reduced);
	double rotlim = posEnvlpRotlim(CrotMaxL, Reduced);
	double rotrel = (rotlim < TrotPuL) ? rotlim : TrotPuL;
	double rotmp2 = TrotMinL - (1.0-pinchY)*minmom/(rE1n*kn);
	double rotch = rotrel + (rotmp2-rotrel)*pinchX;

	double tmpmo1;
	double tmpmo2;
	
	if(TstrainL> TrotPuL) {
		TtangentL = rE1p*kp;
		TstressL = CstressL + TtangentL*dStrain;
		if(TstressL<=0.0) {
			TstressL = 0.0;
			TtangentL = rE1p*1.0e-9;
		}
	}
	else if (TstrainL <= TrotPuL && TstrainL > rotch) {
		if (TstrainL >= rotrel) {
			TstressL = 0.0;
			TtangentL = rE1n*1.0e-9;
		}
		else {
			TtangentL = minmom*pinchY/(rotch-rotrel);
			tmpmo1 = CstressL + rE1n*kn*dStrain;
			tmpmo2 = (TstrainL-rotrel)*TtangentL;
			if (tmpmo1 > tmpmo2) {
				TstressL = tmpmo1;
				TtangentL = rE1n*kn;
			}
			else
				TstressL = tmpmo2;
		}
	}
	else {
		TtangentL = (1.0-pinchY)*minmom/(TrotMinL-rotch);
		tmpmo1 = CstressL + rE1n*kn*dStrain;
		tmpmo2 = pinchY*minmom + (TstrainL-rotch)*TtangentL;
		if (tmpmo1 > tmpmo2) {
			TstressL = tmpmo1;
			TtangentL = rE1n*kn;
		}
		else
			TstressL = tmpmo2;
	}
}

void
NonlinearBS::negativeIncrementT(double dStrain, double Reduced)
{
	this->setEnvelope(Reduced);
	double RenergyA = 0.5 * (Rrot1p*Rmom1p + (Rrot2p-Rrot1p)*(Rmom2p+Rmom1p) + (Rrot3p-Rrot2p)*(Rmom3p+Rmom2p) +
		Rrot1n*Rmom1n + (Rrot2n-Rrot1n)*(Rmom2n+Rmom1n) + (Rrot3n-Rrot2n)*(Rmom3n+Rmom2n));
	
	double kn = pow(CrotMinT/Rrot1n,beta);
	kn = (kn<1.0) ? 1.0 : 1.0/kn;
	double kp = pow(CrotMaxT/Rrot1p,beta);
	kp = (kp<1.0) ? 1.0 : 1.0/kp;
	
	if (TloadIndicatorT == 1) {
		TloadIndicatorT = 2;
		if(CstressT>=0.0){
			TrotPuT = CstrainT - CstressT/(rE1p*kp);
			double energy = CenergyDT - 0.5*CstressT/(rE1p*kp)*CstressT;
			double damfc = 0.0;
			if(CrotMaxT>Rrot1p){
				damfc = damfc2*energy/RenergyA;
				damfc += damfc1*(CrotMaxT-Rrot1p)/Rrot1p;
			}
			TrotMinT = CrotMinT*(1.0+damfc);
		}
	}
	
	TloadIndicatorT = 2;
	TrotMinT = (TrotMinT < Rrot1n) ? TrotMinT : Rrot1n;
	
	double minmom = negEnvlpStress(TrotMinT, Reduced);
	double rotlim = posEnvlpRotlim(CrotMaxT, Reduced);
	double rotrel = (rotlim < TrotPuT) ? rotlim : TrotPuT;
	double rotmp2 = TrotMinT - (1.0-pinchY)*minmom/(rE1n*kn);
	double rotch = rotrel + (rotmp2-rotrel)*pinchX; // changed on 7/11/2006

	double tmpmo1;
	double tmpmo2;
	
	if(TstrainT > TrotPuT) {
		TtangentT = rE1p*kp;
		TstressT = CstressT + TtangentT*dStrain;
		if (TstressT <= 0.0) {
			TstressT = 0.0;
			TtangentT = rE1p*1.0e-9;
		}
	}
	else if (TstrainT <= TrotPuT && TstrainT > rotch) {
		if (TstrainT >= rotrel) {
			TstressT = 0.0;
			TtangentT = rE1n*1.0e-9;
		}
		else {
			TtangentT = minmom*pinchY/(rotch-rotrel);
			tmpmo1 = CstressT + rE1n*kn*dStrain;
			tmpmo2 = (TstrainT-rotrel)*TtangentT;
			if (tmpmo1 > tmpmo2) {
				TstressT = tmpmo1;
				TtangentT = rE1n*kn;
			}
			else
				TstressT = tmpmo2;
		}
	}
	else {
		TtangentT = (1.0-pinchY)*minmom/(TrotMinT-rotch);
		tmpmo1 = CstressT + rE1n*kn*dStrain;
		tmpmo2 = pinchY*minmom + (TstrainT-rotch)*TtangentT;
		if (tmpmo1 >tmpmo2) {
			TstressT = tmpmo1;
			TtangentT = rE1n*kn;
		}
		else
			TstressT = tmpmo2;
	}
}

Response*
NonlinearBS::setResponse (const char **argv, int argc, Information &matInfo, OPS_Stream &output)
{
	//cout << "Saving stresses and strains" << endln;
	
	if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());
	
	else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
	
	else if (strcmp(argv[0],"tangent") == 0)
		return new MaterialResponse(this, 3, this->getTangent());
	
	return 0;
}

int
NonlinearBS::getResponse (int responseID, Information &matInfo)
{
	switch (responseID) {
	case -1:
		return-1;
	case 1:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStress();
		return 0;
	case 2:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStrain();
		return 0;
	case 3:
		if (matInfo.theMatrix != 0)
			*(matInfo.theMatrix) = getTangent();
		return 0;
	default:
		return -1;
	}
}

void
NonlinearBS::Tens_UnloadingP (double epsc, double epsc2, double epscP, double sigcP, double ecmin, double dept, double friction,
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
				//if ( abs(Sige0) < 10e-12) {
				//	AA = 0.0;
				//	BB = 0.0;
				//	CC = 0.0;
				//} else {
					CC = Sige0;
					BB = ((-friction) + (-sigFF/M)/M - Sige0 + Sige0/M/M)/(eSlip0+eSlip0/M);
					AA = ((-friction) - BB*eSlip0 - Sige0)/eSlip0/eSlip0;
				//}

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
NonlinearBS::Tens_UnloadingN (double epsc, double epsc2, double epscP, double sigcP, double ecmin, double dept, double friction,
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
		if (sigc <=0) {
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
		//if ( abs(Sige0) < 10e-12) {
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
NonlinearBS::Compr_Reloading (double epsc, double epsc2, double epscP, double sigcP, double ecmin, double dept, double friction,
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
	//if ( abs(Sige0) < 10e-12) {
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
		double epn = ept + dept;
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
