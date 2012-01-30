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
// $Date: 2008-11-06 20:49:49 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/uniaxial/CCFT_concMaterial.cpp,v $

// Written: MHS
// Created: Aug 2001

// Modified by : Cenk Tort - tort0008@umn.edu
// University of Minnesota - Civil Engineering Department
// Date: Wed Jul 23 17:40:25 EDT 2003
// Description: This file contains the class implementation for
// cyclic uniaxial stress-strain relationship of concrete for
// RCFT members.
// Further Modified by: Mark Denavit
// University of Illinois at Urbana-Champaign
// for circular CFT (CCFT) members

#include "CCFT_concMaterial.h"
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <fstream>
#include <iomanip>

#define SMALL_STRAIN 1.0e-20
//#define SMALL_STRESS 1.0e-16
#define SMALL_TANGENT 1.0e-7
#define SMALL_STRESS 0.0
//#define SMALL_TANGENT 0.0
#define SMALL_NUMBER 1.0e-10

//using std::ofstream;
//using std::ios;
//using std::endl;
using namespace std;

CCFT_concMaterial::CCFT_concMaterial( int tag, double i1, double i2, double i3, double i4, double i5 )
  :UniaxialMaterial(tag, MAT_TAG_CCFT_conc), 
  Fc(i1), D(i2), t(i3), Fy(i4), Es(i5)
{
	// Make sure Fc is positive
	if ( Fc < 0.0 ) {
		Fc = -Fc;
	}
	
	// Define Backbone Curve Parameters
	//Ec = 40.0 * sqrt (Fc*1000) + 1000.0;
	Ec = 185.0 * pow(1000*Fc,0.375);
	
	// Compressive Branch
	double alphaHoop = 0.134; //calibrated parameter (see thesis)
	double confinementPressure = 2*alphaHoop*Fy/(D/t-1);
	double x = confinementPressure/Fc; //ratio of confiment pressure to unconfined compressive strength
	
	Fc_n = -Fc*(-1.254+2.254*sqrt(1+7.94*x)-2*x); //confined compressive strength

	double k2 = 5*(-Fc_n-Fc)/confinementPressure;
	double ec_unconf = -1*pow(1000*Fc,0.25)/4000; //unconfined strain at peak stress
    ec_n = ec_unconf*(1+k2*x); //confined strain at peak stress
    
	double m = (1000*Fc - 36) / (6*sqrt(1000*Fc));
	Fc_res_n = -Fc*(sqrt(m*x)+x);
	if (Fc_res_n < Fc_n) {
		Fc_res_n = Fc_n;  // Residual strength is not to be larger than confined strength
	}
    
	double K = -1*Fc_n/Fc;
	double n3 = 46.0 / pow(1000*Fc,0.375);
	double r3 = 1000*Fc/750 - 1.9;
	double y3, scratch;
	tsaiEquation(3.0,r3,n3,y3,scratch);	
	double dfc = Fc*(1-y3);
	Ec_sft = K*dfc*( 0.8/pow(K,5) + 0.2 ) / (2*ec_n); //softening slope of the linear decending branch
		
	ec_res_n = ( Fc_res_n - Fc_n + Ec_sft * ec_n ) / Ec_sft; //strain at start of constant post-peak branch
	
	nn = 0.8 + Fc * 1000 / 2500;
	r_n = Fc/0.75 - 1.9;
	//xcr_n = 1.04;
		
	// Tensile Branch
	np = 0.8 + Fc * 1000/ 2500;  
	Fc_p = 6 * sqrt(Fc*1000)/1000; //tensile strength
	ec_p = 1.23 * Fc_p / Ec;
	//ecr_p = 5.3*ec_p;	
	ecr_p = 4.0*ec_p;	
	
		r_p = 4;	
	double n_p = fabs(Ec * ec_p / Fc_p);
	double xcr = ecr_p/ec_p;
	double ycr, zcr;
	tsaiEquation(xcr, r_p, n_p, ycr, zcr);
	ecrk_p = ecr_p - (Fc_p * ycr)/(Ec * zcr);
	
	#ifdef _HAJJAR_COMPOSITE_DEBUG
		// Print some data to a file
		ofstream initConcData;
		initConcData.open("initConcData.dat",ios::app);
		initConcData<<endl<<endl<<"Inital Concrete Data:"<<endl;
		initConcData<<"Ec: "<<Ec<<" Fc: "<<Fc<<endl;
		initConcData<<"ec_unconf: "<<ec_unconf<<" k2: "<<k2<<" x: "<<x<<endl;
		initConcData<<"Fc_n: "<<Fc_n<<" ec_n: "<<ec_n<<" Fc_res_n: "<<Fc_res_n<<" Ec_sft: "<<Ec_sft<<endl; 
		initConcData<<"Fc_p: "<<Fc_p<<" ec_p: "<<ec_p<<" ecr_p: "<<ecr_p<<" ecrk_p: "<<ecrk_p<<endl; 
	#endif
	
	// Initilize Parameters
	isCracked = false;	CisCracked = false;
	Trule = 0; Crule = 0;
	//isCracked = true;	CisCracked = true; // For Testing in Post-Cracking
	//Trule = 0; Crule = 1; // For Testing in Post-Cracking	
    Tstrain = 0.0; Cstrain = 0.0;
	Tstress = 0.0; Cstress = 0.0;
	Ttangent = 0.0; Ctangent = 0.0;
	eo = 0.0; Ceo = 0.0;
	er1 = 0.0; er3 = 0.0; er6 = 0.0; er13 = 0.0; er14 = 0.0;
	Cer1 = 0.0; Cer3 = 0.0; Cer6 = 0.0; Cer13 = 0.0; Cer14 = 0.0;
	eb = 0.0; Ceb = 0.0;
	Esec_n = 0.0; CEsec_n = 0.0;
	Epl_n = 0.0; CEpl_n = 0.0;
	delta_f_n = 0.0; Cdelta_f_n = 0.0;
	delta_e_n = 0.0; Cdelta_e_n = 0.0;
	epl_n = 0.0; Cepl_n = 0.0;
	fnew_n = 0.0; Cfnew_n = 0.0;
	Enew_n = 0.0; CEnew_n = 0.0;
	ere_n = 0.0; Cere_n = 0.0;
	fnew_str_n = 0.0; Cfnew_str_n = 0.0;
	Enew_str_n = 0.0; CEnew_str_n = 0.0;
	ere_str_n = 0.0; Cere_str_n = 0.0;
	f71target = 0.0; Cf71target = 0.0;
	E71target = 0.0; CE71target = 0.0;
	e72target = 0.0; Ce72target = 0.0;
	f3target = 0.0; Cf3target = 0.0;
	f81target = 0.0; Cf81target = 0.0;
	E81target = 0.0; CE81target = 0.0;
	e82target = 0.0; Ce82target = 0.0;
	f4target = 0.0; Cf4target = 0.0;
	
	delta_e_p = 0.0; Cdelta_e_p = 0.0;  
	delta_f_p = 0.0; Cdelta_f_p = 0.0;
	fnew_p = 0.0; Cfnew_p = 0.0;
	Enew_p = 0.0; CEnew_p = 0.0; 
	ere_p = 0.0; Cere_p = 0.0; 
	fnew_str_p = 0.0; Cfnew_str_p = 0.0; 
	Enew_str_p = 0.0; CEnew_str_p = 0.0;
	ere_str_p = 0.0; Cere_str_p = 0.0; 
	epl_p = 0.0; Cepl_p = 0.0; 
	Epl_p = 0.0; CEpl_p = 0.0;
	er2 = 0.0; er4 = 0.0; er9 = 0.0; er10 = 0.0; 
	Cer2 = 0.0; Cer4 = 0.0; Cer9 = 0.0; Cer10 = 0.0;
	Esec_p = 0.0; CEsec_p = 0.0; 
	ea = 0.0; Cea = 0.0;
}

CCFT_concMaterial::CCFT_concMaterial()
  :UniaxialMaterial(0,MAT_TAG_CCFT_conc)
  // @todo declare things as zero
{
	//DOES NOTHING
}

CCFT_concMaterial::~CCFT_concMaterial()
{
	//DOES NOTHING
}

int 
CCFT_concMaterial::setTrialStrain(double strain, double strainRate)
{
	double strain_incr;

	#ifdef _HAJJAR_COMPOSITE_DEBUG
		ofstream trialConcData;
		trialConcData.open("trialConcData.dat",ios::app);
	#endif

	// Return state variables to last commited values
	backToCommitStateVar();		
		
	// Define trial strain and strain increment
	//Tstrain = Cstrain + strain; // For use with Cenk's elements
	Tstrain = strain; // For use with other elements
	strain_incr = Tstrain - Cstrain; // Works for both
	
	//trialConcData<<endl<<" Fc: "<<Fc<<" Tstrain: "<<Tstrain<<" Cstrain: "<<Cstrain<<" strain_incr: "<<strain_incr<<endl;


	if(isCracked == false) { // Pre-Cracking 
		switch ( Crule ) {
		case 0:
			// It should only arrive here before loading has started
			if( fabs(Tstrain) < SMALL_STRAIN ) {
				Trule = 0;
			} else if( Tstrain <= 0.0 ) { 
				// Compressive Strain
				Trule = 1;
			} else if ( Tstrain <= ecrk_p + eo ) { 
				// Tensile Strain less than cracking 
				Trule = 2;
			} else { 
				// Tensile strain greater than cracking
				isCracked = true;
				Trule = 6;
			}
			break;
			
		case 1:
			if ( strain_incr <= SMALL_STRAIN ) { 
				// Continued compressive loading
				Trule = 1;
			} else { 
				// Unloading
				er1 = Cstrain;
				double fr1, scratch;
				negativeEnvelope(er1,fr1,scratch);
				Esec_n = Ec * ( fabs( fr1/(Ec*ec_n) ) + 0.57 ) / ( fabs( er1/ec_n ) + 0.57 );
				Epl_n = 0.1*Ec*exp(-2*er1/ec_n);
				delta_f_n = 0.09*fr1*sqrt(er1/ec_n);
				delta_e_n = er1 / (1.15 + 2.75*fabs(er1/ec_n));
				epl_n = er1 - fr1 / Esec_n;
				fnew_n = fr1 - delta_f_n;
				Enew_n = fnew_n / (er1 - epl_n );
				ere_n = er1 + delta_e_n;
				
				// Start of shifting the tensile branch
				double eo_old = eo;
				
				// if xu_p < xu_n force a new reverasal from tensile branch and reset eo
				if ( fabs((er2-eo)/ec_p) < fabs(er1/ec_n) ) {
					er2 = er1 * ec_p / ec_n;
					eo = 0.0;
					eo_old = 0.0;
					// compute reversal for new er2
					double fr2, scratch;
					positiveEnvelope(er2,fr2,scratch);
					Esec_p = Ec * ( fabs( fr2/(Ec*ec_p) ) + 0.67 ) / ( fabs( (er2-eo)/ec_p ) + 0.67 );
					Epl_p = Ec / ( pow((er2-eo)/ec_p,1.1) + 1.0 );
					delta_f_p = 0.15*fr2;
					delta_e_p = 0.22*(er2-eo);
					epl_p = er2 - fr2 / Esec_p;
					fnew_p = fr2 - delta_f_p;
					Enew_p = fnew_p / (er2 - epl_p );
					ere_p = er2 + delta_e_p;
				}
				
				// compute new eo
				double fr2;
				positiveEnvelope(er2,fr2,scratch);
				double Deo = 2*fr2/(Esec_p+Epl_p);
				eo = epl_n - (er2-eo-Deo);
				
				// shift all necessary strains (some may not be necessary)
				er2 += (eo-eo_old);
				er4 += (eo-eo_old);
				er6 += (eo-eo_old);
				er9 += (eo-eo_old);
				ea += (eo-eo_old);
				e82target += (eo-eo_old);
				ere_p += (eo-eo_old);
				ere_str_p += (eo-eo_old);
				epl_p += (eo-eo_old);
				
				// End of shifting the tensile branch
				
				
				if ( Tstrain <= epl_n ) {
					Trule = 3;
					f3target = Cstress;
				} else if ( Tstrain <= er2 ) {
					Trule = 9;
				} else if ( Tstrain <= ere_p ) {
					Trule = 8;
					f81target = fnew_p;
					E81target = Enew_p;
					e82target = ere_p; // Because arriving on Rule 8 after full unloading
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else { 
					// Tensile strain greater than cracking
					isCracked = true;
					Trule = 6;
				}
			}
			break;
			
		case 2:
			if ( strain_incr >= -SMALL_STRAIN ) { // Continued tensile loading
				if ( Tstrain <= ecrk_p + eo ) { 
					// Before cracking
					Trule = 2;
				} else { 
					// After Cracking
					isCracked = true;
					Trule = 6;
				}
			} else { 
				// Unloading
				er2 = Cstrain;
				double fr2, scratch;
				positiveEnvelope(er2,fr2,scratch);
				Esec_p = Ec * ( fabs( fr2/(Ec*ec_p) ) + 0.67 ) / ( fabs( (er2-eo)/ec_p ) + 0.67 );
				Epl_p = Ec / ( pow((er2-eo)/ec_p,1.1) + 1.0 );
				delta_f_p = 0.15*fr2;
				delta_e_p = 0.22*(er2-eo);
				epl_p = er2 - fr2 / Esec_p;
				fnew_p = fr2 - delta_f_p;
				Enew_p = fnew_p / (er2 - epl_p );
				ere_p = er2 + delta_e_p;
				if ( Tstrain >= epl_p ) {
					Trule = 4;
					f4target = Cstress;
				} else if ( Tstrain >= er1 ) {
					Trule = 10;
				} else if ( Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else {
					Trule = 1;
				}
			}
			break;
			
		case 3:
			if ( strain_incr >= -SMALL_STRAIN ) { // Continued tensile loading
				if ( Tstrain <= epl_n ) { 
					Trule = 3;
					// target remains as it is
				} else if ( Tstrain <= er2 ) {
					Trule = 9;
				} else if ( Tstrain <= ere_p ) {
					Trule = 8;
					f81target = fnew_p;
					E81target = Enew_p;
					e82target = ere_p; // Because arriving on Rule 8 after full unloading
				} else if ( Tstrain <= ecrk_p + eo ) {	
					Trule = 2;
				} else { 
					isCracked = true;
					Trule = 6;
				}
			} else { 
				// Reversal from Rule 3
				er3 = Cstrain;
				double fr1, scratch;
				negativeEnvelope( er1, fr1, scratch );
				double fr3 = Cstress;
				fnew_str_n = fr1 - delta_f_n * ( er1 - er3 )/( er1 - epl_n );
				Enew_str_n = ( fnew_str_n - fr3 ) / ( er1 - er3 );
				ere_str_n = er1 + delta_e_n * ( er1 - er3 ) / ( er1 - epl_n );
				
				if ( Tstrain >= er1 ) {
					Trule = 16;
					f71target = fnew_str_n;
					E71target = Enew_str_n; 
					e72target = ere_str_n; // Because arriving on Rule 16 after partial unloading			
				} else if ( Tstrain >= ere_str_n )  {
					Trule = 7;
					f71target = fnew_str_n;
					E71target = Enew_str_n;
					e72target = ere_str_n; // Because arriving on Rule 7 after partial unloading
				} else {
					Trule = 1;
				}
			}
				
			break;
			
		case 4:
			if ( strain_incr <= SMALL_STRAIN ) { 
				// Continued compressive unloading
				if ( Tstrain >= epl_p ) {
					Trule = 4;
					// target remains as it is
				} else if ( Tstrain >= er1 ) {
					Trule = 10;
				} else if ( Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else {
					Trule = 1;
				}
			} else { 
				// Reversal
				er4 = Cstrain;
				double fr2, scratch;
				positiveEnvelope( er2, fr2, scratch );
				double fr4 = Tstress;
				fnew_str_p = fr2 - delta_f_p * ( er2 - er4 )/( er2 - epl_p );
				Enew_str_p = ( fnew_str_p - fr4 ) / ( er2 - er4 );
				ere_str_p = er2 + delta_e_p * ( er2 - er4 ) / ( er2 - epl_p );
				
				if ( Tstrain <= er2 ) {
					Trule = 17;
					f81target = fnew_str_p;
					E81target = Enew_str_p; 
					e82target = ere_str_p; // Because arriving on Rule 17 after partial unloading
				} else if ( Tstrain <= ere_str_p ) {
					Trule = 8;
					f81target = fnew_str_p;
					E81target = Enew_str_p;
					e82target = ere_str_p; // Because arriving on Rule 8 after partial unloading
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else { 
					isCracked = true;
					Trule = 6;
				}
			}
			break;
			
		case 7:
			if ( strain_incr <= SMALL_STRAIN ) {
				// Continued compressive
				if ( Tstrain >= e72target ) {
					Trule = 7;
				} else {
					Trule = 1;
				}
			} else {
				// Reversal from Second Branch
				er1 = Cstrain;
				double fr1, scratch;
				negativeEnvelope( er1, fr1, scratch );
				Esec_n = Ec * ( fabs( fr1/(Ec*ec_n) ) + 0.57 ) / ( fabs( er1/ec_n ) + 0.57 );
				Epl_n = 0.1*Ec*exp(-2*er1/ec_n);
				delta_f_n = 0.09*fr1*sqrt(er1/ec_n);
				delta_e_n = er1 / (1.15 + 2.75*fabs(er1/ec_n));
				epl_n = er1 - fr1 / Esec_n;
				fnew_n = fr1 - delta_f_n;
				Enew_n = fnew_n / (er1 - epl_n );
				ere_n = er1 + delta_e_n;
				if ( Tstrain <= epl_n ) {
					Trule = 3;
					f3target = Cstress;
				} else if ( Tstrain <= er2 ) {
					Trule = 9;
				} else if ( Tstrain <= ere_p ) {
					Trule = 8;
					f81target = fnew_p;
					E81target = Enew_p;
					e82target = ere_p; // Because arriving on Rule 8 after full unloading
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else { 
					isCracked = true;
					Trule = 6;
				}
			}
			break;
			
		case 8:
			if ( strain_incr >= SMALL_STRAIN ) {
				// Continued tensile
				if ( Tstrain <= e82target ) {
					Trule = 8;
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else { 
					// Tensile strain greater than cracking
					isCracked = true;
					Trule = 6;
				}
			} else {
				// Reversal from Second Branch of Rule 8
				er2 = Cstrain;
				double fr2, scratch;
				positiveEnvelope(er2,fr2,scratch);
				Esec_p = Ec * ( fabs( fr2/(Ec*ec_p) ) + 0.67 ) / ( fabs( (er2-eo)/ec_p ) + 0.67 );
				Epl_p = Ec / ( pow((er2-eo)/ec_p,1.1) + 1.0 );
				delta_f_p = 0.15*fr2;
				delta_e_p = 0.22*er2;
				epl_p = er2 - fr2 / Esec_p;
				fnew_p = fr2 - delta_f_p;
				Enew_p = fnew_p / (er2 - epl_p );
				ere_p = er2 + delta_e_p;
				if ( Tstrain >= epl_p ) {
					Trule = 4;
					f4target = Cstress;
				} else if ( Tstrain >= er1 ) {
					Trule = 10;
				} else if ( Tstrain <= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else { 
					Trule = 1;
				}
			}
			break;
			
		case 9:
			if ( strain_incr >= -SMALL_STRAIN ) { 
				// Continued tensile loading
				if ( Tstrain <= er2 ) {
					Trule = 9;
				} else if ( Tstrain <= ere_p ) {
					Trule = 8;
					f81target = fnew_p;
					E81target = Enew_p;
					e82target = ere_p; // Because arriving on Rule 8 after full unloading
				} else if ( Tstrain <= ecrk_p + eo ) {	
					Trule = 2;
				} else { 
					isCracked = true;
					Trule = 6;
				}
			} else { 
				// Reversal from Rule 9
				er9 = Cstrain;
				eb = er1 - ((er9 - epl_n)/(er2 - epl_n))*(er1-epl_p);
				
				if ( Tstrain >= eb ) {
					Trule = 11;
				} else if ( Tstrain >= er1 ) {
					Trule = 10;
				} else if ( Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else {
					Trule = 1;
				}
			}
			break;
			
		case 10:
			if ( strain_incr <= SMALL_STRAIN ) { 
				// Continued compressive 
				if ( Tstrain >= er1 ) {
					Trule = 10;
				} else if ( Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else {
					Trule = 1;
				}
			} else { 
				// Reversal
				er10 = Cstrain;
				ea = epl_n + ((er1 - er10)/(er1 - epl_p))*(er2-epl_n);
				
				if ( Tstrain <= ea ) {
					Trule = 12;
				} else if ( Tstrain <= er2 ) {
					Trule = 9;
				} else if ( Tstrain < ere_p ) {
					Trule = 8;
					f81target = fnew_p;
					E81target = Enew_p;
					e82target = ere_p; // Because arriving on Rule 8 after full unloading
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else { 
					isCracked = true;
					Trule = 6;
				}
			}
			break;
			
		case 11:
			if ( strain_incr <= SMALL_STRAIN ) { 
				// Continued compressive
				if ( Tstrain >= eb ) {
					Trule = 11;
				} else if ( Tstrain >= er1 ) {
					Trule = 10;
				} else if ( Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else {
					Trule = 1;
				}
			} else { 
				// Reversal
				if ( Tstrain <= er9 ) {
					Trule = 11;
				} else if ( Tstrain <= er2 ) {
					Trule = 9;
				} else if ( Tstrain < ere_p ) {
					Trule = 8;
					f81target = fnew_p;
					E81target = Enew_p;
					e82target = ere_p; // Because arriving on Rule 8 after full unloading
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else { 
					isCracked = true;
					Trule = 6;
				}
			}
			break;
			
		case 12:
			if ( strain_incr >= -SMALL_STRAIN ) { 
				// Continued tensile loading
				if ( Tstrain < ea ) {
					Trule = 12;
				} else if ( Tstrain <= er2 ) {
					Trule = 9;
				} else if ( Tstrain <= ere_p ) {
					Trule = 8;
					f81target = fnew_p;
					E81target = Enew_p;
					e82target = ere_p; // Because arriving on Rule 8 after full unloading
				} else if ( Tstrain <= ecrk_p + eo ) {	
					Trule = 2;
				} else { 
					isCracked = true;
					Trule = 6;
				}
			} else { 
				// Reversal from Rule 12
				if ( Tstrain >= er10 ) {
					Trule = 12;
				} else if ( Tstrain >= er1 ) {
					Trule = 10;
				} else if ( Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else {
					Trule = 1;
				}
			}
			break;

		case 16:
			if ( strain_incr <= SMALL_STRAIN ) {
				if ( Tstrain > er1 ) {
					Trule = 16;
					// Targets stay as is
				} else if ( Tstrain > e72target ) {
					Trule = 7;
					f71target = fnew_str_n;
					E71target = Enew_str_n;
					e72target = ere_str_n; // Because arriving on Rule 7 after partial unloading
				} else {
					Trule = 1;
				}
			} else {
				// Reversal from First Branch of Rule 7, which is renamed Rule 16
				if ( Tstrain <= er3 ) {
					Trule = 16;
					// Targets stay as is
				} else if ( Tstrain <= epl_n ) { 
					Trule = 3;
					// target remains as it is
				} else if ( Tstrain <= er2 ) {
					Trule = 9;
				} else if ( Tstrain <= ere_p ) {
					Trule = 8;
					f81target = fnew_p;
					E81target = Enew_p;
					e82target = ere_p; // Because arriving on Rule 8 after full unloading
				} else if ( Tstrain <= ecrk_p + eo ) {	
					Trule = 2;
				} else { 
					isCracked = true;
					Trule = 6;
				}
			}
			break;					
			
		case 17:
			if ( strain_incr >= SMALL_STRAIN ) {
				// Continued tensile
				if ( Tstrain <= er2 ) {
					Trule = 17;
					// Targets stay as is
				} else if ( Tstrain <= e82target ) {
					Trule = 8;
					f81target = fnew_str_p;
					E81target = Enew_str_p;
					e82target = ere_str_p; // Because arriving on Rule 8 after partial unloading
				} else if ( Tstrain <= ecrk_p + eo ) {
					Trule = 2;
				} else { 
					// Tensile strain greater than cracking
					isCracked = true;
					Trule = 6;
				}
			} else {
				// Reversal from First Branch
				if ( Tstrain >= er4 ) {
					Trule = 17;
					// Targets stay as is
				} else if ( Tstrain >= epl_p ) {
					Trule = 4;
					// target remains as it is
				} else if ( Tstrain >= er1 ) {
					Trule = 10;
				} else if ( Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else { 
					Trule = 1;
				}
			}
			break;
			
		default:
			// It should not arrive here
			break;
			
	} // End of switch for Crule (for pre-cracking region)

	} else { // Post-Cracking
		//trialConcData<<"I got to the post cracking area"<<endl;

		switch ( Crule ) {
		case 1:
			//trialConcData<<"I got to the rule 1 area"<<endl;
			if ( strain_incr <= SMALL_STRAIN ) { 
				// Continued compressive loading
				Trule = 1;
			} else { // Unloading
				er1 = Cstrain;
				double fr1, scratch;
				negativeEnvelope(er1,fr1,scratch);
				Esec_n = Ec * ( fabs( fr1/(Ec*ec_n) ) + 0.57 ) / ( fabs( er1/ec_n ) + 0.57 );
				Epl_n = 0.1*Ec*exp(-2*er1/ec_n);
				delta_f_n = 0.09*fr1*sqrt(er1/ec_n);
				delta_e_n = er1 / (1.15 + 2.75*fabs(er1/ec_n));
				epl_n = er1 - fr1 / Esec_n;
				fnew_n = fr1 - delta_f_n;
				Enew_n = fnew_n / (er1 - epl_n );
				ere_n = er1 + delta_e_n;
				//trialConcData<<"I got to the unloading from rule 1 area"<<endl;
				//trialConcData<<" er1: "<<er1<<" fr1: "<<fr1<<" Esec_n: "<<Esec_n<<" Epl_n: "<<Epl_n<<" delta_f_n: "<<delta_f_n<<endl;
				//trialConcData<<" delta_e_n: "<<delta_e_n<<" epl_n: "<<epl_n<<" fnew_n: "<<fnew_n<<" Enew_n: "<<Enew_n<<" ere_n: "<<ere_n<<endl;
				if ( Tstrain <= epl_n ) {
					Trule = 3;
					f3target = Cstress;
				} else {
					Trule = 6;
				}
			}
			break;
			
		case 3:
			if ( strain_incr >= -SMALL_STRAIN ) { // Continued unloading
				if ( Tstrain <= epl_n ) {
					Trule = 3;
				} else {
					Trule = 6;
				}
			} else { // Reloading (reversal)
				er3 = Cstrain;
				double fr1, scratch;
				negativeEnvelope( er1, fr1, scratch );
				double fr3 = Cstress;
				fnew_str_n = fr1 - delta_f_n * ( er1 - er3 )/( er1 - epl_n );
				Enew_str_n = ( fnew_str_n - fr3 ) / ( er1 - er3 );
				ere_str_n = er1 + delta_e_n * ( er1 - er3 ) / ( er1 - epl_n );
				
				if ( Tstrain >= er1 ) {
					Trule = 16;
					f71target = fnew_str_n;
					E71target = Enew_str_n;	
					e72target = ere_str_n; // Because arriving on Rule 7 after partial unloading				
				} else if ( Tstrain >= ere_str_n ) {
					Trule = 7;
					f71target = fnew_str_n;
					E71target = Enew_str_n;
					e72target = ere_str_n; // Because arriving on Rule 7 after partial unloading
				} else {
					Trule = 1;
				}
			}
			break;
			
		case 6:
			if ( strain_incr >= -SMALL_STRAIN ) { // Continued tensile loading
				Trule = 6;
			} else { // reversal from 6
				er6 = Cstrain;
				if ( Tstrain >= er1 ) {
					Trule = 13;
				} else if (Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else {
					Trule = 1;
					negativeEnvelope(Tstrain,Tstress,Ttangent);
				}
			}
			break;
			
		case 7:
			if ( strain_incr <= SMALL_STRAIN ) {
				if ( Tstrain > e72target ) {
					Trule = 7;
				} else {
					Trule = 1;
				}
			} else {
				// Reversal from Second Branch
				er1 = Cstrain;
				double fr1, scratch;
				negativeEnvelope( er1, fr1, scratch );
				Esec_n = Ec * ( fabs( fr1/(Ec*ec_n) ) + 0.57 ) / ( fabs( er1/ec_n ) + 0.57 );
				Epl_n = 0.1*Ec*exp(-2*er1/ec_n);
				delta_f_n = 0.09*fr1*sqrt(er1/ec_n);
				delta_e_n = er1 / (1.15 + 2.75*fabs(er1/ec_n));
				epl_n = er1 - fr1 / Esec_n;
				fnew_n = fr1 - delta_f_n;
				Enew_n = fnew_n / (er1 - epl_n );
				ere_n = er1 + delta_e_n;
				if ( Tstrain <= epl_n ) {
					Trule = 3;
					f3target = Cstress;
				} else {
					Trule = 6;
				}
			}
			break;
			
		case 13:
			if ( strain_incr <= SMALL_STRAIN ) {
				if ( Tstrain >= er1 ) {
					Trule = 13;
				} else if ( Tstrain >= ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
					//trialConcData << "I just got to Rule 7 from Rule 13, f71target" << f71target << "E71target" << E71target<< "e72target" << e72target << endl;
				} else {
					Trule = 1;
				}
			} else {
				er13 = Cstrain;
				eb = er13 - Cstress/Esec_n;
				if ( Tstrain < eb ) {
					Trule = 14;
				} else {
					Trule = 6;
				}
			}
			break;
			
		case 14:
			if ( strain_incr >= -SMALL_STRAIN ) {
				if ( Tstrain < eb ) {
					Trule = 14;
				} else {
					Trule = 6;
				}
			} else { // Reversal from Rule 14
				er14 = Cstrain;
				if ( Tstrain >= er13 ) {
					Trule = 15;
				} else if ( Tstrain > er1 ) {
					Trule = 13;
				} else if ( Tstrain > ere_n ) {
					Trule = 7;
					// will be the second branch
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else {
					Trule = 1;
				}
			}
			break;
			
		case 15:
			if ( strain_incr <= SMALL_STRAIN ) {
				if ( Tstrain > er13) {
					Trule = 15;
				} else if ( Tstrain > er1 ) {
					Trule = 13;
				} else if ( Tstrain > ere_n ) {
					Trule = 7;
					f71target = fnew_n;
					E71target = Enew_n;
					e72target = ere_n; // Because arriving on Rule 7 after full unloading
				} else {
					Trule = 1;
				}
			} else {
				if ( Tstrain <= er14 ) {
					Trule = 15; //stay in Rule 15
				} else if ( Tstrain <= eb ) {
					Trule = 14;
				} else {
					Trule = 6;
				}
			}
			break;
		
		case 16:
			if ( strain_incr <= SMALL_STRAIN ) {
				if ( Tstrain > er1 ) {
					Trule = 16;
					// Targets stay as is
				} else if ( Tstrain > e72target ) {
					Trule = 7;
					f71target = fnew_str_n;
					E71target = Enew_str_n;
					e72target = ere_str_n; // Because arriving on Rule 7 after partial unloading
				} else {
					Trule = 1;
				}
			} else {
				// Reversal from First Branch of Rule 7, which is renamed Rule 16
				if ( Tstrain <= er3 ) {
					Trule = 16;
					// Targets stay as is
				} else if ( Tstrain <= epl_n ) {
					Trule = 3;
				} else {
					Trule = 6;
				}
			}
			break;			
			
		default:
			// It should not arrive here
			break;
			
	} // End of switch for Crule
	} // End of if statment for isCracked
	
	// Now that we have found the new rule and updated any state variables, we find the new stress and tangent
	switch ( Trule ) {
		case 0:
			// Loading hasn't started yet
			Tstress = SMALL_STRESS;
			Ttangent = Ec;
			break;
		case 1:
			negativeEnvelope(Tstrain,Tstress,Ttangent);
			break;
		case 2:
			positiveEnvelope(Tstrain,Tstress,Ttangent);
			break;
		case 3:
			transitionCurve(Tstrain, Tstress, Ttangent, er1, f3target, Ec, epl_n, 0.0, Epl_n);
			//trialConcData<<"I am in the rule 3 area"<<endl;
			//trialConcData<<" Tstrain: "<<Tstrain<<" er1: "<<er1<<" fi: "<<fi<<" Ec: "<<Ec<<" epl_n: "<<epl_n<<" Epl_n: "<<Epl_n<<endl;
			break;
		case 4:
			transitionCurve(Tstrain, Tstress, Ttangent, er2, f4target, Ec, epl_p, 0.0, Epl_p);
			break;
		case 6:
			Tstress = SMALL_STRESS;
			Ttangent = SMALL_TANGENT;
			break;
		case 7:
			{
				// Second Branch of Rule 7 remains Rule 7
				// First Branch has been renamed Rule 16
				double ff, Ef;
				negativeEnvelope(e72target,ff,Ef); // Calculate target stress and tangent		
				transitionCurve(Tstrain, Tstress, Ttangent, er1, f71target, E71target, e72target, ff, Ef);
			}
			break;
		case 8:
			{
				// Second Branch of RUle 8 remains Rule 8
				// First Branch has been renamed Rule 17
				double ff, Ef;
				positiveEnvelope(e82target,ff,Ef); // Calculate target stress and tangent		
				//trialConcData<<"    I am in rule 8"<<endl;
				//trialConcData<<"Rule 8 Set Area, er2: "<<er2<<" f81target: "<<f81target<<" E81target: "<<E81target<<" e82target: "<<e82target<<" ff: "<<ff<<" Ef: "<<Ef<<endl;
				transitionCurve(Tstrain, Tstress, Ttangent, er2, f81target, E81target, e82target, ff, Ef);
			}
			break;
		case 9:
			transitionCurve(Tstrain, Tstress, Ttangent, epl_n, 0.0, Epl_n, er2, fnew_p, Enew_p);
			break;
		case 10:
			transitionCurve(Tstrain, Tstress, Ttangent, epl_p, 0.0, Epl_p, er1, fnew_n, Enew_n);
			break;
		case 11:
			{
				double fr9, fb, Eb, scratch;
				//Find the inital stress (using Rule 9)
				transitionCurve(er9, fr9, scratch, epl_n, 0.0, Epl_n, er2, fnew_p, Enew_p);
				//Find the final stress and tangent (using Rule 10)
				transitionCurve(eb, fb, Eb, epl_p, 0.0, Epl_p, er1, fnew_n, Enew_n);
				//Set Tstress and Ttangent (Using Rule 11)
				transitionCurve(Tstrain, Tstress, Ttangent, er9, fr9, Ec, eb, fb, Eb);
			}
			break;
		case 12:
			{
				double fr10, fa, Ea, scratch;
				//Find the inital stress (using Rule 10)
				transitionCurve(er10, fr10, scratch, epl_p, 0.0, Epl_p, er1, fnew_n, Enew_n);
				//trialConcData<<"Rule 12 Set Area (Rule 10), epl_p: "<<epl_p<<" Epl_p: "<<Epl_p<<" er1: "<<er1<<" fnew_n: "<<fnew_n<<" Enew_n: "<<Enew_n<<endl;
				//Find the final stress and tangent (using Rule 9)
				transitionCurve(ea, fa, Ea, epl_n, 0.0, Epl_n, er2, fnew_p, Enew_p);
				//trialConcData<<"Rule 12 Set Area (Rule 9), epl_n: "<<epl_n<<" Epl_n: "<<Epl_n<<" er2: "<<er2<<" fnew_p: "<<fnew_p<<" Enew_p: "<<Enew_p<<endl;
				//Set Tstress and Ttangent (Using Rule 12)
				transitionCurve(Tstrain, Tstress, Ttangent, er10, fr10, Ec, ea, fa, Ea);
				//trialConcData<<"Rule 12 Set Area, er10: "<<er10<<" fr10: "<<fr10<<" Ec: "<<Ec<<" ea: "<<ea<<" fa: "<<fa<<" Ea: "<<Ea<<endl;
			}
			break;
		case 13:
			transitionCurve(Tstrain, Tstress, Ttangent, er6, 0.0, 0.0, er1, fnew_n, Enew_n);
			break;
		case 14:
			{
				double fr13, scratch;
				//Find the inital stress (Using Rule 13)
				transitionCurve(er13, fr13, scratch, er6, 0.0, 0.0, er1, fnew_n, Enew_n);
				//Set Tstress and Ttangent (Using Rule 14)
				transitionCurve(Tstrain, Tstress, Ttangent, er13, fr13, Ec, eb, 0.0, 0.0);
			}
			break;
		case 15:
			{
				double fr13, fr14, Er13, scratch;
				//Find fr13 and Er13 (Using Rule 13)
				transitionCurve(er13, fr13, Er13, er6, 0.0, 0.0, er1, fnew_n, Enew_n);
				//Find fr14 (Using Rule 14)
				transitionCurve(er14, fr14, scratch, er13, fr13, Ec, eb, 0.0, 0.0);
				//Set Tstress and Ttangent (Using Rule 15)
				transitionCurve(Tstrain, Tstress, Ttangent, er14, fr14, Ec, er13, fr13, Er13);
			}
			break;
		case 16:
			{
				// First Branch of Rule 7 (per Chang and Mander) is renamed Rule 16
				double fr3, scratch;
				transitionCurve(er3, fr3, scratch, er1, f3target, Ec, epl_n, 0.0, Epl_n);
				//trialConcData<<"Rule 16 Set Area (Rule 3), er1: "<<er1<<" f3target: "<<f3target<<" Ec: "<<Ec<<" epl_n: "<<epl_n<<" Epl_n: "<<Epl_n<<endl;
					// starting stress for Rule 16 (first branch) comes from Rule 3
				//trialConcData<<"Rule 16 Set Area, er3: "<<er3<<" fr3: "<<fr3<<" Ec: "<<Ec<<" er1: "<<er1<<" f71target: "<<f71target<<" E71target: "<<E71target<<endl;
				transitionCurve(Tstrain, Tstress, Ttangent, er3, fr3, Ec, er1, f71target, E71target);
			}
			break;
		case 17:
			{
				// First Branch of Rule 8 (per Chang and Mander) is renamed Rule 17
				double fr4, scratch;
				transitionCurve(er4, fr4, scratch, er2, f4target, Ec, epl_p, 0.0, Epl_p);
					// starting stress for Rule 8 (first branch) comes from Rule 4
				transitionCurve(Tstrain, Tstress, Ttangent, er4, fr4, Ec, er2, f81target, E81target);
			}
			break;
		default:
			// It should not arrive here
			break;			
	}
	
	//trialConcData<<Tstrain<<"   "<<Tstress<<"    "<<Ttangent<<"    "<<Trule<<"   "<<Cstrain<<"  "<<Cstress<<"  "<<Ctangent<<"  "<<Crule<<endl;
	return 0;
}

double 
CCFT_concMaterial::getStress(void)
{
  return Tstress;
}

double 
CCFT_concMaterial::getTangent(void)
{
  return Ttangent;
}

double 
CCFT_concMaterial::getInitialTangent(void)
{
  return Ec;
}

double 
CCFT_concMaterial::getStrain(void)
{
  return Tstrain;
}

int 
CCFT_concMaterial::commitState(void)
{
	#ifdef _HAJJAR_COMPOSITE_DEBUG
		ofstream commitConcData;
		commitConcData.open("commitConcData.dat",ios::app);  
		commitConcData<<Tstrain<<"   "<<Tstress<<"    "<<Ttangent<<"    "<<Trule<<"   "<<Cstrain<<"  "<<Cstress<<"  "<<Ctangent<<"  "<<Crule<<endl;
	#endif

	this->commitStateVar();
	return 0;
}

int 
CCFT_concMaterial::revertToLastCommit(void)
{
	this->backToCommitStateVar();
	return 0;
}

int 
CCFT_concMaterial::revertToStart(void)
{
	// @todo code CCFT_concMaterial::revertToStart(void)
	return Trule;
}


void 
CCFT_concMaterial::negativeEnvelope(double strain, double &stress, double &tangent)
{
	#ifdef _HAJJAR_COMPOSITE_DEBUG
		//ofstream trialConcData;
		//trialConcData.open("trialConcData.dat",ios::app);
	#endif

	if( strain <= ec_res_n ) {
		// constant post-peak branch
		stress = Fc_res_n;
		tangent = SMALL_TANGENT;
		// Cenk had a small slope here: strs = Fc_res_n + (0.000001) * ( strn - ec_res_n );
	} else if ( strain <= ec_n ) {
		// linear post-peak softening branch
		stress = Fc_n + Ec_sft * ( strain - ec_n );
		tangent = Ec_sft;
	} else if ( strain <= 0.0 ) {
		// nonlinear ascending curve defined by Chang and Mander
		double x_n, n_n, y, z;
		x_n = strain/ec_n;
		n_n = fabs(Ec * ec_n / Fc_n);
		tsaiEquation(x_n, r_n, n_n, y, z);
		stress = Fc_n * y;
		tangent = Ec * z;
		//trialConcData<<"I am in the nonlinear section"<<endl;
		//trialConcData<<"Fc: "<<Fc<<" y: "<<y<<" x_n: "<<x_n<<" n_n: "<<n_n<<" strs: "<<strs<<" strn: "<<strn<<" y: "<<y<<endl;
	} else {
		// it should not get here but if it does return zero
		stress = 0.0;
		tangent = 0.0;
	}
	
//	if ( strain <= -4.0/4000.0 ) {
//		stress = -4.0;
//		tangent = SMALL_TANGENT;
//	} else if ( strain <= 0.0 ) {
//		stress = strain*4000.0;
//		tangent = 4000.0;
//	} else {
//		stress = 0.0;
//		tangent = SMALL_TANGENT;
//	}
	
	return;
}

void 
CCFT_concMaterial::positiveEnvelope(double strain, double &stress, double &tangent)
{
	double n_p, x_p;
	n_p = fabs(Ec * ec_p / Fc_p);
	x_p = (strain - eo)/ec_p;
	 
	if( strain - eo <= 0 ) {
		// it should not get here but if it does return negative small stress
		stress = -1*SMALL_STRESS;
		tangent = 0.0;
		  
	} else if ( strain - eo <= ecr_p ) {
		// nonlinear pre-critical branch
		double y, z;
		tsaiEquation(x_p, r_p, n_p, y, z);
		stress = Fc_p * y;
		tangent = Ec * z;
		  
	} else if ( strain - eo <= ecrk_p ) {
		// linear softening branch between critical and cracked
		double x_cr, y_cr, z_cr;
		x_cr = ecr_p/ec_p;
		tsaiEquation(x_cr, r_p, n_p, y_cr, z_cr);
		stress = Fc_p * (y_cr + n_p * z_cr * (x_p - x_cr));
		tangent = Ec * z_cr;
		  
	} else {
		// cracked, return positive small stress
		stress = SMALL_STRESS;
		tangent = SMALL_TANGENT;
		//Cenk had this line: crack = 1;
		  
	}
	return;	
}


void 
CCFT_concMaterial::transitionCurve(double Tstrain, double &Tstress, double &Ttangent, double ei, double fi, double Ei, double ef, double ff, double Ef)
{
	#ifdef _HAJJAR_COMPOSITE_DEBUG
		ofstream trialConcData;
		trialConcData.open("trialConcData.dat",ios::app);
	#endif	
	//trialConcData<<"Start of Transition Curv  ff: "<<ff<<" fi: "<<fi<<" ef: "<<ef<<" ei: "<<ei<<endl;

	double R, A, Esec;
	if ( fabs (ei/ef - 1) < 1.0e-4) {
		// Inital and final strain are too close together, linear function
		Ttangent = (ff - fi) / (ef - ei);
		Tstress = fi + Ttangent * (Tstrain - ei);
		//trialConcData<<"transistionCurve -> first linear option"<<endl;
		return;
	} 
	
	Esec = ( ff - fi ) / ( ef - ei );
	//trialConcData<<"  ff: "<<ff<<" fi: "<<fi<<" ef: "<<ef<<" ei: "<<ei<<" Esec: "<<Esec<<endl;
	
	if( fabs(Esec-Ei) <= SMALL_NUMBER ) {
		R = 0.0;
	} else {
		R = ( Ef - Esec ) / ( Esec - Ei );
	}
	
	if(R < 0.0 || pow( fabs(ef - ei), R ) < SMALL_NUMBER ) {
		// Smooth transition curve without change of curvature not possible, linear function
		Ttangent = (ff - fi) / (ef - ei);
		Tstress = fi + Ttangent * (Tstrain - ei);
		//trialConcData<<"transistionCurve -> second linear option"<<endl;
		return;
	}
	
	if ( R > 100 ) {
		//trialConcData<<"Maybe the power fcn is bad, R too large, R: "<<R<<endl;
		Ttangent = (ff - fi) / (ef - ei);
		Tstress = fi + Ttangent * (Tstrain - ei);
		//trialConcData<<"transistionCurve -> third linear option"<<endl;
		return;
	}
	
	if ( pow( fabs(ef - ei), R ) < SMALL_NUMBER ) {
		//trialConcData<<"Maybe the power fcn is bad, pow( fabs(ef - ei), R ) too small: "<<pow( fabs(ef - ei), R )<<endl;
		Ttangent = (ff - fi) / (ef - ei);
		Tstress = fi + Ttangent * (Tstrain - ei);
		//trialConcData<<"transistionCurve -> fourth linear option"<<endl;
		return;
	}
	
	if ( fabs( Ef/Esec - 1 ) < 0.001 || fabs( Ei/Esec - 1) < 0.001 ) {
		//trialConcData<<"Maybe the power fcn is bad, Ef or Ei is too small close to Esec, Ei: "<<Ei<<" Ef: "<<Ef<<" Esec: "<<Esec<<endl;
		Ttangent = (ff - fi) / (ef - ei);
		Tstress = fi + Ttangent * (Tstrain - ei);
		//trialConcData<<"transistionCurve -> fifth linear option"<<endl;
		return;
	}
//		// Power function not ideal, multilinear function
//		// Intersection point of tangent lines from inital and final points
//		double ej = (Ei*ei - Ef*ef - fi + ff)/(Ei-Ef);
//		// Points to make the turns at
//		double eti = ei + (2.0/3.0)*(ej-ei);
//		double fti = fi + Ei * (eti - ei); 
//		double etf = ef + (2.0/3.0)*(ej-ef);
//		double ftf = ff + Ef * (etf - ef);
//		//trialConcData<<"transistionCurve -> multilinear option"<<endl;
//		if ( Tstrain <= max(ei,eti) && Tstrain >= min(ei,eti) ) {
//			// We are between ei and eti, now linearly interpolate
//			Ttangent = (fti - fi) / (eti - ei);
//			Tstress = fi + Ttangent * (Tstrain - ei);
//			return;
//		} else if ( Tstrain <= max(eti,etf) && Tstrain >= min(eti,etf) ) {
//			// We are between eti and etf, now linearly interpolate
//			Ttangent = (ftf - fti) / (etf - eti);
//			Tstress = fti + Ttangent * (Tstrain - eti);
//			return;		
//		} else if ( Tstrain <= max(etf,ef) && Tstrain >= min(etf,ef) ) {
//			// We are between etf and ef, now linearly interpolate
//			Ttangent = (ff - ftf) / (ef - etf);
//			Tstress = ftf + Ttangent * (Tstrain - etf);
//			return;		
//		} else {
//			// We should not be here
//			Tstress = 0.0;
//			Ttangent = 0.0;
//			//trialConcData<<"transistionCurve -> should not be here of the multi linear option"<<endl;
//			return;
//		}
//	} else {
		// Power function should be okay, use it
		//trialConcData<<"transistionCurve -> power option"<<endl;
		A = ( Esec - Ei ) / pow( fabs(ef - ei), R );
		Tstress = fi + (Tstrain-ei) * ( Ei + A * pow( fabs(Tstrain-ei), R ) );
		Ttangent = Ei + A * (R+1) * pow( fabs(Tstrain-ei), R );
		//trialConcData<<"  A: "<<A<<" R: "<<R<<" Esec: "<<Esec<<endl;
		return;
//	}  
}


void 
CCFT_concMaterial::tsaiEquation(double x, double r, double n, double &y, double &z)
{
	double D;
	if( r == 1 ) {
	    D = 1 + ( n - 1 + log( x ) ) * x;
	} else {
	    D = 1 + ( n - r / ( r - 1 ) ) * x + pow ( x , r ) / ( r - 1 );
	}
	y = ( n * x / D );
	z = ( 1 - pow( x, r ) ) / ( D * D );
	return;	
}


UniaxialMaterial *
CCFT_concMaterial::getCopy(void)
{
	CCFT_concMaterial *theCopy = new CCFT_concMaterial(this->getTag(), Fc, D, t, Fy, Es);
	theCopy->Ttangent = Ttangent;
	theCopy->Tstrain = Tstrain;
	theCopy->Tstress = Tstress;
	theCopy->Trule = Trule;
	// @todo code CCFT_concMaterial::getCopy(void)
	return theCopy;
}

int 
CCFT_concMaterial::sendSelf(int cTag, Channel &theChannel)
{
	// @todo code CCFT_concMaterial::sendSelf(int cTag, Channel &theChannel)
	return 0;
}

int 
CCFT_concMaterial::recvSelf(int cTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
	// @todo code CCFT_concMaterial::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
	return 0;	
}

void 
CCFT_concMaterial::Print(OPS_Stream &s, int flag)
{
  s<<"CCFT_conc, tag:"<<this->getTag()<<endln;
  s<<" Fcn:"<<Fc_n<<endln;
  s<<" Fcp:"<<Fc_p<<endln;
  s<<" Ec: "<<Ec<<endln;
  s<<" eo:"<<eo<<endln;
  s<<" Ec_sft: "<<Ec_sft<<endln;
  // @todo put better stuff in to CCFT_concMaterial::Print(OPS_Stream &s, int flag)
  return;
}

void
CCFT_concMaterial::backToCommitStateVar(void){
	isCracked = CisCracked;
	Trule = Crule;	
	Tstrain = Cstrain;
	Tstress = Cstress;
	Ttangent = Ctangent;
	eo = Ceo;
	er1 = Cer1; er3 = Cer3; er6 = Cer6; er13 = Cer13; er14 = Cer14;
	eb = Ceb;
	Esec_n = CEsec_n;
	Epl_n = CEpl_n;
	delta_f_n = Cdelta_f_n;
	delta_e_n = Cdelta_e_n;
	epl_n = Cepl_n;
	fnew_n = Cfnew_n;
	Enew_n = CEnew_n;
	ere_n = Cere_n;
	fnew_str_n = Cfnew_str_n;
	Enew_str_n = CEnew_str_n;
	ere_str_n = Cere_str_n;
	f71target = Cf71target;
	E71target = CE71target;
	e72target = Ce72target;
	f3target = Cf3target;
	f81target = Cf81target;
	E81target = CE81target;
	e82target = Ce82target;
	f4target = Cf4target;
	delta_e_p = Cdelta_e_p;  
	delta_f_p = Cdelta_f_p;
	fnew_p = Cfnew_p;
	Enew_p = CEnew_p; 
	ere_p = Cere_p; 
	fnew_str_p = Cfnew_str_p; 
	Enew_str_p = CEnew_str_p;
	ere_str_p = Cere_str_p; 
	epl_p = Cepl_p; 
	Epl_p = CEpl_p;
	er2 = Cer2;
	er4 = Cer4; 
	er9 = Cer9; 
	Esec_p = CEsec_p; 
	er10 = Cer10;
	ea = Cea;
}

void
CCFT_concMaterial::commitStateVar(void){
	CisCracked = isCracked;
	Crule = Trule;	
	Cstrain = Tstrain;
	Cstress = Tstress;
	Ctangent = Ttangent;
	Ceo = eo;
	Cer1 = er1; Cer3 = er3; Cer6 = er6; Cer13 = er13; Cer14 = er14;
	Ceb = eb;
	CEsec_n = Esec_n;
	CEpl_n = Epl_n;
	Cdelta_f_n = delta_f_n;
	Cdelta_e_n = delta_e_n;
	Cepl_n = epl_n;
	Cfnew_n = fnew_n;
	CEnew_n = Enew_n;
	Cere_n = ere_n;
	Cfnew_str_n = fnew_str_n;
	CEnew_str_n = Enew_str_n;
	Cere_str_n = ere_str_n;
	Cf71target = f71target;
	CE71target = E71target;
	Ce72target = e72target;  
	Cf3target = f3target;
	Cf81target = f81target;
	CE81target = E81target;
	Ce82target = e82target;
	Cf4target = f4target;
	Cdelta_e_p = delta_e_p;  
	Cdelta_f_p = delta_f_p;
	Cfnew_p = fnew_p;
	CEnew_p = Enew_p; 
	Cere_p = ere_p; 
	Cfnew_str_p = fnew_str_p; 
	CEnew_str_p = Enew_str_p;
	Cere_str_p = ere_str_p; 
	Cepl_p = epl_p; 
	CEpl_p = Epl_p;
	Cer2 = er2;
	Cer4 = er4; 
	Cer9 = er9; 
	CEsec_p = Esec_p; 
	Cer10 = er10;
	Cea = ea;
}

