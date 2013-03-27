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

// $Revision: 1.2 $
// $Date: 2008-07-03 18:08:38 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/uniaxial/RCFT_stlMaterial.cpp,v $

// Written: MHS
// Created: Aug 2001

// Modified by : Cenk Tort - tort0008@umn.edu
// University of Minnesota - Civil Engineering Department
// Date: Wed Jul 23 17:40:25 EDT 2003
// Description: This file contains the class implementation for
// cyclic uniaxial stress-strain relationship of steel for
// RCFT members.

#include "RCFT_stlMaterial.h"
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <fstream>

//using std::ofstream;
//using std::ios;
//using std::endl;
using namespace std;


RCFT_stlMaterial::RCFT_stlMaterial( int tag, double Fy, double Fu, double es, double depth, double thickness, double Eplo )
  :UniaxialMaterial(tag, MAT_TAG_RCFT_stl), fy(Fy), fu(Fu), Ee(es), D(depth), t(thickness), epo(Eplo) 
{
       if(fy > 52){
         Rbso =  1.06 * fy;
         //Rbso = 1.15 * fy;
	 //Epoi = 0.490 * 1.02 * 0.01 * Ee;
         Epoi = 0.0001 * Ee;
         alfa = 0.175;
	 a = -0.553;
	 bb = 6.47;
	 c = 34.8;
	 w = 2.67 / fy;
	 ksi = 8.04 * 0.001;
	 e = 700;
	 fE = 0.361;
       }else if( (fy > 40) && (fy <= 52) ){
         Rbso =  1.13 * fy;
         //Epoi = 0.411 * 3.41 * 0.01 * Ee;
         Epoi = 0.411 * 3.41 * 0.01 * Ee;
         //Epoi = 0.001 * Ee;
         alfa = 0.217;
         a = -0.0528;
         bb = 1.88;
         c = 18.7;
         w = 4.0 / fy;
         ksi = 1.52 * 0.001;
         //ksi = 10 * 0.001;
         //ksi = 20 * 0.001;
         //ksi = 5 * 0.001;
         e = 316;
         fE = 0.484;
       }else if(fy <= 40){
         Rbso =  1.15 * fy;
	 Epoi = 0.361 * 2.49 * 0.01 * Ee;
	 alfa = 0.191;
	 a = -0.505;
	 bb = 2.17;
	 c = 14.4;
	 w = 3.08 / fy;
	 ksi = 9.89 * 0.0001;
         e = 500;
         fE = 0.300;
       }

       if ( (D/t) * (fy/Ee) <= 0.08 ){
         Ksft = 0.001;
         //Ksft = - Epoi;
       }
       else if ( (D/t) * (fy/Ee) > 0.08 ){
         Ksft = ( 1033.50 * (D/t) * (fy/Ee) - 86.32 ) * fy;
       }

       if ( (D/t) * sqrt(fy/Ee) <= 0.92 ){
         lbW = 0.005;
       }
       else if ( ( (D/t) * sqrt(fy/Ee) > 0.92 ) && ( (D/t) * sqrt(fy/Ee) < 1.45 ) ){
         lbW = 0.0284 * ( (D/t) * sqrt(fy/Ee) ) - 0.0212;
       }
       else if ( (D/t) * sqrt(fy/Ee) >= 1.45 ){
         lbW = 0.02;   
       }

       //Ksft = 2846.34;
       //eplbf = 0.0032 - 0.0027;
       eplbf = ( 3.14 * ( pow ( (D/t) * sqrt(fy/Ee), -1.48 ) ) * (fy/Ee + 0.002) - (fy/Ee + 0.002) );
       //eplbf = 10;
       if ( Ksft < 0.01 ){
         fres_n = fy * 0.9999;
       }
       else{
         fres_n = (-7.31 * (D/t) * (fy/Ee) + 1.58) * fy;
       }
       //fres_n = 57.0;
       eresp_n = -( ( fres_n - fy ) / Ksft ) + 3.14 * ( pow ( (D/t) * sqrt(fy/Ee), -1.48 ) ) * (fy/Ee + 0.002) - (fy/Ee + 0.002);
       //eresp_n = 0.011 - 0.0027;

       //Rlso = CRlso = fy;
       //Rls = CRls = 0.0;
       Rlso = CRlso = fy * ( alfa - a * exp( -bb * epo * 100 ) - ( alfa - a - 1 ) * exp( -c * epo * 100 ) );
       Rls = CRls = fy * ( alfa - a * exp( -bb * epo * 100 ) - ( alfa - a - 1 ) * exp( -c * epo * 100 ) );
       Tls_p = Cls_p =  fy * ( alfa - a * exp( -bb * epo * 100 ) - ( alfa - a - 1 ) * exp( -c * epo * 100 ) );
       Tls_n = Cls_n = -fy * ( alfa - a * exp( -bb * epo * 100 ) - ( alfa - a - 1 ) * exp( -c * epo * 100 ) );
       //Elastic Local Buckling
       if(eplbf < 0.0){
         eplbf = 0.000001;
         fres_n = Tls_p * ( fres_n / fy );
         eresp_n = -( ( fres_n - fy ) / Ksft ) + 3.14 * ( pow ( (D/t) * sqrt(fy/Ee), -1.48 ) ) * (fy/Ee + 0.002) - (fy/Ee + 0.002);
       }
       Tbs_p = Cbs_p =  Rbso;
       Tbs_n = Cbs_n = -Rbso;
       Tmem_p = Cmem_p = Tls_p;
       Tmem_n = Cmem_n = Tls_n;
       Tvbs_p = Cvbs_p = 0.0;
       Tvbs_n = Cvbs_n = 0.0;
       elastic = Celastic = 0;
       plastic = Cplastic = 0;
       memory = Cmemory = 0;
       lb = Clb = 0;
       elb = Celb = 0;
       ep = Cep = 0.0;
       ep_ref = Cep_ref = 0.0;
       elb_ref = Celb_ref = 0.0;
       epmin = Cepmin = 0.0;
       epmax = Cepmax = 0.0;
       ebar_p = Cebar_p = 0.0;
       Ep = CEp = 0.0;
       Epo = CEpo = 0.0;
       W = CW = 0.0;
       delta = Cdelta = 0.0;
       delta_p = Cdelta_p = 0.0;
       delta_y = Cdelta_y = 0.0;
       delta_in = Cdelta_in = 0.0;
       deltap_in = Cdeltap_in = 0.0;
       h = Ch = 0.0;
       cbs = Ccbs = 0.0;
       elbf = Celbf = 0.0;
       slbf = Cslbf = 0.0;
       Trule = Crule = 0;
       Ttangent = Ctangent = Ee;
       Tstress = Cstress = 0.0;
       Tstrain = Cstrain = 0.0;
       strs = 0.0;
       tgnt = 0.0;
       Cstrs = 0.0;
       Ctgnt = 0.0;
       strain_inc = 0.0;
       repeat = 0;
       ld = Cld = 0;
       eunld = Ceunld = 0.0;
       funld = Cfunld = 0.0;
       Cmax_strs = 0.0;
       Tmax_strs = 0.0;
       buckled = 0;
       lbstrn = 0.0;
       lbstrs = 0.0;
}

RCFT_stlMaterial::RCFT_stlMaterial()
  :UniaxialMaterial(0,MAT_TAG_RCFT_stl),
  fy(0.0), fu(0.0), Ee(0.0), D(0.0), t(0.0), Ksft(0.0), eplbf(0.0), Rbso(0.0),
  Epoi(0.0), alfa(0.0), a(0.0), bb(0.0), c(0.0), w(0.0), ksi(0.0), e(0.0), fE(0.0),
  Rlso(0.0), CRlso(0.0), Rls(0.0), CRls(0.0), Tls_p(0.0), Cls_p(0.0), Tls_n(0.0), Cls_n(0.0),
  Tbs_p(0.0), Cbs_p(0.0), Tbs_n(0.0), Cbs_n(0.0), Tmem_p(0.0), Cmem_p(0.0), Tmem_n(0.0), Cmem_n(0.0),
  Tvbs_p(0.0), Cvbs_p(0.0), Tvbs_n(0.0), Cvbs_n(0.0), elastic(0), Celastic(0), plastic(0), Cplastic(0),
  memory(0), Cmemory(0), lb(0), Clb(0), elb(0), Celb(0), ep(0.0), Cep(0.0), ep_ref(0.0), Cep_ref(0.0), 
  elb_ref(0.0), Celb_ref(0.0), epmin(0.0), Cepmin(0.0), epmax(0.0), Cepmax(0.0), ebar_p(0.0), Cebar_p(0.0),
  Ep(0.0), CEp(0.0), Epo(0.0), CEpo(0.0), W(0.0), CW(0.0), delta(0.0),  Cdelta(0.0), delta_p(0.0),
  Cdelta_p(0.0), delta_y(0.0), Cdelta_y(0.0), delta_in(0.0), Cdelta_in(0.0), deltap_in(0.0), 
  Cdeltap_in(0.0), h(0.0), Ch(0.0), cbs(0.0), Ccbs(0.0), elbf(0.0), Celbf(0.0), slbf(0.0), Cslbf(0.0),
  Trule(0), Crule(0), Ttangent(0.0), Ctangent(0.0), Tstress(0.0), Cstress(0.0), Tstrain(0.0), Cstrain(0.0),
  strs(0.0), tgnt(0.0), strain_inc(0.0), repeat(0), ld(0), Cld(0), eunld(0.0), Ceunld(0.0), funld(0.0), Cfunld(0.0),
  fres_n(0.0), eresp_n(0.0), epo(0), Cstrs(0.0), Ctgnt(0.0), Cmax_strs(0.0), Tmax_strs(0.0), buckled(0), lbstrn(0.0), lbstrs(0.0)  
{
       //DOES NOTHING
}

RCFT_stlMaterial::~RCFT_stlMaterial()
{

}

int 
RCFT_stlMaterial::setTrialStrain(double strain, double strainRate)
{
#ifdef COMPOSITE_DEBUG
   ofstream output;
   output.open("pldata.dat",ios::app);
#endif
   /* INITIALIZE STATE VARIBALES TO THEIR VALUES AT THE MOST CONVERGED CONFIGURATION */	
   backtocommitStatevar();	
   /************************************************************************/
   /* CALCULATE STRAIN INCREMENT AND STRESS INCREMENT FOR THE CURRENT 	   */
   /* LOAD INCREMENT							   */
   /************************************************************************/
   /* CAUTION : UPDATING THE CURRENT STRAIN LEVEL IS DIFFERENT FOR OTHER MATERIAL MODELS IN OPENSEES */
   /* Tstrain = strain IF YOU ARE USING THIS MATERIAL MODEL WITH OTHER ELEMENTS IN OPENSEES          */
   Tstrain = Cstrain + strain;
   //Tstrain = strain;
   double stress_inc, Pstress, dep;
   strain_inc = Tstrain - Cstrain;
   /*CALCULATE TRIAL STRES */
   Pstress = Cstress + Ctangent * strain_inc;
   /* CHECK WHETHER THE STRESS POINT IS INSIDE THE LOADING SURFACE ( ELASTIC LOADING ) */
   if( ( Pstress - Tls_p < 0.0 ) && ( Pstress - Tls_n > 0.0 ) ){
	   /* ELASTIC LOADING */  
	   elastic = 1;
	   plastic = 0;
	   Trule = 1;
	   /* LOADING IN THE SOFTENING REGION IS PLASTIC */ 
	   if( ( lb == 1 ) && ( strain_inc <= 0.0 ) ){
		   plastic = 1;
		   elastic = 0;
		   Trule = 2;
	   }
   }
   else{
	   /* PLASTIC LOADING */ 
	   elastic = 0;
	   plastic = 1;
	   Trule = 2;
	   /* REVERSE LOADING IN THE SOFTENING REGION IS ELASTIC */
	   if( ( lb == 1 ) && ( strain_inc >= 0.0 ) ){
		   elastic = 1;
		   plastic = 0;
		   Trule = 1;
	   }
   }
   /* ELASTIC LOADING ( Stress point is within loading surfaces ) */
   if ( elastic == 1 ){
	   /* RESET LOCAL BUCKLING IN THE CASE OF ELASTIC LOADING */
	   lb = 0;
	   /* MAKE SURE REFERENCE PLASTIC STRAIN TO DETECT LOCAL BUCKLING REMAINS CORRECT DURING ELASTIC LOADING */
	   if( Cld != ld ){
		   ep_ref = ep;
	   }
	   else{
		   ep_ref = Cep_ref;
	   }	 
	   Ttangent = Ee;
	   stress_inc = strain_inc * Ttangent;
	   Tstress = Cstress + stress_inc;
	   if ( Tstress <=   Epo * ep + Tbs_n ){
		   Tstress =   Epo * ep + Tbs_n;
	   } 
	   if ( Tstress >=   Epo * ep + Tbs_p ){
		   Tstress =  Epo * ep + Tbs_p;
	   }   
	   /* IF ELASTIC LOADING TAKES PLACE RIGHT AFTER PLASTIC LOADING UPDATE STATE VARIABLES */
	   if( Crule == 2 ){
		   /* IF STRESS LEVEL DOES NOT REACH MEMORY LINE CALCULATE delta_y */	  
		   if( ( strain_inc < 0.0 ) && ( Cstress < Epo * ep + Tmem_p ) ){
			   delta_y = fabs( -Epo * ep - Tmem_p + Tls_p );
			   /* MEMORY SURFACE IS NOT BREACHED */
			   memory = 1;
		   }
		   /* IF STRESS LEVEL REACHES MEMORY LINE UPDATE MEMORY LINES */
		   else if( ( strain_inc < 0.0 ) && ( Cstress >= Epo * ep + Tmem_p ) ){
			   //Tmem_p =  - Epo * ep + Tls_p;
			   //Tmem_n = Tmem_p - 2 * ( Tmem_p - cbs );
			   Tmem_p = Epo * ep + Cmax_strs;
			   Tmem_n = - Tmem_p; 
			   /* MEMORY SURFACE IS BREACHED */
			   memory = 0;
		   }
		   /* IF STRESS LEVEL DOES NOT REACH MEMORY LINE CALCULATE delta_y */
		   if( ( strain_inc > 0.0 ) && ( Cstress > Epo * ep + Tmem_n ) ){
			   delta_y = fabs( -Epo * ep - Tmem_n + Tls_n );
			   /* MEMORY SURFACE IS NOT BREACHED */
			   memory = 1;
		   }
		   /* IF STRESS LEVEL REACHES MEMORY LINE UPDATE MEMORY LINES */
		   else if( ( strain_inc > 0.0 ) && ( Cstress <=  Epo * ep + Tmem_n ) ){
			   //Tmem_n = - Epo * ep + Tls_n;
			   //Tmem_p = Tmem_n - 2 * ( Tmem_n - cbs );
			   Tmem_n = -Epoi * ep - Cmax_strs;
			   Tmem_p = -Tmem_n; 
			   /* MEMORY SURFACE IS BREACHED */
			   memory = 0;
		   }
		   /* IF MEMORY SURFACE IS BREACHED DETERMINE BOUNDING SURFACES */
		   if ( ( strain_inc < 0.0 ) && ( memory == 0 ) ){
			   Tbs_n = - ( fu + ( Rbso - fu ) * exp ( - pow( fy / Ee, -2 )* ksi  * pow ( 0.5 * ebar_p, 2 ) ) );
			   if( buckled == 1 ){
				   Tbs_n = lbstrs;
				   if( lbstrs <  - ( -lbW * W + 1.0 ) * Rbso ){
					   Tbs_n = - ( -lbW * W + 1.0 ) * Rbso;
				   }
				   if( Tbs_n > -0.50 * Rbso ){
					   Tbs_n = -0.50 * Rbso;
				   }
			   }
			   if( Tbs_n <= -fu ){
				   Tbs_n = -fu;
			   }
			   cbs = Tbs_n + ( ( fabs( Tbs_n ) + fabs( Tbs_p ) ) / 2 );
			   cbs = 0.0; 
		   }
		   /* IF MEMORY SURFACE IS NOT BREACHED DETERMINE VIRTUAL BOUNDING SURFACES */
		   if ( ( strain_inc < 0.0 ) && ( memory == 1 ) ){
			   Tbs_n = - ( fu + ( Rbso - fu ) * exp ( - pow( fy / Ee, -2 )* ksi * pow ( 0.5 * ebar_p, 2 ) ) );
			   if( buckled == 1 ){
				   Tbs_n = lbstrs;
				   if( lbstrs < - ( -lbW * W + 1.0 ) * Rbso ){
					   Tbs_n = - ( -lbW * W + 1.0 ) * Rbso;
				   }
				   if( Tbs_n > -0.50 * Rbso ){
					   Tbs_n = -0.50 * Rbso;
				   }
			   }
			   if( Tbs_n <= -fu ){
				   Tbs_n = -fu;
			   }
			   Tvbs_n = Tbs_n - delta_y;
			   cbs = Tvbs_n + ( ( fabs( Tvbs_n ) + fabs( Tbs_p ) ) / 2 );
			   cbs = 0.0;
		   }
		   /* IF MEMORY SURFACE IS BREACHED DETERMINE BOUNDING SURFACES */
		   if ( ( strain_inc > 0.0 ) && ( memory == 0 ) ){
			   Tbs_p = fu + ( Rbso - fu ) * exp ( - pow( fy / Ee, -2 )* ksi * pow ( 0.5 * ebar_p, 2 ) );
			   if( buckled == 1 ){
				   Tbs_p = ( -lbW * W + 1.0 ) * Rbso; 
			   } 
			   if( Tbs_p < 0.50 * Rbso ){
				   Tbs_p = 0.50 * Rbso;
			   }
			   if( Tbs_p >= fu ){
				   Tbs_p = fu;
			   }
			   cbs = Tbs_p - ( ( fabs( Tbs_n ) + fabs( Tbs_p ) ) / 2 );
			   cbs = 0.0;
		   }
		   /* IF MEMORY SURFACE IS NOT BREACHED DETERMINE VIRTUAL BOUNDING SURFACES */ 
		   if ( ( strain_inc > 0.0 ) && ( memory == 1 ) ){
			   Tbs_p = fu + ( Rbso - fu ) * exp ( - pow( fy / Ee, -2 )* ksi * pow ( 0.5 * ebar_p, 2 ) );
			   if( buckled == 1 ){
				   Tbs_p = ( -lbW * W + 1.0 ) * Rbso; 
			   }
			   if( Tbs_p < 0.50 * Rbso ){
				   Tbs_p = 0.50 * Rbso;
			   }
			   if( Tbs_p >= fu ){
				   Tbs_p = fu;
			   }
			   Tvbs_p = Tbs_p + delta_y;
			   cbs = Tbs_n + ( ( fabs( Tbs_n ) + fabs( Tvbs_p ) ) / 2 );
			   cbs = 0.0;
		   }
		   /* UPDATE UNLOADING STRESS AND STRAIN */
		   funld = Cstress;
		   eunld = Cstrain;
	   }
   }
   /* PLASCTIC LOADING ( Stress point is outside of loading surfaces ) */
   if ( plastic == 1 ){
	   /* IMMEDIATELY AFTER PLASTIC LOADING, MEMORY SURFACE IS BREACHED, NO LOCAL BUCKLING OCCURED */  
	   if( ( Crule == 1 ) && ( memory == 0 ) && ( lb == 0 ) ){
		   if( strain_inc <= 0.0 ){
			   Epo = Epoi / ( 1 + w * W );
			   delta_in = fabs( Epo * ep + Tbs_n - Tls_n);
			   h = e * delta_in + fE * Ee;
			   Ep = Epo + h * ( delta_in ) / 0.000000000001;
			   Ttangent = Ee * Ep / ( Ee + Ep );
			   /* CORRECT THE OVERSHOOTING OF THE YIELD SURFACE RATHER THAN HAVING Tstress = Cstress + Ttangent * strain_incr */
			   Tstress = Tls_n + ( ( Pstress - Tls_n ) / Ee ) * Ttangent;
			   /* MAKE SURE THAT LOADING POINT DOES NOT BREACH BOUNDING SURFACE */
			   if ( Tstress <= Epo * ep + Tbs_n ){
				   Tstress = Epo * ep + Tbs_n;
			   } 
			   /* CALCULATE INCREMENTAL PLASTIC STRAIN */ 
			   dep = ( ( Pstress - Tls_n ) / Ee ) * Ttangent / Ep;
			   /* DETECT DIRECTION OF LOADING */
			   if(dep > 0.0){
				   ld = 1;
			   }
			   else if(dep < 0.0){
				   ld = -1;
			   }
			   /* IF LOADING DIRECTION IS NOT REVERSED DO NOT UPDATE THE STATE VRAIABLES */
			   if( ld == Cld ){
				   delta_in = Cdelta_in;
				   h = Ch;
				   delta = fabs( Epo * (ep+dep) + Tbs_n - Tls_n);
				   Ep = Epo + h * ( delta ) / (delta_in-delta);
				   Ttangent = Ee * Ep / ( Ee + Ep );
				   stress_inc = (Tstrain - eunld) * Ttangent;
				   Tstress = funld + stress_inc;
				   dep = stress_inc / Ep;
				   ep_ref = Cep_ref;
			   }
			   /* IF LOADING DIRECTION IS NOT REVERSED DO NOT UPDATE THE REFERENCE PLASTIC STRAIN FOR LOCAL BUCKLING*/
			   if( ld != Cld ){
				   ep_ref = ep;
				   //Tls_n = Tstress;
				   //Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 )
				   //  - ( alfa - a - 1 ) * exp( -c * ebar_p * 100 ) );
				   //Tls_p = Tstress + 2 * Rls;
			   }
			   ep = ep + dep;
			   /* CHECK LOCAL BUCKLING */
			   //if ( ( ep - ep_ref <= - eplbf ) && ( Tstress < -fy  ) && ( ep < 0.0 )  ){
			   //if ( ( ep - ep_ref <= - eplbf ) && ( ( Tstress < -fy  ) | ( Tstress < Tls_n ) ) ){
			   //     lb = 1;
			   //     buckled = 1;
			   //     Tbs_n = - fy;
			   //}
			   /* CHECK LOCAL BUCKLING */
			   if ( ( ep - ep_ref <= - eplbf ) && ( ( Tstrain < lbstrn ) && ( Tstress < 0.99 * lbstrs ) ) && ( buckled == 1 ) ){
				   lb = 1;
				   //buckled = 1;
			   }
			   if ( ( ep - ep_ref <= - eplbf ) && ( Tstress < Tls_n ) && ( buckled == 0 ) && ( Tstress < -fres_n ) ){
				   lb = 1;
				   buckled = 1;
			   }
			   /* UPDATE AND CALCULATE STATE VARIABLES */
			   if( ep >= epmax ){
				   epmax = ep;
			   }
			   if( ep <= epmin ){ 
				   epmin = ep;
			   }
			   ebar_p = fabs( epmin ) + fabs ( epmax );
			   Tls_n = Tstress;
			   Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 )
				   - ( alfa - a - 1 ) * exp( -c * ebar_p * 100 ) );
			   Tls_p = Tstress + 2 * Rls;
			   W = W + Tstress * dep;
		   } 
		   if( strain_inc >= 0.0 ){
			   Epo = Epoi / ( 1 + w * W );
			   delta_in = fabs( Epo * ep + Tbs_p - Tls_p) ;
			   h = e * delta_in + fE * Ee;
			   Ep = Epo + h * ( delta_in ) / 0.000000000001;
			   Ttangent = Ee * Ep / ( Ee + Ep );
			   /* CORRECT THE OVERSHOOTING OF THE YIELD SURFACE RATHER THAN HAVING Tstress = Cstress + Ttangent * strain_incr */
			   Tstress = Tls_p + ( ( Pstress - Tls_p ) / Ee ) * Ttangent;
			   /* MAKE SURE THAT LOADING POINT DOES BREACH BOUNDING SURFACE */
			   if ( Tstress >= Epo * ep + Tbs_p ){
				   Tstress = Epo * ep + Tbs_p;
			   }
			   /* CALCULATE INCREMENTAL PLASTIC STRAIN */
			   dep = ( ( Pstress - Tls_p ) / Ee ) * Ttangent / Ep;
			   /* DETECT DIRECTION OF LOADING */
			   if(dep > 0.0){
				   ld = 1;
			   }
			   else if(dep < 0.0){
				   ld = -1;
			   }
			   /* IF LOADING DIRECTION IS NOT REVERSED DO NOT UPDATE THE STATE VRAIABLES */
			   if( ld == Cld ){
				   delta_in = Cdelta_in;
				   h = Ch;
				   delta = fabs( Epo * (ep+dep) + Tbs_p - Tls_p);
				   Ep = Epo + h * ( delta ) / (delta_in-delta);
				   Ttangent = Ee * Ep / ( Ee + Ep );
				   stress_inc = (Tstrain - eunld) * Ttangent;
				   Tstress = funld + stress_inc;
				   dep = stress_inc / Ep;
				   ep_ref = Cep_ref;                   
			   }
			   /* IF LOADING DIRECTION IS NOT REVERSED DO NOT UPDATE THE REFERENCE PLASTIC STRAIN FOR LOCAL BUCKLING*/
			   if( ld != Cld ){
				   ep_ref = ep;
				   //Tls_p = Tstress;
				   //Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 )
				   // - ( alfa - a -1 ) * exp( -c * ebar_p * 100 ) );
				   //Tls_n = Tstress - 2 * Rls;
			   }
			   ep = ep + dep;
			   /* UPDATE AND CALCULATE STATE VARIABLES */
			   if( ep >= epmax ){
				   epmax = ep;
			   }
			   if( ep <= epmin ){
				   epmin = ep;
			   }
			   ebar_p = fabs( epmin ) + fabs ( epmax );
			   Tls_p = Tstress;
			   Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 )
				   - ( alfa - a -1 ) * exp( -c * ebar_p * 100 ) );
			   Tls_n = Tstress - 2 * Rls;
			   W = W + Tstress * dep;
		   } 
	   } /* if( (Crule == 1 ) && ( memory == 0 ) ) 		   */
	   /* IMMEDIATELY AFTER PLASTIC LOADING, MEMORY SURFACE IS NOT BREACHED, NO LOCAL BUCKLING OCCURED */
	   else if ( ( Crule == 1 ) && ( memory == 1 ) && ( lb == 0 ) ){
		   if( strain_inc <= -0.0 ){
			   Epo = Epoi / ( 1 + w * W );
			   /* CALCULATE delta_in with respect to VIRTUAL BOUNDING SURFACE */
			   deltap_in = fabs( Epo * ep + Tvbs_n - Tls_n);
			   delta_in = fabs ( Epo * ep + Tbs_n - Tls_n);
			   h = e * deltap_in + fE * Ee;
			   Ep = Epo + h * ( delta_in ) / 1.0E-30;
			   Ttangent = Ee * Ep / ( Ee + Ep );
			   /* CORRECT THE OVERSHOOTING OF THE YIELD SURFACE RATHER THAN HAVING Tstress = Cstress + Ttangent * strain_incr */
			   //Tstress = Tls_n + ( ( Pstress - Tls_n ) / Ee )* Ttangent;
			   Tstress = Cstress + Ttangent * strain_inc;
			   /* MAKE SURE THAT LOADING POINT DOES BREACH BOUNDING SURFACE */
			   if ( Tstress <= Epo * ep + Tbs_n ){
				   Tstress = Epo * ep + Tbs_n;
			   }
			   /* CALCULATE INCREMENTAL PLASTIC STRAIN */
			   dep = ( ( Pstress - Tls_n ) / Ee ) * Ttangent / Ep;
			   /* DETECT DIRECTION OF LOADING */
			   if(dep > 0.0){
				   ld = 1;
			   }
			   else if(dep < 0.0){
				   ld = -1;
			   }
			   /* IF LOADING DIRECTION IS NOT REVERSED DO NOT UPDATE THE STATE VRAIABLES */
			   if( ld == Cld ){
				   deltap_in = Cdeltap_in;
				   delta_in = Cdelta_in;
				   h = Ch;
				   delta = fabs( Epo * (ep+dep) + Tbs_n - Tls_n);
				   memory = Cmemory;
				   if( fabs( delta_in-delta ) < 1.0E-30 ){
					   Ep =  Epo + h * ( delta ) / 1.0E-30;
				   }
				   else{
					   Ep = Epo + h * ( delta ) / (delta_in-delta);
				   }
				   Ttangent = Ee * Ep / ( Ee + Ep );
				   stress_inc = ( Tstrain - eunld ) * Ttangent;
				   Tstress = funld + stress_inc;
				   dep = stress_inc / Ep;
			   }
			   /* IF LOADING DIRECTION IS NOT REVERSED DO NOT UPDATE THE REFERENCE PLASTIC STRAIN FOR LOCAL BUCKLING*/
			   if( ld == Cld ){
				   ep_ref = Cep_ref;
			   }
			   else{
				   ep_ref = ep;
			   }
			   ep = ep + dep;
			   ///* CHECK LOCAL BUCKLING */
			   //if ( ( ep - ep_ref <= -eplbf ) && ( Tstress < -fy ) && ( Tstrain < 0.0 )){
			   //if ( ( ep - ep_ref <= -eplbf ) && ( ( Tstress < -fy ) | ( Tstress < Tbs_n ) ) ){
			   //   lb = 1;
			   //   buckled = 1; 
			   //   Tbs_n = -fy;
			   //}
			   /* CHECK LOCAL BUCKLING */
			   if ( ( ep - ep_ref <= - eplbf ) && ( ( Tstrain < lbstrn ) && ( Tstress < 0.99 * lbstrs ) ) && ( buckled == 1 ) ){
				   lb = 1;
				   //buckled = 1;
			   }
			   if ( ( ep - ep_ref <= - eplbf ) && ( Tstress < Tls_n ) && ( buckled == 0 ) && ( Tstress < -fres_n ) ){
				   lb = 1;
				   buckled = 1;
			   }
			   /* UPDATE AND CALCULATE STATE VARIABLES */
			   if( ep >= epmax ){
				   epmax = ep;
			   }
			   if( ep <= epmin ){
				   epmin = ep;
			   }
			   ebar_p = fabs( epmin ) + fabs ( epmax );
			   Tls_n = Tstress;
			   Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 )
				   - ( alfa - a -1 ) * exp( -c * ebar_p * 100 ) );
			   Tls_p = Tstress + 2 * Rls;
			   W = W + Tstress * dep;
		   } 
		   if( strain_inc >= 0.0 ){
			   Epo = Epoi / ( 1 + w * W );
			   /* CALCULATE delta_in with respect to VIRTUAL BOUNDING SURFACE */
			   deltap_in = fabs( Epo * ep + Tvbs_p - Tls_p) ;
			   delta_in = fabs( Epo * ep + Tbs_p - Tls_p) ;
			   h = e * deltap_in + fE * Ee;
			   Ep = Epo + h * ( deltap_in ) / 1.0E-30;
			   Ttangent = Ee * Ep / ( Ee + Ep );
			   /* CORRECT THE OVERSHOOTING OF THE YIELD SURFACE RATHER THAN HAVING Tstress = Cstress + Ttangent * strain_incr */
			   //Tstress = Tls_p + ( ( Pstress - Tls_p ) / Ee )* Ttangent;
			   Tstress = Cstress + Ttangent * strain_inc;
			   /* MAKE SURE THAT LOADING POINT DOES BREACH BOUNDING SURFACE */
			   if ( Tstress >= Epo * ep + Tbs_p ){
				   Tstress = Epo * ep + Tbs_p;
			   }
			   /* CALCULATE INCREMENTAL PLASTIC STRAIN */
			   dep = ( ( Pstress - Tls_p ) / Ee ) * Ttangent / Ep;
			   /* DETECT DIRECTION OF LOADING */
			   if(dep > 0.0){
				   ld = 1;
			   }
			   else if(dep < 0.0){
				   ld = -1;
			   }
			   /* IF LOADING DIRECTION IS NOT REVERSED DO NOT UPDATE THE STATE VARIABLES */
			   if( ld == Cld ){
				   deltap_in = Cdeltap_in;
				   delta_in = Cdelta_in;
				   h = Ch;
				   delta = fabs( Epo * (ep+dep) + Tbs_p - Tls_p);
				   memory = Cmemory;
				   if( fabs( delta_in-delta ) < 1.0E-30 ){
					   Ep =  Epo + h * ( delta ) / 1.0E-30;
				   }
				   else{
					   Ep = Epo + h * ( delta ) / (delta_in-delta);
				   }
				   Ttangent = Ee * Ep / ( Ee + Ep );
				   stress_inc = ( Tstrain - eunld) * Ttangent;
				   Tstress = funld + stress_inc;
				   dep = stress_inc / Ep;
			   }
			   /* IF LOADING DIRECTION IS NOT REVERSED DO NOT UPDATE THE REFERENCE PLASTIC STRAIN FOR LOCAL BUCKLING*/
			   if( ld == Cld ){
				   ep_ref = Cep_ref;
			   }
			   else{
				   ep_ref = ep;
			   }
			   /* UPDATE AND CALCULATE STATE VARIABLES */
			   ep = ep + dep;
			   if( ep >= epmax ){
				   epmax = ep;
			   }
			   if( ep <= epmin ){
				   epmin = ep;
			   }
			   ebar_p = fabs( epmin ) + fabs ( epmax );
			   Tls_p = Tstress;
			   Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 )
				   - ( alfa - a -1 ) * exp( -c * ebar_p * 100 ) );
			   Tls_n = Tstress - 2 * Rls;
			   W = W + Tstress * dep;
		   } 
	   } /* else if ( ( Crule == 1 ) && ( memory == 1 ) )*/   
	   /* COUNTINUED PLASTIC LOADING, MEMORY SURFACE BREACHED, NO LOCAL BUCKLING */ 
	   else if  ( ( Crule == 2 ) && ( memory == 0 ) && ( lb == 0 ) ){
		   if( Cstress <= Tls_n ){
			   delta = fabs( Epo * ep + Tbs_n - Cstress) ;
			   if( ( Epo * ep + Tbs_n - Cstress ) > 0.0 ){ 
				   delta = 0.0;
			   }
		   } 
		   if( Cstress >= Tls_p ){
			   delta = fabs( Epo * ep + Tbs_p - Cstress);
			   if( ( Epo * ep + Tbs_p - Cstress ) < 0.0 ){
				   delta = 0.0;
			   }
		   } 
		   h = e * delta + fE * Ee;
		   Ep = Epo + h * ( delta ) / ( delta_in - delta );
		   Ttangent = Ee * Ep / ( Ee + Ep );
		   stress_inc = strain_inc * Ttangent;
		   Tstress = Cstress + stress_inc;
		   if ( Tstress >=  Epo * ep + Tbs_p ){
			   Tstress =  Epo * ep + Tbs_p;
		   }
		   if ( Tstress <=  Epo * ep + Tbs_n ){
			   Tstress = Epo * ep + Tbs_n;
		   }
		   dep = stress_inc / Ep;
		   ep = ep + dep;
		   if( strain_inc < 0 ){
			   //if ( ( ep - ep_ref <= -eplbf ) && ( Tstress < -fy ) && ( Tstrain < 0.0 ) ){
			   //if ( ( ep - ep_ref <= -eplbf ) && ( ( Tstress < -fy ) | ( Tstress < Tbs_n ) ) ) {
			   //   lb = 1;
			   //   buckled = 1;
			   //   Tbs_n = -fy;
			   //}
			   if ( ( ep - ep_ref <= - eplbf ) && ( ( Tstrain < lbstrn ) && ( Tstress < 0.99 * lbstrs ) ) && ( buckled == 1 ) ){
				   lb = 1;
				   //buckled = 1;
			   }
			   if ( ( ep - ep_ref <= - eplbf ) && ( Tstress < Tls_n ) && ( buckled == 0 ) && ( Tstress < -fres_n )){
				   lb = 1;
				   buckled = 1;
			   }
		   }
		   if( ep >= epmax ){
			   epmax = ep;
		   } 
		   if( ep <= epmin ){
			   epmin = ep;
		   } 
		   ebar_p = fabs( epmin ) + fabs ( epmax );
		   if( strain_inc > 1.0E-30 ){
			   Tls_p = Tstress;
			   Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 )
				   - ( alfa - a -1 ) * exp( -c * ebar_p * 100 ) );
			   Tls_n = Tstress - 2 * Rls;
		   } 
		   if( strain_inc < -1.0E-30 ){
			   Tls_n = Tstress;
			   Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 )
				   - ( alfa - a -1 ) * exp( -c * ebar_p * 100 ) );
			   Tls_p = Tstress + 2 * Rls;
		   } 
		   W =W + Tstress * dep;
		   elb_ref = Tstrain;
	   } /* else if ( Crule == 2 ) && ( memory == 0 ) */
	   else if ( ( Crule == 2 ) && ( memory == 1 ) && ( lb == 0 ) ){
		   if( ( Cstress >=  Epo * ep + Tmem_p ) | ( Cstress <= Epo * ep + Tmem_n ) )
		   {
			   memory = 0;
			   if( Cstress <= Tls_n ){
				   delta = fabs( Epo * ep + Tbs_n - Cstress) ;
				   if( ( Epo * ep + Tbs_n - Cstress ) > 0.0 ){ 
					   delta = 0.0;
				   }
			   }
			   if( Cstress >= Tls_p ){
				   delta = fabs( Epo * ep + Tbs_p - Cstress);
				   if( ( Epo * ep + Tbs_p - Cstress ) < 0.0 ){
					   delta = 0.0;
				   }
			   }
			   h = e * delta + fE * Ee;
			   Ep = Epo + h * ( delta ) / ( delta_in - delta );
			   Ttangent = Ee * Ep / ( Ee + Ep );
			   stress_inc = strain_inc * Ttangent;
			   Tstress = Cstress + stress_inc;
			   if ( Tstress >=  Epo * ep + Tbs_p ){
				   Tstress =  Epo * ep + Tbs_p;
			   }
			   if ( Tstress <=  Epo * ep + Tbs_n ){
				   Tstress =  Epo * ep + Tbs_n;
			   } 
			   dep = stress_inc / Ep;
			   ep = ep + dep;
			   if( strain_inc < 0.0 ){
				   //if ( ( ep <= -eplbf ) && ( Tstress < -fy ) && ( Tstrain < 0.0 )){
				   //if ( ( ep <= -eplbf ) && ( ( Tstress < -fy ) | ( Tstress < Tbs_n ) ) ){
				   //   lb = 1;
				   //   buckled = 1;
				   //   Tbs_n = -fy;
				   //}
				   if ( ( ep - ep_ref <= - eplbf ) && ( ( Tstrain < lbstrn ) && ( Tstress < 0.99 * lbstrs ) ) && ( buckled == 1 ) ){
					   lb = 1;
					   //buckled = 1;
				   }
				   if ( ( ep - ep_ref <= - eplbf ) && ( Tstress < Tls_n ) && ( buckled == 0 ) && ( Tstress < -fres_n ) ){
					   lb = 1;
					   buckled = 1;
				   }
			   }
			   if( ep >= epmax ){
				   epmax = ep;
			   }
			   if( ep <= epmin ){
				   epmin = ep;
			   }
			   ebar_p = fabs( epmin ) + fabs ( epmax );
			   if( strain_inc > 1.0E-20 ){
				   Tls_p = Tstress;
				   Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 )
					   - ( alfa - a -1 ) * exp( -c * ebar_p * 100 ) );
				   Tls_n = Tstress - 2 * Rls;
			   }
			   if( strain_inc < -1.0E-20 ){
				   Tls_n = Tstress;
				   Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 )
					   - ( alfa - a -1 ) * exp( -c * ebar_p * 100 ) );
				   Tls_p = Tstress + 2 * Rls;
			   }
			   W = W + Tstress * dep;
		   } /* if( ( Cstress <= Tmem_p ) && ( Cstress >= Tmem_n ) ) */
		   else if( ( Cstress < Epo * ep + Tmem_p ) && ( Cstress > Epo * ep + Tmem_n ) ){
			   if( Cstress <= Tls_n ){
				   delta = fabs( Epo * ep + Tbs_n - Cstress) ;
				   delta_p = fabs( Epo * ep + Tvbs_n - Cstress) ;
				   if( ( Epo * ep + Tvbs_n - Tstress ) > 0.0 ){
					   delta = 0.0;
				   }
			   }
			   if( Cstress >= Tls_p ){
				   delta = fabs( Epo * ep + Tbs_p - Cstress);
				   delta_p = fabs( Epo * ep + Tvbs_p - Cstress) ;
				   if( ( Epo * ep + Tvbs_p - Cstress ) < 0.0 ){
					   delta = 0.0;
				   }
			   }
			   h = e * delta + fE * Ee;
			   Ep = Epo + h * ( delta + delta_y ) / ( delta_in - delta );
			   Ttangent = Ee * Ep / ( Ee + Ep );
			   stress_inc = strain_inc * Ttangent;
			   Tstress = Cstress + stress_inc;
			   if ( Tstress >=  Epo * ep + Tbs_p ){
				   Tstress =  Epo * ep + Tbs_p;
			   }
			   if ( Tstress <=  Epo * ep + Tbs_n ){
				   Tstress =  Epo * ep + Tbs_n;
			   }
			   dep = stress_inc / Ep;
			   ep = ep + dep;
			   if( strain_inc < 0.0 ){
				   //if( ( ep - ep_ref <= -eplbf ) && ( Tstress < -fy ) && ( Tstrain < 0.0 ) ){
				   //if( ( ep - ep_ref <= -eplbf ) && ( ( Tstress < -fy ) | ( Tstress < Tbs_n ) ) ){
				   //  lb = 1;
				   //  buckled = 1;
				   //  Tbs_n = -fy;
				   //}
				   if ( ( ep - ep_ref <= - eplbf ) && ( ( Tstrain < lbstrn ) && ( Tstress < 0.99 * lbstrs ) ) && ( buckled == 1 ) ){
					   lb = 1;
					   //buckled = 1;
				   }
				   if ( ( ep - ep_ref <= - eplbf ) && ( Tstress < Tls_n ) && ( buckled == 0 ) && ( Tstress < -fres_n  ) ){
					   lb = 1;
					   buckled = 1;
				   }
			   }
			   if( ep >= epmax ){
				   epmax = ep;
			   }
			   if( ep <= epmin ){
				   epmin = ep;
			   }
			   ebar_p = fabs( epmin ) + fabs ( epmax );
			   if( strain_inc > 1.0E-20 ){
				   Tls_p = Tstress;
				   Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 )
					   - ( alfa - a -1 ) * exp( -c * ebar_p * 100 ) );
				   Tls_n = Tstress - 2 * Rls;
			   }
			   if( strain_inc < -1.0E-20 ){
				   Tls_n = Tstress;
				   Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 )
					   - ( alfa - a -1 ) * exp( -c * ebar_p * 100 ) );
				   Tls_p = Tstress + 2 * Rls;
			   }
			   W = W + Tstress * dep;
			   elb_ref = Tstrain;
		   } /* else if( ( Cstress <= Tmem_p ) && ( Cstress >= Tmem_n ) ) */
	   } /* else if ( ( Crule == 2 ) && ( memory == 1 ) ) */
	   else if ( lb == 1 ){
		   if( ep - ep_ref > -eresp_n ){
			   Ttangent = -Ksft;
			   Ep = Ksft * Ee / ( Ee + Ksft );
			   stress_inc = strain_inc * Ttangent;
			   Tstress = Cstress + stress_inc;
			   if ( Tstrain < lbstrn ){
				   lbstrn = Tstrain;
				   lbstrs = Tstress;
			   }
			   dep = -(stress_inc / Ep);
			   ep = ep + dep;
		   }
		   else if( ep - ep_ref < -eresp_n ){
			   Ttangent = 0.1;
			   Ep = 0.1;
			   Tstress = Cstress;
			   if ( Tstrain < lbstrn ){
				   lbstrn = Tstrain;
				   lbstrs = Tstress;
			   }
			   dep = strain_inc;
			   ep = ep + dep;
		   }
		   if( ep >= epmax ){
			   epmax = ep;
		   }
		   if( ep <= epmin ){
			   epmin = ep;
		   }
		   W = W + Tstress * dep;
		   ebar_p = fabs( epmin ) + fabs ( epmax );
		   if( strain_inc < 0.0 ){
			   Tls_n = Tstress;
			   Rls = Rlso * ( alfa - a * exp( -bb * ebar_p * 100 )
				   - ( alfa - a -1 ) * exp( -c * ebar_p * 100 ) );
			   Tls_p = Tstress + 2 * Rls;
		   }
	   } /* else if ( lb == 1) */
   } /* if ( plastic == 1 ) */
   strs = Tstress;
   tgnt = Ttangent;
   if( Tstress > Cmax_strs ){
	   Tmax_strs = Tstress;
   }
   return 0;
}
double 
RCFT_stlMaterial::getStress(void)
{
  return strs;
}

double 
RCFT_stlMaterial::getTangent(void)
{
  return tgnt;
}

double 
RCFT_stlMaterial::getInitialTangent(void)
{
  return Ee;
}

double 
RCFT_stlMaterial::getStrain(void)
{
  return Tstrain;
}

int 
RCFT_stlMaterial::commitState(void)
{
#ifdef COMPOSITE_DEBUG
  ofstream output;
  output.open("pldata.dat",ios::app);
#endif
  Cstrain = Tstrain;
  Cstress = Tstress;
  Ctangent = Ttangent;
  Crule = Trule;
  Cstrs = strs;
  Ctgnt = tgnt;
  /* COMMIT STATE VARIABLES */
  commitStatevar();
  return 0;
}

int 
RCFT_stlMaterial::revertToLastCommit(void)
{
#ifdef COMPOSITE_DEBUG
  ofstream output;
  output.open("pldata.dat",ios::app);
#endif
  Tstrain = Cstrain;
  Tstress = Cstress;
  Ttangent = Ctangent;
  Trule = Crule;
  strs = Cstrs;
  tgnt = Ctgnt;
#ifdef COMPOSITE_DEBUG
  output<<Tstrain<<"   "<<Tstress<<"   "<<Ttangent<<"   "<<Cstrain<<"   "<<Cstress<<"   "<<Ctangent<<"  "<<Ksft<<"  "<<(D/t) * (fy/Ee)<<"  "<<fres_n<<endl;
#endif
  this->backtocommitStatevar();
  return 0;
}

int 
RCFT_stlMaterial::revertToStart(void)
{
  return lb;
}

void
RCFT_stlMaterial::commitStatevar(void)
{
   CRlso = Rlso;
   CRls = Rls;
   Cls_p = Tls_p;
   Cls_n = Tls_n;
   Cbs_p = Tbs_p;
   Cbs_n = Tbs_n;
   Cmem_p = Tmem_p;
   Cmem_n = Tmem_n;
   Cvbs_p = Tvbs_p;
   Cvbs_n = Tvbs_n;
   Celastic = elastic;
   Cplastic = plastic;
   Cmemory = memory;
   Clb = lb;
   Celb = elb;
   Cep = ep;
   Cep_ref = ep_ref;
   Celb_ref = elb_ref;
   Cepmin = epmin;
   Cepmax = epmax;
   Cebar_p = ebar_p;
   CEp = Ep;
   CEpo = Epo;
   CW = W;
   Cdelta = delta;
   Cdelta_p = delta_p;
   Cdelta_y = delta_y;
   Cdelta_in = delta_in;
   Cdeltap_in = deltap_in;
   Cld = ld;
   Ch = h;
   Ccbs = cbs;
   Ceunld = eunld;
   Cfunld = funld; 
}
void
RCFT_stlMaterial::backtocommitStatevar(void)
{
   Rlso = CRlso;
   Rls = CRls;
   Tls_p = Cls_p;
   Tls_n = Cls_n;
   Tbs_p = Cbs_p;
   Tbs_n = Cbs_n;
   Tmem_p = Cmem_p;
   Tmem_n = Cmem_n;
   Tvbs_p = Cvbs_p;
   Tvbs_n = Cvbs_n;
   elastic = Celastic;
   plastic = Cplastic;
   memory = Cmemory;
   lb = Clb;
   elb = Celb;
   ep = Cep;
   ep_ref = Cep_ref;
   elb_ref = Celb_ref;
   epmin = Cepmin;
   epmax = Cepmax;
   ebar_p = Cebar_p;
   Ep = CEp;
   Epo = CEpo;
   W = CW;
   delta = Cdelta;
   delta_p = Cdelta_p;
   delta_y = Cdelta_y;
   delta_in = Cdelta_in;
   deltap_in = Cdeltap_in;
   ld = Cld;
   h = Ch;
   cbs = Ccbs;
   eunld = Ceunld;
   funld = Cfunld;
}	

UniaxialMaterial *
RCFT_stlMaterial::getCopy(void)
{
  RCFT_stlMaterial *theCopy = new RCFT_stlMaterial(this->getTag(), fy, fu, Ee, D, t, epo);
  theCopy->fy = fy; theCopy->fu = fu; theCopy->Ee = Ee; theCopy->D = D; theCopy->t = t;
  theCopy->Ksft = Ksft; theCopy->eplbf = eplbf; theCopy->Rbso = Rbso; theCopy->Epoi = Epoi; theCopy->alfa = alfa; 
  theCopy->a = a; theCopy->bb = bb; theCopy->c = c; theCopy->w = w; theCopy->ksi = ksi; theCopy->e = e; theCopy->fE = fE;
  theCopy->Rlso = Rlso; theCopy->CRlso = CRlso; theCopy->Rls = Rls; theCopy->CRls = CRls; theCopy->Tls_p = Tls_p; theCopy->Cls_p = Cls_p;
  theCopy->Tls_n = Tls_n; theCopy->Cls_n = Cls_n; theCopy->Tbs_p = Tbs_p; theCopy->Cbs_p = Cbs_p;
  theCopy->Tbs_n = Tbs_n; theCopy->Cbs_n = Cbs_n; theCopy->Tmem_p = Tmem_p; theCopy->Cmem_p = Cmem_p;
  theCopy->Tmem_n = Tmem_n; theCopy->Cmem_n = Cmem_n; theCopy->Tvbs_p = Tvbs_p; theCopy->Cvbs_p = Cvbs_p; 
  theCopy->Tvbs_n = Tvbs_n, theCopy->Cvbs_n = Cvbs_n;
  theCopy->elastic = elastic; theCopy->Celastic = Celastic; theCopy->plastic = plastic; theCopy->Cplastic = Cplastic;
  theCopy->memory = memory; theCopy->Cmemory = Cmemory; theCopy->lb = lb; theCopy->Clb = Clb; theCopy->elb = elb; theCopy->Celb = Celb;
  theCopy->ep = ep; theCopy->Cep = Cep; theCopy->ep_ref = ep_ref; theCopy->Cep_ref = Cep_ref; theCopy->elb_ref = elb_ref;
  theCopy->Celb_ref = Celb_ref; theCopy->epmin = epmin; theCopy->Cepmin = Cepmin; theCopy->epmax = epmax; theCopy->Cepmax = Cepmax;
  theCopy->ebar_p = ebar_p; theCopy->Cebar_p = Cebar_p; theCopy->Ep = Ep; theCopy->CEp = CEp; theCopy->Epo = Epo;
  theCopy->CEpo = CEpo; theCopy->W = W; theCopy->CW = CW; theCopy->delta = delta; theCopy->Cdelta = Cdelta;
  theCopy->delta_p = delta_p; theCopy->Cdelta_p = Cdelta_p;
  theCopy->delta_y = delta_y; theCopy->Cdelta_y = Cdelta_y; theCopy->delta_in = delta_in; 
  theCopy->Cdelta_in = Cdelta_in; theCopy->deltap_in = deltap_in; theCopy->Cdeltap_in = Cdeltap_in;
  theCopy->h = h; theCopy->Ch = Ch; theCopy->cbs = cbs; theCopy->Ccbs = Ccbs;
  theCopy->elbf = elbf; theCopy->Celbf = Celbf; theCopy->slbf = slbf; theCopy->Cslbf = Cslbf;
  theCopy->Trule = Trule; theCopy->Crule = Crule; theCopy->Ttangent = Ttangent; theCopy->Ctangent = Ctangent;
  theCopy->Tstress = Tstress; theCopy->Cstress = Cstress;
  theCopy->Tstrain = Tstrain; theCopy->Cstrain = Cstrain; theCopy->strs = strs; theCopy->strn = strn; theCopy->strain_inc = strain_inc;
  theCopy->repeat = repeat; theCopy->eunld = eunld; theCopy->funld = funld; theCopy->Ceunld = Ceunld; theCopy->Cfunld = Cfunld;
  theCopy->fres_n = fres_n; theCopy->eresp_n = eresp_n; theCopy->epo = epo;
  return theCopy;
}
int 
RCFT_stlMaterial::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(102);
  data(0) = this->getTag();
  data(1) = fy; data(2) = fu; data(3) = Ee; data(4) = D; data(5) = t;
  data(6) = Ksft; data(7) = eplbf; data(8) = Rbso; data(9) = Epoi; data(10) = alfa; 
  data(11) = a; data(12) = bb; data(13) = c; data(14) = w; data(15) = ksi; data(16) = e; data(17) = fE;
  data(18) = Rlso; data(19) = CRlso; data(20) = Rls; data(21) = CRls; data(22) = Tls_p; data(23) = Cls_p;
  data(24) = Tls_n; data(25) = Cls_n; data(26) = Tbs_p; data(27) = Cbs_p;
  data(28) = Tbs_n; data(29) = Cbs_n; data(30) = Tmem_p; data(31) = Cmem_p;
  data(32) = Tmem_n; data(33) = Cmem_n; data(34) = Tvbs_p; data(35) = Cvbs_p; 
  data(36) = Tvbs_n, data(37) = Cvbs_n;
  data(38) = elastic; data(39) = Celastic; data(40) = plastic; data(41) = Cplastic;
  data(42) = memory; data(43) = Cmemory; data(44) = lb; data(45) = Clb; data(46) = elb; data(47) = Celb;
  data(48) = ep; data(49) = Cep; data(50) = ep_ref; data(51) = Cep_ref; data(52) = elb_ref;
  data(53) = Celb_ref; data(54) = epmin; data(55) = Cepmin; data(56) = epmax; data(57) = Cepmax;
  data(54) = ebar_p; data(58) = Cebar_p; data(59) = Ep; data(60) = CEp; data(61) = Epo;
  data(62) = CEpo; data(63) = W; data(64) = CW; data(65) = delta; data(66) = Cdelta;
  data(67) = delta_p; data(68) = Cdelta_p;
  data(69) = delta_y; data(70) = Cdelta_y; data(71) = delta_in; 
  data(72) = Cdelta_in; data(73) = deltap_in; data(74) = Cdeltap_in;
  data(75) = h; data(76) = Ch; data(77) = cbs; data(78) = Ccbs;
  data(79) = elbf; data(80) = Celbf; data(81) = slbf; data(82) = Cslbf;
  data(83) = Trule; data(84) = Crule; data(85) = Ttangent; data(86) = Ctangent;
  data(87) = Tstress; data(88) = Cstress;
  data(89) = Tstrain; data(90) = Cstrain; data(91) = strs; data(92) = strn; data(93) = strain_inc; 
  data(94) = repeat; data(95) = eunld; data(96) = funld; data(97) = Ceunld; data(98) = Cfunld;
  data(99) = fres_n; data(100) = eresp_n; data(101) = epo;
  
  int res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if(res < 0)
  opserr<<"RCFT_stlMaterial::sendSelf() - failed to send data\n";
  return res;
}

int 
RCFT_stlMaterial::recvSelf(int cTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  static Vector data(102);
  int res = theChannel.recvVector(this->getDbTag(),cTag, data);
  if(res<0) {
    opserr<<"RCFT_stlMaterial::recvSelf() - failed to receive data \n";
    this->setTag(0);
  }
  else{
    this->setTag((int)data(0));
    fy = data(0);  fu = data(1); Ee = data(2); D = data(3); t = data(4);
    Ksft = data(5); eplbf = data(6); Rbso = data(7); Epoi = data(8); alfa = data(9); 
    a = data(10); bb = data(11); c = data(12); w = data(13);  ksi = data(14); e = data(15); fE = data(16);
    Rlso = data(17); CRlso = data(18); Rls = data(19); CRls = data(20); Tls_p = data(21); Cls_p = data(22);
    Tls_n = data(23); Cls_n = data(24); Tbs_p = data(25); Cbs_p = data(26);
    Tbs_n = data(27); Cbs_n = data(28); Tmem_p = data(29); Cmem_p = data(30);
    Tmem_n = data(31); Cmem_n = data(32); Tvbs_p = data(33); Cvbs_p = data(34); 
    Tvbs_n = data(35), Cvbs_n = data(36);
    elastic = data(37); Celastic = data(38); plastic = data(39); Cplastic = data(40);
    memory = data(41); Cmemory = data(42); lb = data(43); Clb = data(44); elb = data(45); Celb = data(46);
    ep = data(47); Cep = data(48); ep_ref = data(49); Cep_ref = data(50); elb_ref = data(51);
    Celb_ref = data(52); epmin = data(53); Cepmin = data(54); epmax = data(55);  Cepmax = data(56);
    ebar_p = data(57); Cebar_p = data(58); Ep = data(59); CEp = data(60); Epo = data(61);
    CEpo = data(62); W = data(63); CW = data(64); delta = data(65); Cdelta = data(66);
    delta_p = data(67); Cdelta_p = data(68);
    delta_y = data(69); Cdelta_y = data(70); delta_in = data(71); 
    Cdelta_in = data(72); deltap_in = data(73); Cdeltap_in = data(74);
    h = data(75); Ch = data(76); cbs = data(77); Ccbs = data(78);
    elbf = data(79); Celbf = data(80); slbf = data(81); Cslbf = data(82);
    Trule = data(83); Crule = data(84); Ttangent = data(85); Ctangent = data(86);
    Tstress = data(87); Cstress = data(88);
    Tstrain = data(89); Cstrain = data(90); strs = data(91); strn = data(92); strain_inc = data(93);
    repeat = data(94); eunld = data(95); funld = data(96); Ceunld = data(97); Cfunld = data(98); 
    fres_n = data(99); eresp_n = data(100); epo = data(101);
    
    //Set the trial state variables
    revertToLastCommit();
  }
  return res;
}

void 
RCFT_stlMaterial::Print(OPS_Stream &s, int flag){
  s<<"RCFT_stl, tag:"<<this->getTag()<<endln;
  s<<" fy:"<<fy<<endln;
  s<<" fu:"<<fu<<endln;
  s<<" Ee: "<<Ee<<endln;
  s<<" Ksft:"<<Ksft<<endln;
  return;
}
