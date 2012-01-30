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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/uniaxial/RCFT_concMaterial.cpp,v $

// Written: MHS
// Created: Aug 2001

// Modified by : Cenk Tort - tort0008@umn.edu
// University of Minnesota - Civil Engineering Department
// Date: Wed Jul 23 17:40:25 EDT 2003
// Description: This file contains the class implementation for
// cyclic uniaxial stress-strain relationship of concrete for
// RCFT members.

#include "RCFT_concMaterial.h"
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <fstream>
#include <iomanip>

//using std::ofstream;
//using std::ios;
//using std::endl;
using namespace std;

RCFT_concMaterial::RCFT_concMaterial( int tag, double fc, double depth, double thickness, double fy, double es )
  :UniaxialMaterial(tag, MAT_TAG_RCFT_conc), cp(0), Fc(fc), D(depth), t(thickness), Fy(fy), Es(es)
   , Cstrain(0.0), Cstress(0.0), Crule(0), Ctangent(0.0), Tstrain(0.0), Tstress(0.0), Ttangent(0.0) 
{
       	Fc_n = -Fc;
	Fc_p = 6 * sqrt(Fc*1000)/1000;
        nn = 0.8 + Fc * 1000 / 2500;
        np = 0.8 + Fc * 1000/ 2500;   
        Ec = 40.0 * sqrt (Fc*1000) + 1000.0; 
        //Ec = Ec * 0.77;
        //Ec = 2799.1;
        //Ec = 3154.096;
        //Ec = 4500.0;
        //Ec = 3530.31;
        //Ec = 1351.95 * sqrt(Fc);
        //Ec = 5322.70;
        //Ec = 2799.1; 
        //Ec = 5003.63;
        //Ec = 4931.1; 
        //Ec = 5076.14;
        //Ec = 4365.5;
        ec_n = ( Fc_n / Ec ) * nn / ( nn - 1.0 );
        ec_p = 1.23 * Fc_p / Ec;
	double R = (D/t) * sqrt(Fy/Es) * Fc/Fy;
	Ec_sft = - ( 332 * R - 9.60 ) * Fc;
        Fc_res_n = - ( Fc * 0.32 * pow(R,-0.5) );
        if((pow(R,-0.5)*0.32)>1.0){
	  Fc_res_n = - Fc * 0.85;
        }
        if(Fc_res_n < Fc ){
          Fc_res_n = - Fc * 0.85;
        }
        //Fc_res_n = 0.5 * Fc_res_n;
        ec_res_n = ( Fc_res_n - Fc_n + Ec_sft * ec_n ) / Ec_sft;
        eunld_n = 0.0;
        eunld_p = 0.0;
        eo_p = 0.0;
        eo_n = 0.0;
        crack = 0;
        r_n = 7;
        r_p = 4;
        xcr_n = 1.04;
        xcr_p = 5.3;

	funld_n = 0.0; funld_p = 0.0; ere_n = 0.0; 
	Ere_n = 0.0; Ttangent = Ec; Tstrain = 0.0; Tstress = 0.0; Esec_n = 0.0;
	Epl_n = 0.0; df_n = 0.0; de_n = 0.0; epl_n = 0.0; fnew_n = 0.0;
	Enew_n = 0.0; fre_n = 0.0; fnew_str_n = 0.0; Enew_str_n = 0.0;
	ere_str_n = 0.0; fre_str_n = 0.0; Ere_str_n = 0.0; deo_n = 0.0;
	Esec_p = 0.0; Epl_p = 0.0; df_p =0.0; de_p = 0.0; epl_p = 0.0;
        fnew_p = 0.0; Enew_p = 0.0; ere_p = 0.0; fre_p = 0.0; Ere_p = 0.0;
        xu_n = 0.0; xu_p = 0.0; deo_p = 0.0; ea = 0.0; fa =0.0; Ea = 0.0;
	ecr = 0.0; er = 0.0; fr = 0.0; 
	eb = 0.0; fb = 0.0; Eb = 0.0; fnew_str_p = 0.0; Enew_str_p = 0.0;
	ere_str_p = 0.0; fre_str_p = 0.0; Ere_str_p = 0.0; Cstrain = 0.0;
        Cstress = 0.0; Cstart_flag = 1; Tstart_flag = 1; direction = 1; cp = 1; history_n = 1;
        Ceunld_n = 0.0; Cfunld_n = 0.0;
        history_p = 1; crack = 0; Trule = 0; Frule = 0; Crule = 0;	

	Ccrack = 0; Cdirection = 1; Ceo_n = 0.0; Ceo_p = 0.0;
	CEsec_n = 0.0; CFrule = 0;
	CEpl_n = 0.0; Cdf_n = 0.0; Cde_n = 0.0; Cepl_n = 0.0;
	Cfnew_n = 0.0; CEnew_n = 0.0; Cere_n = 0.0; Cfre_n = 0.0;
	CEre_n = 0.0;  Cfnew_str_n = 0.0; CEnew_str_n = 0.0;
	Cere_str_n = 0.0; Cfre_str_n = 0.0; CEre_str_n = 0.0;
	Cdeo_n = 0.0; Cfunld_p = 0.0; Ceunld_p = 0.0; CEsec_p = 0.0;
	CEpl_p = 0.0; Cdf_p = 0.0; Cde_p = 0.0; Cepl_p = 0.0; Cfnew_p = 0.0;
	CEnew_p = 0.0; Cere_p = 0.0; Cfre_p = 0.0; CEre_p = 0.0; Cxu_n = 0.0;
	Cxu_p = 0.0; Cdeo_p = 0.0; Cea = 0.0; Cfa = 0.0; CEa = 0.0; Cer = 0.0; Cfr = 0.0;
	Er = CEr = 0.0;
	Ceb = 0.0; Cfb = 0.0; CEb = 0.0; Cfnew_str_p = 0.0; CEnew_str_p = 0.0;
	Cere_str_p = 0.0; Cfre_str_p = 0.0; CEre_str_p = 0.0;
	Cea1 = 0.0; ea1 = 0.0; Ceb1 = 0.0; eb1 = 0.0;
	Cfa1 = 0.0; fa1 = 0.0; Cfb1 = 0.0; fb1 = 0.0;
	CEa1 = 0.0; Ea1 = 0.0; CEb1 = 0.0; Eb1 = 0.0;
        Cer1 = 0.0; er1 = 0.0; 
	Cfr1 = 0.0; fr1 = 0.0;
	CEr1 = 0.0; Er1 = 0.0;
	CUrule = 0; Urule = 0;

	eunld_n7 = 0.0; funld_n7 = 0.0; Eunld_n7 = 0.0;
        eunld_p8 = 0.0; funld_p8 = 0.0; Eunld_p8 = 0.0;
        Ceunld_n7 = 0.0; Cfunld_n7 = 0.0; CEunld_n7 = 0.0;
        Ceunld_p8 = 0.0; Cfunld_p8 = 0.0; CEunld_p8 = 0.0;
	
	Chistory_n = 1; Chistory_p = 1;

        //cout<<Fc_n<<"  "<<ec_n<<"  "<<Fc_res_n<<"  "<<ec_res_n<<"  "<<1.0<<"  "<<Fc_p<<"  "<<Ec_sft<<endl;

}

RCFT_concMaterial::RCFT_concMaterial()
  :UniaxialMaterial(0,MAT_TAG_RCFT_conc),
   Cstart_flag(1), Tstart_flag(1), crack(0), Fc(0.0), D(0.0), t(0.0), Fy(0.0), Es(0.0), eo_n(0.0), eo_p(0.0) 
   , Cstrain(0.0), Cstress(0.0), Crule(0), Ctangent(0.0), Tstrain(0.0), Tstress(0.0), Ttangent(0.0)
{
       //DOES NOTHING
}

RCFT_concMaterial::~RCFT_concMaterial()
{

}

int 
RCFT_concMaterial::setTrialStrain(double strain, double strainRate)
{ 	
  double strain_incr, stress_incr;
  const double SMALL = 0.00000000000000000001;
  ofstream concrete;
  concrete.open("concrete.dat",ios::app);
  //ofstream conc;
  //conc.open("conc.dat",ios::app);

  Tstrain = Cstrain + strain;
  //Tstrain = strain;
  strain_incr = Tstrain - Cstrain;
  if( cp == 1 ){
     backtocommitStatevar();	   
  }
  if( fabs(Tstrain) < SMALL )
  {
     Tstress = 0.00000000000000001;
     Ttangent = Ec;
  }
  // DETERMINE THE DIRECTION OF LOADING
  if( ( Cstart_flag == 1 ) && ( Tstrain < -SMALL ) )
  {
     Trule = 1;
     Tstart_flag = 0;
     direction = 1;
     history_n = 0;
     cp = 0;
  }
  if( ( Cstart_flag == 1 ) && ( Tstrain > SMALL ) )
  {
     Trule = 2;
     Tstart_flag = 0;
     direction = 2;
     history_p = 0;
     cp = 0;
  }  
  if( ( Crule == 7 ) && ( strain_incr < -SMALL ) && ( Tstrain < ere_str_n ) && ( Frule == 3 ) &&( cp == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Crule == 7 ) && ( strain_incr < -SMALL ) && ( Tstrain < ere_n ) && ( cp == 1 ) && ( history_n == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Crule == 3 ) && ( strain_incr < -SMALL ) && ( Tstrain < eunld_n ) && ( cp == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Crule == 7 ) && ( strain_incr < -SMALL ) && ( Tstrain < ere_n ) && ( Frule == 10 ) && ( cp == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Crule == 7 ) && ( strain_incr < -SMALL ) && ( Tstrain < ere_n ) && ( Frule == 13 ) && ( cp == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Crule == 7 ) && ( strain_incr < -SMALL ) && ( Tstrain < ere_n ) && ( Frule == 11 ) && ( cp == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Crule == 7 ) && ( strain_incr < -SMALL ) && ( Tstrain < ere_n ) && ( Frule == 5 ) && ( cp == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Crule == 7 ) && ( strain_incr < -SMALL ) && ( Tstrain < ere_n ) && ( Frule == 55 ) && ( cp == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Crule == 11 ) && ( strain_incr < -SMALL ) && ( Tstrain < ere_n ) && ( cp == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Crule == 55 ) && ( strain_incr < -SMALL ) && ( Tstrain < ere_n ) && ( cp == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Crule == 5 ) && ( strain_incr < -SMALL ) && ( Tstrain < ere_n ) && ( cp == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Crule == 10 ) && ( strain_incr < -SMALL ) && ( Tstrain < 0.0 ) && ( cp == 1 ) && ( history_n == 1 ) )
  {
     Trule = 1;
     cp = 0;
     history_n = 0;
  }
  if( ( Crule == 4 ) && ( strain_incr < -SMALL ) && ( Tstrain < 0.0 ) && ( cp == 1 ) && ( history_n == 1 ) )
  {
     Trule = 1;
     cp = 0;
     history_n = 0;
  }
  if( ( Crule == 4 ) && ( strain_incr < -SMALL ) && ( Tstrain < ere_n ) && ( cp == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Crule == 10 ) && ( strain_incr < -SMALL ) && ( Tstrain < ere_n ) && ( cp == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Crule == 13 ) && ( strain_incr < -SMALL ) && ( Tstrain < 0.0 ) && ( cp == 1 ) && ( history_n == 1 ) )
  {
     Trule = 1;
     cp = 0;
     history_n = 0;
  }
  if( ( Crule == 13 ) && ( strain_incr < -SMALL ) && ( Tstrain < ere_n ) && ( cp == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Crule == 9 ) && ( strain_incr < -SMALL ) && ( Tstrain < ere_n ) && ( cp == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Crule == 1 ) && ( strain_incr < -SMALL ) && ( cp == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Crule == 2 ) && ( strain_incr < -SMALL ) && ( Tstrain < 0.0 ) && ( cp == 1 ) )
  {
     Trule = 1;
     cp = 0;
  }
  if( ( Trule == 1 ) && ( cp == 0 ) )
  {
     Tstress = stress_envlp_n ( Tstrain );
     Ttangent = tangent_envlp_n ( Tstrain );
  }
  /**********************************************RULE 2*****************************************************/
  if( ( Crule == 9 ) && ( strain_incr > SMALL ) && ( Tstrain > ere_p ) && ( cp == 1 ) )
  {
     Trule = 2;
     cp = 0;
  }
  if( ( Crule == 8 ) && ( strain_incr > SMALL ) && ( Tstrain > ere_p ) && ( Frule == 9 ) && ( cp == 1 ) )
  {
     Trule = 2;
     cp = 0;
  }
  if( ( Crule == 8 ) && ( strain_incr > SMALL ) && ( Tstrain > ere_p ) && ( Frule == 12 ) && ( cp == 1 ) )
  {
     Trule = 2;
     cp = 0;
  }
  if( ( Crule == 8 ) && ( strain_incr > SMALL ) && ( Tstrain > ere_str_p ) && ( Frule == 4 ) && ( cp == 1 ) )
  {
     Trule = 2;
     cp = 0;
  }
  if( ( Crule == 8 ) && ( strain_incr > SMALL ) && ( Tstrain > ere_p ) && ( Frule == 65 ) && ( cp == 1 ) )
  {
     Trule = 2;
     cp = 0;
  }
  if( ( Crule == 8 ) && ( strain_incr > SMALL ) && ( Tstrain > ere_p ) && ( Frule == 6 ) && ( cp == 1 ) )
  {
     Trule = 2;
     cp = 0;
  }
  if( ( Crule == 4 ) && ( strain_incr > SMALL ) && ( Tstrain > ere_p ) && ( cp == 1 ) )
  {
     Trule = 2;
     cp = 0;
  }
  if( ( Crule == 8 ) && ( strain_incr > SMALL ) && ( Tstrain > ere_str_p ) && ( Frule == 4 ) && ( cp == 1 ) )
  {
     Trule = 2;
     cp = 0;
  }
  if( ( Crule == 2 ) && ( strain_incr > SMALL ) && ( cp == 1 ) && ( crack == 0 ) )
  {
     Trule = 2;
     cp = 0;
  }
  if( ( Crule == 12 ) && ( strain_incr > SMALL ) && ( Tstrain > ere_p ) && ( cp == 1 ) )
  {
     Trule = 2;
     cp = 0;
  }
  if( ( Crule == 10 ) && ( strain_incr > SMALL ) && ( Tstrain > eunld_p ) && ( cp == 1 ) )
  {
     Trule = 2;
     cp = 0;
  }
  if( ( Crule == 65 ) && ( strain_incr > SMALL ) && ( Tstrain > ere_p ) && ( cp == 1 ) )
  {
     Trule = 2;
     cp = 0;
  }
  if( ( Crule == 6 ) && ( strain_incr > SMALL ) && ( Tstrain > ere_p ) && ( cp == 1 ) )
  {
     Trule = 2;
     cp = 0;
  }
  if( ( Crule == 1 ) && ( strain_incr > SMALL ) && ( Tstrain > 0.0 ) && ( cp == 1 ) )
  {
     Trule = 2;
     cp = 0;
  }
  if( ( Trule == 2 ) && ( cp == 0 ) )
  {
     if( history_n == 1 )
     {
	 Esec_p = Ec * ( fabs ( funld_p / ( Ec * ec_p ) ) + 0.67 ) / ( fabs( ( eunld_p - eo_p ) / ec_p ) + 0.67 );
	 epl_p = eunld_p - funld_p/Esec_p;
     }
     Tstress = stress_envlp_p ( Tstrain - eo_p );
     Ttangent = tangent_envlp_p ( Tstrain - eo_p );
  }
  /********************************************RULE 3*****************************************************/
  if( ( Crule == 1  ) && ( strain_incr > SMALL ) && ( cp == 1 ) )
  {
     Trule = 3;
     Frule = 1;
     funld_n = Cstress;
     eunld_n = Cstrain;
     if( Cstrain > ec_res_n ){
     Esec_n = Ec * ( fabs ( funld_n / ( Ec * ec_n ) ) + 0.57 ) / ( fabs( eunld_n / ec_n ) + 0.57 );
     Epl_n = 0.1 * Ec * pow( 2.71828 , -2 * fabs ( eunld_n / ec_n ) );
     df_n = 0.09 * funld_n * sqrt ( fabs( eunld_n / ec_n ) );
     de_n = ( eunld_n ) / ( 1.15 + 2.75 * fabs ( eunld_n / ec_n ) );
     epl_n = eunld_n - funld_n/Esec_n;
     fnew_n = funld_n - df_n;
     Enew_n = fnew_n / ( eunld_n - epl_n );
     ere_n = eunld_n + de_n;
     fre_n = stress_envlp_n ( ere_n );
     Ere_n = tangent_envlp_n ( ere_n );
     }
     if( Cstrain < ec_res_n ){
     Esec_n = Ec * ( fabs ( Fc_res_n / ( Ec * ec_n ) ) + 0.57 ) / ( fabs( eunld_n / ec_n ) + 0.57 );
     Epl_n = 0.1 * Ec * pow( 2.71828 , -2 * fabs ( eunld_n / ec_n ) );
     df_n = 0.09 * Fc_res_n * sqrt ( fabs( ec_res_n / ec_n ) );
     de_n = ( eunld_n ) / ( 1.15 + 2.75 * fabs ( eunld_n / ec_n ) );
     epl_n = eunld_n - Fc_res_n/Esec_n;
     fnew_n = Fc_res_n - df_n;
     Enew_n = fnew_n / ( eunld_n - epl_n );
     ere_n = eunld_n + de_n;
     fre_n = stress_envlp_n ( ere_n );
     Ere_n = tangent_envlp_n ( ere_n );
     }
     er = Cstrain;
     fr = Cstress;
     cp = 0;
  }
  if( ( Crule == 7 ) && ( strain_incr > SMALL ) && ( Cstrain > eunld_n )&& ( Tstrain < epl_n ) && ( cp == 1 ) )
  {
     Trule = 3;
     Frule = 7;
     ea = Cstrain;
     fa = Cstress;
     Ea = Ctangent;
     cp = 0;
  }
  if( ( Crule == 3 ) && ( strain_incr > SMALL ) && ( Tstrain < epl_n ) && ( cp == 1 ) )
  {
     Trule = 3;
     cp = 0;
  }
  if( ( Trule == 3 ) && ( Frule == 1 ) && ( cp == 0 ) )
  {
     Tstress = stress_tran(Tstrain, eunld_n, funld_n, Ec, epl_n, 0.0, Epl_n);
     Ttangent = tangent_tran(Tstrain, eunld_n, funld_n, Ec, epl_n, 0.0, Epl_n);
  }
  if( ( Trule == 3 ) && ( Frule == 7 ) && ( Tstrain < er ) && ( Tstrain > ea ) && ( cp == 0 ) )
  {
     Tstress = stress_tran(Tstrain, ea, fa, Ec, er, fr, Er);
     Ttangent = tangent_tran(Tstrain, ea, fa, Ec, er, fr, Er);
  }
  if( ( Trule == 3 ) && ( Frule == 7 ) && ( Tstrain > er ) && ( cp == 0 ) )
  {
     Tstress = stress_tran(Tstrain, eunld_n, funld_n, Ec, epl_n, 0.0, Epl_n);
     Ttangent = tangent_tran(Tstrain, eunld_n, funld_n, Ec, epl_n, 0.0, Epl_n);
  }
  /********************************************RULE 5*****************************************************/
  if( ( Crule == 7 ) && ( Frule == 3 ) && ( strain_incr > SMALL ) && ( Cstrain < eunld_n ) && ( Cstrain > ere_str_n ) && ( cp == 1 ) )
  {
     Trule = 5;
     Frule = 7;
     eunld_n7 = Cstrain;
     funld_n7 = Cstress;
     Eunld_n7 = Ctangent;
     er = Cstrain;
     fr = Cstress;
     Er = Ctangent;
     cp = 0;
  }
  if( ( Crule == 7 ) && ( Frule == 7 ) && ( strain_incr > SMALL ) && ( Cstrain < eunld_n ) && ( Cstrain > ere_str_n ) && ( cp == 1 ) )
  {
     Trule = 5;
     Frule = 7;
     eunld_n7 = Cstrain;
     funld_n7 = Cstress;
     Eunld_n7 = Ctangent;
     er = Cstrain;
     fr = Cstress;
     Er = Ctangent;
     cp = 0;
  }
  if( ( Crule == 7 ) && ( Frule == 5 ) && ( strain_incr > SMALL ) && ( Cstrain < eunld_n ) && ( Cstrain > ere_n ) && ( cp == 1 ) )
  {
     Trule = 5;
     Frule = 7;
     eunld_n7 = Cstrain;
     funld_n7 = Cstress;
     Eunld_n7 = Ctangent;
     er = Cstrain;
     fr = Cstress;
     Er = Ctangent;
     cp = 0;
  }
  if( ( Crule == 7 ) && ( Frule == 10 ) && ( strain_incr > SMALL ) && ( Cstrain < eunld_n ) && ( Cstrain > ere_n ) && ( cp == 1 ) )
  {
     Trule = 5;
     Frule = 7;
     eunld_n7 = Cstrain;
     funld_n7 = Cstress;
     Eunld_n7 = Ctangent;
     er = Cstrain;
     fr = Cstress;
     Er = Ctangent;
     cp = 0;
  }
  if( ( Crule == 7 ) && ( Frule == 11 ) && ( strain_incr > SMALL ) && ( Cstrain < eunld_n ) && ( Cstrain > ere_n ) && ( cp == 1 ) )
  {
     Trule = 5;
     Frule = 7;
     eunld_n7 = Cstrain;
     funld_n7 = Cstress;
     Eunld_n7 = Ctangent;
     er = Cstrain;
     fr = Cstress;
     Er = Ctangent;
     cp = 0;
  }
  if( ( Crule == 7 ) && ( Frule == 13 ) && ( strain_incr > SMALL ) && ( Cstrain < eunld_n ) && ( Cstrain > ere_n ) && ( cp == 1 ) )
  {
     Trule = 5;
     Frule = 7;
     eunld_n7 = Cstrain;
     funld_n7 = Cstress;
     Eunld_n7 = Ctangent;
     er = Cstrain;
     fr = Cstress;
     Er = Ctangent;
     cp = 0;
  }
  if( ( Crule == 7 ) && ( Frule == 55 ) && ( strain_incr > SMALL ) && ( Cstrain < eunld_n ) && ( Cstrain > ere_n ) && ( cp == 1 ) )
  {
     Trule = 5;
     Frule = 7;
     eunld_n7 = Cstrain;
     funld_n7 = Cstress;
     Eunld_n7 = Ctangent;
     er = Cstrain;
     fr = Cstress;
     Er = Ctangent;
     cp = 0;
  }
  if( ( Crule == 55 ) && ( strain_incr > SMALL ) && ( Cstrain > ere_n ) && ( cp == 1 ) )
  {
     Trule = 5;
     Frule = 55;
     ea = Cstrain;
     fa = Cstress;
     Ea = Ctangent;
     cp = 0;
  }	  
  if( ( Crule == 5 ) && ( strain_incr > SMALL ) && ( Tstrain < epl_n ) && ( cp == 1 ) )
  {
    Trule = 5;
    cp = 0;
  }	  
  if( ( Trule == 5 ) && ( Frule == 55 ) && ( Tstrain < er ) && ( Tstrain > ea ) && ( cp == 0 ) )
  {
     Tstress = stress_tran(Tstrain, ea, fa, Ec, er, fr, Er);
     Ttangent = tangent_tran(Tstrain, ea, fa, Ec, er, fr, Er);
  }
  if( ( Trule == 5 ) && ( Frule == 55 ) && ( Tstrain > er ) && ( cp == 0 ) )
  {
     Tstress = stress_tran(Tstrain, eunld_n7, funld_n7, Ec, epl_n, 0.0, Epl_n);
     Ttangent = tangent_tran(Tstrain, eunld_n7, funld_n7, Ec, epl_n, 0.0, Epl_n);
  }
  if( ( Trule == 5 ) && ( Frule == 7 ) && ( cp == 0 ) )
  {
     Tstress = stress_tran(Tstrain, eunld_n7, funld_n7, Ec, epl_n, 0.0, Epl_n);
     Ttangent = tangent_tran(Tstrain, eunld_n7, funld_n7, Ec, epl_n, 0.0, Epl_n);
  }
  /********************************************RULE 55*****************************************************/
  if( ( Crule == 5 ) && ( strain_incr < -SMALL ) && ( Tstrain < epl_n ) && ( Tstrain > eunld_n7 ) && ( cp == 1 ) )
  {
     Trule = 55;
     if( Cstrain > er ){
     er = Cstrain;
     fr = Cstress;
     Er = Ctangent;
     }
     er1 = Cstrain;
     fr1 = Cstress;
     Er1 = Ctangent;
     cp = 0;
  } 
  if( ( Crule == 55 ) && ( strain_incr < - SMALL ) && ( Tstrain > eunld_n7 ) && ( cp == 1 ) )
  {
     Trule = 55;
     cp = 0;
  }	  
  if( ( Trule == 55 ) && ( cp == 0 ) )
  { 
     Tstress = stress_tran(Tstrain, er1, fr1, Ec, eunld_n7, funld_n7, Eunld_n7);
     Ttangent = tangent_tran(Tstrain, er1, fr1, Ec, eunld_n7, funld_n7, Eunld_n7);      
  }	  
  /*******************************************RULE 7*******************************************************/
  if( ( Crule == 3  ) && ( strain_incr < -SMALL ) && ( cp == 1 ) )
  {
     Trule = 7;
     Frule = 3;
     if( (Cstrain > er ) | ( Frule == 1 ) ){
     er = Cstrain;
     fr = Cstress;
     Er = Ctangent; 
     }
     er1 = Cstrain;
     fr1 = Cstress;
     Er1 = Ctangent;
     fnew_str_n = funld_n - df_n * ( eunld_n - Cstrain ) / ( eunld_n - epl_n );
     Enew_str_n = ( fnew_str_n - Cstress ) / ( eunld_n - Cstrain );
     ere_str_n = eunld_n + de_n * ( eunld_n - Cstrain ) / ( eunld_n - epl_n );
     fre_str_n = stress_envlp_n ( ere_str_n );
     Ere_str_n = tangent_envlp_n (ere_str_n );
     cp = 0;
  }
  if( ( Crule == 10 ) && ( strain_incr < -SMALL ) && ( Tstrain < eunld_n ) && ( Tstrain > ere_n )&& ( cp == 1 ) && ( history_n == 0 ) )
  {
     Trule = 7;
     Frule = 10;
     cp = 0;
  }
  if( ( Crule == 13 ) && ( strain_incr < -SMALL ) && ( Tstrain < eunld_n ) && ( Tstrain > ere_n ) && ( cp == 1 ) && ( history_n == 0 ) )
  {
     Trule = 7;
     Frule = 13;
     cp = 0;
  }
  if( ( Crule == 11 ) && ( strain_incr < -SMALL ) && ( Tstrain < eunld_n ) && ( cp == 1 ) )
  {
     Trule = 7;
     Frule = 11;
     cp = 0;
  }
  if( ( Crule == 5 ) && ( strain_incr < -SMALL ) && ( Tstrain < eunld_n7 ) && ( cp == 1 ) )
  {
     Trule = 7;
     Frule = 5;
     cp = 0;
  }
  if( ( Crule == 55 ) && ( strain_incr < -SMALL ) && ( Tstrain < eunld_n7 ) && ( cp == 1 ) )
  {
     Trule = 7;
     Frule = 55;
     cp = 0;
  }
  if( ( Crule == 7 ) && ( Frule == 3 ) && ( strain_incr < -SMALL ) && ( Tstrain > ere_str_n ) && ( cp == 1 ) )
  {
     Trule = 7;
     cp = 0;
  }
  if( ( Crule == 7 ) && ( ( Frule == 10 ) | ( Frule == 13 ) | ( Frule == 11 ) | ( Frule == 5 ) | ( Frule == 55 ) ) && ( strain_incr < -SMALL ) && ( Tstrain > ere_n ) && ( cp == 1 ) )
  {
     Trule = 7;
     cp = 0;
  }
  if ( ( Trule == 7 ) && ( Frule == 10 ) && ( cp == 0 ))
  {
     Tstress = stress_tran( Tstrain, eunld_n, fnew_n, Enew_n, ere_n, fre_n, Ere_n );
     Ttangent = tangent_tran( Tstrain, eunld_n, fnew_n, Enew_n, ere_n, fre_n, Ere_n );
  }
  if ( ( Trule == 7 ) && ( Frule == 13 ) && ( cp == 0 ))
  {
     Tstress = stress_tran( Tstrain, eunld_n, fnew_n, Enew_n, ere_n, fre_n, Ere_n );
     Ttangent = tangent_tran( Tstrain, eunld_n, fnew_n, Enew_n, ere_n, fre_n, Ere_n );
  }
  if( ( Trule == 7 ) && ( Frule == 3 ) && ( Tstrain > eunld_n ) && ( cp == 0 ))
  {
     Tstress = stress_tran ( Tstrain, er1, fr1, Ec, eunld_n, fnew_str_n, Enew_str_n );
     Ttangent = tangent_tran ( Tstrain, er1, fr1, Ec, eunld_n, fnew_str_n, Enew_str_n );
  }
  if( ( Trule == 7 ) && ( Frule == 3 ) && ( Tstrain < eunld_n ) && ( Tstrain > ere_str_n ) && ( cp == 0 ))
  {
     Tstress = stress_tran( Tstrain, eunld_n, fnew_str_n, Enew_str_n, ere_str_n, fre_str_n, Ere_str_n );
     Ttangent = tangent_tran( Tstrain, eunld_n, fnew_str_n, Enew_str_n, ere_str_n, fre_str_n, Ere_str_n );  
  }
  if( ( Trule == 7 ) && ( Frule == 11 ) && ( Tstrain < eunld_n ) && ( Tstrain > ere_n ) && ( cp == 0 ))
  {
     Tstress = stress_tran( Tstrain, eunld_n, fnew_n, Enew_n, ere_n, fre_n, Ere_n );
     Ttangent = tangent_tran( Tstrain, eunld_n, fnew_n, Enew_n, ere_n, fre_n, Ere_n );
  }
  if( ( Trule == 7 ) && ( Frule == 5 ) && ( Tstrain < eunld_n ) && ( Tstrain > ere_n ) && ( cp == 0 ))
  {
     Tstress = stress_tran( Tstrain, eunld_n7, funld_n7, Eunld_n7, ere_n, fre_n, Ere_n );
     Ttangent = tangent_tran( Tstrain, eunld_n7, funld_n7, Eunld_n7, ere_n, fre_n, Ere_n );
  }
  if( ( Trule == 7 ) && ( Frule == 55 ) && ( Tstrain < eunld_n ) && ( Tstrain > ere_n ) && ( cp == 0 ))
  {
     Tstress = stress_tran( Tstrain, eunld_n7, funld_n7, Eunld_n7, ere_n, fre_n, Ere_n );
     Ttangent = tangent_tran( Tstrain, eunld_n7, funld_n7, Eunld_n7, ere_n, fre_n, Ere_n );
  }
  if( ( Trule == 7 ) && ( Frule == 3 ) && ( Tstrain < ere_str_n ) && ( cp == 0 ))
  {
     Tstress = stress_envlp_n( Tstrain );
     Ttangent = tangent_envlp_n( Tstrain );
  }
  if( ( Trule == 7 ) && ( Frule == 3 ) && ( Tstrain < ere_n ) && ( cp == 0 ))
  {
     Tstress = stress_envlp_n( Tstrain );
     Ttangent = tangent_envlp_n( Tstrain );
  }
  if( ( Trule == 7 ) && ( Frule == 11 ) && ( Tstrain < ere_n ) && ( cp == 0 ))
  {
     Tstress = stress_envlp_n( Tstrain );
     Ttangent = tangent_envlp_n( Tstrain );
  }
  /*******************************************RULE 9********************************************************/
  if( ( Crule == 3 ) && ( strain_incr > SMALL ) && ( Tstrain > epl_n ) && ( cp == 1 ) && ( crack == 0 ) )
  {
     Trule = 9;
     Frule = 3;
     cp = 0;
     xu_n = fabs( eunld_n  / ec_n );
     xu_p = fabs( ( eunld_p - eo_p ) / ec_p );
     if( xu_p <= xu_n )
     {
        xu_p = xu_n;
        eo_p = 0;
        eunld_p = xu_p * ec_p;
        funld_p = stress_envlp_p ( eunld_p - eo_p );
     }
     Esec_p = Ec * ( fabs ( funld_p / ( Ec * ec_p ) ) + 0.67 ) / ( fabs( (eunld_p - eo_p) / ec_p ) + 0.67 );
     deo_p = 2 * funld_p / ( Esec_p + Epl_n);
     eo_p = epl_n + deo_p - xu_p * ec_p;
     eunld_p = xu_p * ec_p + eo_p;
     funld_p = stress_envlp_p ( eunld_p - eo_p );
     Epl_p = Ec / ( pow( fabs( ( eunld_p - eo_p ) / ( ec_p ) ), 1.1 ) + 1 );
     df_p = 0.15 * funld_p;
     de_p = 0.22 * eunld_p;
     epl_p = eunld_p - funld_p/Esec_p;
     fnew_p = funld_p - df_p;
     Enew_p = fnew_p / ( eunld_p - epl_p );
     ere_p = eunld_p + fabs( de_p );
     fre_p = stress_envlp_p ( fabs( ere_p - eo_p ) );
     Ere_p = tangent_envlp_p ( fabs( ere_p - eo_p ) );
  }
  if( ( Crule == 7 ) && ( strain_incr > SMALL ) && ( Tstrain > epl_n ) && ( cp == 1 ) && ( crack == 0 ) )
  {
     Trule = 9;
     Frule = 7;
     cp = 0;
     xu_n = fabs( eunld_n  / ec_n );
     xu_p = fabs( ( eunld_p - eo_p ) / ec_p );
     if( xu_p <= xu_n )
     {
       xu_p = xu_n;
       eo_p = 0;
       eunld_p = xu_p * ec_p;
       funld_p = stress_envlp_p ( eunld_p - eo_p );
     }
     Esec_p = Ec * ( fabs ( funld_p / ( Ec * ec_p ) ) + 0.67 ) / ( fabs( (eunld_p - eo_p) / ec_p ) + 0.67 );
     deo_p = 2 * funld_p / ( Esec_p + Epl_n);
     eo_p = epl_n + deo_p - xu_p * ec_p;
     eunld_p = xu_p * ec_p + eo_p;
     funld_p = stress_envlp_p ( eunld_p - eo_p );
     Epl_p = Ec / ( pow( fabs( ( eunld_p - eo_p ) / ( ec_p ) ), 1.1 ) + 1 );
     df_p = 0.15 * funld_p;
     de_p = 0.22 * eunld_p;
     epl_p = eunld_p - funld_p/Esec_p;
     fnew_p = funld_p - df_p;
     Enew_p = fnew_p / ( eunld_p - epl_p );
     ere_p = eunld_p + fabs( de_p );
     fre_p = stress_envlp_p ( fabs( ere_p - eo_p ) );
     Ere_p = tangent_envlp_p ( fabs( ere_p - eo_p ) );
  }
  if( ( Crule == 5 ) && ( strain_incr > SMALL ) && ( Tstrain > epl_n ) && ( cp == 1 ) && ( crack == 0 ) )
  {
     Trule = 9;
     Frule = 5;
     cp = 0;
     xu_n = fabs( eunld_n  / ec_n );
     xu_p = fabs( ( eunld_p - eo_p ) / ec_p );
     if( xu_p <= xu_n )
     {
       xu_p = xu_n;
       eo_p = 0;
       eunld_p = xu_p * ec_p;
       funld_p = stress_envlp_p ( eunld_p - eo_p );
     }
     Esec_p = Ec * ( fabs ( funld_p / ( Ec * ec_p ) ) + 0.67 ) / ( fabs( (eunld_p - eo_p) / ec_p ) + 0.67 );
     deo_p = 2 * funld_p / ( Esec_p + Epl_n);
     eo_p = epl_n + deo_p - xu_p * ec_p;
     eunld_p = xu_p * ec_p + eo_p;
     funld_p = stress_envlp_p ( eunld_p - eo_p );
     Epl_p = Ec / ( pow( fabs( ( eunld_p - eo_p ) / ( ec_p ) ), 1.1 ) + 1 );
     df_p = 0.15 * funld_p;
     de_p = 0.22 * eunld_p;
     epl_p = eunld_p - funld_p/Esec_p;
     fnew_p = funld_p - df_p;
     Enew_p = fnew_p / ( eunld_p - epl_p );
     ere_p = eunld_p + fabs( de_p );
     fre_p = stress_envlp_p ( fabs( ere_p - eo_p ) );
     Ere_p = tangent_envlp_p ( fabs( ere_p - eo_p ) );
  }
  if( ( Crule == 12 ) && ( Urule == 9 ) && ( strain_incr > SMALL ) && ( Tstrain > ea ) && ( cp == 1 ) )
  {
     Trule = 9;
     cp = 0;
     xu_n = fabs( eunld_n  / ec_n );
     xu_p = fabs( ( eunld_p - eo_p ) / ec_p );
     if( xu_p <= xu_n )
     {
        xu_p = xu_n;
        eo_p = 0;
        eunld_p = xu_p * ec_p;
        funld_p = stress_envlp_p ( eunld_p - eo_p );
     }
     Esec_p = Ec * ( fabs ( funld_p / ( Ec * ec_p ) ) + 0.67 ) / ( fabs( (eunld_p - eo_p) / ec_p ) + 0.67 );
     deo_p = 2 * funld_p / ( Esec_p + Epl_n);
     eo_p = epl_n + deo_p - xu_p * ec_p;
     eunld_p = xu_p * ec_p + eo_p;
     funld_p = stress_envlp_p ( eunld_p - eo_p );
     Epl_p = Ec / ( pow( fabs( ( eunld_p - eo_p ) / ( ec_p ) ), 1.1 ) + 1 );
     df_p = 0.15 * funld_p;
     de_p = 0.22 * eunld_p;
     epl_p = eunld_p - funld_p/Esec_p;
     fnew_p = funld_p - df_p;
     Enew_p = fnew_p / ( eunld_p - epl_p );
     ere_p = eunld_p + fabs( de_p );
     fre_p = stress_envlp_p ( fabs( ere_p - eo_p ));
     Ere_p = tangent_envlp_p ( fabs( ere_p - eo_p ));
   }
   if( ( Crule == 9 ) && ( strain_incr > SMALL ) && ( Tstrain < eunld_p ) && ( cp == 1 ) )
   {
     Trule = 9;
     cp = 0;
   }
   if( ( Crule == 12 ) && ( strain_incr < -SMALL ) && ( Tstrain > ea ) && ( cp == 1 ) )
   {
     Trule = 9;
     cp = 0;
   }
   if( ( Crule == 11 ) && ( Frule == 9 ) && ( strain_incr > SMALL ) && ( Tstrain > ea ) && ( cp == 1 ) )
   {
     Trule = 9;
     cp = 0;
   }
   if( ( Trule == 9 ) && ( cp == 0 ) )
   {
     //Tstress = stress_tran ( Tstrain, epl_n, 0.0, Epl_n, eunld_p, fnew_p, Enew_p );
     //Ttangent = tangent_tran ( Tstrain, epl_n, 0.0, Epl_n, eunld_p, fnew_p, Enew_p );
     Tstress = ( ( fnew_p ) / ( eunld_p - epl_n ) ) * ( Tstrain - epl_n );
     Ttangent = ( fnew_p ) / ( eunld_p - epl_n );
     if( Tstrain > eunld_p ){
       Tstress = stress_envlp_p ( Tstrain );
       Ttangent = tangent_envlp_p ( Tstrain );
     }	     
   }
   /*******************************************RULE 16**************************************************/
   if( ( Crule == 3 ) && ( strain_incr > SMALL ) && ( Tstrain > epl_n ) && ( crack == 1 ) && ( cp == 1 ) )
   {
     Trule = 16;
     cp = 0;
   }
   if( ( Crule == 5 ) && ( strain_incr > SMALL ) && ( Tstrain > epl_n ) && ( crack == 1 ) && ( cp == 1 ) )
   {
     Trule = 16;
     cp = 0;
   }
   if( ( Crule == 2 ) && ( strain_incr > SMALL ) && ( crack == 1 ) && ( cp == 1 ) )
   {
     Trule = 16;
     cp = 0;
   }
   if( ( Crule == 13 ) && ( strain_incr > SMALL ) && ( Tstrain > eb ) && ( crack == 1 ) && ( cp == 1 ) )
   {
     Trule = 16;
     cp = 0;
   }	   
   if( ( Crule == 14 ) && ( strain_incr > SMALL ) && ( Tstrain > eb ) && ( crack == 1 ) && ( cp == 1 ) )
   {
     Trule = 16;
     cp = 0;
   }
   if( ( Crule == 15 ) && ( strain_incr > SMALL ) && ( Tstrain > eb ) && ( crack == 1 ) && ( cp == 1 ) )
   {
     Trule = 16;
     cp = 0;
   }
   if( ( Crule == 16 ) && ( strain_incr > SMALL ) && ( crack == 1 ) && ( cp == 1 ) )
   {
     Trule = 16;
     cp = 0;
   }
   if( ( Trule == 16 ) && ( cp == 0 ) )
   {
     Tstress = 0.0000001;
     Ttangent = 0.0000001;
   }
   /*******************************************RULE 13************************************************/
   if( ( Crule == 16 ) && ( strain_incr < -SMALL ) && ( crack == 1 ) && ( cp == 1 )  )
   {
     Trule = 13;
     ecr = Cstrain;
     cp = 0;
   }
   if( ( Crule == 14 ) && ( strain_incr < -SMALL ) && ( Tstrain < ea ) && ( crack == 1 ) && ( cp == 1 ) )
   {
     Trule = 13;
     cp = 0;
   }
   if( ( Crule == 15 ) && ( strain_incr < -SMALL ) && ( Tstrain < ea ) && ( crack == 1 ) && ( cp == 1 ) )
   {
     Trule = 13;
     cp = 0;
   }
   if( ( Crule == 13 ) && ( strain_incr < -SMALL ) && ( crack == 1 ) && ( cp == 1 ) )
   {
     Trule = 13;
     eb = Cstrain - Cstress / Esec_n;
     fb = 0;
     Eb = 0;
     cp = 0;
   }
   if( ( Trule == 13 ) && ( cp == 0 ) ) 
   {
     //if( history_n == 0 ){
     //Tstress = stress_tran ( Tstrain , ecr, 0.0, 0.0, eunld_n, fnew_n, Enew_n );
     //Ttangent = tangent_tran ( Tstrain , ecr, 0.0, 0.0, eunld_n, fnew_n, Enew_n );
     //}
     //else if( history_n == 1 ){
     //Tstress = 0.000001;
     //Ttangent = 0.000001; 
     //}
     if( fabs(fnew_n) > 0.0000000000001 ){
     Tstress = stress_tran (Tstrain, ecr, 0.0, 0.0, eunld_n, fnew_n, Enew_n);
     Ttangent = tangent_tran (Tstrain, ecr, 0.0, 0.0, eunld_n, fnew_n, Enew_n);
     }
     else{
     Tstress = 0.0;
     Ttangent = 0.0;
     }
     if( fabs(Tstress) < 1E-10){
        Tstress = 0.0000000000001;
     }
     if( fabs(Ttangent) < 1E-10){
        Ttangent = 0.0000000000001;
     }
   }
   /******************************************RULE 11************************************************/
   if( ( Crule == 9 ) && ( strain_incr < -SMALL ) && ( cp == 1 ) )
   {
     Trule = 11;
     //if( Frule == 3 ){
     //  eb = eunld_n;
     //  fb = fnew_n;
     //  Eb = Enew_n;
     //}
     //if( Frule == 11 ){
     //  eb = eunld_n7;
     //  fb = funld_n7;
     //  Eb = Eunld_n7;
     //}	
     if( fabs( eunld_n7 ) > fabs( eunld_n ) ){
         eb = eunld_n7;
	 fb = funld_n7;
	 Eb = Eunld_n7;
     }	  
     else{
         eb = eunld_n;
	 fb = fnew_n;
	 Eb = Enew_n;
     }	     
     Frule = 9;
     Urule = 9;
     ea = Cstrain;
     fa = Cstress;
     Ea = Ctangent;
     eb1 = ea - fa / Ec;
     fb1 = 0.0;
     Eb1 = ( fb - fb1 ) / ( eb - eb1 );    
     ea1 = Cstrain;
     fa1 = Cstress;
     Ea1 = Ec;
     cp = 0;
   }
   if ( ( Crule == 12 ) && ( strain_incr < -SMALL ) && ( Tstrain > eb ) && ( Tstrain < ea ) && ( cp == 1 ) )
   {
     Trule = 11;
     Frule = 12;
     if( Cstrain < eb1 ){
    	eb1 = Cstrain;
     	fb1 = Cstress;
     	Eb1 = Ec;
     }	     
     if( Cstress > 0.0 ){
	ea1 = Cstrain; 
        fa1 = Cstress;
        Ea1 = Ctangent;	
     	eb1 = Cstrain -Cstress/Ec;
     	fb1 = 0.0;
     	Eb1 = Ec;
     }
     if( Cstress < 0.0 ){
     	ea1 = Cstrain;
        fa1 = Cstress;
     	Ea1 = Ctangent;     
     }
     cp = 0;
   }
   if( ( Crule == 11 ) && ( strain_incr < -SMALL ) && ( Tstrain > eb ) && ( cp == 1 ) )
   {
     Trule = 11;
     cp = 0;
   }
   if ( ( Trule == 11 ) && ( Frule == 12 ) && ( Tstrain < ea1 ) && ( Tstrain > eb1 ) && ( cp == 0 ) )
   {
     Tstress = fa1 + Ec * ( Tstrain - ea1 );
     Ttangent = Ec;
   }	   
   if ( ( Trule == 11 ) && ( Frule == 12 ) && ( Tstrain < eb1 ) && ( cp == 0 ) )
   {
     Tstress = fb + ( ( fb - fb1 ) / ( eb - eb1 ) ) * ( Tstrain - eb );    
     Ttangent = ( fb - fb1 ) / ( eb - eb1 );
   }
   if ( ( Trule == 11 ) && ( Frule == 9 ) && ( Tstrain < ea ) && ( Tstrain > eb1 ) && ( cp == 0 ) )
   { 
     Tstress  = fa + Ec * ( Tstrain - ea );
     Ttangent = Ec;	     
   }	   
   if ( ( Trule == 11 ) && ( Frule == 9 ) && ( Tstrain < eb1 ) && ( Tstrain > eb ) && ( cp == 0 ) )
   {
     Tstress  = fb + ( ( fb - fb1 ) / ( eb - eb1 ) ) * ( Tstrain - eb );
     Ttangent = ( fb - fb1 ) / ( eb - eb1 );
   }
   /*****************************************RULE 12**************************************************/
   if ( ( Crule == 10 ) && ( strain_incr > SMALL ) && ( Tstrain < eunld_p ) && ( cp == 1 ) )
   {
      Trule = 12;
      if( fabs( eunld_p8 ) > fabs( eunld_p ) ){
         ea = eunld_p8;
	 fa = funld_p8;
	 Ea = Eunld_p8;
      }
      else{
	 ea = eunld_p;
	 fa = fnew_p;
	 Ea = Enew_p;
      }
      Frule = 10;
      Urule = 10;
      eb = Cstrain;
      fb = Cstress;
      Eb = Ctangent;
      ea1 = eb - fb / Ec;
      fa1 = 0.0;
      Ea1 = Ec; 
      eb1 = Cstrain;
      fb1 = Cstress;
      Eb1 = Ec; 
      cp = 0;
   }
   if( ( Crule == 11 ) && ( strain_incr > SMALL ) && ( Tstrain > eb ) && ( Tstrain < ea ) && ( cp == 1 ) )
   {
      Trule = 12;
      Frule = 11;
      if( Cstrain > ea1 ){
     	   ea1 = Cstrain;
           fa1 = Cstress;
           Ea1 = Ec;
      } 	      
      if( Cstress < 0.0 ){
	   eb1 = Cstrain;
	   fb1 = Cstress;
	   Eb1 = Ctangent;   
           ea1 = Cstrain - Cstress / Ec;
           fa1 = 0.0;
           Ea1 = Ec;
      }
      else if( Cstress > 0.0 ){
           eb1 = Cstrain;
           fb1 = Cstress;
           Eb1 = Ctangent;      
      }
      cp = 0;
   }
   if( ( Crule == 12 ) && ( strain_incr > SMALL ) && ( Tstrain < ea ) && ( cp == 1 ) )
   {
      Trule = 12;
      cp = 0;
   }
   if( ( Trule == 12 ) && ( Frule == 10 ) && ( Tstrain > eb ) && ( Tstrain < ea1 ) && ( cp == 0 ) )
   { 
      Tstress  = fb + Ec * ( Tstrain - eb );	   
      Ttangent = Ec; 		   
      if( history_n == 1 ){
         Tstress  = fa + ((fa - fa1)/(ea - ea1)) * (Tstrain - ea);
	 Ttangent = (fa - fa1)/(ea - ea1);   
      }
   }
   if( ( Trule == 12 ) && ( Frule == 10 ) && ( Tstrain > ea1 ) && ( Tstrain < ea ) && ( cp == 0 ) )
   {
      Tstress  = fa + ((fa - fa1)/(ea - ea1)) * (Tstrain - ea);
      Ttangent = (fa - fa1)/(ea - ea1);   
   }	   
   if( ( Trule == 12 ) && ( Frule == 11 ) && ( Tstrain > eb1 ) && ( Tstrain < ea1 ) && ( cp == 0 ) )
   {
      Tstress  = fb1 + Ec * (Tstrain - eb1);
      Ttangent = Ec;
      if( (history_n == 1) && ( Tstress > 0.000001 ) ) {
        Tstress  = fa + ((fa - fb1)/(ea - eb1)) * (Tstrain - ea);
        Ttangent = (fa - fb1)/(ea - eb1);   
      }
   }
   if( ( Trule == 12 ) && ( Frule == 11 ) && ( Tstrain > ea1 ) && ( Tstrain < ea ) && ( cp == 0 ) )
   {
      Tstress  = fa1 +  ( ( fa - fa1 ) / ( ea - ea1 ) ) * (Tstrain - ea1);
      Ttangent = ( fa - fa1 ) / ( ea - ea1 );
   }
   /*****************************************RULE 4****************************************************/
   if( ( Crule == 2 ) && ( strain_incr < -SMALL ) && ( cp == 1 ) && ( crack == 0 ) )
   {
     Trule = 4;
     Frule = 2;
     funld_p = Cstress;
     eunld_p = Cstrain;
     ere_p = 1.22 * eunld_p; 
     Esec_p = Ec * ( fabs ( funld_p / ( Ec * ec_p ) ) + 0.67 ) / ( fabs( ( eunld_p - eo_p ) / ec_p ) + 0.67 );
     epl_p = eunld_p - funld_p/Esec_p;
     Epl_p = Ec / ( pow( fabs( ( eunld_p - eo_p ) / ( ec_p ) ), 1.1 ) + 1 );
     df_p = 0.15 * funld_p;
     fnew_p = funld_p - df_p;
     Enew_p = fnew_p / ( eunld_p - epl_p );
     df_p = 0.15 * funld_p;
     de_p = fabs(0.22 * eunld_p);
     fre_p = stress_envlp_p ( fabs( ere_p - eo_p ) );
     Ere_p = tangent_envlp_p ( fabs( ere_p - eo_p ) );
     fnew_str_p = funld_p - df_p * ( eunld_p - Cstrain ) / ( eunld_p - epl_p );
     Enew_str_p = ( fnew_p - Cstress ) / ( eunld_p - Cstrain );
     ere_str_p = eunld_p + de_p * ( eunld_p - Cstrain ) / ( eunld_p - epl_p );
     fre_str_p = stress_envlp_p ( ere_str_p - eo_p );
     Ere_str_p = tangent_envlp_p ( ere_str_p - eo_p );
     er = Cstrain;
     fr = Cstress;
     Er = Ctangent;
     if( Tstrain < epl_p){
        Trule = 10;
     }	     
     cp = 0;
   }
   if( ( Crule == 8 ) && ( Frule == 4 ) && ( strain_incr < -SMALL ) && ( Tstrain < eunld_p ) && ( cp == 1 ) && ( crack == 0 ) )
   {
     Trule = 4;
     Frule = 8;
     ea = Cstrain;
     fa = Cstress;
     Ea = Ctangent;
     cp = 0;
   }
   if( ( Crule == 4 ) && ( strain_incr < -SMALL ) && ( Tstrain > epl_p ) && ( cp == 1 ) && ( crack == 0 ) )
   {
     Trule = 4;
     cp = 0;
   }
   if( ( Trule == 4 ) && ( Frule == 2 ) && ( cp == 0 ) )
   {
     Tstress = stress_tran (Tstrain, eunld_p, funld_p, Ec, epl_p , 0.0, Epl_p);
     Ttangent = tangent_tran (Tstrain, eunld_p, funld_p, Ec, epl_p, 0.0, Epl_p);
   }
   if( ( Trule == 4 ) && ( Frule == 8 ) && ( Tstrain < ea ) && ( Tstrain > er ) && ( cp == 0 ) )
   {
     Tstress = stress_tran(Tstrain, ea, fa, Ec, er, fr, Er);
     Ttangent = tangent_tran(Tstrain, ea, fa, Ec, er, fr, Er);
   }
   if( ( Trule == 4 ) && ( Frule == 8 ) && ( Tstrain < er ) && ( cp == 0 ) )
   {
     Tstress = stress_tran (Tstrain, eunld_p, funld_p, Ec, epl_p , 0.0, Epl_p);
     Ttangent = tangent_tran (Tstrain, eunld_p, funld_p, Ec, epl_p, 0.0, Epl_p);
   }
   /********************************************RULE 6*****************************************************/
   if( ( Crule == 8 ) && ( Frule == 4 ) && ( strain_incr < -SMALL ) && ( Cstrain > eunld_p ) && ( Cstrain < ere_str_p ) && ( cp == 1 ) )
   {
     Trule = 6;
     Frule = 8;
     eunld_p8 = Cstrain;
     funld_p8 = Cstress;
     Eunld_p8 = Ctangent;
     er = Cstrain;
     fr = Cstress;
     Er = Ctangent;
     cp = 0;
   }
   if( ( Crule == 8 ) && ( Frule == 9 ) && ( strain_incr < -SMALL ) && ( Cstrain > eunld_p ) && ( Cstrain < ere_p ) && ( cp == 1 ) )
   {
     Trule = 6;
     Frule = 8;
     eunld_p8 = Cstrain;
     funld_p8 = Cstress;
     Eunld_p8 = Ctangent;
     er = Cstrain;
     fr = Cstress;
     Er = Ctangent;
     cp = 0;
   }
   if( ( Crule == 8 ) && ( Frule == 12 ) && ( strain_incr < -SMALL ) && ( Cstrain > eunld_p ) && ( Cstrain < ere_p ) && ( cp == 1 ) )
   {
     Trule = 6;
     Frule = 8;
     eunld_p8 = Cstrain;
     funld_p8 = Cstress;
     Eunld_p8 = Ctangent;
     er = Cstrain;
     fr = Cstress;
     Er = Ctangent;
     cp = 0;
   }
   if( ( Crule == 8 ) && ( Frule == 65 ) && ( strain_incr < -SMALL ) && ( Cstrain > eunld_p ) && ( Cstrain < ere_p ) && ( cp == 1 ) )
   {
     Trule = 6;
     Frule = 8;
     eunld_p8 = Cstrain;
     funld_p8 = Cstress;
     Eunld_p8 = Ctangent;
     er = Cstrain;
     fr = Cstress;
     Er = Ctangent;
     cp = 0;
   }
   if( ( Crule == 8 ) && ( Frule == 6 ) && ( strain_incr < -SMALL ) && ( Cstrain > eunld_p ) && ( Cstrain < ere_p ) && ( cp == 1 ) ){
     Trule = 6;
     Frule = 8;
     eunld_p8 = Cstrain;
     funld_p8 = Cstress;
     Eunld_p8 = Ctangent;
     er = Cstrain;
     fr = Cstress;
     Er = Ctangent;
     cp = 0;
   }
   if( ( Crule == 65 ) && ( strain_incr < -SMALL ) && ( Tstrain < ere_p ) && ( cp == 1 ) )
   {
     Trule = 6;
     Frule = 65;
     ea = Cstrain;
     fa = Cstress;
     Ea = Ctangent;
     cp = 0;
   }
   if( ( Crule == 6 ) && ( strain_incr < -SMALL ) && ( Tstrain > epl_p ) && ( cp == 1 ) )
   {
     Trule = 6;
     cp = 0;
   }
   if( ( Trule == 6 ) && ( Frule == 65 ) && ( Tstrain  < ea ) && ( Tstrain > er ) && ( cp == 0 ) )
   {
     Tstress = stress_tran(Tstrain, ea, fa, Ec, er, fr, Er);
     Ttangent = tangent_tran(Tstrain, ea, fa, Ec, er, fr, Er);
   }
   if( ( Trule == 6 ) && ( Frule == 65 ) && ( Tstrain < er ) && ( cp == 0 ) )
   {
     Tstress = stress_tran(Tstrain, eunld_p8, funld_p8, Ec, epl_p, 0.0, Epl_p);
     Ttangent = tangent_tran(Tstrain, eunld_p8, funld_p8, Ec, epl_p, 0.0, Epl_p);
   }
   if( ( Trule == 6 ) && ( Frule == 8 ) && ( cp == 0 ) )
   {
     Tstress = stress_tran(Tstrain, eunld_p8, funld_p8, Ec, epl_p, 0.0, Epl_p);
     Ttangent = tangent_tran(Tstrain, eunld_p8, funld_p8, Ec, epl_p, 0.0, Epl_p);
   }
   /********************************************RULE 65*****************************************************/
   if( ( Crule == 6 ) && ( strain_incr > SMALL ) && ( Tstrain > epl_p ) && ( Tstrain < eunld_p8 ) && ( cp == 1 ) )
   {
     Trule = 65;
     er1 = Cstrain;
     fr1 = Cstress;
     Er1 = Ctangent;
     cp = 0;
   }
   if( ( Crule == 65 ) && ( strain_incr > SMALL ) && ( Tstrain < eunld_p8 ) && ( cp == 1 ) )
   {
     Trule = 65;
     cp = 0;
   }
   if( ( Trule == 65 ) && ( cp == 0 ) )
   {
     Tstress = stress_tran(Tstrain, er1, fr1, Ec, eunld_p8, funld_p8, Eunld_p8);
     Ttangent = tangent_tran(Tstrain, er1, fr1, Ec, eunld_p8, funld_p8, Eunld_p8);
   }
   /****************************************RULE 10****************************************************/
   if( ( Crule == 4 ) && ( strain_incr < -SMALL ) && ( Tstrain < epl_p ) && ( cp == 1 ) )
   {
     Trule = 10;
     Frule = 4;
     if( history_n == 1 )
     {
        eunld_n = 0.0;
        fnew_n = 0.0;
        Enew_n = 0.0;
     }
     cp = 0;
   }
   if( ( Crule == 6 ) && ( strain_incr < -SMALL ) && ( Tstrain < epl_p ) && ( cp == 1 ) )
   {
     Trule = 10;
     Frule = 6;
     if( history_n == 1 )
     {
        eunld_n = 0.0;
        fnew_n = 0.0;
        Enew_n = 0.0;
     }
     cp = 0;
   }
   if( ( Crule == 2 ) && ( strain_incr < -SMALL ) && ( Tstrain < epl_p ) && ( cp == 1 ) )
   {
     Trule = 10;
     cp = 0;
   }	   
   if( ( Crule == 8 ) && ( strain_incr < -SMALL ) && ( Tstrain < epl_p ) && ( cp == 1 ) )
   {
     Trule = 10;
     cp = 0;
   }
   if( ( Crule == 11 ) && ( Urule == 10 ) && ( Tstrain < eb ) && ( strain_incr < -SMALL ) && ( cp == 1 ) )
   {
     Trule = 10;
     if( history_n == 1 )
     {
        eunld_n = 0.0;
        fnew_n = 0.0;
        Enew_n = 0.0;
     }
     cp = 0;
   }
   if( ( Crule == 12 ) && ( Urule == 10 ) && ( Tstrain < eb ) && ( strain_incr < -SMALL ) && ( cp == 1 ) )
   {
     Trule = 10;
     if( history_n == 1 )
     {
        eunld_n = 0.0;
        fnew_n = 0.0;
        Enew_n = 0.0;
     }
     cp = 0;
   }
   if( ( Crule == 10 ) && ( strain_incr < -SMALL ) && ( Tstrain > eunld_n ) && ( cp == 1 ) )
   {
     Trule = 10;
     ere_str_n = eunld_n + de_n * ( eunld_n - Cstrain ) / ( eunld_n - epl_n );
     fre_str_n = stress_envlp_n ( ere_str_n );
     cp = 0;
   }
   if( ( Crule == 12 ) && ( Frule == 10 ) && ( strain_incr < -SMALL ) && ( Tstrain < eb ) && ( cp == 1 ) )
   {
     Trule = 10;
     cp = 0;
   }
   if( ( Trule == 10 ) && ( cp == 0 ) )
   {
     if( fabs(fnew_n) > 0.0000000000001 ){
     Tstress = stress_tran (Tstrain, epl_p, 0.0, Epl_p, eunld_n, fnew_n, Enew_n);
     Ttangent = tangent_tran (Tstrain, epl_p, 0.0, Epl_p, eunld_n, fnew_n, Enew_n);
     }
     else{
     Tstress = 0.0;
     Ttangent = 0.0;
     }
     if( fabs(Tstress) < 1E-10){
        Tstress = 0.0000000000001;
     }
     if( fabs(Ttangent) < 1E-10){
        Ttangent = 0.0000000000001;
     }
   }
   /****************************************RULE 8**************************************************/
   if( ( Crule == 4 ) && ( strain_incr > SMALL ) && ( Tstrain > epl_p ) && ( cp == 1 ) )
   {
     Trule = 8;
     Frule = 4;
     if( Cstrain < er ){
     er = Cstrain;
     fr = Cstress;
     Er = Ctangent;
     }
     er1 = Cstrain;
     fr1 = Cstress;
     Er1 = Ctangent;
     cp = 0;
   }
   if( ( Crule == 9 ) && ( strain_incr > SMALL ) && ( Tstrain > eunld_p ) && ( Tstrain < ere_p ) && ( cp == 1 ) )
   {
     Trule = 8;
     Frule = 9;
     cp = 0;
   }
   if( ( Crule == 12 ) && ( strain_incr > SMALL ) && ( Tstrain > eunld_p ) && ( Tstrain < ere_p ) && ( cp == 1 ) )
   {
     Trule = 8;
     Frule = 12;
     cp = 0;
   }
   if( ( Crule == 65 ) && ( strain_incr > SMALL ) && ( Tstrain > eunld_p8 ) && ( cp == 1 ) )
   {
     Trule = 8;
     Frule = 65;
     cp = 0;
   }
   if( ( Crule == 6 ) && ( strain_incr > SMALL ) && ( Tstrain > eunld_p8 ) && ( cp == 1 ) )
   {
     Trule = 8;
     Frule = 6;
     cp = 0;
   }
   if( ( Crule == 8 ) && ( Frule == 9 ) && ( strain_incr > SMALL ) && ( Tstrain < ere_p ) && ( cp == 1 ) )
   {
     Trule = 8;
     cp = 0;
   }
   if( ( Crule == 8 ) && ( Frule == 12 ) && ( strain_incr > SMALL ) && ( Tstrain < ere_p ) && ( cp == 1 ) )
   {
     Trule = 8;
     cp = 0;
   }
   if( ( Crule == 8 ) && ( Frule == 65 ) && ( strain_incr > SMALL ) && ( Tstrain < ere_p ) && ( cp == 1 ) )
   {
     Trule = 8;
     cp = 0;
   }
   if( ( Crule == 8 ) && ( Frule == 6 ) && ( strain_incr > SMALL ) && ( Tstrain < ere_p ) && ( cp == 1 ) )
   {
     Trule = 8;
     cp = 0;
   }
   if( ( Crule == 8 ) && ( Frule == 4 ) && ( strain_incr > SMALL ) && ( Tstrain < ere_str_p ) && ( cp == 1 ) )
   {
     Trule = 8;
     cp = 0;
   }
   if( ( Trule == 8 ) && ( Frule == 4 ) && ( Tstrain < eunld_p ) && ( cp == 0 ) )
   {
     Tstress =  stress_tran (Tstrain, er1, fr1, Ec, ere_str_p, fre_str_p, Ere_str_p);
     Ttangent = tangent_tran (Tstrain, er1, fr1, Ec, ere_str_p, fre_str_p, Ere_str_p);
   }
   if( ( Trule == 8 ) && ( Frule == 4 ) && ( Tstrain > eunld_p ) && ( Tstrain < ere_str_p ) && ( cp == 0 ) )
   {
     Tstress = stress_tran (Tstrain, eunld_p, fnew_str_p, Enew_str_p, ere_str_p, fre_str_p, Ere_str_p);
     Ttangent = tangent_tran (Tstrain, eunld_p, fnew_str_p, Enew_str_p, ere_str_p, fre_str_p, Ere_str_p);
   }
   if( ( Trule == 8 ) && ( Frule == 4 ) && ( Tstrain > ere_str_p ) && ( cp == 0 ) )
   {
     Tstress = stress_envlp_p(Tstrain-eo_p);
     Ttangent = tangent_envlp_p(Tstrain-eo_p);
   }
   if( ( Trule == 8 ) && ( Frule == 9 ) && ( Tstrain > eunld_p ) && ( Tstrain < ere_p ) && ( cp == 0 ) )
   {
     Tstress = stress_tran(Tstrain, eunld_p, fnew_p, Enew_p, ere_p, fre_p, Ere_p);
     Ttangent = tangent_tran(Tstrain, eunld_p, fnew_p, Enew_p, ere_p, fre_p, Ere_p);
   }
   if( ( Trule == 8 ) && ( Frule == 12 ) && ( Tstrain > eunld_p ) && ( Tstrain < ere_p ) && ( cp == 0 ) )
   {
     Tstress = stress_tran(Tstrain, eunld_p, fnew_p, Enew_p, ere_p, fre_p, Ere_p);
     Ttangent = tangent_tran(Tstrain, eunld_p, fnew_p, Enew_p, ere_p, fre_p, Ere_p);
   }
   if( ( Trule == 8 ) && ( Frule == 65 ) && ( Tstrain > eunld_p ) && ( Tstrain < ere_p ) && ( cp == 0 ) )
   {
     Tstress = stress_tran(Tstrain, eunld_p8, funld_p8, Eunld_p8, ere_p, fre_p, Ere_p);
     Ttangent = tangent_tran(Tstrain, eunld_p8, funld_p8, Eunld_p8, ere_p, fre_p, Ere_p);
   }
   if( ( Trule == 8 ) && ( Frule == 6 ) && ( Tstrain > eunld_p ) && ( Tstrain < ere_p ) && ( cp == 0 ) )
   {
     Tstress = stress_tran(Tstrain, eunld_p8, funld_p8, Eunld_p8, ere_p, fre_p, Ere_p);
     Ttangent = tangent_tran(Tstrain, eunld_p8, funld_p8, Eunld_p8, ere_p, fre_p, Ere_p);
   }
   if( ( Trule == 8 ) && ( Frule == 9 ) && ( Tstrain > ere_p ) && ( cp == 0 ) )
   {
     Tstress = stress_envlp_p(Tstrain-eo_p);
     Ttangent = tangent_envlp_p(Tstrain-eo_p);
   }
   if( ( Trule == 8 ) && ( Frule == 12 ) && ( Tstrain > ere_p ) && ( cp == 0 ) )
   {
     Tstress = stress_envlp_p(Tstrain-eo_p);
     Ttangent = tangent_envlp_p(Tstrain-eo_p);
   }
   /*****************************************RULE 14************************************************/
   if( ( Crule == 13 ) && ( strain_incr > SMALL ) && ( Tstrain > eunld_n ) && ( crack == 1 ) && ( cp == 1 ) )
   {
     Trule = 14;
     Frule = 13;
     ea = Cstrain;
     fa = Cstress;
     Ea = Ctangent;
     eb = ea - fa / Esec_n;
     fb = 0;
     Eb = 0;
     if( Tstrain > eb ){
       Trule = 16;
     }	     
     cp = 0;
   }
   if( ( Crule == 15 ) && ( strain_incr > SMALL ) && ( Tstrain < eb ) && ( crack == 1 ) && ( cp == 1 ) )
   {
     Trule = 14;
     Frule = 15;
     er = Cstrain;
     fr = Cstress;
     cp = 0;
   }
   if( ( Crule == 14 ) && ( strain_incr > SMALL ) && ( Tstrain < eb ) && ( crack == 1 ) && ( cp == 1 ) )
   {
     Trule = 14;
     cp = 0;
   }
   if ( (Trule == 14 ) && ( Frule == 13 ) && ( cp == 0 ) )
   {
     if( history_n == 0 ){
     Tstress = ( fa / ( ea - eb ) )  * ( Tstrain - eb ); 
     Ttangent = ( fa / ( ea - eb ) );
     }
     else if( history_n == 1 ){
     Tstress = 0.00000001;
     Ttangent = 0.00000001;
     }
   }
   if ( (Trule == 14 ) && ( Frule == 15 ) && ( cp == 0 ) )
   {
     Tstress = ( fr / ( er - eb ) )  * ( Tstrain - eb );
     Ttangent = ( fr / ( er - eb ) );
   }
   /**********************************************RULE 15*****************************************************/
   if( ( Crule == 14 ) && ( strain_incr < -SMALL ) && ( Tstrain > ea ) && ( Tstrain < eb ) && ( crack == 1 ) && ( cp == 1 ) )
   {
     Trule = 15;
     er = Cstrain;
     fr = Cstress;
     cp = 0;
   }
   if( ( Crule == 15 ) && ( strain_incr < -SMALL ) && ( Tstrain > ea ) && ( Tstrain < eb ) && ( crack == 1 ) && ( cp == 1 ) )
   {
     Trule = 15;
     cp = 0;
   }
   if ( ( Trule == 15 ) && ( cp == 0 ) )
   {
     //Tstress = stress_tran (Tstrain, er, fr, Ec, ea, fa, Ea);
     //Ttangent = tangent_tran (Tstrain, er, fr, Ec, ea, fa, Ea);
     Tstress = fa + ( ( fa - fr ) / ( ea - er ) ) * ( Tstrain - ea ); 
     Ttangent = ( fa - fr ) / ( ea - er );	     
   }
   cp = 1;
   if( ( Tstress < Fc_n ) | ( Tstress > Fc_p ) ) {
   //concrete<<Tstrain<<"   "<<Tstress<<"    "<<Ttangent<<"    "<<Trule<<"   "<<Cstrain<<"  "<<Cstress<<"  "<<Ctangent<<"  "<<Crule<<"  "<<"eunld_n "<<eunld_n<<"   "<<"funld_n "<<funld_n<<"   "<<"fnew_str_n "<<fnew_str_n<<"   "<<"ere_str_n "<<ere_str_n<<" fre_str_n "<<fre_str_n<<"   "<<"epl_n "<<epl_n<<"   "<<"er  "<<er<<"  "<<"fr  "<<fr<<"  "<<"ere_n "<<ere_n<<"  "<<"fre_n "<<fre_n<<"  "<<"ea "<<ea<<"  "<<"fa "<<fa<<"  "<<"epl_p "<<epl_p<<"  "<<" Ere_n "<<Ere_n<<"Epl_p "<<Epl_p<<"   "<<"eb"<<eb<<"  "<<"fb "<<fb<<"  "<<"fnew_n "<<fnew_n<<"  "<<" crack "<<crack<<" eunld_p "<<eunld_p<<" funld_p "<<funld_p<<" de_n "<<de_n<<"  df_n "<<df_n<<" ec_n "<<ec_n<<" Ec "<<Ec<<"  "<<"Enew_n"<<"  "<<Enew_n<<"  "<<"Enew_p"<<"  "<<Enew_p<<"  "<<"Eb"<<"  "<<Eb<<"  "<<"Ea"<<"  "<<Ea<<"  "<<"Epl_n"<<"  "<<Epl_n<<" ere_p "<<ere_p<<" fre_p "<<fre_p<<" Ere_p "<<Ere_p<<"  "<<" ere_str_p "<<ere_str_p<<"  "<<" fre_str_p "<<fre_str_p<<"  "<<" Ere_str_p "<<Ere_str_p<<" fnew_str_p "<<fnew_str_p<<" Esec_n "<<Esec_n<<" ea1 "<<ea1<<" fa1 "<<fa1<<" fb1 "<<fb1<<" eb1 "<<eb1<<" fnew_p "<<fnew_p<<" Frule "<<Frule<<" er1 "<<er1<<" fr1 "<<fr1<<endl;
   }
  return 0;
}
double 
RCFT_concMaterial::getStress(void)
{
  return Tstress;
}

double 
RCFT_concMaterial::getTangent(void)
{
  return Ttangent;
}

double 
RCFT_concMaterial::getInitialTangent(void)
{
  return Ec;
}

double 
RCFT_concMaterial::getStrain(void)
{
  return Tstrain;
}

int 
RCFT_concMaterial::commitState(void)
{
  ofstream conc;
  conc.open("conc.dat",ios::app);  
  //conc<<Tstrain<<"   "<<Tstress<<"    "<<Ttangent<<"    "<<Trule<<"   "<<Cstrain<<"  "<<Cstress<<"  "<<Ctangent<<"  "<<Crule<<endl;
  //conc<<eunld_p8<<"  "<<funld_p8<<"   "<<Eunld_p8<<"   "<<ere_p<<"   "<<fre_p<<"   "<<Ere_p<<endl;
  //conc<<Ec_sft<<"  "<<Fc_res_n<<endl;
  //if( ( Tstress < Fc_n ) | ( Tstress > Fc_p ) ) {
  //conc<<Tstrain<<"   "<<Tstress<<"    "<<Ttangent<<"    "<<Trule<<"   "<<Cstrain<<"  "<<Cstress<<"  "<<Ctangent<<"  "<<Crule<<"  "<<"eunld_n "<<eunld_n<<"   "<<"funld_n "<<funld_n<<"   "<<"fnew_str_n "<<fnew_str_n<<"   "<<"ere_str_n "<<"ere_str_n"<<" fre_str_n "<<fre_str_n<<"   "<<"epl_n "<<epl_n<<"   "<<"er  "<<er<<"  "<<"fr  "<<fr<<"  "<<"ere_n "<<ere_n<<"  "<<"fre_n "<<fre_n<<"  "<<"ea "<<ea<<"  "<<"fa "<<fa<<"  "<<"epl_p "<<epl_p<<"  "<<" Ere_n "<<Ere_n<<"Epl_p "<<Epl_p<<"   "<<"eb"<<eb<<"  "<<"fb "<<fb<<"  "<<"fnew_n "<<fnew_n<<"  "<<" crack "<<crack<<" eunld_p "<<eunld_p<<" funld_p "<<funld_p<<" de_n "<<de_n<<"  df_n "<<df_n<<" ec_n "<<ec_n<<" Ec "<<Ec<<"  "<<"Enew_n"<<"  "<<Enew_n<<"  "<<"Enew_p"<<"  "<<Enew_p<<"  "<<"Eb"<<"  "<<Eb<<"  "<<"Ea"<<"  "<<Ea<<"  "<<"Epl_n"<<"  "<<Epl_n<<" ere_p "<<ere_p<<" fre_p "<<fre_p<<" Ere_p "<<Ere_p<<"  "<<" ere_str_p "<<ere_str_p<<"  "<<" fre_str_p "<<fre_str_p<<"  "<<" Ere_str_p "<<Ere_str_p<<" fnew_str_p "<<fnew_str_p<<" Esec_n "<<Esec_n<<" ea1 "<<ea1<<" fa1 "<<fa1<<" fb1 "<<fb1<<" eb1 "<<eb1<<" fnew_p "<<fnew_p<<" Frule "<<Frule<<" er1 "<<er1<<" fr1 "<<fr1<<" eunld_n7 "<<eunld_n7<<" funld_n7 "<<funld_n7<<"Ec_sft "<<Ec_sft<<" ec_re_n "<<ec_res_n<<" Fc_res_n "<<Fc_res_n<<" ecr "<<ecr<<" Fc_p "<<Fc_p<<"  "<<"ec_p"<<"  "<<ec_p<<"  "<<"r_p"<<"  "<<r_p<<endl;
  //}
  Cstrain = Tstrain;
  Cstress = Tstress;
  Ctangent = Ttangent;
  Crule = Trule;
  Cstart_flag = Tstart_flag;
  Ceunld_n = eunld_n;
  Cfunld_n = funld_n;
  this->commitStatevar();
  return 0;
}

int 
RCFT_concMaterial::revertToLastCommit(void)
{
  Tstrain = Cstrain;
  Tstress = Cstress;
  Ttangent = Ctangent;
  Trule = Crule;
  Tstart_flag = Cstart_flag;
  eunld_n = Ceunld_n;
  funld_n = Cfunld_n;
  this->backtocommitStatevar();
  return 0;
}

int 
RCFT_concMaterial::revertToStart(void)
{
  //Trule = 0;
  //Tstrain = 0;
  //Tstress = 0;
  //Ttangent = 0;
  return Trule;
}

double
RCFT_concMaterial::stress_envlp_n(double strn)
{
  double strs, x_n, n_n, Dcr_n, ycr_n, zcr_n, x_sp, D, y ;
  n_n = fabs(Ec * ec_n / Fc_n);
  if( r_n != 1 )
  {
      Dcr_n = 1 + ( n_n - r_n / ( r_n - 1 ) ) * xcr_n + pow ( xcr_n , r_n ) / ( r_n - 1 );
  }
  else if( r_n == 1 )
  {
      Dcr_n = 1 + ( n_n - 1 + log( xcr_n ) ) * xcr_n;
  }
  ycr_n = n_n * xcr_n / Dcr_n;
  zcr_n = ( 1 - pow( xcr_n, r_n ) ) / ( Dcr_n * Dcr_n);
  x_sp = xcr_n - ycr_n / ( n_n * zcr_n );
  x_n = fabs(strn / ec_n);
  if(r_n != 1)
  {
     D = 1 + ( n_n - r_n / ( r_n - 1 ) ) * x_n + pow ( x_n, r_n ) / ( r_n - 1 );
  }
  if(r_n == 1)
  {
     D = 1 + ( n_n - 1 + log( x_n ) ) * x_n;
  }
  y = n_n * x_n / D;
  if( strn >= ec_n )
  {
     strs = Fc_n * y;
  }
  else if( ( strn < ec_n ) && ( strn >= ec_res_n ) )
  {
     strs = Fc_n + Ec_sft * ( strn - ec_n );
  }
  else if( strn < ec_res_n )
  {
     strs = Fc_res_n + (0.000001) * ( strn - ec_res_n );
  }
  return strs;
}

double 
RCFT_concMaterial::tangent_envlp_n(double strn)
{
  double tgnt, x_n, n_n, Dcr_n, ycr_n, zcr_n, x_sp, D, z ;
  n_n = fabs(Ec * ec_n / Fc_n);
  if( r_n != 1 )
  {
      Dcr_n = 1 + ( n_n - r_n / ( r_n - 1 ) ) * xcr_n + pow ( xcr_n , r_n ) / ( r_n - 1 );
  }
  else if( r_n == 1 )
  {
      Dcr_n = 1 + ( n_n - 1 + log( xcr_n ) ) * xcr_n;
  }
  ycr_n = n_n * xcr_n / Dcr_n;
  zcr_n = ( 1 - pow( xcr_n, r_n ) ) / ( Dcr_n * Dcr_n);
  x_sp = xcr_n - ycr_n / ( n_n * zcr_n );
  x_n = fabs(strn / ec_n);
  if(r_n != 1)
  {
     D = 1 + ( n_n - r_n / ( r_n - 1 ) ) * x_n + pow ( x_n, r_n ) / ( r_n - 1 );
  }
  if(r_n == 1)
  {
     D = 1 + ( n_n - 1 + log( x_n ) ) * x_n;
  }
  z = ( 1 - pow( x_n, r_n ) ) / ( D * D );
  if(  strn >= ec_n )
  {
     tgnt = Ec * z;
  }
  else if( ( strn < ec_n ) && ( strn >= ec_res_n ) )
  {
     tgnt = Ec_sft;
  }
  else if( strn < ec_res_n )
  {
     tgnt = 0.0001;
  }
  return tgnt;
}

double
RCFT_concMaterial::stress_envlp_p(double strn)
{
  double strs, x_p, n_p, Dcr_p, ycr_p, zcr_p, x_crk, D, y ;
  n_p = fabs(Ec * ec_p / Fc_p);
  if( r_p != 1 )
  {
      Dcr_p = 1 + ( n_p - r_p / ( r_p - 1 ) ) * xcr_p + pow ( xcr_p , r_p ) / ( r_p - 1 );
  }
  else if( r_p == 1 )
  {
      Dcr_p = 1 + ( n_p - 1 + log( xcr_p ) ) * xcr_p;
  }
  ycr_p = n_p * xcr_p / Dcr_p;
  zcr_p = ( 1 - pow( xcr_p, r_p ) ) / ( Dcr_p * Dcr_p );
  x_crk = xcr_p - ycr_p / ( n_p * zcr_p );
  x_p = fabs(strn / ec_p);
  if(r_p != 1)
  {
     D = 1 + ( n_p - r_p / ( r_p - 1 ) ) * x_p + pow ( x_p, r_p ) / ( r_p - 1 );
  }
  if(r_p == 1)
  {
     D = 1 + ( n_p - 1 + log( x_p ) ) * x_p;
  }
  y = n_p * x_p / D;
  if( x_p < xcr_p )
  {
     strs = Fc_p * y;
  }
  else if( ( x_p >= xcr_p ) && ( x_p <= x_crk ) )
  {
     strs = Fc_p * ( ycr_p + n_p * zcr_p * ( x_p - xcr_p ) );
  }
  else if( x_p > x_crk )
  {
     strs = 0.000001;
     crack = 1;
  }
  return strs;
}
                                                                                                                             
double
RCFT_concMaterial::tangent_envlp_p(double strn)
{
  double tgnt, x_p, n_p, Dcr_p, ycr_p, zcr_p, x_crk, D, z ;
  n_p = fabs(Ec * ec_p / Fc_p);
  if( r_p != 1 )
  {
      Dcr_p = 1 + ( n_p - r_p / ( r_p - 1 ) ) * xcr_p + pow ( xcr_p , r_p ) / ( r_p - 1 );
  }
  else if( r_p == 1 )
  {
      Dcr_p = 1 + ( n_p - 1 + log( xcr_p ) ) * xcr_p;
  }
  ycr_p = n_p * xcr_p / Dcr_p;
  zcr_p = ( 1 - pow( xcr_p, r_p ) ) / ( Dcr_p * Dcr_p );
  x_crk = xcr_p - ycr_p / ( n_p * zcr_p );
  x_p = fabs(strn / ec_p);
  if(r_p != 1)
  {
     D = 1 + ( n_p - r_p / ( r_p - 1 ) ) * x_p + pow ( x_p, r_p ) / ( r_p - 1 );
  }
  if(r_p == 1)
  {
     D = 1 + ( n_p - 1 + log( x_p ) ) * x_p;
  }
  z = ( 1 - pow( x_p, r_p ) ) / ( D * D );
  if( x_p < xcr_p )
  {
     tgnt = Ec * z;
  }
  else if( ( x_p >= xcr_p ) && ( x_p <= x_crk ) )
  {
     tgnt = Ec * zcr_p;
  }
  else if( x_p > x_crk )
  {
     tgnt = 0.000001;
  }
  return tgnt;
}

double 
RCFT_concMaterial::stress_tran(double strn, double ei, double fi, double Ei, double ef, double ff, double Ef)
{
  ofstream concrete;
  concrete.open("values.dat",ios::app);
  double R, A, Esec, strs;
  const double SMALL = 1.0e-30;
  Esec = ( ff - fi ) / ( ef - ei );
  if( (Trule == 3) | (Trule == 5) | (Trule == 4 ) | ( Trule == 6 ) ){
        if( Ei < 1.2*Esec ){
            Ei = 1.2*Esec;
        }
        if( Ef > Esec){
            Ef = Esec;
        }
  }
  if( (Trule == 7) | (Trule == 55) | ( Trule == 65 ) ){
	if( Ei < 1.2*Esec ){
	    Ei = 1.2*Esec;
	}
	if( Ef > Esec){
	    Ef = Esec;
	}
  }
  if( Trule == 8 ){
  	if( Ei > 0 ){
           if( Ei < 1.2*Esec ){
  	      Ei = 1.2*Esec;
           }
  	   if( Ef > Esec){
  	      Ef = Esec;
  	   }
  	}
  	else if( Ei < 0 ){
  	   if( Ei > 0.8*Esec ){
  	      Ei = 0.8*Esec;
  	   }
  	   if( Ef > Esec){
  	      Ef = Esec;
  	   }
  	}
  } 
  if( (Trule == 13) | (Trule == 10) ){
        if( (Ei > 0.8*Esec) && (Esec>0.0000001) ){
            Ei = 0.8*Esec;
        }
        if( (Ef < Esec) && (Esec>0.0000001) ){
            Ef = Esec;
        }
	if( (Ef > 2.0*Esec) && (Esec>0.0000001) ){
	    Ef = 2.0*Esec;
	}
  }
  R = ( Ef - Esec ) / ( Esec - Ei );
  //concrete<<strn<<"   "<<R<<"  "<<scientific<<pow(fabs(ef-ei),R)<<endl;
  if( ( Ef / Esec < 0.0 ) | ( Trule == 8 ) | ( Trule == 5 ) )
  {	 
     strs = fi + Esec* ( strn - ei ); 
     //concrete<<"linear"<<strs<<endl;
  }
  else 
  { 
     A = ( Esec - Ei ) / pow( fabs(ef - ei), R );
     strs = fi + ( strn - ei ) * ( Ei + A * pow( fabs( strn - ei ), R ) );
     //concrete<<"non-linear"<<strs<<endl;
  }
  //A = ( Esec - Ei ) / pow( fabs(ef - ei), R );
  //strs = fi + ( strn - ei ) * ( Ei + A * pow( fabs( strn - ei ), R ) );
  return strs;
}

double
RCFT_concMaterial::stress_tran3(double strn, double ei, double fi, double Ei, double ef, double ff, double Ef)
{
   double A, B, C, N, strs;
   A = ei;
   B = Ei;
   N = ( Ei - Ef ) * ( ef - ei ) / ( fi + Ei * ( ef - ei ) );
   C = ( Ef - Ei ) / ( N * pow( ef - ei , N - 1 ) );
   strs = A + B * ( strn - ei ) + C * pow( strn - ei, N - 1 );  
   return strs;   
}	

double
RCFT_concMaterial::tangent_tran3(double strn, double ei, double fi, double Ei, double ef, double ff, double Ef)
{
   double A, B, C, N, tgnt;
   A = ei;
   B = Ei;
   N = ( Ei - Ef ) * ( ef - ei ) / ( fi + Ei * ( ef - ei ) );
   C = ( Ef - Ei ) / ( N * pow( ef - ei , N - 1 ) );
   tgnt = B + C * ( N - 1 ) * pow( strn - ei, N - 1 );
   return tgnt;	
}

double
RCFT_concMaterial::tangent_tran(double strn, double ei, double fi, double Ei, double ef, double ff, double Ef)
{ 
  ofstream output;
  output.open("infi.dat",ios::app);  
  ofstream conc;
  conc.open("conc.dat",ios::app);
  double R, A, Esec, tgnt;
  const double SMALL = 1.0e-30;
  Esec = ( ff - fi ) / ( ef - ei );
  if( (Trule == 3) | (Trule == 5) | (Trule == 4 ) | ( Trule == 6 ) ){
     if( Ei < 1.2*Esec ){
        Ei = 1.2*Esec;
     }
     if( Ef > Esec){
        Ef = Esec;
     }
  } 
  if( (Trule == 7) | (Trule == 55) | ( Trule == 65 ) ){
     if( Ei < 1.2*Esec ){
        Ei = 1.2*Esec;
     }
     if( Ef > Esec){
        Ef = Esec;
     }
  }
  if( Trule == 8 ){
     if( Ei > 0 ){	  
       if( Ei < 1.2*Esec ){
          Ei = 1.2*Esec;
       }
       if( Ef > Esec){
          Ef = Esec;
       }
     }
     else if( Ei < 0 ){
       if( Ei > 0.8*Esec ){
          Ei = 0.8*Esec;
       }
       if( Ef > Esec){
          Ef = Esec;
       }
     }	     
  }
  if( (Trule == 13) | (Trule == 10) ){
     if( (Ei > 0.8*Esec) && (Esec>0.0000001) ){
        Ei = 0.8*Esec;
     }
     if( (Ef < Esec) && (Esec>0.0000001) ){
        Ef = Esec;
     }
     if( (Ef > 2.0*Esec) && (Esec>0.0000001) ){
        Ef = 2.0*Esec;
     }
  }
  R = ( Ef - Esec ) / ( Esec - Ei );
  //conc<<"___________________________"<<endl;
  //conc<<strn<<"  "<<R<<"  "<<pow(fabs(ef-ei),R )<<"   "<<( Esec - Ei ) / pow( fabs(ef - ei), R )<<endl;
  if ( (Ef / Esec < 0.0) | ( Trule == 8) | ( Trule == 5 ) )
  {
     tgnt = Esec;
  }
  else
  {
     A = ( Esec - Ei ) / pow( fabs(ef - ei), R );
     tgnt = Ei + A * ( R + 1 ) * pow ( fabs( strn - ei ), R );
  }
  //A = ( Esec - Ei ) / pow( fabs(ef - ei), R );
  //output<<strn<<"  "<<R<<"  "<<A<<"   "<<Esec<<"  "<<Ei<<"  "<<Ef<<endl;
  //tgnt = Ei + A * ( R + 1 ) * pow ( fabs( strn - ei ), R );
  return tgnt;
}

double
RCFT_concMaterial::stress_tran_p(double strn, double ei, double fi, double Ei, double ef, double ff, double Ef)
{
  double aa, bb, cc, strs, T;
  T = ei*ei + ef*ef - 2*ei*ef;
  aa = -fi/T + ff/T + Ei/(ei-ef);
  bb = 2*ei*fi/T - 2*ei*ff/T + (-ei-ef)*Ei/(ei-ef);
  cc = -ef*(-ef+2*ei)*fi/T + ei*ei*ff/T + ei*ef*Ei/(ei-ef);
  strs = aa*strn*strn + bb*strn + cc;  
  return strs;
}

double
RCFT_concMaterial::tangent_tran_p(double strn, double ei, double fi, double Ei, double ef, double ff, double Ef)
{
  double aa, bb, tgnt, T;
  T = ei*ei + ef*ef - 2*ei*ef;
  aa = -fi/T + ff/T + Ei/(ei-ef);
  bb = 2*ei*fi/T - 2*ei*ff/T + (-ei-ef)*Ei/(ei-ef);
  tgnt = 2*aa*strn+bb;
  return tgnt;
}

UniaxialMaterial *
RCFT_concMaterial::getCopy(void)
{
  RCFT_concMaterial *theCopy = new RCFT_concMaterial(this->getTag(), Fc, D, t, Fy, Es);
  theCopy->Ttangent = Ttangent;
  theCopy->Tstrain = Tstrain;
  theCopy->Tstress = Tstress;
  theCopy->Trule = Trule;
  theCopy->Frule = Frule;
  theCopy->Cstart_flag = Cstart_flag;
  theCopy->Tstart_flag = Tstart_flag;
  theCopy->crack= crack;
  theCopy->Fc_n = Fc_n;
  theCopy->Fc_p = Fc_p;
  theCopy->nn = nn;
  theCopy->np = np;
  theCopy->Ec = Ec;
  theCopy->ec_n = ec_n;
  theCopy->ec_p = ec_p;
  theCopy->Ec_sft = Ec_sft;
  theCopy->Fc_res_n = Fc_res_n;
  theCopy->ec_res_n = ec_res_n;
  theCopy->eo_n = eo_n;
  theCopy->eo_p = eo_p;
  theCopy->funld_n = funld_n;
  theCopy->eunld_n = eunld_n;
  theCopy->Esec_n = Esec_n;
  theCopy->Epl_n = Epl_n;
  theCopy->df_n = df_n;
  theCopy->de_n = de_n;
  theCopy->epl_n = epl_n;
  theCopy->fnew_n = fnew_n;
  theCopy->Enew_n = Enew_n;
  theCopy->ere_n = ere_n;  
  theCopy->fre_n = fre_n;
  theCopy->Ere_n = Ere_n;
  theCopy->fnew_str_n = fnew_str_n;
  theCopy->Enew_str_n = Enew_str_n;
  theCopy->ere_str_n = ere_str_n;
  theCopy->fre_str_n = fre_str_n;
  theCopy->Ere_str_n = Ere_str_n;
  theCopy->deo_n = deo_n;
  theCopy->funld_p = funld_p;
  theCopy->eunld_p = eunld_p;
  theCopy->Esec_p = Esec_p;
  theCopy->Epl_p = Epl_p;
  theCopy->df_p = df_p;
  theCopy->de_p = de_p;
  theCopy->epl_p = epl_p;
  theCopy->fnew_p = fnew_p;
  theCopy->Enew_p = Enew_p;
  theCopy->ere_p = ere_p;
  theCopy->fre_p = fre_p;
  theCopy->Ere_p = Ere_p;
  theCopy->xu_n = xu_n;
  theCopy->xu_p = xu_p;
  theCopy->deo_p = deo_p;
  theCopy->ea = ea;
  theCopy->fa = fa;
  theCopy->Ea = Ea;
  theCopy->er = er;
  theCopy->fr = fr;
  theCopy->eb = eb;
  theCopy->fb = fb;
  theCopy->Eb = Eb;
  theCopy->fnew_str_p = fnew_str_p;
  theCopy->Enew_str_p = Enew_str_p;
  theCopy->ere_str_p = ere_str_p;
  theCopy->fre_str_p = fre_str_p;
  theCopy->fre_str_p = fre_str_p;
  theCopy->Ere_str_p = Ere_str_p;
  theCopy->Cstrain = Cstrain;
  theCopy->Cstress = Cstress;
  theCopy->Crule = Crule;
  theCopy->Ctangent = Ctangent;
  theCopy->direction = direction;
  theCopy->r_n = r_n;
  theCopy->r_p = r_p;
  theCopy->xcr_n = xcr_n;
  theCopy->xcr_p = xcr_p;
  theCopy->history_n = history_n;
  theCopy->history_p = history_p;
  theCopy->Ccrack = Ccrack; theCopy->Cdirection = Cdirection; 
  theCopy->Ceo_n = Ceo_n; theCopy->Ceo_p = Ceo_p;
  theCopy->Cfunld_n = Cfunld_n; theCopy->Ceunld_n = Ceunld_n; theCopy->CEsec_n = CEsec_n;
  theCopy->CEpl_n = Cepl_n; theCopy->Cdf_n = Cdf_n; theCopy->Cde_n = Cde_n; theCopy->Cepl_n = Cepl_n;
  theCopy->Cfnew_n = Cfnew_n; theCopy->CEnew_n = CEnew_n; theCopy->Cere_n = Cere_n; theCopy->Cfre_n = Cfre_n;
  theCopy->CEre_n = CEre_n; theCopy->Cfnew_str_n = Cfnew_str_n; theCopy->CEnew_str_n = CEnew_str_n;
  theCopy->Cere_str_n = Cere_str_n; theCopy->Cfre_str_n = Cfre_str_n;  theCopy->CEre_str_n = CEre_str_n;
  theCopy->Cdeo_n = Cdeo_n; theCopy->Cfunld_p = Cfunld_n; theCopy->Ceunld_p = Ceunld_n; theCopy->CEsec_p = CEsec_p;
  theCopy->CEpl_p = CEpl_p; theCopy->Cdf_p = Cdf_p; theCopy->Cde_p = Cde_p; theCopy->Cepl_p = Cepl_p; theCopy->Cfnew_p = Cfnew_p;
  theCopy->CEnew_p = CEnew_p; theCopy->Cere_p = Cere_p; theCopy->Cfre_p = Cfre_p; theCopy->CEre_p = CEre_p; theCopy->Cxu_n = Cxu_n;
  theCopy->Cxu_p = Cxu_p; theCopy->Cdeo_p = Cdeo_p;  theCopy->Cea = Cea; theCopy->Cfa = Cfa; theCopy->CEa = CEa; theCopy->Cer = Cer; 
  theCopy->Cfr = Cfr;
  theCopy->Ceb = Ceb; theCopy->Cfb = Cfb; theCopy->CEb = CEb; theCopy->Cfnew_str_p = Cfnew_str_p; theCopy->CEnew_str_p = CEnew_str_p;
  theCopy->Cere_str_p = Cere_str_p; theCopy->Cfre_str_p = Cfre_str_p; theCopy->CEre_str_p = CEre_str_p;
  theCopy->Chistory_n = Chistory_n; theCopy->Chistory_p = Chistory_p; theCopy->CFrule = CFrule; theCopy->Cstart_flag = Cstart_flag;
  theCopy->Cea1 = Cea1; theCopy->ea1 = ea1; theCopy->Ceb1 = Ceb1; theCopy->eb1 = eb1;
  theCopy->Cfa1 = Cfa1; theCopy->fa1 = fa1; theCopy->Cfb1 = Cfb1; theCopy->fb1 = fb1;
  theCopy->CEa1 = CEa1; theCopy->Ea1 = Ea1; theCopy->CEb1 = CEb1; theCopy->Eb1 = Eb1;
  theCopy->CUrule = CUrule; theCopy->Urule = Urule;
  theCopy->CEr = CEr; theCopy->Er = Er;
  theCopy->eunld_n7 = eunld_n7;
  theCopy->Ceunld_n7 = Ceunld_n7;
  theCopy->funld_n7 = funld_n7;
  theCopy->Cfunld_n7 = Cfunld_n7; 
  theCopy->Eunld_n7 = Eunld_n7;
  theCopy->CEunld_n7 = CEunld_n7;
  theCopy->eunld_p8 = eunld_p8;
  theCopy->Ceunld_p8 = Ceunld_p8; 
  theCopy->funld_p8 = funld_p8;
  theCopy->Cfunld_p8 = Cfunld_p8; 
  theCopy->Eunld_p8 = Eunld_p8;
  theCopy->CEunld_p8 = CEunld_p8;
  theCopy->er1 = er1;
  theCopy->Cer1 = Cer1; 
  theCopy->fr1 = fr1;
  theCopy->Cfr1 = Cfr1;
  theCopy->Er1 = Er1;
  theCopy->CEr1 = CEr1;
  return theCopy;
}
int 
RCFT_concMaterial::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int 
RCFT_concMaterial::recvSelf(int cTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  return 0;	
}

void 
RCFT_concMaterial::Print(OPS_Stream &s, int flag)
{
  s<<"RCFT_conc, tag:"<<this->getTag()<<endln;
  s<<" Fcn:"<<Fc_n<<endln;
  s<<" Fcp:"<<Fc_p<<endln;
  s<<" Ec: "<<Ec<<endln;
  s<<" eo_n:"<<eo_n<<endln;
  s<<" eo_p:"<<eo_p<<endln;
  s<<" Ec_sft: "<<Ec_sft<<endln;
  return;
}

void
RCFT_concMaterial::backtocommitStatevar(void){
  crack = Ccrack; direction = Cdirection; eo_n = Ceo_n; eo_p = Ceo_p;
  funld_n = Cfunld_n; eunld_n = Ceunld_n; Esec_n = CEsec_n;
  Epl_n = CEpl_n; df_n = Cdf_n; de_n = Cde_n; epl_n = Cepl_n;
  fnew_n = Cfnew_n; Enew_n = CEnew_n; ere_n = Cere_n; fre_n = Cfre_n;
  Ere_n = CEre_n; fnew_str_n = Cfnew_str_n; Enew_str_n = CEnew_str_n; 
  ere_str_n = Cere_str_n; fre_str_n = Cfre_str_n; Ere_str_n = CEre_str_n;
  deo_n = Cdeo_n; funld_p = Cfunld_p; eunld_p = Ceunld_p; Esec_p = CEsec_p; 
  Epl_p = CEpl_p; df_p = Cdf_p; de_p = Cde_p; epl_p = Cepl_p; fnew_p = Cfnew_p;
  Enew_p = CEnew_p; ere_p = Cere_p; fre_p = Cfre_p; Ere_p = CEre_p; xu_n = Cxu_n;
  xu_p = Cxu_p; deo_p = Cdeo_p; ea = Cea; fa = Cfa; Ea = CEa; er = Cer; fr = Cfr;
  eb = Ceb; fb = Cfb; Eb = CEb; fnew_str_p = Cfnew_str_p; Enew_str_p = CEnew_str_p; 
  ere_str_p = Cere_str_p; fre_str_p = Cfre_str_p; Ere_str_p = CEre_str_p;
  history_n = Chistory_n; history_p = Chistory_p; Frule = CFrule; Tstart_flag = Cstart_flag;
  ea1 = Cea1; eb1 = Ceb1; fa1 = Cfa1; fb1 = Cfb1;
  Er = CEr;
  Urule == CUrule;
  eunld_n7 = Ceunld_n7; funld_n7 = Cfunld_n7; Eunld_n7 = CEunld_n7;
  eunld_p8 = Ceunld_p8; funld_p8 = Cfunld_p8; Eunld_p8 = CEunld_p8;
  er1 = Cer1; fr1 = Cfr1; Er1 = CEr1;
}

void
RCFT_concMaterial::commitStatevar(void){
  Ccrack = crack; Cdirection = direction; Ceo_n = eo_n; Ceo_p = eo_p;
  Cfunld_n = funld_n; Ceunld_n = eunld_n; CEsec_n = Esec_n;
  CEpl_n = Epl_n; Cdf_n = df_n; Cde_n = de_n; Cepl_n = epl_n;
  Cfnew_n = fnew_n; CEnew_n = Enew_n; Cere_n = ere_n; Cfre_n = fre_n;
  CEre_n = Ere_n; Cfnew_str_n = fnew_str_n; CEnew_str_n = Enew_str_n;
  Cere_str_n = ere_str_n; Cfre_str_n = fre_str_n; CEre_str_n = Ere_str_n;
  Cdeo_n = deo_n; Cfunld_p = funld_p; Ceunld_p = eunld_p; CEsec_p = Esec_p;
  CEpl_p = Epl_p; Cdf_p = df_p; Cde_p = de_p; Cepl_p = epl_p; Cfnew_p = fnew_p;
  CEnew_p = Enew_p; Cere_p = ere_p; Cfre_p = fre_p; CEre_p = Ere_p; Cxu_n = xu_n;
  Cxu_p = xu_p; Cdeo_p = deo_p; Cea = ea; Cfa = fa; CEa = Ea; Cer = er; Cfr = fr;
  Ceb = eb; Cfb = fb; CEb = Eb; Cfnew_str_p = fnew_str_p; CEnew_str_p = Enew_str_p;
  Cere_str_p = ere_str_p; Cfre_str_p = fre_str_p; CEre_str_p = Ere_str_p;
  Chistory_n = history_n; Chistory_p = history_p; CFrule = Frule; Cstart_flag = Tstart_flag;
  Cea1 = ea1; Ceb1 = eb1; Cfa1 = fa1; Cfb1 = fb1; CEr = Er;
  Cer1 = er1; Cfr1 = fr1; CEr1 = Er1;
  CUrule = Urule;
  Ceunld_n7 = eunld_n7; Cfunld_n7 = funld_n7; CEunld_n7 = Eunld_n7;
  Ceunld_p8 = eunld_p8; Cfunld_p8 = funld_p8; CEunld_p8 = Eunld_p8;
}

