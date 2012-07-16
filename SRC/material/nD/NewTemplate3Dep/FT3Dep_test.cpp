///////////////////////////////////////////////////////////////////////////////
//   COPYLEFT (C): Woody's viral GPL-like license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              
// CLASS:             
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           
//
// RETURN:
// VERSION:
// LANGUAGE:          C++
// TARGET OS:         
// DESIGNER:          Zhao Cheng, Boris Jeremic
// PROGRAMMER:        Zhao Cheng, 
// DATE:              Fall 2005
// UPDATE HISTORY:    
//
///////////////////////////////////////////////////////////////////////////////
//


#include <string.h>
#include <math.h>

#include <fstream.h>
#include <iostream.h>

#include "NewTemplate3Dep.h"

#include "ElasticState.h"
#include "Isotropic_Elastic.h"
#include "elnp_Elastic.h"
#include "DM04_Elastic.h"

#include "YieldFunction.h"
#include "DP_YF.h"
#include "VM_YF.h"
#include "RMC_YF.h"
#include "CC_YF.h"
#include "DM04_YF.h"

#include "PlasticFlow.h"
#include "DP_PF.h"
#include "VM_PF.h"
#include "RMC_PF.h"
#include "CC_PF.h"
#include "DM04_PF.h"

#include "ScalarEvolution.h"
#include "Linear_Eeq.h"
#include "CC_Ev.h"

#include "TensorEvolution.h"
#include "Linear_Eij.h"
#include "AF_Eij.h"
#include "DM04_alpha_Eij.h"
#include "DM04_z_Eij.h"

#include <OPS_Globals.h>
#include <G3Globals.h>
#include <ConsoleErrorHandler.h>
#include <OPS_Stream.h>

//Guanzhou Suggested 20mart2007
#include <StandardStream.h>
#include <FileStream.h>
StandardStream sserr;
OPS_Stream &opserr = sserr;


ErrorHandler *g3ErrorHandler =0;
double        ops_Dt = 0.0;
Domain       *ops_TheActiveDomain = 0;
Element      *ops_TheActiveElement = 0;

OPS_Stream *opserrPtr;

int main() {

cout << "*** T E S T ***" << "\n";

ofstream outStress ("Results.txt");

/*
 double rho = 0.0;
 double e0 = 0.8;
 double G0 = 125;
 double nu = 0.05;
 double Pat = 100.0;
 double p_cut = 0.01;
 double Mc = 1.25;
 double c = 0.8;
 double lambda_c = 0.019;
 double xi = 0.7;
 double ec_ref = 0.934;
 double m = 0.01;
 double h0 = 7.05;
 double ch = 0.968;
 double nb = 1.1;
 double A0 = 0.704;
 double nd = 3.5;
 double z_max = 4.0;
 double cz = 600.0;
 stresstensor zeroT1;
 stresstensor zeroT2;
 stresstensor initStress;
 double p0 = 100;
 initStress.val(1,1) = -p0; initStress.val(2,2) = -p0; initStress.val(3,3) = -p0;

 //----------------1    2   3  4  5    6   7      8  9         10  11 12 13   14  15  16  17  18    19
 double MC[19] = {rho, e0, G0, nu, Pat, p_cut, Mc, c, lambda_c, xi, ec_ref, m, h0, ch, nb, A0, nd, z_max, cz};
 stresstensor TS[2] = {zeroT1, zeroT2};
 MaterialParameter matpar(MC,19, TS,2);
 DM04_Elastic le(3, 4, 5, 6, 2, initStress);
 DM04_YF dpy(0,12, 2,1);
 DM04_PF dpf(0, 2, 0,11, 0,9, 0,10, 0,5, 0,12, 0,7, 0,8, 0,16, 0,17, 2,1, 2,2 );
 DM04_alpha_Eij Eij1(2, 11, 9, 10, 5, 12, 7, 8, 15, 13, 14, 3, 1, 2 );
 DM04_z_Eij Eij2(2, 11, 9, 10, 5, 12, 7, 8, 16, 17, 19, 18, 1, 2 );
 TensorEvolution *TE[2] = {&Eij1, &Eij2};
 NewTemplate3Dep FTEP(1, &matpar, &le, &dpy, &dpf, TE, 1);

stresstensor thisE;
stresstensor tStress;
straintensor tStrain;

double E11 = 0.0;    double E12 = 0.0;    double E13 = 0.0;
double E21 = 0.0;    double E22 = 0.0;    double E23 = 0.0;
double E31 = 0.0;    double E32 = 0.0;    double E33 = 0.0;

int Num_cyc = 1;
int I_cyc = 0;

//double ee = 0.0;
double pp = 0.0;
double qq = 0.0;
double tau = 0.0;

//double  q_cut = 40.0;
double ep_cut = 0.01;

double d_g = 1e-5;
 
 do {   
    // undrained tri-axial compression 
    E11 += (-d_g);    E22 += (d_g*0.5);    E33 += (d_g*0.5);
    
    // anisotropic compression shearing
    //E13 += d_g; E31 += d_g;
    
    thisE.val(1,1) = E11; thisE.val(1,2) = E12; thisE.val(1,3) = E13;
    thisE.val(2,1) = E21; thisE.val(2,2) = E22; thisE.val(2,3) = E23;
    thisE.val(3,1) = E31; thisE.val(3,2) = E32; thisE.val(3,3) = E33;

    FTEP.setTrialStrain(thisE);

    tStress = FTEP.getStressTensor();
    tStrain = FTEP.getStrainTensor();
    
    pp = tStress.p_hydrostatic();
    
    qq = - tStress.cval(1,1) + tStress.cval(2,2);
    //tau = tStress.val(1,3);
    //qq = tStress.q_deviatoric();

    outStress << -tStrain.cval(1,1)  <<  "  "  <<  pp << "  " << qq << "\n";
    //outStress << tStrain.val(1,3)  <<  "  "  << -tStress.val(3,3) << "  " << tau << "\n";
    //outStress << tStrain.val(1,3)  <<  "  "  << pp << "  " << qq << "\n";
    
    if (fabs(E11) >= ep_cut) {
    //if (fabs(qq) >= q_cut) {
    //if (fabs(tau) >= q_cut) {
    //if (fabs(E13) >= ep_cut) {
        I_cyc++;
        d_g *= (-1.0);
    }

 }  while ( I_cyc < Num_cyc );

 return 0;

}
*/


 // Test for Von_Mises
 double rho = 0.0;
 double E = 1.0e5;
 double v = 0.25;
 double k = 50.0;
 double H = 2.0e4;
 double Hij = 0.0;
 double ha = 5000.0;
 double Cr = 500.0;
 stresstensor zero;
 stresstensor initStress;
 double p0 = 100;
 initStress.val(1,1) = -p0; initStress.val(2,2) = -p0; initStress.val(3,3) = -p0;
 double MC[7] = {rho, E, v, H, Hij, ha, Cr};
 double IS[1] = {k};
 stresstensor TS[1] = {zero};
 MaterialParameter matpar(MC,7, IS,1, TS,1);
 
 Isotropic_Elastic le(2, 3, initStress);
 VM_YF dpy(1,1, 2,1);
 VM_PF dpf(2,1);
 Linear_Eeq Eep(4);
 Linear_Eij Eij(5);
 //AF_Eij Eij(6, 7, 1);
 ScalarEvolution *SE = {&Eep};
 TensorEvolution *TE = {&Eij};

 NewTemplate3Dep FTEP(1, &matpar, &le, &dpy, &dpf, &SE, &TE,  1);

 int Num_step = 40;
 double d_g = 2.0e-4;

 stresstensor thisE;
 stresstensor tStress;
 straintensor tStrain;
 
 double pp = 0.0;
 double qq = 0.0;

 tensor aStress(2, def_dim_2, 0.0);

 for ( int i = 0; i <= Num_step; i++ )
 {

 //Uniaxial Loading
 double E11 = d_g*double(i);    double E12 = 0.0;    double E13 = 0.0;
 double E21 = 0.0;    double E22 = -d_g*double(i)*0.5;    double E23 = 0.0;
 double E31 = 0.0;    double E32 = 0.0;    double E33 = -d_g*double(i)*0.5;
 //

 thisE.val(1,1) = E11; thisE.val(1,2) = E12; thisE.val(1,3) = E13;
 thisE.val(2,1) = E21; thisE.val(2,2) = E22; thisE.val(2,3) = E23;
 thisE.val(3,1) = E31; thisE.val(3,2) = E32; thisE.val(3,3) = E33;
 // thisE.print("F","F:");

 FTEP.setTrialStrain(thisE);
 
 tStress = FTEP.getStressTensor();
 tStrain = FTEP.getStrainTensor();
    
 pp  = tStress.p_hydrostatic();
 qq = tStress.val(1,1) - tStress.val(2,2);

 outStress << -tStrain.val(1,1)  <<  "  "  <<  pp << "  " << qq << endln;

 }

 return 0;
 
 }



/*
 // Test for DP

 //double rho = 0.0;
 //double E = 1.0e4;
 //double v = 0.25;
 //double M = 1.2;
 //double H1 = 100.0;
 //double H2 = 100.0; 
 //double p0 = 10.0;
 //stresstensor initStress;
 //initStress.val(1,1) = -p0; initStress.val(2,2) = -p0; initStress.val(3,3) = -p0;
 //stresstensor zeroT;
 //double MC[5] = {rho, E, v, H1, H2};
 //double IS[1] = {M};
 //stresstensor TS[1] = {zeroT};
 //MaterialParameter inpar(MC,5, IS,1, TS,1);
 //Isotropic_Elastic le(2, 3, initStress);
 //DP_YF dpy(1,1, 2,1);
 //DP_PF dpf(1,1, 2,1);
 //Linear_Eeq Eeq(4);
 //ScalarEvolution *SE = {&Eeq};
 //Linear_Eij Eij(5);
 //TensorEvolution *TE = {&Eij};

 double rho = 0.0;
 double E = 1.0e4;
 double v = 0.25;
 double M = 0.8;
 double H2 = 10.0;
 double ha = 20.0;
 double Cr = 2.0;  
 double p0 = 100.0;
 stresstensor initStress;
 initStress.val(1,1) = -p0; initStress.val(2,2) = -p0; initStress.val(3,3) = -p0;
 stresstensor zeroT;
 double MC[7] = {rho, E, v, M, H2, ha, Cr};
 stresstensor TS[1] = {zeroT};
 MaterialParameter inpar(MC,7, TS,1);
 Isotropic_Elastic le(2, 3, initStress);
 DP_YF dpy(0,4, 2,1);
 DP_PF dpf(0,4, 2,1);
 //Linear_Eij Eij(5);
 AF_Eij Eij(6,7, 1);
 TensorEvolution *TE = {&Eij};
                                                                            
 NewTemplate3Dep FTEP(1, &inpar, &le, &dpy, &dpf, &TE, 1);

 int Num_step = 40;
 double d_g = -5.0e-4;

 stresstensor thisE;
 stresstensor tStress;
 straintensor tStrain;
 
 double pp = 0.0;
 double qq = 0.0;

 tensor aStress(2, def_dim_2, 0.0);

 for ( int i = 0; i <= Num_step; i++ )
 {

 // Uniaxial Loading
 double E11 = d_g*double(i);    double E12 = 0.0;    double E13 = 0.0;
 double E21 = 0.0;    double E22 = -d_g*double(i)*0.5;    double E23 = 0.0;
 double E31 = 0.0;    double E32 = 0.0;    double E33 = -d_g*double(i)*0.5;
 //

 thisE.val(1,1) = E11; thisE.val(1,2) = E12; thisE.val(1,3) = E13;
 thisE.val(2,1) = E21; thisE.val(2,2) = E22; thisE.val(2,3) = E23;
 thisE.val(3,1) = E31; thisE.val(3,2) = E32; thisE.val(3,3) = E33;

 FTEP.setTrialStrain(thisE);
 
 tStress = FTEP.getStressTensor();
 tStrain = FTEP.getStrainTensor();
    
 pp  = tStress.p_hydrostatic();
 qq = -tStress.val(1,1) + tStress.val(2,2);

 outStress << -tStrain.val(1,1)  <<  "  "  <<  pp << "  " << qq << endln;

 }

 return 0;
}
*/

/*
 // Test for Cam Clay 

 double rho = 0.0;
 double e0 = 0.8;
 double M = 0.8;
 double lambda = 0.15;
 double kappa = 0.05;
 double v = 0.25;
 double Kc = 200.0;
 double p0 = 200.0;
 stresstensor initStress;
 initStress.val(1,1) = -p0/1.1; initStress.val(2,2) = -p0/1.1; initStress.val(3,3) = -p0/1.1;

 stresstensor zero;
 double MC[7] = {rho, e0, M, lambda, kappa, v, Kc};
 double IS[1] = {p0};
 MaterialParameter inpar(MC,7, IS,1);
 elnp_Elastic le(5, 6, 7, 2, initStress);
 CC_YF dpy(0,3, 1,1);
 CC_PF dpf(0,3, 1,1);
 CC_Ev Ev(3, 4, 5, 2, 1);
 ScalarEvolution *SE = {&Ev};
                                                                            
 NewTemplate3Dep FTEP(1, &inpar, &le, &dpy, &dpf, &SE, 1);

 int Num_step = 1200;
 double d_g = -0.0001;

 stresstensor thisE;
 stresstensor tStress;
 straintensor tStrain;
 
 double pp = 0.0;
 double qq = 0.0;

 tensor aStress(2, def_dim_2, 0.0);

 for ( int i = 0; i <= Num_step; i++ )
 {

 //Uniaxial Loading
 double E11 = d_g*double(i);    double E12 = 0.0;    double E13 = 0.0;
 double E21 = 0.0;    double E22 = -d_g*double(i)*0.5;    double E23 = 0.0;
 double E31 = 0.0;    double E32 = 0.0;    double E33 = -d_g*double(i)*0.5;
 //

 thisE.val(1,1) = E11; thisE.val(1,2) = E12; thisE.val(1,3) = E13;
 thisE.val(2,1) = E21; thisE.val(2,2) = E22; thisE.val(2,3) = E23;
 thisE.val(3,1) = E31; thisE.val(3,2) = E32; thisE.val(3,3) = E33;

 FTEP.setTrialStrain(thisE);
 
 tStress = FTEP.getStressTensor();
 tStrain = FTEP.getStrainTensor();
    
 pp  = tStress.p_hydrostatic();
 qq = -tStress.val(1,1) + tStress.val(2,2);

 outStress << -tStrain.val(1,1)  <<  "  "  <<  pp << "  " << qq << endln;

 }

 return 0;
}
*/


//  //*********************   Simple Shear ******************************
//  double E11 = -0.001;    double E12 = 0.0;    double E13 = d_g*double(i);
//  double E21 = 0.0;    double E22 = -0.001;    double E23 = 0.0;
//  double E31 = d_g*double(i);    double E32 = 0.0;    double E33 = -0.001;
//  //*******************************************************************

 
// /***********************************************************************************/
// /*   Uniaxial                                                                      */
// double E11 = d_g*double(i);    double E12 = 0.0;    double E13 = 0.0;
// double E21 = 0.0;    double E22 = d_g*double(i)*(-v);    double E23 = 0.0;
// double E31 = 0.0;    double E32 = 0.0;    double E33 = d_g*double(i)*(-v); 

// /***********************************************************************************/
// /*   Triaxial Compression                                                          */
// double E11 = 0.0;    double E12 = 0.0;    double E13 = 0.0;
// double E21 = 0.0;    double E22 = 0.0;    double E23 = 0.0;
// double E31 = 0.0;    double E32 = 0.0;    double E33 = -d_g*double(i);

////*********************   Uniaxial Loading **************************
//double E11 = d_g*double(i);    double E12 = 0.0;    double E13 = 0.0;
//double E21 = 0.0;    double E22 = d_g*double(i)*(-v);    double E23 = 0.0;
//double E31 = 0.0;    double E32 = 0.0;    double E33 = d_g*double(i)*(-v);
////*******************************************************************

/*
 // Test for Von_Mises, Linear Isotropic, Error Map
 int i, j, k, err;
 double pi = 3.14159265358979323846;
 double twoO3 = 2.0/3.0; 
   
 // AF
 double rho = 0.0;
 double E = 1.0e5;
 double v = 0.25;
 double kk = 50.0;
 double H =  2.0e4;
 double Hij = 0.0;
 stresstensor zero;
 stresstensor initStress;
 double p0 = 100;
 initStress.val(1,1) = -p0; initStress.val(2,2) = -p0; initStress.val(3,3) = -p0;
 double MC[5] = {rho, E, v, H, Hij};
 double IS[1] = {kk};
 stresstensor TS[1] = {zero};
 MaterialParameter matpar(MC,5, IS,1, TS,1);
 Isotropic_Elastic le(2, 3, initStress);
 VM_YF dpy(1,1, 2,1);
 VM_PF dpf(2,1);
 Linear_Eeq Eep(4);
 Linear_Eij Eij(5);
 ScalarEvolution *SE = {&Eep};
 TensorEvolution *TE = {&Eij}; 
 
 stresstensor stressIncr;
 straintensor strainIncr;
 straintensor strainsubIncr;
 straintensor strainThis;
 stresstensor stressConv;
 tensor Ee(4, def_dim_4, 0.0);
 tensor Ce(4, def_dim_4, 0.0);
 
 double p = 0.0;
 double q = 0.0;
 double theta = 0.0;

 int n1 =  20;
 int n2 =  10;
 int n3 =  12;
 
 double Dpq = 9.99;

 double se1, se2, se3, se4, se5, se6;
 double si1, si2, si3, si4, si5, si6;
 double errmap = 0.0;
 double err1 = 0.0;
 double err2 = 0.0;

 NewTemplate3Dep *NTEP0 = new NewTemplate3Dep(1, &matpar, &le, &dpy, &dpf, &SE, &TE, 0);
 Ee = NTEP0->getTangentTensor();
 err += NTEP0->Stiffness2Compliance(Ee, Ce);
 delete NTEP0;

 for (i=0; i<=n1; i++) { 
   p = i*Dpq - Dpq*10;
   for (j=0; j<=n2; j++) {
	 q = j*Dpq;
     for (k=0; k<=n3; k++) {
	   theta = k*pi/3.0/n3;
	   stressIncr.val(1,1) = -p + twoO3 *q *cos(theta);
	   stressIncr.val(2,2) = -p + twoO3 *q *cos(theta - twoO3*pi);
	   stressIncr.val(3,3) = -p + twoO3 *q *cos(theta + twoO3*pi);
	   strainIncr = Ce("ijmn") * stressIncr("mn");
       strainIncr.null_indices();

       // implicit
       NewTemplate3Dep *NTEPi = new NewTemplate3Dep(2, &matpar, &le, &dpy, &dpf, &SE, &TE, 1);
	   NTEPi->setTrialStrain(strainIncr);
	   stressConv = NTEPi->getStressTensor();
       si1 = stressConv.cval(1,1);
       si2 = stressConv.cval(2,2);
       si3 = stressConv.cval(3,3);
       si4 = stressConv.cval(1,2);
       si5 = stressConv.cval(1,3);
       si6 = stressConv.cval(2,3);
       cout << si1 << " " << si2 << " " << si3 << endln;
       delete NTEPi;

       // explicit
       NewTemplate3Dep *NTEPe = new NewTemplate3Dep(3, &matpar, &le, &dpy, &dpf, &SE, &TE, 0, 50);
       NTEPe->setTrialStrain(strainIncr);
       stressConv = NTEPe->getStressTensor();
       se1 = stressConv.cval(1,1);
       se2 = stressConv.cval(2,2);
       se3 = stressConv.cval(3,3);
       se4 = stressConv.cval(1,2);
       se5 = stressConv.cval(1,3);
       se6 = stressConv.cval(2,3);
       cout << se1 << " " << se2 << " " << se3 << endln;
       delete NTEPe;

       err1 = sqrt( pow(se1-si1, 2) + pow(se2-si2, 2) + pow(se3-si3, 2) + 2.0*pow(se4-si4, 2) + 2.0*pow(se5-si5, 2) + 2.0*pow(se5-si5, 2) );
       err2 = sqrt ( pow(se1, 2) + pow(se2, 2) + pow(se3, 2) + 2.0*pow(se4, 2) + 2.0*pow(se5, 2) + 2.0*pow(se6, 2) );
       errmap = err1 / err2;	   
       outStress << p + p0 << " " << q << " " << theta << " " << errmap << " " << err1/p0/sqrt(3.0) << endln;
     }
   }
 } 

 return 0;
}
*/

/* 
 // Test for Von_Mises, AF, Error Map
 int i, j, k, err;
 double pi = 3.14159265358979323846;
 double twoO3 = 2.0/3.0; 
   
 // AF
 double rho = 0.0;
 double E = 1.0e5;
 double v = 0.25;
 double kk = 50.0;
 double H =  0.0;
 double Hij = 0.0;
 double ha = 5.0e4;
 double Cr = 2.5e3;
 stresstensor zero;
 stresstensor initStress;
 double p0 = 100;
 initStress.val(1,1) = -p0; initStress.val(2,2) = -p0; initStress.val(3,3) = -p0;
 double MC[7] = {rho, E, v, H, Hij, ha, Cr};
 double IS[1] = {kk};
 stresstensor TS[1] = {zero};
 MaterialParameter matpar(MC,7, IS,1, TS,1);
 Isotropic_Elastic le(2, 3, initStress);
 VM_YF dpy(1,1, 2,1);
 VM_PF dpf(2,1);
 Linear_Eeq Eep(4);
 //Linear_Eij Eij(5);
 AF_Eij Eij(6, 7, 1);
 ScalarEvolution *SE = {&Eep};
 TensorEvolution *TE = {&Eij}; 
 
 stresstensor stressIncr;
 straintensor strainIncr;
 straintensor strainsubIncr;
 straintensor strainThis;
 stresstensor stressConv;
 tensor Ee(4, def_dim_4, 0.0);
 tensor Ce(4, def_dim_4, 0.0);
 
 double p = 0.0;
 double q = 0.0;
 double theta = 0.0;

 int n1 =  20;
 int n2 =  10;
 int n3 =  12;
 
 double Dpq = 9.99;

 double se1, se2, se3, se4, se5, se6;
 double si1, si2, si3, si4, si5, si6;
 double errmap = 0.0;
 double err1 = 0.0;
 double err2 = 0.0;

 NewTemplate3Dep *NTEP0 = new NewTemplate3Dep(1, &matpar, &le, &dpy, &dpf, &SE, &TE, 0);
 Ee = NTEP0->getTangentTensor();
 err += NTEP0->Stiffness2Compliance(Ee, Ce);
 delete NTEP0;

 for (i=0; i<=n1; i++) { 
   p = i*Dpq - Dpq*10;
   for (j=0; j<=n2; j++) {
	 q = j*Dpq;
     for (k=0; k<=n3; k++) {
	   theta = k*pi/3.0/n3;
	   stressIncr.val(1,1) = -p + twoO3 *q *cos(theta);
	   stressIncr.val(2,2) = -p + twoO3 *q *cos(theta - twoO3*pi);
	   stressIncr.val(3,3) = -p + twoO3 *q *cos(theta + twoO3*pi);
	   strainIncr = Ce("ijmn") * stressIncr("mn");
       strainIncr.null_indices();

       // implicit
       NewTemplate3Dep *NTEPi = new NewTemplate3Dep(2, &matpar, &le, &dpy, &dpf, &SE, &TE, 1);
	   NTEPi->setTrialStrain(strainIncr);
	   stressConv = NTEPi->getStressTensor();
       si1 = stressConv.cval(1,1);
       si2 = stressConv.cval(2,2);
       si3 = stressConv.cval(3,3);
       si4 = stressConv.cval(1,2);
       si5 = stressConv.cval(1,3);
       si6 = stressConv.cval(2,3);
       cout << si1 << " " << si2 << " " << si3 << endln;
       delete NTEPi;

       // explicit
       NewTemplate3Dep *NTEPe = new NewTemplate3Dep(3, &matpar, &le, &dpy, &dpf, &SE, &TE, 0, 50);
       NTEPe->setTrialStrain(strainIncr);
       stressConv = NTEPe->getStressTensor();
       se1 = stressConv.cval(1,1);
       se2 = stressConv.cval(2,2);
       se3 = stressConv.cval(3,3);
       se4 = stressConv.cval(1,2);
       se5 = stressConv.cval(1,3);
       se6 = stressConv.cval(2,3);
       cout << se1 << " " << se2 << " " << se3 << endln;
       delete NTEPe;

       err1 = sqrt( pow(se1-si1, 2) + pow(se2-si2, 2) + pow(se3-si3, 2) + 2.0*pow(se4-si4, 2) + 2.0*pow(se5-si5, 2) + 2.0*pow(se5-si5, 2) );
       err2 = sqrt ( pow(se1, 2) + pow(se2, 2) + pow(se3, 2) + 2.0*pow(se4, 2) + 2.0*pow(se5, 2) + 2.0*pow(se6, 2) );
       errmap = err1 / err2;	   
       outStress << p + p0 << " " << q << " " << theta << " " << errmap << " " << err1/p0/sqrt(3.0) << endln;
     }
   }
 } 

 return 0;
}
*/

/*
 // Test for DP, PP, Error Map
 int i, j, k, err;
 double pi = 3.14159265358979323846;
 double twoO3 = 2.0/3.0;

 double rho = 0.0;
 double E = 1.0e4;
 double v = 0.25;
 double M = 0.8;
 double p0 = 100.0;
 stresstensor initStress;
 initStress.val(1,1) = -p0; initStress.val(2,2) = -p0; initStress.val(3,3) = -p0;
 stresstensor zeroT;
 double MC[4] = {rho, E, v, M};
 stresstensor TS[1] = {zeroT};
 MaterialParameter matpar(MC,4);
 Isotropic_Elastic le(2, 3, initStress);
 DP_YF dpy(0,4);
 DP_PF dpf(0,4);
 

 stresstensor stressIncr;
 straintensor strainIncr;
 straintensor strainsubIncr;
 straintensor strainThis;
 stresstensor stressConv;
 tensor Ee(4, def_dim_4, 0.0);
 tensor Ce(4, def_dim_4, 0.0);
 
 double p = 0.0;
 double q = 0.0;
 double theta = 0.0;

 int n1 =  20;
 int n2 =  10;
 int n3 =  12;
 
 double Dpq = 9.99;

 double se1, se2, se3, se4, se5, se6;
 double si1, si2, si3, si4, si5, si6;
 double errmap = 0.0;
 double err1 = 0.0;
 double err2 = 0.0;

 NewTemplate3Dep *NTEP0 = new NewTemplate3Dep(1, &matpar, &le, &dpy, &dpf, 0);
 Ee = NTEP0->getTangentTensor();
 err += NTEP0->Stiffness2Compliance(Ee, Ce);
 delete NTEP0;

 for (i=0; i<=n1; i++) { 
   p = i*Dpq - Dpq*10;
   for (j=0; j<=n2; j++) {
	 q = j*Dpq;
     for (k=0; k<=n3; k++) {
	   theta = k*pi/3.0/n3;
	   stressIncr.val(1,1) = -p + twoO3 *q *cos(theta);
	   stressIncr.val(2,2) = -p + twoO3 *q *cos(theta - twoO3*pi);
	   stressIncr.val(3,3) = -p + twoO3 *q *cos(theta + twoO3*pi);
	   strainIncr = Ce("ijmn") * stressIncr("mn");
       strainIncr.null_indices();

       // implicit
       NewTemplate3Dep *NTEPi = new NewTemplate3Dep(2, &matpar, &le, &dpy, &dpf, 1);
	   NTEPi->setTrialStrain(strainIncr);
	   stressConv = NTEPi->getStressTensor();
       si1 = stressConv.cval(1,1);
       si2 = stressConv.cval(2,2);
       si3 = stressConv.cval(3,3);
       si4 = stressConv.cval(1,2);
       si5 = stressConv.cval(1,3);
       si6 = stressConv.cval(2,3);
       cout << si1 << " " << si2 << " " << si3 << endln;
       delete NTEPi;

       // explicit
       NewTemplate3Dep *NTEPe = new NewTemplate3Dep(3, &matpar, &le, &dpy, &dpf, 0, 50);
       NTEPe->setTrialStrain(strainIncr);
       stressConv = NTEPe->getStressTensor();
       se1 = stressConv.cval(1,1);
       se2 = stressConv.cval(2,2);
       se3 = stressConv.cval(3,3);
       se4 = stressConv.cval(1,2);
       se5 = stressConv.cval(1,3);
       se6 = stressConv.cval(2,3);
       cout << se1 << " " << se2 << " " << se3 << endln;
       delete NTEPe;

       err1 = sqrt( pow(se1-si1, 2) + pow(se2-si2, 2) + pow(se3-si3, 2) + 2.0*pow(se4-si4, 2) + 2.0*pow(se5-si5, 2) + 2.0*pow(se5-si5, 2) );
       err2 = sqrt ( pow(se1, 2) + pow(se2, 2) + pow(se3, 2) + 2.0*pow(se4, 2) + 2.0*pow(se5, 2) + 2.0*pow(se6, 2) );
       errmap = err1 / err2;	   
       outStress << p + p0 << " " << q << " " << theta << " " << errmap << " " << err1/p0/sqrt(3.0) << endln;
     }
   }
 } 

 return 0; 
}
*/

/*
 // Test for DP, AF, Error Map
 int i, j, k, err;
 double pi = 3.14159265358979323846;
 double twoO3 = 2.0/3.0;

 double rho = 0.0;
 double E = 1.0e4;
 double v = 0.25;
 double M = 0.8;
 double H1 = 0.0;
 double H2 = 10.0; 
 double ha = 20.0;
 double Cr = 2.0;
 double p0 = 100.0;
 stresstensor initStress;
 initStress.val(1,1) = -p0; initStress.val(2,2) = -p0; initStress.val(3,3) = -p0;
 stresstensor zeroT;
 double MC[7] = {rho, E, v, H1, H2, ha, Cr};
 double IS[1] = {M};
 stresstensor TS[1] = {zeroT};
 MaterialParameter matpar(MC,7, IS,1, TS,1);
 Isotropic_Elastic le(2, 3, initStress);
 DP_YF dpy(1,1, 2,1);
 DP_PF dpf(1,1, 2,1);
 Linear_Eeq Eeq(4);
 ScalarEvolution *SE = {&Eeq};
 //Linear_Eij Eij(5);
 AF_Eij Eij(6, 7, 1);
 TensorEvolution *TE = {&Eij};
 

 stresstensor stressIncr;
 straintensor strainIncr;
 straintensor strainsubIncr;
 straintensor strainThis;
 stresstensor stressConv;
 tensor Ee(4, def_dim_4, 0.0);
 tensor Ce(4, def_dim_4, 0.0);
 
 double p = 0.0;
 double q = 0.0;
 double theta = 0.0;

 int n1 =  20;
 int n2 =  10;
 int n3 =  12;
 
 double Dpq = 9.99;

 double se1, se2, se3, se4, se5, se6;
 double si1, si2, si3, si4, si5, si6;
 double errmap = 0.0;
 double err1 = 0.0;
 double err2 = 0.0;

 NewTemplate3Dep *NTEP0 = new NewTemplate3Dep(1, &matpar, &le, &dpy, &dpf, &SE, &TE, 0);
 Ee = NTEP0->getTangentTensor();
 err += NTEP0->Stiffness2Compliance(Ee, Ce);
 delete NTEP0;

 for (i=0; i<=n1; i++) { 
   p = i*Dpq - Dpq*10;
   for (j=0; j<=n2; j++) {
	 q = j*Dpq;
     for (k=0; k<=n3; k++) {
	   theta = k*pi/3.0/n3;
	   stressIncr.val(1,1) = -p + twoO3 *q *cos(theta);
	   stressIncr.val(2,2) = -p + twoO3 *q *cos(theta - twoO3*pi);
	   stressIncr.val(3,3) = -p + twoO3 *q *cos(theta + twoO3*pi);
	   strainIncr = Ce("ijmn") * stressIncr("mn");
       strainIncr.null_indices();

       // implicit
       NewTemplate3Dep *NTEPi = new NewTemplate3Dep(2, &matpar, &le, &dpy, &dpf, &SE, &TE, 1);
	   NTEPi->setTrialStrain(strainIncr);
	   stressConv = NTEPi->getStressTensor();
       si1 = stressConv.cval(1,1);
       si2 = stressConv.cval(2,2);
       si3 = stressConv.cval(3,3);
       si4 = stressConv.cval(1,2);
       si5 = stressConv.cval(1,3);
       si6 = stressConv.cval(2,3);
       cout << si1 << " " << si2 << " " << si3 << endln;
       delete NTEPi;

       // explicit
       NewTemplate3Dep *NTEPe = new NewTemplate3Dep(3, &matpar, &le, &dpy, &dpf, &SE, &TE, 0, 50);
       NTEPe->setTrialStrain(strainIncr);
       stressConv = NTEPe->getStressTensor();
       se1 = stressConv.cval(1,1);
       se2 = stressConv.cval(2,2);
       se3 = stressConv.cval(3,3);
       se4 = stressConv.cval(1,2);
       se5 = stressConv.cval(1,3);
       se6 = stressConv.cval(2,3);
       cout << se1 << " " << se2 << " " << se3 << endln;
       delete NTEPe;

       err1 = sqrt( pow(se1-si1, 2) + pow(se2-si2, 2) + pow(se3-si3, 2) + 2.0*pow(se4-si4, 2) + 2.0*pow(se5-si5, 2) + 2.0*pow(se5-si5, 2) );
       err2 = sqrt ( pow(se1, 2) + pow(se2, 2) + pow(se3, 2) + 2.0*pow(se4, 2) + 2.0*pow(se5, 2) + 2.0*pow(se6, 2) );
       errmap = err1 / err2;	   
       outStress << p + p0 << " " << q << " " << theta << " " << errmap << " " << err1/p0/sqrt(3.0) << endln;
     }
   }
 } 

 return 0; 
}
*/

/*

 // Test for CC, Error Map
 int i, j, k, err;
 double pi = 3.14159265358979323846;
 double twoO3 = 2.0/3.0;

 double rho = 0.0;
 double e0 = 0.8;
 double M = 0.8;
 double lambda = 0.12;
 double kappa = 0.04;
 double v = 0.25;
 double Kc = 1000.0;
 double p0 = 200.0;
 stresstensor initStress;
 initStress.val(1,1) = -100; initStress.val(2,2) = -100; initStress.val(3,3) = -100;

 
 stresstensor zero;
 double MC[7] = {rho, e0, M, lambda, kappa, v, Kc};
 double IS[1] = {p0};
 MaterialParameter matpar(MC,7, IS,1);
 elnp_Elastic le(5, 6, 7, 2, initStress);
 CC_YF dpy(0,3, 1,1);
 CC_PF dpf(0,3, 1,1);
 CC_Ev Ev(3, 4, 5, 2, 1);
 ScalarEvolution *SE = {&Ev}; 

 stresstensor stressIncr;
 straintensor strainIncr;
 straintensor strainsubIncr;
 straintensor strainThis;
 stresstensor stressConv;
 tensor Ee(4, def_dim_4, 0.0);
 tensor Ce(4, def_dim_4, 0.0);
 
 double p = 0.0;
 double q = 0.0;
 double theta = 0.0;

 int n1 =  20;
 int n2 =  10;
 int n3 =  12;
 
 double Dpq = 9.99;

 double se1, se2, se3, se4, se5, se6;
 double si1, si2, si3, si4, si5, si6;
 double errmap = 0.0;
 double err1 = 0.0;
 double err2 = 0.0;

 NewTemplate3Dep *NTEP0 = new NewTemplate3Dep(1, &matpar, &le, &dpy, &dpf, &SE, 0);
 Ee = NTEP0->getTangentTensor();
 err += NTEP0->Stiffness2Compliance(Ee, Ce);
 delete NTEP0;

 for (i=0; i<=n1; i++) { 
   p = i*Dpq - Dpq*10;
   for (j=0; j<=n2; j++) {
	 q = j*Dpq;
     for (k=0; k<=n3; k++) {
	   theta = k*pi/3.0/n3;
	   stressIncr.val(1,1) = -p + twoO3 *q *cos(theta);
	   stressIncr.val(2,2) = -p + twoO3 *q *cos(theta - twoO3*pi);
	   stressIncr.val(3,3) = -p + twoO3 *q *cos(theta + twoO3*pi);
	   strainIncr = Ce("ijmn") * stressIncr("mn");
       strainIncr.null_indices();

       // implicit
       NewTemplate3Dep *NTEPi = new NewTemplate3Dep(2, &matpar, &le, &dpy, &dpf, &SE, 1, 1);
	   NTEPi->setTrialStrain(strainIncr);
	   stressConv = NTEPi->getStressTensor();
       si1 = stressConv.cval(1,1);
       si2 = stressConv.cval(2,2);
       si3 = stressConv.cval(3,3);
       si4 = stressConv.cval(1,2);
       si5 = stressConv.cval(1,3);
       si6 = stressConv.cval(2,3);
       cout << si1 << " " << si2 << " " << si3 << endln;
       delete NTEPi;

       // explicit
       NewTemplate3Dep *NTEPe = new NewTemplate3Dep(3, &matpar, &le, &dpy, &dpf, &SE, 0, 50);
       NTEPe->setTrialStrain(strainIncr);
       stressConv = NTEPe->getStressTensor();
       se1 = stressConv.cval(1,1);
       se2 = stressConv.cval(2,2);
       se3 = stressConv.cval(3,3);
       se4 = stressConv.cval(1,2);
       se5 = stressConv.cval(1,3);
       se6 = stressConv.cval(2,3);
       cout << se1 << " " << se2 << " " << se3 << endln;
       delete NTEPe;

       err1 = sqrt( pow(se1-si1, 2) + pow(se2-si2, 2) + pow(se3-si3, 2) + 2.0*pow(se4-si4, 2) + 2.0*pow(se5-si5, 2) + 2.0*pow(se5-si5, 2) );
       err2 = sqrt ( pow(se1, 2) + pow(se2, 2) + pow(se3, 2) + 2.0*pow(se4, 2) + 2.0*pow(se5, 2) + 2.0*pow(se6, 2) );
       errmap = err1 / err2;	   
       outStress << p + p0 << " " << q << " " << theta << " " << errmap << " " << err1/p0/sqrt(3.0) << endln;
     }
   }
 } 

 return 0; 
}
*/


/*
 // Error Map
 int i, j, k, err;
 double pi = 3.14159265358979323846;
 double twoO3 = 2.0/3.0;

 double rho = 0.0;
 double e0 = 0.8;
 double G0 = 125;
 double nu = 0.05;
 double Pat = 100.0;
 double p_cut = 0.01;
 double Mc = 1.25;
 double c = 0.8;
 double lambda_c = 0.019;
 double xi = 0.7;
 double ec_ref = 0.934;
 double m = 0.01;
 double h0 = 7.05;
 double ch = 0.968;
 double nb = 1.1;
 double A0 = 0.704;
 double nd = 3.5;
 double z_max = 4.0;
 double cz = 600.0;
 stresstensor zeroT1;
 //zeroT1.val(1,1) = -(2.0*nu-1.0)/(1.0+nu);
 //zeroT1.val(2,2) = -(2.0*nu-1.0)/(1.0+nu);
 //zeroT1.val(3,3) = -(2.0-4.0*nu)/(1.0+nu);
 stresstensor zeroT2;
 stresstensor initStress;
 double p0 = 100;
 initStress.val(1,1) = -p0; initStress.val(2,2) = -p0; initStress.val(3,3) = -p0;
 //initStress.val(1,1) = -p*nu/(1.0-nu); initStress.val(2,2) = -p*nu/(1.0-nu); initStress.val(3,3) = -p;

 //DMImplicit NTEPi(1, rho, Pat, G0, nu, p_cut, e0, Mc, c, ec_ref, lambda_c, xi,
 //		   m, h0, ch, nb, A0, nd, z_max, cz, zeroT1, zeroT2, initStress); 


 //----------------1    2   3  4  5    6   7      8  9         10  11 12 13   14  15  16  17  18    19
 double MC[19] = {rho, e0, G0, nu, Pat, p_cut, Mc, c, lambda_c, xi, ec_ref, m, h0, ch, nb, A0, nd, z_max, cz};
 stresstensor TS[2] = {zeroT1, zeroT2};
 MaterialParameter matpar(MC,19, TS,2);
 DM04_Elastic le(3, 4, 5, 6, 2, initStress);
 DM04_YF dpy(0,12, 2,1);
 DM04_PF dpf(0, 2, 0,11, 0,9, 0,10, 0,5, 0,12, 0,7, 0,8, 0,16, 0,17, 2,1, 2,2 );
 DM04_alpha_Eij Eij1(2, 11, 9, 10, 5, 12, 7, 8, 15, 13, 14, 3, 1, 2 );
 DM04_z_Eij Eij2(2, 11, 9, 10, 5, 12, 7, 8, 16, 17, 19, 18, 1, 2 );
 TensorEvolution *TE[2] = {&Eij1, &Eij2};
 //NewTemplate3Dep NTEPe(2, &matpar, &le, &dpy, &dpf, TE);
 
 stresstensor stressIncr;
 straintensor strainIncr;

 straintensor strainsubIncr;
 straintensor strainThis;
 stresstensor stressConv;
 tensor Ee(4, def_dim_4, 0.0);
 tensor Ce(4, def_dim_4, 0.0);

 double p = 0.0;
 double q = 0.0;
 double theta = 0.0;

 // n1, n2, n3 > 1;
 int n1 =  20;
 int n2 =  10;
 int n3 =  12;
 
 double Dpq = 9.99;

 double se1, se2, se3, se4, se5, se6;
 double si1, si2, si3, si4, si5, si6;
 double errmap = 0.0;
 double err1 = 0.0;
 double err2 = 0.0;
 
 NewTemplate3Dep *NTEP0 = new NewTemplate3Dep(1, &matpar, &le, &dpy, &dpf, TE, 0);
 Ee = NTEP0->getTangentTensor();
 err += NTEP0->Stiffness2Compliance(Ee, Ce);
 delete NTEP0;

 for (i=0; i<=n1; i++) { 
   p = i*Dpq - 99.9;
   for (j=0; j<=n2; j++) {
	 q = j*Dpq;
     for (k=0; k<=n3; k++) {
	   theta = k*pi/3.0/n3;
	   stressIncr.val(1,1) = -p + twoO3 *q *cos(theta);
	   stressIncr.val(2,2) = -p + twoO3 *q *cos(theta - twoO3*pi);
	   stressIncr.val(3,3) = -p + twoO3 *q *cos(theta + twoO3*pi);
	   strainIncr = Ce("ijmn") * stressIncr("mn");
       strainIncr.null_indices();

       // implicit
       //DMImplicit NTEPi(1, rho, Pat, G0, nu, p_cut, e0, Mc, c, ec_ref, lambda_c, xi,
	   //	   m, h0, ch, nb, A0, nd, z_max, cz, zeroT1, zeroT2, initStress);
       NewTemplate3Dep *NTEPi = new NewTemplate3Dep(2, &matpar, &le, &dpy, &dpf, TE, 2);
       NTEPi->setTrialStrain(strainIncr);
	   stressConv = NTEPi->getStressTensor();
       si1 = stressConv.cval(1,1);
       si2 = stressConv.cval(2,2);
       si3 = stressConv.cval(3,3);
       si4 = stressConv.cval(1,2);
       si5 = stressConv.cval(1,3);
       si6 = stressConv.cval(2,3);
       cout << si1 << " " << si2 << " " << si3 << endln;
       delete NTEPi;

       // explicit
       NewTemplate3Dep *NTEPe = new NewTemplate3Dep(3, &matpar, &le, &dpy, &dpf, TE, 0, 50);
       NTEPe->setTrialStrain(strainIncr);
       stressConv = NTEPe->getStressTensor();
       se1 = stressConv.cval(1,1);
       se2 = stressConv.cval(2,2);
       se3 = stressConv.cval(3,3);
       se4 = stressConv.cval(1,2);
       se5 = stressConv.cval(1,3);
       se6 = stressConv.cval(2,3);
       cout << se1 << " " << se2 << " " << se3 << endln;
       delete NTEPe;

       err1 = sqrt( pow(se1-si1, 2) + pow(se2-si2, 2) + pow(se3-si3, 2) + 2.0*pow(se4-si4, 2) + 2.0*pow(se5-si5, 2) + 2.0*pow(se5-si5, 2) );
       err2 = sqrt ( pow(se1, 2) + pow(se2, 2) + pow(se3, 2) + 2.0*pow(se4, 2) + 2.0*pow(se5, 2) + 2.0*pow(se6, 2) );
       errmap = err1 / err2;	   
       outStress << p + p0 << " " << q << " " << theta << " " << errmap << " " << err1/p0/sqrt(3.0) << endln;
     }
   }
 } 

 return 0;
 
}
*/
