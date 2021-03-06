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
                                                                        
//Modified by:  Jian Zhang(j.zhang@ed.ac.uk)---------07,2010// 
//              Panagiotis Kotsovinos(P.Kotsovinos@ed.ac.uk)// 



#include <Steel01PThermal.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <math.h>
#include <float.h>

#include <elementAPI.h>
#include <OPS_Globals.h>


void *
OPS_NewSteel01PThermal()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[7];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel01PThermal tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 3 && numData != 7) {
    opserr << "Invalid #args, want: uniaxialMaterial Steel01PThermal " << iData[0] << " fy? E? b? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial Steel01PThermal " << iData[0] << " fy? E? b? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (numData == 3) {
    dData[3] = STEEL_01_DEFAULT_A1;
    dData[4] = STEEL_01_DEFAULT_A2;
    dData[5] = STEEL_01_DEFAULT_A3;
    dData[6] = STEEL_01_DEFAULT_A4;
  }

  // Parsing was successful, allocate the material
  theMaterial = new Steel01PThermal(iData[0], dData[0], dData[1], 
				   dData[2], dData[3], dData[4], 
				   dData[5], dData[6]);
  
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Steel01PThermal Material\n";
    return 0;
  }

  return theMaterial;
}


Steel01PThermal::Steel01PThermal
(int tag, double FY, double E, double B,
 double A1, double A2, double A3, double A4):
   UniaxialMaterial(tag,MAT_TAG_Steel01PThermal),
   fyT(FY), E0T(E), b(B), a1(A1), a2(A2), a3(A3), a4(A4)
{
   // Sets all history and state variables to initial values
   // History variables
   CminStrain = 0.0;
   CmaxStrain = 0.0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;

   TminStrain = 0.0;
   TmaxStrain = 0.0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   //Ctangent = E0;

   Ctangent = E0T;///JZ, 07/10//


   Tstrain = 0.0;
   Tstress = 0.0;
  // Ttangent = E0;

   Ttangent = E0T;///JZ, 07/10//

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////

	  ThermalElongation = 0; //initialize //JZ, 07/10//
	  E0 = E0T;//JZ, 07/10//
	  fy = fyT;//JZ, 07/10//
	  fp = 0;//JZ, 11/10//
	  TemperautreC = 0;
}

Steel01PThermal::Steel01PThermal():UniaxialMaterial(0,MAT_TAG_Steel01PThermal),
 //fy(0.0), E0(0.0), b(0.0), a1(0.0), a2(0.0), a3(0.0), a4(0.0)
 fyT(0.0), E0T(0.0), b(0.0), a1(0.0), a2(0.0), a3(0.0), a4(0.0) //JZ, 07/10//
{

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////

	  ThermalElongation = 0; //initialize //JZ, 07/10//
	  E0 = E0T;//JZ, 07/10//
	  fy = fyT;//JZ, 07/10//
	  fp = 0;//JZ, 11/10//
	  TemperautreC = 0;

}

Steel01PThermal::~Steel01PThermal ()
{
// AddingSensitivity:BEGIN /////////////////////////////////////
	if (SHVs != 0) 
		delete SHVs;
// AddingSensitivity:END //////////////////////////////////////
}

int Steel01PThermal::setTrialStrain(double strain, double FiberTemperature, double strainRate)
{

  Temp = FiberTemperature;


   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;
   
   if (fabs(dStrain) > DBL_EPSILON || dStrain == 0) {
     // Set trial strain
     Tstrain = strain;
     
     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);
     
   }

   return 0;
}

int Steel01PThermal::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;

   if (fabs(dStrain) > DBL_EPSILON) {
     // Set trial strain
     Tstrain = strain;

     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);

   }

   stress = Tstress;
   tangent = Ttangent;

   return 0;
}

void Steel01PThermal::determineTrialState (double dStrain)
{
      double fyOneMinusB = fy * (1.0 - b);

      double Esh = b*E0;
      double epsy = fy/E0;
      
      double c1 = Esh*Tstrain;
      
      double c2 = TshiftN*fyOneMinusB;

      double c3 = TshiftP*fyOneMinusB;

      double c = Cstress + E0*dStrain;
	  if (E0!=E0T) 
		   c = E0*(Cstrain+dStrain);
	  //Modified by Liming for considering decreasing E0.

      /**********************************************************
         removal of the following lines due to problems with
	 optimization may be required (e.g. on gnucc compiler
         with optimization turned on & -ffloat-store option not
         used) .. replace them with line that follows but which 
         now requires 2 function calls to achieve same result !!
      ************************************************************/

      double c1c3 = c1 + c3;

      if (c1c3 < c)
	Tstress = c1c3;
      else
	Tstress = c;

      double c1c2 = c1-c2;

      if (c1c2 > Tstress)
	Tstress = c1c2;

      /* ***********************************************************
      and replace them with:

      Tstress = fmax((c1-c2), fmin((c1+c3),c));
      **************************************************************/

      if (fabs(Tstress-c) < DBL_EPSILON)
	  Ttangent = E0;
      else
	Ttangent = Esh;

      //************JZ 11/10 S 
	  double EpsiPT = fp/E0;
	  double EpsiYT = 0.02;
	  double EpsiT = 0.2;
	  double EpsiU = 0.2;
	  double CT = (fy-fp)*(fy-fp)/((EpsiYT-EpsiPT)*E0-2*(fy - fp));
	  double BT = pow(CT*(EpsiYT-EpsiPT)*E0+CT*CT, 0.5);
	  double AT = pow((EpsiYT-EpsiPT)*(EpsiYT-EpsiPT+CT/E0),0.5);
	  
      if (Tloading == 0 && dStrain != 0.0) {
	  if (dStrain > 0.0)
	    Tloading = 1;
	  else
	    Tloading = -1;
      }

      // Transition from loading to unloading, i.e. positive strain increment
      // to negative strain increment
      if (Tloading == 1 && dStrain < 0.0) {
	  Tloading = -1;
	  if (Cstrain > TmaxStrain)
	    TmaxStrain = Cstrain;
	  TshiftN = 1 + a1*pow((TmaxStrain-TminStrain)/(2.0*a2*epsy),0.8);
      }

      // Transition from unloading to loading, i.e. negative strain increment
      // to positive strain increment
      if (Tloading == -1 && dStrain > 0.0) {
	  Tloading = 1;
	  if (Cstrain < TminStrain)
	    TminStrain = Cstrain;
	  TshiftP = 1 + a3*pow((TmaxStrain-TminStrain)/(2.0*a4*epsy),0.8);
      }

}

void Steel01PThermal::detectLoadReversal (double dStrain)
{
   // Determine initial loading condition
   if (Tloading == 0 && dStrain != 0.0)
   {
      if (dStrain > 0.0) 
         Tloading = 1;
      else
         Tloading = -1;
   }

   double epsy = fy/E0;

   // Transition from loading to unloading, i.e. positive strain increment
   // to negative strain increment
   if (Tloading == 1 && dStrain < 0.0)
   {
      Tloading = -1;
      if (Cstrain > TmaxStrain)
         TmaxStrain = Cstrain;
      TshiftN = 1 + a1*pow((TmaxStrain-TminStrain)/(2.0*a2*epsy),0.8);
   }

   // Transition from unloading to loading, i.e. negative strain increment
   // to positive strain increment
   if (Tloading == -1 && dStrain > 0.0)
   {
      Tloading = 1;
      if (Cstrain < TminStrain)
         TminStrain = Cstrain;
      TshiftP = 1 + a3*pow((TmaxStrain-TminStrain)/(2.0*a4*epsy),0.8);
   }
}

double Steel01PThermal::getStrain ()
{
   return Tstrain;
}

double Steel01PThermal::getStress ()
{
   return Tstress;
}

double Steel01PThermal::getTangent ()
{
   return Ttangent;
}

double 
Steel01PThermal::getThermalElongation(void) //***JZ
{
  return ThermalElongation;
}

//JZ 07/10 /////////////////////////////////////////////////////////////start
double 
Steel01PThermal::getElongTangent(double TempT, double &ET, double &Elong, double TempTmax) //PK add to include max temp
{
  //JZ updated, from rebar to C steel
  
  // EN 1992 pt 1-2-1. Class N hot rolled  reinforcing steel at elevated temperatures
  if (TempT <= 80) {
    fy = fyT;
    E0 = E0T*(1-TempT*(1-0.98)/80);
    
    //b=TempT*0.00325/80;
    
    fp = fyT*(1-TempT*(1-0.68)/80);
  }
  else if (TempT <= 180) {
    fy = fyT*(1 - (TempT-80)*(1-0.87)/100);
    E0 = E0T*(0.98 - (TempT - 80)*(0.98-0.95)/100);
    
    //b=0.00325+(TempT - 80)*0.00325/100;
    
    fp = fyT*(0.68 - (TempT - 80)*(0.68-0.51)/100);
    
  }
  else if (TempT <= 280) {
    fy = fyT*(0.87 - (TempT - 180)*(0.87-0.70)/100);
    E0 = E0T*(0.95 - (TempT - 180)*(0.95-0.88)/100);
    
    //b=0.0065+(TempT - 180)*0.00325/100;
    
    fp = fyT*(0.51 - (TempT - 180)*(0.51-0.32)/100);
  }
  else if (TempT <= 380) {
    fy = fyT*(0.70 - (TempT - 280)*(0.70-0.50)/100);
    E0 = E0T*(0.88 - (TempT - 280)*(0.88 - 0.81)/100);
    
    //b=0.00975+(TempT - 280)*0.00355/100;
    
    fp = fyT*(0.32 - (TempT - 280)*(0.32 - 0.13)/100);
  }
  else if (TempT <= 480) {
    fy = fyT*(0.50 - (TempT - 380)*(0.50-0.30)/100);
    E0 = E0T*(0.81 - (TempT - 380)*(0.81-0.54)/100);
    
    //b=0.0133+(TempT - 380)*0.0133/100;
    
    fp = fyT*(0.13 - (TempT - 380)*(0.13 - 0.07)/100);
  }
  else if (TempT <= 580) {
    fy = fyT*(0.30 - (TempT - 480)*(0.30-0.14)/100);
    E0 = E0T*(0.54 - (TempT - 480)*(0.54-0.41)/100);
    
    //b=0.0266+(TempT - 480)*0.0136/100;
    
    fp = fyT*(0.07 - (TempT - 480)*(0.07 - 0.05)/100);
  }
  else if (TempT <= 680) {
    fy = fyT*(0.14 - (TempT - 580)*(0.14-0.06)/100);
    E0 = E0T*(0.41 - (TempT - 580)*(0.41-0.10)/100);
    
    // b=0.0402-(TempT - 580)*0.0067/100;
    
    fp = fyT*(0.05 - (TempT - 580)*(0.05 - 0.03)/100);
  }
  else if (TempT <= 780) {
    fy = fyT*(0.06 - (TempT - 680)*(0.06-0.04)/100);
    E0 = E0T*(0.10 - (TempT - 680)*(0.10-0.07)/100);
    
    // b=0.0335-(TempT - 680)*0.0067/100;
    
    fp = fyT*(0.03 - (TempT - 680)*(0.03 - 0.02)/100);
  }
  else if (TempT <= 880) {
    fy = fyT*(0.04 - (TempT - 780)*(0.04-0.02)/100);
    E0 = E0T*(0.07 - (TempT - 780)*(0.07-0.03)/100);
    
    //  b=0.0268-(TempT - 780)*0.0067/100;
    
    fp = fyT*(0.02 - (TempT - 780)*(0.02 - 0.01)/100);
  }
  else if (TempT <= 980) {
    fy = fyT*(0.02 - (TempT - 880)*0.02/100);
    E0 = E0T*(0.03 - (TempT - 880)*0.03/100);
    
    //  b=0.0201-(TempT - 880)*0.0067/100;
    
    fp = fyT*(0.01- (TempT - 880)*0.01/100);
  }
  else if (TempT <= 1080) {
    fy = fyT*0;
    E0 = E0T*0;
    
    // b=0.0134-(TempT - 980)*0.0067/100;
    
    fp = fyT*0;
  }
  else if (TempT <= 1180) {
    fy = fyT*0;
    E0 = E0T*0;
    
    //  b=0.0067-(TempT - 980)*0.0067/100;
    
    fp = fyT*0;
  }
  else  {
    opserr << "the temperature is invalid\n"; 
  } 

  // caculation of thermal elongation of reinforcing steel. JZ
///*	
	if (TempT <= 1) {
		  ThermalElongation = TempT * 1.2164e-5;
	  }
  else if (TempT <= 730) {
      ThermalElongation = -2.416e-4 + 1.2e-5 *(TempT+20) + 0.4e-8 *(TempT+20)*(TempT+20);
  }
  else if (TempT <= 840) {
      ThermalElongation = 11e-3;
  }
  else if (TempT <= 1180) {
      ThermalElongation = -6.2e-3 + 2e-5*TempT;
  }
  else {
	  opserr << "the temperature is invalid\n";
  }
//*/
  //double Taf=(TempT+40)/2;
  //double alpha = (11.09+0.0062*Taf)/1000000;
  //ThermalElongation = TempT * alpha;


/*
// fp,fy,and E0 is taken from the Paper of A. Rubert
  if (TempT <= 80){
	  E0 = E0T;
  }
E0 = E0T*(1-0.001*(TempT-80));
 
if (TempT <= 380){
	fy = fyT;
  }
else if (TempT <= 680){
	fy = fyT*(1-0.8/300*(TempT-380));
  }
else if (TempT <= 980){
	fy = fyT*(0.2-0.2/300*(TempT-680));
  }
else { fy = 0;}



if (TempT <= 80){
	fp = fyT ;
   
  }
else if (TempT <= 180){
	fp = fyT*(1-0.001*(TempT-80));
  }
else if (TempT <= 280){
	fp = fyT*(0.9-0.003*(TempT-180));
  }
else if (TempT <= 480){
	fp = fyT*(0.6-0.001*(TempT-280));
  }
else if (TempT <= 630){
	fp = fyT*(0.4-4/1500*(TempT-480));
  }
else if (TempT <= 980){
	fp = fyT*(0.1-1/3500*(TempT-630));
  }
else {fp = 0;}
*/
  //ThermalElongation = ThermalElongation*(1e-6) ;   //debug  Liming
  ET = E0;   
  
  Elong = ThermalElongation;
  TemperautreC = TempT;

  //opserr << "\getelongation: " << ET << "\ temp:" << TemperautreC <<endln; //PK Check


  return 0;
}
//JZ 07/10 /////////////////////////////////////////////////////////////end 


int Steel01PThermal::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CmaxStrain = TmaxStrain;
   CshiftP = TshiftP;
   CshiftN = TshiftN;
   Cloading = Tloading;

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int Steel01PThermal::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;

   // Reset trial state variables to last committed state
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}

int Steel01PThermal::revertToStart ()
{
   // History variables
   CminStrain = 0.0;
   CmaxStrain = 0.0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;

   TminStrain = 0.0;
   TmaxStrain = 0.0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = E0;

   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = E0;

// AddingSensitivity:BEGIN /////////////////////////////////
	if (SHVs != 0) 
		SHVs->Zero();
// AddingSensitivity:END //////////////////////////////////

   return 0;
}

UniaxialMaterial* Steel01PThermal::getCopy ()
{
   Steel01PThermal* theCopy = new Steel01PThermal(this->getTag(), fy, E0, b,
				  a1, a2, a3, a4);

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CmaxStrain = CmaxStrain;
   theCopy->CshiftP = CshiftP;
   theCopy->CshiftN = CshiftN;
   theCopy->Cloading = Cloading;

   // Trial history variables
   theCopy->TminStrain = TminStrain;
   theCopy->TmaxStrain = TmaxStrain;
   theCopy->TshiftP = TshiftP;
   theCopy->TshiftN = TshiftN;
   theCopy->Tloading = Tloading;

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

int Steel01PThermal::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(16);
   data(0) = this->getTag();

   // Material properties
   data(1) = fy;
   data(2) = E0;
   data(3) = b;
   data(4) = a1;
   data(5) = a2;
   data(6) = a3;
   data(7) = a4;

   // History variables from last converged state
   data(8) = CminStrain;
   data(9) = CmaxStrain;
   data(10) = CshiftP;
   data(11) = CshiftN;
   data(12) = Cloading;

   // State variables from last converged state
   data(13) = Cstrain;
   data(14) = Cstress;
   data(15) = Ctangent;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "Steel01PThermal::sendSelf() - failed to send data\n";

   return res;
}

int Steel01PThermal::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(16);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
   if (res < 0) {
      opserr << "Steel01PThermal::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties
      fy = data(1);
      E0 = data(2);
      b = data(3);
      a1 = data(4);
      a2 = data(5);
      a3 = data(6);
      a4 = data(7);

      // History variables from last converged state
      CminStrain = data(8);
      CmaxStrain = data(9);
      CshiftP = data(10);
      CshiftN = data(11);
      Cloading = int(data(12));

      // Copy converged history values into trial values since data is only
      // sent (received) after convergence
      TminStrain = CminStrain;
      TmaxStrain = CmaxStrain;
      TshiftP = CshiftP;
      TshiftN = CshiftN;
      Tloading = Cloading;

      // State variables from last converged state
      Cstrain = data(13);
      Cstress = data(14);
      Ctangent = data(15);      

      // Copy converged state values into trial values
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }
    
   return res;
}

void Steel01PThermal::Print (OPS_Stream& s, int flag)
{
   s << "Steel01PThermal tag: " << this->getTag() << endln;
   s << "  fy: " << fy << " ";
   s << "  E0: " << E0 << " ";
   s << "  b:  " << b << " ";
   s << "  a1: " << a1 << " ";
   s << "  a2: " << a2 << " ";
   s << "  a3: " << a3 << " ";
   s << "  a4: " << a4 << " ";
}




// AddingSensitivity:BEGIN ///////////////////////////////////
int
Steel01PThermal::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"sigmaY") == 0 || strcmp(argv[0],"fy") == 0)
    return param.addObject(1, this);
  
  if (strcmp(argv[0],"E") == 0)
    return param.addObject(2, this);
  
  if (strcmp(argv[0],"b") == 0)
    return param.addObject(3, this);
  
  if (strcmp(argv[0],"a1") == 0)
    return param.addObject(4, this);
  
  if (strcmp(argv[0],"a2") == 0)
    return param.addObject(5, this);
  
  if (strcmp(argv[0],"a3") == 0)
    return param.addObject(6, this);
  
  if (strcmp(argv[0],"a4") == 0)
    return param.addObject(7, this);

  return -1;
}



int
Steel01PThermal::updateParameter(int parameterID, Information &info)
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
	case 4:
		this->a1 = info.theDouble;
		break;
	case 5:
		this->a2 = info.theDouble;
		break;
	case 6:
		this->a3 = info.theDouble;
		break;
	case 7:
		this->a4 = info.theDouble;
		break;
	default:
		return -1;
	}

	Ttangent = E0;          // Initial stiffness

	return 0;
}




int
Steel01PThermal::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}



double
Steel01PThermal::getStressSensitivity(int gradIndex, bool conditional)
{
	// Initialize return value
	double gradient = 0.0;


	// Pick up sensitivity history variables
	double CstrainSensitivity = 0.0;
	double CstressSensitivity = 0.0;
	if (SHVs != 0) {
		CstrainSensitivity = (*SHVs)(0,gradIndex);
		CstressSensitivity = (*SHVs)(1,gradIndex);
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


	// Compute min and max stress
	double Tstress;
	double dStrain = Tstrain-Cstrain;
	double sigmaElastic = Cstress + E0*dStrain;
	double fyOneMinusB = fy * (1.0 - b);
	double Esh = b*E0;
	double c1 = Esh*Tstrain;
	double c2 = TshiftN*fyOneMinusB;
	double c3 = TshiftP*fyOneMinusB;
	double sigmaMax = c1+c3;
	double sigmaMin = c1-c2;


	// Evaluate stress sensitivity 
	if ( (sigmaMax < sigmaElastic) && (fabs(sigmaMax-sigmaElastic)>1e-5) ) {
		Tstress = sigmaMax;
		gradient = E0Sensitivity*b*Tstrain 
				 + E0*bSensitivity*Tstrain
				 + TshiftP*(fySensitivity*(1-b)-fy*bSensitivity);
	}
	else {
		Tstress = sigmaElastic;
		gradient = CstressSensitivity 
			     + E0Sensitivity*(Tstrain-Cstrain)
				 - E0*CstrainSensitivity;
	}
	if (sigmaMin > Tstress) {
		gradient = E0Sensitivity*b*Tstrain
			     + E0*bSensitivity*Tstrain
				 - TshiftN*(fySensitivity*(1-b)-fy*bSensitivity);
	}

	return gradient;
}




double
Steel01PThermal::getInitialTangentSensitivity(int gradIndex)
{
	// For now, assume that this is only called for initial stiffness 
	if (parameterID == 2) {
		return 1.0; 
	}
	else {
		return 0.0;
	}
}


int
Steel01PThermal::commitSensitivity(double TstrainSensitivity, int gradIndex, int numGrads)
{
	if (SHVs == 0) {
		SHVs = new Matrix(2,numGrads);
	}


	// Initialize unconditaional stress sensitivity
	double gradient = 0.0;


	// Pick up sensitivity history variables
	double CstrainSensitivity = 0.0;
	double CstressSensitivity	 = 0.0;
	if (SHVs != 0) {
		CstrainSensitivity = (*SHVs)(0,gradIndex);
		CstressSensitivity = (*SHVs)(1,gradIndex);
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


	// Compute min and max stress
	double Tstress;
	double dStrain = Tstrain-Cstrain;
	double sigmaElastic = Cstress + E0*dStrain;
	double fyOneMinusB = fy * (1.0 - b);
	double Esh = b*E0;
	double c1 = Esh*Tstrain;
	double c2 = TshiftN*fyOneMinusB;
	double c3 = TshiftP*fyOneMinusB;
	double sigmaMax = c1+c3;
	double sigmaMin = c1-c2;


	// Evaluate stress sensitivity ('gradient')
	if ( (sigmaMax < sigmaElastic) && (fabs(sigmaMax-sigmaElastic)>1e-5) ) {
		Tstress = sigmaMax;
		gradient = E0Sensitivity*b*Tstrain 
				 + E0*bSensitivity*Tstrain
				 + E0*b*TstrainSensitivity
				 + TshiftP*(fySensitivity*(1-b)-fy*bSensitivity);
	}
	else {
		Tstress = sigmaElastic;
		gradient = CstressSensitivity 
			     + E0Sensitivity*(Tstrain-Cstrain)
				 + E0*(TstrainSensitivity-CstrainSensitivity);
	}
	if (sigmaMin > Tstress) {
		gradient = E0Sensitivity*b*Tstrain
			     + E0*bSensitivity*Tstrain
			     + E0*b*TstrainSensitivity
				 - TshiftN*(fySensitivity*(1-b)-fy*bSensitivity);
	}


	// Commit history variables
	(*SHVs)(0,gradIndex) = TstrainSensitivity;
	(*SHVs)(1,gradIndex) = gradient;

	return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////


//this function is no use, just for the definiation of pure virtual function.
int Steel01PThermal::setTrialStrain (double strain, double strainRate)
{
  opserr << "Steel01PThermal::setTrialStrain (double strain, double strainRate) - should never be called\n";
  return 0;
}


int 
Steel01PThermal::getVariable(const char *variable, Information &info)
{
  if (strcmp(variable,"ThermalElongation") == 0) {
    info.theDouble = ThermalElongation;    
    return 0;
  } else if (strcmp(variable,"ElongTangent") == 0) {
    Vector *theVector = info.theVector;
    if (theVector != 0) {
      double tempT, ET, Elong, TempTmax;
      tempT = (*theVector)(0);
	  ET = (*theVector)(1);
	  Elong = (*theVector)(2);
      TempTmax = (*theVector)(3);
      this->getElongTangent(tempT, ET, Elong, TempTmax);
	  (*theVector)(0) = tempT;
      (*theVector)(1) = ET;
      (*theVector)(2) = Elong;
	  (*theVector)(3) = TempTmax;
    }
    return 0;
  }
  return -1;
}


