// File: MCFTConcrete01.cpp
// Written: Li Ning (Neallee@tju.edu.cn)
// Created: 2011.7
// Description: This file contains the class implementation for
// MCFTConcrete01

#include "MCFTConcrete01.h"
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include <MaterialResponse.h>
#include <elementAPI.h>

using namespace std;

#define OPS_Export 

OPS_Export void *
OPS_NewMCFTConcrete01Material(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[2];
  int numData = 1;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs < 3) {
    opserr << "Want: uniaxialMaterial MCFTConcrete01 tag? fpc? epsc0?" << endln;
    return 0;	
  }

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial MCFTConcrete01 tag" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial MCFTConcrete01 tag? fpc? epsc0?" << endln;
    return 0;	
  }

  // Parsing was successful, allocate the material
  theMaterial = new MCFTConcrete01(iData[0], dData[0], dData[1]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type MCFTConcrete01\n";
    return 0;
  }

  return theMaterial;
}

MCFTConcrete01::MCFTConcrete01
(int tag, double FPC, double EPSC0): 
UniaxialMaterial(tag, MAT_TAG_MCFTConcrete01), fpc(FPC), epsc0(EPSC0), D(1.0),
	epsCM(0.0), sigCM(0.0), epsTM(0.0), sigTM(0.0), epsTM1(0.0), tensionOccur(0)
{
  // Make all concrete parameters negative
  if (fpc > 0.0)
    fpc = -fpc;
  
  if (epsc0 > 0.0)
    epsc0 = -epsc0;
  
  // Set trial values
  this->revertToStart();
 }

MCFTConcrete01::MCFTConcrete01():
UniaxialMaterial(0, MAT_TAG_MCFTConcrete01), fpc(0.0), epsc0(0.0), D(1.0),
	epsCM(0.0), sigCM(0.0), epsTM(0.0), sigTM(0.0), epsTM1(0.0), tensionOccur(0)
{
	// Set trial values
	this->revertToStart();
 }

MCFTConcrete01::~MCFTConcrete01 ()
{
   // Does nothing
}

int MCFTConcrete01::setTrialStrain (double strain, double strainRate)
{
  if (epslonTP > 0.0) {
    // add K into zeta, K is delta
    zeta =  (K) * 5.8 / sqrt( -fpc * ( 1.0 + 400.0 * epslonTP / itap ) ); 
    if ( zeta >= 0.9 )
      zeta = 0.9;
    if ( zeta <= 0.25 ) //min zeta
      zeta = 0.25;
  } else {
    zeta = 1.0;
  }
 
  // Reset history variables to last converged state
  TloadingState = CloadingState;  

  // Set trial strain
  Tstrain = strain;

  // Determine change in strain from last converged state
  double dStrain = Tstrain - Cstrain;

  if ( Tstrain < -0.02 || Tstrain > 0.0025) {
	Tstress = 0.0;
	Ttangent = 10e-10;
	return 0;
  }

  // Calculate the trial state given the trial strain
  if (fabs(dStrain) > DBL_EPSILON)   
    determineTrialState (dStrain); 
  
  return 0;
}

int MCFTConcrete01::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
    return 0;
}

double MCFTConcrete01::getStrain ()
{
   return Tstrain;
}

double MCFTConcrete01::getStress ()
{
   return Tstress;
}

double MCFTConcrete01::getTangent ()
{
   return Ttangent;
}

double MCFTConcrete01::getSecant ()
{
	if ( Tstrain == 0.0 )
	{
		return Ec0;
	}
    else
	{
		return Tstress/Tstrain;
	}
}

double MCFTConcrete01::getZeta ()
{
	return zeta;
}

int MCFTConcrete01::commitState ()
{
   // History variables
   CloadingState = TloadingState;

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int MCFTConcrete01::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TloadingState = CloadingState;

   // Reset trial state variables to last committed state
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}

int MCFTConcrete01::revertToStart ()
{
   // History variables
   
   TloadingState = 0;
   CloadingState = 0;
   
   reloadPath = 0;
   unloadPath = 0;

   reverseFromOneStrain = 0.0;
   reverseFromOneStress = 0.0;
   reverseFromTwoStrain = 0.0; 
   reverseFromTwoStress = 0.0;
   reverseFromFourStrain = 0.0; 
   reverseFromFourStress = 0.0;   

   interFiveSevenStrain = 0.0;

   approachFiveToComStrain = 0.0;
   approachSixToComStrain = 0.0;
 
   // State variables
   zeta     = 1.0;
   itap     = 1.0;
   epslonTP = 0.0;
   Ec0      = 2.0*fpc/epsc0;
   D        = 1.0;
   epscp    = 0.0;

   Cstrain  = 0.0;
   Cstress  = 0.0;
   Ctangent = Ec0;

   Tstrain  = 0.0;
   Tstress  = 0.0;
   Ttangent = Ec0;

   return 0;
}

void MCFTConcrete01::determineTrialState (double dStrain)
{
	switch (TloadingState)
	{
	case 0:
		envelope();
		break;

	case 1:
		if ( dStrain < 0 ) //Continues on envelope 1
		{
		    envelope();
		}
		else //Reverse occurs at ascending 1 branch, go to path 5
		{
		    updateCT();
			reverseFromOneStrain = Cstrain;
			reverseFromOneStress = Cstress;

			TloadingState = 5;
			unloadPath = 1; // reverse from ascending compression branch
            if (epsTM > 0.0) reloadPath = 2; // reverse from tension to ...
			else reloadPath = 1; // reverse from compression to ...

        	pathFive();

			if ( Tstress > 0 ) {// Reach path 7, further break through zero stress from path 5
				TloadingState = 7;
				pathSeven();
			}
			epscpCP = epscp;
		}
		opserr << "TloadState = " << TloadingState << endln;
		break;
	
	case 2:
		updateCT();
		
		if ( dStrain < 0 ) { //Continues on envelope 2
			envelope();
		} else { // reverse from branch 2 to path 5
			
			reverseFromTwoStrain = Cstrain;
			reverseFromTwoStress = Cstress;

			TloadingState = 5;
			unloadPath = 2; // reverse from descending compression branch
			if (epsTM > 0.0) reloadPath = 2; // reverse from tension to ...
			else reloadPath = 1; // reverse from compression to ...

			//interFiveSevenStrain = reverseFromTwoStrain - reverseFromTwoStress/(0.8*Ec0);
            //getApproachFiveToComStrain ( )
			pathFive();

			if ( Tstress > 0 ) { // Reach path 7, further break through zero stress from path 5
				TloadingState = 7;
				pathSeven();
			}
			epscpCP = epscp;
		}
		opserr << "TloadState = " << TloadingState << endln;
		break;

	case 3:
	    //determineCompEpscp(epsTM);
	    envelope();
        updateCT();
	    opserr << "TloadState = " << TloadingState << endln;
	    break;

	case 4 :
		if ( dStrain >= 0) { //Continues on envelope
			envelope ();
		} else {
			updateCT();

			reverseFromFourStrain = Cstrain;
			reverseFromFourStress = Cstress;

			TloadingState = 6;

			unloadPath = 2;   //reverse from descending tensional branch
			reloadPath = 2;   //reverse from tension to..

			pathSix();

			epscpTP = epscp;
		}
		opserr << "TloadState = " << TloadingState << endln;
		break;

	case 5 :
		//getApproachFiveToComStrain( );
        updateCT();

		if (dStrain >= 0.0) {
		  pathFive();
		} else {
		  TloadingState = 8;

		  reverseFromOneStrain = Cstrain;
		  reverseFromOneStress = Cstress;
		  reverseFromTwoStrain = Cstrain;
		  reverseFromTwoStress = Cstress;

		  pathEight();
		}
		
		if ( Tstress > -1.0*DBL_EPSILON ) {  //Reach path 7
		  TloadingState = 7;
		  pathSeven();
		}
		break;

	case 6:
	  //reverseFromFourStrain = Cstrain;
	  //reverseFromFourStress = Cstress;
	  updateCT();
	  reloadPath = 2;
	  //  getApproachSixToComStrain ( );
	  if (dStrain <= 0.0) // unloading on path 6
	  {
	    if ( Tstrain < 0 ) {  //Reach path 8
	      TloadingState = 8;
	      pathEight();
	  
	    } else { // going on the unloading path 6
	      pathSix();
	  
	    }
	  }

	  else   //reloading from path 6 to path 7
	  {
	    TloadingState = 7;
	    pathSeven();
	  }
	  
	  if (Tstress < 0.0) {
	  	TloadingState = 8;
	    pathEight();
	  }
	  break;

	case 7:
	  updateCT();

	  reverseFromFourStrain = Cstrain;
	  reverseFromFourStress = Cstress;

	  if (dStrain >= 0.0) { //continue tension
	    pathSeven();
	    if (Tstrain >= epsTM + DBL_EPSILON) { // tension to descending branch
	      TloadingState = 4;
	      envelope(epscpTP); //
		  reloadPath = 2;
	    }
	  }
	  else {             //reverse to compression
		epsTM = -epscp;
	    if (Tstrain >=0.0) {
	      TloadingState = 6;
	      pathSix();
	    } else {
	      TloadingState = 8;
	      pathEight();
	    }
	  }
	  break; // if TloadingState ==7
	
	case 8:
	    updateCT();

		if (dStrain < 0.0) {
			if ( Tstrain < epsCM ){
				TloadingState = 1;
				envelope();
			} else {
				pathEight();
			}
		} else {
			if (Tstrain >= 0.0) {
				TloadingState = 7;
			    pathSeven();
			} else { 
				TloadingState = 5;
				pathFive();
			}
		}
		break; // if TloadingState ==8
	
	default:
		opserr << " MCFTConcrete01::determineTrialState -- impropter TloadingState: " 
		       << TloadingState << "\n";
		return;
	}

} // end trialState

void inline
MCFTConcrete01::updateCT()
{
  if (Cstrain < epsCM) {
	epsCM = Cstrain;
	sigCM = Cstress;
  }
  if (Cstrain > epsTM) {
    epsTM = Cstrain;
    sigTM = Cstress;
  }
}

void
MCFTConcrete01::determineCompEpscp(double eps)
{
  epscp = eps - zeta * epsc0 * ( 0.868*(eps/zeta/epsc0)
	    - 0.166 * pow( eps/zeta/epsc0, 2.0) );
}

void
MCFTConcrete01::determineTensEpscp(double eps)
{
  epscp = 146.0 * pow(eps, 2.0) + 0.523 * eps;
}

void
MCFTConcrete01::pathFive()
{
  
  double slop2 = Ec0;
  double slop3 = 0.071 * Ec0;
  
  if ( unloadPath == 1 )
  {
    determineCompEpscp(reverseFromOneStrain);

	if (reverseFromOneStrain <= 0.62*epsc0) // need further revision ?!?!
	{
      double N = (slop2-slop3)*(epscp-reverseFromOneStrain)
	           / (reverseFromOneStress+slop2*(epscp-reverseFromOneStrain));
	  
      double dStrain = Tstrain - reverseFromOneStrain;
	  
      Tstress = reverseFromOneStress+slop2*dStrain+(slop3-slop2)*pow(dStrain, N)
	  	    / N / pow(epscp-reverseFromOneStrain, N-1);
	  
      Ttangent = slop2 + (slop3-slop2)*pow(dStrain/(epscp-reverseFromOneStrain), N-1);
	
	} else {

      Ttangent = reverseFromOneStress/(reverseFromOneStrain-epscp);
	  Tstress = Ttangent*(Tstrain-epscp);

	}
  } // if ( unloadPath == 1 )

  else if ( unloadPath == 2 )
  {
    slop2 *= 0.8;
	slop3 *= 0.8;
	determineCompEpscp(reverseFromTwoStrain);

    double N = (slop2-slop3)*(epscp-reverseFromTwoStrain)
	         / (reverseFromTwoStress+slop2*(epscp-reverseFromTwoStrain));

    double dStrain = Tstrain - reverseFromTwoStrain;

    Tstress = reverseFromTwoStress+slop2*dStrain+(slop3-slop2)*pow(dStrain, N)
	        / N / pow(epscp-reverseFromTwoStrain, N-1);

    Ttangent = slop2 + (slop3-slop2)*pow(dStrain/(epscp-reverseFromTwoStrain), N-1);

  } // if ( unloadPath == 2 )
  
  else
  {
  	opserr << " MCFTConcrete01::pathFive -- improper unloadPath : " 
  		 << unloadPath << "\n";
  }
  //opserr << "TloadingState = " << TloadingState << ". Path 5." << endln;
}

void
MCFTConcrete01::pathSix()
{  
  double slop5 = Ec0;
  double slop6 = 0.071 * Ec0*(0.001/reverseFromFourStrain);

  if (reverseFromFourStrain> 0.001) slop6 *= 0.75;

  determineTensEpscp(reverseFromFourStrain);

  double N = (slop5-slop6)*(reverseFromFourStrain-epscp)
	         / (slop5*(reverseFromFourStrain-epscp)-reverseFromFourStress);

  double dStrain = reverseFromFourStrain-Tstrain;

  Tstress = reverseFromFourStress-slop5*dStrain+(slop5-slop6)*pow(dStrain, N)
		    / N / pow(reverseFromFourStrain-epscp, N-1);

  if (epscp >= reverseFromFourStrain)
    Ttangent = 0.5 * (slop5 + slop6);
  else 
    Ttangent = slop5+(slop5-slop6)*pow(dStrain/(reverseFromFourStrain-epscp), N-1);
	 
  //opserr << "TloadingState = " << TloadingState << ". Path 6." << endln;
}

void
MCFTConcrete01::pathSeven()
{
    // from compression to tension region

	double ecr = 0.00008; //crack strain
	//double fcr = 0.31 * sqrt( -fpc ); //crack strength

	determineCompEpscp(epsCM);

    // Approach to reverse point at descending branch of tensile envelope
    if ( epsTM > ecr ) // has tension to crack
	{ 
	  double slope1 = reverseFromFourStress / ( reverseFromFourStrain - epscp );
	  Tstress = slope1 * ( Tstrain - epscp );
	  Ttangent = slope1;
	}
	// Approach to peak point of the tensile envelope
	else if ( epsTM < DBL_EPSILON ) // epsTM = 0 , no tension occur
	{
	    envelope(epscpCP);
	} 
	
	else // tension not to crack, path 3 with epscp
	{
	  //envelope(epscp);
	  Tstress = 0.0;
	  Ttangent = 10e-10;
	}
	//opserr << "TloadingState = " << TloadingState << ". Path 7." << endln;
}

void
MCFTConcrete01::pathEight()
{
	//double ecr = 0.00008; //crack strain
	//double fcr = 0.31 * sqrt( -fpc ); //crack strength

    determineCompEpscp(epsCM);
	//  // reverse from tension region
	if (reverseFromFourStrain >= 0.0 ) {
	  if (Tstrain >= 0.0) {
	    Tstress = 0.0;
	    Ttangent = 10e-10;
	  } else {
	    if (reloadPath == 2) { // reverse from tension to..
	    	double slope1 = sigCM / epsCM;
	    	Tstress = slope1 * Tstrain;
	    	Ttangent = slope1;
	    } else {               // reverse from compression to ..
	 	    double slope1 = sigCM / (epsCM-epscp);
	 	    Tstress = slope1 * Tstrain;
	 	    Ttangent = slope1;
	    }
	  }
	}
	// Approach to peak point of the compression envelope
	else // tension region has not been achieved ever
	{
	    if (unloadPath != 2) {
			double slope1 = reverseFromOneStress / ( reverseFromOneStrain - epscp );
		    Tstress = slope1 * ( Tstrain - epscp);
		    Ttangent = slope1;
		} else {
			double slope1 = reverseFromTwoStress / ( reverseFromTwoStrain - epscp );
			Tstress = slope1 * ( Tstrain - epscp);
			Ttangent = slope1;
		}
	}

	//opserr << "TloadingState = " << TloadingState << ". Path 8." << endln;
}

void
MCFTConcrete01::envelope( )
{
	//double fcr = 0.5167 * sqrt( -fpc ); //crack strength for hament's test
	double fcr = 0.31 * sqrt( -fpc ); //crack strength
	double ecr = 0.00008; //crack strain
	
	if ( Tstrain >= 0 ) // Tension
	{
		if ( Tstrain <= ecr ) //Ascending branch
		{
			//double Ec = 6458.0 * sqrt( -fpc ); // for hament's test
			//double Ec = 3875.0 * sqrt( -fpc );
			Tstress = Ec0 * Tstrain;
			Ttangent = Ec0;
			TloadingState = 3;
		}
		else //Descending branch
		{
			Tstress = fcr * pow((ecr/Tstrain), 0.4);
			Ttangent = - fcr * 0.4 * pow(ecr, 0.4) * pow(Tstrain, -1.4);
			TloadingState = 4;
		}
	} 
	else //Compression
	{
		if ( Tstrain >= zeta*epsc0 ) //Ascending branch
        {	
			TloadingState = 1;
			double eta = Tstrain/(zeta*epsc0);
			Tstress = D * zeta * fpc * (2*eta-eta*eta);
			Ttangent = D*Ec0*(1.0-eta);

		}
		else //Descending branch
		{
			TloadingState = 2;
			double temp = (Tstrain/(zeta*epsc0)-1.0)/(4.0/zeta-1.0);
			
			// initially descending branch
			Tstress = D * zeta * fpc * ( 1.0 - temp*temp );
			Ttangent = - D*2.0*fpc*temp/epsc0/(4.0/zeta-1.0);

			// change by using X as powder, 2004/08, X is n
            // Tstress = D * zeta * fpc * ( 1.0 - pow(temp, X) );
			// Ttangent = - D*fpc*(X)*pow(temp,-1.0+X)/epsc0/(4.0/zeta-1.0);

			//check if Tstress > 0.2*zeta*fpc, change to platum part
			if ( Tstress > D*0.2*zeta*fpc )
			{
				Tstress = D*0.2*zeta*fpc;
				Ttangent = 1.0e-10;
			}

		}
	}
}

void
MCFTConcrete01::envelope(double offsetEsp)
{
	//double fcr = 0.5167 * sqrt( -fpc ); //crack strength for hament's test
	double fcr = 0.31 * sqrt( -fpc ); //crack strength
	double ecr = 0.00008; //crack strain
	
	if ( Tstrain >= offsetEsp ) // Tension
	{
		if ( Tstrain <= offsetEsp+ecr ) //Ascending branch
		{
			//double Ec = 6458.0 * sqrt( -fpc ); // for hament's test
			//double Ec = 3875.0 * sqrt( -fpc );
			Tstress = Ec0 * (Tstrain-offsetEsp);
			Ttangent = Ec0;
			//TloadingState = 7;
		}
		else //Descending branch
		{
			Tstress = fcr * pow((ecr/(Tstrain-offsetEsp)), 0.4);
			Ttangent = - fcr * 0.4 * pow(ecr, 0.4) * pow(Tstrain-offsetEsp, -1.4);
			//TloadingState = 7;
		}
	} 
	else //Compression
	{
		//TloadingState = 7;
			
		// initially zero
		Tstress = 0.0;
		Ttangent = 10e-10;

		// change by using X as powder, 2004/08, X is n
        // Tstress = D * zeta * fpc * ( 1.0 - pow(temp, X) );
		// Ttangent = - D*fpc*(X)*pow(temp,-1.0+X)/epsc0/(4.0/zeta-1.0);

		//check if Tstress > 0.2*zeta*fpc, change to platum part
		//if ( Tstress > D*0.2*zeta*fpc )
		//{
		//	Tstress = D*0.2*zeta*fpc;
		//	Ttangent = 0.0;
		//}
	}
}

UniaxialMaterial*
MCFTConcrete01::getCopy ()
{
   MCFTConcrete01* theCopy = new MCFTConcrete01(this->getTag(), fpc, epsc0);  
   
   // History variables
   theCopy->TloadingState = TloadingState;
   theCopy->CloadingState = CloadingState;
   
   theCopy->reloadPath = reloadPath;
   theCopy->unloadPath = unloadPath;

   theCopy->reverseFromOneStrain = reverseFromOneStrain;
   theCopy->reverseFromOneStress = reverseFromOneStress;
   theCopy->reverseFromTwoStrain = reverseFromTwoStrain; 
   theCopy->reverseFromTwoStress = reverseFromTwoStress;
   theCopy->reverseFromFourStrain = reverseFromFourStrain; 
   theCopy->reverseFromFourStress = reverseFromFourStress;   

   theCopy->interFiveSevenStrain = interFiveSevenStrain;

   theCopy->approachFiveToComStrain = approachFiveToComStrain;
   theCopy->approachSixToComStrain = approachSixToComStrain;
 
   // State variables
   theCopy->zeta = zeta;
   theCopy->itap = itap;
   theCopy->epslonTP = epslonTP;
   theCopy->Ec0     = Ec0;
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   theCopy->Tstrain = Tstrain;
   theCopy->Tstress = Tstress;
   theCopy->Ttangent = Ttangent;

   theCopy->D = D;
   theCopy->epscp =epscp;

   return theCopy;
}

int
MCFTConcrete01::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(24);
   data(0) = this->getTag();

   // Material properties
   data(1) = fpc;
   data(2) = epsc0; 
   data(3) = zeta;
   data(4) = itap;
   data(5) = epslonTP;

   // History variables from last converged state
   data(6)  = CloadingState;
   data(7)  = reloadPath;
   data(8)  = unloadPath;
   data(9)  = reverseFromOneStrain;
   data(10) = reverseFromOneStress;
   data(11) = reverseFromTwoStrain; 
   data(12) = reverseFromTwoStress;
   data(13) = reverseFromFourStrain; 
   data(14) = reverseFromFourStress;
   data(15) = interFiveSevenStrain;
   data(16) = approachFiveToComStrain;
   data(17) = approachSixToComStrain;

   // State variables from last converged state
   data(18) = Cstrain;
   data(19) = Cstress;
   data(20) = Ctangent;

   data(21) = D;
   data(22) = Ec0;
   data(23) = epscp;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "MCFTConcrete01::sendSelf() - failed to send data\n";

   return res;
}

int
MCFTConcrete01::recvSelf (int commitTag, Channel& theChannel,
                                 FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(24);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);

   if (res < 0) {
      opserr << "MCFTConcrete01::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties 
      fpc = data(1);
      epsc0 = data(2); 
      zeta = data(3);
      itap = data(4);
      epslonTP = data(5);

      // History variables from last converged state
	  CloadingState = int(data(6));
      reloadPath = int(data(7));
	  unloadPath = int(data(8));
	  reverseFromOneStrain = data(9);
      reverseFromOneStress = data(10);
      reverseFromTwoStrain = data(11); 
      reverseFromTwoStress = data(12);
      reverseFromFourStrain = data(13); 
      reverseFromFourStress = data(14);
      interFiveSevenStrain = data(15);
      approachFiveToComStrain = data(16);
      approachSixToComStrain = data(17);

	  // Copy converged history values into trial values since data is only
      // sent (received) after convergence	  
      TloadingState = CloadingState;

      // State variables from last converged state
      Cstrain = data(18);
      Cstress = data(19);
      Ctangent = data(20);

	  D = data(21);
	  Ec0 = data(22);
	  epscp = data(23);

      // Set trial state variables
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }

   return res;
}

Response* 
MCFTConcrete01::setResponse(const char **argv, int argc,
			 OPS_Stream &theOutput)
{
  Response *theResponse = 0;

  if (strcmp(argv[0],"getPD") == 0) {
    double data = 0.0;
    theResponse = new MaterialResponse(this, 100, data);
  } else if (strcmp(argv[0],"setWallVar") == 0) {
    theResponse = new MaterialResponse(this, 101, Vector(6));
  } else
    return this->UniaxialMaterial::setResponse(argv, argc, theOutput);

  return theResponse;
}
 
int 
MCFTConcrete01::getResponse(int responseID, Information &matInfo)
{
  if (responseID == 100) {
    matInfo.theDouble = this->getPD();
  } else if (responseID == 101){
    Vector *theVector = matInfo.theVector;
	(*theVector)(0) = X;
	(*theVector)(1) = K;
	(*theVector)(2) = D;
	(*theVector)(3) = itap;
	(*theVector)(4) = epslonTP;
	(*theVector)(5) = epscp;
  } else
    return this->UniaxialMaterial::getResponse(responseID, matInfo);

  return 0;
}

double
MCFTConcrete01::getPD ()
{
	double PD;
	double tempRatio;

	if ( epslonTP <= 0.0 )
		PD = 0.0;
	else {
		if ( TloadingState == 1 ) {
			// Compression Ascending branch
			tempRatio = Tstrain/(zeta*epsc0);

			PD = - D * 1160.0 * sqrt(-fpc)/itap * pow((1+400.0*epslonTP/itap), -1.5)
				 * pow(tempRatio,2.0); // 1160, 400

		} else if ( TloadingState == 2 ) { 
		    // Compression Descending branch
			if ( Ttangent == 0.0 ) {
				// at the end platum part of descending branch FMK CHANGED FROM = 0.0
				PD = 0.0;
			} else {
				tempRatio = Tstrain/(zeta*epsc0);
				PD = - D * 1160.0 * sqrt(-fpc)/itap * pow((1+400.0*epslonTP/itap), -1.5) 
 	  	           * (1.0 -(tempRatio-1)/pow(4.0/zeta-1.0,3.0)*(1.0-12.0/zeta+(4.0/zeta+1.0)*tempRatio)); 
				// 1160, 400
			}
		} else {
			PD = 0.0;
		}
		if ( (zeta == 0.9) || (zeta == 0.25)) // zeta = max or min value
			PD = 0.0;
	}

    return PD;
}

void 
MCFTConcrete01::Print (OPS_Stream& s, int flag)
{
   s << "MCFTConcrete01, tag: " << this->getTag() << endln;
   //s << "  fpc: " << fpc << endln;
   //s << "  epsc0: " << epsc0 << endln;

   s << " strain: " << this->getStrain() << endln;
   s << " stress: " << this->getStress() << endln;
   s << " tangent: " << this->getTangent() << endln;
   //s << " PD: " << this->getPD() << endln;
   s << " zeta: " << zeta << endln;
   s << " D: " << D << endln;
   s << " TloadingState: " << TloadingState << endln;
   s << " reverseFromFourStrain: " << reverseFromFourStrain << endln;   

}


