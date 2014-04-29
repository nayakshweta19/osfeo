#ifndef MCFTConcrete01_h
#define MCFTConcrete01_h

// File MCFTConcrete01.h
// Written: Li Ning (Neallee@tju.edu.cn)
// Created: 2011.7

#include <UniaxialMaterial.h>

class MCFTConcrete01 : public UniaxialMaterial
{
  public:
  MCFTConcrete01(int tag, double fpc, double eco);
  MCFTConcrete01();
  ~MCFTConcrete01();
  
  int setTrialStrain(double strain, double strainRate = 0.0); 	
  int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
  
  double getStrain(void);              
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void) {return 2.0*fpc/epsc0;};

  double getSecant(void);
  double getZeta(void);
  double getPD(void);

  int commitState(void);
  int revertToLastCommit(void);    
  int revertToStart(void);   
  
  UniaxialMaterial *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    

  Response *setResponse (const char **argv, int argc, 
			 OPS_Stream &theOutputStream);
  int getResponse (int responseID, Information &matInformation);   
  void Print(OPS_Stream &s, int flag =0);
  
 protected:
  
 private:
  
  /*** Material Properties ***/        
  double fpc;      // Compressive strength
  double epsc0;    // Strain at compressive strength
  double epscpCP;  // pre compression epscp
  double epscpTP;  // pre tension epscp
  double epscp;    // recent epscp
  double Ec0;
  double zeta;     // Softening effect
  double itap;     // Parameter needed for calculating zeta
  double epslonTP; // Strain in the perpendicular direction, needed to get the zeta
  double D;        // Damage factor for strength, get from parameter
  double X;
  double K;
  bool   tensionOccur;

  double epsCM;    // maximum compression strain before reverse occur 
  double epsTM;    // maximum tension strain before reverse occur 
  double epsTM1;   // considering offset of tension strain
  double sigCM;
  double sigTM;

  /*** History Variables ***/
  int CloadingState; // Flag for loading state
  // 1 = ascending branch of envelope in compression
  // 2 = descending branch of envelope in compression
  // 3 = ascending branch of envelope in tension
  // 4 = descending branch of envelope in tension
  // 5 = unloading at compression envelope 1 or 2
  // 6 = unloading from tension to compression zone ( from 4 to 1 or 2 )
  // 7 = reloading from compression to tension (reload from 5 approach to tension zone)
  // 8 = reloading from tension to compression (reload from 6 approach to compression zone)
  // 
  
  int reloadPath;   // Flag for reversing state
  // = 0 initially
  // = 1 if reverse from compression to ..
  // = 2 if reverse from tension to..
  
  int unloadPath;   // Flag for reversing state
  // = 0 initially
  // = 1 if reverse from ascending compression branch
  // = 2 if reverse from descending compression branch

  double reverseFromOneStrain; // Strain and stress of the reversed point from path 1
  double reverseFromOneStress;
  double reverseFromTwoStrain; // Strain and stress of the reversed point from path 2
  double reverseFromTwoStress;
  double reverseFromFourStrain; // Strain and stress of the reversed point from path 4
  double reverseFromFourStress;
  
  /*** Current Key Points According To State Variables ***/
  double interFiveSevenStrain;  // Strain of start point of path 7 (end point of path 5)
  
  double approachFiveToComStrain;    // Intersection point of path 5 to compressive envelope
  double approachSixToComStrain;     // Intersection point of path 6 to compressive envelope
  
  
  /*** CONVERGED State Variables ***/    
  double Cstrain;
  double Cstress;
  double Ctangent;    
  
  /*** TRIAL History Variables ***/
  int TloadingState;	
  
  /*** TRIAL State Variables ***/
  double Tstrain;
  double Tstress;
  double Ttangent; // Not really a state variable, but declared here
  // for convenience
  
  
  // Calculates the trial state variables based on the trial strain
  void determineTrialState(double dStrain);
  
  // Calculate the plastic strain offsets
  void determineCompEpscp(double eps);
  void determineTensEpscp(double eps);

  //Envelope curve in beginning
  void envelope( );
  void envelope(double offsetEsp);
  void inline updateCT();

  void pathFive( );
  void pathSix ( );
  void pathSeven();
  void pathEight();

};

#endif
