#ifndef MCFTSteel01_h
#define MCFTSteel01_h

// MCFTSteel01.h
// Written: Li Ning (Neallee@tju.edu.cn)
// Description: This file contains the class definition for MCFTSteel01.h

#include <UniaxialMaterial.h>

// Default values for reloading/unloading path parameters AC, RC
#define STEEL_Z01_DEFAULT_AC        1.9
#define STEEL_Z01_DEFAULT_RC       10.0
#define LOOP_NUM_LIMIT               30
#define SIZE                         30 //limit of array number

class MCFTSteel01 : public UniaxialMaterial
{
  public:
    MCFTSteel01(int tag,
	    double fy, double E0, double b,
	    double R0, double cR1, double cR2,
	    double a1, double a2, double a3, double a4, double sigInit =0.0);
    
    // Constructor for no isotropic hardening
    MCFTSteel01(int tag,
	    double fy, double E0, double b,
	    double R0, double cR1, double cR2);
    
    // Constructor for no isotropic hardening
    // Also provides default values for R0, cR1, and cR2
    MCFTSteel01(int tag, double fy, double E0, double b);
	    
    MCFTSteel01(void);
    ~MCFTSteel01();

    const char *getClassType(void) const {return "MCFTSteel01";};

    double getInitialTangent(void);
    UniaxialMaterial *getCopy(void);

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);      
    double getStress(void);
    double getTangent(void);
    double getSecant(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

	Response *setResponse (const char **argv, int argc, 
		OPS_Stream &theOutputStream);
	int getResponse (int responseID, Information &matInformation);          

    void Print(OPS_Stream &s, int flag =0);

  protected:
    
  private:
    // matpar : STEEL FIXED PROPERTIES
    double Fy;  //  = matpar(1)  : yield stress
    double E0;  //  = matpar(2)  : initial stiffness
    double b;   //  = matpar(3)  : hardening ratio (Esh/E0)
    double R0;  //  = matpar(4)  : exp transition elastic-plastic
    double cR1; //  = matpar(5)  : coefficient for changing R0 to R
    double cR2; //  = matpar(6)  : coefficient for changing R0 to R
    double a1;  //  = matpar(7)  : coefficient for isotropic hardening in compression
    double a2;  //  = matpar(8)  : coefficient for isotropic hardening in compression
    double a3;  //  = matpar(9)  : coefficient for isotropic hardening in tension
    double a4;  //  = matpar(10) : coefficient for isotropic hardening in tension
    double sigini; // initial 
    // hstvP : STEEL HISTORY VARIABLES
    double epsminP; //  = hstvP(1) : max eps in compression
    double epsmaxP; //  = hstvP(2) : max eps in tension
    double epsplP;  //  = hstvP(3) : plastic excursion
    double epss0P;  //  = hstvP(4) : eps at asymptotes intersection
    double sigs0P;  //  = hstvP(5) : sig at asymptotes intersection
    double epssrP;  //  = hstvP(6) : eps at last inversion point
    double sigsrP;  //  = hstvP(7) : sig at last inversion point
    int    konP;    //  = hstvP(8) : index for loading/unloading
    // hstv : STEEL HISTORY VARIABLES   
    double epsP;  //  = strain at previous converged step
    double sigP;  //  = stress at previous converged step
    double eP;    //   stiffness modulus at last converged step;

    double epsmin; 
    double epsmax; 
    double epspl;  
    double epss0;  
    double sigs0; 
    double epsr;  
    double sigr;  
    int    kon;    
    double sig;   
    double e;     
    double eps;   //  = strain at current step

  
    void determineTrialState (double dStrain);  
	// Calculates the trial state variables based on the trial strain

    void initialEnvelope( ); 
	void tensionEnvelope( );
	void compressionEnvelope( );
	
	void determineTrialLoop(double dStrain); // Calculate the trail state variables in the loop
	void determineDownPathPoint( ); //determine key points of down and up path of hysterics loop
	void determineUpPathPoint( );

	void downPath( ); 
	void upPath( );

	void reverseFromTenEnvelope( );
	void reverseFromComEnvelope( );

	void reverseLoopSetZero( ); // Loop variables set to zero when merge into envelope
	
	double tt1; // for check
	double tt2;
	double ttStrain;

};

#endif
