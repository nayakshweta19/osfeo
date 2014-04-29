#ifndef MCFTSteel02_h
#define MCFTSteel02_h

// MCFTSteel02.h
// Written: Li Ning (Neallee@tju.edu.cn)
// Created: 2003.7
// Description: This file contains the class definition for MCFTSteel02.h

#include <UniaxialMaterial.h>

class MCFTSteel02 : public UniaxialMaterial
{
  public:
    MCFTSteel02(int tag,
	    double fy, double E0, double b,
	    double r, double cR1, double cR2, double a1, double a2);
	    
    MCFTSteel02(void);
    ~MCFTSteel02();

    const char *getClassType(void) const {return "MCFTSteel02";};

	double getInitialTangent(void) {return E0;};
    UniaxialMaterial *getCopy(void);

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
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

// AddingSensitivity:BEGIN //////////////////////////////////////////
	int    setParameter             (const char **argv, int argc, Parameter &param);
	int    updateParameter          (int parameterID, Information &info);
	int    activateParameter        (int parameterID);
	double getStressSensitivity     (int gradNumber, bool conditional);
	double getStrainSensitivity		(int gradNumber);
	double getInitialTangentSensitivity(int gradNumber);
	int    commitSensitivity        (double strainGradient, int gradNumber, int numGrads);
// AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    /*** Material Properties ***/
    double fy;  // Yield stress
    double E0;  // stiffness
    double b;   // Hardening ratio (b = Esh/E0)
    double r0;	// radius of rounded corners  20
    
	double coeffR1;  //18.5
    double coeffR2;  //0.15
    double a1;            //0
    double a2;  // 0

    /*** CONVERGED History Variables ***/
    double CminStrain;  // Minimum strain in compression
    double CmaxStrain;  // Maximum strain in tension

    int Cloading;       // Flag for loading/unloading
                        // 1 = loading (positive strain increment)
                        // -1 = unloading (negative strain increment)
                        // 0 initially
    double CYieldStrain;
    double CYieldStress;
    double CReverStrain;
    double CReverStress;
    double CPlasticExcursion;

    /*** CONVERGED State Variables ***/    
    double Cstrain;
    double Cstress;
    double Ctangent;

    /*** TRIAL History Variables ***/
    double TminStrain;   // abs of minimum strain
    double TmaxStrain;   // abs of minimum strain

    int Tloading;
    double TYieldStrain;
    double TYieldStress;
    double TReverStrain;
    double TReverStress;
    double TPlasticExcursion;
    
    /*** TRIAL State Variables ***/
    double Tstrain;
    double Tstress;
    double Ttangent; // Not really a state variable, but declared here
                     // for convenience

    // Calculates the trial state variables based on the trial strain
    void determineTrialState (double dStrain);
    double getR (double x_in);

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
	Matrix *SHVs;
// AddingSensitivity:END ///////////////////////////////////////////

};

#endif
