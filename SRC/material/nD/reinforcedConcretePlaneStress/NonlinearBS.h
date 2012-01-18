//# REF: (C) Copyright 1999, The Regents of the University of California
#ifndef NonlinearBS_h
#define NonlinearBS_h

#define POS_INF_STRAIN 1.0e16
#define NEG_INF_STRAIN -1.0e16

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Vector.h>
#include <Matrix.h>
#include <NDMaterial.h>

class NonlinearBS : public NDMaterial
{
  public:
	//full constructor
	NonlinearBS(int tag, double fc, double epsc0, double fcu,
				double epscu, double rat, double ft, double Ets,
				double reinfD1, double reinfR1, double reinf_dmm1, double reinfD2, double reinfR2, double reinf_dmm2,
				double mom1p, double rot1p, double mom2p, double rot2p,	double mom3p, double rot3p,
				double mom1n, double rot1n, double mom2n, double rot2n, double mom3n, double rot3n,
				double pinchX, double pinchY,
				double damfc1 = 0.0, double damfc2 = 0.0,
				double beta = 0.0);
	
	NonlinearBS(int tag,double fc, double epsc0, double fcu,
				double epscu, double rat, double ft, double Ets,
				double reinfD1, double reinfR1, double reinf_dmm1, double reinfD2, double reinfR2, double reinf_dmm2,
				double mom1p, double rot1p, double mom2p, double rot2p,
				double mom1n, double rot1n, double mom2n, double rot2n,
				double pinchX, double pinchY,
				double damfc1 = 0.0, double damfc2 = 0.0,
				double beta = 0.0);
	
	//null constructor
	NonlinearBS(void);
	
	//destructor
	~NonlinearBS();
	
	const char *getClassType(void) const {return "NonlinearBS";};
	
	NDMaterial *getCopy (void);
	NDMaterial *getCopy (const char *type);
	const char *getType (void) const;
	
	int setTrialStrain(const Vector &strain_from_element);
	
	
	const Vector &getStrain(void);
	const Vector &getStress(void);
	const Matrix &getTangent(void);
	
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);
	
	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
	
	Response *setResponse (const char **argv, int argc, Information &matlnfo, OPS_Stream &s);
	int getResponse (int responseID, Information &matInformation);
	void Print(OPS_Stream &s, int flag =0);
	
  protected:
  
  private:
  
	void Tens_Envlp (double epsc, double epscP, double Alpha, double &sigc, double &Ect);
	void Compr_Envlp (double epsc, double epsc2, double &sigc, double &Ect);
	void Tens_UnloadingP (double epsc, double epsc2, double epscP, double sigcP,
			double ecmin, double dept, double friction, double Alpha, double &sigc, double &Ect);
	
	void Tens_UnloadingN (double epsc, double epsc2, double epscP, double sigcP,
			double ecmin, double dept, double friction, double Alpha, double &sigc, double &Ect);
	
	void Compr_Reloading (double epsc, double epsc2, double epscP, double sigcP,
			double ecmin, double dept, double friction, double Alpha, double &sigc, double &Ect);
	
	// matpar : Concrete FIXED PROPERTIES
	double fc; // concrete compression strength : mp(1)
	double epsc0; // strain at compression strength : mp(2)
	double fcu; // stress at ultimate (crushing) strain : mp(3)
	double epscu; // ultimate (crushing) strain : mp(4)
	double rat; // ratio between unloading slope at epscu and original slope : mp(5)
	double ft; // concrete tensile strength : mp(6)
	double Ets; // tension stiffening slope : mp(7)
	
	double reinfD1;
	double reinfR1;
	double reinf_dmm1;
	double reinfD2;
	double reinfR2;
	double reinf_dmm2;
	
	// hstvPX : Concerete HISTORY VARIABLES last committed step
	double ecminPX; // hstP(1)
	double deptPX; // hstP(2)
	double epsPX; // strain at previous converged step
	double sigPX; // stress at previous converged step
	double ePX; // stiffness modulus at last converged step;
	
	// hstvPY : Concerete HISTORY VARIABLES last committed step
	double ecminPY; // hstP(1)
	double deptPY; // hstP(2)
	double epsPY; // strain at previous converged step
	double sigPY; // stress at previous converged step
	double ePY; // stiffness modulus at last converged step;
	
	//shear strain
	double gammaXY;
	
	// hstvX: Concerete HISTORY VARIABLES current step
	double ecminX;
	double deptX;
	double sigX;
	double eX;
	double epsX;
	
	// friction
	double Tens_Alpha1;
	double Tens_Alpha2;
	
	double frictionX;
	double frictionY;
	
	double ReducedL;
	double ReducedT;
	
	double slip;
	double M;
	double Yint;
	
	//Shear Retention
	double ShearBeta;
	
	// hstvY: Concerete HISTORY VARIABLES current step
	double ecminY;
	double deptY;
	double sigY;
	double eY;
	double epsY;
	
	double tauXY;
	double thetaSP;
	double betaR;
	
	double eOXX;
	double eOYY;
	double gOXY;
	
	double sigOXX;
	double sigOYY;
	double tauOXY;
	
	double epsL;
	double epsT;
	double gammaLT;
	
	double sigL;
	double sigT;
	double tauLT;
	
	double dsigXdeX;
	double dsigXdeY;
	double dsigYdeX;
	double dsigYdeY;
	double dtauXYdgXY;
	
	double dsigLdeL;
	double dsigTdeT;
	
	double thetaC;
	double shear;
	
	static Vector sig_p;
	static Vector eps_p;
	static Vector strain_vec;
	static Vector eSlip;
	
	static Vector sig_XY;
	static Vector stress_vec;
	static Vector sigSlip;
	static Matrix tangent_matrix;
	static Matrix initialtangent_matrix;
	
	static Matrix iDs;
	static Matrix iDc;
	static Matrix Ds;
	static Matrix Dc;
	static Matrix T;
	static Matrix Ts1;
	static Matrix Ts2;
	static Matrix Dxy;
	
	// Pinching parameters
	double pinchX;   // Deformation pinching
	double pinchY;   // Force pinching
	
	// Damage parameters
	double damfc1;   // Deformation
	double damfc2;   // Energy
	
	// Unloading parameter
	double beta;
	
	// Trial history variables
	double TrotMaxL;
	double TrotMinL;
	double TrotPuL;
	double TrotNuL;
	double TenergyDL;
	int TloadIndicatorL;
	
	// Trial history variables
	double TrotMaxT;
	double TrotMinT;
	double TrotPuT;
	double TrotNuT;
	double TenergyDT;
	int TloadIndicatorT;
	
	// Trial state variables
	double TtangentL;
	double TstressL;
	double TstrainL;
	
	// Trial state variables
	double TtangentT;
	double TstressT;
	double TstrainT;
	
	// Converged history variables
	double CrotMaxL;
	double CrotMinL;
	double CrotPuL;
	double CrotNuL;
	double CenergyDL;
	int CloadIndicatorL;
	
	// Converged history variables
	double CrotMaxT;
	double CrotMinT;
	double CrotPuT;
	double CrotNuT;
	double CenergyDT;
	int CloadIndicatorT;
	
	// Converged state variables
	double CstressL;
	double CstrainL;
	
	// Converged state variables
	double CstressT;
	double CstrainT;
	
	// Backbone parameters
	double mom1p, rot1p;
	double mom2p, rot2p;
	double mom3p, rot3p;
	double mom1n, rot1n;
	double mom2n, rot2n;
	double mom3n, rot3n;
	
	double Rmom1p, Rrot1p;
	double Rmom2p, Rrot2p;
	double Rmom3p, Rrot3p;
	double Rmom1n, Rrot1n;
	double Rmom2n, Rrot2n;
	double Rmom3n, Rrot3n;
	
	double E1p, E1n;
	double E2p, E2n;
	double E3p, E3n;
	
	double rE1p, rE1n;
	double rE2p, rE2n;
	double rE3p, rE3n;
	
	double energyA;
	
	void setEnv(void);
	
	double setEnvelope(double Reduced);
	
	double posEnvlpStress(double strain, double Reduced);
	double negEnvlpStress(double strain, double Reduced);
	
	double posEnvlpTangent(double strain, double Reduced);
	double negEnvlpTangent(double strain, double Reduced);
	
	double posEnvlpRotlim(double strain, double Reduced);
	double negEnvlpRotlim(double strain, double Reduced);
	
	void positiveIncrementL(double dStrain, double Reduced);
	void negativeIncrementL(double dStrain, double Reduced);
	
	void positiveIncrementT(double dStrain, double Reduced);
	void negativeIncrementT(double dStrain, double Reduced);
};
#endif
	