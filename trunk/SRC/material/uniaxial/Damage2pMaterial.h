#ifndef DAMAGE2PMATERIAL_h
#define DAMAGE2PMATERIAL_h

#include <UniaxialMaterial.h>
#include <Matrix.h>

class Damage2p : public UniaxialMaterial
{
  public:
    Damage2p(int tag,
	 double Ec,
	 double nu,
	 double r0,
	 double f0_t,
	 double f0_c,
	 double beta,
	 double Gf,
	 double A_c,
	 double B_c);
	
    Damage2p(int tag, double fc);

    ~Damage2p();

    const char *getClassType(void) const {return "Damage2p";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        
    UniaxialMaterial *getCopy(void);
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    void Print(OPS_Stream &s, int flag =0);
    
    double getInitialTangent(void);

  protected:
    
  private:

    // Material parameters
    double E;
    double nu;
    double r0;
    double f0_t;
    double f0_c;
    double beta;
    double Gf;
    double A_c, A_t;
    double B_c, B_t;
    
    double fc;

    double G, K; // shear modulus and volume modulus
    // History variables (trial and committed)
    double lambda; //multiplier
    double Tstrain, Cstrain;
    double damage_t, damage_c; //state 1,2
    double r_t, r_c;           //state 3,4
    double Cdamage_t, Cdamage_c; //committed state 1,2
    double Cr_t, Cr_c;           //committed state 3,4
    double r0_t, r0_c;
    // Other variables
    double Tstress, TeStress, tangent;
    double Cstress, CeStress;

};


#endif



