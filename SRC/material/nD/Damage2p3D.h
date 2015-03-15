/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.5 $
// $Date: 2008/08/26 16:22:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/Damage2p3D.h,v $


//
//

#ifndef Damage2p3DMaterial_h
#define Damage2p3DMaterial_h

#include <NDMaterial.h>
#include <Matrix.h>

class Damage2p3DMaterial : public NDMaterial
{
  public:
    Damage2p3DMaterial(int tag, 
		    double E,
		    double nu,
		    double r0,
		    double f0_t,
		    double f0_c,
		    double beta,
		    double Gf,
		    double A_c,
		    double B_c);
    Damage2p3DMaterial();
    ~Damage2p3DMaterial();

    const char *  getClassType(void) const {return "Damage2p3D";};

    int           setTrialStrain(const Vector &v);
	int           setTrialStrain(const Vector &v, const Vector &r);
    int           setTrialStrain(const Tensor &v);
	int           setTrialStrain(const Tensor &v, const Tensor &r);
    int           setTrialStrainIncr(const Vector &v);
    int           setTrialStrainIncr(const Vector &v, const Vector &r);

    const Vector& getStrain(void);          
    const Vector& getStress(void);
    const Matrix& getTangent(void);
    //const Tensor& getTangentTensor(void);
    double        signum(double);

    int           commitState(void);
    int           revertToLastCommit(void);
    int           revertToStart(void);

    NDMaterial    *getCopy(void);
    NDMaterial    *getCopy(const char *type);
    int           sendSelf(int commitTag, Channel &theChannel);  
    int           recvSelf(int commitTag, Channel &theChannel, 
		                   FEM_ObjectBroker &theBroker);    
    void          Print(OPS_Stream &s, int flag =0);

    const char *getType(void) const;
    //int getOrder(void) const {return 0;};

    //Response *setResponse (const char **argv, int argc, 
	//			   OPS_Stream &s);
    //int getResponse (int responseID, Information &matInformation);

	// Reliability and sensitivity stuff
    const Matrix& getInitialTangent        (void);
    int           setParameter (const char **argv, int argc, Parameter &param);
    int           updateParameter          (int parameterID, Information &info);
	int           activateParameter        (int parameterID);
	// AddingSensitivity:BEGIN //////////////////////////////////////////
    const Vector & getStressSensitivity     (int gradIndex, bool conditional);
    const Vector & getStrainSensitivity     (int gradIndex);
    //const Matrix & getTangentSensitivity    (int gradIndex);
    const Matrix & getInitialTangentSensitivity    (int gradIndex);
    //const Matrix & getDampTangentSensitivity(int gradIndex);
    //double         getRhoSensitivity        (int gradIndex);
    int            commitSensitivity        (Vector & strainGradient, int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

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
    
    double K, G, lambda; //, volmod elastic constants
    // History variables (trial and commited)
    double damage_t, damage_c; //state 1,2
    double r_t, r_c;           //state 3,4
    double Cdamage_t, Cdamage_c; //committed state 1,2
    double Cr_t, Cr_c;           //committed state 3,4
    double r0_t, r0_c;

	// Ohter variables
    Vector Tstrain, Cstrain;
    Vector Tstress, TeStress, Cstress, CeStress;
    Matrix tangent;

	double tolerance;
	int maxNumIter;

	// Sensitivit stuff
    int parameterID;
	Matrix *SHVs;
};


#endif

