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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/BoucWen3DMaterial.h,v $


//
//

#ifndef BoucWen3DMaterial_h
#define BoucWen3DMaterial_h

#include <NDMaterial.h>
#include <Matrix.h>

class BoucWen3DMaterial : public NDMaterial
{
  public:
    BoucWen3DMaterial(int tag, 
		    double alpha,
		    double ko,
		    double n,
		    double gamma,
		    double beta,
		    double Ao,
		    double deltaA,
		    double deltaNu,
		    double deltaEta,
		    double tolerance,
		    int maxNumIter);
    BoucWen3DMaterial();
    ~BoucWen3DMaterial();

    const char *  getClassType(void) const {return "BoucWen3DMaterial";};

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
    double alpha;
    double ko;
	double n;
    double gamma;
    double beta;
    double Ao;
    double deltaA;
    double deltaNu;
    double deltaEta;

    // History variables (trial and commited)
    Vector Tstrain, Cstrain;
	Vector Tz, Cz;
	Vector Te, Ce;

	// Ohter variables
	Vector Tstress;
	Matrix Ttangent;

	double tolerance;
	int maxNumIter;

	// Sensitivit stuff
    int parameterID;
	Matrix *SHVs;
};


#endif

