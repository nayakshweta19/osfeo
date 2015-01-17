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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ConcreteDPM1.h,v $


//
//

#ifndef ConcreteDPM1_h
#define ConcreteDPM1_h

#include <NDMaterial.h>
#include <Matrix.h>
/**
* This class contains the combination of a local plasticity model for concrete with a local isotropic damage model.
* The yield surface of the plasticity model is based on the extension of the Menetrey and Willam yield criterion.
* The flow rule is nonassociated. The evolution laws of the hardening variables depend on the stress state.
* The plasticity model describes only hardening and perfect plasticity. It is based on the effective stress.
* The damage parameter of the isotropic damage model is based on the total volumetric strain.
* An exponential softening law is implemented.
*/
#define SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT

class ConcreteDPM1 : public NDMaterial
{
  public:
    ConcreteDPM1(int tag,
                  double fc,
                  double ft,
                  double epsc0,
                  double ecc,
                  double nu,
                  double AHard,
                  double BHard,
                  double CHard,
                  double DHard,
                  double yieldHardInitial,
                  double ASoft,
                  double helem,
                  double href);
    ConcreteDPM1();
    ~ConcreteDPM1();

    const char *  getClassType(void) const {return "ConcreteDPM1";};

    int           setTrialStrain(const Vector &v);
	int           setTrialStrain(const Vector &v, const Vector &r);
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

    virtual void giveRealStressVector(Vector &answer, const Vector &reducedStrain);//GaussPoint *gp, TimeStep *tStep, 

    /**
     * @param gp Gauss point.
     * @param strain Strain vector of this Gauss point.
     */
    void performPlasticityReturn(Vector &strain);//GaussPoint *gp, 

    /**
     * Check if the trial stress state falls within the vertex region of the plasticity model at the apex of triaxial extension or triaxial compression.
     * @returns true for vertex case and false if regular stress return can be used.
     * @param answer The volumetric apex stress.
     * @param sig The volumetric stress.
     * @param tempKappa The hardening variable.
     */

    bool checkForVertexCase(double &answer, const double sig, const double tempKappa);

    /**
     * Perform regular stress return for the plasticity model, i.e. if the trial stress state does not lie in the vertex region.
     * @param stress Stress vector which is computed.
     * @param gp Gauss point.
     */
    void performRegularReturn(Vector &stress); //GaussPoint *gp

    /**
     * Perform stress return for vertex case of the plasticity model, i.e. if the trial stress state lies within the vertex region.
     * @param stress Stress vector of this Gauss point.
     * @param apexStress Volumetric stress at the apex of the yield surface.
     * @param gp Gauss point.
     */
    void performVertexReturn(Vector &stress, double apexStress);//GaussPoint *gp

    /**
     * Compute the yield value based on stress and hardening variable.
     * @param sig Volumetric stress.
     * @param rho Length of the deviatoric stress.
     * @param theta Lode angle of the stress state.
     * @param tempKappa Hardening variable.
     * @return Yield value.
     */
    double computeYieldValue(const double sig, const double rho, const double theta, const double tempKappa) const;

    /**
     * Compute the value of the hardening function based on the hardening
     * variable.
     * @param tempKappa Hardening variable.
     * @return Value of the hardening function.
     */
    double computeHardeningOne(const double tempKappa) const;

    /**
     * Compute the derivative of the hardening function based on the
     * hardening parameter.
     * @param tempKappa Hardening variable.
     * @return The derivative of the hardening function.
     */
    double computeHardeningOnePrime(const double tempKappa) const;

    /**
     * Compute the derivative of the yield surface with respect to the hardening
     * variable based on the stress state and the hardening variable.
     * @param sig Volumetric stress.
     * @param rho Deviatoric length.
     * @param tempKappa Hardening variable.
     * @return The derivative of the yield surface.
     */
    double computeDFDKappa(const double sig, const double rho, const double tempKappa);

    /**
     * Compute the derivative of kappa with respect of delta
     * lambda based on the stress state and the hardening variable.
     * @param sig Volumetric stress.
     * @param rho Length of the deviatoric stress.
     * @param tempKappa Hardening variable.
     * @return Derivative of kappa with respect to delta lambda.
     */
    double computeDKappaDDeltaLambda(const double sig, const double rho, const double tempKappa);

    /**
     * Compute the ductility measure based on the stress state.
     * @param sig Volumetric stress.
     * @param rho Length of the deviatoric strength.
     * @param theta Lode angle of stress state.
     * @returns Ductility measure.
     */
    virtual double computeDuctilityMeasure(const double sig, const double rho, const double theta);

    /**
     * Computes the first derivative of the ductility measure with respect to
     * the invariants sig and rho based on the stress state and the hardening parameter.
     */
    void computeDDuctilityMeasureDInv(Vector &answer, const double sig, const double rho, const double tempKappa);

    /**
     * This matrix is the core of the closest point projection and collects
     * the derivative of the flow rule and the hardening parameters.
     */
    void computeAMatrix(Matrix &answer, const double sig, const double rho, const double tempKappa);

    /**
     * Here, the first derivative of the plastic potential with respect
     * to the invariants sig and rho are computed.
     */
    void computeDGDInv(Vector &answer, const double sig, const double rho, const double tempKappa);

    /**
     * This function computes the ratio of the volumetric and deviatoric component
     * of the flow direction. It is used within the vertex return to check,
     * if the vertex return is admissible.
     */
    double computeRatioPotential(const double sig, const double tempKappa);

    /**
     * Here, the second derivative of the plastic potential with respect to the
     * invariants sig and rho are computed.
     */
    void computeDDGDDInv(Matrix &answer, const double sig, const double rho, const double tempKappa);

    /**
     * Here, the mixed derivative of the plastic potential with respect
     * to the invariants and the hardening parameter are determined.
     */
    void computeDDGDInvDKappa(Vector &answer, const double sig, const double rho, const double tempKappa);

    /**
     * Computes the mixed derivative of the hardening parameter kappa with
     * respect to the plastic multiplier delta Lambda and the invariants sig
     * and rho.
     */
    void computeDDKappaDDeltaLambdaDInv(Vector &answer, const double sig, const double rho, const double tempKappa);

    /**
     * Computes the derivative of the evolution law of the hardening
     * parameter kappa with respect to the hardening variable kappa.
     */
    double computeDDKappaDDeltaLambdaDKappa(const double sig, const double rho, const double tempKappa);
    
    /**
     * Computes the derivative of the yield surface with respect to the
     * invariants sig and rho.
     */
    void computeDFDInv(Vector &answer, const double sig, const double rho, const double tempKappa) const;

    /**
     * Compute temporary kappa.
     */
    double computeTempKappa(const double kappaInitial, const double sigTrial, const double rhoTrial, const double sig);

    /**
     * Perform stress return for the damage model, i.e. if the trial stress state does not violate the plasticity surface.
     * @param strain Strain.
     * @param gp Gauss point.
     * @param tStep Time step.
     * @return Damage.
     */
    double computeDamage(const Vector &strain);// GaussPoint *gp, TimeStep *tStep

    /// Compute damage parameter.
    virtual double computeDamageParam(double kappa);//GaussPoint *gp

    /// Compute the damage-driving variable from given damage.
    double computeInverseDamage(double dam);//GaussPoint *gp

    /// Compute equivalent strain value.
    virtual void computeEquivalentStrain(double &kappaD, const Vector &elasticStrain);// GaussPoint *gp, TimeStep *tStep

    /// Compute the ductility measure for the damage model.
    double computeDuctilityMeasureDamage(const Vector &strain);//GaussPoint *gp

    /**
     * Initialize the characteristic length, if damage is not yet activated.
     */
    void initDamaged(double kappa, const Vector &elasticStrain);//GaussPoint *gp

    /// Compute the trial coordinates.
    void computeTrialCoordinates(const Vector &stress);//GaussPoint *gp

    /// Assign state flag.
    void assignStateFlag();//GaussPoint *gp

    /// Computes the derivative of rho with respect to the stress.
    void computeDRhoDStress(Vector &answer, const Vector &stress) const;

    /// Computes the derivative of sig with respect to the stress.
    void computeDSigDStress(Vector &answer) const;

    /// Computes the second derivative of rho with the respect to the stress.
    void computeDDRhoDDStress(Matrix &answer, const Vector &stress) const;

    /// Computes the derivative of costheta with respect to the stress.
    void computeDCosThetaDStress(Vector &answer, const Vector &stress) const;

    /// Compute the derivative of R with respect to costheta.
    double computeDRDCosTheta(const double theta, const double ecc) const;

    virtual void give3dMaterialStiffnessMatrix(Matrix &answer);//MatResponseMode mode, GaussPoint *gp, TimeStep *tStep


    virtual bool isCharacteristicMtrxSymmetric() { return false; }// MatResponseMode rMode

    //virtual int setIPValue(const Matrix &value); //GaussPoint *gp, InternalStateType type

    //virtual int giveIPValue(Matrix &answer); //GaussPoint *gp,InternalStateType type, TimeStep *tStep

    /**
    * Get the plastic strain deviator from the material status.
    * @return Plastic strain deviator.
    */
    const Vector &givePlasticStrain() const { return plasticStrain; }
    /**
    * Get the temp value of the full plastic strain vector from the material status.
    * @return Temp value of plastic strain vector.
    */
    const Vector &giveTempPlasticStrain() const { return tempPlasticStrain; }
    /**
    * Assign the temp value of the hardening variable of the damage model.
    * @param v New temp value of the hardening variable.
    */
    //void letTempKappaDBe(double v) { tempKappaD = v; }
    double ConcreteDPM1::computeMeanSize();

    void ConcreteDPM1::computePrincipalValues(Matrix &answer);

    double giveVolumetricPlasticStrain() const {
      return 1. / 3. * (plasticStrain(0) + plasticStrain(1) + plasticStrain(2));
    }

    /**
    * Get the deviatoric plastic strain norm from the material status.
    * @return Deviatoric plasticStrainNorm.
    */
    double giveDeviatoricPlasticStrainNorm() {
      Vector deviatoricPlasticStrain(plasticStrain.Size());
      double volumetricPlasticStrain;
      plasticStrain.computeDeviatoricVolumetricSplit(deviatoricPlasticStrain,
        volumetricPlasticStrain);
      return deviatoricPlasticStrain.computeStrainNorm();
    }

  protected:
    /**
    * Parameters of the yield surface of the plasticity model. 
    fc is the uniaxial compressive strength, 
    ft the uniaxial tensile strength and 
    ecc controls the out of roundness of the deviatoric section.
    */
    enum Concrete_VertexType { VT_Regular, VT_Tension, VT_Compression };
    Concrete_VertexType vertexType;

    double fc, ft, epsc0, E, ecc, href;

    /// Parameter of the ductilityMeasure of the plasticity model.
    double AHard;
    double BHard;
    double CHard;
    double DHard;

    double ASoft;/// Parameter of the ductilityMeasure of the damage model.
    double yieldHardInitial;/// Parameter of the hardening law of the plasticity model.
    double dilationConst;/// Control parameter for the volumetric plastic flow of the plastic potential

    double deltaLambda;/// Plastic multiplier of the plasticity model.

    double sig;/// the volumetric stress.
    double rho;/// The length of the deviatoric stress.

    double thetaTrial;/// The lode angle of the trial stress.

    double m;/// The friction parameter of the yield surface.

    double mQ;/// The dilation parameter of the plastic potential.

    double helem;/// Element size (to be used in fracture energy approach (crack band).

    //LinearElasticMaterial *linearElasticMaterial;/// Pointer for linear elastic material

    double eM;/// Elastic Young's modulus.
    double gM;/// Elastic shear modulus.
    double kM;/// Elastic bulk modulus.
    double nu;/// Elastic Poisson's ration.

    /// Hardening variable of plasticity model.
    double kappaP;
    double tempKappaP;

    /// Hardening variable of damage model.
    double kappaD;
    double tempKappaD;

    /// Damage variable of damage model.
    double damage;
    double tempDamage;

    /// Control parameter for the exponential softening law.
    double ef;

    /// Yield tolerance for the plasticity model.
    double yieldTol;

    /// Maximum number of iterations for stress return.
    int newtonIter;

    /// Stress and its deviatoric part.
    Vector effectiveStress;

  private:
    // Material parameters
    Vector plasticStrain;
    Vector tempPlasticStrain;
    double tempVolumetricPlasticStrain;

    double le;

    /// @name History variables of the damage model
    //@{
    double equivStrain;
    double tempEquivStrain;

    double deltaEquivStrain;
    //@}

    /// @name Indicates the state (i.e. elastic, unloading, plastic, damage, vertex) of the Gauss point
    //@{
    int state_flag;
    int temp_state_flag;
    //@}

#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    ///  @name History variable of the modified size-dependent adjustment
    /// (indicating value of omega*ft/E+kappaD at the onset of localization)
    double epsloc;
    double tempEpsloc;
#endif
    // History variables (trial and committed)
    Vector Tstrain, Cstrain;
	Vector Tz, Cz;
	Vector Te, Ce;

	// Other variables
	Vector Tstress;
	Matrix Ttangent;

	double tolerance;
	int maxNumIter;

	// Sensitivity stuff
    int parameterID;
	Matrix *SHVs;
};


#endif

