///////////////////////////////////////////////////////////////////////////////
//   COPYLEFT (C): Woody's viral GPL-like license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              
// CLASS:             
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           
//
// RETURN:
// VERSION:
// LANGUAGE:          C++
// TARGET OS:         
// DESIGNER:          Zhao Cheng, Boris Jeremic
// PROGRAMMER:        Zhao Cheng, 
// DATE:              Fall 2005
// UPDATE HISTORY:    06/2006, add functions for matrix based elements, CZ
//                    10/2006, add various more algorithms, CZ
//                    Guanzhou Jie updated for parallel Dec 2006
//                    Feb2008 adding ScaledExplicit, Mahdi's idea
//                    Nima Tafazzoli updated for API (Feb 2009)
/////////////////////////////////////////////////////////////////////////////
//


//! 
//! @mainpage NewTemplate3Dep Documentation
//! \n
//! This documentation covers the public APIs provided by the NewTemplate3dep libraries. It is intended mainly for programmers, those working on NewTemplate3Dep itself,
//! as well as developers extending functionality and/or material models intending to use these APIs. For more information about theory background 
//! see the <A HREF="http://sokocalo.engr.ucdavis.edu/~jeremic/research/index.html#LectureNotes" target="_blank">Lecture Notes
//! \n
//! \n
//!




#ifndef NewTemplate3Dep_H
#define NewTemplate3Dep_H

#include <stresst.h>
#include <straint.h>
#include <BJtensor.h>
#include <BJmatrix.h>
#include <BJvector.h>
#include <Matrix.h>
#include <Vector.h>

#include <NDMaterial.h>

#include "MaterialParameter.h"
#include "ElasticState.h"
#include "YieldFunction.h"
#include "PlasticFlow.h"
#include "ScalarEvolution.h"
#include "TensorEvolution.h"
//cout
#include <OPS_Globals.h>

// Boris Jeremic 17Nov2008
#include <iostream>
using namespace std;



#include <Channel.h>

using namespace std;



class NewTemplate3Dep : public NDMaterial
{

public:

//!
//! The NewTemplate3Dep class is a sub-class of NDMaterial and a container holding the components of material parameters, 
//! yield functions, plastic flows, scalar and tensor evolution laws
//! - tag                                            ---> Material tag number
//! - ptr_material_parameter_in  ---> pointer to MaterialParameter
//! - ptr_elastic_state_in              ---> pointer to ElasticState
//! - ptr_yield_function_in           ---> pointer to YieldFunction
//! - ptr_plastic_flow_in               ---> pointer to PlasticFlow
//! - ptr_scalar_evolution_in        ---> double pointer to ScalarEvolution
//! - ptr_tensor_evolution_in       ---> double pointer to TensorEvolution
//! - caseIndex_in                         ---> Integration algorithm which 0 (default) is explicit, 1 is implicit, 
//! - subStep_in                              ---> Substepping number (default = 1)
//! The multiple constructor functions enable the material model with both scalar and tensor evolution, 
//! with only scalar or tensor evolution, or without any evolution (perfectly-plasticity) applicable in this framework.
//! 
//! MaterialParameter: A class to contain the material parameters. The material parameters includes three type,
//! namely, material constants, initial scalar internal variables and initial tensor internal variables. The numbers of
//! these three types of material variables should be also stored in the class.
//!
//! ElasticState: Virtual base class for the elastic part of the material models (the most important function
//! is getElasticStiffness, which returns the elastic stiffness tensor)
//!
//! YieldFunction: Virtual base class for the properties of the yield function.
//! - YieldFunctionValue function is to evaluate the value of the yield function
//! - StressDerivative function gives the derivative of the yield function to the stresses
//! - InScalarDerivative function returns the derivative of the yield function to the scalar internal variable
//! - IntensorDerivative function returns the derivative of the yield function to the tensor internal variable
//!
//! PlasticFlow class is the virtual base class mainly to return the plastic flow tensor.
//!
//! ScalarEvolution class is the virtual base class to evaluate the evolution modulus of the scalar evolution H,
//! which may be the function of the stress, strain, and material parameters.
//!
//! TensorEvolution class is the virtual base class to evaluate the evolution modulus of the tensor evolution Hij,
//! which may be the function of the plastic flow, as well as the stresses, strain, and material parameters.
//!

 


    NewTemplate3Dep( int tag,
                  MaterialParameter *ptr_material_parameter_in,
                  ElasticState      *ptr_elastic_state_in,
                  YieldFunction     *ptr_yield_function_in,        
                  PlasticFlow       *ptr_plastic_flow_in,
                  ScalarEvolution   **ptr_scalar_evolution_in,
                  TensorEvolution   **ptr_tensor_evolution_in,
                  int caseIndex_in, int subStep_in = 1);


//!
//! The NewTemplate3Dep class is a sub-class of NDMaterial and a container holding the components of material parameters, 
//! yield functions, plastic flows, scalar and tensor evolution laws
//! - tag                                            ---> Material tag number
//! - ptr_material_parameter_in  ---> pointer to MaterialParameter
//! - ptr_elastic_state_in              ---> pointer to ElasticState
//! - ptr_yield_function_in           ---> pointer to YieldFunction
//! - ptr_plastic_flow_in               ---> pointer to PlasticFlow
//! - ptr_scalar_evolution_in        ---> double pointer to ScalarEvolution
//! - caseIndex_in                         ---> Integration algorithm which 0 (default) is explicit, 1 is implicit, 
//! - subStep_in                              ---> Substepping number (default = 1)
//! The multiple constructor functions enable the material model with both scalar and tensor evolution, 
//! with only scalar or tensor evolution, or without any evolution (perfectly-plasticity) applicable in this framework.
//! 
//! MaterialParameter: A class to contain the material parameters. The material parameters includes three type,
//! namely, material constants, initial scalar internal variables and initial tensor internal variables. The numbers of
//! these three types of material variables should be also stored in the class.
//!
//! ElasticState: Virtual base class for the elastic part of the material models (the most important function
//! is getElasticStiffness, which returns the elastic stiffness tensor)
//!
//! YieldFunction: Virtual base class for the properties of the yield function.
//! - YieldFunctionValue function is to evaluate the value of the yield function
//! - StressDerivative function gives the derivative of the yield function to the stresses
//! - InScalarDerivative function returns the derivative of the yield function to the scalar internal variable
//! - IntensorDerivative function returns the derivative of the yield function to the tensor internal variable
//!
//! PlasticFlow class is the virtual base class mainly to return the plastic flow tensor.
//!
//! ScalarEvolution class is the virtual base class to evaluate the evolution modulus of the scalar evolution H,
//! which may be the function of the stress, strain, and material parameters.
//!
//!

     NewTemplate3Dep( int tag,
                  MaterialParameter *ptr_material_parameter_in,
                  ElasticState      *ptr_elastic_state_in,
                  YieldFunction     *ptr_yield_function_in,        
                  PlasticFlow       *ptr_plastic_flow_in,
                  ScalarEvolution   **ptr_scalar_evolution_in,
                  int caseIndex_in, int subStep_in = 1);


//!
//! The NewTemplate3Dep class is a sub-class of NDMaterial and a container holding the components of material parameters, 
//! yield functions, plastic flows, scalar and tensor evolution laws
//! - tag                                            ---> Material tag number
//! - ptr_material_parameter_in  ---> pointer to MaterialParameter
//! - ptr_elastic_state_in              ---> pointer to ElasticState
//! - ptr_yield_function_in           ---> pointer to YieldFunction
//! - ptr_plastic_flow_in               ---> pointer to PlasticFlow
//! - ptr_tensor_evolution_in       ---> double pointer to TensorEvolution
//! - caseIndex_in                         ---> Integration algorithm which 0 (default) is explicit, 1 is implicit, 
//! - subStep_in                              ---> Substepping number (default = 1)
//! The multiple constructor functions enable the material model with both scalar and tensor evolution, 
//! with only scalar or tensor evolution, or without any evolution (perfectly-plasticity) applicable in this framework.
//! 
//! MaterialParameter: A class to contain the material parameters. The material parameters includes three type,
//! namely, material constants, initial scalar internal variables and initial tensor internal variables. The numbers of
//! these three types of material variables should be also stored in the class.
//!
//! ElasticState: Virtual base class for the elastic part of the material models (the most important function
//! is getElasticStiffness, which returns the elastic stiffness tensor)
//!
//! YieldFunction: Virtual base class for the properties of the yield function.
//! - YieldFunctionValue function is to evaluate the value of the yield function
//! - StressDerivative function gives the derivative of the yield function to the stresses
//! - InScalarDerivative function returns the derivative of the yield function to the scalar internal variable
//! - IntensorDerivative function returns the derivative of the yield function to the tensor internal variable
//!
//! PlasticFlow class is the virtual base class mainly to return the plastic flow tensor.
//!
//! TensorEvolution class is the virtual base class to evaluate the evolution modulus of the tensor evolution Hij,
//! which may be the function of the plastic flow, as well as the stresses, strain, and material parameters.
//!

    
     NewTemplate3Dep( int tag,
                  MaterialParameter *ptr_material_parameter_in,
                  ElasticState      *ptr_elastic_state_in,
                  YieldFunction     *ptr_yield_function_in,        
                  PlasticFlow       *ptr_plastic_flow_in,
                  TensorEvolution   **ptr_tensor_evolution_in,
                  int caseIndex_in, int subStep_in = 1);


//!
//! The NewTemplate3Dep class is a sub-class of NDMaterial and a container holding the components of material parameters, 
//! yield functions, plastic flows, scalar and tensor evolution laws
//! - tag                                            ---> Material tag number
//! - ptr_material_parameter_in  ---> pointer to MaterialParameter
//! - ptr_elastic_state_in              ---> pointer to ElasticState
//! - ptr_yield_function_in           ---> pointer to YieldFunction
//! - ptr_plastic_flow_in               ---> pointer to PlasticFlow
//! - caseIndex_in                         ---> Integration algorithm which 0 (default) is explicit, 1 is implicit, 
//! - subStep_in                              ---> Substepping number (default = 1)
//! The multiple constructor functions enable the material model with both scalar and tensor evolution, 
//! with only scalar or tensor evolution, or without any evolution (perfectly-plasticity) applicable in this framework.
//! 
//! MaterialParameter: A class to contain the material parameters. The material parameters includes three type,
//! namely, material constants, initial scalar internal variables and initial tensor internal variables. The numbers of
//! these three types of material variables should be also stored in the class.
//!
//! ElasticState: Virtual base class for the elastic part of the material models (the most important function
//! is getElasticStiffness, which returns the elastic stiffness tensor)
//!
//! YieldFunction: Virtual base class for the properties of the yield function.
//! - YieldFunctionValue function is to evaluate the value of the yield function
//! - StressDerivative function gives the derivative of the yield function to the stresses
//! - InScalarDerivative function returns the derivative of the yield function to the scalar internal variable
//! - IntensorDerivative function returns the derivative of the yield function to the tensor internal variable
//!
//! PlasticFlow class is the virtual base class mainly to return the plastic flow tensor.
//!



     NewTemplate3Dep( int tag,
                  MaterialParameter *ptr_material_parameter_in,
                  ElasticState      *ptr_elastic_state_in,
                  YieldFunction     *ptr_yield_function_in,        
                  PlasticFlow       *ptr_plastic_flow_in,
                  int caseIndex_in, int subStep_in = 1);
    
    NewTemplate3Dep(void);

    
    ~NewTemplate3Dep(void);

    const char *getClassType(void) const {return "NewTemplate3Dep";};

	

//! setTrialStrain and setTrialStrainIncr are the functions to drive material point state
//! from the current step to the next using specific integration algorithm.

	
    // methods to set and retrieve state using the Matrix class, added 06/2006
    int setTrialStrain (const Vector &v);
    int setTrialStrain (const Vector &v, const Vector &r);
    int setTrialStrainIncr (const Vector &v);
    int setTrialStrainIncr (const Vector &v, const Vector &r);


//! getTangent:

    const Matrix &getTangent (void);


//! getTangentBJmatrix:

    const BJmatrix& getTangentBJmatrix(void);


//! getStress: to obtain the current stress in form of vector

    const Vector &getStress (void);


//! getStrain: to obtain the current strain in form of vector

    const Vector &getStrain (void);

    // methods to set and retrieve state using the Tensor class    
    int setTrialStrain(const Tensor& v);
    int setTrialStrain(const Tensor& v, const Tensor& r);    
    int setTrialStrainIncr(const Tensor& v);
    int setTrialStrainIncr(const Tensor& v, const Tensor& r);


//! getTangentTensor: 

    const BJtensor& getTangentTensor(void);


//! getStressTensor: to obtain the current stress in form of tensor

    const stresstensor& getStressTensor(void);


//! getStrainTensor: to obtain the current strain in form of tensor

    const straintensor& getStrainTensor(void);


//! getPlasticStrainTensor: to obtain the current plastic strain in form of tensor

    const straintensor& getPlasticStrainTensor(void);

// Nima
    const stresstensor& getPredictorStressTensor(void);

// Nima
    const stresstensor& getintersection_stress(void);

// Nima
    const straintensor& getintersection_strain(void);

// Nima
    double getEnergy(double);
    double resetEnergy(void);


// moved out BJ 8Jan2008
//    straintensor BardetConstraint(int , double);


    double getRho();

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    
    NDMaterial* getCopy(void);
    NDMaterial* getCopy(const char *code);

    const char *getType(void) const;

    //Guanzhou
	int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream& s, int flag =0);
  
private:


//! Explicit: function to use the explicit integration method

    int Explicit(const straintensor& strain_incr, int NumStep_in = 1);


//! Implicit: function to use the implicit integration method

    int Implicit(const straintensor& strain_incr, int NumStep_in = 1);      


//! ImplicitLineSearch: function to use the implicit linesearch integration method

    int ImplicitLineSearch(const straintensor& strain_incr);


//! ScaledExplicit: function to use the scaled explicit integration method

    int ScaledExplicit(const straintensor& strain_incr, int NumStep_in = 1);


//! PredictorEPState:

    int PredictorEPState(const straintensor& strain_incr);


//! yield_surface_cross: trying to find intersection point

    stresstensor yield_surface_cross(const stresstensor& start_stress, 
                                     const stresstensor& end_stress,
				     double a);


//! zbrentstress: Routine used by yield_surface_cross to find the stresstensor at cross point

    double zbrentstress(const stresstensor& start_stress,
                        const stresstensor& end_stress,
                        double x1, double x2, double tol) const;


//! func: 

    double func( const stresstensor& start_stress,
                 const stresstensor& end_stress,
                 const MaterialParameter& ptr_material_parameter, 
                 double alfa ) const;



//! Tensor2MatrixSysR4:

    int Tensor2MatrixSysR4(const tensor& T, Matrix& M);    


//! Tensor2VectorSysR2:

    int Tensor2VectorSysR2(const tensor& T, Vector& V);


//! Matrix2TensorSysR4:

    int Matrix2TensorSysR4(const Matrix& M, tensor& T);


//! Vector2TensorSysR2:

    int Vector2TensorSysR2(const Vector& V, tensor& T, int Num0 = 0);

public: // temp, for test only, to be deleted


//! Stiffness2Compliance:

    int Stiffness2Compliance(const tensor& S, tensor& C);
 
private:



//! TrialStrain: a strain tensor contains the current strain

    straintensor TrialStrain;


//! TrialStress: a stress tensor contains the current stress

    stresstensor TrialStress;


//! TrialPlastic_Strain: a strain tensor contains the current plastic strain

    straintensor TrialPlastic_Strain;

// Nima
    stresstensor TrialPredictor_Stress;
    stresstensor Trialintersection_stress;
    straintensor Trialintersection_strain;

// Nima : For Energy
    stresstensor tStress;
    straintensor tStrain;
    stresstensor pStress;
    straintensor pStrain;
    straintensor incStrain;
    stresstensor incStress;
    stresstensor TotalMultiply;
    double TotalEnergy;
//    MatPoint3D ** matpoint;
    NDMaterial * pMat;


//! CommitStress:
    
    stresstensor CommitStress;


//! CommitStrain:

    straintensor CommitStrain; 


//! CommitPlastic_Strain:

    straintensor CommitPlastic_Strain;


//! Stiffness:
  
    BJtensor Stiffness;
                
    MaterialParameter *ptr_material_parameter;
    ElasticState      *ptr_elastic_state;           
    YieldFunction     *ptr_yield_function;
    PlasticFlow       *ptr_plastic_flow;
    ScalarEvolution   **ptr_scalar_evolution;
    TensorEvolution   **ptr_tensor_evolution;
    
    int caseIndex;
    int subStep;
    int divergeOrnot;
    
    static const  straintensor ZeroStrain;
    static const  stresstensor ZeroStress;
    static const  BJtensor ZeroI4;
    static const int ISMAX;
    static const int ITMAX;
    static const double TOL;
    static const double FTOL;

    // For Matrix based elements
    static Vector sigma; 
    static Matrix D;
//BJmatrix NewTemplate3Dep::TangentMatrix( 6, 6, 0.0);
    static BJmatrix TangentMatrix;
    static Vector epsilon;
};

#endif
