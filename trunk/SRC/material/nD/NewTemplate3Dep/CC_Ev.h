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
// UPDATE HISTORY:    Guanzhou Jie updated for parallel, Dec 2006
//                    Nima Tafazzoli updated for API (Feb 2009)
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef CC_Ev_H
#define CC_Ev_H 

#include "ScalarEvolution.h"
#include "ElasticState.h"
#define SE_TAG_CC 121001

class CC_Ev : public ScalarEvolution
{
  public:
  
//! Special volumetric isotropic hardening/softening law for Cam-Clay model
//! Evolution function for Modified Cam-Clay: 
//! inputs:
//! - M_index_in: to locate the position in the defined (constant) MaterialParameter command for M (critial state ratio)
//! - lambda_index_in: to locate the position in the defined (constant) MaterialParameter command for lambda (normal consolidation line slope)
//! - kappa_index_in: to locate the position in the defined (constant) MaterialParameter command for kappa (unloading-reloading line slope)
//! - e0_index_in: to locate the position in the defined (constant) MaterialParameter command for e0 (initial void ratio)
//! - p0_index_in: to locate the position in the defined (constant) MaterialParameter command for p0 (initial internal scalar variable)
//!

    CC_Ev(int M_index_in,
         int lambda_index_in,
         int kappa_index_in,
         int e0_index_in,
         int p0_index_in);

    CC_Ev() : ScalarEvolution(SE_TAG_CC) {};

    ScalarEvolution* newObj();

//! Evolution law for the internal variable p0
//! H: \dot{p0} = \dot{plastic multiplier} * H
//! inputs:
//! -PlasticFlow& plastic_flow:
//! - stresstensor &Stre: the tensor of stress
//! - straintensor &Stra: the tensor of strain
//! - MaterialParameter &MaterialParameter_in: the class of material parameter

    double H(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
             const straintensor& Stra, const MaterialParameter& material_parameter);

//! derivation of the evolution law for the internal variable p0 (H) with respect to stress
//! DH_Ds:
//! PlasticFlow& plastic_flow:
//! - stresstensor &Stre: the tensor of stress
//! - straintensor &Stra: the tensor of strain
//! - MaterialParameter &MaterialParameter_in: the class of material parameter

    const tensor& DH_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                   const straintensor& Stra, const MaterialParameter& material_parameter);

//! derivation of the evolution law for the internal variable p0 (H) with respect to p0
//! DH_Diso:
//! PlasticFlow& plastic_flow:
//! - stresstensor &Stre: the tensor of stress
//! - straintensor &Stra: the tensor of strain
//! - MaterialParameter &MaterialParameter_in: the class of material parameter

    double DH_Diso(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                   const straintensor& Stra, const MaterialParameter& material_parameter);
 
    //Guanzhou added for parallel
//! sendSelf: for parallel computing
    int sendSelf(int commitTag, Channel &theChannel);  
//! recvSelf: for parallel computing
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
  
      
  private:
    
//! getM: to get M (critical state ratio)
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 

    double getM(const MaterialParameter& material_parameter) const;

//! getlambda: to get lambda (normal consolidation line slope)
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 

    double getlambda(const MaterialParameter& material_parameter) const;

//! getkappa: to get kappa (unloading-reloading line slope)
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 

    double getkappa(const MaterialParameter& material_parameter) const;

//! gete0: to get e0 (initial viod ratio)
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 

    double gete0(const MaterialParameter& material_parameter) const;

//! getp0: to get p0 (initial internal scalar variable)
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 

    double getp0(const MaterialParameter& material_parameter) const; 
    
  private:
  
    int M_index;
    int lambda_index;
    int kappa_index;
    int e0_index;
    int p0_index;
};

#endif






