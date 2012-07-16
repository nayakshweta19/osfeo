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
// UPDATE HISTORY:    Guanzhou Jie updated for parallel Dec 2006
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef Isotropic_Elastic_H
#define Isotropic_Elastic_H

#include "ElasticState.h"
#define ES_TAG_IsotropicElastic 123003

//stresstensor zerostress;
//straintensor zerostrain;

class Isotropic_Elastic : public ElasticState
{  
  public:
                    
    Isotropic_Elastic() : ElasticState(ES_TAG_IsotropicElastic) {}; //Guanzhou added for parallel processing

    Isotropic_Elastic(int E_in, 
                   int v_in,
                   const stresstensor& initialStress = zerostress, 
                   const straintensor& initialStrain = zerostrain);

// Nima Tafazzoli added for new material models (January 2010)    
    Isotropic_Elastic(const stresstensor& initialStress, 
                                     const straintensor& initialStrain);

    ElasticState *newObj();
    
    const BJtensor& getElasticStiffness(const MaterialParameter &MaterialParameter_in) const;
    
    //Guanzhou added for parallel
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

  private:
    
    double getE(const MaterialParameter &MaterialParameter_in) const;
    double getv(const MaterialParameter &MaterialParameter_in) const;
  
  private:

    static BJtensor ElasticStiffness;
    
    int E_index;   
    int v_index;
};


#endif

