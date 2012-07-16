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

#ifndef elnp_Elastic_H
#define elnp_Elastic_H

#include "ElasticState.h"
#define ES_TAG_elnpElastic 123002
//stresstensor zerostress;
//straintensor zerostrain;

class elnp_Elastic : public ElasticState
{  
  public:
  
    elnp_Elastic() : ElasticState(ES_TAG_elnpElastic) {}; //Guanzhou added for parallel processing

	elnp_Elastic(int kappa_in, 
                 int v_in, 
                 int K_c_in,
                 int e0_in,
                 const stresstensor& initialStress = zerostress, 
                 const straintensor& initialStrain = zerostrain);

// Nima Tafazzoli added for new material models (Febuary 2010)    
    elnp_Elastic(const stresstensor& initialStress, 
                                     const straintensor& initialStrain);
                                       
    ElasticState* newObj();
    
    const BJtensor& getElasticStiffness(const MaterialParameter &MaterialParameter_in) const;
    
    //Guanzhou added for parallel
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

  private:
  
    double getkappa(const MaterialParameter &MaterialParameter_in) const;
    double getv(const MaterialParameter &MaterialParameter_in) const;
    double getK_c(const MaterialParameter &MaterialParameter_in) const;
    double gete0(const MaterialParameter &MaterialParameter_in) const;
     
  private:
    
    static BJtensor ElasticStiffness;
    
    int kappa_index;
    int v_index;
    int K_c_index;
    int e0_index;
};


#endif

