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
// DESIGNER:          Mahdi Taiebat, Boris Jeremic
// PROGRAMMER:        Mahdi Taiebat, Boris Jeremic 
// Note:              
// DATE:              Spring 2007
// UPDATE HISTORY:    Nima Tafazzoli updated for API (Mar 2009)
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef SANISAND_Elastic_H
#define SANISAND_Elastic_H

#include "ElasticState.h"
#define ES_TAG_SANISAND 123005

class SANISAND_Elastic : public ElasticState
{  
  public:
  
    SANISAND_Elastic() : ElasticState(ES_TAG_SANISAND) {}; //Guanzhou added for parallel processing


//! Reference: Taiebat & Dafalias (2008)
//! Elastic function for SANISAND:
//! inputs:
//! - G0_in:                         Reference elastic shear modulus (same unit as stress);
//! - K0_in:                         Reference elastic bulk modulus (same unit as stress); 
//! - Pat_in:                        Atmospheric pressure;
//! - k_c_in:                        cut-off factor;                
//!                                  for p < k_c*Pat, use p = k_c*Pat for calculation of G;
//!                                  (a default value of k_c = 0.01 should work fine)
//! - e0_in:                         initial void ratio;
//! - stresstensor& initialStress:   stress tensor;
//! - straintensor& initialStrain:   strain tensor;

    SANISAND_Elastic(int G0_in, 
                     int K0_in,
                     int Pat_in,
                     int k_c_in,
                     int e0_in,
                     const stresstensor& initialStress = zerostress, 
                     const straintensor& initialStrain = zerostrain);                    
    
    ElasticState* newObj();

//! Reference: Taiebat & Dafalias (2008)
//! getElasticStiffness: 
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
    const BJtensor& getElasticStiffness (const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! getStress: 
    const stresstensor& getStress() const;    
    
     //Guanzhou added for parallel
//! sendSelf: for parallel computing
    int sendSelf(int commitTag, Channel &theChannel);  
//! recvSelf: for parallel computing
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

  private:
  
//! Reference: Taiebat & Dafalias (2008)
//! getG0: to get 
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getG0(const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! getK0: to get 
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getK0(const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! getPat: to get 
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getPat(const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! getk_c: to get 
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getk_c(const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! gete0: to get 
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double gete0(const MaterialParameter &MaterialParameter_in) const;    
     
  private:

    static BJtensor ElasticStiffness;
  
    int G0_index;   
    int K0_index;
    int Pat_index;
    int k_c_index;
    int e0_index;    
};

#endif

