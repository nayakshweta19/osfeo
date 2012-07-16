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
// UPDATE HISTORY:    Guanzhou Jie updated for parallel, Dec. 2006
//                    Nima Tafazzoli updated for API (Feb 2009)
//
///////////////////////////////////////////////////////////////////////////////
//
#ifndef CC_PF_H
#define CC_PF_H

#include "PlasticFlow.h"
#include <math.h>

class CC_PF : public PlasticFlow
{
  public:   

//!
//! Yield function for Modified Cam-Clay: Q=q^2-M^2*p*(p-p0)
//! inputs:
//! - M_which_in: to define the material parameter type (critial state ratio: M), 
//!               0 for the (constant) material parameter (set 0 for Cam-Clay)
//!               1 for the (initial) internal scalar (N/A)
//!               2 for the (initial) internal tensor (N/A)
//! - index_M_in: to locate the position in the defined (constant) MaterialParameter command for M
//!
//! - p0_which_in: to define the material parameter type (initial internal variable: p0), 
//!               0 for the (constant) material parameter (N/A)
//!               1 for the (initial) internal scalar (set 1 for Cam-Clay)
//!               2 for the (initial) internal tensor (N/A)
//! - index_p0_in: to locate the position in the defined (scalar internal variable) MaterialParameter command for p0
//!

    CC_PF(int M_which_in = -1, int index_M_in = 0, 
          int p0_which_in = -1, int index_p0_in = 0);
          
//! Deconstructor for Modified Cam-Clay plastic function
    ~CC_PF();
         
    PlasticFlow* newObj();

    const char *getPlasticFlowType(void) const {return "CC_PF";};    


//! PlasticFlowTensor: to define the plastic flow
//! inputs: 
//! - stresstensor &Stre: the tensor of stress
//! - straintensor &Stra: the tensor of strain
//! - MaterialParameter &MaterialParameter_in: the class of material parameter

    const straintensor& PlasticFlowTensor(const stresstensor &Stre, 
                                          const straintensor &Stra, 
                                          const MaterialParameter &MaterialParameter_in) const;


//! Dm_Ds: to obtain the derivative of the plastic flow respect to the scalar internal varible
//! inputs: 
//! - stresstensor &Stre: the tensor of stress
//! - straintensor &Stra: the tensor of strain
//! - MaterialParameter &MaterialParameter_in: the class of material parameter

    const tensor& Dm_Ds(const stresstensor &Stre, 
                        const straintensor &Stra, 
                        const MaterialParameter &MaterialParameter_in) const;


//! Dm_Diso: to obtain the derivative of the plastic flow respect to the tensor internal varible
//! inputs: 
//! - stresstensor &Stre: the tensor of stress
//! - straintensor &Stra: the tensor of strain
//! - MaterialParameter &MaterialParameter_in: the class of material parameter

    const tensor& Dm_Diso(const stresstensor &Stre, 
                          const straintensor &Stra, 
                          const MaterialParameter &MaterialParameter_in) const;

    //Guanzhou added for parallel
//! sendSelf: for parallel computing
    int sendSelf(int commitTag, Channel &theChannel);  
//! recvSelf: for parallel computing
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

  private:

//! getM: to get M
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getM(const MaterialParameter &MaterialParameter_in) const;

//! getp0: to get p0
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getP0(const MaterialParameter &MaterialParameter_in) const;

  private:
    
    int M_which;
    int index_M;   
    int p0_which;
    int index_p0; 
        
    static straintensor CCm;
};


#endif


