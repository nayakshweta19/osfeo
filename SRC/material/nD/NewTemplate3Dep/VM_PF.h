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
//                    Nima Tafazzoli updated for API (Feb 2009)
//
///////////////////////////////////////////////////////////////////////////////
//
#ifndef VM_PF_H
#define VM_PF_H

#include "PlasticFlow.h"
#include <math.h>

class VM_PF : public PlasticFlow
{
  public:   
  
//!
//! Plastic Flow function for von-Mises:
//! inputs:
//! - alpha_which: to define the material parameter type (initial back stress), 
//!              - 0 for the (constant) material parameter, 
//!              - 1 for the (initial) internal scalar (N/A), and 
//!              - 2 for the (initial) internal tensor
//! - index_alpha_in: to locate the position in the defined MaterialParameter command
//!

    VM_PF(int alpha_which = -1, int index_alpha_in = 0);
    
//! Deconstructor for von-Mises plastic flow function	
    ~VM_PF();
         
    PlasticFlow* newObj();
    
    const char *getPlasticFlowType(void) const {return "VM_PF";};

    
//! PlasticFlowTensor: to define the plastic flow
//! inputs: 
//! - stresstensor &Stre: the tensor of stress
//! - straintensor &Stra: the tensor of strain
//! - MaterialParameter &MaterialParameter_in: the class of material parameter

    const straintensor& PlasticFlowTensor(const stresstensor& Stre, 
                                          const straintensor& Stra, 
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

//! Dm_Dkin:
//! inputs:
//! - stresstensor &Stre: the tensor of stress
//! - straintensor &Stra: the tensor of strain
//! - MaterialParameter &MaterialParameter_in: the class of material parameter

    const tensor& Dm_Dkin(const stresstensor &Stre, 
                          const straintensor &Stra, 
                          const MaterialParameter &MaterialParameter_in) const;

    //Guanzhou added for parallel
//! sendSelf: for parallel computing
    int sendSelf(int commitTag, Channel &theChannel);  
//! recvSelf: for parallel computing
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

  private:

//! getalpha: to get alpha - the back stress
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    const stresstensor& getalpha(const MaterialParameter &MaterialParameter_in) const; 

  private:
    
    int alpha_which;
    int index_alpha;   
    
    static straintensor VMm;
    static stresstensor VMb;
};


#endif


