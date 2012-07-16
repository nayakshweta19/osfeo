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


//! 
//! @mainpage
//! \n
//! \n
//! This documentation is the API for von-Mises material model. For more information about theory background 
//! see the <A HREF="http://sokocalo.engr.ucdavis.edu/~jeremic/CG/LN.pdf#page=30" target="_blank">Lecture Notes
//! \n
//! \n
//!



#ifndef VM_YF_H
#define VM_YF_H

#include "YieldFunction.h"
#include <math.h>



class VM_YF : public YieldFunction
{
  public:

//!
//! Yield function for von-Mises:
//! inputs:
//! - k_which_in: to define the material parameter type (unconfined compression strength), 
//!              - 0 for the (constant) material parameter, 
//!              - 1 for the (initial) internal scalar, and 
//!              - 2 for the (initial) internal tensor (N/A)
//! - index_k_in: to locate the position in the defined MaterialParameter command
//! - alpha_which_in: to define the material parameter type (initial back stress), 
//!              - 0 for the (constant) material parameter, 
//!              - 1 for the (initial) internal scalar (N/A), and 
//!              - 2 for the (initial) internal tensor
//! - index_alpha_in: to locate the position in the defined MaterialParameter command
//!
  
    VM_YF(int k_which_in = -1, int index_k_in = 0, 
          int alpha_which_in = -1, int index_alpha_in = 0);

//! Deconstructor for von-Mises yield function
    ~VM_YF();
    
    const char *getYieldFunctionType(void) const {return "VM_YF";}; 
     
    YieldFunction *newObj();


//! YieldFunctionValue: to define the yield function, f = [(sij-aij)*(sij-aij)]^0.5 - sqrt(2.0/3)*k
//! inputs:  
//! - stresstensor &Stre: the tensor of stress
//! - MaterialParameter &MaterialParameter_in: the class of material parameter

    double YieldFunctionValue(const stresstensor &Stre, 
                              const MaterialParameter &MaterialParameter_in) const;

//! StressDerivative: to obtain the derivative of the yield function respect to the stress
//! inputs: 
//! - stresstensor &Stre: the tensor of stress
//! - MaterialParameter &MaterialParameter_in: the class of material parameter

    const stresstensor& StressDerivative(const stresstensor &Stre, 
                                         const MaterialParameter &MaterialParameter_in) const;
    
//! InScalarDerivative: to obtain the derivative of the yield function respect to the scalar (k, if appliable)
//! inputs: 
//! - stresstensor &Stre: the tensor of stress
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
//! - which: 0 if k is a constant material scalar, 1 if k is a changing internal scalar variable

    double InScalarDerivative(const stresstensor &Stre, 
                              const MaterialParameter &MaterialParameter_in, 
                              int which) const;

//! InTensorDerivative: to obtain the derivative of the yield function respect to the scalar (alpha_ij, if appliable)
//! inputs: 
//! - stresstensor &Stre: the tensor of stress
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
//! - which: 0 if alpha is a constant material tensor, 1 if alpha is a changing internal tensor variable
                           
    const stresstensor& InTensorDerivative(const stresstensor &Stre, 
                                           const MaterialParameter &MaterialParameter_in, 
                                           int which) const;

//! getNumInternalScalar: to get the number of scalar internal varibles
    int getNumInternalScalar() const;

//! getNumInternalTensor: to get the number of tensor internal varibles
    int getNumInternalTensor() const;

//! getYieldFunctionRank: to get the rank of the yield function to the unit of stress
    int getYieldFunctionRank() const;

	//Guanzhou added for parallel
//! sendSelf: for parallel computing
    int sendSelf(int commitTag, Channel &theChannel);  

//! recvSelf: for parallel computing
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

	
  private:
    
//! getk: to get k - the unconfined compression strength
//! inputs: 
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
    double getk(const MaterialParameter &MaterialParameter_in) const;
    
//! getbackstress: to get alpha - the back stress
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
    const stresstensor& getbackstress(const MaterialParameter &MaterialParameter_in) const;
    
  private:
      
    int k_which;
    int index_k;
    int alpha_which;
    int index_alpha;
        
    static stresstensor VMst;
    static stresstensor VMa;
};


#endif


