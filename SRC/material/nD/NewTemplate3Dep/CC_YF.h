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

//!
//! @mainpage
//! \n
//! \n
//! This documentation is the API for Modified Cam-Clay material model. For more information about theory background 
//! see the <A HREF="http://sokocalo.engr.ucdavis.edu/~jeremic/CG/LN.pdf#page=30" target="_blank">Lecture Notes
//! \n
//! \n
//!


#ifndef CC_YF_H
#define CC_YF_H

#include "YieldFunction.h"
#include <math.h>

class CC_YF : public YieldFunction
{
  public:

//!
//! Yield function for Modified Cam-Clay: f=q^2-M^2*p*(p-p0)
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
//! - index_p0_in: to locate the position in the defined (scalar internal variable) MaterialParameter command for P0
//!

    CC_YF(int M_which_in = -1, int index_M_in = 0, 
          int p0_which_in = -1, int index_p0_in = 0);

//! Deconstructor for Modified Cam-Clay yield function
    ~CC_YF();
    
    const char *getYieldFunctionType(void) const {return "CC_YF";};
      
    YieldFunction *newObj();


//! YieldFunctionValue: to define the yield function, f=q^2-M^2*p*(p-p0)
//! inputs:
//! - stresstensor &Stre: the tensor of stress
//! - MaterialParameter &MaterialParameter_in: the class of material parameter

    double YieldFunctionValue(const stresstensor& Stre, 
                              const MaterialParameter &MaterialParameter_in) const;

//! StressDerivative: to obtain the derivative of the yield function respect to the stress
//! inputs:
//! - stresstensor &Stre: the tensor of stress
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
                      
    const stresstensor& StressDerivative(const stresstensor& Stre, 
                                         const MaterialParameter &MaterialParameter_in) const;

//! InScalarDerivative: to obtain the derivative of the yield function respect to the scalar (k, if appliable)
//! inputs: 
//! - stresstensor &Stre: the tensor of stress
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
//! - which: 0 if k is a constant material scalar, 1 if k is a changing internal scalar variable


    double InScalarDerivative(const stresstensor& Stre, 
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

//! getM: to get M
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter

    double getM(const MaterialParameter &MaterialParameter_in) const;

//! getP0:to get p0
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
    double getP0(const MaterialParameter &MaterialParameter_in) const;
        
  private:    
    int M_which;
    int index_M;
    int p0_which;
    int index_p0;
    
  private:
    static stresstensor CCst;
};

#endif


