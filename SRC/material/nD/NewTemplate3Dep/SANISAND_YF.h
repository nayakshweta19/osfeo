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

//!
//! @mainpage
//! \n
//! \n
//! This documentation is the API for SANISAND material model. For more information about theory background 
//! see the <A HREF="http://sokocalo.engr.ucdavis.edu/~jeremic/CG/LN.pdf#page=30" target="_blank">Lecture Notes
//! \n
//! \n
//!


#ifndef SANISAND_YF_H
#define SANISAND_YF_H

#include "YieldFunction.h"
#include <math.h>
#define YF_TAG_SANISAND 124005

class SANISAND_YF : public YieldFunction
{
  public:

//! Reference: Taiebat & Dafalias (2008)
//! Yield function for SANISAND:
//! inputs:
//! - m_which_in:      type of parameter m (yield surface opening) --> 0
//! - index_m_in:      id of parameter   m 
//! - p0_which_in:     type of parameter p0 (yield surface size) --> 1
//! - index_p0_in:     id of parameter   p0
//! - alpha_which_in:  type of parameter alpha (back stress ratio) --> 2 
//! - index_alpha_in:  id of parameter   alpha
//!
    SANISAND_YF(int m_which_in,      int index_m_in, 
                int p0_which_in,     int index_p0_in,
                int alpha_which_in,  int index_alpha_in);

//! Deconstructor for SANISAND yield function
    ~SANISAND_YF(); 
    
    SANISAND_YF() : YieldFunction(YF_TAG_SANISAND) {};

    YieldFunction *newObj();
        
    const char *getYieldFunctionType(void) const {return "SANISAND_YF";};   


//! YieldFunctionValue: to define the yield function
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
//! - which: 0 if M is a constant material scalar, 1 if M is a changing internal scalar variable
      
    double InScalarDerivative(const stresstensor &Stre, 
                              const MaterialParameter &MaterialParameter_in, 
                              int which) const; 

//! InTensorDerivative: to obtain the derivative of the yield function respect to the scalar (alpha_ij, if appliable)
//! inputs: 
//! - stresstensor &Stre: the tensor of stress
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
//! - which: 0 if alpha is a constant material tensor, 1 if alpha is a changing internal tensor variable

    const stresstensor& InTensorDerivative(const stresstensor& Stre, 
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
//! getm: to get m (constant in yield function)
//! inputs: 
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
    double getm(const MaterialParameter &MaterialParameter_in) const;

//! getp0: to get p0 
//! inputs: 
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
    double getp0(const MaterialParameter &MaterialParameter_in) const;

//! getalpha: to get alpha
//! inputs: 
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
    const stresstensor& getalpha(const MaterialParameter &MaterialParameter_in) const;
    
  private:    
    int m_which; 
    int index_m;
    int p0_which; 
    int index_p0;
    int alpha_which; 
    int index_alpha;
        
    static stresstensor SANISANDst;
};


#endif

