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
// UPDATE HISTORY:    
//
///////////////////////////////////////////////////////////////////////////////
//
#ifndef YieldFunction_H
#define YieldFunction_H

#include <stresst.h>
#include "MaterialParameter.h"
#include <FEM_ObjectBroker.h>
#include <Channel.h>

class YieldFunction  : public MovableObject
{
  public:    

    YieldFunction(int clsTag) : MovableObject(clsTag) {};//Guanzhou changed for parallel

    virtual ~YieldFunction() {};
    
    virtual YieldFunction *newObj() = 0;

    virtual const char *getYieldFunctionType(void) const = 0;  

    virtual double YieldFunctionValue(const stresstensor &Stre, 
                                      const MaterialParameter &MaterialParameter_in) const = 0;
                     
    virtual const stresstensor& StressDerivative(const stresstensor &Stre, 
                                                 const MaterialParameter &MaterialParameter_in) const = 0;
    
    virtual double InScalarDerivative(const stresstensor& Stre, 
                                      const MaterialParameter &MaterialParameter_in, 
                                      int which) const;
                                   
    virtual const stresstensor& InTensorDerivative(const stresstensor& Stre, 
                                                   const MaterialParameter &MaterialParameter_in, 
                                                   int which) const;

    //Guanzhou added for parallel, pure virtual, subclasses must override
    virtual int sendSelf(int commitTag, Channel &theChannel) = 0;  
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) = 0;    

	virtual int getNumInternalScalar() const = 0;
    virtual int getNumInternalTensor() const = 0;
    virtual int getYieldFunctionRank() const = 0;
    
  private:
    
    static stresstensor stressYF;    
};


#endif

