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
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef TensorEvolution_H
#define TensorEvolution_H 

#include <stresst.h>
#include <straint.h>
#include <BJtensor.h>

#include "MaterialParameter.h"
#include "ElasticState.h"
#include "PlasticFlow.h"

class TensorEvolution : public MovableObject
{
  public:
    
    TensorEvolution(int clsTag) : MovableObject(clsTag) {};//Guanzhou changed for parallel
    
    virtual ~TensorEvolution(){};
    
    virtual TensorEvolution *newObj() = 0;

    virtual const straintensor& Hij(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                                    const straintensor& Stra, const MaterialParameter& material_parameter);

    virtual const tensor& DHij_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                                    const straintensor& Stra, const MaterialParameter& material_parameter);

    virtual const tensor& DHij_Diso(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                                    const straintensor& Stra, const MaterialParameter& material_parameter);

    virtual const tensor& DHij_Dkin(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                                    const straintensor& Stra, const MaterialParameter& material_parameter);

    // for D-M model 
    virtual const tensor& DHij_Dkin2(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                                    const straintensor& Stra, const MaterialParameter& material_parameter);

    //Guanzhou added for parallel, pure virtual, subclasses must override
    virtual int sendSelf(int commitTag, Channel &theChannel) = 0;  
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) = 0;    

  protected:
    
    static straintensor TensorEvolutionHij;
    static tensor TE_tensorR2;
    static tensor TE_tensorR4;  
};


#endif

