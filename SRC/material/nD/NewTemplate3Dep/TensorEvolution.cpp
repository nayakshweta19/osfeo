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

#ifndef TensorEvolution_CPP
#define TensorEvolution_CPP

#include "TensorEvolution.h"

straintensor TensorEvolution::TensorEvolutionHij;
tensor TensorEvolution::TE_tensorR2(2, def_dim_2, 0.0);
tensor TensorEvolution::TE_tensorR4(4, def_dim_4, 0.0);

///////////////////////////////////////////////////////////////////////////////
//Guanzhou TensorEvolution::TensorEvolution()
//Guanzhou {
//Guanzhou 
//Guanzhou }

///////////////////////////////////////////////////////////////////////////////
const straintensor& TensorEvolution::Hij(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                                         const straintensor& Stra, const MaterialParameter& material_parameter)
{   
    straintensor Z2;
    TensorEvolution::TensorEvolutionHij.Initialize(Z2);

    return TensorEvolution::TensorEvolutionHij;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& TensorEvolution::DHij_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                                         const straintensor& Stra, const MaterialParameter& material_parameter)
{   
    tensor TE_R4(4, def_dim_4, 0.0);
    TensorEvolution::TE_tensorR4.Initialize(TE_R4);
    
    return TensorEvolution::TE_tensorR4;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& TensorEvolution::DHij_Diso(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                                         const straintensor& Stra, const MaterialParameter& material_parameter)
{   
    tensor TE_R2(2, def_dim_2, 0.0);
    TensorEvolution::TE_tensorR2.Initialize(TE_R2);
    
    return TensorEvolution::TE_tensorR2;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& TensorEvolution::DHij_Dkin(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                                         const straintensor& Stra, const MaterialParameter& material_parameter)
{   
    tensor TE_R4(4, def_dim_4, 0.0);
    TensorEvolution::TE_tensorR4.Initialize(TE_R4);
    
    return TensorEvolution::TE_tensorR4;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& TensorEvolution::DHij_Dkin2(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                                         const straintensor& Stra, const MaterialParameter& material_parameter)
{   
    tensor TE_R4(4, def_dim_4, 0.0);
    TensorEvolution::TE_tensorR4.Initialize(TE_R4);
    
    return TensorEvolution::TE_tensorR4;
}

#endif

