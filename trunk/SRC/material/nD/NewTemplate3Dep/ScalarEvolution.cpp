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

#ifndef ScalarEvolution_CPP
#define ScalarEvolution_CPP

#include "ScalarEvolution.h"

tensor ScalarEvolution::SE_tensorR2(2, def_dim_2, 0.0);

///////////////////////////////////////////////////////////////////////////////
//Guanzhou ScalarEvolution::ScalarEvolution()
//Guanzhou {
//Guanzhou 
//Guanzhou }

///////////////////////////////////////////////////////////////////////////////
double ScalarEvolution::H(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
    return 0.0;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& ScalarEvolution::DH_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
    tensor SE_R2(2, def_dim_2, 0.0);
    ScalarEvolution::SE_tensorR2.Initialize(SE_R2);

    return ScalarEvolution::SE_tensorR2;
}

///////////////////////////////////////////////////////////////////////////////
double ScalarEvolution::DH_Diso(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
    return 0.0;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& ScalarEvolution::DH_Dkin(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
    tensor SE_R2(2, def_dim_2, 0.0);
    ScalarEvolution::SE_tensorR2.Initialize(SE_R2);

    return ScalarEvolution::SE_tensorR2;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& ScalarEvolution::DH_Dkin2(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
    tensor SE_R2(2, def_dim_2, 0.0);
    ScalarEvolution::SE_tensorR2.Initialize(SE_R2);

    return ScalarEvolution::SE_tensorR2;
}

#endif

