///////////////////////////////////////////////////////////////////////////////
//   COPYLEFT (C): Woody's viral GPL-like license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ours, cause we don't give
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

#ifndef PlasticFlow_CPP
#define PlasticFlow_CPP

#include "PlasticFlow.h"

tensor PlasticFlow::PF_tensorR2(2, def_dim_2, 0.0);
tensor PlasticFlow::PF_tensorR4(4, def_dim_4, 0.0);

///////////////////////////////////////////////////////////////////////////////
const tensor& PlasticFlow::Dm_Ds(const stresstensor &Stre, 
                                 const straintensor &Stra, 
                                 const MaterialParameter &MaterialParameter_in) const 
{
    tensor PF_R4(4, def_dim_4, 0.0);
    PlasticFlow::PF_tensorR4.Initialize(PF_R4);
    
    return PlasticFlow::PF_tensorR4;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& PlasticFlow::Dm_Diso(const stresstensor &Stre, 
                                 const straintensor &Stra, 
                                 const MaterialParameter &MaterialParameter_in) const 
{
    tensor PF_R2(2, def_dim_2, 0.0);
    PlasticFlow::PF_tensorR2.Initialize(PF_R2);
    
    return PlasticFlow::PF_tensorR2;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& PlasticFlow::Dm_Dkin(const stresstensor &Stre, 
                                 const straintensor &Stra, 
                                 const MaterialParameter &MaterialParameter_in) const 
{
    tensor PF_R4(4, def_dim_4, 0.0);
    PlasticFlow::PF_tensorR4.Initialize(PF_R4);
    
    return PlasticFlow::PF_tensorR4;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& PlasticFlow::Dm_Dkin2(const stresstensor &Stre, 
                                 const straintensor &Stra, 
                                 const MaterialParameter &MaterialParameter_in) const 
{
    tensor PF_R4(4, def_dim_4, 0.0);
    PlasticFlow::PF_tensorR4.Initialize(PF_R4);
    
    return PlasticFlow::PF_tensorR4;
}

#endif

