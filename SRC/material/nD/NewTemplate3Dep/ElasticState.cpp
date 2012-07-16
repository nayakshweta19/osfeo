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

#ifndef ElasticState_CPP
#define ElasticState_CPP

#include "ElasticState.h"

BJtensor ElasticState::ElasticCompliance(4, def_dim_4, 0.0);
const stresstensor ElasticState::zerostress(0.0);
const straintensor ElasticState::zerostrain(0.0);

////////////////////////////////////////////////////////////////
ElasticState::ElasticState(int clsTag, const stresstensor &initialStress, const straintensor &initialStrain)
:MovableObject(clsTag)
{
    Stress.Initialize(initialStress);
    Strain.Initialize(initialStrain);
}

////////////////////////////////////////////////////////////////
ElasticState::ElasticState(int clsTag, const stresstensor &initialStress)
:MovableObject(clsTag)
{
    Stress.Initialize(initialStress);

    straintensor ZeroStra;
    Strain.Initialize(ZeroStra);
}

////////////////////////////////////////////////////////////////
ElasticState::ElasticState(int clsTag)
:MovableObject(clsTag)
{
    stresstensor ZeroStre;
    Stress.Initialize(ZeroStre);
    
    straintensor ZeroStra;
    Strain.Initialize(ZeroStra);
}
                                     

////////////////////////////////////////////////////////////////
const stresstensor& ElasticState::getStress() const 
{ 
    return this->Stress;
}

////////////////////////////////////////////////////////////////
const straintensor& ElasticState::getStrain() const 
{ 
    return this->Strain;
}

/////////////////////////////////////////////////////////////////
int ElasticState::setStress(const stresstensor &Stre_in) 
{
    this->Stress.Initialize(Stre_in);
    
    return 0;
}

/////////////////////////////////////////////////////////////////
int ElasticState::setStrain(const straintensor &Stra_in) 
{
    this->Strain.Initialize(Stra_in);
    
    return 0;
}

#endif
