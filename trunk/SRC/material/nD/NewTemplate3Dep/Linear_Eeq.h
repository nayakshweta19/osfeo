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

#ifndef Linear_Eeq_H
#define Linear_Eeq_H 

#include "ScalarEvolution.h"
#define SE_TAG_Linear_Eeq 121003

class Linear_Eeq : public ScalarEvolution
{
  public:
  
    Linear_Eeq() : ScalarEvolution(SE_TAG_Linear_Eeq) {}; //Guanzhou added for parallel processing

    Linear_Eeq(int LinearFactor_index_in);
    
    ScalarEvolution* newObj();

    double H(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
             const straintensor& Stra, const MaterialParameter& material_parameter);

    const tensor& DH_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                   const straintensor& Stra, const MaterialParameter& material_parameter);

    double DH_Diso(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                   const straintensor& Stra, const MaterialParameter& material_parameter);

    const tensor& DH_Dkin(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                   const straintensor& Stra, const MaterialParameter& material_parameter);
    
	//Guanzhou added for parallel
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
     
  private:
  
    double getLinearFactor(const MaterialParameter& material_parameter) const;  
    
  private:
  
    int LinearFactor_index;
    
};


#endif

