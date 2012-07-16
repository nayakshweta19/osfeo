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

#ifndef Isotropic_Elastic_CPP
#define Isotropic_Elastic_CPP

#include "Isotropic_Elastic.h"

BJtensor Isotropic_Elastic::ElasticStiffness(4, def_dim_4, 0.0);

Isotropic_Elastic::Isotropic_Elastic(int E_in,
                               int v_in,
                               const stresstensor& initialStress, 
                               const straintensor& initialStrain)
: ElasticState(ES_TAG_IsotropicElastic, initialStress, initialStrain),
  E_index(E_in),
  v_index(v_in)
{

}

// Nima Tafazzoli added for new Material Models
Isotropic_Elastic::Isotropic_Elastic(const stresstensor& initialStress, 
                                     const straintensor& initialStrain)
: ElasticState(ES_TAG_IsotropicElastic, initialStress, initialStrain)
{

}


// Create a new 
ElasticState* Isotropic_Elastic::newObj() 
{
    ElasticState *Els = new  Isotropic_Elastic(this->E_index, 
                                               this->v_index,
                                               this->Stress,
                                               this->Strain);
    return Els;
}

// Get Stiffness Tensor
const BJtensor& Isotropic_Elastic::getElasticStiffness(const MaterialParameter &MaterialParameter_in) const
{
    // Kronecker delta tensor
    BJtensor I2("I", 2, def_dim_2);

    BJtensor I_ijkl = I2("ij")*I2("kl");
    I_ijkl.null_indices();
    BJtensor I_ikjl = I_ijkl.transpose0110();
    BJtensor I_iljk = I_ijkl.transpose0111();
    BJtensor I4s = (I_ikjl+I_iljk)*0.5;
    
    double E = getE(MaterialParameter_in);
    double v = getv(MaterialParameter_in);
    
    if (E< 0.0 || v < -1.0 || v >= 0.5) {
      opserr << "Isotropic_Elastic: Invalid Input. " << endln;
      exit (1);
    }
        
    // Building elasticity tensor
    Isotropic_Elastic::ElasticStiffness = I_ijkl*( E*v / ( (1.0+v)*(1.0 - 2.0*v) ) ) + I4s*( E / (1.0 + v) );

    return Isotropic_Elastic::ElasticStiffness;
}

// Get Young's modulus
double Isotropic_Elastic::getE(const MaterialParameter &MaterialParameter_in) const
{
    if ( E_index > MaterialParameter_in.getNum_Material_Constant() || E_index < 1) {
      opserr << "Isotropic_Elastic: Invalid Input. " << endln;
      exit (1);
    }
    else    
      return MaterialParameter_in.getMaterial_Constant(E_index - 1);
}


// Get Poisson's ratio
double Isotropic_Elastic::getv(const MaterialParameter &MaterialParameter_in) const
{
    if ( v_index > MaterialParameter_in.getNum_Material_Constant() || v_index  < 1) { 
      opserr << "Isotropic_Elastic: Invalid Input. " << endln;
      exit (1);
    }
    else
      return MaterialParameter_in.getMaterial_Constant(v_index - 1);  
}

//Guanzhou added for parallel
int Isotropic_Elastic::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(2);
    idData.Zero();

    idData(0) = E_index;
    idData(1) = v_index;  
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "Isotropic_Elastic::sendSelf -- failed to send ID\n";
   	return -1;
    }

    if (theChannel.sendnDarray(dataTag, commitTag, this->Stress) < 0) {
    	opserr << "Isotropic_Elastic::sendSelf() -  failed to send nDarray Stress\n";
    	return -1;
    }
   
    if (theChannel.sendnDarray(dataTag, commitTag, this->Strain) < 0) {
    	opserr << "Isotropic_Elastic::sendSelf() -  failed to send nDarray Strain\n";
    	return -1;
    }
    
    return 0;
}

//Guanzhou added for parallel
int Isotropic_Elastic::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(2);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "Isotropic_Elastic::recvSelf -- failed to recv ID\n";
	return -1;
    }

    E_index     = idData(0);
    v_index     = idData(1);
     
    if (theChannel.recvnDarray(dataTag, commitTag, this->Stress) < 0) {
    	opserr << "Isotropic_Elastic::recvSelf() -  failed to recv nDarray Stress\n";
    	return -1;
    }

    if (theChannel.recvnDarray(dataTag, commitTag, this->Strain) < 0) {
    	opserr << "Isotropic_Elastic::recvSelf() -  failed to recv nDarray Strain\n";
    	return -1;
    }

    return 0;
}

#endif

