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

#ifndef Linear_Eij_CPP
#define Linear_Eij_CPP

#include "Linear_Eij.h"

Linear_Eij::Linear_Eij(int LinearFactor_index_in)
: TensorEvolution(TE_TAG_Linear_Eij), LinearFactor_index(LinearFactor_index_in)
{

}

///////////////////////////////////////////////////////////////////////////////
TensorEvolution* Linear_Eij::newObj()
{
   TensorEvolution* nObj = new Linear_Eij(this->LinearFactor_index);

   return nObj;
}

///////////////////////////////////////////////////////////////////////////////
const straintensor& Linear_Eij::Hij(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                            const straintensor& Stra, const MaterialParameter& material_parameter)
{
    // mainly for D-P "rotational" model
    straintensor PF = plastic_flow.PlasticFlowTensor(Stre, Stra, material_parameter);    
    double L = getLinearFactor(material_parameter);
    
    TensorEvolution::TensorEvolutionHij = PF.deviator() * L; 
    
    return TensorEvolution::TensorEvolutionHij;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& Linear_Eij::DHij_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
    double L = getLinearFactor(material_parameter);
    tensor I2("I", 2, def_dim_2);
    tensor I_ijkl = I2("ij")*I2("kl");
    I_ijkl.null_indices();
    
    tensor dmOds = plastic_flow.Dm_Ds(Stre, Stra, material_parameter);
    tensor tensor1 = I_ijkl("pqij")*dmOds("pqmn");
    tensor1.null_indices();
    dmOds = dmOds - tensor1 *(1.0/3.0);
    
    TensorEvolution::TE_tensorR4 = dmOds * L;
    
    return TensorEvolution::TE_tensorR4;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& Linear_Eij::DHij_Diso(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
    double L = getLinearFactor(material_parameter);
    tensor I2("I", 2, def_dim_2);
    tensor I_ijkl = I2("ij")*I2("kl");
    I_ijkl.null_indices();
    
    tensor dmOdq = plastic_flow.Dm_Diso(Stre, Stra, material_parameter);
    tensor tensor1 = I_ijkl("pqij")*dmOdq("pq");
    tensor1.null_indices();
    dmOdq = dmOdq - tensor1 *(1.0/3.0);
    
    TensorEvolution::TE_tensorR2 = dmOdq * L;
    
    return TensorEvolution::TE_tensorR2;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& Linear_Eij::DHij_Dkin(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
    double L = getLinearFactor(material_parameter);
    tensor I2("I", 2, def_dim_2);
    tensor I_ijkl = I2("ij")*I2("kl");
    I_ijkl.null_indices();
    
    tensor dmOda = plastic_flow.Dm_Dkin(Stre, Stra, material_parameter);
    tensor tensor1 = I_ijkl("pqij")*dmOda("pqmn");
    tensor1.null_indices();
    dmOda = dmOda - tensor1 *(1.0/3.0);
    
    TensorEvolution::TE_tensorR4 = dmOda * L;
    
    return TensorEvolution::TE_tensorR4;
}

///////////////////////////////////////////////////////////////////////////////
double Linear_Eij::getLinearFactor(const MaterialParameter& material_parameter) const
{
    if ( LinearFactor_index <= material_parameter.getNum_Material_Constant() && LinearFactor_index > 0)
        return material_parameter.getMaterial_Constant(LinearFactor_index -1);
    else {
        opserr << "Linear_Eij: Invalid Input. " << endln;
        exit (1);
    }  
}

//Guanzhou added for parallel
int Linear_Eij::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(1);
    idData.Zero();

    idData(0) = LinearFactor_index; 
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "Linear_Eij::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int Linear_Eij::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(1);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "Linear_Eij::recvSelf -- failed to recv ID\n";
	return -1;
    }
    
    LinearFactor_index = idData(0);

    return 0;
}


#endif

