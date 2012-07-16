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
// UPDATE HISTORY:    Guanzhou Jie, changed this to the parallel interface, Dec 2006
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef AF_Eij_CPP
#define AF_Eij_CPP
#include <OPS_Globals.h>

#include <iostream>
using namespace std;

#include "AF_Eij.h"
	
stresstensor AF_Eij::AFal;

AF_Eij::AF_Eij(int ha_index_in, 
               int Cr_index_in,
               int alpha_index_in)
:TensorEvolution(TE_TAG_AF), ha_index(ha_index_in), 
  Cr_index(Cr_index_in),
  alpha_index(alpha_index_in)
{

}

TensorEvolution* AF_Eij::newObj()
{
    TensorEvolution* nObj = new AF_Eij(this->ha_index,
                                       this->Cr_index,
                                       this->alpha_index);
                                       
    return nObj;
}

///////////////////////////////////////////////////////////////////////////////
const straintensor& AF_Eij::Hij(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                                const straintensor& Stra, const MaterialParameter& material_parameter)
{
    straintensor PF = plastic_flow.PlasticFlowTensor(Stre, Stra, material_parameter); 

    double ha = getha(material_parameter);
    double Cr = getCr(material_parameter);
    stresstensor a = getalpha(material_parameter);

    TensorEvolution::TensorEvolutionHij = PF * (2.0*ha/3.0) - (a * PF.equivalent() *Cr); 
    
    return TensorEvolution::TensorEvolutionHij;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& AF_Eij::DHij_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
    double twoOthree = 2.0/3.0;
    double ha = getha(material_parameter);
    double Cr = getCr(material_parameter);
    stresstensor a = getalpha(material_parameter);

    straintensor PF = plastic_flow.PlasticFlowTensor(Stre, Stra, material_parameter);
    tensor dmOds = plastic_flow.Dm_Ds(Stre, Stra, material_parameter);
        
    TensorEvolution::TE_tensorR4 = dmOds * (twoOthree * ha);

    double m_eq =  sqrt(twoOthree * (PF("mn")*PF("mn")).trace());
    tensor tensor1 = PF("ij") * dmOds("ijmn");
    tensor1.null_indices();
    tensor tensor2 = tensor1("ij") * a("mn");
    tensor2.null_indices();
        
    if (m_eq != 0.0)
       TensorEvolution::TE_tensorR4 -= tensor2 * (twoOthree * Cr / m_eq);
    
    return TensorEvolution::TE_tensorR4;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& AF_Eij::DHij_Diso(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
    double twoOthree = 2.0/3.0;
    double ha = getha(material_parameter);
    double Cr = getCr(material_parameter);
    stresstensor a = getalpha(material_parameter);

    straintensor PF = plastic_flow.PlasticFlowTensor(Stre, Stra, material_parameter);
    tensor dmOdq = plastic_flow.Dm_Diso(Stre, Stra, material_parameter);
        
    TensorEvolution::TE_tensorR2 = dmOdq * (twoOthree * ha);

    double m_eq =  sqrt(twoOthree * (PF("mn")*PF("mn")).trace());
    double scalar1 = (PF("mn")*dmOdq("mn")).trace();

    if (m_eq != 0.0)     
      TensorEvolution::TE_tensorR2 -= a * (-scalar1 * twoOthree * Cr / m_eq);

    return TensorEvolution::TE_tensorR2;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& AF_Eij::DHij_Dkin(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
    double twoOthree = 2.0/3.0;
    double ha = getha(material_parameter);
    double Cr = getCr(material_parameter);
    stresstensor a = getalpha(material_parameter);

    straintensor PF = plastic_flow.PlasticFlowTensor(Stre, Stra, material_parameter);
    tensor dmOda = plastic_flow.Dm_Dkin(Stre, Stra, material_parameter);
        
    TensorEvolution::TE_tensorR4 = dmOda * (twoOthree * ha);

    double m_eq =  sqrt(twoOthree * (PF("mn")*PF("mn")).trace());
    tensor tensor1 = PF("ij") * dmOda("ijmn");
    tensor1.null_indices();
    tensor tensor2 = tensor1("ij") * a("mn");
    tensor2.null_indices();
        
    if (m_eq != 0.0)
       TensorEvolution::TE_tensorR4 -= tensor2 * (twoOthree * Cr / m_eq);

    tensor I2("I", 2, def_dim_2);
    tensor I_ijkl = I2("ij")*I2("kl");
    I_ijkl.null_indices();
    tensor I_ikjl = I_ijkl.transpose0110();
    tensor I_iljk = I_ijkl.transpose0111();
    tensor I4s = (I_ikjl+I_iljk)*0.5;
    
    TensorEvolution::TE_tensorR4 -= I4s * (m_eq * Cr);
    
    return TensorEvolution::TE_tensorR4;
}

///////////////////////////////////////////////////////////////////////////////
double AF_Eij::getha(const MaterialParameter& material_parameter) const
{
    if ( ha_index <= material_parameter.getNum_Material_Constant() && ha_index > 0)
        return material_parameter.getMaterial_Constant(ha_index -1);
    else {
        opserr << "AF_Eij: Invalid Input of " << ha_index << endln;
        exit (1);
    }
}

///////////////////////////////////////////////////////////////////////////////
double AF_Eij::getCr(const MaterialParameter& material_parameter) const
{
    if ( Cr_index <= material_parameter.getNum_Material_Constant() && Cr_index > 0)
        return material_parameter.getMaterial_Constant(Cr_index -1);
    else {
        opserr << "AF_Eij: Invalid Input of " << Cr_index << endln;
        exit (1);
    }
}

///////////////////////////////////////////////////////////////////////////////
const stresstensor& AF_Eij::getalpha(const MaterialParameter& material_parameter) const
{
    if ( alpha_index <= material_parameter.getNum_Internal_Tensor() && alpha_index > 0) {
        AF_Eij::AFal = material_parameter.getInternal_Tensor(alpha_index -1);
        return AF_Eij::AFal;
    }
    else {
        opserr << "AF_Eij: Invalid Input of " << alpha_index << endln;
        exit (1);
    }
}

//Guanzhou added for parallel
int AF_Eij::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();

    static ID idData(3);
    idData.Zero();

    idData(0) = ha_index;
    idData(1) = Cr_index;
    idData(2) = alpha_index;

    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "AF_Eij::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int AF_Eij::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(3);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "AF_Eij::recvSelf -- failed to recv ID\n";
	return -1;
    }

    ha_index    = idData(0);
    Cr_index    = idData(1);
    alpha_index = idData(2);
    
    return 0;
}

#endif

