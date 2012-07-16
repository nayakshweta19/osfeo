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

#ifndef Linear_Eeq_CPP
#define Linear_Eeq_CPP

#include "Linear_Eeq.h"

///////////////////////////////////////////////////////////////////////////////
Linear_Eeq::Linear_Eeq(int LinearFactor_index_in)
: ScalarEvolution(SE_TAG_Linear_Eeq), LinearFactor_index(LinearFactor_index_in)
{

}

///////////////////////////////////////////////////////////////////////////////
ScalarEvolution* Linear_Eeq::newObj()
{
    ScalarEvolution* nObj = new Linear_Eeq(this->LinearFactor_index);
    
    return nObj;
}

///////////////////////////////////////////////////////////////////////////////
double Linear_Eeq::H(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                     const straintensor& Stra, const MaterialParameter& material_parameter)
{
    
    straintensor PF = plastic_flow.PlasticFlowTensor(Stre, Stra, material_parameter);
//    PF.report("PF");

    double L = getLinearFactor(material_parameter);
    
    double hh = PF.equivalent() * L;
    
    return  hh;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& Linear_Eeq::DH_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
    double twoOthree = 2.0/3.0;
    double L = getLinearFactor(material_parameter);

    straintensor PF = plastic_flow.PlasticFlowTensor(Stre, Stra, material_parameter);
    tensor dmOds = plastic_flow.Dm_Ds(Stre, Stra, material_parameter);
    double m_eq =  sqrt(twoOthree * (PF("mn")*PF("mn")).trace());
    tensor tensor1 = PF("ij") * dmOds("ijmn");
    tensor1.null_indices();
    tensor tensor2 = tensor1 * (twoOthree * L);
    
    if (m_eq != 0.0)
       ScalarEvolution::SE_tensorR2 = tensor2 *(1.0/m_eq);
    else 
       ScalarEvolution::SE_tensorR2 = tensor2 *0.0;
    
    return ScalarEvolution::SE_tensorR2;
}

///////////////////////////////////////////////////////////////////////////////
double Linear_Eeq::DH_Diso(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
    double twoOthree = 2.0/3.0;
    double L = getLinearFactor(material_parameter);

    straintensor PF = plastic_flow.PlasticFlowTensor(Stre, Stra, material_parameter);
    tensor dmOdq = plastic_flow.Dm_Diso(Stre, Stra, material_parameter);
    double m_eq =  sqrt(twoOthree * (PF("mn")*PF("mn")).trace());
    double scalar1 = (PF("mn")*dmOdq("mn")).trace();
    double scalar2 = scalar1 * twoOthree * L;
    if (m_eq != 0.0)
       scalar2 *= (1.0/m_eq);
    else 
       scalar2 = 0.0;
    
    return scalar2;
}

/////////////////////////////////////////////////////////////////////////////////
const tensor& Linear_Eeq::DH_Dkin(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
    double twoOthree = 2.0/3.0;
    double L = getLinearFactor(material_parameter);

    straintensor PF = plastic_flow.PlasticFlowTensor(Stre, Stra, material_parameter);
    tensor dmOda = plastic_flow.Dm_Dkin(Stre, Stra, material_parameter);
    double m_eq =  sqrt(twoOthree * (PF("mn")*PF("mn")).trace());
    tensor tensor1 = PF("ij") * dmOda("ijmn");
    tensor1.null_indices();
    tensor tensor2 = tensor1 * (twoOthree * L);
    
    if (m_eq != 0.0)
       ScalarEvolution::SE_tensorR2 = tensor2 *(1.0/m_eq);
    else 
       ScalarEvolution::SE_tensorR2 = tensor2 *0.0;
    
    return ScalarEvolution::SE_tensorR2;
}

///////////////////////////////////////////////////////////////////////////////
double Linear_Eeq::getLinearFactor(const MaterialParameter& material_parameter) const
{
    if ( LinearFactor_index <= material_parameter.getNum_Material_Constant() && LinearFactor_index > 0)
        return material_parameter.getMaterial_Constant(LinearFactor_index-1);
    else {
        opserr << "Linear_Eeq: Invalid Input. " << endln;
        exit (1);
    }
} 

//Guanzhou added for parallel
int Linear_Eeq::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(1);
    idData.Zero();

    idData(0) = LinearFactor_index; 
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "Linear_Eeq::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int Linear_Eeq::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(1);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "Linear_Eeq::recvSelf -- failed to recv ID\n";
	return -1;
    }
    
    LinearFactor_index = idData(0);

    return 0;
}


#endif

