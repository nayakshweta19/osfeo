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
//
//                    Mahdi and Boris correcting derivatives 
//                    as per Mahdi's derivations, 27April2007
//                        
///////////////////////////////////////////////////////////////////////////////
//

#ifndef CC_Ev_CPP
#define CC_Ev_CPP

#include "CC_Ev.h"

CC_Ev::CC_Ev(int M_index_in,
             int lambda_index_in, 
             int kappa_index_in,
             int e0_index_in,
             int p0_index_in)
: ScalarEvolution(SE_TAG_CC), M_index(M_index_in),
  lambda_index(lambda_index_in), 
  kappa_index(kappa_index_in),
  e0_index(e0_index_in),
  p0_index(p0_index_in)
{

}

ScalarEvolution* CC_Ev::newObj()
{
    ScalarEvolution* nObj = new CC_Ev(this->M_index,
                                      this->lambda_index,
                                      this->kappa_index,
                                      this->e0_index,
                                      this->p0_index);
    return nObj;
}

///////////////////////////////////////////////////////////////////////////////
double CC_Ev::H(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                const straintensor& Stra, const MaterialParameter& material_parameter)
{    
// this is \bar{p_0}
    // Make sure f = Q = q*q - M*M*p*(po - p) = 0
    
    double p0 = getp0(material_parameter);
    double lambda = getlambda(material_parameter);
    double kappa = getkappa(material_parameter);
    double e0 = gete0(material_parameter);
    
    double e = e0 + (1.0 + e0) *Stra.Iinvariant1();

    //// One way
    //straintensor PF = plastic_flow.PlasticFlowTensor(Stre, Stra, material_parameter);
    //double d_Ev = - PF.Iinvariant1(); // note "minus"

    // Another way
    double M = getM(material_parameter);
    double p = Stre.p_hydrostatic();
//    double d_Ev = -M*M*(p0-2.0*p); // Zhao's version
//    double d_Ev = (-1.0)*M*M*(p0 - 2.0*p); // Mahdi and Boris 27April2007
   double d_Ev = M*M*(2.0*p - p0);
   
    return (1.0 + e) * p0 * d_Ev / (lambda - kappa);
}

///////////////////////////////////////////////////////////////////////////////
const tensor& CC_Ev::DH_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
// this is d \bar{p_0} / d sigma_ij

    tensor I2("I", 2, def_dim_2);

    double M = getM(material_parameter);
    double p0 = getp0(material_parameter);
    double lambda = getlambda(material_parameter);
    double kappa = getkappa(material_parameter);
    double e0 = gete0(material_parameter);
    
    double e = e0 + (1.0 + e0) *Stra.Iinvariant1();

    double scalar1 = (1.0+e)*p0*M*M*(-2.0/3.0)/(lambda-kappa);

    ScalarEvolution::SE_tensorR2 = I2 * scalar1;
    
    return ScalarEvolution::SE_tensorR2;
}

///////////////////////////////////////////////////////////////////////////////
double CC_Ev::DH_Diso(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
// this is d \bar{p_0} / d p_0


    double M = getM(material_parameter);
    double p0 = getp0(material_parameter);
    double lambda = getlambda(material_parameter);
    double kappa = getkappa(material_parameter);
    double e0 = gete0(material_parameter);
    
    double e = e0 + (1.0 + e0) *Stra.Iinvariant1();
    double p = Stre.p_hydrostatic();

//    double scalar1 = (1.0+e)*p0*M*M*(p0-p)/(lambda-kappa); // Zhao's version
    double scalar1 = (-2.0)*(1.0+e)*p0*M*M*(p0-p)/(lambda-kappa); // Mahdi
    
    return scalar1;
}

// Get M
double CC_Ev::getM(const MaterialParameter& material_parameter) const
{
    if ( M_index <= material_parameter.getNum_Material_Constant() && lambda_index > 0)
        return material_parameter.getMaterial_Constant(M_index-1);
    else {
        opserr << "CC_Ev: Invalid Input of " << lambda_index << endln;
        exit (1);
    }
}

// Get lambda
double CC_Ev::getlambda(const MaterialParameter& material_parameter) const
{
    if ( lambda_index <= material_parameter.getNum_Material_Constant() && lambda_index > 0)
        return material_parameter.getMaterial_Constant(lambda_index-1);
    else {
        opserr << "CC_Ev: Invalid Input of " << lambda_index << endln;
        exit (1);
    }
}

// Get kappa
double CC_Ev::getkappa(const MaterialParameter& material_parameter) const
{
    if ( kappa_index <= material_parameter.getNum_Material_Constant() && kappa_index > 0)
        return material_parameter.getMaterial_Constant(kappa_index-1);
    else {
        opserr << "CC_Ev: Invalid Input of " << kappa_index << endln;
        exit (1);
    }
}

// Get e0
double CC_Ev::gete0(const MaterialParameter& material_parameter) const
{
    if ( e0_index <= material_parameter.getNum_Material_Constant() && e0_index > 0)
        return material_parameter.getMaterial_Constant(e0_index-1);
    else {
        opserr << "CC_Ev: Invalid Input of " << e0_index << endln;
        exit (1);
    }
}

// Get p0
double CC_Ev::getp0(const MaterialParameter& material_parameter) const
{
    if ( p0_index <= material_parameter.getNum_Internal_Scalar() && p0_index > 0)
        return material_parameter.getInternal_Scalar(p0_index-1);
    else {
        opserr << "CC_Ev: Invalid Input of " << p0_index << endln;
        exit (1);
    }
}

//Guanzhou added for parallel
int CC_Ev::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(5);
    idData.Zero();

    idData(0) = lambda_index;
    idData(1) = kappa_index;
    idData(2) = e0_index;
    idData(3) = p0_index;
    idData(4) = M_index;
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "CC_Ev::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int CC_Ev::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(5);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "CC_Ev::recvSelf -- failed to recv ID\n";
	return -1;
    }

    lambda_index   = idData(0);
    kappa_index    = idData(1);
    e0_index = idData(2);
    p0_index = idData(3);
    M_index = idData(4);   
    
    return 0;
}

#endif

