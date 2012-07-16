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
//
//                    Mahdi and Boris correcting derivatives 
//                    as per Mahdi's derivations, 27April2007
//
//                    Nima Tafazzoli updated for PlasticFlowTensor  (Feb 2010)
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef CC_PF_CPP
#define CC_PF_CPP

#include "CC_PF.h"
#include <OPS_Globals.h>
#define PF_TAG_CC 120002
straintensor CC_PF::CCm;

//================================================================================
CC_PF::CC_PF(int M_which_in, int index_M_in, 
             int p0_which_in, int index_p0_in)
: PlasticFlow(PF_TAG_CC), M_which(M_which_in), index_M(index_M_in),
  p0_which(p0_which_in), index_p0(index_p0_in)
{

}

//================================================================================
CC_PF::~CC_PF() 
{  

}

//================================================================================
PlasticFlow* CC_PF::newObj() 
{  
     PlasticFlow  *new_PF = new CC_PF(this->M_which, 
                                      this->index_M, 
                                      this->p0_which, 
                                      this->index_p0);
     
     return new_PF;
}

//================================================================================
const straintensor& CC_PF::PlasticFlowTensor(const stresstensor &Stre, 
                                             const straintensor &Stra, 
                                             const MaterialParameter &MaterialParameter_in) const
{
	// Q = q*q - M*M*p*(po - p) = 0

    tensor I2("I", 2, def_dim_2);
    stresstensor sij = Stre.deviator();
    
    double M = getM(MaterialParameter_in);
    double p0 = getP0(MaterialParameter_in);
    double p = Stre.p_hydrostatic();
//	double p0 = (-1.0)*getP0(MaterialParameter_in);//Mahdi - II
//	double p = (-1.0)*Stre.p_hydrostatic();        //Mahdi - II

    double scalar1 = M*M*(p0-2.0*p)*(1.0/3.0);
    	
//   	CC_PF::CCm = (I2 * scalar1) + (sij * 3.0); // Zhao's version
   	CC_PF::CCm =  (sij * 3.0) + (I2 * scalar1); // Nima
//   	CC_PF::CCm =  (sij * 3.0) - (I2 * scalar1); // Mahdi and Boris 27April2007 --- II
    
    return CC_PF::CCm;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& CC_PF::Dm_Ds(const stresstensor &Stre, 
                           const straintensor &Stra, 
                           const MaterialParameter &MaterialParameter_in) const 
{
// out 27April2007 Mahdi and Boris   Zhao's version has error
// out 27April2007 Mahdi and Boris
// out 27April2007 Mahdi and Boris     tensor I2("I", 2, def_dim_2);
// out 27April2007 Mahdi and Boris      tensor I_ijkl = I2("ij")*I2("kl");
// out 27April2007 Mahdi and Boris      I_ijkl.null_indices();
// out 27April2007 Mahdi and Boris      tensor I_ikjl = I_ijkl.transpose0110();
// out 27April2007 Mahdi and Boris      tensor I_iljk = I_ijkl.transpose0111();
// out 27April2007 Mahdi and Boris      tensor I4s = (I_ikjl+I_iljk)*0.5;
// out 27April2007 Mahdi and Boris  
// out 27April2007 Mahdi and Boris      double M = getM(MaterialParameter_in);
// out 27April2007 Mahdi and Boris  
// out 27April2007 Mahdi and Boris      PlasticFlow::PF_tensorR4 = I4s *3.0 - I_ijkl *(1.0 + 2.0*M*M/9.0);


// Mahdi and Boris 27April2007
     tensor delta("I", 2, def_dim_2);
     tensor temp2 = delta("ij") * delta("mn");
     tensor temp3 = temp2.transpose_1234_1324(); // this in now delta_im delta_jn
 
 
     double M = getM(MaterialParameter_in);
 
     PlasticFlow::PF_tensorR4 = temp3 *3.0 + temp2 *( 2.0*M*M/9.0 - 1.0 );
    
     return PlasticFlow::PF_tensorR4;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& CC_PF::Dm_Diso(const stresstensor &Stre, 
                             const straintensor &Stra, 
                             const MaterialParameter &MaterialParameter_in) const 
{
    tensor I2("I", 2, def_dim_2);

    double M = getM(MaterialParameter_in);
    
//    PlasticFlow::PF_tensorR2 = I2 *(M*M/3.0); //Zhao's version
    PlasticFlow::PF_tensorR2 = I2 *(M * M * 0.333333333333); // Mahdi and Boris 27April2007
        
    return PlasticFlow::PF_tensorR2;
}

//================================================================================   
double CC_PF::getM(const MaterialParameter &MaterialParameter_in) const
{
	// to get M
	if ( M_which == 0 && index_M <= MaterialParameter_in.getNum_Material_Constant() && index_M > 0)
		return MaterialParameter_in.getMaterial_Constant(index_M-1);
	else {
		opserr << "Warning!! CC_PF: Invalid Input (M). " << endln;
		exit (1);
	}
}

//================================================================================ 
double CC_PF::getP0(const MaterialParameter &MaterialParameter_in) const
{
	//to get P0
	if ( p0_which == 1 && index_p0 <= MaterialParameter_in.getNum_Internal_Scalar() && index_p0 > 0)
		return MaterialParameter_in.getInternal_Scalar(index_p0-1);
	else {
		opserr << "Warning!! CC_PF: Invalid Input (po). " << endln;
		exit (1);
	}
}

//Guanzhou added for parallel
int CC_PF::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(4);
    idData.Zero();

    idData(0) = M_which; 
    idData(1) = index_M; 
    idData(2) = p0_which;
    idData(3) = index_p0;
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "CC_PF::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int CC_PF::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(4);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "CC_PF::recvSelf -- failed to recv ID\n";
	return -1;
    }

    M_which    = idData(0);
    index_M    = idData(1);
    p0_which = idData(2);
    index_p0 = idData(3);

    return 0;
}

#endif

