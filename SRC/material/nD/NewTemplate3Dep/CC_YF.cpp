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
//                    Mahdi and Boris correcting derivatives 
//                    as per Mahdi's derivations, 27April2007
//                    Nima Tafazzoli updated for API (Feb 2009) 
//                    Nima Tafazzoli updated for InScalarDerivative  (Feb 2010)
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef CC_YF_CPP
#define CC_YF_CPP

#include "CC_YF.h"
#define YF_TAG_CC 124001

stresstensor CC_YF::CCst;

//================================================================================
CC_YF::CC_YF(int M_which_in, int index_M_in, 
             int p0_which_in, int index_p0_in)
: YieldFunction(YF_TAG_CC), M_which(M_which_in), index_M(index_M_in), 
  p0_which(p0_which_in), index_p0(index_p0_in)
{

}

//================================================================================
CC_YF::~CC_YF() 
{

}

//================================================================================
YieldFunction* CC_YF::newObj() 
{
	YieldFunction  *new_YF = new CC_YF(this->M_which, 
                                       this->index_M, 
                                       this->p0_which, 
                                       this->index_p0);

	return new_YF;
}

//================================================================================
double CC_YF::YieldFunctionValue( const stresstensor& Stre, 
                                  const MaterialParameter &MaterialParameter_in ) const
{
	// f = q*q - M*M*p*(po - p) = 0
	
	double M = getM(MaterialParameter_in);
	double p0 = getP0(MaterialParameter_in);
	double p = Stre.p_hydrostatic();
	double q = Stre.q_deviatoric();

	return  q*q - M*M*p*(p0 - p);

}

//================================================================================
const stresstensor& CC_YF::StressDerivative(const stresstensor& Stre, 
                                            const MaterialParameter &MaterialParameter_in) const
{
    tensor I2("I", 2, def_dim_2);
    stresstensor sij = Stre.deviator();
    
    double M = getM(MaterialParameter_in);
    double p0 = getP0(MaterialParameter_in);
    double p = Stre.p_hydrostatic();

    double scalar1 = M*M*(p0-2.0*p)*(1.0/3.0);

    CC_YF::CCst = (sij * 3.0) + (I2 * scalar1);

    return CC_YF::CCst;
}

//================================================================================
double CC_YF::InScalarDerivative(const stresstensor& Stre, 
                                 const MaterialParameter &MaterialParameter_in, 
                                 int which) const
{
	if (p0_which == 1 && which == index_p0) {
		double M = getM(MaterialParameter_in);
		double p = Stre.p_hydrostatic();
//		return  (-1.0)*M*M*p;
		return M*M*p; //Nima Tafazzoli removed the (-) sign

	}
	else {
		opserr << "Warning!! CC_YF: Invalid Input Parameter. " << endln;
		exit (1);
	}
}

//================================================================================
int CC_YF::getNumInternalScalar() const
{
	return 1;
}

//================================================================================
int CC_YF::getNumInternalTensor() const
{
	return 0;
}

//================================================================================   
int CC_YF::getYieldFunctionRank() const
{
	return 2;
}

//================================================================================   
double CC_YF::getM(const MaterialParameter &MaterialParameter_in) const
{
	// to get M
	if ( M_which == 0 && index_M <= MaterialParameter_in.getNum_Material_Constant() && index_M > 0 )
		return MaterialParameter_in.getMaterial_Constant(index_M-1);
	else {
		opserr << "Warning!! CC_YF: Invalid Input (M). " << endln;
		exit (1);
	}
}

//================================================================================ 
double CC_YF::getP0(const MaterialParameter &MaterialParameter_in) const
{
	//to get P0
	if ( p0_which == 1 && index_p0 <= MaterialParameter_in.getNum_Internal_Scalar() && index_p0 > 0 )
		return MaterialParameter_in.getInternal_Scalar(index_p0-1);
	else {
		opserr << "Warning!! CC_YF: Invalid Input (po). " << endln;
		exit (1);
	}
}

//Guanzhou added for parallel
int CC_YF::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(4);
    idData.Zero();

    idData(0) = M_which; 
    idData(1) = index_M; 
    idData(2) = p0_which;
    idData(3) = index_p0;
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "CC_YF::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int CC_YF::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(4);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "CC_YF::recvSelf -- failed to recv ID\n";
	return -1;
    }

    M_which    = idData(0);
    index_M    = idData(1);
    p0_which = idData(2);
    index_p0 = idData(3);

    return 0;
}

#endif

