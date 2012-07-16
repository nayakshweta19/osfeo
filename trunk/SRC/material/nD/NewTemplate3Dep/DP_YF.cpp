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

#ifndef DP_YF_CPP
#define DP_YF_CPP

#include "DP_YF.h"
#define YF_TAG_DP 124003

#include <iostream>
using namespace std;

stresstensor DP_YF::DPst;
stresstensor DP_YF::DPa;

//================================================================================
DP_YF::DP_YF(int M_which_in, int index_M_in, 
             int alpha_which_in, int index_alpha_in)
: YieldFunction(YF_TAG_DP), M_which(M_which_in), index_M(index_M_in), 
  alpha_which(alpha_which_in), index_alpha(index_alpha_in)
{

}

//================================================================================
DP_YF::~DP_YF() 
{

}

//================================================================================
YieldFunction* DP_YF::newObj() 
{
	YieldFunction  *new_YF = new DP_YF(this->M_which, this->index_M, 
	                                   this->alpha_which, this->index_alpha);

	return new_YF;
}

//================================================================================
double DP_YF::YieldFunctionValue( const stresstensor& Stre, 
                                 const MaterialParameter &MaterialParameter_in ) const
{

// f = [(sij-p*aij)(sij-p*aij)]^0.5 - sqrt(2/3)*(p*m) = 0


    double M = getM(MaterialParameter_in);
    double p = Stre.p_hydrostatic();
    
    stresstensor alpha;

    if (alpha_which == 2) 
       { 
         alpha = getalpha(MaterialParameter_in);
       }

    stresstensor s_bar = Stre - alpha*p;
    s_bar = s_bar.deviator();
 
    double s_norm = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
    double YF = s_norm - sqrt(2.0/3.0)*M*p;

    
// Nima
/*
    // g = 1.5*(sij-p*aij)(sij-p*aij) - (m*p)^2 = 0
    double M = getM(MaterialParameter_in);
    stresstensor alpha = getalpha(MaterialParameter_in);
    double p = Stre.p_hydrostatic();
    stresstensor s_bar = Stre - alpha*p;
    s_bar = s_bar.deviator();
    double s_norm = (s_bar("ij")*s_bar("ij")).trace();
    double YF = s_norm * (3.0/2.0) + - (M * M * p * p);
*/    
    return YF;
}

//================================================================================
const stresstensor& DP_YF::StressDerivative(const stresstensor& Stre, 
                                            const MaterialParameter &MaterialParameter_in) const
{    

    tensor I2("I", 2, def_dim_2);
    double M = getM(MaterialParameter_in);
    double p = Stre.p_hydrostatic();

    stresstensor alpha;

    if (alpha_which == 2) 
       { 
         alpha = getalpha(MaterialParameter_in);
       }

    stresstensor s_bar = Stre - alpha*p;
    s_bar = s_bar.deviator();

    double s_norm = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
    stresstensor nij;
    if (s_norm != 0.0)
      nij =  s_bar *(1.0/s_norm);
    double a_n = (alpha("ij")*nij("ij")).trace();

    DP_YF::DPst = nij + I2 * (M*sqrt(2.0/27.0) + a_n/3.0);

    return DP_YF::DPst;

// Nima
/*
    tensor I2("I", 2, def_dim_2);
    double M = getM(MaterialParameter_in);
    double p = Stre.p_hydrostatic();
    stresstensor alpha = getalpha(MaterialParameter_in);
    stresstensor s_bar;
    if (alpha_which == 2) 
      { 
        stresstensor alpha = getalpha(MaterialParameter_in);
        s_bar = Stre - (alpha * p);
      }
    else
      {
        s_bar = Stre;
      }

    s_bar = s_bar.deviator();
    double scalar1 = (alpha("ij")*s_bar("ij")).trace();

    DP_YF::DPst = s_bar *3.0 + I2 *(scalar1 + M*M*p*2.0/3.0);
*/


}

//================================================================================
double DP_YF::InScalarDerivative(const stresstensor& Stre, 
                                 const MaterialParameter &MaterialParameter_in, 
                                 int which) const
{

    double scalar1 = 0.0;

	if (M_which == 1 && which == 1) {
      double p = Stre.p_hydrostatic();
      scalar1 = -p * sqrt(2.0/3.0);
	}


// Nima
/*
    double scalar1 = 0.0;

	if (M_which == 1 && which == 1) {
      double p = Stre.p_hydrostatic();
      double M = getM(MaterialParameter_in);
      scalar1 = -2 * p * p * M;
	}
*/

	return scalar1;
}

//================================================================================
const stresstensor& DP_YF::InTensorDerivative(const stresstensor& Stre, 
                                              const MaterialParameter &MaterialParameter_in, 
                                              int which) const
{


    stresstensor nij;
    double p;

    if (alpha_which == 2 || which == 1) {
      stresstensor alpha = getalpha(MaterialParameter_in);
      p = Stre.p_hydrostatic();
      stresstensor s_bar = Stre - (alpha * p);
      s_bar = s_bar.deviator();
      double s_norm = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
      if (s_norm != 0.0)
        nij = s_bar *(1.0/s_norm); 
    }

    DP_YF::DPst = nij *(-p);


// Nima
/*
    stresstensor s_bar;
    double p;

    if (alpha_which == 2 || which == 1)
    {
      stresstensor alpha = getalpha(MaterialParameter_in);
      p = Stre.p_hydrostatic();
      stresstensor s_bar = Stre - (alpha * p);
      s_bar = s_bar.deviator();
    }
    DP_YF::DPst = s_bar *(-3*p);
*/

    return DP_YF::DPst;
}

//================================================================================
int DP_YF::getNumInternalScalar() const
{
	if (M_which == 1)
		return 1;
	else
		return 0;
}

//================================================================================
int DP_YF::getNumInternalTensor() const
{
	if (alpha_which == 2)
		return 1;
	else
		return 0;
}

//================================================================================   
int DP_YF::getYieldFunctionRank() const
{
	return 1;
}

//================================================================================   
double DP_YF::getM(const MaterialParameter &MaterialParameter_in) const
{
	// to get M
	if ( M_which == 0) {
		if ( index_M <= MaterialParameter_in.getNum_Material_Constant() && index_M > 0)
			return MaterialParameter_in.getMaterial_Constant(index_M-1); 
		else {
			opserr << "DP_YF: Invalid Input. " << endln;
			exit (1);
		}
	}
	else if ( M_which == 1) {
		if ( index_M <= MaterialParameter_in.getNum_Internal_Scalar() && index_M > 0)
			return MaterialParameter_in.getInternal_Scalar(index_M-1); 
		else {
			opserr << "DP_YF: Invalid Input. " << endln;
			exit (1);
		}
    }
	else {
		opserr << "DP_YF: Invalid Input. " << endln;
		exit(1);
	}
}

//================================================================================ 
const stresstensor& DP_YF::getalpha(const MaterialParameter &MaterialParameter_in) const
{
	//to get alpha(backstress)
	if ( alpha_which == 2 && index_alpha <= MaterialParameter_in.getNum_Internal_Tensor() && index_alpha > 0) {
		DP_YF::DPa = MaterialParameter_in.getInternal_Tensor(index_alpha-1);
		return DP_YF::DPa;
	}
	else if ( alpha_which == -1 ) {
        DP_YF::DPa = DP_YF::DPa*0.0;
        return DP_YF::DPa;
    }    
    else {
		opserr << "DP_YF: Invalid Input. " << endln;
		exit (1);
	}
}

//Guanzhou added for parallel
int DP_YF::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(4);
    idData.Zero();

    idData(0) = M_which;    
    idData(1) = index_M;    
    idData(2) = alpha_which;
    idData(3) = index_alpha;
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "DP_YF::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int DP_YF::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(6);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "DP_YF::recvSelf -- failed to recv ID\n";
	return -1;
    }

    M_which       = idData(0);
    index_M       = idData(1);
    alpha_which     = idData(2);
    index_alpha     = idData(3);
    
    return 0;
}

#endif

