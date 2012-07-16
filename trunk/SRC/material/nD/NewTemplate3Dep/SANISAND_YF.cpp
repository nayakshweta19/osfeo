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
// DESIGNER:          Mahdi Taiebat, Boris Jeremic
// PROGRAMMER:        Mahdi Taiebat, Boris Jeremic 
// Note:              
// DATE:              Spring 2007
// UPDATE HISTORY:    
//
///////////////////////////////////////////////////////////////////////////////
//

// This is based on:
// f = (3/2) * [(sij-p*aij)(sij-p*aij)] - m^2 * p^2 * [1-(p/p0)^n] = 0;

// Parameters:
// 1- m
// 2- p0
// 3- alpha

// Note: n=20 is used by default

#ifndef SANISAND_YF_CPP
#define SANISAND_YF_CPP

#include "SANISAND_YF.h"

stresstensor SANISAND_YF::SANISANDst;

// Note: n=20 is used by default
static double n = 20.0; // this is the exponent to yield function


//================================================================================
SANISAND_YF::SANISAND_YF(int m_which_in,      int index_m_in, 
                         int p0_which_in,     int index_p0_in,
                         int alpha_which_in,  int index_alpha_in)
: YieldFunction(YF_TAG_SANISAND), 
m_which(m_which_in), index_m(index_m_in), 
p0_which(p0_which_in), index_p0(index_p0_in), 
alpha_which(alpha_which_in), index_alpha(index_alpha_in)
{

}

//================================================================================
SANISAND_YF::~SANISAND_YF() 
{

}

//================================================================================
YieldFunction* SANISAND_YF::newObj() 
{
    if ( this->m_which < 0) 
      {
        ::fprintf(stderr,"\n\n SANISAND_YF: this->m_which <= 0.0 \n\n\n ");
        exit(1);
      }

// CONTINUE with checks

    YieldFunction  *new_YF = new SANISAND_YF(this->m_which,     this->index_m, 
                                             this->p0_which,    this->index_p0,
                                             this->alpha_which, this->index_alpha);
    return new_YF;
}

//================================================================================
double SANISAND_YF::YieldFunctionValue( const stresstensor& Stre, 
                                        const MaterialParameter &MaterialParameter_in ) const
{
//    double n = 20.0;
    double p = Stre.p_hydrostatic();
    double m = getm(MaterialParameter_in);
    double p0 = getp0(MaterialParameter_in);
    stresstensor alpha = getalpha(MaterialParameter_in);
    stresstensor s_bar = Stre.deviator() - (alpha * p);
    double temp1 = 0.0;
    if ( p>0 ) 
      temp1 = pow(p/p0,n); 
    double temp2 = 1.5 * ( s_bar("ij") * s_bar("ij") ).trace();
    double temp3 = m*m*p*p*(1-temp1);
    return temp2 - temp3;
}

//================================================================================
const stresstensor& SANISAND_YF::StressDerivative(const stresstensor& Stre, 
                                                  const MaterialParameter &MaterialParameter_in) const
{
//    double n = 20.0;
    double p = Stre.p_hydrostatic();     
    double m = getm(MaterialParameter_in);
    double p0 = getp0(MaterialParameter_in);
    stresstensor alpha = getalpha(MaterialParameter_in);     
    stresstensor s_bar = Stre.deviator() - (alpha * p);
    stresstensor dfods = s_bar*(3.0);
    double temp1 = 0.0;
    if ( p>0 ) 
      temp1 = pow(p/p0,n); 
    double dfodp = (alpha("mn")*s_bar("mn")).trace() *(-3.0) - 2.0*m*m*p + (2.0+n)*m*m*p*temp1;
    BJtensor I2("I", 2, def_dim_2);
    SANISAND_YF::SANISANDst = dfods - I2*(1.0/3.0)*dfodp;
    return SANISAND_YF::SANISANDst;
}

//================================================================================
double SANISAND_YF::InScalarDerivative(const stresstensor& Stre, 
                                       const MaterialParameter &MaterialParameter_in, 
                                       int which) const
{
    double p = Stre.p_hydrostatic();     
    double m = getm(MaterialParameter_in);
    double p0 = getp0(MaterialParameter_in);
    double temp1 = (-n/p0) * m*m*p*p* pow(p/p0,n); 	    
    return temp1;
}

//================================================================================
const stresstensor& SANISAND_YF::InTensorDerivative(const stresstensor& Stre, 
                                                    const MaterialParameter &MaterialParameter_in, 
                                                    int which) const
{
    stresstensor alpha = getalpha(MaterialParameter_in);
    double p = Stre.p_hydrostatic();	        
    stresstensor s_bar = Stre.deviator() - (alpha * p);	    		
    SANISAND_YF::SANISANDst = s_bar * (-3.0) * p;	    
    return SANISAND_YF::SANISANDst;
}

//================================================================================
int SANISAND_YF::getNumInternalScalar() const
{
    return 1;
}

//================================================================================
int SANISAND_YF::getNumInternalTensor() const
{
    return 1;
}

//================================================================================   
int SANISAND_YF::getYieldFunctionRank() const
{
    return 1;
}

//================================================================================
// to get m
double SANISAND_YF::getm(const MaterialParameter &MaterialParameter_in) const
{
    double m = 0.0;
    m = MaterialParameter_in.getMaterial_Constant(index_m-1);
    return m;
}

//================================================================================
// to get p0
double SANISAND_YF::getp0(const MaterialParameter &MaterialParameter_in) const
{
    double p0 = 0.0;
    p0 = MaterialParameter_in.getInternal_Scalar(index_p0-1);
    return p0;
}

//================================================================================
//to get alpha
const stresstensor& SANISAND_YF::getalpha(const MaterialParameter &MaterialParameter_in) const
{
    SANISAND_YF::SANISANDst = MaterialParameter_in.getInternal_Tensor(index_alpha-1);
    return SANISAND_YF::SANISANDst;
}

//Guanzhou added for parallel
int SANISAND_YF::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(6);
    idData.Zero();

    idData(0) = m_which;    
    idData(1) = index_m;    
    idData(2) = p0_which;    
    idData(3) = index_p0;    
    idData(4) = alpha_which;
    idData(5) = index_alpha;
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) 
    {
   	opserr << "SANISAND_YF::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int SANISAND_YF::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(6);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) 
    {
    	opserr << "SANISAND_YF::recvSelf -- failed to recv ID\n";
	return -1;
    }

    m_which     = idData(0);
    index_m     = idData(1);
    p0_which    = idData(2);
    index_p0    = idData(3);
    alpha_which = idData(4);
    index_alpha = idData(5);

    return 0;
}

#endif

