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
// PROGRAMMER:        Zhao Cheng 
// Note:              Helpful discuss with Mahdi Taiebat and Professor Y.F. Dafalias
// DATE:              Fall 2005
// UPDATE HISTORY:    Guanzhou Jie updated for parallel Dec. 2006
//
///////////////////////////////////////////////////////////////////////////////
//

// Ref: Dafalias and Manzari 2004: J. Eng. Mech. 130(6), pp 622-634
// Constants:
// 1- m:      f = [(sij-p*aij)(sij-p*aij)]^0.5 - sqrt(2/3)*p*m = 0;
// 2- alpha:  f = [(sij-p*aij)(sij-p*aij)]^0.5 - sqrt(2/3)*p*m = 0;

#ifndef DM04_YF_CPP
#define DM04_YF_CPP

#include "DM04_YF.h"

stresstensor DM04_YF::DM04st;

//================================================================================
DM04_YF::DM04_YF(int m_which_in, int index_m_in, 
                 int alpha_which_in, int index_alpha_in)
: YieldFunction(YF_TAG_DM04), m_which(m_which_in), index_m(index_m_in), 
  alpha_which(alpha_which_in), index_alpha(index_alpha_in)
{

}

//================================================================================
DM04_YF::~DM04_YF() 
{

}

//================================================================================
YieldFunction* DM04_YF::newObj() 
{
	YieldFunction  *new_YF = new DM04_YF(this->m_which, this->index_m, 
                                         this->alpha_which, this->index_alpha);

	return new_YF;
}

//================================================================================
double DM04_YF::YieldFunctionValue( const stresstensor& Stre, 
                                    const MaterialParameter &MaterialParameter_in ) const
{
	//// f1 = [(sij-p*aij)(sij-p*aij)] - (2/3)*(p*m)^2 = 0
	
    //double p = Stre.p_hydrostatic();
	//double m = getm(MaterialParameter_in);
	//stresstensor alpha = getalpha(MaterialParameter_in);
    //stresstensor s_bar = Stre.deviator() - (alpha * p);
	//double temp1 = ( s_bar("ij") * s_bar("ij") ).trace();
    //return temp1 - (2.0/3.0)*m*m*p*p;

    // f2 = [(sij-p*aij)(sij-p*aij)]^(0.5) - sqrt(2/3)*(p*m) = 0
    
    double m = getm(MaterialParameter_in);
    stresstensor alpha = getalpha(MaterialParameter_in);
    double p = Stre.p_hydrostatic();
    stresstensor s_bar = Stre - alpha*p;
    s_bar = s_bar.deviator();
    double s_norm = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
    
    return s_norm - sqrt(2.0/3.0)*m*p;
}

//================================================================================
const stresstensor& DM04_YF::StressDerivative(const stresstensor& Stre, 
                                              const MaterialParameter &MaterialParameter_in) const
{
     //// f1
     
     //BJtensor KroneckerI("I", 2, def_dim_2);
     //double p = Stre.p_hydrostatic();     
     //double m = getm(MaterialParameter_in);
     //stresstensor alpha = getalpha(MaterialParameter_in);     
	 //stresstensor s_bar = Stre.deviator() - (alpha * p);
     //double s_bar_alpha = (s_bar("ij")*alpha("ij")).trace();     
     //DM04st = s_bar + ( KroneckerI * ( s_bar_alpha/3.0 + (4.0/9.0)*m*m*p ));     
     
     //return DM04st;

     // f2
     
     tensor I2("I", 2, def_dim_2);
     double m = getm(MaterialParameter_in);
     stresstensor alpha = getalpha(MaterialParameter_in);
     double p = Stre.p_hydrostatic();
     stresstensor s_bar = Stre - alpha*p;
     s_bar = s_bar.deviator();
     double s_norm = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
     stresstensor nij;
     if (s_norm != 0.0)
       nij =  s_bar *(1.0/s_norm);
     double a_n = (alpha("ij")*nij("ij")).trace();
     DM04_YF::DM04st = nij + I2 * (m*sqrt(2.0/27.0) + a_n/3.0);
     
     return DM04_YF::DM04st;
}

//================================================================================
const stresstensor& DM04_YF::InTensorDerivative(const stresstensor& Stre, 
                                                const MaterialParameter &MaterialParameter_in, 
                                                int which) const
{
     //// f1
     
     //if (which == index_alpha) {    
     //    stresstensor alpha = getalpha(MaterialParameter_in);
	 //    double p = Stre.p_hydrostatic();	        
     //    stresstensor s_bar = Stre.deviator() - (alpha *p);	    		
	 //	DM04_YF::DM04st = s_bar *(-p);	    
     //}
	 //else {
	 //	opserr << "DM04_YF: Invalid Input. " << endln;
	 //	exit (1);
	 //}    
	 
     //return DM04_YF::DM04st;

     // f2
      
     if (which == index_alpha) {    
       stresstensor alpha = getalpha(MaterialParameter_in);
       double p = Stre.p_hydrostatic();
       stresstensor s_bar = Stre - alpha*p;
       s_bar = s_bar.deviator();
       double s_norm = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
       stresstensor nij;
       if (s_norm != 0.0)
         nij = s_bar *(1.0/s_norm);
       DM04_YF::DM04st = nij * (-p);
     }
	 else {
	 	opserr << "DM04_YF: Invalid Input. " << endln;
	 	exit (1);
	 }    	 
     
     return DM04_YF::DM04st;
}

//================================================================================
int DM04_YF::getNumInternalScalar() const
{
	return 0;
}

//================================================================================
int DM04_YF::getNumInternalTensor() const
{
	return 1;
}

//================================================================================   
int DM04_YF::getYieldFunctionRank() const
{
	return 1;
}

//================================================================================   
double DM04_YF::getm(const MaterialParameter &MaterialParameter_in) const
{
	// to get m
	double m = 0.0;
    if ( m_which == 0 && index_m <= MaterialParameter_in.getNum_Material_Constant() && index_m > 0 ) {
        m = MaterialParameter_in.getMaterial_Constant(index_m-1);
	    if (m <= 0.0) {
		    opserr << "DM04_YF: Invalid Input, m <= 0.0. " << endln;
		    exit (1);
	    }
		return m;
    }
	else {
		opserr << "Warning!! DM04_YF: Invalid Input. " << endln;
		exit (1);
	}
}

//================================================================================ 
const stresstensor& DM04_YF::getalpha(const MaterialParameter &MaterialParameter_in) const
{
	//to get alpha
	if ( alpha_which == 2 && index_alpha <= MaterialParameter_in.getNum_Internal_Tensor() && index_alpha > 0 ) {
		DM04_YF::DM04st = MaterialParameter_in.getInternal_Tensor(index_alpha-1);
		return DM04_YF::DM04st;
	}
	else {
		opserr << "Warning!! DM04_YF: Invalid Input. " << endln;
		exit (1);
	}
}

//Guanzhou added for parallel
int DM04_YF::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(4);
    idData.Zero();

    idData(0) = m_which;    
    idData(1) = index_m;    
    idData(2) = alpha_which;
    idData(3) = index_alpha;
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "DM04_YF::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int DM04_YF::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(4);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "DM04_YF::recvSelf -- failed to recv ID\n";
	return -1;
    }

    m_which      = idData(0);
    index_m      = idData(1);
    alpha_which = idData(2);
    index_alpha = idData(3);

    return 0;
}

#endif

