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
// UPDATE HISTORY:    Guanzhou Jie updated for parallel Dec 2006s
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef RMC_YF_CPP
#define RMC_YF_CPP

#include "RMC_YF.h"
#define YF_TAG_RMC 124004

stresstensor RMC_YF::RMCst;

//================================================================================
RMC_YF::RMC_YF(int a_which_in, int index_a_in, 
             int k_which_in, int index_k_in, 
             int r_which_in, int index_r_in)
: YieldFunction(YF_TAG_RMC), a_which(a_which_in), index_a(index_a_in), 
  k_which(k_which_in), index_k(index_k_in),
  r_which(r_which_in), index_r(index_r_in)
{

}

//================================================================================
RMC_YF::~RMC_YF() 
{

}

//================================================================================
YieldFunction* RMC_YF::newObj() 
{
	YieldFunction  *new_YF = new RMC_YF(this->a_which, this->index_a, 
	                                   this->k_which, this->index_k, 
	                                   this->r_which, this->index_r);

	return new_YF;
}

//================================================================================
double RMC_YF::YieldFunctionValue( const stresstensor& Stre, 
                                 const MaterialParameter &MaterialParameter_in ) const
{
	// f = a*I1 + (J2D)^0.5/g(theta) - k = 0; or
	// f = -3*a*p + q/(sqrt(3)*g(theta) - k = 0;
	
	double a = geta(MaterialParameter_in);
	double k = getk(MaterialParameter_in);
	double r = getr(MaterialParameter_in);
	
	double theta = Stre.theta();
	double g = RoundedFunction(theta, r);
	
	return Stre.Iinvariant1() *a + sqrt(Stre.Jinvariant2())/g - k;
}

//================================================================================
const stresstensor& RMC_YF::StressDerivative(const stresstensor& Stre, 
                                            const MaterialParameter &MaterialParameter_in) const
{
	double a = geta(MaterialParameter_in);
	double r = getr(MaterialParameter_in);
	
	double q = Stre.q_deviatoric();
	double theta = Stre.theta();
	
	// g = f1/f2, 1/g = f2/f1;
	// d(1/g) = d(f2/f1) = (df2*f1 - df1*f2)/(f1)^2;
	double f1 = RoundedFunctionf1(theta, r);
	double f2 = RoundedFunctionf2(theta, r);
	double df1 = RoundedFunctiondf1(theta, r);
	double df2 = RoundedFunctiondf2(theta, r);
	double dginv = (df2*f1 - df1*f2) /(f1*f1);
	
	double dfodp = -3.0*a;
	double dfodq = f1/(f2*sqrt(3.0));
	double dfodtheta = dginv*q/sqrt(3.0);
	
	RMCst = Stre.dpoverds() *dfodp + Stre.dqoverds() *dfodq + Stre.dthetaoverds() *dfodtheta;

	return RMCst;
}

//================================================================================
double RMC_YF::InScalarDerivative(const stresstensor& Stre, 
                                 const MaterialParameter &MaterialParameter_in, 
                                 int which) const
{
	if (a_which == 1 && which == index_a)
		return Stre.Iinvariant1();
	
	if (k_which == 1 && which == index_k)
		return -1.0;
		
	return 0.0;
}

//================================================================================
int RMC_YF::getNumInternalScalar() const
{
	int Numyf = 0;
	
	if ( a_which == 1)
		Numyf++;
	
	if ( k_which == 1)
		Numyf++;
	
	return Numyf;
}

//================================================================================
int RMC_YF::getNumInternalTensor() const
{
	return 0;
}

//================================================================================   
int RMC_YF::getYieldFunctionRank() const
{
	return 1;
}

//================================================================================   
double RMC_YF::geta(const MaterialParameter &MaterialParameter_in) const
{
	// to get a
	if ( a_which == 0) {
		if ( index_a <= MaterialParameter_in.getNum_Material_Constant() && index_a > 0)
			return MaterialParameter_in.getMaterial_Constant(index_a-1); 
		else {
			opserr << "RMC_YF: Invalid Input. " << endln;
			exit (1);
		}
	}
	else if ( a_which == 1) {
		if ( index_a <= MaterialParameter_in.getNum_Internal_Scalar() && index_a > 0)
			return MaterialParameter_in.getInternal_Scalar(index_a-1); 
		else {
			opserr << "RMC_YF: Invalid Input. " << endln;
			exit (1);
		}
    }
	else {
		opserr << "RMC_YF: Invalid Input. " << endln;
		exit(1);
	}
}

//================================================================================   
double RMC_YF::getk(const MaterialParameter &MaterialParameter_in) const
{
	// to get k
	if ( k_which == 0) {
		if ( index_k <= MaterialParameter_in.getNum_Material_Constant() && index_k > 0)
			return MaterialParameter_in.getMaterial_Constant(index_k-1); 
		else {
			opserr << "RMC_YF: Invalid Input. " << endln;
			exit (1);
		}
	}
	else if ( k_which == 1) {
		if ( index_k <= MaterialParameter_in.getNum_Internal_Scalar() && index_k > 0)
			return MaterialParameter_in.getInternal_Scalar(index_k-1); 
		else {
			opserr << "RMC_YF: Invalid Input. " << endln;
			exit (1);
		}
    }
	else {
		opserr << "RMC_YF: Invalid Input. " << endln;
		exit(1);
	}
}

//================================================================================   
double RMC_YF::getr(const MaterialParameter &MaterialParameter_in) const
{
	// to get r
	if ( r_which == 0) {
		if ( index_r <= MaterialParameter_in.getNum_Material_Constant() && index_r > 0)
			return MaterialParameter_in.getMaterial_Constant(index_r-1); 
		else {
			opserr << "RMC_YF: Invalid Input. " << endln;
			exit (1);
		}
	}
	else if ( r_which == 1) {
		if ( index_r <= MaterialParameter_in.getNum_Internal_Scalar() && index_r > 0)
			return MaterialParameter_in.getInternal_Scalar(index_r-1); 
		else {
			opserr << "RMC_YF: Invalid Input. " << endln;
			exit (1);
		}
    }
	else {
		opserr << "RMC_YF: Invalid Input. " << endln;
		exit(1);
	}
}

//================================================================================   
double RMC_YF::RoundedFunctionf1(double s, double r) const
{
	double t1 = sqrt(-4.0*r + 5.0*r*r + 4.0*(1.0 - r*r)*pow(cos(s),2));
	return 2.0*(1.0 - r*r)*cos(s) + (-1.0 + 2.0*r) * t1;
}

//================================================================================   
double RMC_YF::RoundedFunctionf2(double s, double r) const
{
	return pow(-1.0 + 2.0*r,2) + 4.0*(1.0 - r*r)*pow(cos(s),2);
}

//================================================================================   
double RMC_YF::RoundedFunctiondf1(double s, double r) const
{
	double t1 = sqrt(-4.0*r + 5.0*r*r + 4.0*(1 - r*r)*pow(cos(s),2));	
	return -2.0*(1.0 - r*r)*sin(s) - (4.0*(-1.0 + 2.0*r)*(1 - r*r)*cos(s)*sin(s))/t1;
}

//================================================================================   
double RMC_YF::RoundedFunctiondf2(double s, double r) const
{
	return -8.0*(1 - r*r)*cos(s)*sin(s);
}

//================================================================================   
double RMC_YF::RoundedFunction(double s, double r) const
{
	double f1 = RoundedFunctionf1(s, r);
	double f2 = RoundedFunctionf2(s, r);
	
	return f1/f2;
}

//Guanzhou added for parallel
int RMC_YF::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(6);
    idData.Zero();

    idData(0) = a_which;    
    idData(1) = index_a;    
    idData(2) = k_which;    
    idData(3) = index_k;    
    idData(4) = r_which;
    idData(5) = index_r;
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "DP_YF::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int RMC_YF::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(6);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "DP_YF::recvSelf -- failed to recv ID\n";
	return -1;
    }

    a_which       = idData(0);
    index_a       = idData(1);
    k_which     = idData(2);
    index_k     = idData(3);
    r_which = idData(4);
    index_r = idData(5);
    
    return 0;
}

#endif
