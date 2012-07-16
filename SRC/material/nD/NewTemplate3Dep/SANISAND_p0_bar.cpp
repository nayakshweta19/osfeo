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

#ifndef SANISAND_p0_bar_CPP
#define SANISAND_p0_bar_CPP

#include "SANISAND_p0_bar.h"

stresstensor SANISAND_p0_bar::SANISAND_p0_bar_t;

// Note: V=1000 is used by default
static double V = 1000.0; // this is the coefficient of the exponential function

const double OneOverThree = 0.3333333333;
const double TwoOverThree = 0.6666666667;
const double rt23 = sqrt(TwoOverThree);    
const double rt32 = 1.0/rt23;    

SANISAND_p0_bar::SANISAND_p0_bar(int e0_index_in,      
                                 int Pat_index_in,
                                 int alpha_cc_index_in,
                                 int p_r_index_in,     
                                 int rho_c_index_in,   
                                 int theta_c_index_in,   
                                 int K0_index_in,      
                                 int p0_index_in,      
                                 int alpha_index_in)   
: ScalarEvolution(SE_TAG_SANISAND_p0_bar), 
  e0_index(e0_index_in),
  Pat_index(Pat_index_in),
  alpha_cc_index(alpha_cc_index_in),
  p_r_index(p_r_index_in),
  rho_c_index(rho_c_index_in),
  theta_c_index(theta_c_index_in),
  K0_index(K0_index_in),
  p0_index(p0_index_in),
  alpha_index(alpha_index_in)
{

}

ScalarEvolution* SANISAND_p0_bar::newObj()
{
    ScalarEvolution* nObj = new SANISAND_p0_bar(this->e0_index,      
                                                this->Pat_index,
                                                this->alpha_cc_index,
                                                this->p_r_index,     
                                                this->rho_c_index,   
                                                this->theta_c_index,   
                                                this->K0_index,      
                                                this->p0_index,      
                                                this->alpha_index);  
    return nObj;
}


double SANISAND_p0_bar::H(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
    double e0 = gete0(material_parameter);
    double Pat = getPat(material_parameter);
    double alpha_cc = getalpha_cc(material_parameter);
    double p_r = getp_r(material_parameter);
    double rho_c = getrho_c(material_parameter);
    double theta_c = gettheta_c(material_parameter);
    double K0 = getK0(material_parameter);
    double p0 = getp0(material_parameter);
    stresstensor alpha = getalpha(material_parameter);

    stresstensor s_bar;
    double norm_s = 0.0;
    double r_ef = 0.0;
    double temp_scalar = 0.0;
    double p_b = 0.0;
    double e = e0;
    double delta = 0.0;
    double p0_bar = 0.0;

    double p = Stre.p_hydrostatic();
    stresstensor s = Stre.deviator();

    s_bar = s - (alpha *p);
    norm_s = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
    if (p > 0.0 && norm_s > 0.0)
      r_ef = rt32 * norm_s / p;
    e = e0 + (1.0 + e0) * Stra.Iinvariant1();
    temp_scalar = (alpha("ij")*alpha("ij")).trace();
    p_b = p_r * pow(e,(-1.0/rho_c));   
    if ( p_r > p )
      delta = 1.0 - (p/p_b) * (1.0 + 3*temp_scalar/(alpha_cc*alpha_cc));
    else
      delta = 1.0 - (1.0 + 3*temp_scalar/(alpha_cc*alpha_cc));

    if (delta > 0 )
      p0_bar = (1.0+e) * p0 * exp(-V*r_ef) / (e * (rho_c - pow(p0/Pat,0.3333333333)*(1.0/K0)) * (1.0 - pow(delta,theta_c)) );
    else
      p0_bar = (1.0+e) * p0 * exp(-V*r_ef) / (e * (rho_c - pow(p0/Pat,0.3333333333)*(1.0/K0)) * (1.0 + pow(-delta,theta_c)) );

    return p0_bar;
}
										   



//out- ///////////////////////////////////////////////////////////////////////////////
//out- const tensor& CC_Ev::DH_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
//out-                           const straintensor& Stra, const MaterialParameter& material_parameter)
//out- {
//out- // this is d \bar{p_0} / d sigma_ij
//out- 
//out-     tensor I2("I", 2, def_dim_2);
//out- 
//out-     double M = getM(material_parameter);
//out-     double p0 = getp0(material_parameter);
//out-     double lambda = getlambda(material_parameter);
//out-     double kappa = getkappa(material_parameter);
//out-     double e0 = gete0(material_parameter);
//out-     
//out-     double e = e0 + (1.0 + e0) *Stra.Iinvariant1();
//out- 
//out- //    double scalar1 = (1.0+e)*p0*M*M*(-2.0/3.0)/(lambda-kappa); // Zhao
//out-     double scalar1 = (1.0+e)*p0*M*M*(-2.0/3.0)/(lambda-kappa); // Mahdi and Boris 27April2007
//out- 
//out-     ScalarEvolution::SE_tensorR2 = I2 * scalar1;
//out- // test    ScalarEvolution::SE_tensorR2 = I2 * scalar1 *0.0;
//out- //    ScalarEvolution::SE_tensorR2 = I2 * scalar1 *0.0;
//out-     
//out-     return ScalarEvolution::SE_tensorR2;
//out- }
//out- 
//out- ///////////////////////////////////////////////////////////////////////////////
//out- double CC_Ev::DH_Diso(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
//out-                           const straintensor& Stra, const MaterialParameter& material_parameter)
//out- {
//out- // this is d \bar{p_0} / d p_0
//out- 
//out- 
//out-     double M = getM(material_parameter);
//out-     double p0 = getp0(material_parameter);
//out-     double lambda = getlambda(material_parameter);
//out-     double kappa = getkappa(material_parameter);
//out-     double e0 = gete0(material_parameter);
//out-     
//out-     double e = e0 + (1.0 + e0) *Stra.Iinvariant1();
//out-     double p = Stre.p_hydrostatic();
//out- 
//out- //    double scalar1 = (1.0+e)*p0*M*M*(p0-p)/(lambda-kappa); // Zhao's version
//out-     double scalar1 = (-2.0)*(1.0+e)*p0*M*M*(p0-p)/(lambda-kappa); // Mahdi
//out-     
//out-     return scalar1;
//out- //    return 0.0;
//out- // test      return 0.0;
//out- }




// to get e0
//================================================================================
double SANISAND_p0_bar::gete0(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, e0_index);
}

// to get alpha_cc
//================================================================================
double SANISAND_p0_bar::getalpha_cc(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, alpha_cc_index);
}

// to get Pat
//================================================================================
double SANISAND_p0_bar::getPat(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, Pat_index);
}

// to get p_r
//================================================================================
double SANISAND_p0_bar::getp_r(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, p_r_index);
}

// to get rho_c
//================================================================================
double SANISAND_p0_bar::getrho_c(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, rho_c_index);
}

// to get theta_c
//================================================================================
double SANISAND_p0_bar::gettheta_c(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, theta_c_index);
}

// to get K0
//================================================================================
double SANISAND_p0_bar::getK0(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, K0_index);
}

// to get p0
//================================================================================
double SANISAND_p0_bar::getp0(const MaterialParameter& material_parameter) const
{
    double p0 = 0.0;
    p0 = material_parameter.getInternal_Scalar(p0_index-1);
    return p0;
}

// to get alpha
//================================================================================
const stresstensor& SANISAND_p0_bar::getalpha(const MaterialParameter& material_parameter) const
{
	if ( alpha_index <= material_parameter.getNum_Internal_Tensor() && alpha_index > 0) {
		SANISAND_p0_bar::SANISAND_p0_bar_t = material_parameter.getInternal_Tensor(alpha_index-1);
		return SANISAND_p0_bar::SANISAND_p0_bar_t;
	}
	else {
		opserr << "SANISAND_alpha: Invalid Input (alpha) " << endln;
		exit (1);
	}
}


//================================================================================
double SANISAND_p0_bar::getParameters(const MaterialParameter& material_parameter, int which) const
{
	if ( which <= material_parameter.getNum_Material_Constant() && which > 0)
		return material_parameter.getMaterial_Constant(which-1);
	else {
		opserr << "SANISAND_p0_bar: Invalid Input - #" << which << endln;
		exit (1);
	}
} 


//Guanzhou added for parallel
int SANISAND_p0_bar::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(9);
    idData.Zero();

    idData(0) = e0_index;       
    idData(1) = alpha_cc_index; 
    idData(2) = Pat_index;      
    idData(3) = p_r_index;      
    idData(4) = rho_c_index;    
    idData(5) = theta_c_index;    
    idData(6) = K0_index;       
    idData(7) = p0_index;       
    idData(8) = alpha_index;    
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "SANISAND_p0_bar::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int SANISAND_p0_bar::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(9);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "SANISAND_p0_bar::recvSelf -- failed to recv ID\n";
	return -1;
    }

    e0_index       = idData(0);		     
    alpha_cc_index = idData(1);		     
    Pat_index      = idData(2);		     
    p_r_index      = idData(3);		     
    rho_c_index    = idData(4);		     
    theta_c_index    = idData(5);		     
    K0_index       = idData(6);		     
    p0_index       = idData(7);		     
    alpha_index    = idData(8);		     
    
    return 0;
}

#endif


