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
// UPDATE HISTORY:    Guanzhou Jie updated for parallel, Dec. 2006
//
///////////////////////////////////////////////////////////////////////////////
//

// Ref: Dafalias and Manzari 2004: J. Eng. Mech. 130(6), pp 622-634
// Parameters:
//  1- e0:         initial void ratio at zero strain;
//  2- e_r:       reference void for critical state line, ec = e_r - lambda_c*(pc/Pat)^xi;
//  3- lambda_c:  parameter for critical state line, ec = e_r - lambda_c*(pc/Pat)^xi;
//  4- xi:        parameter for critical state line, ec = e_r - lambda_c*(pc/Pat)^xi;
//  5- Pat:       atmospherics pressure for critical state line, ec = e0 - lambda_c*(pc/Pat)^xi;
//  6- m:         parameter in the yield function;
//  7- M:         critical state stress ration;
//  8- cc:        tension-compression strength ratio;
//  9- A0:        parameter;
// 10- nd         parameter;
// 11- c_z:       parameter
// 12- z_max      parameter
// 13- alpha:     "back-stress" tensor in yield function; (the 1st tensorial internal variable);
// 14- z:         fabric dilation internal tensor (the 2nd tensorial internal variable); 

#ifndef DM04_z_Eij_CPP
#define DM04_z_Eij_CPP

#include "DM04_z_Eij.h"

stresstensor DM04_z_Eij::DM04_z_t;

DM04_z_Eij::DM04_z_Eij(int index_e0_in,
            int index_e_r_in,
            int index_lambda_c_in,
            int index_xi_in,
            int index_Pat_in,
            int index_m_in,
            int index_M_cal_in,
            int index_cc_in,        
            int index_A0_in,
            int index_nd_in,
            int c_z_index_in,
            int z_max_index_in,
            int alpha_index_in,
            int z_index_in)
: TensorEvolution(TE_TAG_DM04_z_Eij), index_e0(index_e0_in),
  index_e_r(index_e_r_in),
  index_lambda_c(index_lambda_c_in),
  index_xi(index_xi_in),
  index_Pat(index_Pat_in),
  index_m(index_m_in),
  index_M_cal(index_M_cal_in),
  index_cc(index_cc_in),        
  index_A0(index_A0_in),
  index_nd(index_nd_in),
  c_z_index(c_z_index_in), 
  z_max_index(z_max_index_in),
  alpha_index(alpha_index_in), 
  z_index(z_index_in)
{

}

TensorEvolution* DM04_z_Eij::newObj()
{
    TensorEvolution* nObj = new DM04_z_Eij(this->index_e0,
                                           this->index_e_r,
                                           this->index_lambda_c,
                                           this->index_xi,
                                           this->index_Pat,
                                           this->index_m,
                                           this->index_M_cal,
                                           this->index_cc,        
                                           this->index_A0,
                                           this->index_nd,
                                           this->c_z_index,
                                           this->z_max_index,
                                           this->alpha_index,
                                           this->z_index);

    return nObj;
}

//================================================================================
const straintensor& DM04_z_Eij::Hij(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                                    const straintensor& Stra, const MaterialParameter& material_parameter)
{
    double c_z = getc_z(material_parameter);
    double z_max = getz_max(material_parameter);        
    stresstensor alpha = getalpha(material_parameter);
    stresstensor z = getz(material_parameter);
    straintensor PF = plastic_flow.PlasticFlowTensor(Stre, Stra, material_parameter);
        
    double p = Stre.p_hydrostatic();
    stresstensor s = Stre.deviator();

    stresstensor n;
    
    stresstensor s_bar = Stre.deviator() - (alpha *p);
    double norm_s = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
    if (p > 0.0 && norm_s > 0.0)
      n = s_bar * (1.0/norm_s);
    
    double D = PF.Iinvariant1();
    if (D < 0.0) 
      D = 0.0;  
   
    TensorEvolution::TensorEvolutionHij = ((n *z_max) +z) *(-c_z*D); 
    
    return TensorEvolution::TensorEvolutionHij;
}


//================================================================================
const tensor& DM04_z_Eij::DHij_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{
    const double oneOver3 = 1.0/3.0;    
    const double rt23 = sqrt(2.0/3.0);

    tensor I2("I", 2, def_dim_2);    
    tensor I4 = I2("ij")*I2("kl");
    tensor I4s = ( I4.transpose0110() + I4.transpose0111() ) *0.5;
    tensor I4dev = I4s - I4 * oneOver3;

    tensor Z4(4, def_dim_4, 0.0);

    double e0 = gete0(material_parameter);
    double e_r = gete_r(material_parameter);
    double lambda_c = getlambda_c(material_parameter);
    double xi = getxi(material_parameter);
    double Pat = getPat(material_parameter);
    //double m = getm(material_parameter);
    double M_cal = getM_cal(material_parameter);
    double cc = getcc(material_parameter);
    double A0 = getA0(material_parameter);
    double nd = getnd(material_parameter); 
    double c_z = getc_z(material_parameter);
    double z_max = getz_max(material_parameter);        
    stresstensor alpha = getalpha(material_parameter);
    stresstensor z = getz(material_parameter);

    double g = 0.0;
    double ec = e_r;
    double stateParameter = 0.0;
    double expnd = 1.0;
    //double ad = 0.0;
    double A_d = 0.0;
    double D = 0.0;
    double D0 = 0.0;
    double epsilon_v = 0.0;
    stresstensor s_bar;
    double e = e0;
    double J3D;
    double cos3theta = 0.0;
    double z_n = 0.0;
    double alpha_n = 0.0;
    double s_n = 0.0;        
        
    double p = Stre.p_hydrostatic();
    stresstensor s = Stre.deviator();

    stresstensor n;
    
    s_bar = Stre.deviator() - (alpha *p);
    double norm_s = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
    if (p > 0.0 && norm_s > 0.0)
        n = s_bar * (1.0/norm_s);

    J3D = n.Jinvariant3();
    cos3theta = -3.0*sqrt(6.0) *J3D;
    
    if (cos3theta > 1.0) 
      cos3theta = 1.0;

    if (cos3theta < -1.0) 
      cos3theta = -1.0;
    
    g = getg(cc, cos3theta);

    if ( (p/Pat) >= 0.0 )
      ec = getec(e_r, lambda_c, xi, Pat, p);
    
    epsilon_v = Stra.Iinvariant1();
    e = e0 + (1.0 + e0) *epsilon_v;
    
    stateParameter = e - ec;
    
    expnd = exp(nd*stateParameter);
    
    alpha_n = (alpha("ij")*n("ij")).trace();
    s_n = (s("ij")*n("ij")).trace();

    // way 1
    //ad = g*M_cal*expnd - m;
    //D0 = rt23 * ad - alpha_n;

    // way 2
    D0 = rt23*g*M_cal*expnd - s_n /p;

    z_n = (z("ij")*n("ij")).trace();
    if (z_n < 0.0) 
      z_n = 0.0;
    A_d = A0 * (1.0 + z_n);

    D = D0 *(-A_d);
    
    if (D <= 0.0) {
       TensorEvolution::TE_tensorR4.Initialize(Z4);
       return TensorEvolution::TE_tensorR4; 
    }   
       
    alpha_n = (alpha("ij")*n("ij")).trace();

    tensor n_n = n("ik")*n("kj");
      n_n.null_indices();

    tensor nt_nt = n("ij")*n("kl");
      nt_nt.null_indices();

    tensor alpha_I = alpha("ij")*I2("kl");
      alpha_I.null_indices();

    tensor n_I = n("ij")*I2("kl");
      n_I.null_indices();

    // dn_ds:
    tensor dn_ds = I4dev - nt_nt + alpha_I*oneOver3 - n_I*(alpha_n*oneOver3);
    dn_ds = dn_ds *(1.0/norm_s);

    // dphi_ds:
    tensor dphi_ds = I2 * (-oneOver3*lambda_c*xi*pow(p, xi-1.0)/pow(Pat, xi));

    // dcos3theta_ds:
    tensor dcos3theta_ds = dn_ds("ijmn")*n_n("ji");
      dcos3theta_ds.null_indices();
    dcos3theta_ds = dcos3theta_ds *(-3.0*sqrt(6.0));

    // dg_ds:
    tensor dg_ds = dcos3theta_ds *(g*g*(1.0-cc)/(2.0*cc));

    // dad_ds:
    tensor dad_ds = (dg_ds + dphi_ds *(g*nd)) *(M_cal*expnd);
    
    // dD_ds:

    // way 1
    //tensor tensor1 = alpha("pq")*dn_ds("pqmn");
    //  tensor1.null_indices();
    
    // way 2
    tensor tensor1 = s("pq")*dn_ds("pqmn");
      tensor1.null_indices();
    tensor1 = (tensor1 + n + I2*(s_n/3.0/p)) *(1.0/p);
    
    tensor dD_ds = (dad_ds *rt23 - tensor1) *(-A_d);    
    if (z_n > 0.0) {
      tensor tensor2 = z("pq")*dn_ds("pqmn");
        tensor2.null_indices();
      dD_ds += tensor2 *(-A0*D0);
    }
                  
    tensor tensor4 = n *z_max + z;
    tensor tensor3 = tensor4("ij")*dD_ds("mn");
      tensor3.null_indices();
    
    TensorEvolution::TE_tensorR4 = (tensor3 + dn_ds*(D*z_max)) *(-c_z); 
    
    return TensorEvolution::TE_tensorR4;
}

//================================================================================
const tensor& DM04_z_Eij::DHij_Dkin(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{  
    const double rt23 = sqrt(2.0/3.0);

    tensor I2("I", 2, def_dim_2);    
    tensor I4 = I2("ij")*I2("kl");
    tensor I4s = ( I4.transpose0110() + I4.transpose0111() ) *0.5;

    tensor Z4(4, def_dim_4, 0.0);

    double e0 = gete0(material_parameter);
    double e_r = gete_r(material_parameter);
    double lambda_c = getlambda_c(material_parameter);
    double xi = getxi(material_parameter);
    double Pat = getPat(material_parameter);
    //double m = getm(material_parameter);
    double M_cal = getM_cal(material_parameter);
    double cc = getcc(material_parameter);
    double A0 = getA0(material_parameter);
    double nd = getnd(material_parameter); 
    double c_z = getc_z(material_parameter);
    double z_max = getz_max(material_parameter);        
    stresstensor alpha = getalpha(material_parameter);
    stresstensor z = getz(material_parameter);

    double g = 0.0;
    double ec = e_r;
    double stateParameter = 0.0;
    double expnd = 1.0;
    //double ad = 0.0;
    double A_d = 0.0;
    double D = 0.0;
    double D0 = 0.0;
    double epsilon_v = 0.0;
    stresstensor s_bar;
    double e = e0;
    double J3D;
    double cos3theta = 0.0;
    double z_n = 0.0;
    double alpha_n = 0.0;
    double s_n = 0.0;        
    
    double p = Stre.p_hydrostatic();
    stresstensor s = Stre.deviator();

    stresstensor n;
    
    s_bar = Stre.deviator() - (alpha *p);
    double norm_s = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
    if (p > 0.0 && norm_s > 0.0)
      n = s_bar * (1.0/norm_s);

    J3D = n.Jinvariant3();
    cos3theta = -3.0*sqrt(6.0) *J3D;
    
    if (cos3theta > 1.0) 
      cos3theta = 1.0;

    if (cos3theta < -1.0) 
      cos3theta = -1.0;
    
    g = getg(cc, cos3theta);

    if ( (p/Pat) >= 0.0 )
      ec = getec(e_r, lambda_c, xi, Pat, p);
    
    epsilon_v = Stra.Iinvariant1();
    e = e0 + (1.0 + e0) *epsilon_v;
    
    stateParameter = e - ec;
    
    expnd = exp(nd*stateParameter);
    
    alpha_n = (alpha("ij")*n("ij")).trace();
    s_n = (s("ij")*n("ij")).trace();

    // way 1
    //ad = g*M_cal*expnd - m;
    //D0 = rt23 * ad - alpha_n;

    // way 2
    D0 = rt23*g*M_cal*expnd - s_n /p;

    z_n = (z("ij")*n("ij")).trace();
    if (z_n < 0.0) 
      z_n = 0.0;
    A_d = A0 * (1.0 + z_n);

    D = D0 *(-A_d);

    if (D <= 0.0) {
       TensorEvolution::TE_tensorR4.Initialize(Z4);
       return TensorEvolution::TE_tensorR4; 
    }  

    // dAd_dz:
    tensor dAd_dz(2, def_dim_2, 0.0);
    if (z_n > 0.0)
      dAd_dz = n *A0;
      
    // dD_dz:
    tensor dD_dz = dAd_dz *(-D0);

    tensor tensor1 = n("ij")*dD_dz("mn");    
    TensorEvolution::TE_tensorR4 = (tensor1*z_max + I4s*D) *(-c_z); 
    
    return TensorEvolution::TE_tensorR4;
}

//================================================================================
const tensor& DM04_z_Eij::DHij_Dkin2(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                          const straintensor& Stra, const MaterialParameter& material_parameter)
{ 
    const double rt23 = sqrt(2.0/3.0);

    tensor I2("I", 2, def_dim_2);    
    tensor I4 = I2("ij")*I2("kl");
    tensor I4s = ( I4.transpose0110() + I4.transpose0111() ) *0.5;

    tensor Z4(4, def_dim_4, 0.0);

    double e0 = gete0(material_parameter);
    double e_r = gete_r(material_parameter);
    double lambda_c = getlambda_c(material_parameter);
    double xi = getxi(material_parameter);
    double Pat = getPat(material_parameter);
    //double m = getm(material_parameter);
    double M_cal = getM_cal(material_parameter);
    double cc = getcc(material_parameter);
    double A0 = getA0(material_parameter);
    double nd = getnd(material_parameter); 
    double c_z = getc_z(material_parameter);
    double z_max = getz_max(material_parameter);        
    stresstensor alpha = getalpha(material_parameter);
    stresstensor z = getz(material_parameter);

    double g = 0.0;
    double ec = e_r;
    double stateParameter = 0.0;
    double expnd = 1.0;
    //double ad = 0.0;
    double A_d = 0.0;
    double D = 0.0;
    double D0 = 0.0;
    double epsilon_v = 0.0;
    stresstensor s_bar;
    double e = e0;
    double J3D;
    double cos3theta = 0.0;
    double z_n = 0.0;
    double alpha_n = 0.0;
    double s_n = 0.0;        
        
    double p = Stre.p_hydrostatic();
    stresstensor s = Stre.deviator();

    stresstensor n;
    
    s_bar = Stre.deviator() - (alpha *p);
    double norm_s = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
    if (p > 0.0 && norm_s > 0.0)
      n = s_bar * (1.0/norm_s);

    J3D = n.Jinvariant3();
    cos3theta = -3.0*sqrt(6.0) *J3D;
    
    if (cos3theta > 1.0) 
      cos3theta = 1.0;

    if (cos3theta < -1.0) 
      cos3theta = -1.0;
    
    g = getg(cc, cos3theta);

    if ( (p/Pat) >= 0.0 )
      ec = getec(e_r, lambda_c, xi, Pat, p);
    
    epsilon_v = Stra.Iinvariant1();
    e = e0 + (1.0 + e0) *epsilon_v;
    
    stateParameter = e - ec;
    
    expnd = exp(nd*stateParameter);
    
    alpha_n = (alpha("ij")*n("ij")).trace();
    s_n = (s("ij")*n("ij")).trace();

    // way 1
    //ad = g*M_cal*expnd - m;
    //D0 = rt23 * ad - alpha_n;

    // way 2
    D0 = rt23*g*M_cal*expnd - s_n /p;

    z_n = (z("ij")*n("ij")).trace();
    if (z_n < 0.0) 
      z_n = 0.0;
    A_d = A0 * (1.0 + z_n);

    D = D0 *(-A_d);
    
    if (D <= 0.0) {
       TensorEvolution::TE_tensorR4.Initialize(Z4);
       return TensorEvolution::TE_tensorR4; 
    }   
       
    alpha_n = (alpha("ij")*n("ij")).trace();

    tensor n_n = n("ik")*n("kj");
      n_n.null_indices();

    tensor nt_nt = n("ij")*n("kl");
      nt_nt.null_indices();

    tensor alpha_I = alpha("ij")*I2("kl");
      alpha_I.null_indices();

    tensor n_I = n("ij")*I2("kl");
      n_I.null_indices();

    // dn_dalpha:
    tensor dn_da = nt_nt - I4s;
    dn_da = dn_da *(p/norm_s);

    // dcos3theta_dalpha:
    tensor dcos3theta_da = dn_da("ijmn")*n_n("ji");
      dcos3theta_da.null_indices();
    dcos3theta_da = dcos3theta_da *(-3.0*sqrt(6.0));

    // dg_da:
    tensor dg_da = dcos3theta_da *(g*g*(1.0-cc)/(2.0*cc));

    // dad_da:
    tensor dad_da = dg_da *(M_cal*expnd);

    // dD_da:
    
    // way 1
    //tensor tensor1 = alpha("pq")*dn_da("pqmn");
    //  tensor1.null_indices();
    //tensor dD_da = (dad_da *rt23 - n - tensor1) *(-A_d);
    
    // way 2
    tensor tensor1 = s("pq")*dn_da("pqmn");
      tensor1.null_indices();
    tensor dD_da = (dad_da *rt23 - tensor1 *(1.0/p)) *(-A_d);
    
    if (z_n > 0.0) {
      tensor tensor2 = z("pq")*dn_da("pqmn");
        tensor2.null_indices();
      dD_da += tensor2 *(-A0*D0);
    }
    
    tensor tensor4 = n *z_max + z;              
    tensor tensor3 = tensor4("ij")*dD_da("mn");
      tensor3.null_indices();
    
    TensorEvolution::TE_tensorR4 = (tensor3 + dn_da*(D*z_max)) *(-c_z); 
    
    return TensorEvolution::TE_tensorR4;
}

// to get e0
//================================================================================
double DM04_z_Eij::gete0(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, index_e0);
}

// to get e_r
//================================================================================
double DM04_z_Eij::gete_r(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, index_e_r);
}

// to get lambda_c
//================================================================================
double DM04_z_Eij::getlambda_c(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, index_lambda_c);
}

// to get xi
//================================================================================
double DM04_z_Eij::getxi(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, index_xi);

}

// to get Pat
//================================================================================
double DM04_z_Eij::getPat(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, index_Pat);
}

// to get m
//================================================================================
double DM04_z_Eij::getm(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, index_m);
}

// to get M
//================================================================================
double DM04_z_Eij::getM_cal(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, index_M_cal);
}

// to get c
//================================================================================
double DM04_z_Eij::getcc(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, index_cc);
}

// to get A0
//================================================================================
double DM04_z_Eij::getA0(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, index_A0);

}

// to get n_d
//================================================================================
double DM04_z_Eij::getnd(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, index_nd);
}
// to get c_z
//================================================================================
double DM04_z_Eij::getc_z(const MaterialParameter& material_parameter) const
{
    if ( c_z_index <= material_parameter.getNum_Material_Constant() && c_z_index > 0)
        return material_parameter.getMaterial_Constant(c_z_index-1);
    else {
        opserr << "DM04_alpha: Invalid Input. " << endln;
        exit (1);
    }
}

// to get c
//================================================================================
double DM04_z_Eij::getz_max(const MaterialParameter& material_parameter) const
{
    if ( z_max_index <= material_parameter.getNum_Material_Constant() && z_max_index > 0)
        return material_parameter.getMaterial_Constant(z_max_index-1);
    else {
        opserr << "DM04_alpha: Invalid Input. " << endln;
        exit (1);
    }
}

// to get alpha
//================================================================================
const stresstensor& DM04_z_Eij::getalpha(const MaterialParameter& material_parameter) const
{
    if ( alpha_index <= material_parameter.getNum_Internal_Tensor() && alpha_index > 0) {
        DM04_z_Eij::DM04_z_t = material_parameter.getInternal_Tensor(alpha_index-1);
        return DM04_z_Eij::DM04_z_t;
    }
    else {
        opserr << "DM04_z: Invalid Input. " << endln;
        exit (1);
    }
}

// to get z
//================================================================================
const stresstensor& DM04_z_Eij::getz(const MaterialParameter& material_parameter) const
{
    if ( z_index <= material_parameter.getNum_Internal_Tensor() && z_index > 0) {
        DM04_z_Eij::DM04_z_t = material_parameter.getInternal_Tensor(z_index-1);
        return DM04_z_Eij::DM04_z_t;
    }
    else {
        opserr << "DM04_z: Invalid Input. " << endln;
        exit (1);
    }
}

//================================================================================
double DM04_z_Eij::getParameters(const MaterialParameter &MaterialParameter_in, int which) const
{
	if ( which <= MaterialParameter_in.getNum_Material_Constant() && which > 0)
		return MaterialParameter_in.getMaterial_Constant(which-1);
	else {
		opserr << "DM04_z_Eij: Invalid Input. " << " and " << which << endln;
		exit (1);
	}
} 


//================================================================================
double DM04_z_Eij::getec(double e_r, double lambda_c, double xi, double Pat, double p_c) const
{
    double ee = e_r; 
    
    if ( p_c/Pat >= 0.0 )
      ee = e_r - lambda_c * pow(p_c/Pat, xi);
    else 
      opserr << "Warning: DM04_z_Eij - 'p_c/Pat' less than zero! " << endln;

    return ee;
}

//================================================================================
double DM04_z_Eij::getg(double c, double cos3theta) const
{
    return 2.0 * c / ( (1.0+c) - (1.0-c)*cos3theta );
}

//Guanzhou added for parallel
int DM04_z_Eij::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(14);
    idData.Zero();

		idData(0) = index_e0;
		idData(1) = index_e_r;
		idData(2) = index_lambda_c;
		idData(3) = index_xi;
		idData(4) = index_Pat;
		idData(5) = index_m;
		idData(6) = index_M_cal;
		idData(7) = index_cc;
		idData(8) = index_A0;
		idData(9) = index_nd;
    idData(10) = c_z_index;  
    idData(11) = z_max_index;
    idData(12) = alpha_index;
    idData(13) =	z_index;    
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "DM04_z_Eij::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int DM04_z_Eij::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(14);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "DM04_z_Eij::recvSelf -- failed to recv ID\n";
	return -1;
    }

	index_e0      =     idData(0); 
	index_e_r		=  idData(1); 
	index_lambda_c= 	 idData(2); 
	index_xi		= 	 idData(3); 
	index_Pat		=  idData(4); 
	index_m		= 	 idData(5); 
	index_M_cal	= 	 idData(6); 
	index_cc		= 	 idData(7); 
	index_A0		= 	 idData(8); 
	index_nd		= 	 idData(9); 
     c_z_index  	 =	 idData(10);
     z_max_index	 =	 idData(11);
     alpha_index	 =	 idData(12);
    	z_index    = 	 idData(13);

    return 0;
}


#endif

