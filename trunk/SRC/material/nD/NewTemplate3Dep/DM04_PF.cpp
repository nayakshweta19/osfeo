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
// 11- alpha:     "back-stress" tensor in yield function; (the 1st tensorial internal variable);
// 12- z:         fabric dilation internal tensor (the 2nd tensorial internal variable); 

#ifndef DM04_PF_CPP
#define DM04_PF_CPP

#include "DM04_PF.h"

straintensor DM04_PF::DM04m;
stresstensor DM04_PF::DM04temp;

//================================================================================
DM04_PF::DM04_PF(int e0_which_in, int index_e0_in,
                 int e_r_which_in, int index_e_r_in, 
                 int lambda_c_which_in, int index_lambda_c_in,
                 int xi_which_in, int index_xi_in,
                 int Pat_which_in, int index_Pat_in,
                 int m_which_in, int index_m_in,
                 int M_cal_which_in, int index_M_cal_in,
                 int cc_which_in, int index_cc_in,
                 int A0_which_in, int index_A0_in,
                 int nd_which_in, int index_nd_in,
                 int alpha_which_in, int index_alpha_in,
                 int z_which_in, int index_z_in)
: PlasticFlow(PF_TAG_DM04), e0_which(e0_which_in), index_e0(index_e0_in), 
  e_r_which(e_r_which_in), index_e_r(index_e_r_in), 
  lambda_c_which(lambda_c_which_in), index_lambda_c(index_lambda_c_in),
  xi_which(xi_which_in), index_xi(index_xi_in),
  Pat_which(Pat_which_in), index_Pat(index_Pat_in),
  m_which(m_which_in), index_m(index_m_in),
  M_cal_which(M_cal_which_in), index_M_cal(index_M_cal_in),
  cc_which(cc_which_in), index_cc(index_cc_in),
  A0_which(A0_which_in), index_A0(index_A0_in),
  nd_which(nd_which_in), index_nd(index_nd_in),
  alpha_which(alpha_which_in), index_alpha(index_alpha_in),
  z_which(z_which_in), index_z(index_z_in)
{

}

//================================================================================
DM04_PF::~DM04_PF() 
{  

}

//================================================================================
PlasticFlow* DM04_PF::newObj() 
{  
     PlasticFlow  *new_PF = new DM04_PF(this->e0_which, this->index_e0,
                                        this->e_r_which, this->index_e_r,
                                        this->lambda_c_which, this->index_lambda_c,
                                        this->xi_which, this->index_xi,
                                        this->Pat_which, this->index_Pat,
                                        this->m_which, this->index_m,
                                        this->M_cal_which, this->index_M_cal,
                                        this->cc_which, this->index_cc,
                                        this->A0_which, this->index_A0,
                                        this->nd_which, this->index_nd,
                                        this->alpha_which, this->index_alpha,
                                        this->z_which, this->index_z);
     
     return new_PF;
}

//================================================================================
const straintensor& DM04_PF::PlasticFlowTensor(const stresstensor& Stre, 
                                               const straintensor& Stra, 
                                               const MaterialParameter &MaterialParameter_in) const
{
    const double oneOver3 = 1.0/3.0;
    const double rt23 = sqrt(2.0/3.0);    
    tensor I2("I", 2, def_dim_2);

    double e0 = gete0(MaterialParameter_in);
    double e_r = gete_r(MaterialParameter_in);
    double lambda_c = getlambda_c(MaterialParameter_in);
    double xi = getxi(MaterialParameter_in);
    double Pat = getPat(MaterialParameter_in);
    //double m = getm(MaterialParameter_in);
    double M_cal = getM_cal(MaterialParameter_in);
    double cc = getcc(MaterialParameter_in);
    double A0 = getA0(MaterialParameter_in);
    double nd = getnd(MaterialParameter_in);    

    stresstensor alpha = getalpha(MaterialParameter_in);
    stresstensor z = getz(MaterialParameter_in);

    stresstensor n;
    stresstensor alpha_d;
    stresstensor alpha_d_alpha;
    double g = 0.0;
    double ec = e_r;
    double stateParameter = 0.0;
    double expnd = 1.0;
    //double ad = 0.0;
    double A_d = 0.0;
    double B = 1.0;
    double C = 0.0;
    double D = 0.0;
    double D0 = 0.0;
    stresstensor s_bar;
    double norm_s = 0.0;
    double epsilon_v = 0.0;
    double e = e0;
    double J3D;
    double cos3theta = 0.0;
    double z_n = 0.0;
    double alpha_n = 0.0;
    double s_n = 0.0;

    double p = Stre.p_hydrostatic();
    stresstensor s = Stre.deviator();
            
    s_bar = s - (alpha *p);
    norm_s = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
    if (p > 0.0 && norm_s > 0.0)
      n = s_bar * (1.0/norm_s);
       
    J3D = n.Jinvariant3();
    cos3theta = -3.0*sqrt(6.0) *J3D;

    if (p <= 0.0)
      cos3theta = 1.0;
    
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
    //D0 = rt23 *ad - alpha_n;

    // way 2, better use this when "p" is small
    D0 = rt23*g*M_cal*expnd - s_n /p;

    z_n = (z("ij")*n("ij")).trace();
    if (z_n < 0.0) 
      z_n = 0.0;
    A_d = A0 * (1.0 + z_n);

    D = D0 *(-A_d);

    B = 1.0 + 1.5 *((1.0-cc)/cc) *g *cos3theta;
    C = 3.0 *sqrt(1.5) *((1.0-cc)/cc) *g;
    
    stresstensor n_n = n("ik")*n("kj");
      n_n.null_indices();

    // note different 'positive-negative' since we assume extension (dilation) positive 
    // which is different from the Ref.
    DM04_PF::DM04m = n *B + n_n *C + I2 *((D-C)*oneOver3);
                           
    return DM04_PF::DM04m;
}

//================================================================================
const tensor& DM04_PF::Dm_Ds(const stresstensor &Stre, 
                          const straintensor &Stra, 
                          const MaterialParameter &MaterialParameter_in) const
{
    const double oneOver3 = 1.0/3.0;
    const double rt23 = sqrt(2.0/3.0);

    tensor I2("I", 2, def_dim_2);    
    tensor I4 = I2("ij")*I2("kl");
    tensor I4s = ( I4.transpose0110() + I4.transpose0111() ) *0.5;
    tensor I4dev = I4s - I4 * oneOver3;

    double e0 = gete0(MaterialParameter_in);
    double e_r = gete_r(MaterialParameter_in);
    double lambda_c = getlambda_c(MaterialParameter_in);
    double xi = getxi(MaterialParameter_in);
    double Pat = getPat(MaterialParameter_in);
    //double m = getm(MaterialParameter_in);
    double M_cal = getM_cal(MaterialParameter_in);
    double cc = getcc(MaterialParameter_in);
    double A0 = getA0(MaterialParameter_in);
    double nd = getnd(MaterialParameter_in);    

    stresstensor alpha = getalpha(MaterialParameter_in);
    stresstensor z = getz(MaterialParameter_in);

    stresstensor n;
    stresstensor alpha_d;
    stresstensor alpha_d_alpha;
    double g = 0.0;
    double ec = e_r;
    double stateParameter = 0.0;
    double expnd = 1.0;
    //double ad = 0.0;
    double A_d = 0.0;
    double B = 1.0;
    double C = 0.0;
    double D0 = 0.0;
    stresstensor s_bar;
    double norm_s = 0.0;
    double epsilon_v = 0.0;
    double e = e0;
    double J3D;
    double cos3theta = 0.0;
    double z_n = 0.0;
    double alpha_n = 0.0;
    double s_n = 0.0;

    double p = Stre.p_hydrostatic();
    stresstensor s = Stre.deviator();
            
    s_bar = s - (alpha *p);
    norm_s = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
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

    B = 1.0 + 1.5 *((1.0-cc)/cc) *g *cos3theta;
    C = 3.0 *sqrt(1.5) *((1.0-cc)/cc) *g;
    
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

    // dB_ds:
    tensor dB_ds = (dg_ds*cos3theta + dcos3theta_ds*g) *(1.5*(1.0-cc)/cc);

    // dC_ds:
    tensor dC_ds = dg_ds *(3.0*sqrt(1.5)*(1.0-cc)/cc);

    // dR_ds:
    tensor tensor1 = n("ij")*dB_ds("mn");
      tensor1.null_indices();
    tensor tensor2 = n_n - I2 *oneOver3;
    tensor tensor3 = tensor2("ij")*dC_ds("mn");
      tensor3.null_indices();
    tensor tensor4 = n("kj")*dn_ds("ikmn");
      tensor4.null_indices();
    tensor4.transpose1100(); 
    tensor dR_ds = dn_ds *B + tensor1 + tensor4 *(2.0*C) + tensor3;

    // dad_ds:
    tensor dad_ds = (dg_ds + dphi_ds *(g*nd)) *(M_cal*expnd);

    // dD_ds

    // way 1
    //tensor tensor5 = alpha("pq")*dn_ds("pqmn");
    //  tensor5.null_indices();
    // way 2
    tensor tensor5 = s("pq")*dn_ds("pqmn");
      tensor5.null_indices();
    tensor5 = (tensor5 + n + I2*(s_n/3.0/p)) *(1.0/p);

    tensor dD_ds = (dad_ds *rt23 - tensor5) *A_d;
    if (z_n > 0.0) {
      tensor tensor6 = z("pq")*dn_ds("pqmn");
        tensor6.null_indices();
      dD_ds += tensor6 *(-A0*D0);
    }

    // dm_ds:
    tensor tensor7 = I2("ij")*dD_ds("mn");
      tensor7.null_indices();
    
    PlasticFlow::PF_tensorR4 = dR_ds + tensor7 *oneOver3;   
    
    return PlasticFlow::PF_tensorR4;
}

//================================================================================
const tensor& DM04_PF::Dm_Dkin(const stresstensor &Stre, 
                          const straintensor &Stra, 
                          const MaterialParameter &MaterialParameter_in) const
{
    const double oneOver3 = 1.0/3.0;
    const double rt23 = sqrt(2.0/3.0);

    tensor I2("I", 2, def_dim_2);    
    tensor I4 = I2("ij")*I2("kl");
    tensor I4s = ( I4.transpose0110() + I4.transpose0111() ) *0.5;

    double e0 = gete0(MaterialParameter_in);
    double e_r = gete_r(MaterialParameter_in);
    double lambda_c = getlambda_c(MaterialParameter_in);
    double xi = getxi(MaterialParameter_in);
    double Pat = getPat(MaterialParameter_in);
    //double m = getm(MaterialParameter_in);
    double M_cal = getM_cal(MaterialParameter_in);
    double cc = getcc(MaterialParameter_in);
    double A0 = getA0(MaterialParameter_in);
    double nd = getnd(MaterialParameter_in);    

    stresstensor alpha = getalpha(MaterialParameter_in);
    stresstensor z = getz(MaterialParameter_in);

    stresstensor n;
    stresstensor alpha_d;
    stresstensor alpha_d_alpha;
    double g = 0.0;
    double ec = e_r;
    double stateParameter = 0.0;
    double expnd = 1.0;
    //double ad = 0.0;
    double A_d = 0.0;
    double B = 1.0;
    double C = 0.0;
    double D0 = 0.0;
    stresstensor s_bar;
    double norm_s = 0.0;
    double epsilon_v = 0.0;
    double e = e0;
    double J3D;
    double cos3theta = 0.0;
    double z_n = 0.0;
    double alpha_n = 0.0;
    double s_n = 0.0;

    double p = Stre.p_hydrostatic();
    stresstensor s = Stre.deviator();
            
    s_bar = s - (alpha *p);
    norm_s = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
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

    B = 1.0 + 1.5 *((1.0-cc)/cc) *g *cos3theta;
    C = 3.0 *sqrt(1.5) *((1.0-cc)/cc) *g;
    
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

    // dB_da:
    tensor dB_da = (dg_da*cos3theta + dcos3theta_da*g) *(1.5*(1.0-cc)/cc);

    // dC_ds:
    tensor dC_da = dg_da *(3.0*sqrt(1.5)*(1.0-cc)/cc);

    // dR_da:
    tensor tensor1 = n("ij")*dB_da("mn");
      tensor1.null_indices();
    tensor tensor2 = n_n - I2 *oneOver3;
    tensor tensor3 = tensor2("ij")*dC_da("mn");
      tensor3.null_indices();
    tensor tensor4 = n("kj")*dn_da("ikmn");
      tensor4.null_indices();
    tensor4.transpose1100(); 
    tensor dR_da = dn_da *B + tensor1 + tensor4 *(2.0*C) + tensor3;

    // dad_da:
    tensor dad_da = dg_da *(M_cal*expnd);

    // dD_da:
    
    // way 1
    //tensor tensor5 = alpha("pq")*dn_da("pqmn");
    //  tensor5.null_indices();
    //tensor dD_da = (dad_da *rt23 - n - tensor5) *(-A_d);
    
    // way 2
    tensor tensor5 = s("pq")*dn_da("pqmn");
      tensor5.null_indices();
    tensor dD_da = (dad_da *rt23 - tensor5 *(1.0/p)) *(-A_d);

    if (z_n > 0.0) {
      tensor tensor6 = z("pq")*dn_da("pqmn");
        tensor6.null_indices();
      dD_da += tensor6 *(-A0*D0);
    }

    // dm_da:
    tensor tensor7 = I2("ij")*dD_da("mn");
      tensor7.null_indices();
    
    PlasticFlow::PF_tensorR4 = dR_da + tensor7 *oneOver3;   
    
    return PlasticFlow::PF_tensorR4;
}

//================================================================================
const tensor& DM04_PF::Dm_Dkin2(const stresstensor &Stre, 
                          const straintensor &Stra, 
                          const MaterialParameter &MaterialParameter_in) const
{
    const double oneOver3 = 1.0/3.0;
    const double rt23 = sqrt(2.0/3.0);
    
    tensor I2("I", 2, def_dim_2);    

    double e0 = gete0(MaterialParameter_in);
    double e_r = gete_r(MaterialParameter_in);
    double lambda_c = getlambda_c(MaterialParameter_in);
    double xi = getxi(MaterialParameter_in);
    double Pat = getPat(MaterialParameter_in);
    //double m = getm(MaterialParameter_in);
    double M_cal = getM_cal(MaterialParameter_in);
    double cc = getcc(MaterialParameter_in);
    double A0 = getA0(MaterialParameter_in);
    double nd = getnd(MaterialParameter_in);    

    stresstensor alpha = getalpha(MaterialParameter_in);
    stresstensor z = getz(MaterialParameter_in);

    stresstensor n;
    stresstensor alpha_d;
    stresstensor alpha_d_alpha;
    double g = 0.0;
    double ec = e_r;
    double stateParameter = 0.0;
    double expnd = 1.0;
    //double ad = 0.0;
    double D0 = 0.0;
    stresstensor s_bar;
    double norm_s = 0.0;
    double epsilon_v = 0.0;
    double e = e0;
    double J3D;
    double cos3theta = 0.0;
    double z_n = 0.0;
    double alpha_n = 0.0;
    double s_n = 0.0;

    double p = Stre.p_hydrostatic();
    stresstensor s = Stre.deviator();
            
    s_bar = s - (alpha *p);
    norm_s = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
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

    tensor dD_dz(2, def_dim_2, 0.0);

    // dD_dz:
    if (z_n > 0.0)
      dD_dz = n *A0;

    // dm_da:
    tensor tensor1 = I2("ij")*dD_dz("mn");
      tensor1.null_indices();
    
    PlasticFlow::PF_tensorR4 = tensor1 *(-D0*oneOver3);   
    
    return PlasticFlow::PF_tensorR4;
}

// to get e0
//================================================================================
double DM04_PF::gete0(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, e0_which, index_e0);
}

// to get e_r
//================================================================================
double DM04_PF::gete_r(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, e_r_which, index_e_r);
}

// to get lambda_c
//================================================================================
double DM04_PF::getlambda_c(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, lambda_c_which, index_lambda_c);
}

// to get xi
//================================================================================
double DM04_PF::getxi(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, xi_which, index_xi);

}

// to get Pat
//================================================================================
double DM04_PF::getPat(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, Pat_which, index_Pat);
}

// to get m
//================================================================================
double DM04_PF::getm(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, m_which, index_m);
}

// to get M
//================================================================================
double DM04_PF::getM_cal(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, M_cal_which, index_M_cal);
}

// to get c
//================================================================================
double DM04_PF::getcc(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, cc_which, index_cc);
}

// to get A0
//================================================================================
double DM04_PF::getA0(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, A0_which, index_A0);

}

// to get n_d
//================================================================================
double DM04_PF::getnd(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, nd_which, index_nd);
}

// to get alpha
//================================================================================
const stresstensor& DM04_PF::getalpha(const MaterialParameter &MaterialParameter_in) const
{
	if ( alpha_which == 2 && index_alpha <= MaterialParameter_in.getNum_Internal_Tensor() && index_alpha > 0) {
		DM04_PF::DM04temp = MaterialParameter_in.getInternal_Tensor(index_alpha-1);
		return DM04_PF::DM04temp;
	}
	else {
		opserr << "DM04_PF: Invalid Input. " << endln;
		exit (1);
	}
}

// to get z
//================================================================================
const stresstensor& DM04_PF::getz(const MaterialParameter &MaterialParameter_in) const
{
	if ( z_which == 2 && index_z <= MaterialParameter_in.getNum_Internal_Tensor() && index_z > 0) {
		DM04_PF::DM04temp = MaterialParameter_in.getInternal_Tensor(index_z-1);
		return DM04_PF::DM04temp;
	}
	else {
		opserr << "DM04_PF: Invalid Input. " << endln;
		exit (1);
	}
}


//================================================================================
double DM04_PF::getParameters(const MaterialParameter &MaterialParameter_in, int parIndex, int which) const
{
	if ( parIndex == 0 && which <= MaterialParameter_in.getNum_Material_Constant() && which > 0)
		return MaterialParameter_in.getMaterial_Constant(which-1);
	else if ( parIndex == 1 && which <= MaterialParameter_in.getNum_Internal_Scalar() && which > 0)
		return MaterialParameter_in.getInternal_Scalar(which-1);
	else {
		opserr << "DM04_PF: Invalid Input. " << parIndex << " and " << which << endln;
		exit (1);
	}
} 


//================================================================================
double DM04_PF::getec(double e_r, double lambda_c, double xi, double Pat, double p_c) const
{
    double ee = e_r; 
    
    if ( p_c/Pat >= 0.0 )
      ee = e_r - lambda_c * pow(p_c/Pat, xi);
    else 
      opserr << "Warning: DM04_PF - 'p_c/Pat' less than zero! " << endln;

    return ee;
}

//================================================================================
double DM04_PF::getg(double c, double cos3theta) const
{
    return 2.0 * c / ( (1.0+c) - (1.0-c)*cos3theta );
}

//Guanzhou added for parallel
int DM04_PF::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(24);
    idData.Zero();

    idData(0) = e0_which;      	  idData(12) = index_e0;      
    idData(1) = e_r_which;     	  idData(13) = index_e_r;     
    idData(2) = lambda_c_which;	  idData(14) = index_lambda_c;
    idData(3) = xi_which;      	  idData(15) = index_xi;      
    idData(4) = Pat_which;     	  idData(16) = index_Pat;     
    idData(5) = m_which;       	  idData(17) = index_m;       
    idData(6) = M_cal_which;   	  idData(18) = index_M_cal;   
    idData(7) = cc_which;      	  idData(19) = index_cc;      
    idData(8) = A0_which;      	  idData(20) = index_A0;      
    idData(9) = nd_which;      	  idData(21) = index_nd;      
    		               	  	                      
    idData(10) = alpha_which;     idData(22) = index_alpha;   
    idData(11) = z_which;         idData(23) = index_z;       
    
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "DM04_PF::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int DM04_PF::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(24);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "DM04_PF::recvSelf -- failed to recv ID\n";
	return -1;
    }

    e0_which       = idData(0);    index_e0      = idData(12);
    e_r_which      = idData(1);    index_e_r     = idData(13);
    lambda_c_which = idData(2);    index_lambda_c= idData(14);
    xi_which       = idData(3);    index_xi      = idData(15);
    Pat_which      = idData(4);    index_Pat     = idData(16);
    m_which        = idData(5);    index_m       = idData(17);
    M_cal_which    = idData(6);    index_M_cal   = idData(18);
    cc_which       = idData(7);    index_cc      = idData(19);
    A0_which       = idData(8);    index_A0      = idData(20);
    nd_which       = idData(9);    index_nd      = idData(21);
                      	                            	                 
    alpha_which   = idData(10);   index_alpha   = idData(22);
    z_which       = idData(11);   index_z       = idData(23);

    return 0;
}

#endif




