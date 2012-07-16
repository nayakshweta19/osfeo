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

// Parameters:
//  1- e0:        initial void ratio at zero strain;
//  2- e_r:       reference void for critical state line, ec = e_r - lambda*(pc/Pat)^xi;
//  3- lambda:    parameter for critical state line;
//  4- xi:        parameter for critical state line;
//  5- Pat:       atmospherics pressure for critical state line;
//  6- alpha_cc:  critical state stress ration;
//  7- c:         tension-compression strength ratio;
//  8- A0:        dilatancy parameter;
//  9- nd         dilatancy parameter;
// 10- m:         opening of the yield surface;
// 11- alpha:     "back-stress" tensor in yield function; (the 1st tensorial internal variable);
// 12- z:         fabric dilation internal tensor (the 2nd tensorial internal variable); 
// 13- X:         LCC parameter; 


#ifndef SANISAND_PF_CPP
#define SANISAND_PF_CPP

#include "SANISAND_PF.h"

straintensor SANISAND_PF::SANISANDm;
stresstensor SANISAND_PF::SANISANDtemp;

// Note: V=1000 is used by default
static double V = 1000.0; // this is the coefficient of the exponential function

const double OneOverThree = 0.3333333333;
const double TwoOverThree = 0.6666666667;
const double rt23 = sqrt(TwoOverThree);    
const double rt32 = 1.0/rt23;    

//================================================================================


SANISAND_PF::SANISAND_PF(int e0_which_in,       int index_e0_in,
                         int e_r_which_in,      int index_e_r_in,
                         int lambda_which_in,   int index_lambda_in,
                         int xi_which_in,       int index_xi_in,
                         int Pat_which_in,      int index_Pat_in,
                         int alpha_cc_which_in, int index_alpha_cc_in,
                         int c_which_in,        int index_c_in,        
                         int A0_which_in,       int index_A0_in,
                         int nd_which_in,       int index_nd_in,
                         int m_which_in,        int index_m_in,
                         int alpha_which_in,    int index_alpha_in,
                         int z_which_in,        int index_z_in,
                         int X_which_in,        int index_X_in)
: PlasticFlow(PF_TAG_SANISAND), 
  e0_which(e0_which_in), index_e0(index_e0_in), 
  e_r_which(e_r_which_in), index_e_r(index_e_r_in), 
  lambda_which(lambda_which_in), index_lambda(index_lambda_in),
  xi_which(xi_which_in), index_xi(index_xi_in),
  Pat_which(Pat_which_in), index_Pat(index_Pat_in),
  alpha_cc_which(alpha_cc_which_in), index_alpha_cc(index_alpha_cc_in),
  c_which(c_which_in), index_c(index_c_in),
  A0_which(A0_which_in), index_A0(index_A0_in),
  nd_which(nd_which_in), index_nd(index_nd_in),
  m_which(m_which_in), index_m(index_m_in),
  alpha_which(alpha_which_in), index_alpha(index_alpha_in),
  z_which(z_which_in), index_z(index_z_in),
  X_which(X_which_in), index_X(index_X_in)
{

}

//================================================================================
SANISAND_PF::~SANISAND_PF() 
{  

}

//================================================================================
PlasticFlow* SANISAND_PF::newObj() 
{  
    PlasticFlow  *new_PF = new SANISAND_PF(this->e0_which,       this->index_e0,
                                           this->e_r_which,      this->index_e_r,
                                           this->lambda_which,   this->index_lambda,
                                           this->xi_which,       this->index_xi,
                                           this->Pat_which,      this->index_Pat,
                                           this->alpha_cc_which, this->index_alpha_cc,
                                           this->c_which,        this->index_c,
                                           this->A0_which,       this->index_A0,
                                           this->nd_which,       this->index_nd,
                                           this->m_which,        this->index_m,
                                           this->alpha_which,    this->index_alpha,
                                           this->z_which,        this->index_z,
                                           this->X_which,        this->index_X);
    return new_PF;
}

//================================================================================
const straintensor& SANISAND_PF::PlasticFlowTensor(const stresstensor& Stre, 
                                               const straintensor& Stra, 
                                               const MaterialParameter &MaterialParameter_in) const
{
//    const double OneOverThree = 0.3333333333;
//    const double TwoOverThree = 0.6666666667;
//    const double rt23 = sqrt(TwoOverThree);
//    const double rt32 = 1.0/rt23;

    tensor I2("I", 2, def_dim_2);

    double e0 = gete0(MaterialParameter_in);
    double e_r = gete_r(MaterialParameter_in);
    double lambda = getlambda(MaterialParameter_in);
    double xi = getxi(MaterialParameter_in);
    double Pat = getPat(MaterialParameter_in);
    double alpha_cc = getalpha_cc(MaterialParameter_in);
    double c = getc(MaterialParameter_in);
    double A0 = getA0(MaterialParameter_in);
    double nd = getnd(MaterialParameter_in);    
    double m = getm(MaterialParameter_in);    
    stresstensor alpha = getalpha(MaterialParameter_in);
    stresstensor z = getz(MaterialParameter_in);
    double X = getX(MaterialParameter_in);    

//out       ofstream outAlpha;
//out       outAlpha.open("Alpha.txt",ios::app);
//out  
//out//      BJtensor MahdiAlapha = ptr_plastic_flow->getalphaSANISAND(*ptr_material_parameter); 
//out//       BJtensor MahdiAlapha = FTEP.ptr_plastic_flow.getalpha(matpar); 
//out  
//out
//out
//out       outAlpha << alpha.val(1,1) <<  "   "  
//out                << alpha.val(1,2) <<  "   " 
//out                << alpha.val(1,3) <<  "   "
//out                << alpha.val(2,1) <<  "   "  
//out                << alpha.val(2,2) <<  "   " 
//out                << alpha.val(2,3) <<  "   " 
//out                << alpha.val(3,1) <<  "   "  
//out                << alpha.val(3,2) <<  "   " 
//out                << alpha.val(3,3) << endln;
//out
//out       outAlpha.close();
//out// 


    stresstensor n;
    stresstensor s_bar;
    double norm_s = 0.0;
    double r_ef = 0.0;
    double cos3theta = 0.0;
    double g = 0.0;
    double ec = e_r;
    double e = e0;
    double psi = 0.0;
    double Ad = 0.0;
    double alpha_d_c = 0.0;
    double alpha_d = 0.0;
    double D = 0.0;
    double z_n = 0.0;
    double dQodp1 = 0.0; 
    double dQodp2 = 0.0;
    straintensor dQods1;
    straintensor dQods2;

    double p = Stre.p_hydrostatic();
    stresstensor s = Stre.deviator();
            
    s_bar = s - (alpha *p);
    norm_s = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
    if (p > 0.0 && norm_s > 0.0)
    {
      n = s_bar * (1.0/norm_s);
      r_ef = rt32 * norm_s / p;
      cos3theta = -3.0 * sqrt(6.0) * n.Jinvariant3();
    }   

    if (p <= 0.0)
      cos3theta = 1.0;

    if (cos3theta > 1.0) 
      cos3theta = 1.0;

    if (cos3theta < -1.0) 
      cos3theta = -1.0;
    
    g = getg(c, cos3theta);

    if ( p >= 0.0 )
      ec = getec(e_r, lambda, xi, Pat, p);

    e = e0 + (1.0 + e0) * Stra.Iinvariant1();
    psi = e - ec;    
    alpha_d_c = alpha_cc * exp(nd*psi);
    alpha_d = g*alpha_d_c;

    z_n = (z("ij")*n("ij")).trace();
    if (z_n < 0.0) 
      z_n = 0.0;
    Ad = A0 * (1.0 + z_n);

    //Ad = A0;

    // Method 1
    double temp = rt32*(alpha("ij")*n("ij")).trace();
    D = Ad * (alpha_d - temp);

    //// Method 2, better to use this when "p" is small
    //double temp = (rt32/p)*(s("ij")*n("ij")).trace();
    //D = Ad * (alpha_d + m - temp);

    dQodp1 = D * r_ef;
    dQodp2 = exp(-V*r_ef);
    dQods1 = n * rt32 * r_ef;
    dQods2 = s * (1.5/p) * X * exp(-V*r_ef);

    SANISAND_PF::SANISANDm = (dQods1+dQods2) - I2*(OneOverThree)*(dQodp1+dQodp2);
                           
    return SANISAND_PF::SANISANDm;
}

//out- //================================================================================
//out- const tensor& SANISAND_PF::Dm_Ds(const stresstensor &Stre, 
//out-                           const straintensor &Stra, 
//out-                           const MaterialParameter &MaterialParameter_in) const
//out- {
//out-     const double OneOverThree = 1.0/3.0;
//out-     const double rt23 = sqrt(2.0/3.0);
//out- 
//out-     tensor I2("I", 2, def_dim_2);    
//out-     tensor I4 = I2("ij")*I2("kl");
//out-     tensor I4s = ( I4.transpose0110() + I4.transpose0111() ) *0.5;
//out-     tensor I4dev = I4s - I4 * OneOverThree;
//out- 
//out-     double e0 = gete0(MaterialParameter_in);
//out-     double e_r = gete_r(MaterialParameter_in);
//out-     double lambda = getlambda(MaterialParameter_in);
//out-     double xi = getxi(MaterialParameter_in);
//out-     double Pat = getPat(MaterialParameter_in);
//out-     //double m = getm(MaterialParameter_in);
//out-     double alpha_cc = getalpha_cc(MaterialParameter_in);
//out-     double c = getc(MaterialParameter_in);
//out-     double A0 = getA0(MaterialParameter_in);
//out-     double nd = getnd(MaterialParameter_in);    
//out- 
//out-     stresstensor alpha = getalpha(MaterialParameter_in);
//out-     stresstensor z = getz(MaterialParameter_in);
//out- 
//out-     stresstensor n;
//out-     stresstensor alpha_d;
//out-     stresstensor alpha_d_alpha;
//out-     double g = 0.0;
//out-     double ec = e_r;
//out-     double psi = 0.0;
//out-     double expnd = 1.0;
//out-     //double ad = 0.0;
//out-     double Ad = 0.0;
//out-     double B = 1.0;
//out-     double C = 0.0;
//out-     double D0 = 0.0;
//out-     stresstensor s_bar;
//out-     double norm_s = 0.0;
//out-     double epsilon_v = 0.0;
//out-     double e = e0;
//out-     double J3D;
//out-     double cos3theta = 0.0;
//out-     double z_n = 0.0;
//out-     double alpha_n = 0.0;
//out-     double s_n = 0.0;
//out- 
//out-     double p = Stre.p_hydrostatic();
//out-     stresstensor s = Stre.deviator();
//out-             
//out-     s_bar = s - (alpha *p);
//out-     norm_s = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
//out-     if (p > 0.0 && norm_s > 0.0)
//out-       n = s_bar * (1.0/norm_s);
//out-        
//out-     J3D = n.Jinvariant3();
//out-     cos3theta = -3.0*sqrt(6.0) *J3D;
//out-     
//out-     if (cos3theta > 1.0) 
//out-       cos3theta = 1.0;
//out- 
//out-     if (cos3theta < -1.0) 
//out-       cos3theta = -1.0;
//out-     
//out-     g = getg(c, cos3theta);
//out- 
//out-     if ( (p/Pat) >= 0.0 )
//out-       ec = getec(e_r, lambda, xi, Pat, p);
//out-     
//out-     epsilon_v = Stra.Iinvariant1();
//out-     e = e0 + (1.0 + e0) *epsilon_v;
//out-     
//out-     psi = e - ec;    
//out-     expnd = exp(nd*psi);
//out- 
//out-     alpha_n = (alpha("ij")*n("ij")).trace();
//out-     s_n = (s("ij")*n("ij")).trace();
//out- 
//out-     // way 1
//out-     //ad = g*alpha_cc*expnd - m;
//out-     //D0 = rt23 * ad - alpha_n;
//out- 
//out-     // way 2
//out-     D0 = rt23*g*alpha_cc*expnd - s_n /p;
//out- 
//out-     z_n = (z("ij")*n("ij")).trace();
//out-     if (z_n < 0.0) 
//out-       z_n = 0.0;
//out-     Ad = A0 * (1.0 + z_n);
//out- 
//out-     B = 1.0 + 1.5 *((1.0-c)/c) *g *cos3theta;
//out-     C = 3.0 *sqrt(1.5) *((1.0-c)/c) *g;
//out-     
//out-     tensor n_n = n("ik")*n("kj");
//out-       n_n.null_indices();
//out- 
//out-     tensor nt_nt = n("ij")*n("kl");
//out-       nt_nt.null_indices();
//out- 
//out-     tensor alpha_I = alpha("ij")*I2("kl");
//out-       alpha_I.null_indices();
//out- 
//out-     tensor n_I = n("ij")*I2("kl");
//out-       n_I.null_indices();
//out- 
//out-     // dn_ds:
//out-     tensor dn_ds = I4dev - nt_nt + alpha_I*OneOverThree - n_I*(alpha_n*OneOverThree);
//out-     dn_ds = dn_ds *(1.0/norm_s);
//out- 
//out-     // dphi_ds:
//out-     tensor dphi_ds = I2 * (-OneOverThree*lambda*xi*pow(p, xi-1.0)/pow(Pat, xi));
//out- 
//out-     // dcos3theta_ds:
//out-     tensor dcos3theta_ds = dn_ds("ijmn")*n_n("ji");
//out-       dcos3theta_ds.null_indices();
//out-     dcos3theta_ds = dcos3theta_ds *(-3.0*sqrt(6.0));
//out- 
//out-     // dg_ds:
//out-     tensor dg_ds = dcos3theta_ds *(g*g*(1.0-c)/(2.0*c));
//out- 
//out-     // dB_ds:
//out-     tensor dB_ds = (dg_ds*cos3theta + dcos3theta_ds*g) *(1.5*(1.0-c)/c);
//out- 
//out-     // dC_ds:
//out-     tensor dC_ds = dg_ds *(3.0*sqrt(1.5)*(1.0-c)/c);
//out- 
//out-     // dR_ds:
//out-     tensor tensor1 = n("ij")*dB_ds("mn");
//out-       tensor1.null_indices();
//out-     tensor tensor2 = n_n - I2 *OneOverThree;
//out-     tensor tensor3 = tensor2("ij")*dC_ds("mn");
//out-       tensor3.null_indices();
//out-     tensor tensor4 = n("kj")*dn_ds("ikmn");
//out-       tensor4.null_indices();
//out-     tensor4.transpose1100(); 
//out-     tensor dR_ds = dn_ds *B + tensor1 + tensor4 *(2.0*C) + tensor3;
//out- 
//out-     // dad_ds:
//out-     tensor dad_ds = (dg_ds + dphi_ds *(g*nd)) *(alpha_cc*expnd);
//out- 
//out-     // dD_ds
//out- 
//out-     // way 1
//out-     //tensor tensor5 = alpha("pq")*dn_ds("pqmn");
//out-     //  tensor5.null_indices();
//out-     // way 2
//out-     tensor tensor5 = s("pq")*dn_ds("pqmn");
//out-       tensor5.null_indices();
//out-     tensor5 = (tensor5 + n + I2*(s_n/3.0/p)) *(1.0/p);
//out- 
//out-     tensor dD_ds = (dad_ds *rt23 - tensor5) *Ad;
//out-     if (z_n > 0.0) {
//out-       tensor tensor6 = z("pq")*dn_ds("pqmn");
//out-         tensor6.null_indices();
//out-       dD_ds += tensor6 *(-A0*D0);
//out-     }
//out- 
//out-     // dm_ds:
//out-     tensor tensor7 = I2("ij")*dD_ds("mn");
//out-       tensor7.null_indices();
//out-     
//out-     PlasticFlow::PF_tensorR4 = dR_ds + tensor7 *OneOverThree;   
//out-     
//out-     return PlasticFlow::PF_tensorR4;
//out- }
//out- 
//out- //================================================================================
//out- const tensor& SANISAND_PF::Dm_Dkin(const stresstensor &Stre, 
//out-                           const straintensor &Stra, 
//out-                           const MaterialParameter &MaterialParameter_in) const
//out- {
//out-     const double OneOverThree = 1.0/3.0;
//out-     const double rt23 = sqrt(2.0/3.0);
//out- 
//out-     tensor I2("I", 2, def_dim_2);    
//out-     tensor I4 = I2("ij")*I2("kl");
//out-     tensor I4s = ( I4.transpose0110() + I4.transpose0111() ) *0.5;
//out- 
//out-     double e0 = gete0(MaterialParameter_in);
//out-     double e_r = gete_r(MaterialParameter_in);
//out-     double lambda = getlambda(MaterialParameter_in);
//out-     double xi = getxi(MaterialParameter_in);
//out-     double Pat = getPat(MaterialParameter_in);
//out-     //double m = getm(MaterialParameter_in);
//out-     double alpha_cc = getalpha_cc(MaterialParameter_in);
//out-     double c = getc(MaterialParameter_in);
//out-     double A0 = getA0(MaterialParameter_in);
//out-     double nd = getnd(MaterialParameter_in);    
//out- 
//out-     stresstensor alpha = getalpha(MaterialParameter_in);
//out-     stresstensor z = getz(MaterialParameter_in);
//out- 
//out-     stresstensor n;
//out-     stresstensor alpha_d;
//out-     stresstensor alpha_d_alpha;
//out-     double g = 0.0;
//out-     double ec = e_r;
//out-     double psi = 0.0;
//out-     double expnd = 1.0;
//out-     //double ad = 0.0;
//out-     double Ad = 0.0;
//out-     double B = 1.0;
//out-     double C = 0.0;
//out-     double D0 = 0.0;
//out-     stresstensor s_bar;
//out-     double norm_s = 0.0;
//out-     double epsilon_v = 0.0;
//out-     double e = e0;
//out-     double J3D;
//out-     double cos3theta = 0.0;
//out-     double z_n = 0.0;
//out-     double alpha_n = 0.0;
//out-     double s_n = 0.0;
//out- 
//out-     double p = Stre.p_hydrostatic();
//out-     stresstensor s = Stre.deviator();
//out-             
//out-     s_bar = s - (alpha *p);
//out-     norm_s = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
//out-     if (p > 0.0 && norm_s > 0.0)
//out-         n = s_bar * (1.0/norm_s);
//out-        
//out-     J3D = n.Jinvariant3();
//out-     cos3theta = -3.0*sqrt(6.0) *J3D;
//out-     
//out-     if (cos3theta > 1.0) 
//out-       cos3theta = 1.0;
//out- 
//out-     if (cos3theta < -1.0) 
//out-       cos3theta = -1.0;
//out-     
//out-     g = getg(c, cos3theta);
//out- 
//out-     if ( (p/Pat) >= 0.0 )
//out-       ec = getec(e_r, lambda, xi, Pat, p);
//out-     
//out-     epsilon_v = Stra.Iinvariant1();
//out-     e = e0 + (1.0 + e0) *epsilon_v;
//out-     
//out-     psi = e - ec;    
//out-     expnd = exp(nd*psi);
//out-     
//out-     alpha_n = (alpha("ij")*n("ij")).trace();
//out-     s_n = (s("ij")*n("ij")).trace();
//out- 
//out-     // way 1
//out-     //ad = g*alpha_cc*expnd - m;
//out-     //D0 = rt23 * ad - alpha_n;
//out- 
//out-     // way 2
//out-     D0 = rt23*g*alpha_cc*expnd - s_n /p;
//out- 
//out-     z_n = (z("ij")*n("ij")).trace();
//out-     if (z_n < 0.0) 
//out-       z_n = 0.0;
//out-     Ad = A0 * (1.0 + z_n);
//out- 
//out-     B = 1.0 + 1.5 *((1.0-c)/c) *g *cos3theta;
//out-     C = 3.0 *sqrt(1.5) *((1.0-c)/c) *g;
//out-     
//out-     tensor n_n = n("ik")*n("kj");
//out-       n_n.null_indices();
//out- 
//out-     tensor nt_nt = n("ij")*n("kl");
//out-       nt_nt.null_indices();
//out- 
//out-     tensor alpha_I = alpha("ij")*I2("kl");
//out-       alpha_I.null_indices();
//out- 
//out-     tensor n_I = n("ij")*I2("kl");
//out-       n_I.null_indices();
//out- 
//out-     // dn_dalpha:
//out-     tensor dn_da = nt_nt - I4s;
//out-     dn_da = dn_da *(p/norm_s);
//out- 
//out-     // dcos3theta_dalpha:
//out-     tensor dcos3theta_da = dn_da("ijmn")*n_n("ji");
//out-       dcos3theta_da.null_indices();
//out-     dcos3theta_da = dcos3theta_da *(-3.0*sqrt(6.0));
//out- 
//out-     // dg_da:
//out-     tensor dg_da = dcos3theta_da *(g*g*(1.0-c)/(2.0*c));
//out- 
//out-     // dB_da:
//out-     tensor dB_da = (dg_da*cos3theta + dcos3theta_da*g) *(1.5*(1.0-c)/c);
//out- 
//out-     // dC_ds:
//out-     tensor dC_da = dg_da *(3.0*sqrt(1.5)*(1.0-c)/c);
//out- 
//out-     // dR_da:
//out-     tensor tensor1 = n("ij")*dB_da("mn");
//out-       tensor1.null_indices();
//out-     tensor tensor2 = n_n - I2 *OneOverThree;
//out-     tensor tensor3 = tensor2("ij")*dC_da("mn");
//out-       tensor3.null_indices();
//out-     tensor tensor4 = n("kj")*dn_da("ikmn");
//out-       tensor4.null_indices();
//out-     tensor4.transpose1100(); 
//out-     tensor dR_da = dn_da *B + tensor1 + tensor4 *(2.0*C) + tensor3;
//out- 
//out-     // dad_da:
//out-     tensor dad_da = dg_da *(alpha_cc*expnd);
//out- 
//out-     // dD_da:
//out-     
//out-     // way 1
//out-     //tensor tensor5 = alpha("pq")*dn_da("pqmn");
//out-     //  tensor5.null_indices();
//out-     //tensor dD_da = (dad_da *rt23 - n - tensor5) *(-Ad);
//out-     
//out-     // way 2
//out-     tensor tensor5 = s("pq")*dn_da("pqmn");
//out-       tensor5.null_indices();
//out-     tensor dD_da = (dad_da *rt23 - tensor5 *(1.0/p)) *(-Ad);
//out- 
//out-     if (z_n > 0.0) {
//out-       tensor tensor6 = z("pq")*dn_da("pqmn");
//out-         tensor6.null_indices();
//out-       dD_da += tensor6 *(-A0*D0);
//out-     }
//out- 
//out-     // dm_da:
//out-     tensor tensor7 = I2("ij")*dD_da("mn");
//out-       tensor7.null_indices();
//out-     
//out-     PlasticFlow::PF_tensorR4 = dR_da + tensor7 *OneOverThree;   
//out-     
//out-     return PlasticFlow::PF_tensorR4;
//out- }
//out- 
//out- //================================================================================
//out- const tensor& SANISAND_PF::Dm_Dkin2(const stresstensor &Stre, 
//out-                           const straintensor &Stra, 
//out-                           const MaterialParameter &MaterialParameter_in) const
//out- {
//out-     const double OneOverThree = 1.0/3.0;
//out-     const double rt23 = sqrt(2.0/3.0);
//out-     
//out-     tensor I2("I", 2, def_dim_2);    
//out- 
//out-     double e0 = gete0(MaterialParameter_in);
//out-     double e_r = gete_r(MaterialParameter_in);
//out-     double lambda = getlambda(MaterialParameter_in);
//out-     double xi = getxi(MaterialParameter_in);
//out-     double Pat = getPat(MaterialParameter_in);
//out-     //double m = getm(MaterialParameter_in);
//out-     double alpha_cc = getalpha_cc(MaterialParameter_in);
//out-     double c = getc(MaterialParameter_in);
//out-     double A0 = getA0(MaterialParameter_in);
//out-     double nd = getnd(MaterialParameter_in);    
//out- 
//out-     stresstensor alpha = getalpha(MaterialParameter_in);
//out-     stresstensor z = getz(MaterialParameter_in);
//out- 
//out-     stresstensor n;
//out-     stresstensor alpha_d;
//out-     stresstensor alpha_d_alpha;
//out-     double g = 0.0;
//out-     double ec = e_r;
//out-     double psi = 0.0;
//out-     double expnd = 1.0;
//out-     //double ad = 0.0;
//out-     double D0 = 0.0;
//out-     stresstensor s_bar;
//out-     double norm_s = 0.0;
//out-     double epsilon_v = 0.0;
//out-     double e = e0;
//out-     double J3D;
//out-     double cos3theta = 0.0;
//out-     double z_n = 0.0;
//out-     double alpha_n = 0.0;
//out-     double s_n = 0.0;
//out- 
//out-     double p = Stre.p_hydrostatic();
//out-     stresstensor s = Stre.deviator();
//out-             
//out-     s_bar = s - (alpha *p);
//out-     norm_s = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
//out-     if (p > 0.0 && norm_s > 0.0)
//out-         n = s_bar * (1.0/norm_s);
//out-        
//out-     J3D = n.Jinvariant3();
//out-     cos3theta = -3.0*sqrt(6.0) *J3D;
//out-     
//out-     if (cos3theta > 1.0) 
//out-       cos3theta = 1.0;
//out- 
//out-     if (cos3theta < -1.0) 
//out-       cos3theta = -1.0;
//out-     
//out-     g = getg(c, cos3theta);
//out- 
//out-     if ( (p/Pat) >= 0.0 )
//out-       ec = getec(e_r, lambda, xi, Pat, p);
//out-     
//out-     epsilon_v = Stra.Iinvariant1();
//out-     e = e0 + (1.0 + e0) *epsilon_v;
//out-     
//out-     psi = e - ec;    
//out-     expnd = exp(nd*psi);
//out-     
//out-     alpha_n = (alpha("ij")*n("ij")).trace();
//out-     s_n = (s("ij")*n("ij")).trace();
//out- 
//out-     // way 1
//out-     //ad = g*alpha_cc*expnd - m;
//out-     //D0 = rt23 * ad - alpha_n;
//out- 
//out-     // way 2
//out-     D0 = rt23*g*alpha_cc*expnd - s_n /p;
//out- 
//out-     tensor dD_dz(2, def_dim_2, 0.0);
//out- 
//out-     // dD_dz:
//out-     if (z_n > 0.0)
//out-       dD_dz = n *A0;
//out- 
//out-     // dm_da:
//out-     tensor tensor1 = I2("ij")*dD_dz("mn");
//out-       tensor1.null_indices();
//out-     
//out-     PlasticFlow::PF_tensorR4 = tensor1 *(-D0*OneOverThree);   
//out-     
//out-     return PlasticFlow::PF_tensorR4;
//out- }

// to get e0
//================================================================================
double SANISAND_PF::gete0(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, e0_which, index_e0);
}

// to get e_r
//================================================================================
double SANISAND_PF::gete_r(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, e_r_which, index_e_r);
}

// to get lambda
//================================================================================
double SANISAND_PF::getlambda(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, lambda_which, index_lambda);
}

// to get xi
//================================================================================
double SANISAND_PF::getxi(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, xi_which, index_xi);

}

// to get Pat
//================================================================================
double SANISAND_PF::getPat(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, Pat_which, index_Pat);
}

// to get alpha_cc
//================================================================================
double SANISAND_PF::getalpha_cc(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, alpha_cc_which, index_alpha_cc);
}

// to get c
//================================================================================
double SANISAND_PF::getc(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, c_which, index_c);
}

// to get A0
//================================================================================
double SANISAND_PF::getA0(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, A0_which, index_A0);

}

// to get n_d
//================================================================================
double SANISAND_PF::getnd(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, nd_which, index_nd);
}

// to get m
//================================================================================
double SANISAND_PF::getm(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, m_which, index_m);
}

// to get alpha
//================================================================================
const stresstensor& SANISAND_PF::getalpha(const MaterialParameter &MaterialParameter_in) const
{
	if ( alpha_which == 2 && index_alpha <= MaterialParameter_in.getNum_Internal_Tensor() && index_alpha > 0) {
		SANISAND_PF::SANISANDtemp = MaterialParameter_in.getInternal_Tensor(index_alpha-1);
		return SANISAND_PF::SANISANDtemp;
	}
	else {
		opserr << "SANISAND_PF: Invalid Input. " << endln;
		exit (1);
	}
}

// to get z
//================================================================================
const stresstensor& SANISAND_PF::getz(const MaterialParameter &MaterialParameter_in) const
{
	if ( z_which == 2 && index_z <= MaterialParameter_in.getNum_Internal_Tensor() && index_z > 0) {
		SANISAND_PF::SANISANDtemp = MaterialParameter_in.getInternal_Tensor(index_z-1);
		return SANISAND_PF::SANISANDtemp;
	}
	else {
		opserr << "SANISAND_PF: Invalid Input. " << endln;
		exit (1);
	}
}

// to get X
//================================================================================
double SANISAND_PF::getX(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, X_which, index_X);
}


//================================================================================
double SANISAND_PF::getParameters(const MaterialParameter &MaterialParameter_in, int parIndex, int which) const
{
	if ( parIndex == 0 && which <= MaterialParameter_in.getNum_Material_Constant() && which > 0)
		return MaterialParameter_in.getMaterial_Constant(which-1);
	else if ( parIndex == 1 && which <= MaterialParameter_in.getNum_Internal_Scalar() && which > 0)
		return MaterialParameter_in.getInternal_Scalar(which-1);
	else {
		opserr << "SANISAND_PF: Invalid Input. " << parIndex << " and " << which << endln;
		exit (1);
	}
} 


//================================================================================
double SANISAND_PF::getec(double e_r, double lambda, double xi, double Pat, double p_c) const
{
    double ee = e_r; 
    
    if ( p_c/Pat >= 0.0 )
      ee = e_r - lambda * pow(p_c/Pat, xi);
    else 
      opserr << "Warning: SANISAND_PF - 'p_c/Pat' less than zero! " << endln;

    return ee;
}

//================================================================================
double SANISAND_PF::getg(double c, double cos3theta) const
{
    return 2.0 * c / ( (1.0+c) - (1.0-c)*cos3theta );
}

//Guanzhou added for parallel
int SANISAND_PF::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(26);
    idData.Zero();

    idData(0) = e0_which;      	  idData(13) = index_e0;      
    idData(1) = e_r_which;     	  idData(14) = index_e_r;     
    idData(2) = lambda_which;	  idData(15) = index_lambda;
    idData(3) = xi_which;      	  idData(16) = index_xi;      
    idData(4) = Pat_which;     	  idData(17) = index_Pat;     
    idData(5) = alpha_cc_which;   idData(18) = index_alpha_cc;   
    idData(6) = c_which;      	  idData(19) = index_c;      
    idData(7) = A0_which;      	  idData(20) = index_A0;      
    idData(8) = nd_which;      	  idData(21) = index_nd;      
    idData(9) = m_which;      	  idData(22) = index_m;      
    idData(10) = alpha_which;     idData(23) = index_alpha;   
    idData(11) = z_which;         idData(24) = index_z;       
    idData(12) = X_which;         idData(25) = index_X;       
    
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "SANISAND_PF::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int SANISAND_PF::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(26);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "SANISAND_PF::recvSelf -- failed to recv ID\n";
	return -1;
    }

    e0_which       = idData(0) ;   index_e0       = idData(13);
    e_r_which      = idData(1) ;   index_e_r      = idData(14);
    lambda_which   = idData(2) ;   index_lambda   = idData(15);
    xi_which       = idData(3) ;   index_xi       = idData(16);
    Pat_which      = idData(4) ;   index_Pat      = idData(17);
    alpha_cc_which = idData(5) ;   index_alpha_cc = idData(18);
    c_which        = idData(6) ;   index_c        = idData(19);
    A0_which       = idData(7) ;   index_A0       = idData(20);
    nd_which       = idData(8) ;   index_nd       = idData(21);
    m_which        = idData(9) ;   index_m        = idData(22);
    alpha_which    = idData(10);   index_alpha    = idData(23);
    z_which        = idData(11);   index_z        = idData(24);
    X_which        = idData(12);   index_X        = idData(25);

    return 0;
}

#endif




