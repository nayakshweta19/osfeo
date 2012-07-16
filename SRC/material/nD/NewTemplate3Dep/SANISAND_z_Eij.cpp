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
//  8- nb:        bounding parameter;
//  9- h0:        bounding parameter;
// 10- ch:        bounding parameter;
// 11- G0:        parameter in the elastic part
// 10- m:         opening of the yield surface;
// 11- alpha:     "back-stress" tensor in yield function; (the 1st tensorial internal variable);

#ifndef SANISAND_z_Eij_CPP
#define SANISAND_z_Eij_CPP

#include "SANISAND_z_Eij.h"

stresstensor SANISAND_z_Eij::SANISAND_z_t;

SANISAND_z_Eij::SANISAND_z_Eij(int e0_index_in,	         
   			       int e_r_index_in,         	                                               
   			       int lambda_index_in,      	                                               
   			       int xi_index_in,	         	                                               
   			       int Pat_index_in,         	                                               
   			       int alpha_cc_index_in,    	                                               
   			       int c_index_in,	          	                                               
   			       int nb_index_in,	         	                                               
   			       int h0_index_in,	         	                                               
   			       int ch_index_in,	         	                                               
   			       int G0_index_in,	         	                                               
   			       int m_index_in,	          	                                               
   			       int c_z_index_in,         	                                               
			       int z_max_index_in,       
			       int alpha_index_in,       
			       int z_index_in)           
: TensorEvolution(TE_TAG_SANISAND_z_Eij), 
  e0_index(e0_index_in),
  e_r_index(e_r_index_in), 
  lambda_index(lambda_index_in),
  xi_index(xi_index_in),
  Pat_index(Pat_index_in),
  alpha_cc_index(alpha_cc_index_in),
  c_index(c_index_in),
  nb_index(nb_index_in),
  h0_index(h0_index_in),   
  ch_index(ch_index_in),
  G0_index(G0_index_in),
  m_index(m_index_in),  
  c_z_index(c_z_index_in),  
  z_max_index(z_max_index_in),  
  alpha_index(alpha_index_in),
  z_index(z_index_in)
{

}

TensorEvolution* SANISAND_z_Eij::newObj()
{
    TensorEvolution* nObj = new SANISAND_z_Eij(this->e0_index,       
                                               this->e_r_index,      
                                               this->lambda_index,   
                                               this->xi_index,       
                                               this->Pat_index,      
                                               this->alpha_cc_index, 
                                               this->c_index,        
                                               this->nb_index,       
                                               this->h0_index,       
                                               this->ch_index,       
                                               this->G0_index,       
                                               this->m_index,        
                                               this->c_z_index,        
                                               this->z_max_index,        
                                               this->alpha_index,
                                               this->z_index);   

    return nObj;
}

//================================================================================
const straintensor& SANISAND_z_Eij::Hij(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
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


//out- //================================================================================
//out- const tensor& SANISAND_z_Eij::DHij_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
//out-                           const straintensor& Stra, const MaterialParameter& material_parameter)
//out- {
//out-     const double oneOver3 = 1.0/3.0;    
//out-     const double rt23 = sqrt(2.0/3.0);
//out- 
//out-     tensor I2("I", 2, def_dim_2);    
//out-     tensor I4 = I2("ij")*I2("kl");
//out-     tensor I4s = ( I4.transpose0110() + I4.transpose0111() ) *0.5;
//out-     tensor I4dev = I4s - I4 * oneOver3;
//out- 
//out-     tensor Z4(4, def_dim_4, 0.0);
//out- 
//out-     double e0 = gete0(material_parameter);
//out-     double e_r = gete_r(material_parameter);
//out-     double lambda_c = getlambda_c(material_parameter);
//out-     double xi = getxi(material_parameter);
//out-     double Pat = getPat(material_parameter);
//out-     //double m = getm(material_parameter);
//out-     double M_cal = getM_cal(material_parameter);
//out-     double cc = getcc(material_parameter);
//out-     double A0 = getA0(material_parameter);
//out-     double nd = getnd(material_parameter); 
//out-     double c_z = getc_z(material_parameter);
//out-     double z_max = getz_max(material_parameter);        
//out-     stresstensor alpha = getalpha(material_parameter);
//out-     stresstensor z = getz(material_parameter);
//out- 
//out-     double g = 0.0;
//out-     double ec = e_r;
//out-     double stateParameter = 0.0;
//out-     double expnd = 1.0;
//out-     //double ad = 0.0;
//out-     double A_d = 0.0;
//out-     double D = 0.0;
//out-     double D0 = 0.0;
//out-     double epsilon_v = 0.0;
//out-     stresstensor s_bar;
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
//out-     stresstensor n;
//out-     
//out-     s_bar = Stre.deviator() - (alpha *p);
//out-     double norm_s = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
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
//out-     g = getg(cc, cos3theta);
//out- 
//out-     if ( (p/Pat) >= 0.0 )
//out-       ec = getec(e_r, lambda_c, xi, Pat, p);
//out-     
//out-     epsilon_v = Stra.Iinvariant1();
//out-     e = e0 + (1.0 + e0) *epsilon_v;
//out-     
//out-     stateParameter = e - ec;
//out-     
//out-     expnd = exp(nd*stateParameter);
//out-     
//out-     alpha_n = (alpha("ij")*n("ij")).trace();
//out-     s_n = (s("ij")*n("ij")).trace();
//out- 
//out-     // way 1
//out-     //ad = g*M_cal*expnd - m;
//out-     //D0 = rt23 * ad - alpha_n;
//out- 
//out-     // way 2
//out-     D0 = rt23*g*M_cal*expnd - s_n /p;
//out- 
//out-     z_n = (z("ij")*n("ij")).trace();
//out-     if (z_n < 0.0) 
//out-       z_n = 0.0;
//out-     A_d = A0 * (1.0 + z_n);
//out- 
//out-     D = D0 *(-A_d);
//out-     
//out-     if (D <= 0.0) {
//out-        TensorEvolution::TE_tensorR4.Initialize(Z4);
//out-        return TensorEvolution::TE_tensorR4; 
//out-     }   
//out-        
//out-     alpha_n = (alpha("ij")*n("ij")).trace();
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
//out-     tensor dn_ds = I4dev - nt_nt + alpha_I*oneOver3 - n_I*(alpha_n*oneOver3);
//out-     dn_ds = dn_ds *(1.0/norm_s);
//out- 
//out-     // dphi_ds:
//out-     tensor dphi_ds = I2 * (-oneOver3*lambda_c*xi*pow(p, xi-1.0)/pow(Pat, xi));
//out- 
//out-     // dcos3theta_ds:
//out-     tensor dcos3theta_ds = dn_ds("ijmn")*n_n("ji");
//out-       dcos3theta_ds.null_indices();
//out-     dcos3theta_ds = dcos3theta_ds *(-3.0*sqrt(6.0));
//out- 
//out-     // dg_ds:
//out-     tensor dg_ds = dcos3theta_ds *(g*g*(1.0-cc)/(2.0*cc));
//out- 
//out-     // dad_ds:
//out-     tensor dad_ds = (dg_ds + dphi_ds *(g*nd)) *(M_cal*expnd);
//out-     
//out-     // dD_ds:
//out- 
//out-     // way 1
//out-     //tensor tensor1 = alpha("pq")*dn_ds("pqmn");
//out-     //  tensor1.null_indices();
//out-     
//out-     // way 2
//out-     tensor tensor1 = s("pq")*dn_ds("pqmn");
//out-       tensor1.null_indices();
//out-     tensor1 = (tensor1 + n + I2*(s_n/3.0/p)) *(1.0/p);
//out-     
//out-     tensor dD_ds = (dad_ds *rt23 - tensor1) *(-A_d);    
//out-     if (z_n > 0.0) {
//out-       tensor tensor2 = z("pq")*dn_ds("pqmn");
//out-         tensor2.null_indices();
//out-       dD_ds += tensor2 *(-A0*D0);
//out-     }
//out-                   
//out-     tensor tensor4 = n *z_max + z;
//out-     tensor tensor3 = tensor4("ij")*dD_ds("mn");
//out-       tensor3.null_indices();
//out-     
//out-     TensorEvolution::TE_tensorR4 = (tensor3 + dn_ds*(D*z_max)) *(-c_z); 
//out-     
//out-     return TensorEvolution::TE_tensorR4;
//out- }
//out- 
//out- //================================================================================
//out- const tensor& SANISAND_z_Eij::DHij_Dkin(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
//out-                           const straintensor& Stra, const MaterialParameter& material_parameter)
//out- {  
//out-     const double rt23 = sqrt(2.0/3.0);
//out- 
//out-     tensor I2("I", 2, def_dim_2);    
//out-     tensor I4 = I2("ij")*I2("kl");
//out-     tensor I4s = ( I4.transpose0110() + I4.transpose0111() ) *0.5;
//out- 
//out-     tensor Z4(4, def_dim_4, 0.0);
//out- 
//out-     double e0 = gete0(material_parameter);
//out-     double e_r = gete_r(material_parameter);
//out-     double lambda_c = getlambda_c(material_parameter);
//out-     double xi = getxi(material_parameter);
//out-     double Pat = getPat(material_parameter);
//out-     //double m = getm(material_parameter);
//out-     double M_cal = getM_cal(material_parameter);
//out-     double cc = getcc(material_parameter);
//out-     double A0 = getA0(material_parameter);
//out-     double nd = getnd(material_parameter); 
//out-     double c_z = getc_z(material_parameter);
//out-     double z_max = getz_max(material_parameter);        
//out-     stresstensor alpha = getalpha(material_parameter);
//out-     stresstensor z = getz(material_parameter);
//out- 
//out-     double g = 0.0;
//out-     double ec = e_r;
//out-     double stateParameter = 0.0;
//out-     double expnd = 1.0;
//out-     //double ad = 0.0;
//out-     double A_d = 0.0;
//out-     double D = 0.0;
//out-     double D0 = 0.0;
//out-     double epsilon_v = 0.0;
//out-     stresstensor s_bar;
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
//out-     stresstensor n;
//out-     
//out-     s_bar = Stre.deviator() - (alpha *p);
//out-     double norm_s = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
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
//out-     g = getg(cc, cos3theta);
//out- 
//out-     if ( (p/Pat) >= 0.0 )
//out-       ec = getec(e_r, lambda_c, xi, Pat, p);
//out-     
//out-     epsilon_v = Stra.Iinvariant1();
//out-     e = e0 + (1.0 + e0) *epsilon_v;
//out-     
//out-     stateParameter = e - ec;
//out-     
//out-     expnd = exp(nd*stateParameter);
//out-     
//out-     alpha_n = (alpha("ij")*n("ij")).trace();
//out-     s_n = (s("ij")*n("ij")).trace();
//out- 
//out-     // way 1
//out-     //ad = g*M_cal*expnd - m;
//out-     //D0 = rt23 * ad - alpha_n;
//out- 
//out-     // way 2
//out-     D0 = rt23*g*M_cal*expnd - s_n /p;
//out- 
//out-     z_n = (z("ij")*n("ij")).trace();
//out-     if (z_n < 0.0) 
//out-       z_n = 0.0;
//out-     A_d = A0 * (1.0 + z_n);
//out- 
//out-     D = D0 *(-A_d);
//out- 
//out-     if (D <= 0.0) {
//out-        TensorEvolution::TE_tensorR4.Initialize(Z4);
//out-        return TensorEvolution::TE_tensorR4; 
//out-     }  
//out- 
//out-     // dAd_dz:
//out-     tensor dAd_dz(2, def_dim_2, 0.0);
//out-     if (z_n > 0.0)
//out-       dAd_dz = n *A0;
//out-       
//out-     // dD_dz:
//out-     tensor dD_dz = dAd_dz *(-D0);
//out- 
//out-     tensor tensor1 = n("ij")*dD_dz("mn");    
//out-     TensorEvolution::TE_tensorR4 = (tensor1*z_max + I4s*D) *(-c_z); 
//out-     
//out-     return TensorEvolution::TE_tensorR4;
//out- }
//out- 
//out- //================================================================================
//out- const tensor& SANISAND_z_Eij::DHij_Dkin2(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
//out-                           const straintensor& Stra, const MaterialParameter& material_parameter)
//out- { 
//out-     const double rt23 = sqrt(2.0/3.0);
//out- 
//out-     tensor I2("I", 2, def_dim_2);    
//out-     tensor I4 = I2("ij")*I2("kl");
//out-     tensor I4s = ( I4.transpose0110() + I4.transpose0111() ) *0.5;
//out- 
//out-     tensor Z4(4, def_dim_4, 0.0);
//out- 
//out-     double e0 = gete0(material_parameter);
//out-     double e_r = gete_r(material_parameter);
//out-     double lambda_c = getlambda_c(material_parameter);
//out-     double xi = getxi(material_parameter);
//out-     double Pat = getPat(material_parameter);
//out-     //double m = getm(material_parameter);
//out-     double M_cal = getM_cal(material_parameter);
//out-     double cc = getcc(material_parameter);
//out-     double A0 = getA0(material_parameter);
//out-     double nd = getnd(material_parameter); 
//out-     double c_z = getc_z(material_parameter);
//out-     double z_max = getz_max(material_parameter);        
//out-     stresstensor alpha = getalpha(material_parameter);
//out-     stresstensor z = getz(material_parameter);
//out- 
//out-     double g = 0.0;
//out-     double ec = e_r;
//out-     double stateParameter = 0.0;
//out-     double expnd = 1.0;
//out-     //double ad = 0.0;
//out-     double A_d = 0.0;
//out-     double D = 0.0;
//out-     double D0 = 0.0;
//out-     double epsilon_v = 0.0;
//out-     stresstensor s_bar;
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
//out-     stresstensor n;
//out-     
//out-     s_bar = Stre.deviator() - (alpha *p);
//out-     double norm_s = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
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
//out-     g = getg(cc, cos3theta);
//out- 
//out-     if ( (p/Pat) >= 0.0 )
//out-       ec = getec(e_r, lambda_c, xi, Pat, p);
//out-     
//out-     epsilon_v = Stra.Iinvariant1();
//out-     e = e0 + (1.0 + e0) *epsilon_v;
//out-     
//out-     stateParameter = e - ec;
//out-     
//out-     expnd = exp(nd*stateParameter);
//out-     
//out-     alpha_n = (alpha("ij")*n("ij")).trace();
//out-     s_n = (s("ij")*n("ij")).trace();
//out- 
//out-     // way 1
//out-     //ad = g*M_cal*expnd - m;
//out-     //D0 = rt23 * ad - alpha_n;
//out- 
//out-     // way 2
//out-     D0 = rt23*g*M_cal*expnd - s_n /p;
//out- 
//out-     z_n = (z("ij")*n("ij")).trace();
//out-     if (z_n < 0.0) 
//out-       z_n = 0.0;
//out-     A_d = A0 * (1.0 + z_n);
//out- 
//out-     D = D0 *(-A_d);
//out-     
//out-     if (D <= 0.0) {
//out-        TensorEvolution::TE_tensorR4.Initialize(Z4);
//out-        return TensorEvolution::TE_tensorR4; 
//out-     }   
//out-        
//out-     alpha_n = (alpha("ij")*n("ij")).trace();
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
//out-     tensor dg_da = dcos3theta_da *(g*g*(1.0-cc)/(2.0*cc));
//out- 
//out-     // dad_da:
//out-     tensor dad_da = dg_da *(M_cal*expnd);
//out- 
//out-     // dD_da:
//out-     
//out-     // way 1
//out-     //tensor tensor1 = alpha("pq")*dn_da("pqmn");
//out-     //  tensor1.null_indices();
//out-     //tensor dD_da = (dad_da *rt23 - n - tensor1) *(-A_d);
//out-     
//out-     // way 2
//out-     tensor tensor1 = s("pq")*dn_da("pqmn");
//out-       tensor1.null_indices();
//out-     tensor dD_da = (dad_da *rt23 - tensor1 *(1.0/p)) *(-A_d);
//out-     
//out-     if (z_n > 0.0) {
//out-       tensor tensor2 = z("pq")*dn_da("pqmn");
//out-         tensor2.null_indices();
//out-       dD_da += tensor2 *(-A0*D0);
//out-     }
//out-     
//out-     tensor tensor4 = n *z_max + z;              
//out-     tensor tensor3 = tensor4("ij")*dD_da("mn");
//out-       tensor3.null_indices();
//out-     
//out-     TensorEvolution::TE_tensorR4 = (tensor3 + dn_da*(D*z_max)) *(-c_z); 
//out-     
//out-     return TensorEvolution::TE_tensorR4;
//out- }

// to get e0
//================================================================================
double SANISAND_z_Eij::gete0(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, e0_index);
}

// to get e_r
//================================================================================
double SANISAND_z_Eij::gete_r(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, e_r_index);
}

// to get lambda
//================================================================================
double SANISAND_z_Eij::getlambda(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, lambda_index);
}

// to get xi
//================================================================================
double SANISAND_z_Eij::getxi(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, xi_index);

}

// to get Pat
//================================================================================
double SANISAND_z_Eij::getPat(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, Pat_index);
}

// to get alpha_cc
//================================================================================
double SANISAND_z_Eij::getalpha_cc(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, alpha_cc_index);
}

// to get c
//================================================================================
double SANISAND_z_Eij::getc(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, c_index);
}

// to get n_d
//================================================================================
double SANISAND_z_Eij::getnb(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, nb_index);
}


// to get h0
//================================================================================
double SANISAND_z_Eij::geth0(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, h0_index);

}


// to get ch
//================================================================================
double SANISAND_z_Eij::getch(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, ch_index);

}

// to get G0
//================================================================================
double SANISAND_z_Eij::getG0(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, G0_index);

}

// to get m
//================================================================================
double SANISAND_z_Eij::getm(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, m_index);
}

// to get c_z
//================================================================================
double SANISAND_z_Eij::getc_z(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, c_z_index);
}

// to get z_max
//================================================================================
double SANISAND_z_Eij::getz_max(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, z_max_index);
}


// to get alpha
//================================================================================
const stresstensor& SANISAND_z_Eij::getalpha(const MaterialParameter& material_parameter) const
{
	if ( alpha_index <= material_parameter.getNum_Internal_Tensor() && alpha_index > 0) {
		SANISAND_z_Eij::SANISAND_z_t = material_parameter.getInternal_Tensor(alpha_index-1);
		return SANISAND_z_Eij::SANISAND_z_t;
	}
	else {
		opserr << "SANISAND_alpha: Invalid Input (alpha) " << endln;
		exit (1);
	}
}

// to get z
//================================================================================
const stresstensor& SANISAND_z_Eij::getz(const MaterialParameter& material_parameter) const
{
	if ( z_index <= material_parameter.getNum_Internal_Tensor() && z_index > 0) {
		SANISAND_z_Eij::SANISAND_z_t = material_parameter.getInternal_Tensor(z_index-1);
		return SANISAND_z_Eij::SANISAND_z_t;
	}
	else {
		opserr << "SANISAND_alpha: Invalid Input (z) " << endln;
		exit (1);
	}
}



//================================================================================
double SANISAND_z_Eij::getParameters(const MaterialParameter& material_parameter, int which) const
{
	if ( which <= material_parameter.getNum_Material_Constant() && which > 0)
		return material_parameter.getMaterial_Constant(which-1);
	else {
		opserr << "SANISAND_alpha: Invalid Input - #" << which << endln;
		exit (1);
	}
} 


//================================================================================
double SANISAND_z_Eij::getec(double e_r, double lambda, double xi, double Pat, double p_c) const
{
    double ee = e_r; 
    
    if ( p_c/Pat >= 0.0 )
      ee = e_r - lambda * pow(p_c/Pat, xi);
    else 
      opserr << "Warning: SANISAND_z_Eij - 'p_c/Pat' less than zero! " << endln;

    return ee;
}

//================================================================================
double SANISAND_z_Eij::getg(double c, double cos3theta) const
{
    return 2.0 * c / ( (1.0+c) - (1.0-c)*cos3theta );
}

//Guanzhou added for parallel
int SANISAND_z_Eij::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(16);
    idData.Zero();

    idData(0) = e0_index;      
    idData(1) = e_r_index;     
    idData(2) = lambda_index;
    idData(3) = xi_index;      
    idData(4) = Pat_index;     
    idData(5) = alpha_cc_index;   
    idData(6) = c_index;      
    idData(7) = nb_index;      
    idData(8) = h0_index;      
    idData(9) = ch_index;      
    idData(10) = G0_index;      
    idData(11) = m_index;       
    idData(12) = c_z_index;       
    idData(13) = z_max_index;       
    idData(14) = alpha_index;   
    idData(15) = z_index;   

    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "SANISAND_z_Eij::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int SANISAND_z_Eij::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(16);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "SANISAND_z_Eij::recvSelf -- failed to recv ID\n";
	return -1;
    }

    e0_index       = idData(0);
    e_r_index      = idData(1);
    lambda_index   = idData(2);
    xi_index       = idData(3);
    Pat_index      = idData(4);
    alpha_cc_index = idData(5);
    c_index        = idData(6);
    nb_index       = idData(7);
    h0_index       = idData(8);
    ch_index       = idData(9);      	                 
    G0_index       = idData(10);  
    m_index        = idData(11);  
    c_z_index      = idData(12);  
    z_max_index    = idData(13);  
    alpha_index    = idData(14);  
    z_index        = idData(15);  

    return 0;
}


#endif

