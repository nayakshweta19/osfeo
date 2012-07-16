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

#ifndef SANISAND_alpha_Eij_CPP
#define SANISAND_alpha_Eij_CPP

#include "SANISAND_alpha_Eij.h"

stresstensor SANISAND_alpha_Eij::SANISAND_alpha_t;

const double OneOverThree = 0.3333333333;
const double TwoOverThree = 0.6666666667;
const double rt23 = sqrt(TwoOverThree);    
const double rt32 = 1.0/rt23;    


SANISAND_alpha_Eij::SANISAND_alpha_Eij(int e0_index_in,
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
                                       int alpha_index_in)
: TensorEvolution(TE_TAG_SANISAND_alpha_Eij), 
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
  alpha_index(alpha_index_in)
{
//   stresstensor zT;
//   alpha_in.Initialize(zT);
}

TensorEvolution* SANISAND_alpha_Eij::newObj()
{
    TensorEvolution* nObj = new SANISAND_alpha_Eij(this->e0_index,
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
                                                   this->alpha_index);
    return nObj;
}

const straintensor& SANISAND_alpha_Eij::Hij(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                                            const straintensor& Stra, const MaterialParameter& material_parameter)
{
//    const double rt23 = sqrt(2.0/3.0);    

//    stresstensor a_a_in;
//    double a_in = 0.0;

    double e0 = gete0(material_parameter);
    double e_r = gete_r(material_parameter);
    double lambda = getlambda(material_parameter);
    double xi = getxi(material_parameter);
    double Pat = getPat(material_parameter);
    double alpha_cc = getalpha_cc(material_parameter);
    double c = getc(material_parameter);
    double nb = getnb(material_parameter);
    double h0 = geth0(material_parameter);
    double ch = getch(material_parameter);
    double G0 = getG0(material_parameter);        
    double m = getm(material_parameter);
    stresstensor alpha = getalpha(material_parameter);

    stresstensor n;
    stresstensor s_bar;
    double norm_s = 0.0;
    double r_ef = 0.0;
    double cos3theta = 0.0;
    double g = 0.0;
    double ec = e_r;
    double e = e0;
    double psi = 0.0;
    double alpha_b_c = 0.0;
    stresstensor alpha_b_tensor;
    stresstensor b_ref;
    stresstensor temp_tensor;
    double lower = 0.0;
    double h = G0*h0;

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
    alpha_b_c = alpha_cc * exp(-nb*psi);
    alpha_b_tensor = n * (rt23 * g * alpha_b_c);
    b_ref = n * rt23 * alpha_b_c * (1.0+c);
//    b_ref = n * rt23 * alpha_cc * (1.0+c);

    // Method 1
    temp_tensor = b_ref - (alpha_b_tensor - alpha);

    //// Method 2, better to use this when "p" is small 
    //temp_tensor = b_ref - (alpha_b_tensor - s*(1.0/p));

    lower = rt32*(temp_tensor("ij")*n("ij")).trace();
    if ( lower>0 ) 
      h = G0 * h0 * (1-ch*e) * sqrt(Pat/p) / (lower*lower);
//      h = G0 * h0 * (1-ch*e) * sqrt(Pat/p);
//      h = h0;

    // Method 1
    temp_tensor = alpha_b_tensor - alpha; 
    
    // Method 2
    //temp_tensor = alpha_b_tensor+n*m - s*(1.0/p); 

    TensorEvolution::TensorEvolutionHij = temp_tensor * (h*r_ef);
     
    return TensorEvolution::TensorEvolutionHij;
}


//out- //================================================================================
//out- const tensor& SANISAND_alpha_Eij::DHij_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
//out-                             const straintensor& Stra, const MaterialParameter& material_parameter)
//out- {
//out-     const double oneOver3 = 1.0/3.0;    
//out-     const double rt23 = sqrt(2.0/3.0);
//out- 
//out-     tensor I2("I", 2, def_dim_2);    
//out-     tensor I4 = I2("ij")*I2("kl");
//out-     tensor I4s = ( I4.transpose0110() + I4.transpose0111() ) *0.5;
//out-     tensor I4dev = I4s - I4 * oneOver3;
//out-     
//out-     stresstensor a_a_in;
//out-     double a_in = 0.0;
//out-     double h = 0.0;
//out-     
//out-     double e0 = gete0(material_parameter);
//out-     double e_r = gete_r(material_parameter);
//out-     double lambda = getlambda(material_parameter);
//out-     double xi = getxi(material_parameter);
//out-     double Pat = getPat(material_parameter);
//out-     double m = getm(material_parameter);
//out-     double alpha_cc = getalpha_cc(material_parameter);
//out-     double c = getc(material_parameter);
//out-     double nb = getnb(material_parameter);
//out-     double h0 = geth0(material_parameter);
//out-     double ch = getch(material_parameter);
//out-     double G0 = getG0(material_parameter);        
//out- 
//out-     stresstensor alpha = getalpha(material_parameter);
//out- //    stresstensor z = getz(material_parameter);
//out- 
//out-     double p = Stre.p_hydrostatic();
//out-     stresstensor s = Stre.deviator();
//out- 
//out-     stresstensor n;
//out-     stresstensor alpha_b;
//out-     stresstensor alpha_b_alpha;
//out-     double b0 = 0.0;
//out-     double g = 0.0;
//out-     double ec = e_r;
//out-     double psi = 0.0;
//out-     double expnb = 1.0;
//out-     double ab = 0.0;
//out-     stresstensor s_bar;
//out-     double norm_s = 0.0;
//out- 
//out-     double J3D;
//out-     double cos3theta = 0.0;
//out-     
//out-     double e = e0 + (1.0 + e0) *Stra.Iinvariant1();
//out-     
//out-     s_bar = Stre.deviator() - (alpha *p);
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
//out-     psi = e - ec;
//out- 
//out-     expnb = exp( -nb *psi );
//out-     
//out-     // way 1
//out-     ab = g*alpha_cc*expnb - m;
//out-     alpha_b = n *(rt23*ab);
//out-     alpha_b_alpha = alpha_b - alpha;
//out-     
//out-     // way 2    
//out-     //ab = g *alpha_cc *expnb;
//out-     //alpha_b_alpha = n *(rt23*ab) - s *(1.0/p);
//out-     
//out-     if ( (p/Pat) > 0.0 )
//out-       b0 = G0 *h0 *(1.0-ch*e) *pow(p/Pat, -0.5);
//out-           
//out-     a_a_in = alpha - alpha_in;
//out-     a_in = (a_a_in("ij")*n("ij")).trace();
//out-     if (a_in < 0.0) 
//out-       alpha_in.Initialize(alpha);
//out-     if (a_in < aTOL)
//out-       a_in = aTOL;      
//out-     
//out-     h = b0 / a_in;
//out- 
//out-     double alpha_n = (alpha("ij")*n("ij")).trace();
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
//out-     // db0_ds
//out-     tensor db0_ds = I2 *(b0/6.0/p);
//out- 
//out-     // dphi_ds:
//out-     tensor dphi_ds = I2 * (-oneOver3*lambda*xi*pow(p, xi-1.0)/pow(Pat, xi));
//out- 
//out-     // dcos3theta_ds:
//out-     tensor dcos3theta_ds = dn_ds("ijmn")*n_n("ji");
//out-       dcos3theta_ds.null_indices();
//out-     dcos3theta_ds = dcos3theta_ds *(-3.0*sqrt(6.0));
//out- 
//out-     // dg_ds:
//out-     tensor dg_ds = dcos3theta_ds *(g*g*(1.0-c)/(2.0*c));
//out- 
//out-     // dab_ds:
//out-     tensor dab_ds = ( dg_ds - dphi_ds *(g*nb) )*(alpha_cc*expnb);
//out- 
//out-     // dabn_ds:
//out-     tensor tensor1 = n("ij")*dab_ds("mn");
//out-       tensor1.null_indices();
//out-     tensor dabn_ds = tensor1 + dn_ds *ab;
//out- 
//out-     // dh_ds:
//out-     tensor dh_ds = db0_ds *(1.0/a_in);
//out-     if (a_in > aTOL) {
//out-       tensor tensor2 = a_a_in("pq")*dn_ds("pqmn");
//out-       dh_ds = dh_ds - tensor2 *(h/a_in);
//out-     }
//out- 
//out-     tensor tensor3 = alpha_b_alpha("ij")*dh_ds("mn");
//out-     
//out-     TensorEvolution::TE_tensorR4 = (tensor3 + dabn_ds *(h*rt23)) *(2.0/3.0);
//out-      
//out-     return TensorEvolution::TE_tensorR4;
//out- }
//out- 
//out- 
//out- //================================================================================
//out- const tensor& SANISAND_alpha_Eij::DHij_Dkin(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
//out-                             const straintensor& Stra, const MaterialParameter& material_parameter)
//out- {
//out-     const double oneOver3 = 1.0/3.0;    
//out-     const double rt23 = sqrt(2.0/3.0);
//out- 
//out-     tensor I2("I", 2, def_dim_2);    
//out-     tensor I4 = I2("ij")*I2("kl");
//out-     tensor I4s = ( I4.transpose0110() + I4.transpose0111() ) *0.5;
//out-     tensor I4dev = I4s - I4 * oneOver3;
//out-     
//out-     stresstensor a_a_in;
//out-     double a_in = 0.0;
//out-     double h = 0.0;
//out-     
//out-     double e0 = gete0(material_parameter);
//out-     double e_r = gete_r(material_parameter);
//out-     double lambda = getlambda(material_parameter);
//out-     double xi = getxi(material_parameter);
//out-     double Pat = getPat(material_parameter);
//out-     double m = getm(material_parameter);
//out-     double alpha_cc = getalpha_cc(material_parameter);
//out-     double c = getc(material_parameter);
//out-     double nb = getnb(material_parameter);
//out-     double h0 = geth0(material_parameter);
//out-     double ch = getch(material_parameter);
//out-     double G0 = getG0(material_parameter);        
//out- 
//out-     stresstensor alpha = getalpha(material_parameter);
//out- //    stresstensor z = getz(material_parameter);
//out- 
//out-     double p = Stre.p_hydrostatic();
//out-     stresstensor s = Stre.deviator();
//out- 
//out-     stresstensor n;
//out-     stresstensor alpha_b;
//out-     stresstensor alpha_b_alpha;
//out-     double b0 = 0.0;
//out-     double g = 0.0;
//out-     double ec = e_r;
//out-     double psi = 0.0;
//out-     double expnb = 1.0;
//out-     double ab = 0.0;
//out-     stresstensor s_bar;
//out-     double norm_s = 0.0;
//out- 
//out-     double J3D;
//out-     double cos3theta = 0.0;
//out-     
//out-     double e = e0 + (1.0 + e0) *Stra.Iinvariant1();
//out-     
//out-     s_bar = Stre.deviator() - (alpha *p);
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
//out-     psi = e - ec;
//out- 
//out-     expnb = exp( -nb *psi );
//out-     
//out-     // way 1
//out-     ab = g*alpha_cc*expnb - m;
//out-     alpha_b = n *(rt23*ab);
//out-     alpha_b_alpha = alpha_b - alpha;
//out-     
//out-     // way 2    
//out-     //ab = g *alpha_cc *expnb;
//out-     //alpha_b_alpha = n *(rt23*ab) - s *(1.0/p);
//out-     
//out-     if ( (p/Pat) > 0.0 )
//out-       b0 = G0 *h0 *(1.0-ch*e) *pow(p/Pat, -0.5);
//out-           
//out-     a_a_in = alpha - alpha_in;
//out-     a_in = (a_a_in("ij")*n("ij")).trace();
//out-     if (a_in < 0.0) 
//out-       alpha_in.Initialize(alpha);
//out-     if (a_in < aTOL)
//out-       a_in = aTOL;      
//out-     
//out-     h = b0 / a_in;
//out- 
//out-     tensor n_n = n("ik")*n("kj");
//out-       n_n.null_indices();
//out- 
//out-     tensor nt_nt = n("ij")*n("kl");
//out-       nt_nt.null_indices();
//out-     
//out-     // dn_da:
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
//out-     // dab_da:
//out-     tensor dab_da = dg_da *(alpha_cc*expnb);
//out- 
//out-     // dabn_da:
//out-     tensor tensor1 = n("ij")*dab_da("mn");
//out-       tensor1.null_indices();
//out-     tensor dabn_da = tensor1 + dn_da *ab;
//out-     
//out-     TensorEvolution::TE_tensorR4 = (dabn_da *rt23 - I4s) *(h*2.0/3.0);
//out-     
//out-     // dh_da:
//out-     if (a_in > aTOL) {
//out-       tensor tensor2 = a_a_in("pq")*dn_da("pqmn");
//out-         tensor2.null_indices();
//out-       tensor dh_da = (n + tensor2) *(-h/a_in);
//out- 
//out-       tensor tensor3 = alpha_b_alpha("ij")*dh_da("mn");
//out-         tensor3.null_indices();
//out-       
//out-       TensorEvolution::TE_tensorR4 += tensor3 *(2.0/3.0);
//out-     }
//out-      
//out-     return TensorEvolution::TE_tensorR4;
//out- }


// to get e0
//================================================================================
double SANISAND_alpha_Eij::gete0(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, e0_index);
}

// to get e_r
//================================================================================
double SANISAND_alpha_Eij::gete_r(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, e_r_index);
}

// to get lambda
//================================================================================
double SANISAND_alpha_Eij::getlambda(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, lambda_index);
}

// to get xi
//================================================================================
double SANISAND_alpha_Eij::getxi(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, xi_index);

}

// to get Pat
//================================================================================
double SANISAND_alpha_Eij::getPat(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, Pat_index);
}

// to get alpha_cc
//================================================================================
double SANISAND_alpha_Eij::getalpha_cc(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, alpha_cc_index);
}

// to get c
//================================================================================
double SANISAND_alpha_Eij::getc(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, c_index);
}

// to get n_d
//================================================================================
double SANISAND_alpha_Eij::getnb(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, nb_index);
}


// to get h0
//================================================================================
double SANISAND_alpha_Eij::geth0(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, h0_index);

}


// to get ch
//================================================================================
double SANISAND_alpha_Eij::getch(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, ch_index);

}

// to get G0
//================================================================================
double SANISAND_alpha_Eij::getG0(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, G0_index);

}

// to get m
//================================================================================
double SANISAND_alpha_Eij::getm(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, m_index);
}


// to get alpha
//================================================================================
const stresstensor& SANISAND_alpha_Eij::getalpha(const MaterialParameter& material_parameter) const
{
	if ( alpha_index <= material_parameter.getNum_Internal_Tensor() && alpha_index > 0) {
		SANISAND_alpha_Eij::SANISAND_alpha_t = material_parameter.getInternal_Tensor(alpha_index-1);
		return SANISAND_alpha_Eij::SANISAND_alpha_t;
	}
	else {
		opserr << "SANISAND_alpha: Invalid Input (alpha) " << endln;
		exit (1);
	}
}


//================================================================================
double SANISAND_alpha_Eij::getParameters(const MaterialParameter& material_parameter, int which) const
{
	if ( which <= material_parameter.getNum_Material_Constant() && which > 0)
		return material_parameter.getMaterial_Constant(which-1);
	else {
		opserr << "SANISAND_alpha: Invalid Input - #" << which << endln;
		exit (1);
	}
} 


//================================================================================
double SANISAND_alpha_Eij::getec(double e_r, double lambda, double xi, double Pat, double p_c) const
{
    double ee = e_r; 
    
    if ( p_c/Pat >= 0.0 )
      ee = e_r - lambda * pow(p_c/Pat, xi);
    else 
      opserr << "Warning: SANISAND_alpha_Eij - 'p_c/Pat' less than zero! " << endln;

    return ee;
}

//================================================================================
double SANISAND_alpha_Eij::getg(double c, double cos3theta) const
{
    return 2.0 * c / ( (1.0+c) - (1.0-c)*cos3theta );
}


//Guanzhou added for parallel
int SANISAND_alpha_Eij::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(13);
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
    idData(12) = alpha_index;   
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "SANISAND_alpha_Eij::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int SANISAND_alpha_Eij::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(13);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "SANISAND_alpha_Eij::recvSelf -- failed to recv ID\n";
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
    alpha_index    = idData(12); 
    
    return 0;
}

#endif

