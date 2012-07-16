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
// UPDATE HISTORY:    
//
///////////////////////////////////////////////////////////////////////////////
//

// Ref: Dafalias and Manzari 2004: J. Eng. Mech. 130(6), pp 622-634
// Parameters:
//  1- e0:        initial void ratio;
//  2- e_r:       reference void for critical state line, ec = e_r - lambda_c*(pc/Pat)^xi;
//  3- lambda_c:  parameter for critical state line, ec = e_r - lambda_c*(pc/Pat)^xi;
//  4- xi:        parameter for critical state line, ec = e_r - lambda_c*(pc/Pat)^xi;
//  5- Pat:       atmospherics pressure for critical state line, ec = e0 - lambda_c*(pc/Pat)^xi;
//  6- m:         parameter in the yield function;
//  7- M:         critical state stress ration;
//  8- cc:        tension-compression strength ratio;
//  9- nb:        parameter;
// 10- h0:        parameter;
// 11- ch:        parameter;
// 12- G0:        parameter in the elastic part
// 13- alpha:     "back-stress ratio" tensor in yield function; (the 1st tensorial internal variable);
// 14- z:         fabric dilation internal tensor (the 2nd tensorial internal variable); 

#ifndef DM04_alpha_Eij_CPP
#define DM04_alpha_Eij_CPP

#include "DM04_alpha_Eij.h"

stresstensor DM04_alpha_Eij::DM04_alpha_t;

const double aTOL = 1.0e-12;

DM04_alpha_Eij::DM04_alpha_Eij(int e0_index_in,
                               int e_r_index_in,
                               int lambda_c_index_in,
                               int xi_index_in,
                               int Pat_index_in,
                               int m_index_in,
                               int M_cal_index_in,
                               int cc_index_in,
                               int nb_index_in,
                               int h0_index_in,
                               int ch_index_in,
                               int G0_index_in,
                               int alpha_index_in,
                               int z_index_in)
: TensorEvolution(TE_TAG_DM04_alpha_Eij), e0_index(e0_index_in),
  e_r_index(e_r_index_in), 
  lambda_c_index(lambda_c_index_in),
  xi_index(xi_index_in),
  Pat_index(Pat_index_in),
  m_index(m_index_in),  
  M_cal_index(M_cal_index_in),
  cc_index(cc_index_in),
  nb_index(nb_index_in),
  h0_index(h0_index_in),   
  ch_index(ch_index_in),
  G0_index(G0_index_in),
  alpha_index(alpha_index_in),
  z_index(z_index_in) 
{
   stresstensor zT;
   alpha_in.Initialize(zT);
}

TensorEvolution* DM04_alpha_Eij::newObj()
{
    TensorEvolution* nObj = new DM04_alpha_Eij(this->e0_index,
                                               this->e_r_index,
                                               this->lambda_c_index,
                                               this->xi_index,
                                               this->Pat_index,
                                               this->m_index,
                                               this->M_cal_index,
                                               this->cc_index,
                                               this->nb_index,
                                               this->h0_index,
                                               this->ch_index,
                                               this->G0_index,
                                               this->alpha_index,
                                               this->z_index);
    return nObj;
}

const straintensor& DM04_alpha_Eij::Hij(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                                        const straintensor& Stra, const MaterialParameter& material_parameter)
{
    const double rt23 = sqrt(2.0/3.0);    

    stresstensor a_a_in;
    double a_in = 0.0;
    double h = 0.0;
    
    double e0 = gete0(material_parameter);
    double e_r = gete_r(material_parameter);
    double lambda_c = getlambda_c(material_parameter);
    double xi = getxi(material_parameter);
    double Pat = getPat(material_parameter);
    double m = getm(material_parameter);
    double M_cal = getM_cal(material_parameter);
    double cc = getcc(material_parameter);
    double nb = getnb(material_parameter);
    double h0 = geth0(material_parameter);
    double ch = getch(material_parameter);
    double G0 = getG0(material_parameter);        

    stresstensor alpha = getalpha(material_parameter);
    stresstensor z = getz(material_parameter);

    double p = Stre.p_hydrostatic();
    stresstensor s = Stre.deviator();

    stresstensor n;
    stresstensor alpha_b;
    stresstensor alpha_b_alpha;
    double b0 = 0.0;
    double g = 0.0;
    double ec = e_r;
    double stateParameter = 0.0;
    double expnb = 1.0;
    double ab = 0.0;
    stresstensor s_bar;
    double norm_s = 0.0;

    double J3D;
    double cos3theta = 0.0;
    
    double e = e0 + (1.0 + e0) *Stra.Iinvariant1();
    
    s_bar = Stre.deviator() - (alpha *p);
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

    stateParameter = e - ec;

    expnb = exp( -nb *stateParameter );
    
    // way 1
    ab = g*M_cal*expnb - m;
    alpha_b = n *(rt23*ab);
    alpha_b_alpha = alpha_b - alpha;
    
    // way 2    
    //ab = g *M_cal *expnb;
    //alpha_b_alpha = n *(rt23*ab) - s *(1.0/p);
    
    if ( (p/Pat) > 0.0 )
      b0 = G0 *h0 *(1.0-ch*e) *pow(p/Pat, -0.5);

    //if (a_index == 0) {
    //  alpha_in.Initialize(alpha);
    //  a_in = 0.0;
    //  a_index = 1;
    //  h = b0 / aTOL;
    //}
    //else {
    //  a_a_in = alpha - alpha_in;
    //  a_in = (a_a_in("ij")*n("ij")).trace();
    //  if ( a_in < 0.0 ) {
    //    a_index = 0;
    //    alpha_in.Initialize(alpha);
    //  }
    //  if ( a_in < aTOL )
    //    a_in = aTOL;
    //  h = b0 /a_in;
    //}
    
    a_a_in = alpha - alpha_in;
    a_in = (a_a_in("ij")*n("ij")).trace();

    if (a_in < 0.0) 
      alpha_in.Initialize(alpha);
    if (a_in < aTOL)
      a_in = aTOL;      
    
    h = b0 / a_in;
   
    TensorEvolution::TensorEvolutionHij = alpha_b_alpha *(h*2.0/3.0);
     
    return TensorEvolution::TensorEvolutionHij;
}


//================================================================================
const tensor& DM04_alpha_Eij::DHij_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                            const straintensor& Stra, const MaterialParameter& material_parameter)
{
    const double oneOver3 = 1.0/3.0;    
    const double rt23 = sqrt(2.0/3.0);

    tensor I2("I", 2, def_dim_2);    
    tensor I4 = I2("ij")*I2("kl");
    tensor I4s = ( I4.transpose0110() + I4.transpose0111() ) *0.5;
    tensor I4dev = I4s - I4 * oneOver3;
    
    stresstensor a_a_in;
    double a_in = 0.0;
    double h = 0.0;
    
    double e0 = gete0(material_parameter);
    double e_r = gete_r(material_parameter);
    double lambda_c = getlambda_c(material_parameter);
    double xi = getxi(material_parameter);
    double Pat = getPat(material_parameter);
    double m = getm(material_parameter);
    double M_cal = getM_cal(material_parameter);
    double cc = getcc(material_parameter);
    double nb = getnb(material_parameter);
    double h0 = geth0(material_parameter);
    double ch = getch(material_parameter);
    double G0 = getG0(material_parameter);        

    stresstensor alpha = getalpha(material_parameter);
    stresstensor z = getz(material_parameter);

    double p = Stre.p_hydrostatic();
    stresstensor s = Stre.deviator();

    stresstensor n;
    stresstensor alpha_b;
    stresstensor alpha_b_alpha;
    double b0 = 0.0;
    double g = 0.0;
    double ec = e_r;
    double stateParameter = 0.0;
    double expnb = 1.0;
    double ab = 0.0;
    stresstensor s_bar;
    double norm_s = 0.0;

    double J3D;
    double cos3theta = 0.0;
    
    double e = e0 + (1.0 + e0) *Stra.Iinvariant1();
    
    s_bar = Stre.deviator() - (alpha *p);
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

    stateParameter = e - ec;

    expnb = exp( -nb *stateParameter );
    
    // way 1
    ab = g*M_cal*expnb - m;
    alpha_b = n *(rt23*ab);
    alpha_b_alpha = alpha_b - alpha;
    
    // way 2    
    //ab = g *M_cal *expnb;
    //alpha_b_alpha = n *(rt23*ab) - s *(1.0/p);
    
    if ( (p/Pat) > 0.0 )
      b0 = G0 *h0 *(1.0-ch*e) *pow(p/Pat, -0.5);
          
    a_a_in = alpha - alpha_in;
    a_in = (a_a_in("ij")*n("ij")).trace();
    if (a_in < 0.0) 
      alpha_in.Initialize(alpha);
    if (a_in < aTOL)
      a_in = aTOL;      
    
    h = b0 / a_in;

    double alpha_n = (alpha("ij")*n("ij")).trace();

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
   
    // db0_ds
    tensor db0_ds = I2 *(b0/6.0/p);

    // dphi_ds:
    tensor dphi_ds = I2 * (-oneOver3*lambda_c*xi*pow(p, xi-1.0)/pow(Pat, xi));

    // dcos3theta_ds:
    tensor dcos3theta_ds = dn_ds("ijmn")*n_n("ji");
      dcos3theta_ds.null_indices();
    dcos3theta_ds = dcos3theta_ds *(-3.0*sqrt(6.0));

    // dg_ds:
    tensor dg_ds = dcos3theta_ds *(g*g*(1.0-cc)/(2.0*cc));

    // dab_ds:
    tensor dab_ds = ( dg_ds - dphi_ds *(g*nb) )*(M_cal*expnb);

    // dabn_ds:
    tensor tensor1 = n("ij")*dab_ds("mn");
      tensor1.null_indices();
    tensor dabn_ds = tensor1 + dn_ds *ab;

    // dh_ds:
    tensor dh_ds = db0_ds *(1.0/a_in);
    if (a_in > aTOL) {
      tensor tensor2 = a_a_in("pq")*dn_ds("pqmn");
      dh_ds = dh_ds - tensor2 *(h/a_in);
    }

    tensor tensor3 = alpha_b_alpha("ij")*dh_ds("mn");
    
    TensorEvolution::TE_tensorR4 = (tensor3 + dabn_ds *(h*rt23)) *(2.0/3.0);
     
    return TensorEvolution::TE_tensorR4;
}


//================================================================================
const tensor& DM04_alpha_Eij::DHij_Dkin(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                            const straintensor& Stra, const MaterialParameter& material_parameter)
{
    const double oneOver3 = 1.0/3.0;    
    const double rt23 = sqrt(2.0/3.0);

    tensor I2("I", 2, def_dim_2);    
    tensor I4 = I2("ij")*I2("kl");
    tensor I4s = ( I4.transpose0110() + I4.transpose0111() ) *0.5;
    tensor I4dev = I4s - I4 * oneOver3;
    
    stresstensor a_a_in;
    double a_in = 0.0;
    double h = 0.0;
    
    double e0 = gete0(material_parameter);
    double e_r = gete_r(material_parameter);
    double lambda_c = getlambda_c(material_parameter);
    double xi = getxi(material_parameter);
    double Pat = getPat(material_parameter);
    double m = getm(material_parameter);
    double M_cal = getM_cal(material_parameter);
    double cc = getcc(material_parameter);
    double nb = getnb(material_parameter);
    double h0 = geth0(material_parameter);
    double ch = getch(material_parameter);
    double G0 = getG0(material_parameter);        

    stresstensor alpha = getalpha(material_parameter);
    stresstensor z = getz(material_parameter);

    double p = Stre.p_hydrostatic();
    stresstensor s = Stre.deviator();

    stresstensor n;
    stresstensor alpha_b;
    stresstensor alpha_b_alpha;
    double b0 = 0.0;
    double g = 0.0;
    double ec = e_r;
    double stateParameter = 0.0;
    double expnb = 1.0;
    double ab = 0.0;
    stresstensor s_bar;
    double norm_s = 0.0;

    double J3D;
    double cos3theta = 0.0;
    
    double e = e0 + (1.0 + e0) *Stra.Iinvariant1();
    
    s_bar = Stre.deviator() - (alpha *p);
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

    stateParameter = e - ec;

    expnb = exp( -nb *stateParameter );
    
    // way 1
    ab = g*M_cal*expnb - m;
    alpha_b = n *(rt23*ab);
    alpha_b_alpha = alpha_b - alpha;
    
    // way 2    
    //ab = g *M_cal *expnb;
    //alpha_b_alpha = n *(rt23*ab) - s *(1.0/p);
    
    if ( (p/Pat) > 0.0 )
      b0 = G0 *h0 *(1.0-ch*e) *pow(p/Pat, -0.5);
          
    a_a_in = alpha - alpha_in;
    a_in = (a_a_in("ij")*n("ij")).trace();
    if (a_in < 0.0) 
      alpha_in.Initialize(alpha);
    if (a_in < aTOL)
      a_in = aTOL;      
    
    h = b0 / a_in;

    tensor n_n = n("ik")*n("kj");
      n_n.null_indices();

    tensor nt_nt = n("ij")*n("kl");
      nt_nt.null_indices();
    
    // dn_da:
    tensor dn_da = nt_nt - I4s;
    dn_da = dn_da *(p/norm_s);

    // dcos3theta_dalpha:
    tensor dcos3theta_da = dn_da("ijmn")*n_n("ji");
      dcos3theta_da.null_indices();
    dcos3theta_da = dcos3theta_da *(-3.0*sqrt(6.0));

    // dg_da:
    tensor dg_da = dcos3theta_da *(g*g*(1.0-cc)/(2.0*cc));

    // dab_da:
    tensor dab_da = dg_da *(M_cal*expnb);

    // dabn_da:
    tensor tensor1 = n("ij")*dab_da("mn");
      tensor1.null_indices();
    tensor dabn_da = tensor1 + dn_da *ab;
    
    TensorEvolution::TE_tensorR4 = (dabn_da *rt23 - I4s) *(h*2.0/3.0);
    
    // dh_da:
    if (a_in > aTOL) {
      tensor tensor2 = a_a_in("pq")*dn_da("pqmn");
        tensor2.null_indices();
      tensor dh_da = (n + tensor2) *(-h/a_in);

      tensor tensor3 = alpha_b_alpha("ij")*dh_da("mn");
        tensor3.null_indices();
      
      TensorEvolution::TE_tensorR4 += tensor3 *(2.0/3.0);
    }
     
    return TensorEvolution::TE_tensorR4;
}


// to get e0
//================================================================================
double DM04_alpha_Eij::gete0(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, e0_index);
}

// to get e_r
//================================================================================
double DM04_alpha_Eij::gete_r(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, e_r_index);
}

// to get lambda_c
//================================================================================
double DM04_alpha_Eij::getlambda_c(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, lambda_c_index);
}

// to get xi
//================================================================================
double DM04_alpha_Eij::getxi(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, xi_index);

}

// to get Pat
//================================================================================
double DM04_alpha_Eij::getPat(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, Pat_index);
}

// to get m
//================================================================================
double DM04_alpha_Eij::getm(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, m_index);
}

// to get M
//================================================================================
double DM04_alpha_Eij::getM_cal(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, M_cal_index);
}

// to get c
//================================================================================
double DM04_alpha_Eij::getcc(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, cc_index);
}

// to get n_d
//================================================================================
double DM04_alpha_Eij::getnb(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, nb_index);
}


// to get h0
//================================================================================
double DM04_alpha_Eij::geth0(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, h0_index);

}


// to get ch
//================================================================================
double DM04_alpha_Eij::getch(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, ch_index);

}

// to get G0
//================================================================================
double DM04_alpha_Eij::getG0(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, G0_index);

}

// to get alpha
//================================================================================
const stresstensor& DM04_alpha_Eij::getalpha(const MaterialParameter& material_parameter) const
{
	if ( alpha_index <= material_parameter.getNum_Internal_Tensor() && alpha_index > 0) {
		DM04_alpha_Eij::DM04_alpha_t = material_parameter.getInternal_Tensor(alpha_index-1);
		return DM04_alpha_Eij::DM04_alpha_t;
	}
	else {
		opserr << "DM04_alpha: Invalid Input (alpha) " << endln;
		exit (1);
	}
}

// to get z
//================================================================================
const stresstensor& DM04_alpha_Eij::getz(const MaterialParameter& material_parameter) const
{
    if ( z_index <= material_parameter.getNum_Internal_Tensor() && z_index > 0) {
		DM04_alpha_Eij::DM04_alpha_t = material_parameter.getInternal_Tensor(z_index-1);
		return DM04_alpha_Eij::DM04_alpha_t;
	}
	else {
		opserr << "DM04_alpha: Invalid Input (z) " << endln;
		exit (1);
	}
}

//================================================================================
double DM04_alpha_Eij::getParameters(const MaterialParameter& material_parameter, int which) const
{
	if ( which <= material_parameter.getNum_Material_Constant() && which > 0)
		return material_parameter.getMaterial_Constant(which-1);
	else {
		opserr << "DM04_alpha: Invalid Input - #" << which << endln;
		exit (1);
	}
} 


//================================================================================
double DM04_alpha_Eij::getec(double e_r, double lambda_c, double xi, double Pat, double p_c) const
{
    double ee = e_r; 
    
    if ( p_c/Pat >= 0.0 )
      ee = e_r - lambda_c * pow(p_c/Pat, xi);
    else 
      opserr << "Warning: DM04_alpha_Eij - 'p_c/Pat' less than zero! " << endln;

    return ee;
}

//================================================================================
double DM04_alpha_Eij::getg(double c, double cos3theta) const
{
    return 2.0 * c / ( (1.0+c) - (1.0-c)*cos3theta );
}


//Guanzhou added for parallel
int DM04_alpha_Eij::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(14);
    idData.Zero();

    idData(0) = e0_index;      
    idData(1) = e_r_index;     
    idData(2) = lambda_c_index;
    idData(3) = xi_index;      
    idData(4) = Pat_index;     
    idData(5) = m_index;       
    idData(6) = M_cal_index;   
    idData(7) = cc_index;      
    idData(8) = nb_index;      
    idData(9) = h0_index;      
    idData(10) = ch_index;      
    idData(11) = G0_index;      
    idData(12) = alpha_index;   
    idData(13) = z_index;       
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "DM04_alpha_Eij::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int DM04_alpha_Eij::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(14);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "DM04_alpha_Eij::recvSelf -- failed to recv ID\n";
	return -1;
    }

    e0_index       = idData(0);
    e_r_index      = idData(1);
    lambda_c_index = idData(2);
    xi_index       = idData(3);
    Pat_index      = idData(4);
    m_index        = idData(5);
    M_cal_index    = idData(6);
    cc_index       = idData(7);
    nb_index       = idData(8);
    h0_index       = idData(9);
    ch_index     = idData(10);      	                 
    G0_index     = idData(11); 
    alpha_index  = idData(12); 
    z_index      = idData(13);
    
    return 0;
}

#endif

