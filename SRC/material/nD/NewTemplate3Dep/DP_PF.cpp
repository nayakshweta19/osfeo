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
// UPDATE HISTORY:    Guanzhou Jie updated for parallel, Dec 2006
//                    Nima Tafazzoli added for II, Feb 2010
///////////////////////////////////////////////////////////////////////////////
//

#ifndef DP_PF_CPP
#define DP_PF_CPP

#include "DP_PF.h"
#define PF_TAG_DP 120003

straintensor DP_PF::DPm;
stresstensor DP_PF::DPa;

//================================================================================
DP_PF::DP_PF(int m_which_in, int index_m_in, int alpha_which_in, int index_alpha_in)
: PlasticFlow(PF_TAG_DP), m_which(m_which_in), index_m(index_m_in), 
  alpha_which(alpha_which_in), index_alpha(index_alpha_in)
{

}

//================================================================================
DP_PF::~DP_PF() 
{  

}

//================================================================================
PlasticFlow* DP_PF::newObj() 
{  
     PlasticFlow  *new_PF = new DP_PF(this->m_which, this->index_m, 
                                      this->alpha_which, this->index_alpha);
     
     return new_PF;
}

//================================================================================
const straintensor& DP_PF::PlasticFlowTensor(const stresstensor &Stre, 
                                             const straintensor &Stra, 
                                             const MaterialParameter &MaterialParameter_in)  const
{

// Nima Tafazzoli added form (II) lecture notes PP. 99
// Febuary 2010
// f = [(sij-p*aij)(sij-p*aij)]^0.5 - sqrt(2/3)*(p*m) = 0

    tensor I2("I", 2, def_dim_2);
    double M = getm(MaterialParameter_in);
    double p = Stre.p_hydrostatic();
    stresstensor alpha;
     
    if (alpha_which == 2) 
       { 
         alpha = getalpha(MaterialParameter_in);
       }
 
    stresstensor s_bar = Stre - alpha*p;
    s_bar = s_bar.deviator();

//      s_bar.print("s_bar" , "s_bar");


    double s_norm = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
    stresstensor nij;
    if (s_norm != 0.0)
      nij =  s_bar *(1.0/s_norm);
    double a_n = (alpha("ij")*nij("ij")).trace();

    DP_PF::DPm = nij + I2 * (M*sqrt(2.0/27.0) + a_n/3.0);

    return DP_PF::DPm;

/*
    // g = 1.5*(sij-p*aij)(sij-p*aij) - (m*p)^2 = 0

    tensor I2("I", 2, def_dim_2);
    double m = getm(MaterialParameter_in);
    double p = Stre.p_hydrostatic();
    stresstensor alpha;
    stresstensor sigma_b;
    if (alpha_which == 2) 
      { 
        stresstensor alpha = getalpha(MaterialParameter_in);
        sigma_b = Stre - (alpha * p);
      }
    else 
      {
        sigma_b = Stre;
      }
    stresstensor s_bar = sigma_b.deviator();
    double scalar1 = ( s_bar("ij") * alpha("ij") ).trace();
    DP_PF::DPm = s_bar *3.0 + I2 *(scalar1 + m*m*p*2.0/3.0);

    return DP_PF::DPm;
*/

}

///////////////////////////////////////////////////////////////////////////////
const tensor& DP_PF::Dm_Ds(const stresstensor &Stre, 
                           const straintensor &Stra, 
                           const MaterialParameter &MaterialParameter_in) const  
{

// Nima Tafazzoli added form (II) lecture notes PP. 99
// Febuary 2010
// f = [(sij-p*aij)(sij-p*aij)]^0.5 - sqrt(2/3)*(p*m) = 0

    stresstensor alpha = getalpha(MaterialParameter_in);
    double p = Stre.p_hydrostatic();
    stresstensor s_bar = Stre - alpha*p;
    s_bar = s_bar.deviator();
    double sijsij = (s_bar("ij")*s_bar("ij")).trace();
    double s_norm = sqrt(sijsij);

    tensor I_ijkl;
    I_ijkl = I_ijkl.I_1234();

    tensor I_kilj;
    I_kilj = I_kilj.I_3142();

    tensor I2("I", 2, def_dim_2);
    tensor I_alpha = I2("kl")*alpha("ij");
    I_alpha.null_indices();

    tensor alpha_I = alpha("kl")*I2("ij");
    alpha_I.null_indices();

    double alpha_alpha = (alpha("pq")*alpha("pq")).trace();

    double alpha_delta = alpha.cval(1,1)+alpha.cval(2,2)+alpha.cval(3,3);

    double s_delta = s_bar.cval(1,1)+s_bar.cval(2,2)+s_bar.cval(3,3);
	
    double alpha_s = (alpha("ij")*s_bar("ij")).trace();

	stresstensor temp = s_bar + I2* ((1.0/3.0)*alpha_s);
	
	stresstensor temp2 = s_bar;

    stresstensor temp3 = I2* ((-1.0/3.0)*s_delta);

	stresstensor temp4 = I2*((1.0/3.0)*alpha_s);

	stresstensor temp5 = temp2 + temp3 + temp4;
	
	stresstensor temp6 = temp("ij")*temp5("mn");
	      temp6.null_indices();
		  
    PlasticFlow::PF_tensorR4 = ( I_kilj - I_ijkl*(1.0/3.0) + I_alpha*(1.0/3.0) )*(1.0/s_norm) +
                               ( alpha_I*(1.0/3.0) - I_ijkl*((1.0/9.0)*alpha_delta) + I_ijkl*((1.0/9.0)*alpha_alpha) )* (1.0/s_norm) -
                                temp6*(1.0/(pow(sijsij,1.5)));

    return PlasticFlow::PF_tensorR4;


/*
    tensor I2("I", 2, def_dim_2);
    tensor I_ijkl = I2("ij")*I2("kl");
    I_ijkl.null_indices();
    tensor I_ikjl = I_ijkl.transpose0110();
    tensor I_iljk = I_ijkl.transpose0111();
    tensor I4s = (I_ikjl+I_iljk)*0.5;
    
    double m = getm(MaterialParameter_in);
    stresstensor alpha;
    if (alpha_which == 2)  
      stresstensor alpha = getalpha(MaterialParameter_in);
    double p = Stre.p_hydrostatic();
    stresstensor sigma_b = Stre - (alpha * p);
    stresstensor s_bar = sigma_b.deviator();

    tensor alpha_I = alpha("ij")*I2("mn");
    alpha_I.null_indices();
    tensor I_alpha = I2("ij")*alpha("mn");
    I_alpha.null_indices();
    double alpha_alpha = (alpha("pq")*alpha("pq")).trace();

    PlasticFlow::PF_tensorR4 =  I4s *3.0 + alpha_I + I_alpha + I_ijkl *(-1.0 + alpha_alpha/3.0 - 2.0*m*m/9.0);        
    
    return PlasticFlow::PF_tensorR4;
*/

}

///////////////////////////////////////////////////////////////////////////////
const tensor& DP_PF::Dm_Diso(const stresstensor &Stre, 
                             const straintensor &Stra, 
                             const MaterialParameter &MaterialParameter_in) const 
{

// Nima Tafazzoli added form (II) lecture notes PP. 99
// Febuary 2010
// f = [(sij-p*aij)(sij-p*aij)]^0.5 - sqrt(2/3)*(p*m) = 0

    tensor I2("I", 2, def_dim_2);    

    if (m_which == 1)  {
      PlasticFlow::PF_tensorR2 = I2 *sqrt(2.0/27.0);
    }
    else {
      stresstensor tensor1;
      PlasticFlow::PF_tensorR2.Initialize(tensor1);
    }
    
    return PlasticFlow::PF_tensorR2;

/*
    tensor I2("I", 2, def_dim_2);    

    if (m_which == 1)  {
      double m = getm(MaterialParameter_in);
      double p = Stre.p_hydrostatic();
      PlasticFlow::PF_tensorR2 = I2 *(m*p*4.0/3.0);
    }
    else {
      stresstensor tensor1;
      PlasticFlow::PF_tensorR2.Initialize(tensor1);
    }
    
    return PlasticFlow::PF_tensorR2;
*/

}

///////////////////////////////////////////////////////////////////////////////
const tensor& DP_PF::Dm_Dkin(const stresstensor &Stre, 
                             const straintensor &Stra, 
                             const MaterialParameter &MaterialParameter_in)  const 
{

// Nima Tafazzoli added form (II) lecture notes PP. 99
// Febuary 2010
// f = [(sij-p*aij)(sij-p*aij)]^0.5 - sqrt(2/3)*(p*m) = 0

    stresstensor alpha = getalpha(MaterialParameter_in);
    double p = Stre.p_hydrostatic();
    stresstensor s_bar = Stre - alpha*p;
    s_bar = s_bar.deviator();
    double sijsij = (s_bar("ij")*s_bar("ij")).trace();
    double s_norm = sqrt(sijsij);

    tensor I_kilj;
    I_kilj = I_kilj.I_3142();

    tensor I2("I", 2, def_dim_2);
    tensor alpha_I = alpha("kl")*I2("ij");
    alpha_I.null_indices();

    tensor I_S = I2("ij")*s_bar("kl");
    I_S.null_indices();

    tensor S_S = s_bar("ij")*s_bar("kl");
    S_S.null_indices();

    double alpha_s = (alpha("ij")*s_bar("ij")).trace();

    PlasticFlow::PF_tensorR4 = ( I_kilj *(-p) + I_S *(1.0/3.0) - alpha_I *((1.0/3.0)*p) ) *(1.0/s_norm) -
                      ( S_S *p +  I_S*(1.0/3.0)*p*alpha_s ) * (1.0/(pow(sijsij,1.5)));

    return PlasticFlow::PF_tensorR4;


/*
    if (alpha_which == 2)  {
      tensor I2("I", 2, def_dim_2);
      tensor I_ijkl = I2("ij")*I2("kl");
      I_ijkl.null_indices();
      tensor I_ikjl = I_ijkl.transpose0110();
      tensor I_iljk = I_ijkl.transpose0111();
      tensor I4s = (I_ikjl+I_iljk)*0.5;

      double p = Stre.p_hydrostatic();
	  stresstensor s = Stre.deviator();
      tensor tensor1 = I2("ij") * s("mn");
      tensor1.null_indices();        
    
      PlasticFlow::PF_tensorR4 = I4s *(-3.0*p) + tensor1;
    }
    else {
      tensor tensor2(4, def_dim_4, 0.0);
      PlasticFlow::PF_tensorR4.Initialize(tensor2);
    }
    
    return PlasticFlow::PF_tensorR4;
*/

}

//================================================================================
double DP_PF::getm(const MaterialParameter &MaterialParameter_in) const
{
	// to get m
	if ( m_which == 0) {
		if ( index_m <= MaterialParameter_in.getNum_Material_Constant() && index_m > 0)
			return MaterialParameter_in.getMaterial_Constant(index_m-1); 
		else {
			opserr << "DP_PF: Invalid Input. " << endln;
			exit (1);
		}
	}
	else if ( m_which == 1) {
		if ( index_m <= MaterialParameter_in.getNum_Internal_Scalar() && index_m > 0)
			return MaterialParameter_in.getInternal_Scalar(index_m-1); 
		else {
			opserr << "DP_PF: Invalid Input. " << endln;
			exit (1);
		}
    }
	else {
		opserr << "DP_PF: Invalid Input. " << endln;
		exit(1);
	}
}

//================================================================================ 
const stresstensor& DP_PF::getalpha(const MaterialParameter &MaterialParameter_in) const
{
	//to get alpha(backstress)
	if ( alpha_which == 2 && index_alpha <= MaterialParameter_in.getNum_Internal_Tensor() && index_alpha > 0) {
		DP_PF::DPa = MaterialParameter_in.getInternal_Tensor(index_alpha-1);
		return DP_PF::DPa;
	}
	else if ( alpha_which == -1 ) {
        DP_PF::DPa = DP_PF::DPa*0.0;
        return DP_PF::DPa;
    }    
    else {
		opserr << "DP_PF: Invalid Input. " << endln;
		exit (1);
	}
}

//Guanzhou added for parallel
int DP_PF::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(4);
    idData.Zero();

    idData(0) = m_which;
    idData(1) = index_m;
    idData(2) = alpha_which;   
    idData(3) = index_alpha;   
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "DP_PF::sendSelf -- failed to send ID\n";
   	return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int DP_PF::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(4);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "CC_PF::recvSelf -- failed to recv ID\n";
	return -1;
    }

    m_which   = idData(0);
    index_m   = idData(1);
    alpha_which   = idData(2);
    index_alpha   = idData(3);

    return 0;
}

#endif

