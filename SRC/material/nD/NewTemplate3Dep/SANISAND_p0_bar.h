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
// UPDATE HISTORY:    Nima Tafazzoli updated for API (Mar 2009)
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef SANISAND_p0_bar_H
#define SANISAND_p0_bar_H

#include "ScalarEvolution.h"
#define SE_TAG_SANISAND_p0_bar 121003

class SANISAND_p0_bar : public ScalarEvolution
{
  public:

//! Reference: Taiebat & Dafalias (2008)
//! 
//! Parameters:
//!  1- e0:        initial void ratio at zero strain;
//!  2- Pat:       atmospherics pressure for critical state line;
//!  3- alpha_cc:  critical state stress ration;
//!  4- p_r:       LCC parameter; 
//!  5- rho_c:     LCC parameter; 
//!  6- theta_c    LCC parameter; 
//!  7- K0:        Reference elastic bulk modulus (same unit as stress); 
//!  8- p0:        yield surface size;
//!  9- alpha:     back stress ratio tensor in yield function; (the 1st tensorial internal variable);
//! 
//! GENERAL NOTE
//! material parameter index: 
//! the location of the parameter in the material parameter part list
//! 
//! Scalar Evolution function for SANISAND:
//! inputs:
//! - e0_index_in:         the location of the parameter e0       in the material parameter part list
//! - Pat_index_in:        the location of the parameter Pat      in the material parameter part list
//! - alpha_cc_index_in:   the location of the parameter alpha_cc in the material parameter part list
//! - p_r_index_in:        the location of the parameter p_r      in the material parameter part list
//! - rho_c_index_in:      the location of the parameter rho_c    in the material parameter part list
//! - theta_c_index_in:    the location of the parameter theta_c  in the material parameter part list
//! - K0_index_in:         the location of the parameter K0       in the material parameter part list
//! - p0_index_in:         the location of the parameter p0       in the material parameter part list
//! - alpha_index_in:      the location of the parameter alpha    in the material parameter part list


    SANISAND_p0_bar(int e0_index_in,	 
                    int Pat_index_in,	 
                    int alpha_cc_index_in,	 
                    int p_r_index_in,	 
                    int rho_c_index_in,	 
                    int theta_c_index_in,	 
                    int K0_index_in,	 
                    int p0_index_in,	 
                    int alpha_index_in);

    SANISAND_p0_bar() : ScalarEvolution(SE_TAG_SANISAND_p0_bar) {}; //Guanzhou added for parallel processing

    ScalarEvolution* newObj();

//! H: 
//! inputs: 
//! - PlasticFlow& plastic_flow:
//! - stresstensor &Stre: the tensor of stress
//! - straintensor &Stra: the tensor of strain
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
    double H(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
             const straintensor& Stra, const MaterialParameter& material_parameter);

//out-     const tensor& DH_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
//out-                    const straintensor& Stra, const MaterialParameter& material_parameter);
//out- 
//out-     double DH_Diso(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
//out-                    const straintensor& Stra, const MaterialParameter& material_parameter);
 
    //Guanzhou added for parallel
//! sendSelf: for parallel computing
    int sendSelf(int commitTag, Channel &theChannel);
//! recvSelf: for parallel computing
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
      
  private:
 

    double gete0(const MaterialParameter& material_parameter) const; 
    double getPat(const MaterialParameter& material_parameter) const;
    double getalpha_cc(const MaterialParameter& material_parameter) const;

//! Reference: Taiebat & Dafalias (2008)
//! getp_r: to get p_r
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getp_r(const MaterialParameter& material_parameter) const;

//! Reference: Taiebat & Dafalias (2008)
//! getrho_c: to get rho_c
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getrho_c(const MaterialParameter& material_parameter) const;

//! Reference: Taiebat & Dafalias (2008)
//! gettheta_c: to get theta_c
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double gettheta_c(const MaterialParameter& material_parameter) const;

//! Reference: Taiebat & Dafalias (2008)
//! getK0: to get K0
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getK0(const MaterialParameter& material_parameter) const;

//! Reference: Taiebat & Dafalias (2008)
//! getp0: to get p0
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getp0(const MaterialParameter& material_parameter) const;

    const stresstensor& getalpha(const MaterialParameter& material_parameter) const;
    
    inline double getParameters(const MaterialParameter& material_parameter, int which_in) const;    
    
  private:

    int e0_index;
    int Pat_index;
    int alpha_cc_index;
    int p_r_index;  
    int rho_c_index;
    int theta_c_index;
    int K0_index;   
    int p0_index;
    int alpha_index;

    static stresstensor SANISAND_p0_bar_t;
      
};


#endif

   
	 

    
    
 
