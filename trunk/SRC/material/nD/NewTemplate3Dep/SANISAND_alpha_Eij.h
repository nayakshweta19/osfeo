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

#ifndef SANISAND_alpha_Eij_H
#define SANISAND_alpha_Eij_H 

#include "TensorEvolution.h"
#define TE_TAG_SANISAND_alpha_Eij 121006

class SANISAND_alpha_Eij : public TensorEvolution
{
  public:


//! Reference: Taiebat & Dafalias (2008)
//! 
//! Parameters:
//!  1- e0:        initial void ratio at zero strain;
//!  2- e_r:       reference void for critical state line, ec = e_r - lambda*(pc/Pat)^xi;
//!  3- lambda:    parameter for critical state line;
//!  4- xi:        parameter for critical state line;
//!  5- Pat:       atmospherics pressure for critical state line;
//!  6- alpha_cc:  critical state stress ration;
//!  7- c:         tension-compression strength ratio;
//!  8- nb:        bounding parameter;
//!  9- h0:        bounding parameter;
//! 10- ch:        bounding parameter;
//! 11- G0:        parameter in the elastic part
//! 12- m:         opening of the yield surface;
//! 13- alpha:     back stress ratio tensor in yield function; (the 1st tensorial internal variable);
//! 
//! GENERAL NOTE
//! material parameter index: 
//! the location of the parameter in the material parameter part list
//! 
//! Tensor Evolution function for alpha in SANISAND:
//! inputs:
//! - e0_index_in:         the location of the parameter e0       in the material parameter part list
//! - e_r_index_in:        the location of the parameter e_r      in the material parameter part list
//! - lambda_c_index_in:   the location of the parameter lambda   in the material parameter part list
//! - xi_index_in,:        the location of the parameter xi       in the material parameter part list
//! - Pat_index_in:        the location of the parameter Pat      in the material parameter part list
//! - alpha_cc_index_in    the location of the parameter alpha_cc in the material parameter part list
//! - c_index_in:          the location of the parameter c        in the material parameter part list
//! - nb_index_in:         the location of the parameter nb       in the material parameter part list
//! - h0_index_in:         the location of the parameter h0       in the material parameter part list
//! - ch_index_in:         the location of the parameter ch       in the material parameter part list 
//! - G0_index_in:         the location of the parameter G0       in the material parameter part list
//! - m_index_in:          the location of the parameter m        in the material parameter part list
//! - alpha_index_in:      the location of the parameter alpha    in the material parameter part list



    SANISAND_alpha_Eij(int e0_index_in,	 
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
                       int alpha_index_in);

    SANISAND_alpha_Eij() : TensorEvolution(TE_TAG_SANISAND_alpha_Eij) {}; //Guanzhou added for parallel processing

    TensorEvolution* newObj();


//! Hij: 
//! inputs: 
//! - PlasticFlow& plastic_flow:
//! - stresstensor &Stre: the tensor of stress
//! - straintensor &Stra: the tensor of strain
//! - MaterialParameter &MaterialParameter_in: the class of material parameter

    const straintensor& Hij(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                            const straintensor& Stra, const MaterialParameter& material_parameter);

//out-     const tensor& DHij_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
//out-                           const straintensor& Stra, const MaterialParameter& material_parameter);
//out- 
//out-     const tensor& DHij_Dkin(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
//out-                             const straintensor& Stra, const MaterialParameter& material_parameter);
       
    //Guanzhou added for parallel
//! sendSelf: for parallel computing
    int sendSelf(int commitTag, Channel &theChannel);  
//! recvSelf: for parallel computing
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

  private:

//! gete0: to get e0
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double gete0(const MaterialParameter& material_parameter) const; 

//! gete_r: to get e_r
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double gete_r(const MaterialParameter& material_parameter) const;      

//! getlambda: to get lambda
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getlambda(const MaterialParameter& material_parameter) const;

//! getxi: to get xi
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getxi(const MaterialParameter& material_parameter) const;

//! getPat: to get pat (constant of atmosphereic pressure)
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getPat(const MaterialParameter& material_parameter) const;

//! getalpha_cc: to get alpha_cc
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getalpha_cc(const MaterialParameter& material_parameter) const;

//! getc: to get c
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getc(const MaterialParameter& material_parameter) const;


//! getnb: to get nb
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getnb(const MaterialParameter& material_parameter) const;


//! geth0: to get h0
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double geth0(const MaterialParameter& material_parameter) const;


//! getch: to get ch
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getch(const MaterialParameter& material_parameter) const;

//! "Mahdi"
//! getG0: to get G0
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getG0(const MaterialParameter& material_parameter) const;


    double getm(const MaterialParameter& material_parameter) const;
    const stresstensor& getalpha(const MaterialParameter& material_parameter) const;
    
    inline double getParameters(const MaterialParameter& material_parameter, int which_in) const;    
    inline double getec(double e_r, double lambda, double xi, double Pat, double p_c) const;
    inline double getg(double c, double cos3theta) const;
    
  private:
  
    // things need to be memorized
    //int a_index;
    // the following line should be deleted after correction of the implementation for the implicit method
    //stresstensor alpha_in;

  private:
  
    int e0_index;
    int e_r_index;
    int lambda_index;
    int xi_index;
    int Pat_index;
    int alpha_cc_index;
    int c_index;
    int nb_index;
    int h0_index;    
    int ch_index;
    int G0_index;
    int m_index;
    int alpha_index;
    
    static stresstensor SANISAND_alpha_t;

};


#endif

