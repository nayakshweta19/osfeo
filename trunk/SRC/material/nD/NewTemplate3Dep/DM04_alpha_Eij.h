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

#ifndef DM04_alpha_Eij_H
#define DM04_alpha_Eij_H 

#include "TensorEvolution.h"
#define TE_TAG_DM04_alpha_Eij 122001

class DM04_alpha_Eij : public TensorEvolution
{
  public:
  
//!
//! rf: Dafalias-Manzari 2004
//! inputs:
//! - e0_index_in: to locate the position in the defined (constant) MaterialParameter command for e0 (initial void ratio)
//! - e_r_index_in: to locate the position in the defined (constant) MaterialParameter command for e_r
//! - lambda_c_index_in: to locate the position in the defined (constant) MaterialParameter command for lambda_c
//! - xi_index_in,: to locate the position in the defined (constant) MaterialParameter command for xi
//! - Pat_index_in: to locate the position in the defined (constant) MaterialParameter command for Pat
//! - m_index_in: to locate the position in the defined (constant) MaterialParameter command for m
//! - M_cal_index_in: to locate the position in the defined (constant) MaterialParameter command for M
//! - cc_index_in: to locate the position in the defined (constant) MaterialParameter command for cc
//! - nb_index_in: to locate the position in the defined (constant) MaterialParameter command for nb
//! - h0_index_in: to locate the position in the defined (constant) MaterialParameter command for h0
//! - ch_index_in: to locate the position in the defined (constant) MaterialParameter command for ch
//! - G0_index_in: to locate the position in the defined (constant) MaterialParameter command for G0
//! - alpha_index_in: to locate the position in the defined (tensor internal variable) MaterialParameter command for alpha (back stress)
//! - z_index_in: to locate the position in the defined (tensor internal variable) MaterialParameter command for e0 z (fabric tensor)
//!

    DM04_alpha_Eij(int e0_index_in,
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
                   int z_index_in);

    DM04_alpha_Eij() : TensorEvolution(TE_TAG_DM04_alpha_Eij) {}; //Guanzhou added for parallel processing

    TensorEvolution* newObj();

//! rf: Dafalias-Manzari 2004
//! Hij: Evolution law for the internal variable of alpha (back stress)
//! inputs: 
//! - PlasticFlow& plastic_flow:
//! - stresstensor &Stre: the tensor of stress
//! - straintensor &Stra: the tensor of strain
//! - MaterialParameter &MaterialParameter_in: the class of material parameter

   const straintensor& Hij(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                            const straintensor& Stra, const MaterialParameter& material_parameter);

//! rf: Dafalias-Manzari 2004
//! DHij_Ds: derivation of the evolution law for the internal variable alpha (back stress) with respect to stress
//! inputs: 
//! - PlasticFlow& plastic_flow:
//! - stresstensor &Stre: the tensor of stress
//! - straintensor &Stra: the tensor of strain
//! - MaterialParameter &MaterialParameter_in: the class of material parameter

    const tensor& DHij_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                            const straintensor& Stra, const MaterialParameter& material_parameter);

//! rf: Dafalias-Manzari 2004
//! DHij_Dkin: 
//! inputs: 
//! - PlasticFlow& plastic_flow:
//! - stresstensor &Stre: the tensor of stress
//! - straintensor &Stra: the tensor of strain
//! - MaterialParameter &MaterialParameter_in: the class of material parameter

    const tensor& DHij_Dkin(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
                            const straintensor& Stra, const MaterialParameter& material_parameter);
       
    //Guanzhou added for parallel
//! sendSelf: for parallel computing
    int sendSelf(int commitTag, Channel &theChannel);  

//! recvSelf: for parallel computing
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

  private:
    
//! rf: Dafalias & Manzari 2004
//! gete0: to get e0
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double gete0(const MaterialParameter& material_parameter) const; 

//! rf: Dafalias & Manzari 2004
//! gete_r: to get e_r
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double gete_r(const MaterialParameter& material_parameter) const;      

//! rf: Dafalias & Manzari 2004
//! getlambda_c: to get lambda_c
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getlambda_c(const MaterialParameter& material_parameter) const;

//! rf: Dafalias & Manzari 2004
//! getxi: to get xi
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getxi(const MaterialParameter& material_parameter) const;

//! rf: Dafalias & Manzari 2004
//! getPat: to get pat (constant of atmosphereic pressure)
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getPat(const MaterialParameter& material_parameter) const;

//! rf: Dafalias & Manzari 2004
//! getm: to get m (constant in yield function)
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getm(const MaterialParameter& material_parameter) const;

//! rf: Dafalias & Manzari 2004
//! getM_cal: to get M_cal (critial state line slope)
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getM_cal(const MaterialParameter& material_parameter) const;

//! rf: Dafalias-Manzari 2004
//! getcc: to get cc
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
    double getcc(const MaterialParameter& material_parameter) const;

//! rf: Dafalias-Manzari 2004
//! getnb: to get nb
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
    double getnb(const MaterialParameter& material_parameter) const;

//! rf: Dafalias-Manzari 2004
//! geth0: to get h0
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
    double geth0(const MaterialParameter& material_parameter) const;

//! rf: Dafalias-Manzari 2004
//! getch: to get ch
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
    double getch(const MaterialParameter& material_parameter) const;

//! rf: Dafalias-Manzari 2004
//! getG0: to get G0
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
    double getG0(const MaterialParameter& material_parameter) const;

//! rf: Dafalias-Manzari 2004
//! getalpha: to get alpha
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter       
    const stresstensor& getalpha(const MaterialParameter& material_parameter) const;

//! rf: Dafalias-Manzari 2004
//! getz: to get z
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter       
    const stresstensor& getz(const MaterialParameter& material_parameter) const;
    
//! rf: Dafalias & Manzari 2004
//! getParameters:
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
//! - parIndex_in, which_in
    inline double getParameters(const MaterialParameter& material_parameter, int which_in) const;    

//! rf: Dafalias & Manzari 2004
//! getec: to get e_c
//! inputs:
//! - e_r, lambda_c, xi, Pat,p_c
    inline double getec(double e_r, double lambda_c, double xi, double Pat, double p_c) const;

//! rf: Dafalias & Manzari 2004
//! getg: to get g(theta)
//! inputs:
//! - c , cos3theta
    inline double getg(double c, double cos3theta) const;

    
  private:
  
    // things need to be memorized
    //int a_index;
    stresstensor alpha_in;

  private:
  
    int e0_index;
    int e_r_index;
    int lambda_c_index;
    int xi_index;
    int Pat_index;
    int m_index;
    int M_cal_index;
    int cc_index;
    int nb_index;
    int h0_index;    
    int ch_index;
    int G0_index;
    int alpha_index;
    int z_index;   
    
    static stresstensor DM04_alpha_t;

};


#endif


