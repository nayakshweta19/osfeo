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

#ifndef DM04_z_Eij_H
#define DM04_z_Eij_H 

#include "TensorEvolution.h"
#define TE_TAG_DM04_z_Eij 121004

class DM04_z_Eij : public TensorEvolution
{
  public:
  
    DM04_z_Eij() : TensorEvolution(TE_TAG_DM04_z_Eij) {}; //Guanzhou added for parallel processing

//!
//! rf: Dafalias-Manzari 2004
//! inputs:
//! - index_e0: to locate the position in the defined (constant) MaterialParameter command for e0
//! - index_e_r_in: to locate the position in the defined (constant) MaterialParameter command for e_r
//! - index_lambda_c_in: to locate the position in the defined (constant) MaterialParameter command for c
//! - index_xi_in: to locate the position in the defined (constant) MaterialParameter command for xi
//! - index_Pat_in: to locate the position in the defined (constant) MaterialParameter command for Pat
//! - index_m_in: to locate the position in the defined (constant) MaterialParameter command for m
//! - index_M_cal_in: to locate the position in the defined (constant) MaterialParameter command for M
//! - index_cc_in: to locate the position in the defined (constant) MaterialParameter command for cc
//! - index_A0_in: to locate the position in the defined (constant) MaterialParameter command for A0
//! - index_nd_in: to locate the position in the defined (constant) MaterialParameter command for nd
//! - alpha_index_in: to locate the position in the defined (tensor internal variable) MaterialParameter command for alpha (back stress)
//! - z_index_in: to locate the position in the defined (tensor internal variable) MaterialParameter command for e0 z (fabric tensor)
//!

    DM04_z_Eij(int index_e0,
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
	    int z_index_in);

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

    const tensor& DHij_Dkin2(const PlasticFlow& plastic_flow, const stresstensor& Stre, 
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
    double gete0(const MaterialParameter &MaterialParameter_in) const;

//! rf: Dafalias & Manzari 2004
//! gete_r: to get e_r
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double gete_r(const MaterialParameter &MaterialParameter_in) const;      

//! rf: Dafalias & Manzari 2004
//! getlambda_c: to get lambda_c
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getlambda_c(const MaterialParameter &MaterialParameter_in) const;

//! rf: Dafalias & Manzari 2004
//! getxi: to get xi
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getxi(const MaterialParameter &MaterialParameter_in) const;

//! rf: Dafalias & Manzari 2004
//! getPat: to get pat (constant of atmosphereic pressure)
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getPat(const MaterialParameter &MaterialParameter_in) const;

//! rf: Dafalias & Manzari 2004
//! getm: to get m (constant in yield function)
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getm(const MaterialParameter &MaterialParameter_in) const;

//! rf: Dafalias & Manzari 2004
//! getM_cal: to get M_cal (critial state line slope)
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getM_cal(const MaterialParameter &MaterialParameter_in) const;

//! rf: Dafalias-Manzari 2004
//! getcc: to get cc
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
    double getcc(const MaterialParameter &MaterialParameter_in) const;

//! rf: Dafalias & Manzari 2004
//! getA0: to get A0
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getA0(const MaterialParameter &MaterialParameter_in) const;

//! rf: Dafalias & Manzari 2004
//! getnd: to get nd
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getnd(const MaterialParameter &MaterialParameter_in) const;

    double getc_z(const MaterialParameter& material_parameter) const;
    double getz_max(const MaterialParameter& material_parameter) const;

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

    inline double getParameters(const MaterialParameter &MaterialParameter_in, int which_in) const;    
    inline double getec(double e_r, double lambda_c, double xi, double Pat, double p_c) const;
    inline double getg(double c, double cos3theta) const;

  private:
    int index_e0;
    int index_e_r;
    int index_lambda_c;
    int index_xi;
    int index_Pat;
    int index_m;
    int index_M_cal;
    int index_cc;
    int index_A0;
    int index_nd; 
    int c_z_index;
    int z_max_index;
    int alpha_index;
    int z_index;

    static stresstensor DM04_z_t;
};


#endif

