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

#ifndef SANISAND_PF_H
#define SANISAND_PF_H

#include "PlasticFlow.h"
#include <math.h>
#include <fstream>
#include <iostream>
#define PF_TAG_SANISAND 120005

class SANISAND_PF : public PlasticFlow
{
  public:   

    SANISAND_PF() : PlasticFlow(PF_TAG_SANISAND) {}; //Guanzhou added for parallel processing


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
//!  8- A0:        dilatancy parameter;
//!  9- nd         dilatancy parameter;
//! 10- m:         opening of the yield surface;
//! 11- alpha:     back stress ratio tensor in yield function; (the 1st tensorial internal variable);
//! 12- z:         fabric dilation internal tensor (the 2nd tensorial internal variable); 
//! 13- X:         LCC parameter; 
//! 
//! GENERAL NOTE
//! 
//! material parameter type: 
//! 0 for the (constant) material parameter, 
//! 1 for the (initial) internal scalar, 
//! 2 for the (initial) internal tensor
//! 
//! material parameter id: 
//! to locate the position in the deﬁned material parameter command 
//! 
//! Plastic Flow function for SANISAND:
//! inputs:
//! - e0_which:              type of parameter e0       --> (0)
//! - index_e0:              id of parameter   e0
//! - e_r_which_in:          type of parameter e_r      --> (0) 
//! - index_e_r_in:          id of parameter   e_r
//! - lambda_c_which_in:     type of parameter lambda   --> (0)
//! - index_lambda_c_in:     id of parameter   lambda
//! - xi_which_in:           type of parameter xi       --> (0)
//! - index_xi_in:           id of parameter   xi       
//! - Pat_which_in:          type of parameter Pat      --> (0)
//! - index_Pat_in:          id of parameter   Pat      
//! - alpha_cc_which_in:     type of parameter alpha_cc --> (0)
//! - index_alpha_cc_in:     id of parameter   alpha_cc
//! - c_which_in:            type of parameter c        --> (0)
//! - index_c_in:            id of parameter   c        
//! - A0_which_in:           type of parameter A0       --> (0)
//! - index_A0_in:           id of parameter   A0       
//! - nd_which_in:           type of parameter nd       --> (0)
//! - index_nd_in:           id of parameter   nd
//! - alpha_which_in:        type of parameter alpha    --> (2)
//! - index_alpha_in:        id of parameter   alpha    
//! - m_which_in:            type of parameter m        --> (0)
//! - index_m_in:            id of parameter   m        
//! - z_which_in:            type of parameter z        --> (2)
//! - index_z_in:            id of parameter   z
//! - X_which_in:            type of parameter X        --> (0)
//! - index_X_in:            id of parameter   X


    SANISAND_PF(int e0_which,          int index_e0,
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
                int X_which_in,        int index_X_in);

//! Deconstructor for SANISAND plastic flow function	
    ~SANISAND_PF();     	    

    PlasticFlow* newObj();

    const char *getPlasticFlowType(void) const {return "SANISAND_PF";}; 


//! PlasticFlowTensor: to define the plastic flow
//! inputs: 
//! - stresstensor &Stre: the tensor of stress
//! - straintensor &Stra: the tensor of strain
//! - MaterialParameter &MaterialParameter_in: the class of material parameter  

    const straintensor& PlasticFlowTensor(const stresstensor &Stre, 
                                          const straintensor &Stra, 
                                          const MaterialParameter &MaterialParameter_in) const;

//out-     const tensor& Dm_Ds(const stresstensor &Stre, 
//out-                         const straintensor &Stra, 
//out-                         const MaterialParameter &MaterialParameter_in) const;
//out- 
//out-     const tensor& Dm_Diso(const stresstensor &Stre, 
//out-                           const straintensor &Stra, 
//out-                           const MaterialParameter &MaterialParameter_in) const;
//out- 
//out-     const tensor& Dm_Dkin(const stresstensor &Stre, 
//out-                           const straintensor &Stra, 
//out-                           const MaterialParameter &MaterialParameter_in) const;
//out- 
//out-     const tensor& Dm_Dkin2(const stresstensor &Stre, 
//out-                            const straintensor &Stra, 
//out-                            const MaterialParameter &MaterialParameter_in) const;

    //Guanzhou added for parallel
//! sendSelf: for parallel computing
    int sendSelf(int commitTag, Channel &theChannel);  
//! recvSelf: for parallel computing
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    


  private: 

//! Reference: Taiebat & Dafalias (2008)
//! gete0: to get e0
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double gete0(const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! gete_r: to get e_r
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double gete_r(const MaterialParameter &MaterialParameter_in) const;      

//! Reference: Taiebat & Dafalias (2008)
//! getlambda: to get lambda
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getlambda(const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! getxi: to get xi
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getxi(const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! getPat: to get pat (constant of atmosphereic pressure)
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getPat(const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! getalpha_cc: to get alpha_cc
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getalpha_cc(const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! getc: to get c
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getc(const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! getA0: to get A0
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getA0(const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! getnd: to get nd
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getnd(const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! getm: to get m (constant in yield function)
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getm(const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! getalpha: to get alpha
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
    const stresstensor& getalpha(const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! getz: to get z
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter
    const stresstensor& getz(const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! getx: to get x
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
    double getX(const MaterialParameter &MaterialParameter_in) const;

//! Reference: Taiebat & Dafalias (2008)
//! getParameters:
//! inputs:
//! - MaterialParameter &MaterialParameter_in: the class of material parameter 
//! - parIndex_in:
//! - which_in:
    inline double getParameters(const MaterialParameter &MaterialParameter_in, int parIndex_in, int which_in) const;    

//! Reference: Taiebat & Dafalias (2008)
//! getec: to get 
//! inputs:
//! - e_r, lambda_c, xi, Pat, p_c

    inline double getec(double e_r, double lambda, double xi, double Pat, double p_c) const;

//! Reference: Taiebat & Dafalias (2008)
//! getg: to get 
//! inputs:
//! - c, cos3theta
    inline double getg(double c, double cos3theta) const;

  private:
    
    int e0_which;         int index_e0;
    int e_r_which;        int index_e_r;
    int lambda_which;     int index_lambda;    
    int xi_which;         int index_xi;
    int Pat_which;        int index_Pat; 
    int alpha_cc_which;   int index_alpha_cc; 
    int c_which;          int index_c;
    int A0_which;         int index_A0;   
    int nd_which;         int index_nd;
    int m_which;          int index_m;
    int alpha_which;      int index_alpha;
    int z_which;          int index_z;
    int X_which;          int index_X;
            
    static straintensor SANISANDm;
    static stresstensor SANISANDtemp;
};


#endif

