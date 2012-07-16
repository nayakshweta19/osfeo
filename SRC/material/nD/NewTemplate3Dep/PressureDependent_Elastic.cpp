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
// UPDATE HISTORY:    Guanzhou Jie updated for parallel Dec 2006
//
///////////////////////////////////////////////////////////////////////////////
//

// This is based on E = E0*(p/p_ref)^m
// 
// Constants:
// 1: E0:    Young's modulus when p at p_ref;
// 2: v:     Poissoin's ratio;
// 3: m:     power factor, usually = 0.5;
// 4: p_ref  rerefence pressure;
// 5: k_cut  cut-off factor, when E = k_c *E0, E = k_c*E0;

#ifndef PressureDependent_Elastic_CPP
#define PressureDependent_Elastic_CPP

#include "PressureDependent_Elastic.h"

BJtensor PressureDependent_Elastic::ElasticStiffness(4, def_dim_4, 0.0);

PressureDependent_Elastic::PressureDependent_Elastic(int E0_in, 
                                                     int v_in,
                                                     int m_in,
                                                     int p_ref_in,
                                                     int k_cut_in,
                                                     const stresstensor& initialStress,
                                                     const straintensor& initialStrain)
: ElasticState(ES_TAG_PressureDependentElastic, initialStress, initialStrain),
  E0_index(E0_in),
  v_index(v_in),
  m_index(m_in),
  p_ref_index(p_ref_in),
  k_cut_index(k_cut_in) 
{

}

// Create a new 
ElasticState* PressureDependent_Elastic::newObj() 
{
	ElasticState * Els = new  PressureDependent_Elastic(this->E0_index, 
	                                                    this->v_index,
	                                                    this->m_index,
	                                                    this->p_ref_index,
	                                                    this->k_cut_index,
	                                                    this->Stress,
	                                                    this->Strain);
     return Els;
}

// Get Stiffness Tensor
const BJtensor& PressureDependent_Elastic::getElasticStiffness (const MaterialParameter &MaterialParameter_in) const
{
    // Kronecker delta tensor
    BJtensor I2("I", 2, def_dim_2);

    BJtensor I_ijkl = I2("ij")*I2("kl");
    I_ijkl.null_indices();
    BJtensor I_ikjl = I_ijkl.transpose0110();
    BJtensor I_iljk = I_ijkl.transpose0111();
    BJtensor I4s = (I_ikjl+I_iljk)*0.5;
    
    double E0 = getE0(MaterialParameter_in);
    double v = getv(MaterialParameter_in);
    if (v >= 0.5 || v < -1.0) {
	    opserr << "Warning!! PressureDependent_Elastic: Invalid possoin's ratio. " << endln;
	    exit (1);
    }
    double m = getm(MaterialParameter_in);
    double p_ref = getp_ref(MaterialParameter_in);
    if ( p_ref <= 0.0) {
	    opserr << "Warning!! PressureDependent_Elastic: Invalid reference pressure. " << endln;
	    exit (1);
    }   
    double k_cut = getk_cut(MaterialParameter_in);
        
    double p = this->getStress().p_hydrostatic();
    
    double E_cut = E0 * k_cut;
    double E_cal = E0 * pow((p/p_ref), m);
    double E = (E_cal > E_cut) ? E_cal : E_cut;
    
    double K = E / ( 3.0*(1.0-2.0*v) ) ;
    double G = E / ( 2.0*(1.0+v) );
       
    // Building elasticity tensor
    PressureDependent_Elastic::ElasticStiffness = I_ijkl *(K - 2.0*G/3.0) + I4s *(2.0*G);

    return PressureDependent_Elastic::ElasticStiffness;
}

// Get E0
double PressureDependent_Elastic::getE0(const MaterialParameter &MaterialParameter_in) const
{
	if ( E0_index > MaterialParameter_in.getNum_Material_Constant() || E0_index < 2) { 
		opserr << "PressureDependent_Elastic: Invalid Input. " << endln;
		exit (1);
	}
	else
		return MaterialParameter_in.getMaterial_Constant(E0_index - 1);
}

// Get v
double PressureDependent_Elastic::getv(const MaterialParameter &MaterialParameter_in) const
{
	if ( v_index > MaterialParameter_in.getNum_Material_Constant() || v_index < 2) { 
		opserr << "PressureDependent_Elastic: Invalid Input. " << endln;
		exit (1);
	}
	else
		return MaterialParameter_in.getMaterial_Constant(v_index - 1); 
}

// Get m
double PressureDependent_Elastic::getm(const MaterialParameter &MaterialParameter_in) const
{
	if ( m_index > MaterialParameter_in.getNum_Material_Constant() || m_index < 2) { 
		opserr << "PressureDependent_Elastic: Invalid Input. " << endln;
		exit (1);
	}
	else
		return MaterialParameter_in.getMaterial_Constant(m_index - 1); 
}

// Get p_ref
double PressureDependent_Elastic::getp_ref(const MaterialParameter &MaterialParameter_in) const
{
	if ( p_ref_index > MaterialParameter_in.getNum_Material_Constant() || p_ref_index < 2) { 
		opserr << "PressureDependent_Elastic: Invalid Input. " << endln;
		exit (1);
	}
	else
		return MaterialParameter_in.getMaterial_Constant(p_ref_index - 1); 
}

// Get p_cut
double PressureDependent_Elastic::getk_cut(const MaterialParameter &MaterialParameter_in) const
{
	if ( k_cut_index > MaterialParameter_in.getNum_Material_Constant() || k_cut_index < 2) { 
		opserr << "PressureDependent_Elastic: Invalid Input. " << endln;
		exit (1);
	}
	else
		return MaterialParameter_in.getMaterial_Constant(k_cut_index - 1); 
}

//Guanzhou added for parallel
int PressureDependent_Elastic::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(5);
    idData.Zero();

    idData(0) = E0_index;   
    idData(1) = v_index;    
    idData(2) = m_index;    
    idData(3) = p_ref_index;
    idData(4) =	k_cut_index;
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "PressureDependent_Elastic::sendSelf -- failed to send ID\n";
   	return -1;
    }

    if (theChannel.sendnDarray(dataTag, commitTag, this->Stress) < 0) {
    	opserr << "PressureDependent_Elastic::sendSelf() -  failed to send nDarray Stress\n";
    	return -1;
    }
   
    if (theChannel.sendnDarray(dataTag, commitTag, this->Strain) < 0) {
    	opserr << "PressureDependent_Elastic::sendSelf() -  failed to send nDarray Strain\n";
    	return -1;
    }
    
    return 0;
}

//Guanzhou added for parallel
int PressureDependent_Elastic::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(5);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "PressureDependent_Elastic::recvSelf -- failed to recv ID\n";
	return -1;
    }

    E0_index    = idData(0);
    v_index     = idData(1);
    m_index = idData(2);
    p_ref_index = idData(3);
    k_cut_index = idData(4);
     
    if (theChannel.recvnDarray(dataTag, commitTag, this->Stress) < 0) {
    	opserr << "PressureDependent_Elastic::recvSelf() -  failed to recv nDarray Stress\n";
    	return -1;
    }

    if (theChannel.recvnDarray(dataTag, commitTag, this->Strain) < 0) {
    	opserr << "PressureDependent_Elastic::recvSelf() -  failed to recv nDarray Strain\n";
    	return -1;
    }

    return 0;
}


#endif

