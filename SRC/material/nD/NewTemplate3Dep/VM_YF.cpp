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
//                    Mahdi and Boris 27April2007 new definitions for f, df/ds...
//                    according to form 2 from lecture notes
//                    Nima Tafazzoli updated for API (Feb 2009)
///////////////////////////////////////////////////////////////////////////////
//

#ifndef VM_YF_CPP
#define VM_YF_CPP

#include "VM_YF.h"
#define YF_TAG_VM 124006

// Nima added (March 2010)
#include <iostream>
using namespace std;

stresstensor VM_YF::VMst;
stresstensor VM_YF::VMa;


//================================================================================
VM_YF::VM_YF(int k_which_in, int index_k_in, int alpha_which_in, int index_alpha_in)
: YieldFunction(YF_TAG_VM), k_which(k_which_in), index_k(index_k_in), 
  alpha_which(alpha_which_in), index_alpha(index_alpha_in)
{

}

//================================================================================
VM_YF::~VM_YF() 
{

}

//================================================================================
YieldFunction* VM_YF::newObj() 
{
 YieldFunction  *new_YF = new VM_YF(this->k_which, 
                                       this->index_k, 
                                       this->alpha_which, 
                                       this->index_alpha);

 return new_YF;
}

//================================================================================


//-out // Mahdi Taiebat 19April2007
//-out //================================================================================
//-out // f1 = 1.5*(sij-aij)*(sij-aij) - k*k = 0
//-out //================================================================================
//-out 
//-out      double VM_YF::YieldFunctionValue( const stresstensor &Stre, 
//-out                                        const MaterialParameter &MaterialParameter_in ) const
//-out      { 
//-out        // f1 = 1.5*(sij-aij)*(sij-aij) - k*k = 0
//-out        double k = getk(MaterialParameter_in);
//-out        stresstensor alpha;
//-out        if (alpha_which == 2)
//-out        {
//-out         alpha = getbackstress(MaterialParameter_in);
//-out        }
//-out        stresstensor s_bar = Stre - alpha;
//-out        s_bar = s_bar.deviator();
//-out        double s_a = (s_bar("ij")*s_bar("ij")).trace();
//-out        double YF = s_a * (3.0/2.0) - k * k;       
//-out        return YF; 
//-out      
//-out      }  
//-out      
//-out      //================================================================================
//-out      const stresstensor& VM_YF::StressDerivative(const stresstensor &Stre, 
//-out                                                  const MaterialParameter &MaterialParameter_in) const
//-out      {     
//-out        // f1 = 1.5*(sij-aij)*(sij-aij) - k*k = 0
//-out        if (alpha_which == -1) 
//-out        {
//-out          VMst = Stre.deviator() *3.0;
//-out          return VMst;
//-out        }
//-out        else if (alpha_which == 2) 
//-out        {
//-out          stresstensor s_back = getbackstress(MaterialParameter_in);
//-out          stresstensor s_bar = Stre.deviator() - s_back;
//-out          VM_YF::VMst = s_bar *3.0;
//-out          return VM_YF::VMst;
//-out        }
//-out        else 
//-out        {
//-out          opserr << "Warning!! VM_YF: Invalid Input Parameter. " << endln;
//-out          exit (1);
//-out        }
//-out      
//-out      }
//-out      
//-out      //================================================================================
//-out      double VM_YF::InScalarDerivative(const stresstensor &Stre, 
//-out                                       const MaterialParameter &MaterialParameter_in, 
//-out                                       int index_) const
//-out      {
//-out        double scalar1 = 0.0;
//-out        double k = getk(MaterialParameter_in);
//-out        
//-out        if (k_which == 1 && index_ == index_k)
//-out         // f1 = 1.5*(sij-aij)*(sij-aij) - k*k = 0
//-out         scalar1 = k * (-2.0);
//-out        return scalar1;
//-out      }
//-out      
//-out      //================================================================================
//-out      const stresstensor& VM_YF::InTensorDerivative(const stresstensor &Stre, 
//-out                                                    const MaterialParameter &MaterialParameter_in, 
//-out                                                    int index_) const
//-out      {
//-out        // this is \partial F / \partial alpha_{ij}
//-out        
//-out        // f1 = 1.5*(sij-aij)*(sij-aij) - k*k = 0
//-out        if (alpha_which == -1) 
//-out        {
//-out          VMst = Stre.deviator() * (-3.0);
//-out          return VMst;
//-out        }
//-out        else if (alpha_which == 2)  
//-out        {
//-out          stresstensor s_back = getbackstress(MaterialParameter_in);
//-out          stresstensor s_bar = Stre.deviator() - s_back;
//-out          VM_YF::VMst = s_bar * (-3.0);
//-out          return VM_YF::VMst;
//-out        }
//-out        else 
//-out        {
//-out          opserr << "Warning!! VM_YF: Invalid Input Parameter. " << endln;
//-out          exit (1);
//-out        }       
//-out      
//-out      }





//================================================================================
// f2 = Sqrt((sij-aij) * (sij-aij)) - Sqrt(2/3) * k = 0
//================================================================================

double VM_YF::YieldFunctionValue( const stresstensor &Stre, 
                                  const MaterialParameter &MaterialParameter_in ) const
{

  // f2 = [(sij-aij)*(sij-aij)]^0.5 - sqrt(2.0/3)*k
  double sqrt23 = sqrt(2.0/3.0);
  double k = getk(MaterialParameter_in);
  stresstensor alpha;
  if (alpha_which == 2)
    alpha = getbackstress(MaterialParameter_in);
  
//  alpha.print("Tinvalue" , "Tinvalue");
  
  stresstensor s_bar = Stre - alpha;
  s_bar = s_bar.deviator();
  double s_norm  = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
  double YF = s_norm - sqrt23*k;
  return YF;

}  

//================================================================================
const stresstensor& VM_YF::StressDerivative(const stresstensor &Stre, 
                                            const MaterialParameter &MaterialParameter_in) const
{

  // f2 = [(sij-aij)*(sij-aij)]^0.5 - sqrt(2.0/3)*k
  stresstensor nij;
  stresstensor alpha;
  if (alpha_which == 2)
    {
      alpha = getbackstress(MaterialParameter_in);  
    }
    
  stresstensor s_bar = Stre - alpha;
  s_bar = s_bar.deviator();
  
  double s_norm  = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
  
  if (s_norm > 0.0) 
    { 
      nij = s_bar *(1.0/s_norm);     
    }
  else 
    {
      fprintf(stderr,"ERROR stresstensor& VM_YF::StressDerivative(  s_norm <= 0 ");
      exit(1);
    }
  VM_YF::VMst.Initialize(nij);

  return VM_YF::VMst;
}

//================================================================================
double VM_YF::InScalarDerivative(const stresstensor &Stre, 
                                 const MaterialParameter &MaterialParameter_in, 
                                 int index_) const
{
  double scalar1 = 0.0;
  double k = getk(MaterialParameter_in);
  
  if (k_which == 1 && index_ == index_k)
   // f2 = [(sij-aij)*(sij-aij)]^0.5 - sqrt(2.0/3)*k
   scalar1 = -sqrt(2.0/3.0);     
  return scalar1;
}

//================================================================================
const stresstensor& VM_YF::InTensorDerivative(const stresstensor &Stre, 
                                              const MaterialParameter &MaterialParameter_in, 
                                              int index_) const
{
  // this is \partial F / \partial alpha_{ij}
  
  // f2 = [(sij-aij)*(sij-aij)]^0.5 - sqrt(2.0/3)*k
  stresstensor nij;  
  stresstensor alpha;
  if (alpha_which == 2 || index_ == index_alpha) 
    {
      alpha = getbackstress(MaterialParameter_in);
    }
//  stresstensor s_bar = Stre.deviator() - alpha;
  stresstensor s_bar = Stre - alpha;
  s_bar = s_bar.deviator();
  double s_norm = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
  if (s_norm > 0.0) 
    { 
      nij = s_bar *(-1.0/s_norm);
    }
  else 
    {
      fprintf(stderr,"ERROR  stresstensor& VM_YF::InTensorDerivative(  s_norm <= 0 ");
      exit(1);
    }

  VM_YF::VMst = nij;
  return VM_YF::VMst;
  
}



//================================================================================
int VM_YF::getNumInternalScalar() const
{
 if ( k_which == 1)
  return 1;
 else
  return 0;
}

//================================================================================
int VM_YF::getNumInternalTensor() const
{
 if (alpha_which == 2)
  return 1;
 else
  return 0;
}

//================================================================================   
int VM_YF::getYieldFunctionRank() const
{
 return 1;
}

//================================================================================   
double VM_YF::getk(const MaterialParameter &MaterialParameter_in) const
{
 // to get k
 if ( k_which == 0 && index_k <= MaterialParameter_in.getNum_Material_Constant() && index_k > 2)
  return MaterialParameter_in.getMaterial_Constant(index_k-1);
 else if( k_which == 1 && index_k <= MaterialParameter_in.getNum_Internal_Scalar() && index_k > 0)
  return MaterialParameter_in.getInternal_Scalar(index_k-1);
 else {
  opserr << "Warning!! VM_YF: Invalid Input. " << endln;
  exit (1);
 }
}

//================================================================================ 
const stresstensor& VM_YF::getbackstress(const MaterialParameter &MaterialParameter_in) const
{
 //to get back-stress
 if ( alpha_which == 2 && index_alpha <= MaterialParameter_in.getNum_Internal_Tensor() && index_alpha > 0) {
  VM_YF::VMa = MaterialParameter_in.getInternal_Tensor(index_alpha-1);
  return VM_YF::VMa;
 }
 else {
  opserr << "Warning!! VM_YF: Invalid Input. " << endln;
  exit (1);
 }
}

//================================================================================ 
//Guanzhou added for parallel
int VM_YF::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(4);
    idData.Zero();

    idData(0) = k_which;    
    idData(1) = index_k;    
    idData(2) = alpha_which;
    idData(3) = index_alpha;
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
    opserr << "VM_YF::sendSelf -- failed to send ID\n";
    return -1;
    }

    return 0;
}

//================================================================================ 
//Guanzhou added for parallel
int VM_YF::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(4);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
     opserr << "VM_YF::recvSelf -- failed to recv ID\n";
 return -1;
    }

    k_which       = idData(0);
    index_k       = idData(1);
    alpha_which = idData(2);
    index_alpha = idData(3);

    return 0;
}

#endif

