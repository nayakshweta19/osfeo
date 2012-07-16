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
//                    Mahdi and Boris 27April2007 new definitions for Q, dQ/ds...
//                    according to form 2 from lecture notes
//
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef VM_PF_CPP
#define VM_PF_CPP

#include "VM_PF.h"
#define PF_TAG_VM 120006

straintensor VM_PF::VMm;
stresstensor VM_PF::VMb;

//================================================================================
VM_PF::VM_PF(int alpha_which_in, int index_alpha_in)
: PlasticFlow(PF_TAG_VM), alpha_which(alpha_which_in), index_alpha(index_alpha_in)
{

}

//================================================================================
VM_PF::~VM_PF() 
{  

}

//================================================================================
PlasticFlow* VM_PF::newObj() 
{  
     PlasticFlow  *new_PF = new VM_PF(this->alpha_which, this->index_alpha);
     
     return new_PF;
}



// out, changing to form 2 Mahdi and Boris 27April2007  
// out, changing to form 2 Mahdi and Boris 27April2007  
// out, changing to form 2 Mahdi and Boris 27April2007  
//================================================================================
// out, Mahdi and Boris 27April2007  const straintensor& VM_PF::PlasticFlowTensor(const stresstensor &Stre, 
// out, Mahdi and Boris 27April2007                                               const straintensor &Stra, 
// out, Mahdi and Boris 27April2007                                               const MaterialParameter &MaterialParameter_in) const
// out, Mahdi and Boris 27April2007  {
// out, Mahdi and Boris 27April2007    if (alpha_which == -1) 
// out, Mahdi and Boris 27April2007      {
// out, Mahdi and Boris 27April2007        VMm = Stre.deviator() *3.0;
// out, Mahdi and Boris 27April2007        return VMm;
// out, Mahdi and Boris 27April2007      } 
// out, Mahdi and Boris 27April2007   else if (alpha_which == 2) {
// out, Mahdi and Boris 27April2007    stresstensor s_back = getalpha(MaterialParameter_in);
// out, Mahdi and Boris 27April2007    stresstensor s_bar = Stre.deviator() - s_back;
// out, Mahdi and Boris 27April2007    VM_PF::VMm = s_bar *3.0;
// out, Mahdi and Boris 27April2007    return VM_PF::VMm;
// out, Mahdi and Boris 27April2007   }
// out, Mahdi and Boris 27April2007   else {
// out, Mahdi and Boris 27April2007    opserr << "Warning!! VM_PF: Invalid Input Parameter. " << endln;
// out, Mahdi and Boris 27April2007    exit (1);
// out, Mahdi and Boris 27April2007   }
// out, Mahdi and Boris 27April2007  }
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007  ///////////////////////////////////////////////////////////////////////////////
// out, Mahdi and Boris 27April2007  const tensor& VM_PF::Dm_Ds(const stresstensor &Stre, 
// out, Mahdi and Boris 27April2007                             const straintensor &Stra, 
// out, Mahdi and Boris 27April2007                             const MaterialParameter &MaterialParameter_in) const 
// out, Mahdi and Boris 27April2007  {
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007  // This is wrong !!!
// out, Mahdi and Boris 27April2007  // Boris Jeremic and Mahdi Taiebat 07April2007
// out, Mahdi and Boris 27April2007  // see lecture notes appendix on useful formulae
// out, Mahdi and Boris 27April2007  // and also see original implementation in BJ master thesis
// out, Mahdi and Boris 27April2007  // and in ZY PhD (Template3DEP)
// out, Mahdi and Boris 27April2007  //out    tensor I2("I", 2, def_dim_2);
// out, Mahdi and Boris 27April2007  //out     tensor I_ijkl = I2("ij")*I2("kl");
// out, Mahdi and Boris 27April2007  //out     I_ijkl.null_indices();
// out, Mahdi and Boris 27April2007  //out     tensor I_ikjl = I_ijkl.transpose0110();
// out, Mahdi and Boris 27April2007  //out     tensor I_iljk = I_ijkl.transpose0111();
// out, Mahdi and Boris 27April2007  //out     tensor I4s = (I_ikjl+I_iljk)*0.5;
// out, Mahdi and Boris 27April2007  //out 
// out, Mahdi and Boris 27April2007  //out     PlasticFlow::PF_tensorR4 = I4s *3.0 - I_ijkl;
// out, Mahdi and Boris 27April2007      
// out, Mahdi and Boris 27April2007      static tensor I("I", 2, def_dim_2);
// out, Mahdi and Boris 27April2007      static tensor temp1 = I("im") * I("jn");
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007      static tensor I2("I", 2, def_dim_2);
// out, Mahdi and Boris 27April2007      static tensor temp2 = I2("mn") * I2("ij") * (1.0/3.0);
// out, Mahdi and Boris 27April2007      
// out, Mahdi and Boris 27April2007      PlasticFlow::PF_tensorR4 = (temp1 - temp2) * 3.0;
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007      return PlasticFlow::PF_tensorR4;
// out, Mahdi and Boris 27April2007  }
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007  ///////////////////////////////////////////////////////////////////////////////
// out, Mahdi and Boris 27April2007  const tensor& VM_PF::Dm_Diso(const stresstensor &Stre, 
// out, Mahdi and Boris 27April2007                               const straintensor &Stra, 
// out, Mahdi and Boris 27April2007                               const MaterialParameter &MaterialParameter_in) const 
// out, Mahdi and Boris 27April2007  {
// out, Mahdi and Boris 27April2007      tensor PF_R2(2, def_dim_2, 0.0);
// out, Mahdi and Boris 27April2007      PlasticFlow::PF_tensorR2.Initialize(PF_R2);
// out, Mahdi and Boris 27April2007      
// out, Mahdi and Boris 27April2007      return PlasticFlow::PF_tensorR2;
// out, Mahdi and Boris 27April2007  }
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007  ///////////////////////////////////////////////////////////////////////////////
// out, Mahdi and Boris 27April2007  const tensor& VM_PF::Dm_Dkin(const stresstensor &Stre, 
// out, Mahdi and Boris 27April2007                               const straintensor &Stra, 
// out, Mahdi and Boris 27April2007                               const MaterialParameter &MaterialParameter_in) const 
// out, Mahdi and Boris 27April2007  {
// out, Mahdi and Boris 27April2007  //// this is also wrong 
// out, Mahdi and Boris 27April2007  // Boris Jeremic and Mahdi Taiebat  07April2007
// out, Mahdi and Boris 27April2007  //out
// out, Mahdi and Boris 27April2007  //out  
// out, Mahdi and Boris 27April2007  //out  
// out, Mahdi and Boris 27April2007  //out      tensor I2("I", 2, def_dim_2);
// out, Mahdi and Boris 27April2007  //out      tensor I_ijkl = I2("ij")*I2("kl");
// out, Mahdi and Boris 27April2007  //out      I_ijkl.null_indices();
// out, Mahdi and Boris 27April2007  //out      tensor I_ikjl = I_ijkl.transpose0110();
// out, Mahdi and Boris 27April2007  //out      tensor I_iljk = I_ijkl.transpose0111();
// out, Mahdi and Boris 27April2007  //out      tensor I4s = (I_ikjl+I_iljk)*0.5;
// out, Mahdi and Boris 27April2007  // similar to upper case for function
// out, Mahdi and Boris 27April2007  // const tensor& VM_PF::Dm_Ds
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007      static tensor I("I", 2, def_dim_2);
// out, Mahdi and Boris 27April2007      static tensor temp1 = I("im") * I("jn");
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007      static tensor I2("I", 2, def_dim_2);
// out, Mahdi and Boris 27April2007      static tensor temp2 = I2("mn") * I2("ij") * (1.0/3.0);
// out, Mahdi and Boris 27April2007      
// out, Mahdi and Boris 27April2007      PlasticFlow::PF_tensorR4 = (temp1 - temp2) * -3.0;
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007      return PlasticFlow::PF_tensorR4;
// out, Mahdi and Boris 27April2007  }
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007  //================================================================================ 
// out, Mahdi and Boris 27April2007  const stresstensor& VM_PF::getalpha(const MaterialParameter &MaterialParameter_in) const
// out, Mahdi and Boris 27April2007  {
// out, Mahdi and Boris 27April2007   //to get alpha (backstress)
// out, Mahdi and Boris 27April2007   if ( alpha_which == 2 && index_alpha <= MaterialParameter_in.getNum_Internal_Tensor() && index_alpha > 0) {
// out, Mahdi and Boris 27April2007    VM_PF::VMb = MaterialParameter_in.getInternal_Tensor(index_alpha-1);
// out, Mahdi and Boris 27April2007    return VM_PF::VMb;
// out, Mahdi and Boris 27April2007   }
// out, Mahdi and Boris 27April2007   
// out, Mahdi and Boris 27April2007   opserr << "Warning!! VM_PF: Invalid Input. " << endln;
// out, Mahdi and Boris 27April2007   exit (1);
// out, Mahdi and Boris 27April2007  }
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007  
// out, Mahdi and Boris 27April2007  
// out, changing to form 2 Mahdi and Boris 27April2007







//================================================================================
const straintensor& VM_PF::PlasticFlowTensor(const stresstensor &Stre, 
                                             const straintensor &Stra, 
                                             const MaterialParameter &MaterialParameter_in) const
{

       // Q2 = [(sij-aij)*(sij-aij)]^0.5 - sqrt(2.0/3)*k
       stresstensor alpha;
       stresstensor nij;
       if (alpha_which == 2)
         {
           alpha = getalpha(MaterialParameter_in);  
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
           fprintf(stderr,"ERROR stresstensor& VM_PF::StressDerivative(  s_norm <= 0 ");
           exit(1);
         }
       VM_PF::VMm.Initialize(nij);  
       return VM_PF::VMm;
     
     }
     


///////////////////////////////////////////////////////////////////////////////
const tensor& VM_PF::Dm_Ds(const stresstensor &Stre, 
                           const straintensor &Stra, 
                           const MaterialParameter &MaterialParameter_in) const 
{

// This is wrong !!!
// Boris Jeremic and Mahdi Taiebat 07April2007
// see lecture notes appendix on useful formulae
// and also see original implementation in BJ master thesis
// and in ZY PhD (Template3DEP)
//out    tensor I2("I", 2, def_dim_2);
//out     tensor I_ijkl = I2("ij")*I2("kl");
//out     I_ijkl.null_indices();
//out     tensor I_ikjl = I_ijkl.transpose0110();
//out     tensor I_iljk = I_ijkl.transpose0111();
//out     tensor I4s = (I_ikjl+I_iljk)*0.5;
//out 
//out     PlasticFlow::PF_tensorR4 = I4s *3.0 - I_ijkl;
    


    tensor I2("I", 2, def_dim_2);
    tensor I_ijmn = I2("ij")*I2("mn");
    I_ijmn.null_indices();
    tensor I_imjn = I_ijmn.transpose_1234_1324();
    tensor I4 = I_imjn - I_ijmn * 0.33333333333;


    stresstensor alpha;
    if (alpha_which == 2)
      {
        alpha = getalpha(MaterialParameter_in);  
      }
    stresstensor s_bar = Stre - alpha;
    
    s_bar = s_bar.deviator();
    double s_norm_square  = (s_bar("ij")*s_bar("ij")).trace();
    if (s_norm_square <= 0.0) 
      { 
        fprintf(stderr,"ERROR stresstensor& VM_PF::StressDerivative(  s_norm <= 0 ");
        exit(1);
      }
    double s_norm = sqrt(s_norm_square);
    double s_norm_1_5 = sqrt(s_norm_square*s_norm_square*s_norm_square);

    

    tensor temp1  = s_bar("ij") * s_bar("mn") * (1.0/s_norm_1_5);

    tensor temp2 = I4 * (1.0/s_norm);

    PlasticFlow::PF_tensorR4 = temp2 - temp1;

    return PlasticFlow::PF_tensorR4;

}

///////////////////////////////////////////////////////////////////////////////
const tensor& VM_PF::Dm_Diso(const stresstensor &Stre, 
                             const straintensor &Stra, 
                             const MaterialParameter &MaterialParameter_in) const 
{
    tensor PF_R2(2, def_dim_2, 0.0);
    PlasticFlow::PF_tensorR2.Initialize(PF_R2);
    
    return PlasticFlow::PF_tensorR2;
}

///////////////////////////////////////////////////////////////////////////////
const tensor& VM_PF::Dm_Dkin(const stresstensor &Stre, 
                             const straintensor &Stra, 
                             const MaterialParameter &MaterialParameter_in) const 
{
//// this is also wrong 
// Boris Jeremic and Mahdi Taiebat  07April2007
//out
//out  
//out  
//out      tensor I2("I", 2, def_dim_2);
//out      tensor I_ijkl = I2("ij")*I2("kl");
//out      I_ijkl.null_indices();
//out      tensor I_ikjl = I_ijkl.transpose0110();
//out      tensor I_iljk = I_ijkl.transpose0111();
//out      tensor I4s = (I_ikjl+I_iljk)*0.5;
// similar to upper case for function
// const tensor& VM_PF::Dm_Ds

       stresstensor alpha;
       if (alpha_which == 2)
         {
           alpha = getalpha(MaterialParameter_in);  
         }
       stresstensor s_bar = Stre - alpha;
       
       s_bar = s_bar.deviator();
       double s_norm_square  = (s_bar("ij")*s_bar("ij")).trace();
       if (s_norm_square <= 0.0) 
         { 
           fprintf(stderr,"ERROR stresstensor& VM_PF::StressDerivative(  s_norm <= 0 ");
           exit(1);
         }
       double s_norm = sqrt(s_norm_square);
       double s_norm_1_5 = sqrt(s_norm_square*s_norm_square*s_norm_square);

       tensor temp1  = s_bar("ij") * s_bar("mn") * (1.0/s_norm_1_5);
       tensor delta("I", 2, def_dim_2);
       tensor temp2 = delta("im") * delta("jn");
       tensor temp3 = temp2.transpose_1234_1324();
       temp3 = temp3 * (1.0/s_norm);
       
       PlasticFlow::PF_tensorR4 = temp1  -  temp3;
       return PlasticFlow::PF_tensorR4;

}

//================================================================================ 
const stresstensor& VM_PF::getalpha(const MaterialParameter &MaterialParameter_in) const
{
 //to get alpha (backstress)
 if ( alpha_which == 2 && index_alpha <= MaterialParameter_in.getNum_Internal_Tensor() && index_alpha > 0) 
  {
    VM_PF::VMb = MaterialParameter_in.getInternal_Tensor(index_alpha-1);
    return VM_PF::VMb;
  }
 
 opserr << "Warning!! VM_PF: Invalid Input. " << endln;
 exit (1);
}


//Guanzhou added for parallel
int VM_PF::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(2);
    idData.Zero();

    idData(0) = alpha_which;
    idData(1) = index_alpha;
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
    opserr << "VM_PF::sendSelf -- failed to send ID\n";
    return -1;
    }

    return 0;
}

//Guanzhou added for parallel
int VM_PF::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(2);
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
     opserr << "VM_PF::recvSelf -- failed to recv ID\n";
 return -1;
    }

    alpha_which  = idData(0);
    index_alpha  = idData(1);

    return 0;
}


#endif

