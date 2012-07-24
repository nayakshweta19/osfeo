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
// DATE:              Fall 2005
// UPDATE HISTORY:    06/2006, add functions for matrix based elements, CZ
//                    10/2006, add various more algorithms, CZ
//                    Guanzhou Jie updated for parallel Dec 2006
//                    Mahdi Taiebat & Boris Jeremic debugged the code for the
//                    implicit method, April2007
//                    Feb2008 adding ScaledExplicit, Mahdi's idea
//                    Nima Tafazzoli updated for API (Feb 2009)
/////////////////////////////////////////////////////////////////////////////


#ifndef NewTemplate3Dep_CPP
#define NewTemplate3Dep_CPP

#include "NewTemplate3Dep.h"


const  straintensor NewTemplate3Dep::ZeroStrain;
const  stresstensor NewTemplate3Dep::ZeroStress;
const  BJtensor NewTemplate3Dep::ZeroI4(4, def_dim_4, 0.0);
const int NewTemplate3Dep::ISMAX = 32;
const int NewTemplate3Dep::ITMAX = 40;
const double NewTemplate3Dep::TOL = 1.0e-7;
const double NewTemplate3Dep::FTOL = 1.0e-7;

// For Matrix based elements
Matrix NewTemplate3Dep::D(6,6);
Vector NewTemplate3Dep::sigma(6);
Vector NewTemplate3Dep::epsilon(6);

BJmatrix NewTemplate3Dep::TangentMatrix( 6, 6, 0.0);

// Constructor
//================================================================================
NewTemplate3Dep::NewTemplate3Dep( int tag,
                                  MaterialParameter *ptr_material_parameter_in,
                                  ElasticState      *ptr_elastic_state_in,
                                  YieldFunction     *ptr_yield_function_in ,
                                  PlasticFlow       *ptr_plastic_flow_in,
                                  ScalarEvolution  **ptr_scalar_evolution_in,
                                  TensorEvolution  **ptr_tensor_evolution_in,
                                  int caseIndex_in, int subStep_in)
:NDMaterial(tag, ND_TAG_NewTemplate3Dep), caseIndex(caseIndex_in), subStep(subStep_in), divergeOrnot(0)
{
    if ( ptr_material_parameter_in )
      ptr_material_parameter = ptr_material_parameter_in->newObj();
    else {
      cout << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the input parameters. " << endln;
      exit(1);
    }

    if ( ptr_elastic_state_in )
      ptr_elastic_state = ptr_elastic_state_in->newObj();
    else{
      cout << "NewTemplate3Dep:: NewTemplate3Dep failed to get copy of elastic material. " << endln;
      exit(1);
    }

    if ( ptr_yield_function_in )
      ptr_yield_function = ptr_yield_function_in->newObj();
    else {
      cout << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the yield function. " << endln;
      exit(1);
    }

    if ( ptr_plastic_flow_in )
      ptr_plastic_flow = ptr_plastic_flow_in->newObj();
    else {
      cout << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the plastic flow. " << endln;
      exit(1);
    }

    // Scalar (isotropic) Evolution
    if ( ptr_material_parameter_in->getNum_Internal_Scalar() > 0 ) {
      ptr_scalar_evolution = new ScalarEvolution* [ptr_material_parameter_in->getNum_Internal_Scalar()];
      for (int i = 0; i < ptr_material_parameter_in->getNum_Internal_Scalar(); i++)
        ptr_scalar_evolution[i] = (ptr_scalar_evolution_in[i])->newObj();
    }
    else
      ptr_scalar_evolution = NULL;

    // Tensor (kinematic) Evolution
    if ( ptr_material_parameter_in->getNum_Internal_Tensor() > 0 ) {
      ptr_tensor_evolution = new TensorEvolution* [ptr_material_parameter_in->getNum_Internal_Tensor()];
      for (int i = 0; i < ptr_material_parameter_in->getNum_Internal_Tensor(); i++)
        ptr_tensor_evolution[i] = (ptr_tensor_evolution_in[i])->newObj();
    }
    else
      ptr_tensor_evolution = NULL;

    int err;
    err = this->revertToStart();
}

// Constructor
//================================================================================
NewTemplate3Dep::NewTemplate3Dep( int tag,
                                  MaterialParameter *ptr_material_parameter_in,
                                  ElasticState      *ptr_elastic_state_in,
                                  YieldFunction     *ptr_yield_function_in ,
                                  PlasticFlow       *ptr_plastic_flow_in,
                                  ScalarEvolution  **ptr_scalar_evolution_in,
                                  int caseIndex_in, int subStep_in)
:NDMaterial(tag, ND_TAG_NewTemplate3Dep), caseIndex(caseIndex_in), subStep(subStep_in), divergeOrnot(0)
{
    if ( ptr_material_parameter_in )
      ptr_material_parameter = ptr_material_parameter_in->newObj();
    else {
      cout << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the input parameters. " << endln;
      exit(1);
    }

    if ( ptr_elastic_state_in )
      ptr_elastic_state = ptr_elastic_state_in->newObj();
    else{
      cout << "NewTemplate3Dep:: NewTemplate3Dep failed to get copy of elastic material. " << endln;
      exit(1);
    }

    if ( ptr_yield_function_in )
      ptr_yield_function = ptr_yield_function_in->newObj();
    else {
      cout << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the yield function. " << endln;
      exit(1);
    }

    if ( ptr_plastic_flow_in )
      ptr_plastic_flow = ptr_plastic_flow_in->newObj();
    else {
      cout << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the plastic flow. " << endln;
      exit(1);
    }

    // Scalar (isotropic) Evolution
    if ( ptr_material_parameter_in->getNum_Internal_Scalar() > 0 ) {
      ptr_scalar_evolution = new ScalarEvolution* [ptr_material_parameter_in->getNum_Internal_Scalar()];
      for (int i = 0; i < ptr_material_parameter_in->getNum_Internal_Scalar(); i++)
        ptr_scalar_evolution[i] = (ptr_scalar_evolution_in[i])->newObj();
    }
    else
      ptr_scalar_evolution = NULL;

    // Tensor (kinematic) Evolution
    ptr_tensor_evolution = NULL;

    int err;
    err = this->revertToStart();
}

// Constructor
//================================================================================
NewTemplate3Dep::NewTemplate3Dep( int tag,
                                  MaterialParameter *ptr_material_parameter_in,
                                  ElasticState      *ptr_elastic_state_in,
                                  YieldFunction     *ptr_yield_function_in ,
                                  PlasticFlow       *ptr_plastic_flow_in,
                                  TensorEvolution **ptr_tensor_evolution_in,
                                  int caseIndex_in, int subStep_in)
:NDMaterial(tag, ND_TAG_NewTemplate3Dep), caseIndex(caseIndex_in), subStep(subStep_in), divergeOrnot(0)
{
    if ( ptr_material_parameter_in )
      ptr_material_parameter = ptr_material_parameter_in->newObj();
    else {
      cout << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the input parameters. " << endln;
      exit(1);
    }

    if ( ptr_elastic_state_in )
      ptr_elastic_state = ptr_elastic_state_in->newObj();
    else{
      cout << "NewTemplate3Dep:: NewTemplate3Dep failed to get copy of elastic material. " << endln;
      exit(1);
    }

    if ( ptr_yield_function_in )
      ptr_yield_function = ptr_yield_function_in->newObj();
    else {
      cout << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the yield function. " << endln;
      exit(1);
    }

    if ( ptr_plastic_flow_in )
       ptr_plastic_flow = ptr_plastic_flow_in->newObj();
    else {
      cout << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the plastic flow. " << endln;
      exit(1);
    }

    // Scalar (isotropic) Evolution
    ptr_scalar_evolution = NULL;

    // Tensor (kinematic) Evolution
    if ( ptr_material_parameter_in->getNum_Internal_Tensor() > 0 ) {
      ptr_tensor_evolution = new TensorEvolution* [ptr_material_parameter_in->getNum_Internal_Tensor()];
      for (int i = 0; i < ptr_material_parameter_in->getNum_Internal_Tensor(); i++)
        ptr_tensor_evolution[i] = ptr_tensor_evolution_in[i]->newObj();
    }
    else
      ptr_tensor_evolution = NULL;

    int err;
    err = this->revertToStart();
}

// Constructor
//================================================================================
NewTemplate3Dep::NewTemplate3Dep( int tag,
                                  MaterialParameter *ptr_material_parameter_in,
                                  ElasticState      *ptr_elastic_state_in,
                                  YieldFunction     *ptr_yield_function_in ,
                                  PlasticFlow       *ptr_plastic_flow_in,
                                  int caseIndex_in, int subStep_in)
:NDMaterial(tag, ND_TAG_NewTemplate3Dep), caseIndex(caseIndex_in), subStep(subStep_in), divergeOrnot(0)
{
    if ( ptr_material_parameter_in )
      ptr_material_parameter = ptr_material_parameter_in->newObj();
    else {
      cout << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the input parameters. " << endln;
      exit(1);
    }

    if ( ptr_elastic_state_in )
      ptr_elastic_state = ptr_elastic_state_in->newObj();
    else{
      cout << "NewTemplate3Dep:: NewTemplate3Dep failed to get copy of elastic material. " << endln;
      exit(1);
    }

    if ( ptr_yield_function_in )
      ptr_yield_function = ptr_yield_function_in->newObj();
    else {
      cout << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the yield function. " << endln;
      exit(1);
    }

    if ( ptr_plastic_flow_in )
       ptr_plastic_flow = ptr_plastic_flow_in->newObj();
    else {
      cout << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the plastic flow. " << endln;
      exit(1);
    }

    // Scalar (isotropic) Evolution
    ptr_scalar_evolution = NULL;

    // Tensor (kinematic) Evolution
    ptr_tensor_evolution = NULL;

    int err;
    err = this->revertToStart();
}

// Constructor
//================================================================================
NewTemplate3Dep::NewTemplate3Dep()
: NDMaterial(0, ND_TAG_NewTemplate3Dep),
  ptr_material_parameter(NULL),
  ptr_elastic_state(NULL),
  ptr_yield_function(NULL),
  ptr_plastic_flow(NULL),
  ptr_scalar_evolution(NULL),
  ptr_tensor_evolution(NULL),
  caseIndex(0), subStep(1), divergeOrnot(0)
{
    int err;
    err = this->revertToStart();
}

// Destructor
//================================================================================
NewTemplate3Dep::~NewTemplate3Dep()
{
    if (ptr_elastic_state)
      delete ptr_elastic_state;

    for (int i = 0; i < ptr_material_parameter->getNum_Internal_Scalar(); i++) {
      if (ptr_scalar_evolution[i])
        delete ptr_scalar_evolution[i];
    }
    if (ptr_scalar_evolution)
      delete [] ptr_scalar_evolution;

    for (int j = 0; j < ptr_material_parameter->getNum_Internal_Tensor(); j++) {
      if (ptr_tensor_evolution[j])
        delete ptr_tensor_evolution[j];
    }
    if (ptr_tensor_evolution)
       delete [] ptr_tensor_evolution;

    if (ptr_yield_function)
       delete ptr_yield_function;

    if (ptr_plastic_flow)
       delete ptr_plastic_flow;

    if (ptr_material_parameter)
       delete ptr_material_parameter;
}

//=================================================================================
// For Matrix based elements
int NewTemplate3Dep::setTrialStrain (const Vector &v)
{
 straintensor temp;

    temp.val(1,1) = v(0);
    temp.val(2,2) = v(1);
    temp.val(3,3) = v(2);
    temp.val(1,2) = 0.5 * v(3);
    temp.val(2,1) = 0.5 * v(3);
    temp.val(3,1) = 0.5 * v(4);
    temp.val(1,3) = 0.5 * v(4);
    temp.val(2,3) = 0.5 * v(5);
    temp.val(3,2) = 0.5 * v(5);

 return this->setTrialStrainIncr(temp - getStrainTensor());
}

//=================================================================================
// For Matrix based elements
int NewTemplate3Dep::setTrialStrain (const Vector &v, const Vector &r)
{
 return this->setTrialStrainIncr(v);
}

//=================================================================================
// For Matrix based elements
int NewTemplate3Dep::setTrialStrainIncr (const Vector &v)
{
 straintensor temp;

    temp.val(1,1) = v(0);
    temp.val(2,2) = v(1);
    temp.val(3,3) = v(2);
    temp.val(1,2) = 0.5 * v(3);
    temp.val(2,1) = 0.5 * v(3);
    temp.val(3,1) = 0.5 * v(4);
    temp.val(1,3) = 0.5 * v(4);
    temp.val(2,3) = 0.5 * v(5);
    temp.val(3,2) = 0.5 * v(5);

    return this->setTrialStrainIncr(temp);
}

//=================================================================================
// For Matrix based elements
int NewTemplate3Dep::setTrialStrainIncr (const Vector &v, const Vector &r)
{
 return this->setTrialStrainIncr(v);
}

//=================================================================================
// For Matrix based elements
const Matrix& NewTemplate3Dep::getTangent (void)
{
   D(0,0) = Stiffness.cval(1,1,1,1);
   D(0,1) = Stiffness.cval(1,1,2,2);
   D(0,2) = Stiffness.cval(1,1,3,3);
   D(0,3) = Stiffness.cval(1,1,1,2) *0.5;
   D(0,4) = Stiffness.cval(1,1,1,3) *0.5;
   D(0,5) = Stiffness.cval(1,1,2,3) *0.5;

   D(1,0) = Stiffness.cval(2,2,1,1);
   D(1,1) = Stiffness.cval(2,2,2,2);
   D(1,2) = Stiffness.cval(2,2,3,3);
   D(1,3) = Stiffness.cval(2,2,1,2) *0.5;
   D(1,4) = Stiffness.cval(2,2,1,3) *0.5;
   D(1,5) = Stiffness.cval(2,2,2,3) *0.5;

   D(2,0) = Stiffness.cval(3,3,1,1);
   D(2,1) = Stiffness.cval(3,3,2,2);
   D(2,2) = Stiffness.cval(3,3,3,3);
   D(2,3) = Stiffness.cval(3,3,1,2) *0.5;
   D(2,4) = Stiffness.cval(3,3,1,3) *0.5;
   D(2,5) = Stiffness.cval(3,3,2,3) *0.5;

   D(3,0) = Stiffness.cval(1,2,1,1);
   D(3,1) = Stiffness.cval(1,2,2,2);
   D(3,2) = Stiffness.cval(1,2,3,3);
   D(3,3) = Stiffness.cval(1,2,1,2) *0.5;
   D(3,4) = Stiffness.cval(1,2,1,3) *0.5;
   D(3,5) = Stiffness.cval(1,2,2,3) *0.5;

   D(4,0) = Stiffness.cval(1,3,1,1);
   D(4,1) = Stiffness.cval(1,3,2,2);
   D(4,2) = Stiffness.cval(1,3,3,3);
   D(4,3) = Stiffness.cval(1,3,1,2) *0.5;
   D(4,4) = Stiffness.cval(1,3,1,3) *0.5;
   D(4,5) = Stiffness.cval(1,3,2,3) *0.5;

   D(5,0) = Stiffness.cval(2,3,1,1);
   D(5,1) = Stiffness.cval(2,3,2,2);
   D(5,2) = Stiffness.cval(2,3,3,3);
   D(5,3) = Stiffness.cval(2,3,1,2) *0.5;
   D(5,4) = Stiffness.cval(2,3,1,3) *0.5;
   D(5,5) = Stiffness.cval(2,3,2,3) *0.5;

   return D;
}



//====================================================================================
// Nima Tafazzoli fixed the starting from (0,0) to (1,1) for BJmatrix (Feb. 2010)
// For Matrix based elements
const BJmatrix& NewTemplate3Dep::getTangentBJmatrix(void)
{
   TangentMatrix.val(1,1) = Stiffness.cval(1,1,1,1);
   TangentMatrix.val(1,2) = Stiffness.cval(1,1,2,2);
   TangentMatrix.val(1,3) = Stiffness.cval(1,1,3,3);
   TangentMatrix.val(1,4) = Stiffness.cval(1,1,1,2) *0.5;
   TangentMatrix.val(1,5) = Stiffness.cval(1,1,1,3) *0.5;
   TangentMatrix.val(1,6) = Stiffness.cval(1,1,2,3) *0.5;
                           
   TangentMatrix.val(2,1) = Stiffness.cval(2,2,1,1);
   TangentMatrix.val(2,2) = Stiffness.cval(2,2,2,2);
   TangentMatrix.val(2,3) = Stiffness.cval(2,2,3,3);
   TangentMatrix.val(2,4) = Stiffness.cval(2,2,1,2) *0.5;
   TangentMatrix.val(2,5) = Stiffness.cval(2,2,1,3) *0.5;
   TangentMatrix.val(2,6) = Stiffness.cval(2,2,2,3) *0.5;
                           
   TangentMatrix.val(3,1) = Stiffness.cval(3,3,1,1);
   TangentMatrix.val(3,2) = Stiffness.cval(3,3,2,2);
   TangentMatrix.val(3,3) = Stiffness.cval(3,3,3,3);
   TangentMatrix.val(3,4) = Stiffness.cval(3,3,1,2) *0.5;
   TangentMatrix.val(3,5) = Stiffness.cval(3,3,1,3) *0.5;
   TangentMatrix.val(3,6) = Stiffness.cval(3,3,2,3) *0.5;
                           
   TangentMatrix.val(4,1) = Stiffness.cval(1,2,1,1);
   TangentMatrix.val(4,2) = Stiffness.cval(1,2,2,2);
   TangentMatrix.val(4,3) = Stiffness.cval(1,2,3,3);
   TangentMatrix.val(4,4) = Stiffness.cval(1,2,1,2) *0.5;
   TangentMatrix.val(4,5) = Stiffness.cval(1,2,1,3) *0.5;
   TangentMatrix.val(4,6) = Stiffness.cval(1,2,2,3) *0.5;
                           
   TangentMatrix.val(5,1) = Stiffness.cval(1,3,1,1);
   TangentMatrix.val(5,2) = Stiffness.cval(1,3,2,2);
   TangentMatrix.val(5,3) = Stiffness.cval(1,3,3,3);
   TangentMatrix.val(5,4) = Stiffness.cval(1,3,1,2) *0.5;
   TangentMatrix.val(5,5) = Stiffness.cval(1,3,1,3) *0.5;
   TangentMatrix.val(5,6) = Stiffness.cval(1,3,2,3) *0.5;
                           
   TangentMatrix.val(6,1) = Stiffness.cval(2,3,1,1);
   TangentMatrix.val(6,2) = Stiffness.cval(2,3,2,2);
   TangentMatrix.val(6,3) = Stiffness.cval(2,3,3,3);
   TangentMatrix.val(6,4) = Stiffness.cval(2,3,1,2) *0.5;
   TangentMatrix.val(6,5) = Stiffness.cval(2,3,1,3) *0.5;
   TangentMatrix.val(6,6) = Stiffness.cval(2,3,2,3) *0.5;


   return TangentMatrix;
}

//=================================================================================
// For Matrix based elements
const Vector& NewTemplate3Dep::getStress (void)
{
   sigma(0) = TrialStress.cval(1,1);
   sigma(1) = TrialStress.cval(2,2);
   sigma(2) = TrialStress.cval(3,3);
   sigma(3) = TrialStress.cval(1,2);
   sigma(4) = TrialStress.cval(1,3);
   sigma(5) = TrialStress.cval(2,3);

   return sigma;
}

//=================================================================================
// For Matrix based elements
const Vector& NewTemplate3Dep::getStrain (void)
{
   epsilon(0) = TrialStrain.cval(1,1);
   epsilon(1) = TrialStrain.cval(2,2);
   epsilon(2) = TrialStrain.cval(3,3);
   epsilon(3) = TrialStrain.cval(1,2) *2.0;
   epsilon(4) = TrialStrain.cval(1,3) *2.0;
   epsilon(5) = TrialStrain.cval(2,3) *2.0;

   return epsilon;
}

//================================================================================
int NewTemplate3Dep::setTrialStrain(const Tensor& v)
{
    return this->setTrialStrainIncr( v - getStrainTensor() );
}


//================================================================================
int NewTemplate3Dep::setTrialStrain(const Tensor& v, const Tensor& r)
{
    return this->setTrialStrainIncr( v - getStrainTensor() );
}

//================================================================================
int NewTemplate3Dep::setTrialStrainIncr(const Tensor& strain_increment)
{
   switch(caseIndex)
     {
       case (0):
         return this->Explicit(strain_increment, this->subStep);

       case (1):
         return this->Implicit(strain_increment, this->subStep);

       case (2):
         return this->ImplicitLineSearch(strain_increment);

       case (3):
         return this->ScaledExplicit(strain_increment, this->subStep);


       default:
         return 0;
     }
}

//================================================================================
int NewTemplate3Dep::setTrialStrainIncr(const Tensor& v, const Tensor& r)
{
    return this->setTrialStrainIncr(v);
}

//================================================================================
double NewTemplate3Dep::getRho(void)
{
    double rho = 0.0;
    if (ptr_material_parameter->getNum_Material_Constant() > 0)
        rho = ptr_material_parameter->getMaterial_Constant(0);
    else {
      cout << "Error!! NewTemplate3Dep:: number of input parameter for material constants less than 1. " << endln;
      cout << "Remind: NewTemplate3Dep:: the 1st material constant is the density. " << endln;
      exit(1);
    }

    return rho;
}

//================================================================================
const BJtensor& NewTemplate3Dep::getTangentTensor(void)
{
    return this->Stiffness;
}

//================================================================================
const stresstensor&  NewTemplate3Dep::getStressTensor(void)
{
    return this->TrialStress;
}


//================================================================================
const straintensor& NewTemplate3Dep::getStrainTensor(void)
{
    return this->TrialStrain;
}

//================================================================================
const straintensor& NewTemplate3Dep::getPlasticStrainTensor(void)
{
    return this->TrialPlastic_Strain;
}

//================================================================================
// Nima
const stresstensor& NewTemplate3Dep::getPredictorStressTensor(void)
{
    return this->TrialPredictor_Stress;
}

//================================================================================
// Nima
const stresstensor& NewTemplate3Dep::getintersection_stress(void)
{
    return this->Trialintersection_stress;
}

//================================================================================
// Nima
const straintensor& NewTemplate3Dep::getintersection_strain(void)
{
    return this->Trialintersection_strain;
}

//================================================================================
// Nima
double NewTemplate3Dep::getEnergy(double ep_cut)
{
	cout << "\n*****************************************************************************\n";
	cout << "\nResults of Direct Calculation: \n\n";


	tStrain = this->getStrainTensor();
	tStress = this->getStressTensor();

	incStrain = tStrain - pStrain;

	incStress = (tStress + pStress) * 0.5;
	TotalMultiply = incStress("ij")*incStrain("ij");
	TotalMultiply.null_indices();


	if (tStrain.val(1,1) == ep_cut) {
	TotalEnergy = 0;
	cout << "\n****** Resetting energy to zero ******" <<endl;
	}


	TotalEnergy = TotalEnergy + TotalMultiply.trace();




//	pStrain.print("pStn" , "pStrain");
	tStrain.print("tStn" , "tStrain");
//	pStress.print("pStr" , "pStress");
//	tStress.print("tStr" , "tStress");

	cout << "\nTotalEnergy= " << TotalEnergy <<endl;
	cout << "\n**************************************************************************\n";


       pStress.val(1,1) = tStress.val(1,1); pStress.val(1,2) = tStress.val(1,2); pStress.val(1,3) = tStress.val(1,3);
       pStress.val(2,1) = tStress.val(2,1); pStress.val(2,2) = tStress.val(2,2); pStress.val(2,3) = tStress.val(2,3);
       pStress.val(3,1) = tStress.val(3,1); pStress.val(3,2) = tStress.val(3,2); pStress.val(3,3) = tStress.val(3,3);


	pStrain = tStrain;

    return TotalEnergy;
}

//================================================================================
// Nima
double NewTemplate3Dep::resetEnergy()
{
//	if (tStrain.val(1,1) == ep_cut) {
//           TotalEnergy = 0;
	cout << "\n****** Resetting energy to zero ******" <<endl;
	return 0;
         }


//================================================================================

// moved out BJ 8Jan2008
// straintensor NewTemplate3Dep::BardetConstraint(int type_of_test, double Increment)
//   {
// //================================================================================
// // for more info look at the paper:
// // @article{Bardet91,
// //  author  = { Bardet, J. P. and Choucair, W. },
// //  title   = { A linearized integration technique for incremental constitutive equations },
// //  journal = { International Journal for Numerical and Analytical Methods in Geomechanics },
// //  year    = { 1991 },
// //    volume   = { 15 },
// //    number   = { 1 },
// //    pages    = { 1-19 },
// //    month    = { },
// //    note     = {  },
// //    napomena     = { local eCM1261 ; 18May2007 },
// // }
// //
// // constant_p_triaxial_loading_strain_control_d_epsilon_11         ->  1
// // drained_triaxial_loading_strain_control_d_epsilon_11            ->  2
// // Undrained_triaxial_loading_strain_control_d_epsilon_11          ->  3
// // Undrained_cyclic_triaxial_loading_stress_control_d_epsilon_11   ->  4
// // Undrained_simple_shear_loading_strain_control_d_epsilon_12      ->  5
// 
//     double K0 = 1.0;  // to be taken care of
// 
//     Matrix CM =  this->getTangent(); 
//     //std::cerr << CM;
// 
//     BJmatrix CurrentStiffness(6,6,0.0);
// 
// 
//    CurrentStiffness.val(1,1) = CM(0,0);
//    CurrentStiffness.val(1,2) = CM(0,1);
//    CurrentStiffness.val(1,3) = CM(0,2);
//    CurrentStiffness.val(1,4) = CM(0,3);
//    CurrentStiffness.val(1,5) = CM(0,4);
//    CurrentStiffness.val(1,6) = CM(0,5);
//                     
//    CurrentStiffness.val(2,1) = CM(1,0);
//    CurrentStiffness.val(2,2) = CM(1,1);
//    CurrentStiffness.val(2,3) = CM(1,2);
//    CurrentStiffness.val(2,4) = CM(1,3);
//    CurrentStiffness.val(2,5) = CM(1,4);
//    CurrentStiffness.val(2,6) = CM(1,5);
//                                       
//    CurrentStiffness.val(3,1) = CM(2,0);
//    CurrentStiffness.val(3,2) = CM(2,1);
//    CurrentStiffness.val(3,3) = CM(2,2);
//    CurrentStiffness.val(3,4) = CM(2,3);
//    CurrentStiffness.val(3,5) = CM(2,4);
//    CurrentStiffness.val(3,6) = CM(2,5);
//                                       
//    CurrentStiffness.val(4,1) = CM(3,0);
//    CurrentStiffness.val(4,2) = CM(3,1);
//    CurrentStiffness.val(4,3) = CM(3,2);
//    CurrentStiffness.val(4,4) = CM(3,3);
//    CurrentStiffness.val(4,5) = CM(3,4);
//    CurrentStiffness.val(4,6) = CM(3,5);
//                                       
//    CurrentStiffness.val(5,1) = CM(4,0);
//    CurrentStiffness.val(5,2) = CM(4,1);
//    CurrentStiffness.val(5,3) = CM(4,2);
//    CurrentStiffness.val(5,4) = CM(4,3);
//    CurrentStiffness.val(5,5) = CM(4,4);
//    CurrentStiffness.val(5,6) = CM(4,5);
//                                       
//    CurrentStiffness.val(6,1) = CM(5,0);
//    CurrentStiffness.val(6,2) = CM(5,1);
//    CurrentStiffness.val(6,3) = CM(5,2);
//    CurrentStiffness.val(6,4) = CM(5,3);
//    CurrentStiffness.val(6,5) = CM(5,4);
//    CurrentStiffness.val(6,6) = CM(5,5);
//                                  
//     //CurrentStiffness.print("CS","CurrentStiffness"); 
// 
// //    CurrentStiffness.print("C","C");
// 
// 
//     double * p_S_values = NULL;
//     double * p_E_values = NULL;
// 
// 
// // constant_p_triaxial_loading_strain_control_d_epsilon_11
//     if  ( type_of_test == 1 )
//       {
//         ////////////////////////////////////////////////////////////////////
//         // constant_p_triaxial_loading_strain_control_d_epsilon_11
//           double S_values[] =
//         {1.0  ,   1.0   ,   1.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   1.0   ,  -1.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   1.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   1.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   1.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0};
//           double E_values[] =
//         {0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          1.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0};
// 
//          p_S_values = S_values;
//          p_E_values = E_values;
// 
//        }
// 
//     else if ( type_of_test == 2 ) 
//       {
//         ////////////////////////////////////////////////////////////////////
//         // drained_triaxial_loading_strain_control_d_epsilon_11
//           double S_values[] =
//         {0.0  ,   1.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   1.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   1.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   1.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   1.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0};
//           double E_values[] =
//         {0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          1.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0};
//          p_S_values = S_values;
//          p_E_values = E_values;
// 
//       }
// 
//     else if ( type_of_test == 3 ) 
//       {
//         ////////////////////////////////////////////////////////////////////
//         // Undrained_triaxial_loading_strain_control_d_epsilon_11
//           double S_values[] =
//         {0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   1.0   ,  -1.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   1.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   1.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   1.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0};
//           double E_values[] =
//         {1.0  ,   1.0   ,   1.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          1.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0};
//          p_S_values = S_values;
//          p_E_values = E_values;
// 
//       }
// 
//     else if ( type_of_test == 4 ) 
//       {
//         ////////////////////////////////////////////////////////////////////
//         // Undrained_cyclic_triaxial_loading_stress_control_d_epsilon_11
//           double S_values[] =
//         {0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   1.0   ,  -1.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   1.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   1.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   1.0,
//          1.0  ,  -1.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0};
//           double E_values[] =
//         {1.0  ,   1.0   ,   1.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0};
//          p_S_values = S_values;
//          p_E_values = E_values;
// 
//       }
// 
//     else if ( type_of_test == 5 ) 
//       {
//         ////////////////////////////////////////////////////////////////////
//         // Undrained_simple_shear_loading_strain_control_d_epsilon_12
//           double S_values[] =
//         {0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          K0   ,  -1.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   1.0   ,  -1.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   1.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   1.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0};
//           double E_values[] =
//         {1.0  ,   1.0   ,   1.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0,
//          0.0  ,   0.0   ,   0.0   ,   1.0   ,   0.0   ,   0.0};
//          p_S_values = S_values;
//          p_E_values = E_values;
// 
//       }
// 
// 
//     BJmatrix S(6,6, p_S_values);
//     
//     //S.print("S","S");
// 
//     BJmatrix E(6,6, p_E_values);
// 
//     //E.print("E","E");
// 
// 
//     double number = 0.0;
// // vector constructor
//     BJvector Y(6, 0.0);
//     Y.val(6) =  Increment;
// 
//     //Y.print("Y","Y");
// //    ::fprintf(stderr," line 1 ");
// 
//     BJmatrix S_times_C_plus_E = S * CurrentStiffness + E;
// 
//     //S_times_C_plus_E.print("SCE","S_times_C_plus_E");
// 
// 
// //    ::fprintf(stderr," line 2 ");
//     BJmatrix Inverse_S_times_C_plus_E = S_times_C_plus_E.inverse();
//     //Inverse_S_times_C_plus_E.print("I_SCE","Inverse_S_times_C_plus_E");
// 
// //    ::fprintf(stderr," line 3 ");
//     BJvector d_epsilon = Inverse_S_times_C_plus_E * Y;
//     //d_epsilon.print("d_e","d_epsilon"); 
// 
// 
//     straintensor d_epsilon_tensor;
//     //d_epsilon_tensor.print("de");
// 
//     double sqrthalf = sqrt(0.5);
// 
//   // Adopted method from Helnwein (2001):
//     d_epsilon_tensor.val(1,1) = d_epsilon.val(1);
//     d_epsilon_tensor.val(2,2) = d_epsilon.val(2);
//     d_epsilon_tensor.val(3,3) = d_epsilon.val(3);
//     d_epsilon_tensor.val(1,2) = d_epsilon_tensor.val(2,1) = d_epsilon.val(4) *sqrthalf;
//     d_epsilon_tensor.val(2,3) = d_epsilon_tensor.val(3,2) = d_epsilon.val(5) *sqrthalf;
//     d_epsilon_tensor.val(1,3) = d_epsilon_tensor.val(3,1) = d_epsilon.val(6) *sqrthalf;
// 
//     //d_epsilon_tensor.print("de");
// 
// //out-  getchar();
// 
//     return d_epsilon_tensor;
//   }

//================================================================================
int NewTemplate3Dep::commitState(void)
{
    int err = 0;

    CommitStress.Initialize(TrialStress);
    CommitStrain.Initialize(TrialStrain);
    CommitPlastic_Strain.Initialize(TrialPlastic_Strain);

    return err;
}

//================================================================================
int NewTemplate3Dep::revertToLastCommit(void)
{
    int err = 0;

    TrialStress.Initialize(CommitStress);
    TrialStrain.Initialize(CommitStrain);
    TrialPlastic_Strain.Initialize(CommitPlastic_Strain);

    return err;
}

//================================================================================
int NewTemplate3Dep::revertToStart(void)
{
    int err = 0;

    CommitStress = ptr_elastic_state->getStress();
    CommitStrain = ptr_elastic_state->getStrain();
    CommitPlastic_Strain.Initialize(ZeroStrain);

    TrialStress.Initialize(CommitStress);
    TrialStrain.Initialize(CommitStrain);
    TrialPlastic_Strain.Initialize(ZeroStrain);

    Stiffness = ptr_elastic_state->getElasticStiffness(*ptr_material_parameter);

    return err;
}

//================================================================================
NDMaterial * NewTemplate3Dep::getCopy(void)
{
   NDMaterial* tmp = new NewTemplate3Dep(this->getTag(),
                                         this->ptr_material_parameter,
                                         this->ptr_elastic_state,
                                         this->ptr_yield_function,
                                         this->ptr_plastic_flow,
                                         this->ptr_scalar_evolution,
                                         this->ptr_tensor_evolution,
                                         this->caseIndex,
                                         this->subStep);
    return tmp;
}


//================================================================================
NDMaterial * NewTemplate3Dep::getCopy(const char *code)
{
    if (strcmp(code,"ThreeDimensional") == 0) {
       NewTemplate3Dep* tmp = new NewTemplate3Dep( this->getTag(),
                                                   this->ptr_material_parameter,
                                                   this->ptr_elastic_state,
                                                   this->ptr_yield_function,
                                                   this->ptr_plastic_flow,
                                                   this->ptr_scalar_evolution,
                                                   this->ptr_tensor_evolution,
                                                   this->caseIndex,
                                                   this->subStep);
       return tmp;
    }
    else {
      cout << "NewTemplate3Dep::getCopy failed to get model: " <<  code << endln;
      exit(1);
    }
    return 0;
}

//================================================================================
const char *NewTemplate3Dep::getType(void) const
{
    return "ThreeDimensional";
}

//================================================================================
int NewTemplate3Dep::sendSelf(int commitTag, Channel &theChannel)
{
    //Guanzhou implemented for parallel processing
# ifdef _PARALLEL_PROCESSING
    static ID idData(13);

    idData.Zero();

    idData(0) = this->getTag();

    if ( ptr_material_parameter != NULL ) idData(1)  = 1;
    if ( ptr_elastic_state != NULL )      idData(2)  = 1;
    if ( ptr_yield_function != NULL )     idData(3)  = 1;
    if ( ptr_plastic_flow != NULL )       idData(4)  = 1;
    if ( ptr_scalar_evolution != NULL )   idData(5)  = 1;
    if ( ptr_tensor_evolution != NULL )   idData(6)  = 1;

    idData(7) = caseIndex;

    if ( ptr_elastic_state != NULL )      idData(8)  = ptr_elastic_state->getClassTag();
    if ( ptr_yield_function != NULL )     idData(9)  = ptr_yield_function ->getClassTag();
    if ( ptr_plastic_flow != NULL )       idData(10) = ptr_plastic_flow->getClassTag();

 idData(11) = subStep;
 idData(12) = divergeOrnot;

    static ID *pID_ScalarTags = NULL;
    if ( ptr_scalar_evolution != NULL ) {
     const int NumScalar = ptr_material_parameter->getNum_Internal_Scalar();
     pID_ScalarTags = new ID(NumScalar);
 for (int i=0; i<NumScalar; i++) (*pID_ScalarTags)(i) = ptr_scalar_evolution[i]->getClassTag();
    }

    static ID *pID_TensorTags = NULL;
    if ( ptr_tensor_evolution != NULL ) {
     const int NumTensor = ptr_material_parameter->getNum_Internal_Tensor();
     pID_TensorTags = new ID(NumTensor);
 for (int i=0; i<NumTensor; i++) (*pID_TensorTags)(i) = ptr_tensor_evolution[i]->getClassTag();
    }

    if (theChannel.sendID(this->getDbTag(), commitTag, idData) < 0) {
    std::cerr << "NewTemplate3Dep::sendSelf -- failed to send ID\n";
    return -1;
    }

    if ( ptr_material_parameter != NULL )
     if (ptr_material_parameter->sendSelf(commitTag, theChannel) < 0) {
         std::cerr << "NewTemplate3Dep::sendSelf -- MaterialParameters failed to send self\n";
     return -1;
        }

    if ( ptr_elastic_state != NULL )
 if (ptr_elastic_state->sendSelf(commitTag, theChannel) < 0) {
         std::cerr << "NewTemplate3Dep::sendSelf -- ElasticState failed to send self\n";
     return -1;
        }

    if ( ptr_yield_function != NULL )
        if (ptr_yield_function->sendSelf(commitTag, theChannel) < 0) {
         std::cerr << "NewTemplate3Dep::sendSelf -- YieldFunction failed to send self\n";
     return -1;
       }

    if ( ptr_plastic_flow != NULL )
        if (ptr_plastic_flow->sendSelf(commitTag, theChannel) < 0) {
         std::cerr << "NewTemplate3Dep::sendSelf -- PlasticFlow failed to send self\n";
     return -1;
        }

    if ( ptr_scalar_evolution != NULL ) {
     if (theChannel.sendID(this->getDbTag(), commitTag, *pID_ScalarTags) < 0) {
        std::cerr << "NewTemplate3Dep::sendSelf -- failed to send pID_ScalarTags\n";
        return -1;
     }

 delete pID_ScalarTags;

 for (int i=0; i<ptr_material_parameter->getNum_Internal_Scalar(); i++) {
     if (ptr_scalar_evolution[i]->sendSelf(commitTag, theChannel) < 0) {
          std::cerr << "NewTemplate3Dep::sendSelf -- ScalarEvolution failed to send self\n";
      return -1;
            }
 }
    }


    if ( ptr_tensor_evolution != NULL ) {
     if (theChannel.sendID(this->getDbTag(), commitTag, *pID_TensorTags) < 0) {
        std::cerr << "NewTemplate3Dep::sendSelf -- failed to send pID_TensorTags\n";
        return -1;
     }
     delete pID_TensorTags;

 for (int i=0; i<ptr_material_parameter->getNum_Internal_Tensor(); i++) {
            if (ptr_tensor_evolution[i]->sendSelf(commitTag, theChannel) < 0) {
          std::cerr << "NewTemplate3Dep::sendSelf -- TensorEvolution failed to send self\n";
      return -1;
            }
 }
    }
# endif
    return 0;
}

//================================================================================
int NewTemplate3Dep::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    //Guanzhou implemented for parallel processing
# ifdef _PARALLEL_PROCESSING
    int dataTag = this->getDbTag();
    static ID idData(13);
    static ID *pID_ScalarTags = NULL, *pID_TensorTags=NULL;
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
     std::cerr << "NewTemplate3Dep::recvSelf -- failed to recv ID\n";
 return -1;
    }

    this->setTag(idData(0));

    if ( idData(1) != 0) {

 //if ( ptr_material_parameter != NULL ) delete ptr_material_parameter;

 ptr_material_parameter = theBroker.getNewMaterialParameterPtr();

 if (ptr_material_parameter->recvSelf(commitTag, theChannel, theBroker) < 0) {
         std::cerr << "NewTemplate3Dep::recvSelf -- MaterialParameter failed to recv self\n";
     return -1;
        }

    }

    if ( idData(2) != 0) {

 //if ( ptr_elastic_state != NULL ) delete ptr_elastic_state;

 ptr_elastic_state = theBroker.getNewElasticStatePtr(idData(8));

 if (ptr_elastic_state->recvSelf(commitTag, theChannel, theBroker) < 0) {
         std::cerr << "NewTemplate3Dep::recvSelf -- ElasticState failed to recv self\n";
     return -1;
        }

    }

    if ( idData(3) != 0) {

 //if ( ptr_yield_function != NULL ) delete ptr_yield_function;

 ptr_yield_function = theBroker.getNewYieldFunctionPtr(idData(9));

 if (ptr_yield_function->recvSelf(commitTag, theChannel, theBroker) < 0) {
         std::cerr << "NewTemplate3Dep::recvSelf -- YieldFunction failed to recv self\n";
     return -1;
        }

    }

    if ( idData(4) != 0) {

 //if ( ptr_yield_function != NULL ) delete ptr_yield_function;

 ptr_plastic_flow = theBroker.getNewPlasticFlowPtr(idData(10));

 if (ptr_plastic_flow->recvSelf(commitTag, theChannel, theBroker) < 0) {
         std::cerr << "NewTemplate3Dep::recvSelf -- PlasticFlow failed to recv self\n";
     return -1;
        }

    }

    if ( idData(5) != 0) {

 const int NumScalar = ptr_material_parameter->getNum_Internal_Scalar();
 //static ID *
 pID_ScalarTags = new ID(NumScalar);
     if (theChannel.recvID(dataTag, commitTag, *pID_ScalarTags) < 0) {
         std::cerr << "NewTemplate3Dep::recvSelf -- failed to recv pID_ScalarTags\n";
     return -1;
     }

       ptr_scalar_evolution = new ScalarEvolution* [NumScalar];
        for (int i = 0; i < NumScalar; i++) {
            ptr_scalar_evolution[i] = theBroker.getNewScalarEvolutionPtr((*pID_ScalarTags)(i));
     if (ptr_scalar_evolution[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
          std::cerr << "NewTemplate3Dep::recvSelf -- ScalarEvolution failed to recv self\n";
      return -1;
            }
 }

 delete pID_ScalarTags;

    }

    if ( idData(6) != 0) {

 const int NumTensor = ptr_material_parameter->getNum_Internal_Tensor();
 //static ID *
 pID_TensorTags = new ID(NumTensor);
     if (theChannel.recvID(dataTag, commitTag, *pID_TensorTags) < 0) {
         std::cerr << "NewTemplate3Dep::recvSelf -- failed to recv pID_TensorTags\n";
     return -1;
     }

       ptr_tensor_evolution = new TensorEvolution* [NumTensor];
        for (int i = 0; i < NumTensor; i++) {
            ptr_tensor_evolution[i] = theBroker.getNewTensorEvolutionPtr((*pID_TensorTags)(i));
     if (ptr_tensor_evolution[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
          std::cerr << "NewTemplate3Dep::recvSelf -- TensorEvolution failed to recv self\n";
      return -1;
            }
 }

 delete pID_TensorTags;

    }

    caseIndex = idData(7);
 subStep = idData(11);
 divergeOrnot = idData(12);

 this->revertToStart();

# endif
    return 0;

}

//================================================================================
void NewTemplate3Dep::Print(OPS_Stream& s, int flag)
{
     s << (*this);
}

//================================================================================
int NewTemplate3Dep::Explicit(const straintensor& strain_incr, int NumStep_in)
{
    int err = 0;
    int i = 0;
    int iStep = 0;
    double f_start = 0.0;
    double f_pred  = 0.0;

    BJtensor Ee(4, def_dim_4, 0.0);
    BJtensor Ep(4, def_dim_4, 0.0);

    straintensor intersection_strain;
    stresstensor intersection_stress;
    double intersection_factor = 0.0;

    int Num_internal_scalar = 0;
    int Num_internal_tensor = 0;
    int Num_internal_scalar_YF = 0;
    int Num_internal_tensor_YF = 0;

    double upper = 0.0;
    double lower = 0.0;
    double Delta_lambda = 0.0;
    double hardMod  = 0.0;
    double h_s = 0.0;
    double xi_s = 0.0;
    stresstensor dFods;
    straintensor dQods;
    BJtensor Hq(2, def_dim_2, 0.0);
    BJtensor Hf(2, def_dim_2, 0.0);

    stresstensor h_t;
    stresstensor xi_t;

    straintensor incr_strain;
    stresstensor incr_stress;
    straintensor incr_Pstrain;
    stresstensor ep_stress;
    stresstensor predicted_stress;

    stresstensor start_stress;
    stresstensor start_strain;
    straintensor start_Pstrain;

  for (iStep = 1; iStep <= NumStep_in; iStep++) {
    start_stress = getStressTensor();
    start_strain = getStrainTensor();
    start_Pstrain = getPlasticStrainTensor();

    intersection_stress.Initialize(start_stress);

    err += ptr_elastic_state->setStress(start_stress);
    err += ptr_elastic_state->setStrain(start_strain);
    Ee = ptr_elastic_state->getElasticStiffness(*ptr_material_parameter);

//      Ee.print("Ee" , "Ee");

    incr_strain = strain_incr *(1.0/(double)(NumStep_in));

//      incr_strain.print("incr_strain" , "incr_strain");


    incr_stress = Ee("ijpq") * incr_strain("pq");
       incr_stress.null_indices();
       
//       incr_stress.print("incr_stress" , "incr_stress");


    predicted_stress = start_stress + incr_stress;

/////////////////////////////////////////////////////
//        start_stress.print("start_stress" , "start_stress");
//        incr_stress.print("incr_stress" , "incr_stress");
//        predicted_stress.print("predicted_stress" , "predicted_stress");
/////////////////////////////////////////////////////

// Nima
//    TrialPredictor_Stress = predicted_stress;



    f_start = ptr_yield_function->YieldFunctionValue( start_stress, *ptr_material_parameter );
    f_pred =  ptr_yield_function->YieldFunctionValue( predicted_stress, *ptr_material_parameter );

/////////////////////////////////////////////////////
//     cout << "f_start = " << f_start << endl;
//     cout << "f_pred = " << f_pred << endl;
/////////////////////////////////////////////////////
# ifdef _TEACHING_MODE
//    fprintf(stdout,"explicit f_start = %12.4e",f_start);
//    fprintf(stdout,"explicit f_pred = %12.4e",f_pred);
# endif


    // If Elastic
    if ( (f_start < 0.0 && f_pred <= FTOL) || f_start > f_pred ) {
      TrialStrain = start_strain + incr_strain;
      TrialStress.Initialize(predicted_stress);

// Nima
//      Trialintersection_stress = TrialStress;
//      Trialintersection_strain = TrialStrain;



      if (iStep == NumStep_in)
        Stiffness = Ee;
    }
    else {
      // If Elastic and then Elastic-Plastic
      if ( f_start < 0.0 )  {
        intersection_factor = zbrentstress( start_stress, predicted_stress, 0.0, 1.0, TOL );
        intersection_stress = yield_surface_cross( start_stress, predicted_stress, intersection_factor );
        intersection_strain = start_strain + (incr_strain * intersection_factor);

        incr_stress = predicted_stress - intersection_stress;  // necessary


// Nima
//      Trialintersection_stress = intersection_stress;
//      Trialintersection_strain = intersection_strain;


        err += ptr_elastic_state->setStress(intersection_stress);
        err += ptr_elastic_state->setStrain(intersection_strain);
        Ee = ptr_elastic_state->getElasticStiffness(*ptr_material_parameter);
      }



      // If E-P Response,
      Delta_lambda = 0.0;


      dFods = ptr_yield_function->StressDerivative( intersection_stress, *ptr_material_parameter );
      dQods = ptr_plastic_flow->PlasticFlowTensor( intersection_stress, intersection_strain, *ptr_material_parameter );


//        dFods.print("nij" , "nij");
//        dQods.print("mij" , "mij");

//    stresstensor Tbefore = ptr_material_parameter->getInternal_Tensor(0);
//    Tbefore.print("Tbeginning" , "Tbeginning");
    
      // E_ijkl * R_kl
      Hq = Ee("ijkl") * dQods("kl");
        Hq.null_indices();

      // L_ij * E_ijkl
      Hf = dFods("ij") * Ee("ijkl");
        Hf.null_indices();

//      Hf.print("Hf" , "Hf");

      // L_ij * E_ijkl * d e_kl ( true EP strain increment)
      upper = ( dFods("ij") * incr_stress("ij") ).trace();

      // L_ij * E_ijkl * R_kl
      lower = ( Hf("ij") * dQods("ij") ).trace();

//      cout << "lower = " << lower << endl;


      hardMod  = 0.0;

      // Evolution of scalar (isotropic) internal variables in yield function
      Num_internal_scalar_YF = ptr_yield_function->getNumInternalScalar();
      for (i = 0; i < Num_internal_scalar_YF; i++) {
        h_s = ptr_scalar_evolution[i]->H(*ptr_plastic_flow, intersection_stress, intersection_strain, *ptr_material_parameter);
        xi_s = ptr_yield_function->InScalarDerivative( intersection_stress, *ptr_material_parameter, i+1);
        hardMod += h_s * xi_s;
      }

      // Evolution of tensor (kinematic) internal variables in yield function
      Num_internal_tensor_YF = ptr_yield_function->getNumInternalTensor();
      for (i = 0; i < Num_internal_tensor_YF; i++) {
        h_t = ptr_tensor_evolution[i]->Hij(*ptr_plastic_flow, intersection_stress, intersection_strain, *ptr_material_parameter);
        xi_t = ptr_yield_function->InTensorDerivative( intersection_stress, *ptr_material_parameter, i+1);
        hardMod += ( h_t("mn") * xi_t("mn") ).trace();
      }

//     cout << "\n upper = " << upper <<endl;
//     cout << "\n hardMod = " << hardMod <<endl;


      lower -= hardMod;

      Delta_lambda = upper / lower;

      if (Delta_lambda < 0.0)
        Delta_lambda = 0.0;


// Nima
/////////////////////////////////////////////////////

//      dFods.print("dFods" , "dFods");
//      dQods.print("dQods" , "dQods");
//      cout << "\nnominator = " << upper << "  " << "denominator = " << lower << endl;
//       cout << "\n upper = " << upper <<endl;
//       cout << " lower = " << lower <<endl;
//       cout << " Delta_lambda = " << Delta_lambda <<endl;
/////////////////////////////////////////////////////

      // Plastic strain increment
      incr_Pstrain = dQods * Delta_lambda;

      ep_stress = predicted_stress - (Hq *Delta_lambda);

      TrialStress.Initialize(ep_stress);
      TrialStrain = start_strain + incr_strain;
      TrialPlastic_Strain = start_Pstrain + incr_Pstrain;

      // Update internal scalar variables
      Num_internal_scalar = ptr_material_parameter->getNum_Internal_Scalar();
      for (i = 0; i < Num_internal_scalar; i++) {
        double dS = ( ptr_scalar_evolution[i]->H(*ptr_plastic_flow, intersection_stress, intersection_strain, *ptr_material_parameter) ) *Delta_lambda;
        double S = ptr_material_parameter->getInternal_Scalar(i);
        err += ptr_material_parameter->setInternal_Scalar(i, S + dS );
	
//         double S1 = ptr_material_parameter->getInternal_Scalar(i);
// 
// 	cout << " k1 = " << S << "  " << " k2 = " << S1 << " dk = " << dS << endl;

      }

      // Update internal tensor variables
      Num_internal_tensor = ptr_material_parameter->getNum_Internal_Tensor();
      for (i = 0; i < Num_internal_tensor; i++) {
        stresstensor dT = ptr_tensor_evolution[i]->Hij(*ptr_plastic_flow, intersection_stress, intersection_strain, *ptr_material_parameter) *Delta_lambda;
        stresstensor T = ptr_material_parameter->getInternal_Tensor(i);
	
// 	dT.print("dT" , "dT");
// 	cout << " Delta_lambda = " << Delta_lambda <<endl;
// 	T.print("Tbeforeset" , "Tbeforeset");
	
        err += ptr_material_parameter->setInternal_Tensor(i, T + dT );
	
//	stresstensor Tafter = ptr_material_parameter->getInternal_Tensor(0);
// 	Tafter.print("Tafterset" , "Tafterset");
//	Tafter.print("Tafterset" , "Tafterset");

      }

      // To obtain Eep, at the last step
      if (iStep == NumStep_in) {
        Ep = Hq("pq") * Hf("mn");
          Ep.null_indices();
        Ep = Ep * (1.0/lower);

        if ( Delta_lambda > 0.0 )
          Stiffness = Ee - Ep;
        else
          Stiffness = Ee;
      }

    }

  }

//     straintensor start_stressend;
//     start_stressend  = getStressTensor();
//     start_stressend.print("start_stressend" , "start_stressend");
    
    
  return err;
}

//================================================================================
int NewTemplate3Dep::Implicit(const straintensor& strain_incr, int NumStep_in)
{
    // Mahdi Taiebat & Boris Jeremic 22April2007
    // these variables should not be defined as static!

    int err = 0;
    int i = 0;
    int j = 0;
    int iStep = 0;
    double YieldFun = 0.0;

    BJtensor Ee(4, def_dim_4, 0.0);
    //BJtensor Ee1(4, def_dim_4, 0.0);
    //BJtensor Ee2(4, def_dim_4, 0.0);
    BJtensor Ce(4, def_dim_4, 0.0);

    int Num_internal_scalar_YF = 0;
    int Num_internal_tensor_YF = 0;
    int Num_internal_scalar = 0;
    int Num_internal_tensor = 0;

    stresstensor start_stress;
    stresstensor start_strain;
    straintensor start_Pstrain;

// Nima
    straintensor incr_Pstrain;

    straintensor incr_strain;
    stresstensor incr_stress;
    stresstensor PredictedStress;

    double start_InScalar = 0.0;
    tensor start_InTensor(2, def_dim_2, 0.0);
    tensor start_InTensor2(2, def_dim_2, 0.0);

    this->divergeOrnot = 0;

  for (iStep = 1; iStep <= NumStep_in; iStep++)
  {
    start_stress = getStressTensor();
    start_strain = getStrainTensor();
    start_Pstrain = getPlasticStrainTensor();

    Num_internal_scalar = ptr_material_parameter->getNum_Internal_Scalar();
    if (Num_internal_scalar > 1)
      {
        cout << "Error, NewTemplate3Dep::Backward Agorithm, scalar internal variables more than 1!" << endln;
        exit(1);
      }
    else if (Num_internal_scalar == 1)
      {
        start_InScalar = ptr_material_parameter->getInternal_Scalar(0);
      }

    Num_internal_tensor = ptr_material_parameter->getNum_Internal_Tensor();
    if (Num_internal_tensor > 2)
      {
        cout << "Error, NewTemplate3Dep::Backward Agorithm, scalar internal variables more than 1!" << endln;
        exit(1);
      }
    else if (Num_internal_tensor == 1)
      {
       start_InTensor = ptr_material_parameter->getInternal_Tensor(0);
      }
    else if (Num_internal_tensor == 2)
      {
        start_InTensor = ptr_material_parameter->getInternal_Tensor(0);
        start_InTensor2 = ptr_material_parameter->getInternal_Tensor(1);
      }

    err += ptr_elastic_state->setStress(start_stress);
    err += ptr_elastic_state->setStrain(start_strain);
    // Ee1
    Ee = ptr_elastic_state->getElasticStiffness(*ptr_material_parameter);

    incr_strain = strain_incr *(1.0/(double)(NumStep_in));

    // Ee1
    incr_stress = Ee("ijpq") * incr_strain("pq");
      incr_stress.null_indices();
    PredictedStress = start_stress + incr_stress;

    TrialStrain = start_strain + incr_strain;

    //err += ptr_elastic_state->setStress(PredictedStress);
    //err += ptr_elastic_state->setStrain(TrialStrain);
    //Ee2 = ptr_elastic_state->getElasticStiffness(*ptr_material_parameter);

    //Ee = (Ee1 + Ee2) *0.5;
    //err += Stiffness2Compliance(Ee, Ce);

    incr_stress = Ee("ijpq") * incr_strain("pq");
      incr_stress.null_indices();
    PredictedStress = start_stress + incr_stress;

# ifdef _TEACHING_MODE
//    PredictedStress.report(" NewTemplate3Dep::Implicit(  PredictedStress");
#endif

    TrialStress.Initialize(PredictedStress);
    TrialPlastic_Strain.Initialize(start_Pstrain);

    YieldFun = ptr_yield_function->YieldFunctionValue(TrialStress, *ptr_material_parameter);
    //cout << "YieldFun = " << YieldFun << endln;

    if ( YieldFun <= FTOL )
// ELASTIC  ELASTIC  ELASTIC   ELASTIC  ELASTIC ELASTIC
      {
        if (iStep == NumStep_in)
          {
            Stiffness = Ee;
          }
// ELASTIC  ELASTIC  ELASTIC   ELASTIC  ELASTIC ELASTIC
      }
    else
// ELASTIC-PLASTIC  ELASTIC-PLASTIC  ELASTIC-PLASTIC  ELASTIC-PLASTIC  ELASTIC-PLASTIC  ELASTIC-PLASTIC
      {
        // Mahdi Taiebat & Boris Jeremic 20April2007
	// these variables should not be defined as static (e.g. it cause problem for resetting iter_counter=0.0)
        double YieldFun1 = 0.0;
        double YieldFun2 = 0.0;
        double YieldFun_relative = 0.0;

        double Delta_lambda  = 0.0;
        double d2_lambda  = 0.0;
        double RNorm = 0.0;
        double RNorm_relative = 0.0;
        int iter_counter = 0;
        double upper = 0.0;
        double lower = 0.0;

        stresstensor dFods;
        double dFodq = 0.0;
        stresstensor dFoda;
        stresstensor dFoda2;

// Nima
	straintensor dQods;

        straintensor m;
        double hi = 0.0;
        straintensor hk;
        straintensor hk2;

        stresstensor Dstress;
        double Dq = 0.0;
        tensor Da(2, def_dim_2, 0.0);
        tensor Da2(2, def_dim_2, 0.0);

        stresstensor Rstress;
        double Rq = 0.0;
        tensor Ra(2, def_dim_2, 0.0);
        tensor Ra2(2, def_dim_2, 0.0);

        tensor dmods(4, def_dim_4, 0.0);

        straintensor DPlastic_Strain;

        Matrix M66(6, 6);  // intermediate
        Vector V6(6);      // intermediate
        Vector T6(6);      // intermediate

        double S = 0.0; //  scalar internal variable (temp) Mahdi Taiebat & Boris Jeremic 20April2007
        stresstensor T; // first tensor internal variable  Mahdi Taiebat & Boris Jeremic 20April2007
        stresstensor T2; // second tensor internal variable  Mahdi Taiebat & Boris Jeremic 20April2007

        int numMV = 6 + Num_internal_scalar + Num_internal_tensor*6;

        Matrix CC(numMV, numMV);  // BIG matrix
        Matrix CI(numMV, numMV);  // inverse of CI
        Vector R19(numMV);        // Residual
        Vector N19(numMV);        // n, q, a
        Vector M19(numMV);        // m, hi, hk
        Vector D19(numMV);        // incremental
        Vector T19(numMV);        // intermediate

        BJtensor Ttemp(4, def_dim_4, 0.0);

//        YieldFun_relative = YieldFun*FTOL > TOL ? YieldFun*FTOL : TOL;
        YieldFun_relative = YieldFun*FTOL > TOL ? YieldFun*FTOL : TOL;
        YieldFun1 = YieldFun;

      // ################## Beginning of do-while ########################
        do
          {
         //cout << "iter_counter = " << iter_counter << endln;

        // M
            M19.Zero();

            m = ptr_plastic_flow->PlasticFlowTensor(TrialStress, TrialStrain, *ptr_material_parameter);

            m =  Ee("ijpq") * m("pq");
              m.null_indices();

            err += Tensor2VectorSysR2(m, V6);
            for (i=0; i<6; i++)
              M19(i) = V6(i);

            if (Num_internal_scalar == 1)
              {
                hi =
                ptr_scalar_evolution[0]->H(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);

                M19(6) = -hi;
               }

            if (Num_internal_tensor >= 1)
              {
                 hk =
                 ptr_tensor_evolution[0]->Hij(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);

                 err += Tensor2VectorSysR2(hk, V6);

                 for (i=0; i<6; i++)
                   M19(6+Num_internal_scalar+i) = -V6(i);
               }

             if (Num_internal_tensor == 2)
               {
                 hk2 =
                 ptr_tensor_evolution[1]->Hij(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);

                 err += Tensor2VectorSysR2(hk2, V6);

                 for (i=0; i<6; i++)
                   M19(12+Num_internal_scalar+i) = -V6(i);
               }

           // N
              N19.Zero();

              dFods = ptr_yield_function->StressDerivative(TrialStress, *ptr_material_parameter);

   //           dFods.print("dFods","\n\n\n\n dFods on line 1228 of NewTemplate3dEP \n ");

              err += Tensor2VectorSysR2(dFods, V6);

              for (i=0; i<6; i++)
                 N19(i) = V6(i);

              Num_internal_scalar_YF = ptr_yield_function->getNumInternalScalar();
              if (Num_internal_scalar_YF > 1)
                {
                   cout << "Error, NewTemplate3Dep::Backward Agorithm, scalar internal variables more than 1!" << endln;
                   exit(1);
                }

              else if (Num_internal_scalar_YF == 1)
                {
                  dFodq = ptr_yield_function->InScalarDerivative(TrialStress, *ptr_material_parameter, 1);
                  N19(6) = dFodq;
                }

              Num_internal_tensor_YF = ptr_yield_function->getNumInternalTensor();

              if (Num_internal_tensor_YF > 2)
                {
                   cout << "Error, NewTemplate3Dep::Backward Agorithm, tensor internal variables more than 1!" << endln;
                   exit(1);
                 }
              else if (Num_internal_tensor_YF == 1)
                {
                  dFoda = ptr_yield_function->InTensorDerivative(TrialStress, *ptr_material_parameter, 1);
                  err += Tensor2VectorSysR2(dFoda, V6);
                  for (i=0; i<6; i++)
                    N19(6+Num_internal_scalar+i) = V6(i);
                }
           else if (Num_internal_tensor_YF == 2)
             {
                 dFoda = ptr_yield_function->InTensorDerivative(TrialStress, *ptr_material_parameter, 1);
                 err += Tensor2VectorSysR2(dFoda, V6);
                 for (i=0; i<6; i++)
                   N19(6+Num_internal_scalar+i) = V6(i);

                 dFoda2 = ptr_yield_function->InTensorDerivative(TrialStress, *ptr_material_parameter, 2);
                 err += Tensor2VectorSysR2(dFoda2, V6);
                 for (i=0; i<6; i++)
                   N19(12+Num_internal_scalar+i) = V6(i);
             }

           // R & CI
             if (iter_counter == 0)
               {
                 R19.Zero();
                 CI.Zero();
   // BJ optitmize!!!!
                 for (i=0; i<numMV; i++)
                   CI(i, i) = 1.0;
               }

           // d2_lambda
           upper = 0.0;
           lower = 0.0;

   // Mahdi Taiebat & Boris Jeremic 17April2007
   // MAHDI GOT THIS ONE it is a bug 17April2007
   //        if (iter_counter == 0)
   //          {
   //            for (i=0; i<numMV; i++)
   //              lower += N19(i) * M19(i);
   //          }

           T19.addMatrixVector(0.0, CI, R19, 1.0);
           for (i=0; i<numMV; i++)
             upper += N19(i) * T19(i);   // N19 -> [n_mn ; \ksi_B]^T
                                         // T19 -> CI * R19
                                         // R19 -> vector of residuals

           T19.addMatrixVector(0.0, CI, M19, 1.0);
           for (i=0; i<numMV; i++)
             lower += N19(i) * T19(i);

           d2_lambda = (YieldFun - upper) / lower;

        // update:
           iter_counter++;

           Delta_lambda += d2_lambda;

           //if ( Delta_lambda < 0.0 )
        //  Delta_lambda = 0.0;

           T19 = R19 + M19*d2_lambda;  // T19 -> R^old + \d \Delta  \lambda * M19
                                       // M19 -> E_ijkl * m_kl and it is a vector
           D19.addMatrixVector(0.0, CI, T19, -1.0);

           err += Vector2TensorSysR2(D19, Dstress, 0);

# ifdef _TEACHING_MODE
//           TrialStress.report(" NewTemplate3Dep::Implicit(   TrialStress");
#endif
           TrialStress += Dstress;

   //        Dstress.report("Dstress");
   //        TrialStress.report("TrialStress");


           if (Num_internal_scalar == 1)
             {

               Dq = D19(6);
               // Mahdi Taiebat & Boris Jeremic 20April2007
               // redefining S with double causes a problem here
               // double S = 0.0; has been added in the beginning of ElastoPlastic part...

               // double S = ptr_material_parameter->getInternal_Scalar(0);
               S = ptr_material_parameter->getInternal_Scalar(0);
               err += ptr_material_parameter->setInternal_Scalar(0, S + Dq );
             }

           if (Num_internal_tensor >= 1)
             {
               err += Vector2TensorSysR2(D19, Da, 6+Num_internal_scalar);
               T = ptr_material_parameter->getInternal_Tensor(0);
               err += ptr_material_parameter->setInternal_Tensor(0, T + Da );
             }

           if (Num_internal_tensor == 2)
             {
               err += Vector2TensorSysR2(D19, Da2, 12+Num_internal_scalar);
               T2 = ptr_material_parameter->getInternal_Tensor(1);
               err += ptr_material_parameter->setInternal_Tensor(1, T2 + Da2 );
             }


//	Ce.print("Ce","\n Ce");
           DPlastic_Strain = Ce("ijkl")*Dstress("kl");
             DPlastic_Strain.null_indices();

//           TrialPlastic_Strain -= DPlastic_Strain;

           YieldFun = ptr_yield_function->YieldFunctionValue(TrialStress, *ptr_material_parameter);
           //cout << "YieldFun = " << YieldFun << endln;

           YieldFun2 = YieldFun;
           // If not convergent, for Line Search Algorithm
           if ( (fabs(YieldFun2) > fabs(YieldFun1) || iter_counter == ITMAX) && this->caseIndex == 2 )
             {
               this->divergeOrnot = 1;

               cout << "Error, NewTemplate3Dep::Backward Agorithm, this->divergeOrnot = 1 !" << endln;

               return 0;
             }
           //
           YieldFun1 = YieldFun;

           // R
           R19.Zero();


       Rstress = TrialStress - PredictedStress + (m *Delta_lambda);



// Nima
// Nima
//	dQods.print("dQods" , "dQods");
//	cout << "\n Delta_lambda = " << Delta_lambda <<endl;

	incr_Pstrain = dQods * Delta_lambda;
	TrialPlastic_Strain = start_Pstrain + incr_Pstrain;

//	cout << "\n Delta_lambda = " << Delta_lambda <<endl;
//	m.print("m","\n m");
//	TrialPlastic_Strain.print("m","\n m");
//	Rstress.print("Rstress","Rstress");


   //     TrialStress.print("TrialStress","\n TrialStress on line 1334");
   //     PredictedStress.print("PredictedStress","\n PredictedStress on line 1334");
   //     m.print("m","\n m on line 1334");
   //     Rstress.print("Rstress","Rstress on line 1334");

           err += Tensor2VectorSysR2(Rstress, V6);
           for (i=0; i<6; i++)
             R19(i) = V6(i);

           if (Num_internal_scalar == 1)
             {
             // Mahdi Taiebat & Boris Jeremic 20April2007
             // redefining S with double causes a problem here
             // double S = 0.0; has been added in the beginning of ElastoPlastic part...

             // double S = ptr_material_parameter->getInternal_Scalar(0);
               S = ptr_material_parameter->getInternal_Scalar(0);
               Rq = S - start_InScalar - (hi *Delta_lambda);
               R19(6) = Rq;
             }

           if (Num_internal_tensor >= 1)
             {
               T = ptr_material_parameter->getInternal_Tensor(0);
               Ra = T - start_InTensor - (hk *Delta_lambda);
               err += Tensor2VectorSysR2(Ra, V6);
               for (i=0; i<6; i++)
                 R19(6+Num_internal_scalar+i) = V6(i);
             }

           if (Num_internal_tensor == 2)
             {
               T2 = ptr_material_parameter->getInternal_Tensor(1);
               Ra2 = T2 - start_InTensor2 - (hk2 *Delta_lambda);
               err += Tensor2VectorSysR2(Ra2, V6);
               for (i=0; i<6; i++)
                 R19(12+Num_internal_scalar+i) = V6(i);
             }

             RNorm = R19.Norm();
           //cout << "RNorm = " << RNorm << endln;

           // CC and CI
           CC.Zero();

           dmods = ptr_plastic_flow->Dm_Ds(TrialStress, TrialStrain, *ptr_material_parameter);
           dmods =  Ee("ijpq") * dmods("pqmn");
             dmods.null_indices();

           // if 0 internal scalar, and 0 internal tensor
           err += Tensor2MatrixSysR4(dmods, M66);
           for (i=0; i<6; i++) {
          for (j=0; j<6; j++) {
               CC(i, j) = M66(i, j);
             }
           }

           // if 1 internal scalar, and 0 internal tensor
           if (Num_internal_scalar == 1 && Num_internal_tensor == 0)
           {

             tensor dmodq(2, def_dim_2, 0.0);

             tensor dhiods(2, def_dim_2, 0.0);
             double dhiodq = 0.0;

             dmodq = ptr_plastic_flow->Dm_Diso(TrialStress, TrialStrain, *ptr_material_parameter);
             dmodq =  Ee("ijpq") * dmodq("pq");
               dmodq.null_indices();

             err += Tensor2VectorSysR2(dmodq, V6);
             for (i=0; i<6; i++)
               CC(i, 6) = V6(i);

             dhiods = ptr_scalar_evolution[0]->DH_Ds(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhiodq = ptr_scalar_evolution[0]->DH_Diso(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);

             err += Tensor2VectorSysR2(dhiods, V6);
             for (i=0; i<6; i++)
               CC(6, i) = -V6(i);

             CC(6, 6) = -dhiodq;
           }

           // if 1 internal scalar, and 1 internal tensor
           if (Num_internal_scalar == 1 && Num_internal_tensor == 1)
           {
             tensor dmodq(2, def_dim_2, 0.0);
             tensor dmoda(4, def_dim_4, 0.0);

             tensor dhiods(2, def_dim_2, 0.0);
             double dhiodq = 0.0;
             tensor dhioda(2, def_dim_2, 0.0);

             tensor dhkods(4, def_dim_4, 0.0);
             tensor dhkodq(2, def_dim_2, 0.0);
             tensor dhkoda(4, def_dim_4, 0.0);

             dmodq = ptr_plastic_flow->Dm_Diso(TrialStress, TrialStrain, *ptr_material_parameter);
             dmoda = ptr_plastic_flow->Dm_Dkin(TrialStress, TrialStrain, *ptr_material_parameter);
             dmodq =  Ee("ijpq") * dmodq("pq");
               dmodq.null_indices();
             dmoda =  Ee("ijpq") * dmoda("pqmn");
               dmoda.null_indices();

             err += Tensor2VectorSysR2(dmodq, V6);
             for (i=0; i<6; i++)
               CC(i, 6) = V6(i);

             err += Tensor2MatrixSysR4(dmoda, M66);
             for (i=0; i<6; i++)
               {
                 for (j=0; j<6; j++)
                   {
                     CC(i, 7+j) = M66(i,j);
                   }
                }

             dhiods = ptr_scalar_evolution[0]->DH_Ds(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhiodq = ptr_scalar_evolution[0]->DH_Diso(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhioda = ptr_scalar_evolution[0]->DH_Dkin(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);

             err += Tensor2VectorSysR2(dhiods, V6);
             for (i=0; i<6; i++)
               CC(6, i) = -V6(i);

             CC(6, 6) = -dhiodq;

             err += Tensor2VectorSysR2(dhioda, V6);
             for (i=0; i<6; i++)
               CC(6, 7+i) = -V6(i);

             dhkods = ptr_tensor_evolution[0]->DHij_Ds(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhkodq = ptr_tensor_evolution[0]->DHij_Diso(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhkoda = ptr_tensor_evolution[0]->DHij_Dkin(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);

             err += Tensor2MatrixSysR4(dhkods, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(7+i, j) = -M66(i,j);
               }
             }

             err += Tensor2VectorSysR2(dhkodq, V6);
             for (i=0; i<6; i++)
               CC(7+i, 6) = -V6(i);

             err += Tensor2MatrixSysR4(dhkoda, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(7+i, 7+j) = -M66(i,j);
               }
             }

           }

           // if 1 internal scalar, and 2 internal tensor
           if (Num_internal_scalar == 1 && Num_internal_tensor == 2)
           {
             tensor dmodq(2, def_dim_2, 0.0);
             tensor dmoda(4, def_dim_4, 0.0);
             tensor dmoda2(4, def_dim_4, 0.0);

             tensor dhiods(2, def_dim_2, 0.0);
             double dhiodq = 0.0;
             tensor dhioda(2, def_dim_2, 0.0);
             tensor dhioda2(2, def_dim_2, 0.0);

             tensor dhkods(4, def_dim_4, 0.0);
             tensor dhkodq(2, def_dim_2, 0.0);
             tensor dhkoda(4, def_dim_4, 0.0);
             tensor dhkoda2(4, def_dim_4, 0.0);

             tensor dhk2ods(4, def_dim_4, 0.0);
             tensor dhk2odq(2, def_dim_2, 0.0);
             tensor dhk2oda(4, def_dim_4, 0.0);
             tensor dhk2oda2(4, def_dim_4, 0.0);

             dmodq = ptr_plastic_flow->Dm_Diso(TrialStress, TrialStrain, *ptr_material_parameter);
             dmoda = ptr_plastic_flow->Dm_Dkin(TrialStress, TrialStrain, *ptr_material_parameter);
             dmoda2 = ptr_plastic_flow->Dm_Dkin2(TrialStress, TrialStrain, *ptr_material_parameter);
             dmodq =  Ee("ijpq") * dmodq("pq");
               dmodq.null_indices();
             dmoda =  Ee("ijpq") * dmoda("pqmn");
               dmoda.null_indices();
             dmoda2 =  Ee("ijpq") * dmoda2("pqmn");
               dmoda2.null_indices();

             err += Tensor2VectorSysR2(dmodq, V6);
             for (i=0; i<6; i++)
               CC(i, 6) = V6(i);
               // CC(0,6)=V6(0);
               // CC(1,6)=V6(1);
               // CC(2,6)=V6(2);
               // CC(3,6)=V6(3);
               // CC(4,6)=V6(4);
               // CC(5,6)=V6(5);

             err += Tensor2MatrixSysR4(dmoda, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(i, 7+j) = M66(i,j);
               }
             }

             err += Tensor2MatrixSysR4(dmoda2, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(i, 13+j) = M66(i,j);
               }
             }

             dhiods = ptr_scalar_evolution[0]->DH_Ds(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhiodq = ptr_scalar_evolution[0]->DH_Diso(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhioda = ptr_scalar_evolution[0]->DH_Dkin(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhioda2 = ptr_scalar_evolution[0]->DH_Dkin2(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);

             err += Tensor2VectorSysR2(dhiods, V6);
             for (i=0; i<6; i++)
               CC(6, i) = -V6(i);

             CC(6, 6) = -dhiodq;

             err += Tensor2VectorSysR2(dhioda, V6);
             for (i=0; i<6; i++)
               CC(6, 7+i) = -V6(i);

             err += Tensor2VectorSysR2(dhioda2, V6);
             for (i=0; i<6; i++)
               CC(6, 13+i) = -V6(i);

             dhkods = ptr_tensor_evolution[0]->DHij_Ds(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhkodq = ptr_tensor_evolution[0]->DHij_Diso(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhkoda = ptr_tensor_evolution[0]->DHij_Dkin(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhkoda2 = ptr_tensor_evolution[0]->DHij_Dkin2(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);

             err += Tensor2MatrixSysR4(dhkods, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(7+i, j) = -M66(i,j);
               }
             }

             err += Tensor2VectorSysR2(dhkodq, V6);
             for (i=0; i<6; i++)
               CC(7+i, 6) = -V6(i);

             err += Tensor2MatrixSysR4(dhkoda, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(7+i, 7+j) = -M66(i,j);
               }
             }

             err += Tensor2MatrixSysR4(dhkoda2, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(7+i, 13+j) = -M66(i,j);
               }
             }

             dhk2ods = ptr_tensor_evolution[1]->DHij_Ds(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhk2odq = ptr_tensor_evolution[1]->DHij_Diso(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhk2oda = ptr_tensor_evolution[1]->DHij_Dkin(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhk2oda2 = ptr_tensor_evolution[1]->DHij_Dkin2(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);

             err += Tensor2MatrixSysR4(dhk2ods, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(13+i, j) = -M66(i,j);
               }
             }

             err += Tensor2VectorSysR2(dhk2odq, V6);
             for (i=0; i<6; i++)
               CC(13+i, 6) = -V6(i);

             err += Tensor2MatrixSysR4(dhk2oda, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(13+i, 7+j) = -M66(i,j);
               }
             }

             err += Tensor2MatrixSysR4(dhk2oda2, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(13+i, 13+j) = -M66(i,j);
               }
             }

           }

           // if 0 internal scalar, and 1 internal tensor
           if (Num_internal_scalar == 0 && Num_internal_tensor == 1)
           {
             tensor dmoda(4, def_dim_4, 0.0);

             tensor dhkods(4, def_dim_4, 0.0);
             tensor dhkoda(4, def_dim_4, 0.0);

             dmoda = ptr_plastic_flow->Dm_Dkin(TrialStress, TrialStrain, *ptr_material_parameter);
             dmoda =  Ee("ijpq") * dmoda("pqmn");
               dmoda.null_indices();

             err += Tensor2MatrixSysR4(dmoda, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(i, 6+j) = M66(i,j);
               }
             }

             dhkods = ptr_tensor_evolution[0]->DHij_Ds(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhkoda = ptr_tensor_evolution[0]->DHij_Dkin(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);

             err += Tensor2MatrixSysR4(dhkods, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(6+i, j) = -M66(i,j);
               }
             }

             err += Tensor2MatrixSysR4(dhkoda, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(6+i, 6+j) = -M66(i,j);
               }
             }

           }

           // if 0 internal scalar, and 2 internal tensor
           if (Num_internal_scalar == 0 && Num_internal_tensor == 2)
           {
             tensor dmoda(4, def_dim_4, 0.0);
             tensor dmoda2(4, def_dim_4, 0.0);

             tensor dhkods(4, def_dim_4, 0.0);
             tensor dhkoda(4, def_dim_4, 0.0);
             tensor dhkoda2(4, def_dim_4, 0.0);

             tensor dhk2ods(4, def_dim_4, 0.0);
             tensor dhk2oda(4, def_dim_4, 0.0);
             tensor dhk2oda2(4, def_dim_4, 0.0);

             dmoda = ptr_plastic_flow->Dm_Dkin(TrialStress, TrialStrain, *ptr_material_parameter);
             dmoda2 = ptr_plastic_flow->Dm_Dkin2(TrialStress, TrialStrain, *ptr_material_parameter);
             dmoda =  Ee("ijpq") * dmoda("pqmn");
               dmoda.null_indices();
             dmoda2 =  Ee("ijpq") * dmoda2("pqmn");
               dmoda2.null_indices();

             err += Tensor2MatrixSysR4(dmoda, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(i, 6+j) = M66(i,j);
               }
             }

             err += Tensor2MatrixSysR4(dmoda2, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(i, 12+j) = M66(i,j);
               }
             }

             dhkods = ptr_tensor_evolution[0]->DHij_Ds(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhkoda = ptr_tensor_evolution[0]->DHij_Dkin(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhkoda2 = ptr_tensor_evolution[0]->DHij_Dkin2(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);

             err += Tensor2MatrixSysR4(dhkods, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(6+i, j) = -M66(i,j);
               }
             }

             err += Tensor2MatrixSysR4(dhkoda, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(6+i, 6+j) = -M66(i,j);
               }
             }

             err += Tensor2MatrixSysR4(dhkoda2, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(6+i, 12+j) = -M66(i,j);
               }
             }

             dhk2ods = ptr_tensor_evolution[1]->DHij_Ds(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhk2oda = ptr_tensor_evolution[1]->DHij_Dkin(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
             dhk2oda2 = ptr_tensor_evolution[1]->DHij_Dkin2(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);

             err += Tensor2MatrixSysR4(dhk2ods, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(12+i, j) = -M66(i,j);
               }
             }

             err += Tensor2MatrixSysR4(dhk2oda, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(12+i, 6+j) = -M66(i,j);
               }
             }

             err += Tensor2MatrixSysR4(dhk2oda2, M66);
             for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                 CC(12+i, 12+j) = -M66(i,j);
               }
             }

           }

           CC *= Delta_lambda;

           for (i=0; i<numMV; i++)
             CC(i, i) = 1.0 + CC(i, i);

           err += CC.Invert(CI);
         //} while ( fabs(YieldFun) > YieldFun_relative && iter_counter < ITMAX);







          }
        while ( (fabs(YieldFun) > YieldFun_relative || RNorm > TOL) && iter_counter < ITMAX);




      // Mahdi Taiebat 22April2007
      //cout <<  "********** iter_counter = " << iter_counter << "\n";





      // ################## End of do-while ########################

      // To get Algorithmic Tangent Stiffness

      // Updated M & N, beginning

      // M
      M19.Zero();

      m = ptr_plastic_flow->PlasticFlowTensor(TrialStress, TrialStrain, *ptr_material_parameter);
      m =  Ee("ijpq") * m("pq");
        m.null_indices();
      err += Tensor2VectorSysR2(m, V6);
      for (i=0; i<6; i++)
        M19(i) = V6(i);

   if (Num_internal_scalar == 1)
     {
       hi = ptr_scalar_evolution[0]->H(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
          M19(6) = -hi;
     }

   if (Num_internal_tensor >= 1)
     {
       hk = ptr_tensor_evolution[0]->Hij(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
       err += Tensor2VectorSysR2(hk, V6);
       for (i=0; i<6; i++)
         M19(6+Num_internal_scalar+i) = -V6(i);
     }

   if (Num_internal_tensor == 2)
     {
       hk2 = ptr_tensor_evolution[1]->Hij(*ptr_plastic_flow, TrialStress, TrialStrain, *ptr_material_parameter);
       err += Tensor2VectorSysR2(hk2, V6);
       for (i=0; i<6; i++)
         M19(12+Num_internal_scalar+i) = -V6(i);
     }

      // N
      N19.Zero();

      dFods = ptr_yield_function->StressDerivative(TrialStress, *ptr_material_parameter);
      err += Tensor2VectorSysR2(dFods, V6);

      for (i=0; i<6; i++)
        N19(i) = V6(i);

      if (Num_internal_scalar_YF == 1)
        {
          dFodq = ptr_yield_function->InScalarDerivative(TrialStress, *ptr_material_parameter, 1);
          N19(6) = dFodq;
        }

      if (Num_internal_tensor_YF == 1)
        {
          dFoda = ptr_yield_function->InTensorDerivative(TrialStress, *ptr_material_parameter, 1);
          err += Tensor2VectorSysR2(dFoda, V6);
          for (i=0; i<6; i++)
            N19(6+Num_internal_scalar+i) = V6(i);
        }

      else if (Num_internal_tensor_YF == 2)
        {
          dFoda = ptr_yield_function->InTensorDerivative(TrialStress, *ptr_material_parameter, 1);
          err += Tensor2VectorSysR2(dFoda, V6);
          for (i=0; i<6; i++)
            N19(6+Num_internal_scalar+i) = V6(i);

          dFoda2 = ptr_yield_function->InTensorDerivative(TrialStress, *ptr_material_parameter, 2);
          err += Tensor2VectorSysR2(dFoda2, V6);
          for (i=0; i<6; i++)
            N19(12+Num_internal_scalar+i) = V6(i);
        }

      // Updated M & N, ending

      if (iStep == NumStep_in)
        {
          T19.Zero();
          R19.Zero();
          for (i=0; i<numMV; i++)
            {
              for (j=0; j<numMV; j++)
                {
                  T19(i) += CI(j, i) * N19(j);
                  R19(i) += CI(i, j) * M19(j);
                }
            }

          lower = 0.0;
          for (i=0; i<numMV; i++)
            lower += N19(i) * R19(i);

          M66.Zero();
          for (i=0; i<6; i++)
            {
              for (j=0; j<6; j++)
                {
                  M66(i, j) = CI(i, j) - R19(i) * T19(j) / lower;
                }
            }

          err += Matrix2TensorSysR4(M66, Ttemp);
          Stiffness = Ttemp("ijkl") * Ee("klmn");
          Stiffness.null_indices();
        }
// ELASTIC-PLASTIC  ELASTIC-PLASTIC  ELASTIC-PLASTIC  ELASTIC-PLASTIC  ELASTIC-PLASTIC  ELASTIC-PLASTIC
    }

  }


    return err;
}


//================================================================================
int NewTemplate3Dep::ImplicitLineSearch(const straintensor& strain_incr)
{
    int err = 0;

    int Num_internal_scalar = 0;
    int Num_internal_tensor = 0;

    straintensor s_stress;
    stresstensor s_strain;
    straintensor s_Pstrain;

    s_stress = getStressTensor();
    s_strain = getStrainTensor();
    s_Pstrain = getPlasticStrainTensor();

    double s_InScalar = 0.0;
    tensor s_InTensor(2, def_dim_2, 0.0);
    tensor s_InTensor2(2, def_dim_2, 0.0);

    Num_internal_scalar = ptr_material_parameter->getNum_Internal_Scalar();
    if (Num_internal_scalar > 1) {
      cout << "Error, NewTemplate3Dep::Backward Agorithm, scalar internal variables more than 1!" << endln;
      exit(1);
    }
    else if (Num_internal_scalar == 1) {
   s_InScalar = ptr_material_parameter->getInternal_Scalar(0);
    }

    Num_internal_tensor = ptr_material_parameter->getNum_Internal_Tensor();
    if (Num_internal_tensor > 2) {
      cout << "Error, NewTemplate3Dep::Backward Agorithm, scalar internal variables more than 1!" << endln;
      exit(1);
    }
    else if (Num_internal_tensor == 1) {
   s_InTensor = ptr_material_parameter->getInternal_Tensor(0);
 }
    else if (Num_internal_tensor == 2) {
      s_InTensor = ptr_material_parameter->getInternal_Tensor(0);
      s_InTensor2 = ptr_material_parameter->getInternal_Tensor(1);
    }

    this->subStep = 1;  // To avoid confusing
    int num_step = 1;

    // Do-While Loop Begins
    do
      {
        err = Implicit(strain_incr, num_step);

        if (divergeOrnot == 1)
          {

            // Return to the start state:
            TrialStress.Initialize(s_stress);
            TrialStrain.Initialize(s_strain);
            TrialPlastic_Strain.Initialize(s_Pstrain);
            if (Num_internal_scalar == 1)
              err += ptr_material_parameter->setInternal_Scalar(0, s_InScalar );
            if (Num_internal_tensor == 1)
              err += ptr_material_parameter->setInternal_Tensor(0, s_InTensor );
            else if (Num_internal_tensor == 2)
              {
                err += ptr_material_parameter->setInternal_Tensor(0, s_InTensor );
                err += ptr_material_parameter->setInternal_Tensor(1, s_InTensor2 );
              }

         // Half the step:
            num_step *= 2;
          }

      } while (divergeOrnot == 1);
    // Do-While Loop Ends

    return err;

}





// Mahdi's idea
//================================================================================
int NewTemplate3Dep::ScaledExplicit(const straintensor& strain_incr, int NumStep_in)
{
    int err = 0;
    int i = 0;
    int iStep = 0;
    double f_start = 0.0;
    double f_pred  = 0.0;
// scaling idea
    double f_final  = 0.0;

    BJtensor Ee(4, def_dim_4, 0.0);
    BJtensor Ep(4, def_dim_4, 0.0);

    straintensor intersection_strain;
    stresstensor intersection_stress;
    double intersection_factor = 0.0;

    int Num_internal_scalar = 0;
    int Num_internal_tensor = 0;
    int Num_internal_scalar_YF = 0;
    int Num_internal_tensor_YF = 0;

    double upper = 0.0;
    double lower = 0.0;
    double Delta_lambda = 0.0;
    double hardMod  = 0.0;
    double h_s = 0.0;
    double xi_s = 0.0;
    stresstensor dFods;
    straintensor dQods;
    BJtensor Hq(2, def_dim_2, 0.0);
    BJtensor Hf(2, def_dim_2, 0.0);

    stresstensor h_t;
    stresstensor xi_t;

    straintensor incr_strain;
    stresstensor incr_stress;
    straintensor incr_Pstrain;
    stresstensor ep_stress;
    stresstensor predicted_stress;

    straintensor start_stress;
    stresstensor start_strain;
    straintensor start_Pstrain;

// scaling, new part
    stresstensor scaling_of_explicit_stress;

  for (iStep = 1; iStep <= NumStep_in; iStep++) {
    start_stress = getStressTensor();
    start_strain = getStrainTensor();
    start_Pstrain = getPlasticStrainTensor();

    intersection_stress.Initialize(start_stress);

    err += ptr_elastic_state->setStress(start_stress);
    err += ptr_elastic_state->setStrain(start_strain);
    Ee = ptr_elastic_state->getElasticStiffness(*ptr_material_parameter);

    incr_strain = strain_incr *(1.0/(double)(NumStep_in));
    incr_stress = Ee("ijpq") * incr_strain("pq");
    incr_stress.null_indices();

    predicted_stress = start_stress + incr_stress;


/////////////////////////////////////////////////////
//       start_stress.print("start_stress" , "start_stress");
//       incr_stress.print("incr_stress" , "incr_stress");
//       predicted_stress.print("predicted_stress" , "predicted_stress");
/////////////////////////////////////////////////////


    f_start = ptr_yield_function->YieldFunctionValue( start_stress, *ptr_material_parameter );
    f_pred =  ptr_yield_function->YieldFunctionValue( predicted_stress, *ptr_material_parameter );


# ifdef _TEACHING_MODE
      fprintf(stdout,"explicit f_start = %12.4e",f_start);
      fprintf(stdout,"explicit f_pred = %12.4e",f_pred);
# endif


    // If Elastic
    if ( (f_start < 0.0 && f_pred <= FTOL) || f_start > f_pred ) 
      {
        TrialStrain = start_strain + incr_strain;
        TrialStress.Initialize(predicted_stress);

        if (iStep == NumStep_in)
          Stiffness = Ee;
      }
    else 
      {
        // If Elastic and then Elastic-Plastic
        if ( f_start < 0.0 )  
          {
            intersection_factor = zbrentstress( start_stress, predicted_stress, 0.0, 1.0, TOL );
            intersection_stress = yield_surface_cross( start_stress, predicted_stress, intersection_factor );
            intersection_strain = start_strain + (incr_strain * intersection_factor);
            
            incr_stress = predicted_stress - intersection_stress;  // necessary
            
            err += ptr_elastic_state->setStress(intersection_stress);
            err += ptr_elastic_state->setStrain(intersection_strain);
            Ee = ptr_elastic_state->getElasticStiffness(*ptr_material_parameter);
          }

      // If E-P Response,
      Delta_lambda = 0.0;

      dFods = ptr_yield_function->StressDerivative( intersection_stress, *ptr_material_parameter );
      dQods = ptr_plastic_flow->PlasticFlowTensor( intersection_stress, intersection_strain, *ptr_material_parameter );

      // E_ijkl * R_kl
      Hq = Ee("ijkl") * dQods("kl");
        Hq.null_indices();

      // L_ij * E_ijkl
      Hf = dFods("ij") * Ee("ijkl");
        Hf.null_indices();

      // L_ij * E_ijkl * d e_kl ( true EP strain increment)
      upper = ( dFods("ij") * incr_stress("ij") ).trace();

      // L_ij * E_ijkl * R_kl
      lower = ( Hf("ij") * dQods("ij") ).trace();

      hardMod  = 0.0;

      // Evolution of scalar (isotropic) internal variables in yield function
      Num_internal_scalar_YF = ptr_yield_function->getNumInternalScalar();
      for (i = 0; i < Num_internal_scalar_YF; i++) 
        {
          h_s = ptr_scalar_evolution[i]->H(*ptr_plastic_flow, intersection_stress, intersection_strain, *ptr_material_parameter);
          xi_s = ptr_yield_function->InScalarDerivative( intersection_stress, *ptr_material_parameter, i+1);
          hardMod += h_s * xi_s;
        }

      // Evolution of tensor (kinematic) internal variables in yield function
      Num_internal_tensor_YF = ptr_yield_function->getNumInternalTensor();
      for (i = 0; i < Num_internal_tensor_YF; i++) 
        {
          h_t = ptr_tensor_evolution[i]->Hij(*ptr_plastic_flow, intersection_stress, intersection_strain, *ptr_material_parameter);
          xi_t = ptr_yield_function->InTensorDerivative( intersection_stress, *ptr_material_parameter, i+1);
          hardMod += ( h_t("mn") * xi_t("mn") ).trace();
        }

      lower -= hardMod;

      Delta_lambda = upper / lower;

      if (Delta_lambda < 0.0)
        Delta_lambda = 0.0;

      // Plastic strain increment
      incr_Pstrain = dQods * Delta_lambda;
      ep_stress = predicted_stress - (Hq *Delta_lambda);


      TrialStress.Initialize(ep_stress);
      TrialStrain = start_strain + incr_strain;
      TrialPlastic_Strain = start_Pstrain + incr_Pstrain;

      // Update internal scalar variables
      Num_internal_scalar = ptr_material_parameter->getNum_Internal_Scalar();
      for (i = 0; i < Num_internal_scalar; i++) 
      {
        double dS = ( ptr_scalar_evolution[i]->H(*ptr_plastic_flow, intersection_stress, intersection_strain, *ptr_material_parameter) ) *Delta_lambda;
        double S = ptr_material_parameter->getInternal_Scalar(i);
        err += ptr_material_parameter->setInternal_Scalar(i, S + dS );
      }

      // Update internal tensor variables
      Num_internal_tensor = ptr_material_parameter->getNum_Internal_Tensor();
      for (i = 0; i < Num_internal_tensor; i++) 
      {
        stresstensor dT = ptr_tensor_evolution[i]->Hij(*ptr_plastic_flow, intersection_stress, intersection_strain, *ptr_material_parameter) *Delta_lambda;
        stresstensor T = ptr_material_parameter->getInternal_Tensor(i);
        err += ptr_material_parameter->setInternal_Tensor(i, T + dT );
      }


// apply scaling down to yield surface according to Mahdi's idea (reference:...)
      lower = (dFods("ij") * dFods("ij")).trace();
      f_final =  ptr_yield_function->YieldFunctionValue( ep_stress, *ptr_material_parameter ); 
      scaling_of_explicit_stress =  dFods("ij") * (f_final/lower);

      ep_stress = ep_stress - scaling_of_explicit_stress;
      TrialStress.Initialize(ep_stress);
// apply scaling down to yield surface according to Mahdi's idea (reference:...)



      // To obtain Eep, at the last step
      if (iStep == NumStep_in) 
      {
        Ep = Hq("pq") * Hf("mn");
          Ep.null_indices();
        Ep = Ep * (1.0/lower);

        if ( Delta_lambda > 0.0 )
          Stiffness = Ee - Ep;
        else
          Stiffness = Ee;
      }

    }
// end of else loop for el-pl



  }

  return err;
}







// Trying to find intersection point according to M. Crisfield's book
// "Non-linear Finite Element Analysis of Solids and Structures "  Chp 6.6.1 pp168.
//================================================================================
stresstensor NewTemplate3Dep::yield_surface_cross(const stresstensor & start_stress,
                                                  const stresstensor & end_stress, double a)
{
    stresstensor delta_stress = end_stress - start_stress;
    stresstensor intersection_stress = start_stress + delta_stress * a;

    return intersection_stress;
}

// Routine used by yield_surface_cross to find the stresstensor at cross point
//================================================================================
double NewTemplate3Dep::zbrentstress(const stresstensor& start_stress,
                                  const stresstensor& end_stress,
                                  double x1, double x2, double tol) const
{
  double EPS = d_macheps();

  int iter;
  double a = x1;
  double b = x2;
  double c = 0.0;
  double d = 0.0;
  double e = 0.0;
  double min1 = 0.0;
  double min2 = 0.0;
  double fc = 0.0;
  double p = 0.0;
  double q = 0.0;
  double r = 0.0;
  double s = 0.0;
  double tol1 = 0.0;
  double xm = 0.0;

  double fa = func(start_stress, end_stress, *ptr_material_parameter, a);
  double fb = func(start_stress, end_stress, *ptr_material_parameter, b);

  if ( (fb * fa) > 0.0) {
      cout << "\a\n Root must be bracketed in ZBRENTstress " << endln;
      exit(1);
  }

  fc = fb;
  for ( iter = 1; iter <= ISMAX; iter++ ) {
      if ( (fb * fc) > 0.0) {
          c = a;
          fc = fa;
          e = d = b - a;
      }
      if ( fabs(fc) < fabs(fb) ) {
          a = b;
          b = c;
          c = a;
          fa = fb;
          fb = fc;
          fc = fa;
      }
      tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;
      xm = 0.5 * (c - b);
      if ( fabs(xm) <= tol1 || fb == 0.0 )
        return b;

      if ( fabs(e) >= tol1 && fabs(fa) > fabs(fb) ) {
          s = fb / fa;
          if (a == c) {
              p = 2.0 * xm * s;
              q = 1.0 - s;
          }
          else {
              q = fa / fc;
              r = fb / fc;
              p = s * ( 2.0 * xm * q * (q - r) - (b - a) * (r - 1.0) );
              q = (q - 1.0) * (r - 1.0) * (s - 1.0);
          }
          if (p > 0.0)
            q = -q;
          p = fabs(p);
          min1 = 3.0 * xm * q - fabs(tol1*q);
          min2 = fabs(e*q);
          if (2.0*p < (min1 < min2 ? min1 : min2)) {
              e = d;
              d = p/q;
          }
          else {
              d = xm;
              e = d;
          }
        }
        else {
          d = xm;
          e = d;
        }
      a = b;
      fa = fb;
      if (fabs(d) > tol1)
        b += d;
      else
        b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
      fb = func(start_stress, end_stress, *ptr_material_parameter, b);
  }

  return 0.0;
}

//================================================================================
double NewTemplate3Dep::func(const stresstensor& start_stress,
                          const stresstensor& end_stress,
                          const MaterialParameter& ptr_material_parameter,
                          double alfa ) const
{
    stresstensor alfa_stress = ( start_stress * (1.0 - alfa) ) + ( end_stress * alfa );

    double f = ptr_yield_function->YieldFunctionValue( alfa_stress, ptr_material_parameter );

    return f;
}

//================================================================================
int NewTemplate3Dep::Tensor2MatrixSysR4(const tensor& T, Matrix& M)
{
  int rank = T.rank();
  if (rank != 4) {
    cout << "NewTemplate3Dep::Tensor2MatrixSysR4 - tensor must be of rank 4" << endln;
    return 1;
  }

  int nr = M.noRows();
  int nc = M.noCols();
  if (nr != 6 || nc != 6) {
    cout << "NewTemplate3Dep::Tensor2MatrixSysR4 - matrix must be of (6, 6)" << endln;;
    return 1;
  }

  double sqrt2 = sqrt(2.0);
  double two = 2.0;

  // Adopt method from Helnwein (2001):

  M(0,0) = T.cval(1,1,1,1);
  M(0,1) = T.cval(1,1,2,2);
  M(0,2) = T.cval(1,1,3,3);
  M(0,3) = T.cval(1,1,1,2) *sqrt2;
  M(0,4) = T.cval(1,1,2,3) *sqrt2;
  M(0,5) = T.cval(1,1,1,3) *sqrt2;

  M(1,0) = T.cval(2,2,1,1);
  M(1,1) = T.cval(2,2,2,2);
  M(1,2) = T.cval(2,2,3,3);
  M(1,3) = T.cval(2,2,1,2) *sqrt2;
  M(1,4) = T.cval(2,2,2,3) *sqrt2;
  M(1,5) = T.cval(2,2,1,3) *sqrt2;

  M(2,0) = T.cval(3,3,1,1);
  M(2,1) = T.cval(3,3,2,2);
  M(2,2) = T.cval(3,3,3,3);
  M(2,3) = T.cval(3,3,1,2) *sqrt2;
  M(2,4) = T.cval(3,3,2,3) *sqrt2;
  M(2,5) = T.cval(3,3,1,3) *sqrt2;

  M(3,0) = T.cval(1,2,1,1) *sqrt2;
  M(3,1) = T.cval(1,2,2,2) *sqrt2;
  M(3,2) = T.cval(1,2,3,3) *sqrt2;
  M(3,3) = T.cval(1,2,1,2) *two;
  M(3,4) = T.cval(1,2,2,3) *two;
  M(3,5) = T.cval(1,2,1,3) *two;

  M(4,0) = T.cval(2,3,1,1) *sqrt2;
  M(4,1) = T.cval(2,3,2,2) *sqrt2;
  M(4,2) = T.cval(2,3,3,3) *sqrt2;
  M(4,3) = T.cval(2,3,1,2) *two;
  M(4,4) = T.cval(2,3,2,3) *two;
  M(4,5) = T.cval(2,3,1,3) *two;

  M(5,0) = T.cval(1,3,1,1) *sqrt2;
  M(5,1) = T.cval(1,3,2,2) *sqrt2;
  M(5,2) = T.cval(1,3,3,3) *sqrt2;
  M(5,3) = T.cval(1,3,1,2) *two;
  M(5,4) = T.cval(1,3,2,3) *two;
  M(5,5) = T.cval(1,3,1,3) *two;

    return 0;
}

//================================================================================
int NewTemplate3Dep::Tensor2VectorSysR2(const tensor& T, Vector& V)
{
  int rank = T.rank();
  if (rank != 2) {
    cout << "NewTemplate3Dep::Tensor2VectorSysR2 - tensor must be of rank 4" << endln;
    return 1;
  }

  int nsize = V.Size();
  if (nsize != 6) {
    cout << "NewTemplate3Dep::Tensor2VectorSysR2 - matrix must be of size 6" << endln;
    return 1;
  }

  double sqrt2 = sqrt(2.0);

  // Adopt method from Helnwein (2001):

  V(0) = T.cval(1,1);
  V(1) = T.cval(2,2);
  V(2) = T.cval(3,3);
  V(3) = T.cval(1,2) *sqrt2;
  V(4) = T.cval(2,3) *sqrt2;
  V(5) = T.cval(1,3) *sqrt2;

  return 0;
}

//================================================================================
int NewTemplate3Dep::Matrix2TensorSysR4(const Matrix& M, tensor& T)
{
  int rank = T.rank();
  if (rank != 4) {
    cout << "NewTemplate3Dep::Matrix2TensorSysR4 - tensor must be of rank 4" << endln;
    return 1;
  }

  int nr = M.noRows();
  int nc = M.noCols();
  if (nr < 6 || nc < 6) {
    cout << "NewTemplate3Dep::Matrix2TensorSysR4 - matrix must be no less than (6, 6)" << endln;
    return 1;
  }

  double sqrthalf = sqrt(0.5);
  double half = 0.5;

  // Adopt method from Helnwein (2001):

  T.val(1,1,1,1) = M(0,0);
  T.val(1,1,2,2) = M(0,1);
  T.val(1,1,3,3) = M(0,2);
  T.val(1,1,1,2) = T.val(1,1,2,1) = M(0,3) *sqrthalf;
  T.val(1,1,2,3) = T.val(1,1,3,2) = M(0,4) *sqrthalf;
  T.val(1,1,1,3) = T.val(1,1,3,1) = M(0,5) *sqrthalf;

  T.val(2,2,1,1) = M(1,0);
  T.val(2,2,2,2) = M(1,1);
  T.val(2,2,3,3) = M(1,2);
  T.val(2,2,1,2) = T.val(2,2,2,1) = M(1,3) *sqrthalf;
  T.val(2,2,2,3) = T.val(2,2,3,2) = M(1,4) *sqrthalf;
  T.val(2,2,1,3) = T.val(2,2,3,1) = M(1,5) *sqrthalf;

  T.val(3,3,1,1) = M(2,0);
  T.val(3,3,2,2) = M(2,1);
  T.val(3,3,3,3) = M(2,2);
  T.val(3,3,1,2) = T.val(3,3,2,1) = M(2,3) *sqrthalf;
  T.val(3,3,2,3) = T.val(3,3,3,2) = M(2,4) *sqrthalf;
  T.val(3,3,1,3) = T.val(3,3,3,1) = M(2,5) *sqrthalf;

  T.val(1,2,1,1) = T.val(2,1,1,1) = M(3,0) *sqrthalf;
  T.val(1,2,2,2) = T.val(2,1,2,2) = M(3,1) *sqrthalf;
  T.val(1,2,3,3) = T.val(2,1,3,3) = M(3,2) *sqrthalf;
  T.val(1,2,1,2) = T.val(2,1,1,2) = T.val(1,2,2,1) = T.val(2,1,2,1) = M(3,3) *half;
  T.val(1,2,2,3) = T.val(2,1,2,3) = T.val(1,2,3,2) = T.val(2,1,3,2) = M(3,4) *half;
  T.val(1,2,1,3) = T.val(2,1,1,3) = T.val(1,2,3,1) = T.val(2,1,3,1) = M(3,5) *half;

  T.val(2,3,1,1) = T.val(3,2,1,1) = M(4,0) *sqrthalf;
  T.val(2,3,2,2) = T.val(3,2,2,2) = M(4,1) *sqrthalf;
  T.val(2,3,3,3) = T.val(3,2,3,3) = M(4,2) *sqrthalf;
  T.val(2,3,1,2) = T.val(3,2,1,2) = T.val(2,3,2,1) = T.val(3,2,2,1) = M(4,3) *half;
  T.val(2,3,2,3) = T.val(3,2,2,3) = T.val(2,3,3,2) = T.val(3,2,3,2) = M(4,4) *half;
  T.val(2,3,1,3) = T.val(3,2,1,3) = T.val(2,3,3,1) = T.val(3,2,3,1) = M(4,5) *half;

  T.val(1,3,1,1) = T.val(3,1,1,1) = M(5,0) *sqrthalf;
  T.val(1,3,2,2) = T.val(3,1,2,2) = M(5,1) *sqrthalf;
  T.val(1,3,3,3) = T.val(3,1,3,3) = M(5,2) *sqrthalf;
  T.val(1,3,1,2) = T.val(3,1,1,2) = T.val(1,3,2,1) = T.val(3,1,2,1) = M(5,3) *half;
  T.val(1,3,2,3) = T.val(3,1,2,3) = T.val(1,3,3,2) = T.val(3,1,3,2) = M(5,4) *half;
  T.val(1,3,1,3) = T.val(3,1,1,3) = T.val(1,3,3,1) = T.val(3,1,3,1) = M(5,5) *half;

  return 0;
}

//================================================================================
int NewTemplate3Dep::Vector2TensorSysR2(const Vector& V, tensor& T, int Num0)
{
  int rank = T.rank();
  if (rank != 2) {
    cout << "NewTemplate3Dep::Vector2TensorSysR2 - tensor must be of rank 4" << endln;
    return 1;
  }

  int nsize = V.Size();
  if (nsize < 6) {
    cout << "NewTemplate3Dep::Vector2TensorSysR2 - matrix less than size 6" << endln;
    return 1;
  }

  double sqrthalf = sqrt(0.5);

  // Adopted method from Helnwein (2001):

  T.val(1,1) = V(Num0+0);
  T.val(2,2) = V(Num0+1);
  T.val(3,3) = V(Num0+2);
  T.val(1,2) = T.val(2,1) = V(Num0+3) *sqrthalf;
  T.val(2,3) = T.val(3,2) = V(Num0+4) *sqrthalf;
  T.val(1,3) = T.val(3,1) = V(Num0+5) *sqrthalf;

  return 0;
}


//================================================================================
int NewTemplate3Dep::Stiffness2Compliance(const tensor& S, tensor& C)
{
  int rank = 0;
  int err = 0;

  rank = S.rank();
  if (rank != 4) {
    cout << "NewTemplate3Dep::Stiffness2Compliance - tensor must be of rank 4" << endln;
    return 1;
  }

  rank = C.rank();
  if (rank != 4) {
    cout << "NewTemplate3Dep::Stiffness2Compliance - tensor must be of rank 4" << endln;
    return 1;
  }

  Matrix S66(6,6);
  Matrix C66(6,6);

  err += Tensor2MatrixSysR4(S, S66);
  err += S66.Invert(C66);
  err += Matrix2TensorSysR4(C66, C);

  return err;
}


#endif

