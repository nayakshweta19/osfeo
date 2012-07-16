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

// This is the container to store material constants:
// Material_Constant:       to store fixed scalar material constants;
//                          Note: the first one should be the density, input 0.0 if not involved;
// Num_Material_Constant:   the number of Material_Constant, should not less than 2, 
//                          given the 1st is the density and the 2nd is void ratio;
// Internal_Scalar:         to store initial evolution scalar variables first time, 
//                          and scalar variables thereafter;
// Num_Internal_Scalar:     the number of Internal_Scalar;
// Internal_Tensor:         to store initial evolution tensor variables first time (usually zero tensor), 
//                          and tensor variables thereafter;
// Num_Internal_Tensor:     the number of Internal_Tensor;


#ifndef MaterialParameter_CPP
#define MaterialParameter_CPP

#include "MaterialParameter.h"
#define MP_TAG 120001

// Constructor 1
MaterialParameter::MaterialParameter(const double *Material_Constant_in, 
                                    int Num_Material_Constant_in, 
                                    const double *Internal_Scalar_in, 
                                    int Num_Internal_Scalar_in, 
                                    const stresstensor *Internal_Tensor_in, 
                                    int Num_Internal_Tensor_in) : MovableObject(MP_TAG)
{
    // Material Constants
    if (Num_Material_Constant_in > 0) {
        Material_Constant = new double [Num_Material_Constant = Num_Material_Constant_in];
        if (!Material_Constant) {
            opserr << "MaterialConstant::Insufficient memory! " << endln;
            exit (1);
        }        
        for (int i = 0; i < Num_Material_Constant; i++)
            Material_Constant[i] = Material_Constant_in[i];
    }
    else {
      Num_Material_Constant = 0;
      Material_Constant = NULL;
    }

    // Scalar (Isotropic) Internal Variables
    if (Num_Internal_Scalar_in > 0) {
        Internal_Scalar = new double [Num_Internal_Scalar = Num_Internal_Scalar_in];
        if (!Internal_Scalar) {
            opserr << "MaterialParameter::Insufficient memory! " << endln;
            exit (1);
        }        
        for (int i = 0; i < Num_Internal_Scalar; i++)
            Internal_Scalar[i] = Internal_Scalar_in[i];
    }
    else {
      Num_Internal_Scalar = 0;
      Internal_Scalar = NULL;
    }

    // Tensor (Kinematic) Internal Variables
    if (Num_Internal_Tensor_in > 0) {
        Internal_Tensor = new stresstensor[Num_Internal_Tensor = Num_Internal_Tensor_in];
        if (!Internal_Tensor) {
            opserr << "MaterialParameter::Insufficient memory! " << endln;
            exit (1);
        }        
        for (int i = 0; i < Num_Internal_Tensor; i++)
            Internal_Tensor[i] = Internal_Tensor_in[i];
    }
    else {
      Num_Internal_Tensor = 0;
      Internal_Tensor = NULL;
    }
}

// Constructor 2
MaterialParameter::MaterialParameter(const double *Material_Constant_in, 
                                     int Num_Material_Constant_in, 
                                     const stresstensor *Internal_Tensor_in, 
                                     int Num_Internal_Tensor_in)
: MovableObject(MP_TAG), Internal_Scalar(NULL), Num_Internal_Scalar(0)
{
    // Material Constants
    if (Num_Material_Constant_in > 0) {
        Material_Constant = new double [Num_Material_Constant = Num_Material_Constant_in];
        if (!Material_Constant) {
            opserr << "MaterialParameter::Insufficient memory! " << endln;
            exit (1);
        }        
        for (int i = 0; i < Num_Material_Constant; i++)
            Material_Constant[i] = Material_Constant_in[i];
    }
    else {
      Num_Material_Constant = 0;
      Material_Constant = NULL;
    }

    // Tensor (Kinematic) Internal Variables
    if (Num_Internal_Tensor_in > 0) {
        Internal_Tensor = new stresstensor[Num_Internal_Tensor = Num_Internal_Tensor_in];
        if (!Internal_Tensor) {
            opserr << "MaterialParameter::Insufficient memory! " << endln;
            exit (1);
        }        
        for (int i = 0; i < Num_Internal_Tensor; i++)
            Internal_Tensor[i] = Internal_Tensor_in[i];
    }
    else {
      Num_Internal_Tensor = 0;
      Internal_Tensor = NULL;
    }
}

//Guanzhou // Constructor 3
//Guanzhou MaterialParameter::MaterialParameter( )
//Guanzhou : Material_Constant(NULL), Num_Material_Constant(0), 
//Guanzhou   Internal_Scalar(NULL), Num_Internal_Scalar(0), 
//Guanzhou   Internal_Tensor(NULL), Num_Internal_Tensor(0)
//Guanzhou {
//Guanzhou 
//Guanzhou }

// Destructor
MaterialParameter::~MaterialParameter( )
{
	if (Material_Constant)
		delete [] Material_Constant;

	if (Internal_Scalar)
		delete [] Internal_Scalar;

	if (Internal_Tensor)
		delete [] Internal_Tensor;
}

// Copy constructor
MaterialParameter::MaterialParameter(const MaterialParameter &refer_MaterialParameter )
: MovableObject(MP_TAG)
{
    // Material Constants
    if (refer_MaterialParameter.getNum_Material_Constant() > 0) {
        Material_Constant = new double [Num_Material_Constant = refer_MaterialParameter.getNum_Material_Constant()];
        if (!Material_Constant) {
            opserr << "MaterialParameter::Insufficient memory! " << endln;
            exit (1);
        }        
        for (int i = 0; i < Num_Material_Constant; i++)
            Material_Constant[i] = refer_MaterialParameter.getMaterial_Constant(i);
    }
    else {
      Num_Material_Constant = 0;
      Material_Constant = NULL;
    }

    // Scalar (Isotropic) Internal Variables
    if (refer_MaterialParameter.getNum_Internal_Scalar() > 0) {
        Internal_Scalar = new double [Num_Internal_Scalar = refer_MaterialParameter.getNum_Internal_Scalar()];
        if (!Internal_Scalar) {
            opserr << "MaterialParameter::Insufficient memory! " << endln;
            exit (1);
        }        
        for (int i = 0; i < Num_Internal_Scalar; i++)
            Internal_Scalar[i] = refer_MaterialParameter.getInternal_Scalar(i);
    }
    else {
      Num_Internal_Scalar = 0;
      Internal_Scalar = NULL;
    }

    // Tensor (Kinematic) Internal Variables
    if (refer_MaterialParameter.getNum_Internal_Tensor() > 0) {
        Internal_Tensor = new stresstensor[Num_Internal_Tensor = refer_MaterialParameter.getNum_Internal_Tensor()];
        if (!Internal_Tensor) {
            opserr << "MaterialParameter::Insufficient memory! " << endln;
            exit (1);
        }        
        for (int i = 0; i < Num_Internal_Tensor; i++)
            Internal_Tensor[i] = refer_MaterialParameter.getInternal_Tensor(i);
    }
    else {
      Num_Internal_Tensor = 0;
      Internal_Tensor = NULL;
    }

}

// Create a new class pointer
MaterialParameter* MaterialParameter::newObj() {
	MaterialParameter *ptr_MaterialParameter = new  MaterialParameter ( this->Material_Constant,
                                                                        this->Num_Material_Constant,
                                                                        this->Internal_Scalar,
                                                                        this->Num_Internal_Scalar,
                                                                        this->Internal_Tensor,
                                                                        this->Num_Internal_Tensor );
    return ptr_MaterialParameter;	
}

// getNum_Material_Parameter
int MaterialParameter::getNum_Material_Constant() const
{
	return Num_Material_Constant;
}

// getNum_Internal_Scalar
int MaterialParameter::getNum_Internal_Scalar() const
{
	return Num_Internal_Scalar;
}

// getNum_Internal_Tensor
int MaterialParameter::getNum_Internal_Tensor() const
{
	return Num_Internal_Tensor;
}

// getMaterial_Parameter
double MaterialParameter::getMaterial_Constant(int which) const
{
	if (which < 0 || which > Num_Material_Constant) {
		opserr << "Error! MaterialParameter::getMaterial_Constant - Invalid index of material constants. " << endln;
		exit(1);
	}	

	return Material_Constant[which];
}

// getInternal_Scalar
double MaterialParameter::getInternal_Scalar(int which) const
{
	if (which < 0 || which > Num_Internal_Scalar) {
		opserr << "Error! MaterialParameter::getInternal_Scalar - Invalid index of internal scalars. " << endln;
		exit(1);
	}	

	return Internal_Scalar[which];
}

// getInternal_Tensor
const stresstensor& MaterialParameter::getInternal_Tensor(int which) const
{
	if (which < 0 || which > Num_Internal_Tensor) {
		opserr << "Error! MaterialParameter::getInternal_Tensor - Invalid index of internal tensors. " << endln;
		exit (1);
	}	

	return Internal_Tensor[which];
}

// setMaterial_Parameter
int MaterialParameter::setMaterial_Constant(int which, double newMaterial_Constant) 
{
	if (which < 0 || which > Num_Material_Constant) {
		opserr << "Error! MaterialParameter::setMaterial_Constant - Invalid index of material constants. " << endln;
		return (1);
	}	
	
	Material_Constant[which] = newMaterial_Constant;
	
	return 0;
}

// setInternal_Scalar
int MaterialParameter::setInternal_Scalar(int which, double newInternal_Scalar) 
{
	if (which < 0 || which > Num_Internal_Scalar) {
		opserr << "Error! MaterialParameter::setInternal_Scalar - Invalid index of internal scalars. " << endln;
		exit (1);
	}	
	
	Internal_Scalar[which] = newInternal_Scalar;
	
	return 0;
}

// setInternal_Tensor
int MaterialParameter::setInternal_Tensor(int which, const stresstensor &newInternal_Tensor) 
{
	if (which < 0 || which > Num_Internal_Tensor) {
		opserr << "Error! MaterialParameter::setInternal_Tensor - Invalid index of internal tensors. " << endln;
		exit (1);
	}	
		
	Internal_Tensor[which] = newInternal_Tensor;

// Nima Tafazzoli
//        newInternal_Tensor.print("alpha after dfods" , "alpha after dfods");
				
	return 0;
}


//Guanzhou added for parallel
int MaterialParameter::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    
    static ID idData(3);
    idData.Zero();

    idData(0) = Num_Material_Constant;
    idData(1) = Num_Internal_Scalar;
    idData(2) = Num_Internal_Tensor;
    
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
   	opserr << "MaterialParameter::sendSelf -- failed to send ID\n";
   	return -1;
    }

    if ( Num_Material_Constant > 0 ) {
        static Vector data(Material_Constant, Num_Material_Constant);

        if (theChannel.sendVector(dataTag, commitTag, data) < 0) {
   	    opserr << "MaterialParameter::sendSelf -- failed to send Vector Material_Constant\n";
   	    return -1;
        }
    }

    if ( Num_Internal_Scalar > 0 ) {
        static Vector data(Internal_Scalar, Num_Internal_Scalar);

        if (theChannel.sendVector(dataTag, commitTag, data) < 0) {
   	    opserr << "MaterialParameter::sendSelf -- failed to send Vector Internal_Scalar\n";
   	    return -1;
        }
    }
    
    if ( Num_Internal_Tensor > 0 ) {
	for (int i = 0; i < Num_Internal_Tensor; i++) {
            if (theChannel.sendnDarray(dataTag, commitTag, Internal_Tensor[i]) < 0) {
   	    	opserr << "MaterialParameter::sendSelf -- failed to send Internal_Tensors\n";
   	    	return -1;
            }
	}
    }
    
    return 0;

}

//Guanzhou added for parallel
int MaterialParameter::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    
    static ID idData(3);
    static Vector *data = NULL;
    idData.Zero();

    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
    	opserr << "MaterialParameter::recvSelf -- failed to recv ID\n";
	return -1;
    }

    Num_Material_Constant = idData(0);
    Num_Internal_Scalar   = idData(1);
    Num_Internal_Tensor   = idData(2);
    
    if ( Num_Material_Constant > 0 ) {
    	Material_Constant = new double [Num_Material_Constant];
        if (Material_Constant == NULL) {
            opserr << "MaterialParameter::Insufficient memory! " << endln;
            exit (1);
        }        
	//static Vector 
	data = new Vector(Material_Constant, Num_Material_Constant);
	if ( theChannel.recvVector(dataTag, commitTag, *data) < 0 ) {
    	    opserr << "MaterialParameter::recvSelf -- failed to recv Material_Constant\n";
	    return -1;
	}
	delete data;
    }
        
    if ( Num_Internal_Scalar > 0 ) {
    	Internal_Scalar = new double [Num_Internal_Scalar];
        if (Internal_Scalar == NULL) {
            opserr << "Internal_Scalar::Insufficient memory! " << endln;
            exit (1);
        }        
	static Vector data(Internal_Scalar, Num_Internal_Scalar);
	if ( theChannel.recvVector(dataTag, commitTag, data) < 0 ) {
    	    opserr << "MaterialParameter::recvSelf -- failed to recv Internal_Scalar\n";
	    return -1;
	}
    }
    
    if ( Num_Internal_Tensor > 0 ) {
    	Internal_Tensor = new stresstensor[Num_Internal_Tensor];
        if (Internal_Tensor == NULL) {
            opserr << "Internal_Tensor::Insufficient memory! " << endln;
            exit (1);
        }        
	for (int i = 0; i < Num_Internal_Tensor; i++) {
	    if ( theChannel.recvnDarray(dataTag, commitTag, Internal_Tensor[i]) < 0 ) {
    	    	opserr << "MaterialParameter::recvSelf -- failed to recv Internal_Tensor\n";
	    	return -1;
	    }
	}
    }
    
    return 0;
}


#endif

