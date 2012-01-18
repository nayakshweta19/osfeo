/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.0 $
// $Date: 2011/08/30 23:00:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DSFMRCPlaneStress.cpp,v $
// File: DSFMRCPlaneStress.cpp
//
// Written: Lining
// Created: 2011.8
//
// Description: This file contains the class definition for
// DSFMRCPlaneStress
// For Detailed explanation of the model, please refer to the journal paper
// entitled "The modified compression-filed theory for reinforced concrete element
// subjected to shear, ACI Journal. s83-22. 1986. pp. 219-231"

#include "DSFMRCPlaneStress.h"
#include <UniaxialMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <math.h>
#include <float.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <stdlib.h>

#include <Tensor.h>

#include <DummyStream.h>
#include <MaterialResponse.h>
#include <elementAPI.h>

using namespace std;

#define OPS_Export 

static int numDSFMRCPlaneStressMaterials = 0;

OPS_Export void *
OPS_NewDSFMRCPlaneStressMaterial()
{
  if (numDSFMRCPlaneStressMaterials == 0) {
    numDSFMRCPlaneStressMaterials++;
    //OPS_Error("DSFMRCPlaneStress material - Written by Lining - Copyright@2011\n", 1);
  }

  // Pointer to a uniaxial material that will be returned
  NDMaterial *theMaterial = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 16) {
    opserr << "Invalid Args want: NDMaterial DSFMRCPlaneStress matTag? rho? UniaxiaMatTag1? UniaxiaMatTag2? UniaxiaMatTag3? UniaxiaMatTag4? \
               angle1? angle2? roux? rouy? db1? db2? fpc? fy? E0? epsc0? aggsize? xd? yd?\n";
    return 0;	
  }

  int tag;
  double rho;
  int    iData[4];
  double dData[14];
  int numData = 0;

  numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid nDMaterial DSFMRCPlaneStress tag" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDouble(&numData, &rho) != 0) {
    opserr << "WARNING invalid parameter rho DSFMRCPlaneStress tag:" << tag << endln;
    return 0;	
  }

  numData = 4;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid nDMaterial DSFMRCPlaneStress tag: " << tag << endln;
    return 0;
  }

  numData = 14;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data DSFMRCPlaneStress tag: " << tag << endln;
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial1 = OPS_GetUniaxialMaterial(iData[0]);
    
  if (theUniaxialMaterial1 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[0];
    opserr << "\nDSFMRCPlaneStress tag: " << tag << endln;
    return 0;
  }
  
  UniaxialMaterial *theUniaxialMaterial2 = OPS_GetUniaxialMaterial(iData[1]);

  if (theUniaxialMaterial2 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[1];
    opserr << "\nDSFMRCPlaneStress tag: " << tag << endln;
    return 0;
  }
  
  UniaxialMaterial *theUniaxialMaterial3 = OPS_GetUniaxialMaterial(iData[2]);
  if (theUniaxialMaterial3 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[2];
    opserr << "\nDSFMRCPlaneStress tag: " << tag << endln;
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial4 = OPS_GetUniaxialMaterial(iData[3]);  
  if (theUniaxialMaterial4 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[3];
    opserr << "\nDSFMRCPlaneStress tag: " << tag << endln;
    return 0;
  }

  //now create the DSFMRCPlaneStress
  theMaterial = new DSFMRCPlaneStress (tag, rho,
						     theUniaxialMaterial1, theUniaxialMaterial2, 
						     theUniaxialMaterial3, theUniaxialMaterial4,
						     dData[0], dData[1], dData[2], dData[3], dData[4],
							 dData[5], dData[6], dData[7], dData[8], dData[9],
							 dData[10],dData[11],dData[12],dData[13]);
       
  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory creating material\n";
    opserr << "DSFMRCPlaneStress: " << tag << endln;
    return 0;
  }

  return theMaterial;
}
 
DSFMRCPlaneStress ::DSFMRCPlaneStress (int tag, 
								   double RHO,
								   UniaxialMaterial *s1,
								   UniaxialMaterial *s2,
								   UniaxialMaterial *c1,
								   UniaxialMaterial *c2,
								   double   ANGLE1,
								   double   ANGLE2,
								   double   ROU1,
								   double   ROU2,
								   double   DB1,
								   double   DB2,
								   double   FPC,
								   double   FYX,
								   double   FYY,
								   double   E,
								   double   EPSC0,
								   double   AGGR,
								   double   XD,
								   double   YD):
  NDMaterial(tag, ND_TAG_DSFMRCPlaneStress), 
  rho(RHO),angle1(ANGLE1),angle2(ANGLE2),roux(ROU1),rouy(ROU2),db1(DB1),db2(DB2),
  fpc(FPC), fyx(FYX), fyy(FYY), E0(E), epsc0(EPSC0), aggr(AGGR), xd(XD), yd(YD),
  lastStress(3), Tstress(3), stress_vec(3), stress0_vec(3), tangent_matrix(3,3), strain_vec(3), 
  epsC_vec(3), epsC12p_prevec(2), epsC12p_nowvec(2), epsSlip_vec(3), epsC0_vec(3), epsS0_vec(3), epsCp_vec(3),
  epsC12cm_vec(2), epsCcm_vec(3), epsC12tm_vec(2), epsCtm_vec(3)
{
    
    lastStress.Zero();  // add at 7.28    
    
    if ( fpc > 0.0 ) { fpc = -fpc; } // set fpc < 0
    if ( epsc0 > 0.0 ) { epsc0 = -epsc0; } // set fpc < 0
    theMaterial = 0;
    
    // Allocate pointers to theSteel1
    theMaterial = new UniaxialMaterial *[4];
    
    if ( theMaterial == 0 ) {
      opserr << " DSFMRCPlaneStress::DSFMRCPlaneStress - failed allocate material array\n";
      exit(-1);
    }
    
    // Get the copy for the Steel1
    theMaterial[0] = s1->getCopy();
    // Check allocation    
    if ( theMaterial[0] == 0 ) {
      opserr << " DSFMRCPlaneStress::DSFMRCPlaneStress - failed to get a copy for steel1\n";
      exit(-1);
    }
    
    // Get the copy for the Steel2
    theMaterial[1] = s2->getCopy();	
    // Check allocation    
    if ( theMaterial[1] == 0 ) {
      opserr << " DSFMRCPlaneStress::DSFMRCPlaneStress - failed to get a copy for steel2\n";
      exit(-1);
    }
    
    // Get the copy for the Concrete1
    theMaterial[2] = c1->getCopy();	
    // Check allocation    
    if ( theMaterial[2] == 0 ) {
      opserr << " DSFMRCPlaneStress::DSFMRCPlaneStress - failed to get a copy for concrete1\n";
      exit(-1);
    }
    
    // Get the copy for the Concrete2
    theMaterial[3] = c2->getCopy();	
    // Check allocation    
    if ( theMaterial[3] == 0 ) {
      opserr << " DSFMRCPlaneStress::DSFMRCPlaneStress - failed to get a copy for concrete2\n";
      exit(-1);
    }
    
    /* FMK */
    theResponses = new Response *[6];  
    
    if ( theResponses == 0) {
      opserr << " DSFMRCPlaneStress::DSFMRCPlaneStress - failed allocate responses array\n";
      exit(-1);
    }
    
    OPS_Stream *theDummyStream = new DummyStream();
    
    const char **argv = new const char *[1];

    argv[0] = "getCommittedStrain";
    theResponses[0] = theMaterial[0]->setResponse(argv, 1, *theDummyStream);
    theResponses[1] = theMaterial[1]->setResponse(argv, 1, *theDummyStream);
    argv[0] = "setWallVar";
    theResponses[2] = theMaterial[2]->setResponse(argv, 1, *theDummyStream);
    theResponses[3] = theMaterial[3]->setResponse(argv, 1, *theDummyStream);
    argv[0] = "getPD";
    theResponses[4] = theMaterial[2]->setResponse(argv, 1, *theDummyStream);
    theResponses[5] = theMaterial[3]->setResponse(argv, 1, *theDummyStream);    

    if ((theResponses[0] == 0) || (theResponses[1] == 0) ||
	(theResponses[2] == 0) || (theResponses[3] == 0) ||
	(theResponses[4] == 0) || (theResponses[5] == 0)) {
      opserr << " DSFMRCPlaneStress::DSFMRCPlaneStress - failed to set appropriate materials tag: " << tag << "\n";
      exit(-1);
    }
    
    delete theDummyStream;
	strain_vec.Zero();
    this->revertToStart();
	determineTrialStress(strain_vec);
}

DSFMRCPlaneStress::DSFMRCPlaneStress()
  :NDMaterial(0, ND_TAG_DSFMRCPlaneStress),
  lastStress(3), Tstress(3), stress_vec(3), stress0_vec(3), tangent_matrix(3,3),strain_vec(3), 
  epsC12p_prevec(2), epsC12p_nowvec(2), epsC_vec(3), epsSlip_vec(3), epsC0_vec(3), epsS0_vec(3), epsCp_vec(3),
  epsC12cm_vec(2), epsCcm_vec(3), epsC12tm_vec(2), epsCtm_vec(3)
{
  theMaterial = 0;
  theResponses = 0;
  strain_vec.Zero();
  this->revertToStart();
  determineTrialStress(strain_vec);
}

DSFMRCPlaneStress::~DSFMRCPlaneStress()
{
  // Delete the pointers
  if (theMaterial != 0) {
    for (int i=0; i<4; i++) {
	  if (theMaterial[i])
	    delete theMaterial[i];
    }
    delete [] theMaterial;
  }
  if (theResponses != 0) {
    for (int j=0; j<6; j++) {
	  if (theResponses[j] != 0)
	    delete theResponses[j];
    }
    delete [] theResponses;
  }
}

double 
DSFMRCPlaneStress::getRho(void)
{
	return rho;
}

// really used one
int 
DSFMRCPlaneStress::setTrialStrain(const Vector &v)
{

  // Set values for strain_vec
  strain_vec = v;
  
  // Set initial values for Tstress
  Tstress.Zero();
  
  Tstress = determineTrialStress(strain_vec);
  
  //determineTangent();

  return 0;
}

int 
DSFMRCPlaneStress::setTrialStrain(const Vector &v, const Vector &r)
{
  opserr << "error: DSFMRCPlaneStress::setTrialStrain(&v, &r) -- not really responsibility" << endln;
  return 0;
}

int 
DSFMRCPlaneStress::setTrialStrainIncr(const Vector &v)
{
  opserr << "error:DSFMRCPlaneStress::setTrialStrainIncr(&v) -- not really responsibility" << endln;
  return 0;
}

int
DSFMRCPlaneStress::setTrialStrainIncr(const Vector &v, const Vector &r)
{
  opserr << "error:DSFMRCPlaneStress::setTrialStrainIncr(&v, &r) -- not really responsibility" << endln;
  return 0;
}

const Matrix& 
DSFMRCPlaneStress::getTangent(void)
{

  determineTrialStress(strain_vec);

  //determineTangent();

  return tangent_matrix;
}

const Vector& 
DSFMRCPlaneStress::getStress(void)
{
  return stress_vec;
}

const Vector& 
DSFMRCPlaneStress::getStrain()
{
  return strain_vec;
}

const Vector& 
DSFMRCPlaneStress::getCommittedStress(void)
{
  return stress_vec;
}

const Vector& 
DSFMRCPlaneStress::getCommittedStrain(void)
{
  return strain_vec;
}

int
DSFMRCPlaneStress::commitState(void)
{
  for (int i=0; i<4; i++) {
    theMaterial[i]->commitState();
  }
  
  lastStress = stress_vec;
  
  return 0;
}

int
DSFMRCPlaneStress::revertToLastCommit(void)
{
  for (int i=0; i<4; i++) {
    theMaterial[i]->revertToLastCommit();
  }

  return 0;
}

int
DSFMRCPlaneStress::revertToStart(void)
{
  for (int i=0; i<4; i++) {
    theMaterial[i]->revertToStart();
  }
  
  Tstress.Zero();
  strain_vec.Zero();
  epsC_vec.Zero();
  epsC12p_prevec.Zero();
  epsC12p_nowvec.Zero();
  epsSlip_vec.Zero();
  epsC0_vec.Zero();
  epsS0_vec.Zero();
  epsCp_vec.Zero();
  stress_vec.Zero();
  stress0_vec.Zero();
  epsC12cm_vec.Zero();
  epsCcm_vec.Zero();
  epsC12tm_vec.Zero();
  epsCtm_vec.Zero();
  
  return 0;
}

NDMaterial* DSFMRCPlaneStress::getCopy(void)
{
  DSFMRCPlaneStress* theCopy =
    new DSFMRCPlaneStress( this->getTag(), 
					 rho, 
					 theMaterial[0], 
					 theMaterial[1], 
					 theMaterial[2], 
					 theMaterial[3], 
					 angle1, angle2, 
					 roux, rouy, 
					 db1, db2,
					 fpc, 
					 fyx, fyy,
					 E0, 
					 epsc0,
					 aggr,
					 xd, yd);
  theCopy->strain_vec = strain_vec;
  return theCopy;
}

NDMaterial* DSFMRCPlaneStress::getCopy(const char *type)
{
	DSFMRCPlaneStress* theModel =
		new DSFMRCPlaneStress( this->getTag(), 
           rho, theMaterial[0], theMaterial[1], theMaterial[2], theMaterial[3],
		   angle1, angle2, roux, rouy, db1, db2, fpc, fyx, fyy, E0, epsc0, aggr, xd, yd );
	theModel->strain_vec = strain_vec;
	theModel->stress_vec = stress_vec;
	return theModel;
}

Response*
DSFMRCPlaneStress::setResponse (const char **argv, int argc, OPS_Stream &output)
{

#ifdef DEBUG
	opserr << "DSFMRCPlaneStress::setResponse(...)" << endln;
#endif

	//Response *theResponse =0;
	const char *matType = this->getType();

	output.tag("NdMaterialOutput");
	output.attr("matType",this->getClassType());
	output.attr("matTag",this->getTag());

	if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());
	else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
	//else if (strcmp(argv[0], "state") == 0)
		//return new MaterialResponse(this, 3, this->getState());
	else
		return 0;
}

int
DSFMRCPlaneStress::getResponse (int responseID, Information &matInfo)
{
#ifdef DEBUG
	opserr << "DSFMRCPlaneStress::getResponse(...)" << endln;
#endif

	switch (responseID) {
	case -1:
		return -1;
	case 1:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStress();
		return 0;
	case 2:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStrain();
		return 0;
	case 3:
		//if (matInfo.theVector != 0)
		//	*(matInfo.theVector) = getState();
		return 0;
	default:
		return -1;
	}
}
//end by LN

void
DSFMRCPlaneStress::Print(OPS_Stream &s, int flag )
{
	s << "\n\tDSFMRCPlaneStress, material id: " << this->getTag() << endln;
    /*
	s << "\tRho: " << rho << endln;
	s << "\tangle1: " << angle1 << endln;
	s << "\tangle2: " << angle2 << endln;
	s << "\trou1: " << roux << endln;
	s << "\trou2: " << rouy << endln;
	s << "\tfpc: " << fpc << endln;
	s << "\tfy: " << fy << endln;
	s << "\tE0: " << E0 << endln;
	//*/

	//s <<  "Principal Strain: citaStrain = "<< citaStrain/3.14159*180.0 << endln;
	//s <<  "Principal Stress: citaStress = "<< citaStress/3.14159*180.0 << endln;
	//s << " v12 = " << miu12 << " v21 = " << miu21 << endln;	
	//s << " steelStatus " << steelStatus << endln;
	//s << " dirStatus " << dirStatus << endln;
	//s << " Damage DOne = " << DDOne << endln;
	//s << " Damage DTwo = " << DDTwo << endln;
	
	//s << " G12 = " << G12 << endln;

	//int i, j;

	s << "\t call the material print() function : "<< endln;
	
	s << "\t the steel 1 information is : " << endln;
	theMaterial[0]->Print(s,flag);
	s << "\t the steel 2 information is : " << endln;
	theMaterial[1]->Print(s,flag);
	s << "\t the concrete 1 information is : " << endln;
	theMaterial[2]->Print(s,flag);
	s << "\t the concrete 2 information is : " << endln;
	theMaterial[3]->Print(s,flag);

	//s << "\tStrain and stress of the uniaxial materials:"<<endln;
	//for ( i=0; i<4; i++)
	//{
	//	s<< "Uniaxial Material "<<i+1<<" :"<<theMaterial[i]->getStrain()<<"   "<<theMaterial[i]->getStress()<< endln;
	//}
}

int
DSFMRCPlaneStress::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  int dataTag = this->getDbTag();
  
  // Packs its data into a Vector and sends this to theChannel
  static Vector data(16);
  data(0) = this->getTag();
  data(1) = rho;
  data(2) = angle1;
  data(3) = angle2;
  data(4) = roux;
  data(5) = rouy;
  data(6) = db1;
  data(7) = db2;
  data(8) = fpc;
  data(9) = fyx;
  data(10) = fyy;
  data(11) = E0;
  data(12) = epsc0;
  data(13) = aggr;
  data(14) = xd;
  data(15) = yd;

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING DSFMRCPlaneStress::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }	      
  
  // Now sends the IDs of its materials
  int matDbTag;
  
  static ID idData(8);
  
  // NOTE: to ensure that the material has a database
  // tag if sending to a database channel.
  
  int i;
  for (i=0; i<4; i++) {
    idData(i) = theMaterial[i]->getClassTag();
    matDbTag = theMaterial[i]->getDbTag();
    if (matDbTag == 0) {
	  matDbTag = theChannel.getDbTag();
	  if (matDbTag != 0)
	    theMaterial[i]->setDbTag(matDbTag);
    }
    idData(i+4) = matDbTag;
  }
  
  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING DSFMRCPlaneStress::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }
  
  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 4; i++) {
    res += theMaterial[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "DSFMRCPlaneStress::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }	
  
  return res;
}

int
DSFMRCPlaneStress::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();
  
  // Quad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(16);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING DSFMRCPlaneStress::recvSelf() - failed to receive Vector\n";
    return res;
  }
  
  this->setTag((int)data(0));
  rho = data(1);
  angle1 = data(2);
  angle2 = data(3);
  roux   = data(4);
  rouy   = data(5);
  db1    = data(6);
  db2    = data(7);
  fpc    = data(8);
  fyx    = data(9);
  fyy    = data(10);
  E0     = data(11);
  epsc0  = data(12);
  aggr   = data(13);
  xd     = data(14);
  yd     = data(15);

  static ID idData(8);
  
  // now receives the tags of its materials
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING DSFMRCPlaneStress::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  if (theMaterial == 0) {
    // Allocate new materials
    theMaterial = new UniaxialMaterial *[4];
    if (theMaterial == 0) {
      opserr << "DSFMRCPlaneStress::recvSelf() - Could not allocate UniaxialMaterial* array\n";
      return -1;
    }
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Allocate new material with the sent class tag
      theMaterial[i] = theBroker.getNewUniaxialMaterial(matClassTag);
      if (theMaterial[i] == 0) {
	     opserr << "DSFMRCPlaneStress::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
	     return -1;
      }
      // Now receive materials into the newly allocated space
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
         opserr << "DSFMRCPlaneStress::recvSelf() - material " << i << "failed to recv itself\n";
	     return res;
      }
    }
  }

  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (theMaterial[i]->getClassTag() != matClassTag) {
	     delete theMaterial[i];
	     theMaterial[i] = theBroker.getNewUniaxialMaterial(matClassTag);
	     if (theMaterial[i] == 0) {
            opserr << "DSFMRCPlaneStress::recvSelf() - material " << i << "failed to create\n";
	        return -1;
		 }
	  }
      // Receive the material
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "DSFMRCPlaneStress::recvSelf() - material " << i << "failed to recv itself\n";
	    return res;
	  }
    }
  }

  return res;
}

Vector
DSFMRCPlaneStress::determineTrialStress(Vector strain)
{ 
  // Get Principal strain direction first
  Vector Tstrain(3);     //epslonx,epslony,gammaxy

  // Get strain values from strain of element
  Tstrain(0) = strain(0);
  Tstrain(1) = strain(1);
  Tstrain(2) = strain(2); // gamma

  epsC_vec = Tstrain; // initial epsC vector assumptions

  double epsC1, epsC2, epsSx, epsSy, epsts;
  double fC1, fC1a, fC1b, fC1c, fC2, fSx, fSy;
  double eC1p, eC2p, eSlip1, eSlip2, eC01, eC02, eC1m, eC2m, eT1m, eT2m;
  double vci, vcimax, s_cita, w, deltaS, gammaS;
  double citaS, citaE, temp_cita, citan1, citan2, citaLag, dCitaE, dCitaS, citaIC; // principal strain direction
  double Cs, Cd, betaD;
  double Ecx, Ecy;

  double temp, tolNorm = 1.0e-8;
  double fcr   = 0.65 * pow(-fpc, 0.33); //0.31*sqrt(fpc);      ----ftp
  double Ec    = 2.0 * fpc/epsc0;
  double epscr = fcr/Ec;

  static Matrix Tc(3,3);     // T(cita)
  static Matrix Ts1(3,3);    // T(angle1)
  static Matrix Ts2(3,3);    // T(angle2)
  static Matrix Dc(3,3);     //
  static Matrix Ds1(3,3);    //
  static Matrix Ds2(3,3);    //
  static Matrix Dcp(3,3);    //
  static Matrix Ds1p(3,3);	 //
  static Matrix Ds2p(3,3);	 //

  // Get citaE based on Tstrain, eq. i-9...
  if ( fabs(Tstrain(0)-Tstrain(1)) < 1e-9 ) {
    citaE = 0.25*PI;	
  } else {  // Tstrain(0) != Tstrain(1) 
    temp_cita = 0.5 * atan(fabs(1.0e6*Tstrain(2)/(1.0e6*Tstrain(0)-1.0e6*Tstrain(1)))); 
    if ( fabs(Tstrain(2)) < 1e-7 ) {
      citaE = 0;
    } else if ( (Tstrain(0) > Tstrain(1)) && ( Tstrain(2) > 0) ) {
      citaE = temp_cita;
    } else if ( (Tstrain(0) > Tstrain(1)) && ( Tstrain(2) < 0) ) {
      citaE = PI - temp_cita;
    } else if ( (Tstrain(0) < Tstrain(1)) && ( Tstrain(2) > 0) ) {
      citaE = 0.5*PI - temp_cita;
    } else if ( (Tstrain(0) < Tstrain(1)) && ( Tstrain(2) < 0) ) {
      citaE = 0.5*PI + temp_cita;
    } else {
      opserr << "CSMMRCPlaneStress::determineTrialStress: Failure to calculate citaE\n";
      opserr << " Tstrain(0) = " << Tstrain(0) << endln;
      opserr << " Tstrain(1) = " << Tstrain(1) << endln;
      opserr << " Tstrain(2) = " << Tstrain(2) << endln;
    }
  }
  citaIC = citaE;

  // vCi need to be satisfied several relationships in the concrete crack region
  static Vector errorVec(3); // epsC_vec converge test need a vector to maintain Error
  
  double errNorm;
  bool vCiconverged = false;

  while (vCiconverged == false) {

    // Get epsC1,epsC2 and citaS based on Tstrain, eq.i-10
    epsC1 = 0.5*(epsC_vec(0)+epsC_vec(1)) 
		  + 0.5*sqrt(pow(epsC_vec(0)-epsC_vec(1), 2.0)+pow(epsC_vec(2), 2.0));

    epsC2 = 0.5*(epsC_vec(0)+epsC_vec(1)) 
		  - 0.5*sqrt(pow(epsC_vec(0)-epsC_vec(1), 2.0)+pow(epsC_vec(2), 2.0));
  
    // Get citaS based on epsC_vec, eq.i-11
    if ( fabs(epsC_vec(0)-epsC_vec(1)) < 1e-7 ) {
      citaS = 0.25*PI;	
    } else {  // epsC_vec(0) != epsC_vec(1) 
      temp_cita = 0.5 * atan(fabs(1.0e6*epsC_vec(2)/(1.0e6*epsC_vec(0)-1.0e6*epsC_vec(1)))); 
      if ( fabs(epsC_vec(2)) < 1e-7 ) {
  	    citaS = 0;
  	  } else if ( (epsC_vec(0) > epsC_vec(1)) && ( epsC_vec(2) > 0) ) {
  	    citaS = temp_cita;
  	  } else if ( (epsC_vec(0) > epsC_vec(1)) && ( epsC_vec(2) < 0) ) {
  	    citaS = PI - temp_cita;
  	  } else if ( (epsC_vec(0) < epsC_vec(1)) && ( epsC_vec(2) > 0) ) {
  	    citaS = 0.5*PI - temp_cita;
  	  } else if ( (epsC_vec(0) < epsC_vec(1)) && ( epsC_vec(2) < 0) ) {
  	    citaS = 0.5*PI + temp_cita;
  	  } else {
  	    opserr << "DSFMRCPlaneStress::determineTrialStress: Failure to calculate citaS\n";
  	    opserr << " epsC_vec(0) = " << epsC_vec(0) << endln;
  	    opserr << " epsC_vec(1) = " << epsC_vec(1) << endln;
  	    opserr << " epsC_vec(2) = " << epsC_vec(2) << endln;
  	  }
    }
    
    //cita=citaS;  //eq.i-11    PI/2.0-citaS = citaCrack
    citan1 = citaS - angle1; //angle1-(PI/2.0-citaS); //eq.i-6
    citan2 = citaS - angle2; //angle2-(PI/2.0-citaS);
    
	Tc(0,0) = pow(cos(citaE), 2.0);// T(cita) cita=citaE // 
	Tc(0,1) = pow(sin(citaE), 2.0);
	Tc(0,2) = cos(citaE)*sin(citaE);
	Tc(1,0) = pow(sin(citaE), 2.0);
	Tc(1,1) = pow(cos(citaE), 2.0);
	Tc(1,2) = -cos(citaE)*sin(citaE);
	Tc(2,0) = -2.0*cos(citaE)*sin(citaE);
	Tc(2,1) = 2.0*cos(citaE)*sin(citaE);
	Tc(2,2) = pow(cos(citaE), 2.0)-pow(sin(citaE), 2.0);

	Ts1(0,0) = pow(cos(angle1), 2.0);// T(angle1)
	Ts1(0,1) = pow(sin(angle1), 2.0);
	Ts1(0,2) = cos(angle1)*sin(angle1);
	Ts1(1,0) = pow(sin(angle1), 2.0);
	Ts1(1,1) = pow(cos(angle1), 2.0);
	Ts1(1,2) = -cos(angle1)*sin(angle1);
	Ts1(2,0) = -2.0*cos(angle1)*sin(angle1);
	Ts1(2,1) = 2.0*cos(angle1)*sin(angle1);
	Ts1(2,2) = pow(cos(angle1), 2.0)-pow(sin(angle1), 2.0);

	Ts2(0,0) = pow(cos(angle2), 2.0);// T(angle2)
	Ts2(0,1) = pow(sin(angle2), 2.0);
	Ts2(0,2) = cos(angle2)*sin(angle2);
	Ts2(1,0) = pow(sin(angle2), 2.0);
	Ts2(1,1) = pow(cos(angle2), 2.0);
	Ts2(1,2) = -cos(angle2)*sin(angle2);
	Ts2(2,0) = -2.0*cos(angle2)*sin(angle2);
	Ts2(2,1) = 2.0*cos(angle2)*sin(angle2);
	Ts2(2,2) = pow(cos(angle2), 2.0)-pow(sin(angle2), 2.0);
  
    //if (strain_vec.Norm() < DBL_EPSILON) {
	//  // EbarC1 EbarC2 GbarC, eq.i-5
    //  Dcp.Zero();
    //  Ds1p.Zero();
    //  Ds2p.Zero();
	//  double Ecx = theMaterial[2]->getInitialTangent();
	//  double Ecy = theMaterial[3]->getInitialTangent();
    //  Dcp(0,0) = Ecx;
    //  Dcp(1,1) = Ecy;
    //  Dcp(2,2) = Ecx*Ecy /(Ecx+Ecy);
    //  // Es1_bar Es2_bar, eq.i-5
    //  Ds1p(0,0) = theMaterial[0]->getInitialTangent();
    //  Ds2p(0,0) = theMaterial[1]->getInitialTangent();
    //  
	//  citaS = 0.0;
	//  Tc(0,0) = pow(cos(citaS),2);// T(cita)
	//  Tc(0,1) = pow(sin(citaS),2);
	//  Tc(0,2) = cos(citaS)*sin(citaS);
	//  Tc(1,0) = pow(sin(citaS),2);
	//  Tc(1,1) = pow(cos(citaS),2);
	//  Tc(1,2) = -cos(citaS)*sin(citaS);
	//  Tc(2,0) = -2.0*cos(citaS)*sin(citaS);
	//  Tc(2,1) = 2.0*cos(citaS)*sin(citaS);
	//  Tc(2,2) = pow(cos(citaS),2)-pow(sin(citaS),2);
	//  
    //  Dc.addMatrixTripleProduct(0.0, Tc, Dcp, 1.0);
    //  Ds1.addMatrixTripleProduct(0.0, Ts1, Ds1p, 1.0);
    //  Ds2.addMatrixTripleProduct(0.0, Ts2, Ds2p, 1.0);
    //  
    //  tangent_matrix = Dc + Ds1 + Ds2;
    //  
    //  // Calculate total stress from steel and concrete
    //  Tstress.Zero();
	//  
	//  return Tstress;
    //}
  
    // Get average strains in reinforcement, eq. i-19
    epsSx = 0.5 * (Tstrain(0)+Tstrain(1))+0.5*(Tstrain(0)-Tstrain(1))*cos(2.0*angle1)
  	      + 0.5 * Tstrain(2)*sin(2.0*angle1);  //+epsS0_vec(0)
  
    epsSy = 0.5 * (Tstrain(0)+Tstrain(1))+0.5*(Tstrain(0)-Tstrain(1))*cos(2.0*angle2)
  	      + 0.5 * Tstrain(2)*sin(2.0*angle2);  //+epsS0_vec(1);
  
    // Calculate Steel stress  
    theMaterial[0]->setTrialStrain(epsSx);  // s1
    fSx = theMaterial[0]->getStress();
  
    theMaterial[1]->setTrialStrain(epsSy);  // s2
    fSy = theMaterial[1]->getStress();
        
	// calculate betaD
    Cs = 0.55;
  
    if (epsC1/epsC2 <= -0.28)
  	  Cd = 0.35*pow(-epsC1/epsC2-0.28, 0.8);
    else
  	  Cd = 0.0;
  
    betaD = 1.0/(1.0+Cs*Cd);
    if (betaD >1.0)
  	  betaD = 1.0;
  
    // Calculate average concrete stress
    // C2 --compression for fc2, epsC2
  
    theMaterial[3]->setTrialStrain(epsC2); //need further revision
    fC2 = betaD * theMaterial[3]->getStress(); // w==5 ??
  
    theResponses[3]->getResponse();
    Information &theInfoC03 = theResponses[3]->getInformation();
    epsC12p_nowvec(1) = (*theInfoC03.theVector)(5); // ecp2 
    
    // C1 -- tensile for fc1, epsC1
    //  case 1: eq.i-33,34
    epsts = 2.0 * 7.5e-3 / fcr / 750.; //eq. i-33
    temp  = (epsC1 - epscr) / (epsts - epscr);

    if      (temp <= 0.0)
  	  fC1a = 0.0;
    else if (temp >  1.0)
  	  fC1a = fcr;
    else
  	  fC1a = fcr * (1 - temp);  //eq. i-34
  
    //  case 2: eq. i-35, 36
    temp = 4.0*(roux/db1*abs(cos(citan1))+rouy/db2*abs(cos(citan2)));
    fC1b = fcr/(1+sqrt(2.2/temp*epsC1));
  
    temp = max(fC1a,fC1b); // eq.i-38 /fC1
  
    // eq.i=5   // for positive value
    fC1c = max(roux* (fyx - fSx) * pow(cos(citan1), 2.0), 0.0)
  	     + max(rouy* (fyy - fSy) * pow(cos(citan2), 2.0), 0.0);
    
    fC1 = min(temp, fC1c);
    
	// 
    theMaterial[2]->setTrialStrain(epsC1);
    theResponses[2]->getResponse();
    Information &theInfoC02 = theResponses[2]->getInformation();
    epsC12p_nowvec(0) = (*theInfoC02.theVector)(5); // ecp1 
    if (epsC12p_nowvec(0) > 0) epsC12p_nowvec(0) = 0.0;

    // deltaEspCp(2),  P. Ceresa 2009, Eq. 4
	double deltaEspC1p = epsC12p_nowvec(0) - epsC12p_prevec(0);
	double deltaEspC2p = epsC12p_nowvec(1) - epsC12p_prevec(1);
  
    epsCp_vec(0) += (0.5 * deltaEspC1p * (1+cos(2.0*citaS)) + 0.5 * deltaEspC2p * (1-cos(2.0*citaS)));
    epsCp_vec(1) += (0.5 * deltaEspC1p * (1-cos(2.0*citaS)) + 0.5 * deltaEspC2p * (1+cos(2.0*citaS)));
    epsCp_vec(2) += deltaEspC1p *sin(2.0*citaS) - deltaEspC2p *sin(2.0*citaS);
  
    eC1p = 0.5*(epsCp_vec(0)+epsCp_vec(1)) + 0.5*sqrt(pow(epsCp_vec(0)-epsCp_vec(1),2.0)+pow(epsCp_vec(2),2.0));
    eC2p = 0.5*(epsCp_vec(0)+epsCp_vec(1)) - 0.5*sqrt(pow(epsCp_vec(0)-epsCp_vec(1),2.0)+pow(epsCp_vec(2),2.0));
    
	// Vecchio: Towards Cyclic Load Modeling of Reinforced Concrete
	// C max
	double deltaEspCm1 = (epsC1 > epsC12cm_vec(0) ? 0 : epsC1-epsC12cm_vec(0));
	double deltaEspCm2 = (epsC2 > epsC12cm_vec(1) ? 0 : epsC2-epsC12cm_vec(1));
  
	epsCcm_vec(0) += (0.5 * deltaEspCm1 * (1+cos(2.0*citaS)) + 0.5 * deltaEspCm2 * (1-cos(2.0*citaS)));
	epsCcm_vec(1) += (0.5 * deltaEspCm1 * (1-cos(2.0*citaS)) + 0.5 * deltaEspCm2 * (1+cos(2.0*citaS)));
	epsCcm_vec(2) += deltaEspCm1 *sin(2.0*citaS) - deltaEspCm2 *sin(2.0*citaS);
  
    eC1m = 0.5*(epsCcm_vec(0)+epsCcm_vec(1)) + 0.5*sqrt(pow(epsCcm_vec(0)-epsCcm_vec(1),2.0)+pow(epsCcm_vec(2),2.0));
    eC2m = 0.5*(epsCcm_vec(0)+epsCcm_vec(1)) - 0.5*sqrt(pow(epsCcm_vec(0)-epsCcm_vec(1),2.0)+pow(epsCcm_vec(2),2.0));

	// T max
	double deltaEspTm1 = (epsC1 < epsC12tm_vec(0) ? 0 : epsC1-epsC12tm_vec(0));
	double deltaEspTm2 = (epsC2 < epsC12tm_vec(1) ? 0 : epsC2-epsC12tm_vec(1));
  
	epsCtm_vec(0) += (0.5*deltaEspTm1*(1+cos(2.0*citaS)) + 0.5*deltaEspTm2*(1-cos(2.0*citaS)));
	epsCtm_vec(1) += (0.5*deltaEspTm1*(1-cos(2.0*citaS)) + 0.5*deltaEspTm2*(1+cos(2.0*citaS)));
	epsCtm_vec(2) += deltaEspTm1 *sin(2.0*citaS) - deltaEspTm2 *sin(2.0*citaS);
  
    eT1m = 0.5*(epsCtm_vec(0)+epsCtm_vec(1))+0.5*sqrt(pow(epsCtm_vec(0)-epsCtm_vec(1),2.0)+pow(epsCtm_vec(2),2.0));
	eT2m = 0.5*(epsCtm_vec(0)+epsCtm_vec(1))-0.5*sqrt(pow(epsCtm_vec(0)-epsCtm_vec(1),2.0)+pow(epsCtm_vec(2),2.0));

    // Calculate local stress at cracks
    int status            = 0; // status to check if iteration satisfied eq.i-7
    int iteration_counter = 0;
    double tolerance      = 1.0e-8; // tolerance for iteration
    bool fC1converged     = false;
    double error;
    double fScrx, fScry;
    double epsScrx        = epsSx;
    double epsScry        = epsSy;
    double epsIncr        = 0.0;
  
    while ( fC1converged == false ) {
      // eq.i-20
  	  epsScrx = epsSx + epsIncr * pow(cos(citan1), 2.0);
  	  epsScry = epsSy + epsIncr * pow(cos(citan2), 2.0);
  	  
  	  theMaterial[0]->setTrialStrain(epsScrx);
  	  fScrx = theMaterial[0]->getStress();
  	  
  	  theMaterial[1]->setTrialStrain(epsScry);
  	  fScry = theMaterial[1]->getStress();
  	  
  	  // eq.i-7
  	  temp = roux * (fScrx - fSx) * pow(cos(citan1), 2.0) 
           + rouy * (fScry - fSy) * pow(cos(citan2), 2.0);
  	  
  	  error = fC1 - temp;
  	  
  	  if (abs(error) <= tolerance) {
  	    fC1converged = true;
	  } else {
  	  // temp is function of citan1 citan2 and fSx and fSy,
  	  // the iterative parameter should include citan1/citan2 angle chosen. 
      // like: if (citan1 > PI/2 && citan2 > PI/2)...
  	  // which would improve the converge ratio
  	    temp = 0.618 * (theMaterial[0]->getSecant() + theMaterial[1]->getSecant());
  	    epsIncr += error/temp;
	  }
  	  //opserr << "error for the fc1 and steel in cracks is: " << error << endln;
      //iteration_counter++;
    }
    // eq. i-8
    vci = roux * (fScrx - fSx) * cos(citan1) * sin(citan1)
  	    + rouy * (fScry - fSy) * cos(citan2) * sin(citan2);
    //vci = (vci < 0 ? 0 : vci);

    // Crack slips determine
    // crack spacing // w   
    s_cita = 1.0/(sin(citaS)/xd+cos(citaS)/yd); // eq.i-21,22
    w = epsC1 * s_cita;
  
    // optional 1
    vcimax = sqrt(-fpc)/(0.31+24*w/(aggr+16)); // shear parameter
  
    // for shear slip strain
    if (w <= 0.0) deltaS = 0.0;           //eq. i-40. 12 
    else          deltaS = vci/(1.8*pow(w,-0.8)+(0.234*pow(w,-0.707)-0.2)*(-fpc*1.25)); // fpc change for fcc = 1.25*fpc
    
    double gammaS1  = deltaS/s_cita;
  
    //need further revision for the initial crack direction, citaIC
    if (roux > 0.0 && rouy > 0.0) citaLag = 5./180.*PI;    // eq.i-41~43 
    else                          citaLag = 7.5/180.*PI;
  
	dCitaE = citaE - citaIC;
	dCitaS = (abs(dCitaE) > citaLag ? dCitaE-citaLag : dCitaE);
	citaS  = citaIC + dCitaS;

	double gammaS2 = Tstrain(2)*cos(2.0*citaS) + (Tstrain(1)-Tstrain(0))*sin(2.0*citaS);  // eq. i-18

	gammaS = max(gammaS1, gammaS2);

    epsSlip_vec(0) = -0.5*gammaS*sin(2.0*citaS);   // eq.i-13~15 
    epsSlip_vec(1) =  0.5*gammaS*sin(2.0*citaS);
    epsSlip_vec(2) =  gammaS*cos(2.0*citaS);
  
	// eSlip1,2
    eSlip1 = 0.5*(epsSlip_vec(0)+epsSlip_vec(1))
  	       + 0.5*sqrt(pow(epsSlip_vec(0)-epsSlip_vec(1), 2.0)+pow(epsSlip_vec(2), 2.0));
    eSlip2 = 0.5*(epsSlip_vec(0)+epsSlip_vec(1))
  	       - 0.5*sqrt(pow(epsSlip_vec(0)-epsSlip_vec(1), 2.0)+pow(epsSlip_vec(2), 2.0));
  
	// eC01,2
	eC01 = 0.5*(epsC0_vec(0)+epsC0_vec(1))
		 + 0.5*sqrt(pow(epsC0_vec(0)-epsC0_vec(1), 2.0)+pow(epsC0_vec(2), 2.0));
    eC02 = 0.5*(epsC0_vec(0)+epsC0_vec(1))
		 - 0.5*sqrt(pow(epsC0_vec(0)-epsC0_vec(1), 2.0)+pow(epsC0_vec(2), 2.0));

    // EbarC1 EbarC2 GbarC, eq.i-5
    Dcp.Zero();  // local stiffness matrix
    Ds1p.Zero();
    Ds2p.Zero();
  
    if (epsC1 == 0.0) {
	  Ecx = theMaterial[2]->getSecant();
	} else if (epsC1 < eC1p+eSlip1+eC01) {
      temp = epsC1 - (eC1p+eSlip1+eC01);
	  epsC_vec(0) += 0.5 * temp * (1+cos(2.0*citaS));
      epsC_vec(1) += 0.5 * temp * (1-cos(2.0*citaS));
	  epsC_vec(2) += temp *sin(2.0*citaS);
	  continue;
	} else {
	  //Ecx = fC1/(epsC1-eC1p-eSlip1-eC01);
	  Ecx = fC1/epsC1;
	}
  
    if (epsC2 == 0.0) {
	  Ecy = theMaterial[3]->getSecant();
	} else if (epsC2 > eC2p+eSlip2+eC02) {
	  temp = epsC2 - (eC2p+eSlip2+eC02);
	  epsC_vec(0) += 0.5 * temp * (1-cos(2.0*citaS));
      epsC_vec(1) += 0.5 * temp * (1+cos(2.0*citaS));
	  epsC_vec(2) -= temp *sin(2.0*citaS);
	  continue;
	} else {
	  Ecy = fC2/epsC2;
	}
  
    Dcp(0,0) = Ecx;
    Dcp(1,1) = Ecy;
    Dcp(2,2) = Ecx * Ecy / ( Ecx + Ecy );
  
    // EbarS1, EbarS2, eq.i-5
    
    if (epsSx == 0.0) Ds1p(0,0) = roux * (theMaterial[0]->getSecant());
    else              Ds1p(0,0) = roux * fSx / epsSx;
    
    if (epsSy == 0.0) Ds2p(0,0) = rouy * (theMaterial[1]->getSecant());
    else              Ds2p(0,0) = rouy * fSy / epsSy;
  
    Dc.addMatrixTripleProduct(0.0, Tc, Dcp, 1.0);
    Ds1.addMatrixTripleProduct(0.0, Ts1, Ds1p, 1.0);
    Ds2.addMatrixTripleProduct(0.0, Ts2, Ds2p, 1.0);

	//tangent_matrix = Dc + Ds1 + Ds2;

	//stress0_vec.addMatrixVector(0.0, Dc, epsCp_vec + epsSlip_vec + epsC0_vec, 1.0);
	//stress0_vec.addMatrixVector(1.0, Ds1, epsS0_vec, 1.0); 
	//stress0_vec.addMatrixVector(1.0, Ds2, epsS0_vec, 1.0);  // (3) in vecchio 2001 

	//Tstress.addMatrixVector(0.0, tangent_matrix, strain_vec, 1.0);

	//tempVec = Tstress + stress0_vec;

	errorVec = epsC_vec - (Tstrain - epsCp_vec - epsSlip_vec - epsC0_vec);

    errNorm = errorVec.pNorm(-1);
	// determine the norm of matrixnorm = tempMat.Norm();
	// errNorm = norm - tempNorm;
	
  	if ( errNorm <= tolNorm) {
	  vCiconverged = true;
	} else {
	  //tempNorm = tempMat.Norm();
	  epsC_vec = Tstrain - epsCp_vec - epsSlip_vec - epsC0_vec;
	  epsC12p_prevec = epsC12p_nowvec;
	  iteration_counter += 1 ;
	  if (iteration_counter > 65) {
		iteration_counter = 0;
		epsC_vec = Tstrain - 0.618 * (epsCp_vec + epsSlip_vec + epsC0_vec);
	  }
	}
  }

  tangent_matrix = Dc + Ds1 + Ds2;

  // Calculate total stress from steel and concrete
  Tstress.addMatrixVector(0.0, tangent_matrix, strain_vec, 1.0);

  stress0_vec.addMatrixVector(0.0, Dc, epsCp_vec + epsSlip_vec + epsC0_vec, 1.0);
  //stress0_vec.addMatrixVector(1.0, Ds1, epsS0_vec, 1.0); 
  //stress0_vec.addMatrixVector(1.0, Ds2, epsS0_vec, 1.0);  // (3) in vecchio 2001 

  // determine internal vectors

  return Tstress;
}

double
DSFMRCPlaneStress::kupferenvelop(double Tstrain, double sig_p, double eps_p)
{
  double sig = 0.0;
  if (Tstrain > eps_p) {
    double eta = Tstrain/eps_p;
    sig = sig_p * (2 * eta - eta * eta);
    //double Ec0 = 2.0*sig_p/eps_p;
    //Ttangent = Ec0*(1.0-eta);
  }
  else if (Tstrain > 2.0 * epsc0) {
	double eta = (Tstrain-eps_p)/(2.0*epsc0-eps_p);
    //Ttangent = (sig_p-fpc)/(eps_p-epscu);
    sig = sig_p *(1.0-eta*eta);
  }
  else {
    sig = 0.2*fpc;
    //Ttangent = 0.0;
  }
  return sig;
}

int
DSFMRCPlaneStress::determineTangent()
{
  Vector strain = strain_vec;
  tangent_matrix.Zero();
  if (strain_vec.Norm() < DBL_EPSILON) {
	// Dc
	double Ec0  = 2.0*fpc/epsc0;
	double Gc   = Ec0/6.0; 
	tangent_matrix(0,0) = Ec0+roux*E0;
	tangent_matrix(1,1) = Ec0+rouy*E0;
	tangent_matrix(2,2) = Gc;
  } else {
    // Dc
	for (int i = 0; i<3; i++) {
	  double temp = strain(i);
	  double deltaE = DBL_EPSILON * 1.e10;
	  strain(i) = temp + deltaE;
	  Vector stress1 = determineTrialStress(strain);
	  strain(i) = temp - deltaE;
	  Vector stress2 = determineTrialStress(strain);

	  for (int j = 0; j<3; j++) {
	    tangent_matrix(i,j) = (stress1(j)-stress2(j))/(2.0*deltaE);
	  }
    }
    
    // Ds
    if (fabs(strain(0))<0.02) {
	  tangent_matrix(0,0) = tangent_matrix(0,0)+ roux*E0;
    } 
    if (fabs(strain(1))<0.02) {
	  tangent_matrix(1,1) = tangent_matrix(1,1)+ rouy*E0;
    } 
  }

  return 0;
}