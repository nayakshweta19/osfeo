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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/MCFTRCPlaneStress02.cpp,v $
// File: MCFTRCPlaneStress02.cpp
//
// Written: Lining
// Created: 2011.8
//
// Description: This file contains the class definition for
// MCFTRCPlaneStress02
// For Detailed explanation of the model, please refer to the journal paper
// entitled "The modified compression-filed theory for reinforced concrete element
// subjected to shear, ACI Journal. s83-22. 1986. pp. 219-231"


#include "MCFTRCPlaneStress02.h"
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

static int numMCFTRCPlaneStress02Materials = 0;

OPS_Export void *
OPS_NewMCFTRCPlaneStress02Material()
{
  if (numMCFTRCPlaneStress02Materials == 0) {
    numMCFTRCPlaneStress02Materials++;
    //OPS_Error("MCFTRCPlaneStress02 material - Written by Lining - Copyright@2011\n", 1);
  }

  // Pointer to a uniaxial material that will be returned
  NDMaterial *theMaterial = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 16) {
    opserr << "Invalid Args want: NDMaterial MCFTRCPlaneStress02 matTag? rho? UniaxiaMatTag1? UniaxiaMatTag2? UniaxiaMatTag3? UniaxiaMatTag4? \
               angle1? angle2? rou1? rou2? db1? db2? fpc? fy? E0? epsc0? aggsize? xd? yd?\n";
    return 0;	
  }

  int tag;
  double rho;
  int    iData[4];
  double dData[15];
  int numData = 0;

  numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid nDMaterial MCFTRCPlaneStress02 tag" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDouble(&numData, &rho) != 0) {
    opserr << "WARNING invalid parameter rho MCFTRCPlaneStress02 tag:" << tag << endln;
    return 0;	
  }

  numData = 4;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid nDMaterial MCFTRCPlaneStress02 tag: " << tag << endln;
    return 0;
  }

  numData = 13;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data MCFTRCPlaneStress02 tag: " << tag << endln;
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial1 = OPS_GetUniaxialMaterial(iData[0]);
    
  if (theUniaxialMaterial1 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[0];
    opserr << "\nMCFTRCPlaneStress02 tag: " << tag << endln;
    return 0;
  }
  
  UniaxialMaterial *theUniaxialMaterial2 = OPS_GetUniaxialMaterial(iData[1]);

  if (theUniaxialMaterial2 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[1];
    opserr << "\nMCFTRCPlaneStress02 tag: " << tag << endln;
    return 0;
  }
  
  UniaxialMaterial *theUniaxialMaterial3 = OPS_GetUniaxialMaterial(iData[2]);
  if (theUniaxialMaterial3 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[2];
    opserr << "\nMCFTRCPlaneStress02 tag: " << tag << endln;
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial4 = OPS_GetUniaxialMaterial(iData[3]);  
  if (theUniaxialMaterial4 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[3];
    opserr << "\nMCFTRCPlaneStress02 tag: " << tag << endln;
    return 0;
  }

  //now create the MCFTRCPlaneStress02
  theMaterial = new MCFTRCPlaneStress02 (tag, rho,
						     theUniaxialMaterial1, theUniaxialMaterial2, 
						     theUniaxialMaterial3, theUniaxialMaterial4,
						     dData[0], dData[1], dData[2], dData[3], dData[4],
							 dData[5], dData[6], dData[7], dData[8], dData[9],
							 dData[10],dData[11],dData[12]);
       
  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory creating material\n";
    opserr << "MCFTRCPlaneStress02: " << tag << endln;
    return 0;
  }

  return theMaterial;
}
 
MCFTRCPlaneStress02 ::MCFTRCPlaneStress02 (int tag, 
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
								   double   FY,
								   double   E,
								   double   EPSC0,
								   double   AGGR,
								   double   XD,
								   double   YD):
  NDMaterial(tag, ND_TAG_MCFTRCPlaneStress02), 
  rho(RHO),angle1(ANGLE1),angle2(ANGLE2),rou1(ROU1),rou2(ROU2),db1(DB1),db2(DB2),
  fpc(FPC), fy(FY), E0(E), epsc0(EPSC0), aggr(AGGR), xd(XD), yd(YD),
  lastStress(3),Tstress(3),
  strain_vec(3),strainC_vec(3),strainSlip_vec(3),strainC0_vec(3),strainS0_vec(3),
  strainCp_vec(3),
  stress_vec(3),tangent_matrix(3,3)
{
    steelStatus = 0;
    dirStatus = 0;
    G12 = 0;
    citaStrain = 10;
    citaStress = 10;
    
    TOneReverseStatus = 0;         
    TOneNowMaxComStrain = 0.0;
    TOneLastMaxComStrain = 0.0;
    
    TTwoReverseStatus = 0;         
    TTwoNowMaxComStrain = 0.0;
    TTwoLastMaxComStrain = 0.0;
    
    COneReverseStatus = 0;         
    COneNowMaxComStrain = 0.0;
    COneLastMaxComStrain = 0.0;
    
    CTwoReverseStatus = 0;         
    CTwoNowMaxComStrain = 0.0;
    CTwoLastMaxComStrain = 0.0;
    
    lastStress.Zero();  // add at 7.28    
    
    if ( fpc > 0.0 ) { fpc = -fpc; } // set fpc < 0
    if ( epsc0 > 0.0 ) { epsc0 = -epsc0; } // set fpc < 0
    theMaterial = 0;
    
    // Allocate pointers to theSteel1
    theMaterial = new UniaxialMaterial *[4];
    
    if ( theMaterial == 0 ) {
      opserr << " MCFTRCPlaneStress02::MCFTRCPlaneStress02 - failed allocate material array\n";
      exit(-1);
    }
    
    // Get the copy for the Steel1
    theMaterial[0] = s1->getCopy();
    // Check allocation    
    if ( theMaterial[0] == 0 ) {
      opserr << " MCFTRCPlaneStress02::MCFTRCPlaneStress02 - failed to get a copy for steel1\n";
      exit(-1);
    }
    
    // Get the copy for the Steel2
    theMaterial[1] = s2->getCopy();	
    // Check allocation    
    if ( theMaterial[1] == 0 ) {
      opserr << " MCFTRCPlaneStress02::MCFTRCPlaneStress02 - failed to get a copy for steel2\n";
      exit(-1);
    }
    
    // Get the copy for the Concrete1
    theMaterial[2] = c1->getCopy();	
    // Check allocation    
    if ( theMaterial[2] == 0 ) {
      opserr << " MCFTRCPlaneStress02::MCFTRCPlaneStress02 - failed to get a copy for concrete1\n";
      exit(-1);
    }
    
    // Get the copy for the Concrete2
    theMaterial[3] = c2->getCopy();	
    // Check allocation    
    if ( theMaterial[3] == 0 ) {
      opserr << " MCFTRCPlaneStress02::MCFTRCPlaneStress02 - failed to get a copy for concrete2\n";
      exit(-1);
    }
    
    /* FMK */
    theResponses = new Response *[6];  
    
    if ( theResponses == 0) {
      opserr << " MCFTRCPlaneStress02::MCFTRCPlaneStress02 - failed allocate responses array\n";
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
      opserr << " MCFTRCPlaneStress02::MCFTRCPlaneStress02 - failed to set appropriate materials tag: " << tag << "\n";
      exit(-1);
    }
    
    delete theDummyStream;
	strain_vec.Zero();
    determineTrialStress(strain_vec);
    this->revertToStart();
}

MCFTRCPlaneStress02::MCFTRCPlaneStress02()
 :NDMaterial(0, ND_TAG_MCFTRCPlaneStress02),lastStress(3),Tstress(3),strain_vec(3),
  strainC_vec(3),strainSlip_vec(3),strainC0_vec(3),strainS0_vec(3),strainCp_vec(3),
  stress_vec(3),tangent_matrix(3,3)
{
  theMaterial = 0;
  theResponses = 0;
  strain_vec.Zero();
  determineTrialStress(strain_vec);
  this->revertToStart();
}

MCFTRCPlaneStress02::~MCFTRCPlaneStress02()
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
MCFTRCPlaneStress02::getRho(void)
{
	return rho;
}

// really used one
int 
MCFTRCPlaneStress02::setTrialStrain(const Vector &v)
{

  // Set values for strain_vec
  strain_vec = v;
  
  // Set initial values for Tstress
  Tstress.Zero();
  
  TOneReverseStatus = COneReverseStatus;         
  TOneNowMaxComStrain = COneNowMaxComStrain;
  TOneLastMaxComStrain = COneLastMaxComStrain;
  
  TTwoReverseStatus = CTwoReverseStatus;         
  TTwoNowMaxComStrain = CTwoNowMaxComStrain;
  TTwoLastMaxComStrain = CTwoLastMaxComStrain;
  
  Tstress = determineTrialStress(strain_vec);
  
  //determineTangent();

  return 0;
}

int 
MCFTRCPlaneStress02::setTrialStrain(const Vector &v, const Vector &r)
{
  opserr << "error:MCFTRCPlaneStress02::setTrialStrain(&v, &r) -- not really responsibility" << endln;
  return 0;
}

int 
MCFTRCPlaneStress02::setTrialStrainIncr(const Vector &v)
{
  opserr << "error:MCFTRCPlaneStress02::setTrialStrainIncr(&v) -- not really responsibility" << endln;
  return 0;
}

int
MCFTRCPlaneStress02::setTrialStrainIncr(const Vector &v, const Vector &r)
{
  opserr << "error:MCFTRCPlaneStress02::setTrialStrainIncr(&v, &r) -- not really responsibility" << endln;
  return 0;
}

const Matrix& 
MCFTRCPlaneStress02::getTangent(void)
{

  determineTrialStress(strain_vec);

  //determineTangent();

  return tangent_matrix;
}

const Vector& 
MCFTRCPlaneStress02::getStress(void)
{
  return stress_vec;
}

const Vector& 
MCFTRCPlaneStress02 :: getStrain()
{
  return strain_vec;
}

const Vector& 
MCFTRCPlaneStress02::getCommittedStress(void)
{
  return stress_vec;
}

const Vector& 
MCFTRCPlaneStress02::getCommittedStrain(void)
{
    return strain_vec;
}

int
MCFTRCPlaneStress02::commitState(void)
{
  for (int i=0; i<4; i++) {
    theMaterial[i]->commitState();
  }
  
  COneReverseStatus = TOneReverseStatus;         
  COneNowMaxComStrain = TOneNowMaxComStrain;
  COneLastMaxComStrain = TOneLastMaxComStrain;
  
  CTwoReverseStatus = TTwoReverseStatus;         
  CTwoNowMaxComStrain = TTwoNowMaxComStrain;
  CTwoLastMaxComStrain = TTwoLastMaxComStrain;
  
  lastStress = stress_vec;
  
  return 0;
}

int
MCFTRCPlaneStress02::revertToLastCommit(void)
{
  for (int i=0; i<4; i++) {
    theMaterial[i]->revertToLastCommit();
  }
  
  TOneReverseStatus = COneReverseStatus;         
  TOneNowMaxComStrain = COneNowMaxComStrain;
  TOneLastMaxComStrain = COneLastMaxComStrain;
  
  TTwoReverseStatus = CTwoReverseStatus;         
  TTwoNowMaxComStrain = CTwoNowMaxComStrain;
  TTwoLastMaxComStrain = CTwoLastMaxComStrain;
  
  return 0;
}

int
MCFTRCPlaneStress02::revertToStart(void)
{
  for (int i=0; i<4; i++) {
    theMaterial[i]->revertToStart();
  }
  
  Tstress.Zero();
  strain_vec.Zero();
  strainC_vec.Zero();
  strainSlip_vec.Zero();
  strainC0_vec.Zero();
  strainS0_vec.Zero();
  strainCp_vec.Zero();
  stress_vec.Zero();
  
  steelStatus = 0;
  dirStatus = 0;
  G12 = 0;
  
  TOneReverseStatus = 0;         
  TOneNowMaxComStrain = 0.0;
  TOneLastMaxComStrain = 0.0;
  
  TTwoReverseStatus = 0;         
  TTwoNowMaxComStrain = 0.0;
  TTwoLastMaxComStrain = 0.0;
  
  COneReverseStatus = 0;         
  COneNowMaxComStrain = 0.0;
  COneLastMaxComStrain = 0.0;
  
  CTwoReverseStatus = 0;         
  CTwoNowMaxComStrain = 0.0;
  CTwoLastMaxComStrain = 0.0;
  
  return 0;
}

NDMaterial* MCFTRCPlaneStress02::getCopy(void)
{
  MCFTRCPlaneStress02* theCopy =
    new MCFTRCPlaneStress02( this->getTag(), 
					 rho, 
					 theMaterial[0], 
					 theMaterial[1], 
					 theMaterial[2], 
					 theMaterial[3], 
					 angle1, angle2, 
					 rou1, rou2, 
					 db1, db2,
					 fpc, 
					 fy, 
					 E0, 
					 epsc0,
					 aggr,
					 xd, yd);
  theCopy->strain_vec = strain_vec;
  return theCopy;
}

NDMaterial* MCFTRCPlaneStress02::getCopy(const char *type)
{
	MCFTRCPlaneStress02* theModel =
		new MCFTRCPlaneStress02( this->getTag(), 
           rho, theMaterial[0], theMaterial[1], theMaterial[2], theMaterial[3],
		   angle1, angle2, rou1, rou2, db1, db2, fpc, fy, E0, epsc0, aggr, xd, yd );
	theModel->strain_vec = strain_vec;
	theModel->stress_vec = stress_vec;
	return theModel;
}

//added by Ln
Response*
MCFTRCPlaneStress02::setResponse (const char **argv, int argc, OPS_Stream &output)
{

#ifdef DEBUG
	opserr << "MCFTRCPlaneStress02::setResponse(...)" << endln;
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
MCFTRCPlaneStress02::getResponse (int responseID, Information &matInfo)
{
#ifdef DEBUG
	opserr << "MCFTRCPlaneStress02::getResponse(...)" << endln;
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
MCFTRCPlaneStress02::Print(OPS_Stream &s, int flag )
{
	s << "\n\tMCFTRCPlaneStress02, material id: " << this->getTag() << endln;
    /*
	s << "\tRho: " << rho << endln;
	s << "\tangle1: " << angle1 << endln;
	s << "\tangle2: " << angle2 << endln;
	s << "\trou1: " << rou1 << endln;
	s << "\trou2: " << rou2 << endln;
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
MCFTRCPlaneStress02::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  int dataTag = this->getDbTag();
  
  // Packs its data into a Vector and sends this to theChannel
  static Vector data(15);
  data(0) = this->getTag();
  data(1) = rho;
  data(2) = angle1;
  data(3) = angle2;
  data(4) = rou1;
  data(5) = rou2;
  data(6) = db1;
  data(7) = db2;
  data(8) = fpc;
  data(9) = fy;
  data(10) = E0;
  data(11) = epsc0;
  data(12) = aggr;
  data(13) = xd;
  data(14) = yd;

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING MCFTRCPlaneStress02::sendSelf() - " << this->getTag() << " failed to send Vector\n";
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
    opserr << "WARNING MCFTRCPlaneStress02::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }
  
  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 4; i++) {
    res += theMaterial[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "MCFTRCPlaneStress02::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }	
  
  return res;
}

int
MCFTRCPlaneStress02::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();
  
  // Quad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(15);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING MCFTRCPlaneStress02::recvSelf() - failed to receive Vector\n";
    return res;
  }
  
  this->setTag((int)data(0));
  rho = data(1);
  angle1 = data(2);
  angle2 = data(3);
  rou1   = data(4);
  rou2   = data(5);
  db1    = data(6);
  db2    = data(7);
  fpc    = data(8);
  fy     = data(9);
  E0     = data(10);
  epsc0  = data(11);
  aggr   = data(12);
  xd     = data(13);
  yd     = data(14);

  static ID idData(8);
  
  // now receives the tags of its materials
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING MCFTRCPlaneStress02::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  if (theMaterial == 0) {
    // Allocate new materials
    theMaterial = new UniaxialMaterial *[4];
    if (theMaterial == 0) {
      opserr << "MCFTRCPlaneStress02::recvSelf() - Could not allocate UniaxialMaterial* array\n";
      return -1;
    }
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Allocate new material with the sent class tag
      theMaterial[i] = theBroker.getNewUniaxialMaterial(matClassTag);
      if (theMaterial[i] == 0) {
	     opserr << "MCFTRCPlaneStress02::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
	     return -1;
      }
      // Now receive materials into the newly allocated space
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
         opserr << "MCFTRCPlaneStress02::recvSelf() - material " << i << "failed to recv itself\n";
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
            opserr << "MCFTRCPlaneStress02::recvSelf() - material " << i << "failed to create\n";
	        return -1;
		 }
	  }
      // Receive the material
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "MCFTRCPlaneStress02::recvSelf() - material " << i << "failed to recv itself\n";
	    return res;
	  }
    }
  }

  return res;
}

Vector
MCFTRCPlaneStress02::determineTrialStress(Vector strain)
{ 
  // Get Principal strain direction first
  
  Vector Tstrain(3);     //epslonx,epslony,0.5*gammaxy

  // Get strain values from strain of element
  Tstrain(0) = strain(0);
  Tstrain(1) = strain(1);
  Tstrain(2) = 0.5*strain(2);
  
  double temp;
  double epsC1, epsC2, citaS, temp_citaS, citaE, temp_citaE; // principal strain direction
  
  // Get epsC1,epsC2 and citaS based on Tstrain, eq.i-10

  epsC1 = 0.5*(Tstrain(0)+Tstrain(1))+0.5*sqrt(pow(Tstrain(0)-Tstrain(1),2.)+pow(Tstrain(2),2.));
  epsC2 = 0.5*(Tstrain(0)+Tstrain(1))-0.5*sqrt(pow(Tstrain(0)-Tstrain(1),2.)+pow(Tstrain(2),2.));

  // Get citaS based on strainC_vec, eq.i-11
  if ( fabs(strainC_vec(0)-strainC_vec(1)) < 1e-7 ) {
    citaS = 0.25*PI;	
  } else {  // strainC_vec(0) != strainC_vec(1) 
    temp_citaS = 0.5 * atan(fabs(2.0e6*strainC_vec(2)/(1.0e6*strainC_vec(0)-1.0e6*strainC_vec(1)))); 
    if ( fabs(strainC_vec(2)) < 1e-7 ) {
	  citaS = 0;
	} else if ( (strainC_vec(0) > strainC_vec(1)) && ( strainC_vec(2) > 0) ) {
	  citaS = temp_citaS;
	} else if ( (strainC_vec(0) > strainC_vec(1)) && ( strainC_vec(2) < 0) ) {
	  citaS = PI - temp_citaS;
	} else if ( (strainC_vec(0) < strainC_vec(1)) && ( strainC_vec(2) > 0) ) {
	  citaS = 0.5*PI - temp_citaS;
	} else if ( (strainC_vec(0) < strainC_vec(1)) && ( strainC_vec(2) < 0) ) {
	  citaS = 0.5*PI + temp_citaS;
	} else {
	  opserr << "MCFTRCPlaneStress02::determineTrialStress: Failure to calculate citaS\n";
	  opserr << " strainC_vec(0) = " << strainC_vec(0) << endln;
	  opserr << " strainC_vec(1) = " << strainC_vec(1) << endln;
	  opserr << " strainC_vec(2) = " << strainC_vec(2) << endln;
	}
  }
  
  //cita=citaS;  //eq.i-11
  double citan1 = citaS-angle1; //eq.i-6
  double citan2 = citaS-angle2;

  // Get citaE based on Tstrain, eq. i-9
  if ( fabs(Tstrain(0)-Tstrain(1)) < 1e-7 ) {
    citaE = 0.25*PI;	
  } else {  // Tstrain(0) != Tstrain(1) 
    temp_citaE = 0.5 * atan(fabs(2.0e6*Tstrain(2)/(1.0e6*Tstrain(0)-1.0e6*Tstrain(1)))); 
    if ( fabs(Tstrain(2)) < 1e-7 ) {
	  citaE = 0;
	} else if ( (Tstrain(0) > Tstrain(1)) && ( Tstrain(2) > 0) ) {
	  citaE = temp_citaE;
	} else if ( (Tstrain(0) > Tstrain(1)) && ( Tstrain(2) < 0) ) {
	  citaE = PI - temp_citaE;
	} else if ( (Tstrain(0) < Tstrain(1)) && ( Tstrain(2) > 0) ) {
	  citaE = 0.5*PI - temp_citaE;
	} else if ( (Tstrain(0) < Tstrain(1)) && ( Tstrain(2) < 0) ) {
	  citaE = 0.5*PI + temp_citaE;
	} else {
	  opserr << "CSMMRCPlaneStress::determineTrialStress: Failure to calculate citaE\n";
	  opserr << " Tstrain(0) = " << Tstrain(0) << endln;
	  opserr << " Tstrain(1) = " << Tstrain(1) << endln;
	  opserr << " Tstrain(2) = " << Tstrain(2) << endln;
	}
  }

  static Matrix Tc(3,3);     // T(cita)
  static Matrix Ts1(3,3);    // T(angle1)
  static Matrix Ts2(3,3);    // T(angle2)
  static Matrix Dc(3,3);     //
  static Matrix Ds1(3,3);    //
  static Matrix Ds2(3,3);    //
  static Matrix Dcp(3,3);    //
  static Matrix Ds1p(3,3);	 //
  static Matrix Ds2p(3,3);	 //

  Tc(0,0) = pow(cos(citaS),2);// T(cita)
  Tc(0,1) = pow(sin(citaS),2);
  Tc(0,2) = cos(citaS)*sin(citaS);
  Tc(1,0) = pow(sin(citaS),2);
  Tc(1,1) = pow(cos(citaS),2);
  Tc(1,2) = -cos(citaS)*sin(citaS);
  Tc(2,0) = -2.0*cos(citaS)*sin(citaS);
  Tc(2,1) = 2.0*cos(citaS)*sin(citaS);
  Tc(2,2) = pow(cos(citaS),2)-pow(sin(citaS),2);

  Ts1(0,0) = pow(cos(angle1),2);// T(angle1)
  Ts1(0,1) = pow(sin(angle1),2);
  Ts1(0,2) = cos(angle1)*sin(angle1);
  Ts1(1,0) = pow(sin(angle1),2);
  Ts1(1,1) = pow(cos(angle1),2);
  Ts1(1,2) = -cos(angle1)*sin(angle1);
  Ts1(2,0) = -2.0*cos(angle1)*sin(angle1);
  Ts1(2,1) = 2.0*cos(angle1)*sin(angle1);
  Ts1(2,2) = pow(cos(angle1),2)-pow(sin(angle1),2);

  Ts2(0,0) = pow(cos(angle2),2);// T(angle2)
  Ts2(0,1) = pow(sin(angle2),2);
  Ts2(0,2) = cos(angle2)*sin(angle2);
  Ts2(1,0) = pow(sin(angle2),2);
  Ts2(1,1) = pow(cos(angle2),2);
  Ts2(1,2) = -cos(angle2)*sin(angle2);
  Ts2(2,0) = -2.0*cos(angle2)*sin(angle2);
  Ts2(2,1) = 2.0*cos(angle2)*sin(angle2);
  Ts2(2,2) = pow(cos(angle2),2)-pow(sin(angle2),2);

  if (strain_vec.Norm() < DBL_EPSILON) {
	// EbarC1 EbarC2 GbarC, eq.i-5
    Dcp.Zero();
    Ds1p.Zero();
    Ds2p.Zero();
	double Ecx = theMaterial[2]->getInitialTangent();
	double Ecy = theMaterial[3]->getInitialTangent();
    Dcp(0,0) = Ecx;
    Dcp(1,1) = Ecy;
    Dcp(2,2) = Ecx*Ecy /(Ecx+Ecy);
    // Es1_bar Es2_bar, eq.i-5
    Ds1p(0,0) = theMaterial[0]->getInitialTangent();
    Ds2p(1,1) = theMaterial[1]->getInitialTangent();
        
    Dc.addMatrixTripleProduct(0.0, Tc, Dcp, 1.0);
    Ds1.addMatrixTripleProduct(0.0, Ts1, Ds1p, 1.0);
    Ds2.addMatrixTripleProduct(0.0, Ts2, Ds2p, 1.0);
    
    tangent_matrix = Dc + Ds1 + Ds2;
    
    // Calculate total stress from steel and concrete
    Tstress.Zero();

	return Tstress;
  }

  // Get average strains in reinforcement, eq. i-19
  double epsS1 = 0.5*(Tstrain(0)+Tstrain(1))+0.5*(Tstrain(0)-Tstrain(1))*cos(2.0*angle1)
	    + 0.5*Tstrain(2)*sin(angle1*2.0)+strainS0_vec(0);

  double epsS2 = 0.5*(Tstrain(0)+Tstrain(1))+0.5*(Tstrain(0)-Tstrain(1))*cos(2.0*angle2)
	    + 0.5*Tstrain(2)*sin(angle2*2.0)+strainS0_vec(1);

  // calculate betaD
  double Cs = 0.55, Cd, betaD;
  temp =-epsC1/epsC2-0.28;
  if (temp >= 0.0) {
	Cd = 0.35*pow(temp, 0.8);
  } else {
	Cd = 0.0;
  }
  betaD = 1.0/(1.0+Cs*Cd);
  if (betaD >1.0) betaD = 1.0;

  // Calculate average concrete stress
  // C2 --compression for fc2, epsC2

  theMaterial[3]->setTrialStrain(betaD*epsC2); //need further revision
  double fc2 = theMaterial[3]->getStress(); // w==5 ??

  // C1 -- tensile for fc1, epsC1
  //  case 1: eq.i-33,34

  double fcr = 0.65*pow(-fpc,0.33); //0.31*sqrt(fpc);      ----ftp
  double Ec  = 2.0*fpc/epsc0;
  double epscr = fcr/Ec;
  double epsts = 2.0*75.0e-6/fcr/0.5; //eq. i-33
  temp = (epsC1-epscr)/(epsts-epscr);

  double fc1, fc1a, fc1b;

  if (temp <= 0.0) {
	fc1a = 0.0;
  } else if (temp>1.0) {
	fc1a = fcr;
  } else {
	fc1a = fcr * (1 - temp);
  }

  //  case 2: eq. i-35, 36
  temp = 4.0*(rou1/db1*abs(cos(angle1))+rou2/db2*abs(cos(angle2)));
  fc1b = fcr/(1+sqrt(2.2/temp*epsC1));

  fc1 = max(fc1a,fc1b); // eq.i-38

  // Calculate Steel stress  
  // s1
  theMaterial[0]->setTrialStrain(epsS1);
  double fs1 = theMaterial[0]->getStress();
  
  // s2
  theMaterial[1]->setTrialStrain(epsS2);
  double fs2 = theMaterial[1]->getStress();
  
  // Calculate local stress at cracks
  int status            = 0; // status to check if iteration satisfied eq.i-7
  int iteration_counter = 0;
  bool converged        = false;
  double tolerance      = 1.0e-8; // tolerance for iteration
  double error;
  double fscr1, fscr2; 
  double epsScr1 = epsS1;
  double epsScr2 = epsS2;
  double epsIncr = 0.0;

  while ( converged == false ) {
    // eq.i-20
	epsScr1 = epsS1 + epsIncr*pow(cos(citan1),2.);
	epsScr2 = epsS2 + epsIncr*pow(cos(citan2),2.);

	theMaterial[0]->setTrialStrain(epsScr1);
	fscr1 = theMaterial[0]->getStress();

	theMaterial[1]->setTrialStrain(epsScr2);
	fscr2 = theMaterial[1]->getStress();
	
	// eq.i-7
	temp = rou1 * (fscr1 - fs1) * pow(cos(citan1), 2.0) 
		  +rou2 * (fscr2 - fs2) * pow(cos(citan2), 2.0);
	error = fc1 - temp;
	if (abs(error) <= tolerance)
	  converged = true;

	temp = 0.618 * (theMaterial[0]->getTangent() + theMaterial[1]->getTangent());
	epsIncr += error/temp;
          
    iteration_counter++;
  }

  double vcimax = 0.0; // shear parameter
  // eq. i-8
  double vci = rou1 * (fscr1 - fs1) * cos(citan1) *sin(citan1)
	          +rou2 * (fscr2 - fs2) * cos(citan2) *sin(citan2);

  // Crack slips determine
  //crack spacing // w   //eq.i-21,22
  double s_cita = 1.0/(sin(citaS)/xd+cos(citaS)/yd);
  double w = epsC1 * s_cita;

  //eq. i-40. 12
  double deltas, gammaS;
  if (w <= 0.0) {
	deltas = 0.0;
  } else {
    deltas = vci/(1.8*pow(w,-0.8)+(0.234*pow(w,-0.707)-0.2)*fpc*1.25); // fpc change for fcc = 1.25*fpc
  }
  gammaS  = deltas/s_cita;

  // eq.i-41~43
  double citaLag;
  if (rou1 > 0.0 && rou2 > 0.0) {
	citaLag = 5./180.*PI;
  } else {
	citaLag = 7.5/180.*PI;
  }

  // eq.i-13~15
  double epsSlipx = -0.5*gammaS*sin(2.0*citaS);
  double epsSlipy = 0.5*gammaS*sin(2.0*citaS);
  double gammaSlipxy = gammaS*cos(2.0*citaS);

  // EbarC1 EbarC2 GbarC, eq.i-5
  Dcp.Zero();
  Ds1p.Zero();
  Ds2p.Zero();

  double Ecx, Ecy;
  if (epsC1 == 0.0) Ecx = theMaterial[2]->getTangent();
  else Ecx = fc1/epsC1;

  if (epsC2 == 0.0) Ecy = theMaterial[3]->getTangent();
  else Ecy = fc2/epsC2;

  Dcp(0,0) = Ecx;
  Dcp(1,1) = Ecy;
  Dcp(2,2) = Ecx * Ecy / ( Ecx + Ecy );

  // Es1_bar Es2_bar, eq.i-5
  // double Esx, Esy;
  if (epsS1 == 0.0) Ds1p(0,0) = rou1 * (theMaterial[0]->getTangent());
  else Ds1p(0,0) = rou1*fs1/epsS1;
  
  if (epsS2 == 0.0) Ds2p(1,1) = rou2 * (theMaterial[1]->getTangent());
  else Ds2p(1,1) = rou2*fs2/epsS2;

  Dc.addMatrixTripleProduct(0.0, Tc, Dcp, 1.0);
  Ds1.addMatrixTripleProduct(0.0, Ts1, Ds1p, 1.0);
  Ds2.addMatrixTripleProduct(0.0, Ts2, Ds2p, 1.0);

  tangent_matrix = Dc + Ds1 + Ds2;

  // Calculate total stress from steel and concrete
  Tstress.addMatrixVector(0.0,tangent_matrix,strain_vec,1.0);

  // determine internal vectors

  return Tstress;
}

double
MCFTRCPlaneStress02::kupferEnvelop(double Tstrain, double sig_p, double eps_p)
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
MCFTRCPlaneStress02::determineTangent()
{
  Vector strain = strain_vec;
  tangent_matrix.Zero();
  if (strain_vec.Norm() < DBL_EPSILON) {
	// Dc
	double Ec0  = 2.0*fpc/epsc0;
	double Gc   = Ec0/6.0; 
	tangent_matrix(0,0) = Ec0+rou1*E0;
	tangent_matrix(1,1) = Ec0+rou2*E0;
	tangent_matrix(2,2) = Gc;
  } else {
    // Dc
	for (int i = 0; i<3; i++) {
	  double temp = strain(i);
	  double deltaE = DBL_EPSILON * 1.e10;
	  strain(i) = temp + deltaE;
	  Vector stress1 = determineTrialStress(strain);
	  strain(i) = temp - deltaE * 2.0;
	  Vector stress2 = determineTrialStress(strain);

	  for (int j = 0; j<3; j++) {
	    tangent_matrix(i,j) = (stress1(j)-stress2(j))/(2.0*deltaE);
	  }
    }
    
    // Ds
    if (fabs(strain(0))<0.02) {
	  tangent_matrix(0,0) = tangent_matrix(0,0)+ rou1*E0;
    } 
    if (fabs(strain(1))<0.02) {
	  tangent_matrix(1,1) = tangent_matrix(1,1)+ rou2*E0;
    } 
  }

  return 0;
}