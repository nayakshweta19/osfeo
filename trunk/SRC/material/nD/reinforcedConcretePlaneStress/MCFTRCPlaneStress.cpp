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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/MCFTRCPlaneStress.cpp,v $
// File: MCFTRCPlaneStress.cpp
//
// Written: Lining
// Created: 2011.8
//
// Description: This file contains the class definition for
// MCFTRCPlaneStress
// For Detailed explanation of the model, please refer to the journal paper
// entitled "The modified compression-filed theory for reinforced concrete element
// subjected to shear, ACI Journal. s83-22. 1986. pp. 219-231"

#include "MCFTRCPlaneStress.h"
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

static int numMCFTRCPlaneStressMaterials = 0;

OPS_Export void *
OPS_NewMCFTRCPlaneStressMaterial()
{
  if (numMCFTRCPlaneStressMaterials == 0) {
    numMCFTRCPlaneStressMaterials++;
    //OPS_Error("MCFTRCPlaneStress material - Written by Lining - Copyright@2011\n", 1);
  }

  // Pointer to a uniaxial material that will be returned
  NDMaterial *theMaterial = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 16) {
    opserr << "Invalid Args want: NDMaterial MCFTRCPlaneStress matTag? rho? UniaxiaMatTag1? UniaxiaMatTag2? UniaxiaMatTag3? UniaxiaMatTag4? \
               angle1? angle2? roux? rouy? dbx? dby? fpc? fyx? fyy? Ec0? epsc0? aggsize? xd? yd?\n";
    return 0;	
  }

  int tag;
  double rho;
  int    iData[4];
  double dData[14];
  int numData = 0;

  numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid nDMaterial MCFTRCPlaneStress tag" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDouble(&numData, &rho) != 0) {
    opserr << "WARNING invalid parameter rho MCFTRCPlaneStress tag:" << tag << endln;
    return 0;	
  }

  numData = 4;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid nDMaterial MCFTRCPlaneStress tag: " << tag << endln;
    return 0;
  }

  numData = 14;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data MCFTRCPlaneStress tag: " << tag << endln;
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial1 = OPS_GetUniaxialMaterial(iData[0]);
    
  if (theUniaxialMaterial1 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[0];
    opserr << "\nMCFTRCPlaneStress tag: " << tag << endln;
    return 0;
  }
  
  UniaxialMaterial *theUniaxialMaterial2 = OPS_GetUniaxialMaterial(iData[1]);

  if (theUniaxialMaterial2 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[1];
    opserr << "\nMCFTRCPlaneStress tag: " << tag << endln;
    return 0;
  }
  
  UniaxialMaterial *theUniaxialMaterial3 = OPS_GetUniaxialMaterial(iData[2]);
  if (theUniaxialMaterial3 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[2];
    opserr << "\nMCFTRCPlaneStress tag: " << tag << endln;
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial4 = OPS_GetUniaxialMaterial(iData[3]);  
  if (theUniaxialMaterial4 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[3];
    opserr << "\nMCFTRCPlaneStress tag: " << tag << endln;
    return 0;
  }

  //now create the MCFTRCPlaneStress
  theMaterial = new MCFTRCPlaneStress (tag, rho,
						     theUniaxialMaterial1, theUniaxialMaterial2, 
						     theUniaxialMaterial3, theUniaxialMaterial4,
						     dData[0], dData[1], dData[2], dData[3], dData[4],
							 dData[5], dData[6], dData[7], dData[8], dData[9],
							 dData[10],dData[11],dData[12],dData[13]);
       
  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory creating material\n";
    opserr << "MCFTRCPlaneStress: " << tag << endln;
    return 0;
  }

  return theMaterial;
}
 
MCFTRCPlaneStress ::MCFTRCPlaneStress (int tag, 
								   double RHO,
								   UniaxialMaterial *s1,
								   UniaxialMaterial *s2,
								   UniaxialMaterial *c1,
								   UniaxialMaterial *c2,
								   double   ANGLE1, double   ANGLE2,
								   double   ROU1,   double   ROU2,
								   double   DB1,    double   DB2,
								   double   FPC,
								   double   FYX,    double   FYY,
								   double   E,
								   double   EPSC0,
								   double   AGGR,
								   double   XD,     double   YD):
  NDMaterial(tag, ND_TAG_MCFTRCPlaneStress), 
  rho(RHO), angle1(ANGLE1), angle2(ANGLE2), roux(ROU1), rouy(ROU2), db1(DB1), db2(DB2),
  fpc(FPC), fyx(FYX), fyy(FYY), Es(E), epsc0(EPSC0), aggr(AGGR), xd(XD), yd(YD),
  lastStress(3), Tstress(3), stress_vec(3), stress0_vec(3), tangent_matrix(3,3), strain_vec(3), 
  epsC_vec(3), epsC12p_prevec(2), epsC12p_nowvec(2), epsCp_vec(3), 
  epsC12cm_vec(2), epsCcm_vec(3), epsC12tm_vec(2), epsCtm_vec(3),
  CepsC_vec(3), CepsCp_vec(3), CepsC12cm_vec(2), CepsCcm_vec(3), CepsC12tm_vec(2), CepsCtm_vec(3)
{
    CepsC_vec.Zero();
	CepsCp_vec.Zero();

	CepsCcm_vec.Zero();
	CepsCtm_vec.Zero();
    CepsC12cm_vec.Zero();
	CepsC12tm_vec.Zero();

    lastStress.Zero();  // add at 7.28    
    
    if ( fpc > 0.0 ) { fpc = -fpc; } // set fpc < 0
    if ( epsc0 > 0.0 ) { epsc0 = -epsc0; } // set fpc < 0
    theMaterial = 0;
    
    // Allocate pointers to theSteel1
    theMaterial = new UniaxialMaterial *[4];
    
    if ( theMaterial == 0 ) {
      opserr << " MCFTRCPlaneStress::MCFTRCPlaneStress - failed allocate material array\n";
      exit(-1);
    }
    
    // Get the copy for the Steel1
    theMaterial[0] = s1->getCopy();
    // Check allocation    
    if ( theMaterial[0] == 0 ) {
      opserr << " MCFTRCPlaneStress::MCFTRCPlaneStress - failed to get a copy for steel1\n";
      exit(-1);
    }
    
    // Get the copy for the Steel2
    theMaterial[1] = s2->getCopy();	
    // Check allocation    
    if ( theMaterial[1] == 0 ) {
      opserr << " MCFTRCPlaneStress::MCFTRCPlaneStress - failed to get a copy for steel2\n";
      exit(-1);
    }
    
    // Get the copy for the Concrete1
    theMaterial[2] = c1->getCopy();	
    // Check allocation    
    if ( theMaterial[2] == 0 ) {
      opserr << " MCFTRCPlaneStress::MCFTRCPlaneStress - failed to get a copy for concrete1\n";
      exit(-1);
    }
    
    // Get the copy for the Concrete2
    theMaterial[3] = c2->getCopy();	
    // Check allocation    
    if ( theMaterial[3] == 0 ) {
      opserr << " MCFTRCPlaneStress::MCFTRCPlaneStress - failed to get a copy for concrete2\n";
      exit(-1);
    }
    
    /* FMK */
    theResponses = new Response *[8];  
    
    if ( theResponses == 0) {
      opserr << " MCFTRCPlaneStress::MCFTRCPlaneStress - failed allocate responses array\n";
      exit(-1);
    }
    
    OPS_Stream *theDummyStream = new DummyStream();
    
    const char **argv = new const char *[1];

    argv[0] = "setVar";
    theResponses[0] = theMaterial[0]->setResponse(argv, 1, *theDummyStream);  //s1
    theResponses[1] = theMaterial[1]->setResponse(argv, 1, *theDummyStream);  //s2
	argv[0] = "getVar";
	theResponses[2] = theMaterial[0]->setResponse(argv, 1, *theDummyStream);  //s1
	theResponses[3] = theMaterial[1]->setResponse(argv, 1, *theDummyStream);  //s2
    argv[0] = "setVar";
    theResponses[4] = theMaterial[2]->setResponse(argv, 1, *theDummyStream);  //c1
    theResponses[5] = theMaterial[3]->setResponse(argv, 1, *theDummyStream);  //c2
    argv[0] = "getVar";
    theResponses[6] = theMaterial[2]->setResponse(argv, 1, *theDummyStream);  //c1
    theResponses[7] = theMaterial[3]->setResponse(argv, 1, *theDummyStream);  //c2

    if ((theResponses[0] == 0) || (theResponses[1] == 0) ||
	(theResponses[2] == 0) || (theResponses[3] == 0) ||
	(theResponses[4] == 0) || (theResponses[5] == 0) ||
	(theResponses[6] == 0) || (theResponses[7] == 0) ) {
      opserr << " MCFTRCPlaneStress::MCFTRCPlaneStress - failed to set appropriate materials tag: " << tag << "\n";
      exit(-1);
    }
    
    delete theDummyStream;
	strain_vec.Zero();
    this->revertToStart();
	determineTrialStress(strain_vec);
    
}

MCFTRCPlaneStress::MCFTRCPlaneStress()
  :NDMaterial(0, ND_TAG_MCFTRCPlaneStress),
  lastStress(3), Tstress(3), stress_vec(3), stress0_vec(3), tangent_matrix(3,3), strain_vec(3), 
  epsC12p_prevec(2), epsC12p_nowvec(2), epsC_vec(3), epsCp_vec(3), 
  epsC12cm_vec(2), epsCcm_vec(3), epsC12tm_vec(2), epsCtm_vec(3)
{
  CepsC_vec.Zero();
  CepsCp_vec.Zero();
  
  CepsCcm_vec.Zero();
  CepsCtm_vec.Zero();
  CepsC12cm_vec.Zero();
  CepsC12tm_vec.Zero();

  theMaterial = 0;
  theResponses = 0;

  strain_vec.Zero();

  this->revertToStart();
  determineTrialStress(strain_vec);

}

MCFTRCPlaneStress::~MCFTRCPlaneStress()
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
    for (int j=0; j<8; j++) {
	  if (theResponses[j] != 0)
	    delete theResponses[j];
    }
    delete [] theResponses;
  }
}

double 
MCFTRCPlaneStress::getRho(void)
{
	return rho;
}

// really used one
int 
MCFTRCPlaneStress::setTrialStrain(const Vector &v)
{

  // Set values for strain_vec
  strain_vec = v;
  
  // Set initial values for Tstress
  //Tstress.Zero();
  
  Tstress = determineTrialStress(strain_vec);
  
  //determineTangent();

  return 0;
}

int 
MCFTRCPlaneStress::setTrialStrain(const Vector &v, const Vector &r)
{
  opserr << "error: MCFTRCPlaneStress::setTrialStrain(&v, &r) -- not really responsibility" << endln;
  return 0;
}

int 
MCFTRCPlaneStress::setTrialStrainIncr(const Vector &v)
{
  opserr << "error:MCFTRCPlaneStress::setTrialStrainIncr(&v) -- not really responsibility" << endln;
  return 0;
}

int
MCFTRCPlaneStress::setTrialStrainIncr(const Vector &v, const Vector &r)
{
  opserr << "error:MCFTRCPlaneStress::setTrialStrainIncr(&v, &r) -- not really responsibility" << endln;
  return 0;
}

const Matrix& 
MCFTRCPlaneStress::getTangent(void)
{
  //determineTrialStress(strain_vec);
  //determineTangent();

  return tangent_matrix;
}

const Vector& 
MCFTRCPlaneStress::getStress(void)
{
  return stress_vec;
}

const Vector& 
MCFTRCPlaneStress::getStrain()
{
  return strain_vec;
}

const Vector& 
MCFTRCPlaneStress::getCommittedStress(void)
{
  return stress_vec;
}

const Vector& 
MCFTRCPlaneStress::getCommittedStrain(void)
{
  return strain_vec;
}

int
MCFTRCPlaneStress::commitState(void)
{
  for (int i=0; i<4; i++) {
    theMaterial[i]->commitState();
  }

  CepsC_vec  = epsC_vec; 
  CepsCp_vec = epsCp_vec;   

  CepsCcm_vec  = epsCcm_vec;  
  CepsCtm_vec  = epsCtm_vec;  
  CepsC12cm_vec = epsC12cm_vec;
  CepsC12tm_vec = epsC12tm_vec;

  lastStress = stress_vec;
  
  return 0;
}

int
MCFTRCPlaneStress::revertToLastCommit(void)
{
  for (int i=0; i<4; i++) {
    theMaterial[i]->revertToLastCommit();
  }

  epsC_vec  = CepsC_vec; 
  epsCp_vec = CepsCp_vec;   

  epsCcm_vec = CepsCcm_vec;  
  epsCtm_vec = CepsCtm_vec;  
  epsC12cm_vec = CepsC12cm_vec;
  epsC12tm_vec = CepsC12tm_vec;

  return 0;
}

int
MCFTRCPlaneStress::revertToStart(void)
{
  for (int i=0; i<4; i++) {
    theMaterial[i]->revertToStart();
  }
  
  Tstress.Zero();
  strain_vec.Zero();
  epsC_vec.Zero();
  epsC12p_prevec.Zero();
  epsC12p_nowvec.Zero();
  epsC12cm_vec.Zero();
  epsC12tm_vec.Zero();
  epsCcm_vec.Zero();
  epsCtm_vec.Zero();
  epsCp_vec.Zero();
  stress_vec.Zero();

  CepsC_vec.Zero();
  CepsCp_vec.Zero();

  CepsCcm_vec.Zero();
  CepsCtm_vec.Zero();
  CepsC12cm_vec.Zero();
  CepsC12tm_vec.Zero();
  
  return 0;
}

NDMaterial* 
MCFTRCPlaneStress::getCopy(void)
{
  MCFTRCPlaneStress* theCopy =
    new MCFTRCPlaneStress( this->getTag(), 
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
					 Es, 
					 epsc0,
					 aggr,
					 xd, yd);
  theCopy->strain_vec = strain_vec;
  return theCopy;
}

NDMaterial* 
MCFTRCPlaneStress::getCopy(const char *type)
{
	MCFTRCPlaneStress* theModel =
		new MCFTRCPlaneStress( this->getTag(), 
           rho, theMaterial[0], theMaterial[1], theMaterial[2], theMaterial[3],
		   angle1, angle2, roux, rouy, db1, db2, fpc, fyx, fyy, Es, epsc0, aggr, xd, yd );
	theModel->strain_vec = strain_vec;
	theModel->stress_vec = stress_vec;
	return theModel;
}

Response*
MCFTRCPlaneStress::setResponse (const char **argv, int argc, OPS_Stream &output)
{

#ifdef DEBUG
	opserr << "MCFTRCPlaneStress::setResponse(...)" << endln;
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
MCFTRCPlaneStress::getResponse (int responseID, Information &matInfo)
{
#ifdef DEBUG
	opserr << "MCFTRCPlaneStress::getResponse(...)" << endln;
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
MCFTRCPlaneStress::Print(OPS_Stream &s, int flag )
{
	s << "\n\tMCFTRCPlaneStress, material id: " << this->getTag() << endln;

	s << "\tRho: " << rho << endln;
	s << "\tangle1: " << angle1 << "\tangle2: " << angle2 << endln;
	s << "\trou1: " << roux << "\trou2: " << rouy << endln;
	s << "\tfpc: " << fpc << "\tepsc0: " << "\taggr: " << aggr <<endln;
	s << "\tfyx: " << fyx << "\tfyy: " << fyy << endln;
	s << "\tEs: " << Es << endln;

	s << "\tEnd of MCFTPSMatreial print()." << endln << endln;

	s << "\t call the material print() function : "<< endln;
	
	s << "\t the steel 1 information is : " << endln;
	theMaterial[0]->Print(s,flag);
	s << "\t the steel 2 information is : " << endln;
	theMaterial[1]->Print(s,flag);
	s << "\t the concrete 1 information is : " << endln;
	theMaterial[2]->Print(s,flag);
	s << "\t the concrete 2 information is : " << endln;
	theMaterial[3]->Print(s,flag);

	s << "\tStrain and stress of the uniaxial materials:"<<endln;
	for (int i=0; i<4; i++)
	{
	  s << "\tUniaxial Material " << i+1 << " :" << theMaterial[i]->getStrain() << "   " << theMaterial[i]->getStress() << endln;
	  s << "\t\t 1DMat Tangent : " << theMaterial[i]->getTangent() << endln;
	}
}

int
MCFTRCPlaneStress::sendSelf(int commitTag, Channel &theChannel)
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
  data(11) = Es;
  data(12) = epsc0;
  data(13) = aggr;
  data(14) = xd;
  data(15) = yd;

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING MCFTRCPlaneStress::sendSelf() - " << this->getTag() << " failed to send Vector\n";
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
    opserr << "WARNING MCFTRCPlaneStress::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }
  
  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 4; i++) {
    res += theMaterial[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "MCFTRCPlaneStress::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }	
  
  return res;
}

int
MCFTRCPlaneStress::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();
  
  // Quad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(16);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING MCFTRCPlaneStress::recvSelf() - failed to receive Vector\n";
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
  Es     = data(11);
  epsc0  = data(12);
  aggr   = data(13);
  xd     = data(14);
  yd     = data(15);

  static ID idData(8);
  
  // now receives the tags of its materials
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING MCFTRCPlaneStress::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  if (theMaterial == 0) {
    // Allocate new materials
    theMaterial = new UniaxialMaterial *[4];
    if (theMaterial == 0) {
      opserr << "MCFTRCPlaneStress::recvSelf() - Could not allocate UniaxialMaterial* array\n";
      return -1;
    }
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Allocate new material with the sent class tag
      theMaterial[i] = theBroker.getNewUniaxialMaterial(matClassTag);
      if (theMaterial[i] == 0) {
	     opserr << "MCFTRCPlaneStress::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
	     return -1;
      }
      // Now receive materials into the newly allocated space
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
         opserr << "MCFTRCPlaneStress::recvSelf() - material " << i << "failed to recv itself\n";
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
            opserr << "MCFTRCPlaneStress::recvSelf() - material " << i << "failed to create\n";
	        return -1;
		 }
	  }
      // Receive the material
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "MCFTRCPlaneStress::recvSelf() - material " << i << "failed to recv itself\n";
	    return res;
	  }
    }
  }

  return res;
}

Vector
MCFTRCPlaneStress::determineTrialStress(Vector strain)
{ 
  // Get Principal strain direction first
  Vector Tstrain(3);     //epslon_x, epslon_y, gamma_xy

  // Get strain values from strain of element
  Tstrain(0) = strain(0);
  Tstrain(1) = strain(1);
  Tstrain(2) = strain(2); // gamma

  double cita, citaS, citaE, temp_cita, citan1, citan2, citaLag, dCitaE, dCitaS, citaIC; // principal strain direction
  double epsC1, epsC2, epsSx, epsSy, epsts;
  double eC1p=0.0, eC2p=0.0, eC1m=0.0, eC2m=0.0, eT1m=0.0, eT2m=0.0;
  double fC1, fC1a, fC1b, fC2, fSx, fSy; //, fC1c

  double fcr   = 0.31*sqrt(-fpc); //0.65 * pow(-fpc, 0.33);      ----ftp
  double Ec    = 2.0 * fpc/epsc0;
  double epscr = fcr/Ec;
  double vci, vcimax, s_cita, w; //delta_S, gamma_S;
  double Cs, Cd, betaD;
  double Ecx, Ecy;

  epsC_vec  = Tstrain; //CepsC_vec; // initial epsC vector assumptions equal to previous committed value
  epsCp_vec = CepsCp_vec;

  // Get citaE based on Tstrain, eq. i-9...
  if ( fabs(Tstrain(0)-Tstrain(1)) < 1e-12 ) {
    citaE = 0.25*PI;
  } else {  // Tstrain(0) != Tstrain(1) 
    temp_cita = 0.5 * atan(fabs(1.0e6*Tstrain(2)/(1.0e6*Tstrain(0)-1.0e6*Tstrain(1)))); 
    if ( fabs(Tstrain(2)) < 1e-9 ) {
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

  citaIC = citaE;  // initial cracked angle

  // Get average strains in reinforcement, eq. i-19
  epsSx = 0.5 * (Tstrain(0)+Tstrain(1))+0.5*(Tstrain(0)-Tstrain(1))*cos(2.0*angle1)
	    + 0.5 * Tstrain(2)*sin(2.0*angle1);

  epsSy = 0.5 * (Tstrain(0)+Tstrain(1))+0.5*(Tstrain(0)-Tstrain(1))*cos(2.0*angle2)
	    + 0.5 * Tstrain(2)*sin(2.0*angle2);

  // Calculate Steel stress incorporate steel offset strains

  //theInfoS01.setVector(theData);
  //theInfoS02.setVector(theData);
  //theResponses[0]->getResponse();  // setVar
  //theResponses[1]->getResponse();  // setVar

  if (theMaterial[0]->setTrialStrain(epsSx) != 0) {
    opserr << "MCFTRCPlaneStress::determineTrialStress(), S1 setTrialStrain() Error" << endln;  // s1
	return -1;
  }
  if (theMaterial[1]->setTrialStrain(epsSy) != 0) {
    opserr << "MCFTRCPlaneStress::determineTrialStress(), S2 setTrialStrain() Error" << endln;  // s2
	return -1;
  }

  fSx = theMaterial[0]->getStress();
  fSy = theMaterial[1]->getStress();
  
  static Vector theData(7); // data for materials passed in

  static Matrix Tc(3,3);     // T(cita)
  static Matrix Ts1(3,3);    // T(angle1)
  static Matrix Ts2(3,3);    // T(angle2)
  static Matrix Dc(3,3);     //
  static Matrix Ds1(3,3);    //
  static Matrix Ds2(3,3);    //
  static Matrix Dcp(3,3);    //
  static Matrix Ds1p(3,3);	 //
  static Matrix Ds2p(3,3);	 //

  // EbarC1 EbarC2 GbarC, eq.i-5
  Dcp.Zero();  // local stiffness matrix
  Ds1p.Zero();
  Ds2p.Zero();

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

  Information &theInfoS01 = theResponses[0]->getInformation(); // S1 setVar
  Information &theInfoS02 = theResponses[1]->getInformation(); // S2 setVar
  Information &theInfoS03 = theResponses[2]->getInformation(); // S1 getVar
  Information &theInfoS04 = theResponses[3]->getInformation(); // S2 getVar

  Information &theInfoC01 = theResponses[4]->getInformation(); // C1 setVar
  Information &theInfoC02 = theResponses[5]->getInformation(); // C2 setVar
  Information &theInfoC03 = theResponses[6]->getInformation(); // C1 getVar
  Information &theInfoC04 = theResponses[7]->getInformation(); // C2 getVar

  //for small strain return quickly
  if (Tstrain.pNorm(-1) <= DBL_EPSILON) {
	double v = 0.2;
	double temp = Ec/(1.0-v*v);

	Dcp(0,0) = temp;
	Dcp(0,1) = v*temp;
	Dcp(1,0) = v*temp;
	Dcp(2,2) = temp*(1.0-v)/2.0;

	Ds1p(0,0)= Es;
	Ds2p(0,0)= Es;

	Dc.addMatrixTripleProduct( 0.0, Tc,  Dcp,  1.0);
	Ds1.addMatrixTripleProduct(0.0, Ts1, Ds1p, 1.0);
	Ds2.addMatrixTripleProduct(0.0, Ts2, Ds2p, 1.0);

	tangent_matrix = Dc + Ds1 + Ds2;

	Tstress.addMatrixVector(0.0, tangent_matrix, strain_vec, 1.0);

	return Tstress;
  }

  // vCi need to be satisfied several relationships in the concrete crack region
  static Vector errorVec(3); // epsC_vec converge test need a vector to maintain Error
  
  double errNorm;
  double temp, tolNorm = 1.0e-9;
  bool vCiconverged = false;
  int iteration_counter = 0;

  while (vCiconverged == false) {

    // Get epsC1,epsC2 and citaS based on Tstrain, eq.i-10
    epsC1 = 0.5*(epsC_vec(0)+epsC_vec(1)) 
		  + 0.5*sqrt(pow(epsC_vec(0)-epsC_vec(1), 2.0)+pow(epsC_vec(2), 2.0));

    epsC2 = 0.5*(epsC_vec(0)+epsC_vec(1)) 
		  - 0.5*sqrt(pow(epsC_vec(0)-epsC_vec(1), 2.0)+pow(epsC_vec(2), 2.0));
  
    // Get citaS based on epsC_vec, eq.i-11
    if ( fabs(epsC_vec(0)-epsC_vec(1)) < DBL_EPSILON) {
      citaS = 0.25*PI;	
    } else {  // epsC_vec(0) != epsC_vec(1) 
      temp_cita = 0.5 * atan(fabs(1.0e6*epsC_vec(2)/(1.0e6*epsC_vec(0)-1.0e6*epsC_vec(1)))); 
      if ( fabs(epsC_vec(2)) < 1e-9 ) {
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
	//if (citan1 < 0.0) citan1 += 2.0*PI;
    citan2 = citaS - angle2; //angle2-(PI/2.0-citaS);
    //if (citan2 < 0.0) citan2 += 2.0*PI;

	// calculate betaD
    Cs = 0.55;
    temp = epsC1/epsC2+0.28;
    if (temp <= 0.0) Cd = 0.35*pow(-temp, 0.8);
    else             Cd = 0.0;
  
    betaD = 1.0/(1.0+Cs*Cd);
    if (betaD > 1.0) betaD = 1.0;
  
    // Calculate average concrete stress
    // C2 --compression for fc2, epsC2 incorporate steel offset strains

	theData(0) = eC2m;  // 
	theData(1) = eT2m;  // 
	theData(5) = eC2p;
	//theData(2) = epsC2;  // 

	int ret = 0;
	ret += theInfoC02.setVector(theData);  // Information for C compression, material [3], response[5]
	ret += theResponses[5]->getResponse(); // C2 setVar
	ret += theMaterial[3]->setTrialStrain(epsC2); //need further revision
    if (ret != 0) {
	  opserr << "MCFTRCPlaneStress::determineTrialStress(), C2 setTrialStrain() Error" << endln;  // C2
	  return ret;
	}
	fC2 = betaD * theMaterial[3]->getStress(); // w==5 ??
  
	theResponses[7]->getResponse(); // C2 getVar to theInfoC04
    epsC12p_nowvec(1) = (*theInfoC04.theVector)(5); // ecp2 
    
    // C1 -- tensile for fc1, epsC1
    //  case 1: eq.i-33,34
    epsts = 2.0 * 7.5e-3 / fcr / 750.; //eq. i-33
    temp  = (epsC1 - epscr) / (epsts - epscr);

    if      (temp <= 0.0) fC1a = 0.0;
    else if (temp >  1.0) fC1a = fcr;
    else                  fC1a = fcr * (1 - temp);  //eq. i-34
  
    //  case 2: eq. i-35, 36
    temp = 4.0*(roux/db1*abs(cos(citan1))+rouy/db2*abs(cos(citan2)));
    fC1b = fcr/(1+sqrt(2.2/temp*epsC1));
  
    //temp = max(fC1a,fC1b); // eq.i-38 /fC1
  
    // eq.i=5   // for positive value
    //fC1c = max(roux* (fyx - fSx) * pow(cos(citan1), 2.0), 0.0)
  	//     + max(rouy* (fyy - fSy) * pow(cos(citan2), 2.0), 0.0);
    
    //fC1 = min(temp, fC1c);
    
	double fC1temp = max(fC1a, fC1b);

	theData(0) = eC1m;  // 
	theData(1) = eT1m;  // 
	theData(5) = eC1p;
	//theData(2) = epsC1;  // 

	ret = 0;
	ret += theInfoC01.setVector(theData);  // Information for C tension, material [2], response[4]
	ret += theResponses[4]->getResponse(); // C1 setVar 
    ret += theMaterial[2]->setTrialStrain(epsC1);
	if (ret != 0) {
	  opserr << "MCFTRCPlaneStress::determineTrialStress(), C1 setTrialStrain() Error" << endln;  // C1
	  return ret;
	}
	fC1 = theMaterial[2]->getStress();

	theResponses[6]->getResponse();  // C1 getVar to theInfoC03
    epsC12p_nowvec(0) = (*theInfoC03.theVector)(5); // ecp1 

    if (epsC12p_nowvec(0) > 0) epsC12p_nowvec(0) = 0.0;

	// MCFT cita corresponding strain field, DSFM cita corresponding stress field
    if ( fabs(epsC_vec(0)-epsC_vec(1)) < 1e-12 ) {
      cita = 0.25*PI;	
    } else {  // epsC_vec(0) != epsC_vec(1) 
      temp_cita = 0.5 * atan(fabs(1.0e6*epsC_vec(2)/(1.0e6*epsC_vec(0)-1.0e6*epsC_vec(1)))); 
      if ( fabs(epsC_vec(2)) < 1e-9 ) {
  	    cita = 0;
  	  } else if ( (epsC_vec(0) > epsC_vec(1)) && ( epsC_vec(2) > 0) ) {
  	    cita = temp_cita;
  	  } else if ( (epsC_vec(0) > epsC_vec(1)) && ( epsC_vec(2) < 0) ) {
  	    cita = PI - temp_cita;
  	  } else if ( (epsC_vec(0) < epsC_vec(1)) && ( epsC_vec(2) > 0) ) {
  	    cita = 0.5*PI - temp_cita;
  	  } else if ( (epsC_vec(0) < epsC_vec(1)) && ( epsC_vec(2) < 0) ) {
  	    cita = 0.5*PI + temp_cita;
  	  } else {
  	    opserr << "DSFMRCPlaneStress::determineTrialStress: Failure to calculate citaS\n";
  	    opserr << " epsC_vec(0) = " << epsC_vec(0) << endln;
  	    opserr << " epsC_vec(1) = " << epsC_vec(1) << endln;
  	    opserr << " epsC_vec(2) = " << epsC_vec(2) << endln;
  	  }
    }

	// Get cita based on epsC_vec+epsCp_vec=Tstrain
    // epsCp_vec = Tstrain - epsC_vec;
	epsC12p_prevec(0) = 0.5*(epsCp_vec(0)+epsCp_vec(1)) 
		+ 0.5*sqrt(pow(epsCp_vec(0)-epsCp_vec(1), 2.0)+pow(epsCp_vec(2), 2.0));

	epsC12p_prevec(1) = 0.5*(epsCp_vec(0)+epsCp_vec(1)) 
		- 0.5*sqrt(pow(epsCp_vec(0)-epsCp_vec(1), 2.0)+pow(epsCp_vec(2), 2.0));

    // deltaEspCp(2),  P. Ceresa 2009, Eq. 4
	double deltaEspC1p = epsC12p_nowvec(0) - epsC12p_prevec(0);
	double deltaEspC2p = epsC12p_nowvec(1) - epsC12p_prevec(1);
  
    epsCp_vec(0) += (0.5 * deltaEspC1p * (1+cos(2.0*cita)) + 0.5 * deltaEspC2p * (1-cos(2.0*cita)));
    epsCp_vec(1) += (0.5 * deltaEspC1p * (1-cos(2.0*cita)) + 0.5 * deltaEspC2p * (1+cos(2.0*cita)));
    epsCp_vec(2) += (deltaEspC1p - deltaEspC2p) * sin(2.0*cita);
  
    eC1p = 0.5 * (epsCp_vec(0)+epsCp_vec(1)) 
		 + 0.5 * ((epsCp_vec(0)-epsCp_vec(1))*cos(2.0*cita)+epsCp_vec(2)*sin(2.0*cita));
    eC2p = 0.5 * (epsCp_vec(0)+epsCp_vec(1)) 
		 - 0.5 * ((epsCp_vec(0)-epsCp_vec(1))*cos(2.0*cita)+epsCp_vec(2)*sin(2.0*cita));
    
	// Vecchio: Towards Cyclic Load Modeling of Reinforced Concrete
	// C max
	double deltaEspCm1 = (epsC1 > epsC12cm_vec(0) ? 0 : epsC1-epsC12cm_vec(0));
	double deltaEspCm2 = (epsC2 > epsC12cm_vec(1) ? 0 : epsC2-epsC12cm_vec(1));
  
	epsCcm_vec(0) += (0.5 * deltaEspCm1 * (1+cos(2.0*cita)) + 0.5 * deltaEspCm2 * (1-cos(2.0*cita)));
	epsCcm_vec(1) += (0.5 * deltaEspCm1 * (1-cos(2.0*cita)) + 0.5 * deltaEspCm2 * (1+cos(2.0*cita)));
	epsCcm_vec(2) += (deltaEspCm1 - deltaEspCm2) * sin(2.0*cita);
  
    eC1m = 0.5*(epsCcm_vec(0)+epsCcm_vec(1))+0.5*((epsCcm_vec(0)-epsCcm_vec(1))*cos(2.0*cita)+epsCcm_vec(2)*sin(2.0*cita));
    eC2m = 0.5*(epsCcm_vec(0)+epsCcm_vec(1))-0.5*((epsCcm_vec(0)-epsCcm_vec(1))*cos(2.0*cita)+epsCcm_vec(2)*sin(2.0*cita));

	// T max
	double deltaEspTm1 = (epsC1 < epsC12tm_vec(0) ? 0 : epsC1-epsC12tm_vec(0));
	double deltaEspTm2 = (epsC2 < epsC12tm_vec(1) ? 0 : epsC2-epsC12tm_vec(1));
  
	epsCtm_vec(0) += (0.5*deltaEspTm1*(1+cos(2.0*cita)) + 0.5*deltaEspTm2*(1-cos(2.0*cita)));
	epsCtm_vec(1) += (0.5*deltaEspTm1*(1-cos(2.0*cita)) + 0.5*deltaEspTm2*(1+cos(2.0*cita)));
	epsCtm_vec(2) += (deltaEspTm1 - deltaEspTm2) * sin(2.0*cita);
  
	eT1m = 0.5*(epsCtm_vec(0)+epsCtm_vec(1))+0.5*((epsCtm_vec(0)-epsCtm_vec(1))*cos(2.0*cita)+epsCtm_vec(2)*sin(2.0*cita));
	eT2m = 0.5*(epsCtm_vec(0)+epsCtm_vec(1))-0.5*((epsCtm_vec(0)-epsCtm_vec(1))*cos(2.0*cita)+epsCtm_vec(2)*sin(2.0*cita));

    // Calculate local stress at cracks
    int status            = 0; // status to check if iteration satisfied eq.i-7

    double tolerance      = 1.0e-9; // tolerance for iteration
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
  	  
  	  if( theMaterial[0]->setTrialStrain(epsScrx) != 0) {
		opserr << "DSFMRCPlaneStress::determineTrialStress(): fail to determine epsScrx" << endln;
		return -1;
	  }
  	  fScrx = theMaterial[0]->getStress();
  	  
  	  if( theMaterial[1]->setTrialStrain(epsScry) != 0) {
		opserr << "DSFMRCPlaneStress::determineTrialStress(): fail to determine epsScry" << endln;
		return -1;
	  }
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
  	    temp = 0.501 * (theMaterial[0]->getSecant() + theMaterial[1]->getSecant());
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
    w      = epsC1 * s_cita;
  
    // optional 1
    vcimax = sqrt(-fpc)/(0.31+24*w/(aggr+16)); // shear parameter

	vci = vci * min(abs(vcimax), abs(vci)) /abs(vci);

    //need further revision for the initial crack direction, citaIC
    if (roux > 0.0 && rouy > 0.0) citaLag = 5./180.*PI;    // eq.i-41~43 
    else                          citaLag = 7.5/180.*PI;
  
	dCitaE = citaE - citaIC;
	dCitaS = (abs(dCitaE) > citaLag ? dCitaE-citaLag : dCitaE);
	citaS  = citaIC + dCitaS;
  
    if (abs(epsC1) < DBL_EPSILON) {
	  Ecx = theMaterial[2]->getSecant();
	//} else if (epsC1 < eC1p) {
    //  temp = epsC1 - (eC1p);
	//  epsC_vec(0) += 0.5 * temp * (1+cos(2.0*citaS));
    //  epsC_vec(1) += 0.5 * temp * (1-cos(2.0*citaS));
	//  epsC_vec(2) += temp *sin(2.0*citaS);
	//  continue;
	} else {
	  //Ecx = fC1/(epsC1-eC1p-eSlip1-eC01);
	  Ecx = fC1/epsC1;
	  if (Ecx > Es) Ecx = Es;
	}
  
    if (abs(epsC2) < DBL_EPSILON) {
	  Ecy = theMaterial[3]->getSecant();
	//} else if (epsC2 > eC2p) {
	//  temp = epsC2 - (eC2p);
	//  epsC_vec(0) += 0.5 * temp * (1-cos(2.0*citaS));
    //  epsC_vec(1) += 0.5 * temp * (1+cos(2.0*citaS));
	//  epsC_vec(2) -= temp *sin(2.0*citaS);
	//  continue;
	} else {
	  Ecy = fC2/epsC2;
	  if (Ecy > Es) Ecy = Es;
	}
  
    Dcp(0,0) = Ecx;
    Dcp(1,1) = Ecy;
    Dcp(2,2) = Ecx * Ecy / ( Ecx + Ecy );
  
    // EbarS1, EbarS2, eq.i-5
    
    if (abs(epsSx) < DBL_EPSILON) Ds1p(0,0) = roux * (theMaterial[0]->getSecant());
    else              Ds1p(0,0) = roux * fSx / epsSx;
    
    if (abs(epsSy) < DBL_EPSILON) Ds2p(0,0) = rouy * (theMaterial[1]->getSecant());
    else              Ds2p(0,0) = rouy * fSy / epsSy;
  
    Dc.addMatrixTripleProduct( 0.0, Tc,  Dcp,  1.0);
    Ds1.addMatrixTripleProduct(0.0, Ts1, Ds1p, 1.0);
    Ds2.addMatrixTripleProduct(0.0, Ts2, Ds2p, 1.0);

	tangent_matrix = Dc + Ds1 + Ds2;

    // Calculate pseudo-prestress vector

    stress0_vec.addMatrixVector(0.0, Dc, epsCp_vec, 1.0);

    static Matrix DInv(3,3);  // D invert
    tangent_matrix.Invert(DInv);

	Vector tempVec(3);
	tempVec.addMatrixVector(0.0, DInv, stress0_vec, 1.0);

	epsCp_vec -= tempVec;

	errorVec = Tstrain - epsCp_vec - epsC_vec;

    errNorm = errorVec.pNorm(-1);
	// determine the norm of matrixnorm = tempMat.Norm();
	// errNorm = norm - tempNorm;
	
  	if ( errNorm <= tolNorm) {
	  vCiconverged = true;
	  OPS_Stream *output = &opserr;
	  this->Print(*output,1);
	} else {
	  //tempNorm = tempMat.Norm();
	  //epsC_vec += 0.5 * errorVec; //Tstrain - epsCp_vec;
	  epsC_vec = Tstrain - epsCp_vec;
	  epsC12p_prevec = epsC12p_nowvec;
	  iteration_counter += 1 ;
	  if (iteration_counter > 20) {
		iteration_counter = 0;
		//epsC_vec = 0.83 *(Tstrain - epsCp_vec);
		epsC12p_prevec *= 0.5;
		epsC12p_nowvec *= 0.5;
	  }
	}
  }

  // Calculate the tangent matrix.
  tangent_matrix = Dc + Ds1 + Ds2;

  // Calculate total stress from steel and concrete
  Tstress.addMatrixVector(0.0, tangent_matrix, strain_vec, 1.0);

  // determine internal vectors

  return Tstress;
}

double
MCFTRCPlaneStress::kupferEnvelop(double Tstrain, double sig_p, double eps_p)
{
  double sig;
  if (Tstrain > eps_p) {
    double eta = Tstrain/eps_p;
    sig = sig_p * (2 * eta - eta * eta);
  }
  else if (Tstrain > 2.0 * epsc0) {
	double eta = (Tstrain-eps_p)/(2.0*epsc0-eps_p);
    sig = sig_p * (1.0 - eta * eta);
  }
  else {
    sig = 1.0e-9 * fpc;
  }
  return sig;

}