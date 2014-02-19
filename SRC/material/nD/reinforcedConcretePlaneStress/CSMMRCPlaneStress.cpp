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
// $Date: 2003/02/14 23:00:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/CSMMRCPlaneStress.cpp,v $
// File: CSMMRCPlaneStress.cpp
//
// Written: Lining
// Created: 2010.11
//
// Description: This file contains the class definition for
// CSMMRCPlaneStress
// For Detailed explanation of the model, please refer to the book
// entitled "Unified Theory of Concrete Structures,"
// by Thomas T.C. Hsu and Y.L. Mo, John Wiley & Sons, April 2010.

#include "CSMMRCPlaneStress.h"
#include <UniaxialMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <math.h>
#include <float.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <stdlib.h>

#include <Tensor.h>

Vector CSMMRCPlaneStress :: strain_vec(3);
Vector CSMMRCPlaneStress :: stress_vec(3);
Matrix CSMMRCPlaneStress :: tangent_matrix(3,3);

double CSMMRCPlaneStress :: epslonOne = 0.0;
double CSMMRCPlaneStress :: epslonTwo = 0.0;
double CSMMRCPlaneStress :: halfGammaOneTwo = 0.0;

double CSMMRCPlaneStress :: sigmaOneC = 0.0;
double CSMMRCPlaneStress :: sigmaTwoC = 0.0;

double CSMMRCPlaneStress::citaR = 0.0;      // principal strain direction
double CSMMRCPlaneStress::lastCitaR = 0.0;  // last converged principle strain direction
int    CSMMRCPlaneStress::steelStatus = 0;  // check if steel yield, 0 not yield, 1 yield
int    CSMMRCPlaneStress::dirStatus = 0; // check if principle direction has exceed 90 degree, 1 yes, 0 no
bool   CSMMRCPlaneStress::isSwapped = 0; // counter-clockwise = 0; clockwise = 1;
int    CSMMRCPlaneStress::lastDirStatus = 0; // 0, 1, 2, 3, 4, 5

#include <DummyStream.h>
#include <MaterialResponse.h>
#include <elementAPI.h>

#define OPS_Export 

static int numCSMMRCPlaneStressMaterials = 0;

OPS_Export void *
OPS_NewCSMMRCPlaneStressMaterial()
{
  if (numCSMMRCPlaneStressMaterials == 0) {
    numCSMMRCPlaneStressMaterials++;
    //OPS_Error("CSMMRCPlaneStress uniaxial material - Written by J.Zhong, Thomas T.C. Hsu and Y.L. Mo - Copyright@2009\n", 1);
  }

  // Pointer to a uniaxial material that will be returned
  NDMaterial *theMaterial = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 14) {
    opserr << "Invalid Args want: NDMaterial CSMMRCPlaneStress matTag? rho? UniaxiaMatTag1? UniaxiaMatTag2? UniaxiaMatTag3? UniaxiaMatTag4? \
              angle1? angle2? rou1? rou2? fpc? fy? E0? epsc0?\n";
    return 0;
  }

  int tag;
  double rho;
  int    iData[4];
  double dData[8];
  int numData = 0;

  numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid nDMaterial CSMMRCPlaneStress tag" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDouble(&numData, &rho) != 0) {
    opserr << "WARNING invalid parameter rho CSMMRCPlaneStress tag:" << tag << endln;
    return 0;
  }

  numData = 4;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid nDMaterial CSMMRCPlaneStress tag: " << tag << endln;
    return 0;
  }

  numData = 8;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data CSMMRCPlaneStress tag: " << tag << endln;
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial1 = OPS_GetUniaxialMaterial(iData[0]);

  if (theUniaxialMaterial1 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[0];
    opserr << "\nCSMMRCPlaneStress tag: " << tag << endln;
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial2 = OPS_GetUniaxialMaterial(iData[1]);

  if (theUniaxialMaterial2 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[1];
    opserr << "\nCSMMRCPlaneStress tag: " << tag << endln;
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial3 = OPS_GetUniaxialMaterial(iData[2]);
  if (theUniaxialMaterial3 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[2];
    opserr << "\nCSMMRCPlaneStress tag: " << tag << endln;
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial4 = OPS_GetUniaxialMaterial(iData[3]);
  if (theUniaxialMaterial4 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[3];
    opserr << "\nCSMMRCPlaneStress tag: " << tag << endln;
    return 0;
  }

  //now create the CSMMRCPlaneStress
  theMaterial = new CSMMRCPlaneStress(tag,
    rho,
    theUniaxialMaterial1,
    theUniaxialMaterial2,
    theUniaxialMaterial3,
    theUniaxialMaterial4,
    dData[0],
    dData[1],
    dData[2],
    dData[3],
    dData[4],
    dData[5],
    dData[6],
    dData[7]);

  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory creating material\n";
    opserr << "CSMMRCPlaneStress: " << tag << endln;
    return 0;
  }

  return theMaterial;
}

CSMMRCPlaneStress::CSMMRCPlaneStress(int      tag,
                                     double   RHO,
                                     UniaxialMaterial *s1,
                                     UniaxialMaterial *s2,
                                     UniaxialMaterial *c1,
                                     UniaxialMaterial *c2,
                                     double   ANGLE1,
                                     double   ANGLE2,
                                     double   ROU1,
                                     double   ROU2,
                                     double   FPC,
                                     double   FY,
                                     double   E,
                                     double   EPSC0) :
NDMaterial(tag, ND_TAG_CSMMRCPlaneStress),
  rho(RHO), angle1(ANGLE1), angle2(ANGLE2), rou1(ROU1), rou2(ROU2),
  fpc(FPC), fy(FY), E0(E), epsc0(EPSC0), lastStress(3), Tstress(3), Tstrain(3)
{
  if (fpc < 0.0) { fpc = -fpc; } // set fpc > 0

  citaStrain = 0.0;
  citaStress = 0.0;

  // Allocate pointers to uniaxial materials
  theMaterial = new UniaxialMaterial *[4];

  if (theMaterial == 0) {
    opserr << " CSMMRCPlaneStress::CSMMRCPlaneStress - failed allocate material array\n";
    exit(-1);
  }

  // Get the copy for theSteel1
  theMaterial[0] = s1->getCopy();
  // Check allocation    
  if (theMaterial[0] == 0) {
    opserr << " CSMMRCPlaneStress::CSMMRCPlaneStress - failed to get a copy for steel1\n";
    exit(-1);
  }

  // Get the copy for theSteel2
  theMaterial[1] = s2->getCopy();
  // Check allocation    
  if (theMaterial[1] == 0) {
    opserr << " CSMMRCPlaneStress::CSMMRCPlaneStress - failed to get a copy for steel2\n";
    exit(-1);
  }

  // Get the copy for theConcrete1
  theMaterial[2] = c1->getCopy();
  // Check allocation    
  if (theMaterial[2] == 0) {
    opserr << " CSMMRCPlaneStress::CSMMRCPlaneStress - failed to get a copy for concrete1\n";
    exit(-1);
  }

  // Get the copy for theConcrete2
  theMaterial[3] = c2->getCopy();
  // Check allocation    
  if (theMaterial[3] == 0) {
    opserr << " CSMMRCPlaneStress::CSMMRCPlaneStress - failed to get a copy for concrete2\n";
    exit(-1);
  }

  /* FMK */
  theResponses = new Response *[6];

  if (theResponses == 0) {
    opserr << " CSMMRCPlaneStress::CSMMRCPlaneStress - failed allocate responses array\n";
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
      opserr << " CSMMRCPlaneStress::CSMMRCPlaneStress - failed to set appropriate materials tag: " << tag << "\n";
      exit(-1);
  }

  delete theDummyStream;
  /* END FMK */
  determineTrialStress();
  this->revertToStart();
}

CSMMRCPlaneStress::CSMMRCPlaneStress()
  :NDMaterial(0, ND_TAG_CSMMRCPlaneStress)
{
  theMaterial = 0;
  theResponses = 0;
  determineTrialStress();
  this->revertToStart();
}


CSMMRCPlaneStress::~CSMMRCPlaneStress()
{
  // Delete the pointers
  if (theMaterial != 0) {
    for (int i = 0; i < 4; i++) {
      if (theMaterial[i] != 0)
        delete theMaterial[i];
    }
    delete[] theMaterial;
  }

  if (theResponses != 0) {
    for (int j = 0; j < 6; j++) {
      if (theResponses[j] != 0)
        delete theResponses[j];
    }
    delete[] theResponses;
  }
}

int
CSMMRCPlaneStress::setTrialStrain(const Vector &v)
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

  determineTrialStress();

  return 0;
}

int
CSMMRCPlaneStress::setTrialStrain(const Vector &v, const Vector &r)
{
  return 0;
}

int
CSMMRCPlaneStress::setTrialStrainIncr(const Vector &v)
{
  return 0;
}

int
CSMMRCPlaneStress::setTrialStrainIncr(const Vector &v, const Vector &r)
{
  return 0;
}

double
CSMMRCPlaneStress::getRho(void)
{
  return rho;
}

const Matrix&
CSMMRCPlaneStress::getTangent(void)
{
  return tangent_matrix;
}

const Vector&
CSMMRCPlaneStress::getStress(void)
{
  return stress_vec;
}

const Vector&
CSMMRCPlaneStress::getStrain()
{
  return strain_vec;
}

const Vector&
CSMMRCPlaneStress::getCommittedStress(void)
{
  return stress_vec;
}

const Vector&
CSMMRCPlaneStress::getCommittedStrain(void)
{
  return strain_vec;
}

int
CSMMRCPlaneStress::commitState(void)
{
  for (int i = 0; i < 4; i++) {
	theMaterial[i]->commitState();
  }

  COneReverseStatus = TOneReverseStatus;
  COneNowMaxComStrain = TOneNowMaxComStrain;
  COneLastMaxComStrain = TOneLastMaxComStrain;

  CTwoReverseStatus = TTwoReverseStatus;
  CTwoNowMaxComStrain = TTwoNowMaxComStrain;
  CTwoLastMaxComStrain = TTwoLastMaxComStrain;

  /*char buffer[200];
  sprintf(buffer, "eS1 = %8.6f, eS2 = %8.6f, eC1 = %8.6f, eC2 = %8.6f; ThetaE = %4.2f, ThetaS = %4.2f; |S(i)|-|S(i-1)|=%8.6e. ",
	  theMaterial[0]->getStrain(), theMaterial[1]->getStrain(), theMaterial[2]->getStrain(), theMaterial[3]->getStrain(),
	  citaStrain / PI * 180, citaStress / PI * 180, fabs(lastStress.Norm() - stress_vec.Norm()));
  opserr << buffer;
  sprintf(buffer, "Te1 = %8.6f, Te2 = %8.6f, Te3 = %8.6f. ", strain_vec(0), strain_vec(1), strain_vec(2));
  opserr << buffer ;//<< endln
  if ((epslonOne > 0.0) && (epslonTwo > 0.0)) {  // both tension
	opserr << "1t2t. e1 = " << epslonOne << ", e2 = " << epslonTwo << ". dir=";
  }
  else if ((epslonOne > 0.0) && (epslonTwo <= 0.0)) {  // one tension, two compression
	opserr << "1t2c. e1 = " << epslonOne << ", e2 = " << epslonTwo << ". dir=";
  }
  else if ((epslonOne <= 0.0) && (epslonTwo > 0.0)) {  // one compression, two tension
	opserr << "1c2t. e1 = " << epslonOne << ", e2 = " << epslonTwo << ". dir=";
  }
  else if ((epslonOne <= 0.0) && (epslonTwo <= 0.0)) {  //both compression
	opserr << "1c2c. e1 = " << epslonOne << ", e2 = " << epslonTwo << ". dir=";
  }
  opserr << dirStatus << ". swap=" << int(isSwapped) << endln;*/

  lastStress = stress_vec;
  lastCitaR = citaR;
  lastDirStatus = dirStatus;

  return 0;
}

int
CSMMRCPlaneStress::revertToLastCommit(void)
{
  for (int i = 0; i < 4; i++)
    theMaterial[i]->revertToLastCommit();

  TOneReverseStatus = COneReverseStatus;
  TOneNowMaxComStrain = COneNowMaxComStrain;
  TOneLastMaxComStrain = COneLastMaxComStrain;

  TTwoReverseStatus = CTwoReverseStatus;
  TTwoNowMaxComStrain = CTwoNowMaxComStrain;
  TTwoLastMaxComStrain = CTwoLastMaxComStrain;

  return 0;
}

int
CSMMRCPlaneStress::revertToStart(void)
{
  for (int i = 0; i < 4; i++)
    theMaterial[i]->revertToStart();

  Tstress.Zero();
  Tstrain.Zero();
  lastStress.Zero();

  strain_vec.Zero();
  stress_vec.Zero();

  steelStatus = 0;
  dirStatus = 0;
  G12 = 0;

  DDOne = 1.0;
  DDTwo = 1.0;

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

NDMaterial*
CSMMRCPlaneStress::getCopy(void)
{
  CSMMRCPlaneStress* theCopy =
    new CSMMRCPlaneStress(this->getTag(),
    rho,
    theMaterial[0],
    theMaterial[1],
    theMaterial[2],
    theMaterial[3],
    angle1,
    angle2,
    rou1,
    rou2,
    fpc,
    fy,
    E0,
    epsc0);
  theCopy->strain_vec = strain_vec;
  theCopy->stress_vec = stress_vec;
  theCopy->lastCitaR = lastCitaR;
  theCopy->lastDirStatus = lastDirStatus;

  return theCopy;
}

NDMaterial*
CSMMRCPlaneStress::getCopy(const char *type)
{
  CSMMRCPlaneStress* theModel =
    new CSMMRCPlaneStress(this->getTag(),
    rho,
    theMaterial[0],
    theMaterial[1],
    theMaterial[2],
    theMaterial[3],
    angle1,
    angle2,
    rou1,
    rou2,
    fpc,
    fy,
    E0,
    epsc0);
  theModel->strain_vec = strain_vec;
  theModel->stress_vec = stress_vec;
  theModel->lastCitaR = lastCitaR;
  theModel->lastDirStatus = lastDirStatus;

  return theModel;
}
//added by Ln
Response*
CSMMRCPlaneStress::setResponse(const char **argv, int argc, OPS_Stream &output)
{

#ifdef DEBUG
  opserr << "CSMMRCPlaneStress::setResponse(...)" << endln;
#endif

  //Response *theResponse =0;
  const char *matType = this->getType();

  output.tag("NdMaterialOutput");
  output.attr("matType", this->getClassType());
  output.attr("matTag", this->getTag());

  if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
    return new MaterialResponse(this, 1, this->getStress());
  else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
    return new MaterialResponse(this, 2, this->getStrain());
  //else if (strcmp(argv[0], "state") == 0)
  //return new MaterialResponse(this, 3, this->getState());
  else
    return 0;
}

int
CSMMRCPlaneStress::getResponse(int responseID, Information &matInfo)
{
#ifdef DEBUG
  opserr << "CSMMRCPlaneStress::getResponse(...)" << endln;
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
CSMMRCPlaneStress::Print(OPS_Stream &s, int flag)
{
  s << "\n\tCSMMRCPlaneStress, material id: " << this->getTag() << endln;

  s << "\tRho: " << rho << endln;
  s << "\tangle1: " << angle1 << endln;
  s << "\tangle2: " << angle2 << endln;
  s << "\trou1: " << rou1 << endln;
  s << "\trou2: " << rou2 << endln;
  s << "\tfpc: " << fpc << endln;
  s << "\tfy: " << fy << endln;
  s << "\tE0: " << E0 << endln;

  s << "Principal Strain: citaStrain = " << citaStrain / 3.14159265453*180.0 << endln;
  s << "Principal Stress: citaStress = " << citaStress / 3.14159265453*180.0 << endln;
  s << " v12 = " << miu12 << " v21 = " << miu21 << endln;
  s << " steelStatus " << steelStatus << endln;
  s << " dirStatus " << dirStatus << endln;
  s << " Damage DOne = " << DDOne << endln;
  s << " Damage DTwo = " << DDTwo << endln;

  s << " G12 = " << G12 << endln;

  int i, j;

  s << "\t call the material print() function : " << endln;

  s << "\t the steel 1 information is : " << endln;
  theMaterial[0]->Print(s, flag);
  s << "\t the steel 2 information is : " << endln;
  theMaterial[1]->Print(s, flag);
  s << "\t the concrete 1 information is : " << endln;
  theMaterial[2]->Print(s, flag);
  s << "\t the concrete 2 information is : " << endln;
  theMaterial[3]->Print(s, flag);

  s << "\tStrain and stress of the uniaxial materials:" << endln;
  for (int i = 0; i < 4; i++) {
    s << "\tUniaxial Material " << i + 1 << " :" << endln;
    s << "\t            Strain : " << theMaterial[i]->getStrain() << endln;
    s << "\t            Stress : " << theMaterial[i]->getStress() << endln;
    s << "\t  Uniaxial Tangent : " << theMaterial[i]->getTangent() << endln;
  }

}

int
CSMMRCPlaneStress::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  int dataTag = this->getDbTag();

  // Packs its data into a Vector and sends this to theChannel
  static Vector data(9);
  data(0) = this->getTag();
  data(1) = rho;
  data(2) = angle1;
  data(3) = angle2;
  data(4) = rou1;
  data(5) = rou2;
  data(6) = fpc;
  data(7) = fy;
  data(8) = E0;

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING CSMMRCPlaneStress::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }

  // Now sends the IDs of its materials
  int matDbTag;

  static ID idData(8);

  // NOTE: to ensure that the material has a database
  // tag if sending to a database channel.

  int i;
  for (i = 0; i < 4; i++) {
    idData(i) = theMaterial[i]->getClassTag();
    matDbTag = theMaterial[i]->getDbTag();
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0)
        theMaterial[i]->setDbTag(matDbTag);
    }
    idData(i + 4) = matDbTag;
  }

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING CSMMRCPlaneStress::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 4; i++) {
    res += theMaterial[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "CSMMRCPlaneStress::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }

  return res;
}

int
CSMMRCPlaneStress::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  // Quad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(9);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING CSMMRCPlaneStress::recvSelf() - failed to receive Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  rho = data(1);
  angle1 = data(2);
  angle2 = data(3);
  rou1 = data(4);
  rou2 = data(5);
  fpc = data(6);
  fy = data(7);
  E0 = data(8);

  static ID idData(8);

  // now receives the tags of its materials
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING CSMMRCPlaneStress::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  if (theMaterial == 0) {
    // Allocate new materials
    theMaterial = new UniaxialMaterial *[4];
    if (theMaterial == 0) {
      opserr << "CSMMRCPlaneStress::recvSelf() - Could not allocate UniaxialMaterial* array\n";
      return -1;
    }
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i + 4);
      // Allocate new material with the sent class tag
      theMaterial[i] = theBroker.getNewUniaxialMaterial(matClassTag);
      if (theMaterial[i] == 0) {
        opserr << "CSMMRCPlaneStress::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
        return -1;
      }
      // Now receive materials into the newly allocated space
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "CSMMRCPlaneStress::recvSelf() - material " << i << "failed to recv itself\n";
        return res;
      }
    }
  }

  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i + 4);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (theMaterial[i]->getClassTag() != matClassTag) {
        delete theMaterial[i];
        theMaterial[i] = theBroker.getNewUniaxialMaterial(matClassTag);
        if (theMaterial[i] == 0) {
          opserr << "CSMMRCPlaneStress::recvSelf() - material " << i << "failed to create\n";
          return -1;
        }
      }
      // Receive the material
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "CSMMRCPlaneStress::recvSelf() - material " << i << "failed to recv itself\n";
        return res;
      }
    }
  }

  return res;
}

int
CSMMRCPlaneStress::determineTrialStress(void)
{
  // Get Principal strain direction first
  Tstrain(0) = strain_vec(0);//epslonx,epslony,0.5*gammaxy	
  Tstrain(1) = strain_vec(1);
  Tstrain(2) = 0.5*strain_vec(2);

  // Get citaR based on Tstrain
  double temp_citaR;
  //double eps = 1e-12;

  if (fabs(Tstrain(0) - Tstrain(1)) < 1e-7)               { citaR = 0.25*PI;             dirStatus = 0; }
  else {// Tstrain[0] != Tstrain[1]
	temp_citaR = 0.5 * atan(fabs(2.0e6*Tstrain(2) / (1.0e6*Tstrain(0) - 1.0e6*Tstrain(1))));
	if (fabs(Tstrain(2)) < 1e-7)                          { citaR = 0;                   dirStatus = 1; }
	else if ((Tstrain(0) > Tstrain(1)) && (Tstrain(2)>0)) { citaR = temp_citaR;          dirStatus = 2; }
	else if ((Tstrain(0) > Tstrain(1)) && (Tstrain(2)<0)) { citaR = PI - temp_citaR;     dirStatus = 3; }
	else if ((Tstrain(0) < Tstrain(1)) && (Tstrain(2)>0)) { citaR = 0.5*PI - temp_citaR; dirStatus = 4; }
	else if ((Tstrain(0) < Tstrain(1)) && (Tstrain(2)<0)) { citaR = 0.5*PI + temp_citaR; dirStatus = 5; }
	else {
	  opserr << "CSMMRCPlaneStress::determineTrialStress: Failure to calculate citaR\n";
	  opserr << " Tstrain(0) = " << Tstrain(0) << endln;
	  opserr << " Tstrain(1) = " << Tstrain(1) << endln;
	  opserr << " Tstrain(2) = " << Tstrain(2) << endln;
	}
  }
  determineConcreteStatus(dirStatus);

  while ((citaR - 0.5*PI) > 1e-8) {
    citaR = citaR - 0.5*PI;
	//dirStatus = 1; 
  }

  citaStrain = citaR;

  int status = 0; // status to check if iteration for principal stress direction
  double tolerance = 0.0088; //PI/360.0; // tolerance for iteration, equals to 0.5 degree
  int iteration_counter = 0;
  double error;

  error = getAngleError(citaR); // first try citaR
  if (error < tolerance)
    status = 1;

  double citaOne = citaR;
  double citaTwo = citaR;
  double minError = error;
  double citaFinal = citaR;
  double citaMinError = citaR;
  //double citaLow = citaR - PI / 2.0; //normally this 45 degree would be enough
  //double citaHigh = citaR + PI / 2.0;

  while ((status == 0) && (citaOne>0-FLT_EPSILON || citaTwo<0.5*PI+FLT_EPSILON)) {
    citaOne = citaOne - PI / 360.0;
    citaTwo = citaTwo + PI / 360.0;

	if (citaOne > 0.0 - FLT_EPSILON) {
      error = getAngleError(citaOne);
      if (minError > error) {
    minError = error;
    citaMinError = citaOne;
      }
      if (error < tolerance) {
    status = 1;
    citaFinal = citaOne;
      }
    }

	if (citaTwo < 0.5*PI + FLT_EPSILON) {
      error = getAngleError(citaTwo);
      if (minError > error) {
    minError = error;
    citaMinError = citaTwo;
      }
      if (error < tolerance) {
	status = 1;
    citaFinal = citaTwo;
      }
    }

    iteration_counter++;
  }

  if (status == 0) {  // does not get converged after iteration
    //citaStress = getPrincipalStressAngle(citaMinError);
	error = getAngleError(citaFinal);
	// if ( minError > 0.05 )
    opserr << "CSMMRCPlaneStress::determineTrialStress(): failure to get converged principal stress. error=" << error/PI*180 << endln;
	//opserr << "citaStrain = " << citaStrain/PI*180. << "бу; citaStress = " << citaStress/PI*180 << "бу." << endln;
  }

  //citaStress = citaFinal;  // assign value for output in the screen

  return 0;
}

double
CSMMRCPlaneStress::getAngleError(double inputCita)
{
  double outputCita = getPrincipalStressAngle(inputCita);

  double error;
  double error1, error2, error3;
  error1 = fabs( inputCita - outputCita);
  error2 = fabs( inputCita - outputCita + 0.5*PI);
  error3 = fabs(-inputCita + outputCita + 0.5*PI);

  if (error1 > error2)
    error = error2;
  else
    error = error1;

  if (error > error3)
    error = error3;

  return error;
}

double
CSMMRCPlaneStress::getPrincipalStressAngle(double inputAngle)
{
  double citaIn = inputAngle; // Trial principal stress direction, obtained from input

  double citaOut = 0.0; // outputAngle, obtained from stresses based on citaIn

  // Define i, j, k, for loop use
  int j;

  // Definition of the variables and matrix
  //Transformation matrix
  Matrix TOne(3, 3);     // T(citaOne)
  Matrix TMOne(3, 3);    // T(minus citaOne)
  Matrix TMSL(3, 3);     // T(minus citaSL)
  Matrix TSL_One(3, 3);  // T(citaSL minus citaOne)
  Matrix TMST(3, 3);     // T(minus citaST)
  Matrix TST_One(3, 3);  // T(citaST minus citaOne)

  Matrix V(3, 3);        //Matrix considering Hsu/Zhu ratios

  //Strain and stress
  Vector tempStrain(3); //temp values of strains
  double stressSL, stressST; // stress of steel layers in L and T direction

  //stiffness of element
  Matrix D(3, 3);      //tangent stiffness matrix
  Matrix DC(3, 3);     //concrete part
  Matrix DSL(3, 3);    //steel part in L direction
  Matrix DST(3, 3);    //steel part in T direction
  Matrix DC_bar(3, 3);   //partial differentiation matrix of concrete part Eq.(49)	
  Matrix tempD(3, 3);  //temp matrix for data transfer

  double epsy = fy / E0;
  double fcr = 0.31*sqrt(fpc);
  double rou = (rou1 < rou2 ? rou1 : rou2);
  //if ( rou < 0.0025 ) rou = 0.0025; //Unified Concrete book, P.273: rou >=0.0015
  if ( rou < 0.0015 ) rou = 0.0015;
  double B = pow((fcr / fy), 1.5) / rou;
  //double fn = fy * (0.91 - 2.0*B) / (0.98 - 0.25*B);  // fn or fs?? P. 258
  double epsn = epsy*(0.91 - 2.0*B) / (0.98 - 0.25*B);

  double citaL = angle1; // angle for direction of steel one
  double citaT = angle2; // angle for direction of steel two

  //Set values for matrix TMSL, TMST
  TMSL(0, 0) = pow(cos(citaL), 2);
  TMSL(0, 1) = pow(sin(citaL), 2);
  TMSL(0, 2) = -2.0*cos(citaL)*sin(citaL);
  TMSL(1, 0) = pow(sin(citaL), 2);
  TMSL(1, 1) = pow(cos(citaL), 2);
  TMSL(1, 2) = 2.0*cos(citaL)*sin(citaL);
  TMSL(2, 0) = cos(citaL)*sin(citaL);
  TMSL(2, 1) = -cos(citaL)*sin(citaL);
  TMSL(2, 2) = pow(cos(citaL), 2) - pow(sin(citaL), 2);

  TMST(0, 0) = pow(cos(citaT), 2);
  TMST(0, 1) = pow(sin(citaT), 2);
  TMST(0, 2) = -2.0*cos(citaT)*sin(citaT);
  TMST(1, 0) = pow(sin(citaT), 2);
  TMST(1, 1) = pow(cos(citaT), 2);
  TMST(1, 2) = 2.0*cos(citaT)*sin(citaT);
  TMST(2, 0) = cos(citaT)*sin(citaT);
  TMST(2, 1) = -cos(citaT)*sin(citaT);
  TMST(2, 2) = pow(cos(citaT), 2) - pow(sin(citaT), 2);

  //Set values for transformation matrix TOne(3,3),TMOne(3,3),TSL_One(3,3),TST_One(3,3)
  TOne(0, 0) = pow(cos(citaIn), 2);
  TOne(0, 1) = pow(sin(citaIn), 2);
  TOne(0, 2) = 2.0*cos(citaIn)*sin(citaIn);
  TOne(1, 0) = pow(sin(citaIn), 2);
  TOne(1, 1) = pow(cos(citaIn), 2);
  TOne(1, 2) = -2.0*cos(citaIn)*sin(citaIn);
  TOne(2, 0) = -cos(citaIn)*sin(citaIn);
  TOne(2, 1) = cos(citaIn)*sin(citaIn);
  TOne(2, 2) = pow(cos(citaIn), 2) - pow(sin(citaIn), 2);

  TMOne(0, 0) = pow(cos(citaIn), 2);
  TMOne(0, 1) = pow(sin(citaIn), 2);
  TMOne(0, 2) = -2.0*cos(citaIn)*sin(citaIn);
  TMOne(1, 0) = pow(sin(citaIn), 2);
  TMOne(1, 1) = pow(cos(citaIn), 2);
  TMOne(1, 2) = 2.0*cos(citaIn)*sin(citaIn);
  TMOne(2, 0) = cos(citaIn)*sin(citaIn);
  TMOne(2, 1) = -cos(citaIn)*sin(citaIn);
  TMOne(2, 2) = pow(cos(citaIn), 2) - pow(sin(citaIn), 2);

  TSL_One(0, 0) = pow(cos(citaL - citaIn), 2);
  TSL_One(0, 1) = pow(sin(citaL - citaIn), 2);
  TSL_One(0, 2) = 2.0*cos(citaL - citaIn)*sin(citaL - citaIn);
  TSL_One(1, 0) = pow(sin(citaL - citaIn), 2);
  TSL_One(1, 1) = pow(cos(citaL - citaIn), 2);
  TSL_One(1, 2) = -2.0*cos(citaL - citaIn)*sin(citaL - citaIn);
  TSL_One(2, 0) = -cos(citaL - citaIn)*sin(citaL - citaIn);
  TSL_One(2, 1) = cos(citaL - citaIn)*sin(citaL - citaIn);
  TSL_One(2, 2) = pow(cos(citaL - citaIn), 2) - pow(sin(citaL - citaIn), 2);

  TST_One(0, 0) = pow(cos(citaT - citaIn), 2);
  TST_One(0, 1) = pow(sin(citaT - citaIn), 2);
  TST_One(0, 2) = 2.0*cos(citaT - citaIn)*sin(citaT - citaIn);
  TST_One(1, 0) = pow(sin(citaT - citaIn), 2);
  TST_One(1, 1) = pow(cos(citaT - citaIn), 2);
  TST_One(1, 2) = -2.0*cos(citaT - citaIn)*sin(citaT - citaIn);
  TST_One(2, 0) = -cos(citaT - citaIn)*sin(citaT - citaIn);
  TST_One(2, 1) = cos(citaT - citaIn)*sin(citaT - citaIn);
  TST_One(2, 2) = pow(cos(citaT - citaIn), 2) - pow(sin(citaT - citaIn), 2);

  // Get strain values from strain of element in x y directions
  /*Tstrain(0) = strain_vec(0);  Tstrain(1) = strain_vec(1);  Tstrain(2) = 0.5*strain_vec(2);*/

  //calculate tempStrain: epslon1,epslon2, 0.5*gamma12 in trial principal stress direction
  tempStrain.addMatrixVector(0.0, TOne, Tstrain, 1.0);

  //double epslonOne, epslonTwo, halfGammaOneTwo;
  epslonOne = tempStrain(0);
  epslonTwo = tempStrain(1);
  halfGammaOneTwo = tempStrain(2);

  //double CSL = theMaterial[0]->getCommittedStrain();
  //double CST = theMaterial[1]->getCommittedStrain();

  theResponses[0]->getResponse();
  theResponses[1]->getResponse();
  Information &theInfoL = theResponses[0]->getInformation();
  Information &theInfoT = theResponses[1]->getInformation();
  double CSL = theInfoL.theDouble;
  double CST = theInfoT.theDouble;

  if ((CSL > epsn) || (CST > epsn)) {
    steelStatus = 1;
  }

  //set v12 and v21 obtained from strain
  double strainSL, strainST; //Biaxial strain of steel in L, T
  double strainSF;           //larger one of strainSL, strainST

  strainSL = pow(cos(citaL), 2)*Tstrain(0) + pow(sin(citaL), 2)*Tstrain(1) + 2.0*sin(citaL)*cos(citaL)*Tstrain(2);
  strainST = pow(cos(citaT), 2)*Tstrain(0) + pow(sin(citaT), 2)*Tstrain(1) + 2.0*sin(citaT)*cos(citaT)*Tstrain(2);

  strainSF = (strainSL > strainST ? strainSL : strainST);

  double v12 = 0.2;
  double v21 = 0.2; //initial values for Hsu/Zhu ratios

  if ((epslonOne > 0.0) && (epslonTwo > 0.0)) {  // both tension
    v12 = 0.0;
    v21 = 0.0;
  }
  else if ((epslonOne > 0.0) && (epslonTwo <= 0.0)) {  // one tension, two compression
    v21 = 0.0;
    if (strainSF > 0.002) {
      v12 = 1.0;
      //v12 = 1.9;
    }
    else if (strainSF < 0) {
      v12 = 0.2;
    }
    else {
      v12 = 0.2 + 400.0*strainSF;
      //v12 = 0.2 + 850.0*strainSF;
    }

    if (steelStatus == 1)
      v12 = 1.0;
    //v12 = 1.9;

  }
  else if ((epslonOne <= 0.0) && (epslonTwo > 0.0)) {  // one compression, two tension
    v12 = 0.0;
    if (strainSF > 0.002) {
      v21 = 1.0;
      //v21 = 1.9;
    }
    else if (strainSF < 0) {
      v21 = 0.2;
    }
    else {
      v21 = 0.2 + 400.0*strainSF;
      //v21 = 0.2 + 850.0*strainSF;
    }

    if (steelStatus == 1)
      v21 = 1.0;
    //v21 = 1.9;

  }
  else if ((epslonOne <= 0.0) && (epslonTwo <= 0.0)) {  //both compression
    if (strainSF > 0.002) {
      v21 = 0.95;
      v12 = 0.95;
      //v21 = 1.9;
      //v12 = 1.9;
    }
    else if (strainSF < 0) {
      v21 = 0.2;
      v12 = 0.2;
    }
    else {
      //v21 = 0.2 + 850.0*strainSF;
      //v12 = 0.2 + 850.0*strainSF;
      v21 = 0.2 + 375.0*strainSF;
      v12 = 0.2 + 375.0*strainSF;
    }

    if (steelStatus == 1) {
      //v21 = 1.9;
      //v12 = 1.9;
      v21 = 0.95;
      v12 = 0.95;
    }

  }
  else {   //error numerical value -- can not happan...
    opserr << "CSMMRCPlaneStress::getPrincipalStressAngle: failure to get Hsu/Zhu ratio!\n";
  }

  miu12 = v12; // record the value for output in screen
  miu21 = v21; // record the value for output in screen

  //set values of matrix V(3,3)
  if (v12*v21 == 1.0) {
    opserr << "CSMMRCPlaneStress::getPrincipalStressAngle: failure to get matrix [V]!\n";
    opserr << "v12= " << v12 << endln;
    opserr << "v21= " << v21 << endln;
    V(0, 0) = 1.0;
	V(0, 1) = 0.0;
	V(0, 2) = 0.0;

	V(1, 0) = 0.0;
    V(1, 1) = 1.0;
	V(1, 2) = 0.0;

	V(2, 0) = 0.0;
	V(2, 1) = 0.0;
    V(2, 2) = 1.0;
  }
  else {
    V(0, 0) = 1.0 / (1.0 - v12*v21);
    V(0, 1) = v12 / (1.0 - v12*v21);
    V(0, 2) = 0.0;
  
    V(1, 0) = v21 / (1.0 - v12*v21);
    V(1, 1) = 1.0 / (1.0 - v12*v21);
    V(1, 2) = 0.0;
  
    V(2, 0) = 0.0;
    V(2, 1) = 0.0;
    V(2, 2) = 1.0;
  }

  //********** get [DC]**************
  DC.addMatrixProduct(0.0, V, TOne, 1.0); // DC = V * TOne;

  //calculate epslon1_bar, epslon2_bar, 0.5*gamma12
  tempStrain(0) = V(0, 0)*epslonOne + V(0, 1)*epslonTwo; //epslon1_bar
  tempStrain(1) = V(1, 0)*epslonOne + V(1, 1)*epslonTwo; //epslon2_bar
  // tempStrain(2) = halfGammaOneTwo;  // hsu-zhu ratio has no effect on this term

  //get stiffness of uniaxial strain of concrete in 12 direction
  //double sigmaOneC; //stress of concrete in 12 direction
  //double sigmaTwoC;
  double GOneTwoC; //shear modulus of concrete in 12 direction

  // get Damage factor: DOne, DTwo
  if (tempStrain(0) < 0) {
    TOneReverseStatus = 0;
    if (TOneNowMaxComStrain > tempStrain(0))
      TOneNowMaxComStrain = tempStrain(0);
  }
  else { // tempStrain(0) > 0
    if (TOneReverseStatus == 0) {     // first reverse from compressive strain
      TOneReverseStatus = 1;
      TOneLastMaxComStrain = COneNowMaxComStrain;
      TOneNowMaxComStrain = 0.0;
    }
  }

  if (tempStrain(1) < 0) {
    TTwoReverseStatus = 0;
    if (TTwoNowMaxComStrain > tempStrain(1))
      TTwoNowMaxComStrain = tempStrain(1);
  }
  else { // tempStrain(1) > 0
    if (TTwoReverseStatus == 0) {      // first reverse from compressive strain
      TTwoReverseStatus = 1;
      TTwoLastMaxComStrain = CTwoNowMaxComStrain;
      TTwoNowMaxComStrain = 0.0;
    }
  }

  double DOne;
  double DTwo;

  if (tempStrain(0) < 0) {
    DOne = 1 - fabs(0.4*TTwoLastMaxComStrain / epsc0);
    if (DOne < 0.2)  DOne = 0.2;
  }

  if (tempStrain(1) < 0) {
    DTwo = 1 - fabs(0.4*TOneLastMaxComStrain / epsc0);
    if (DTwo < 0.2)  DTwo = 0.2;
  }

  //DOne = 1.0; // commented by Ln 
  //DTwo = 1.0;

  DDOne = DOne;  // assign values for screen output
  DDTwo = DTwo;

  //for xx, kk, m, xx=n, kk=delta, keci
  double xx, kk;
//  if (((fabs(lastStress(0)) + fabs(lastStress(1)) + fabs(lastStress(2))) == 0.0) ||
//    (fabs(lastStress(0)) == 0.0)) {
//    xx = 2.0;
//    kk = 1.0;
//  }
//  else {
//    if ((lastStress(0) < 0.0) && (lastStress(1) < 0.0)) {
//      if (lastStress(2) == 0.0) {
//        xx = 2.0;
//        kk = 1.0;
//        //kk = 0;
//      }
//      else {
//        double keci = sqrt((fabs(lastStress(0)) / fabs(lastStress(2)) + 1) * (fabs(lastStress(1)) / fabs(lastStress(2)) + 1));
//        
//        xx = 2 / pow(keci, 3.0);
//        if (xx < 0.6) xx = 0.6;
//        
//        double a, b;
//        if (keci < 1.5) {
//          a = 0;
//          b = 1.0;
//        }
//        else {
//          a = (1.3 - 1.0) / (1.9 - 1.5);
//          b = 1.0 - a*1.5;
//        }
//        kk = a*keci + b;
//        //kk = 0.105*keci+1;
//      }
//    }
//    else if ((lastStress(0) > 0.0) && (lastStress(1) > 0.0))  {
//      kk = 1.0;
//      xx = 2.0;
//    }
//    else { // under tension and compression
//      kk = 1.0;
//      double keciN;
//      if (lastStress(0) < 0) {
//        keciN = sqrt(fabs(lastStress(0)) / fabs(lastStress(2)) + 1);
//      }
//      else {
//        keciN = sqrt(fabs(lastStress(1)) / fabs(lastStress(2)) + 1);
//      }
//      xx = 2 / pow(keciN, 3.0);
//      if (xx < 0.6) xx = 0.6;
//    }
//  }

  xx = 2.0; // for normal cases without axial loads
  kk = 1.0;

  //for Concrete material average responses

  //theMaterial[2]->setTrialStrain(xx, kk, DOne,ita,tempStrain[1],tempStrain[0]); 
  //theMaterial[3]->setTrialStrain(xx, kk, DTwo,ita,tempStrain[0],tempStrain[1]); 

  Information &theInfoC02 = theResponses[2]->getInformation();
  Information &theInfoC03 = theResponses[3]->getInformation();

  double beta = 0.5*atan2(halfGammaOneTwo*2.e6, (epslonOne - epslonTwo)*1.e6); // set beta to 1.0 for simplification, previously
  
  Vector theData(5);
  theData(0) = xx;
  theData(1) = kk;
  theData(3) = fabs(beta)/PI*180.;

  if (isSwapped) {
    theData(2) = DTwo;
    theData(4) = tempStrain(0);   //epslon1_bar
    theInfoC02.setVector(theData);
    theResponses[2]->getResponse();
    
    theData(2) = DOne;
    theData(4) = tempStrain(1);   //epslon2_bar
    theInfoC03.setVector(theData);
    theResponses[3]->getResponse();
    
    theMaterial[2]->setTrialStrain(tempStrain(1), 0.0); //epslon2_bar
    theMaterial[3]->setTrialStrain(tempStrain(0), 0.0); //epslon1_bar
    
    sigmaOneC = theMaterial[3]->getStress();
    sigmaTwoC = theMaterial[2]->getStress();
  }
  else {
    theData(2) = DOne;
	theData(4) = tempStrain(1);   //epslon2_bar
	theInfoC02.setVector(theData);
	theResponses[2]->getResponse();

	theData(2) = DTwo;
	theData(4) = tempStrain(0);   //epslon1_bar
	theInfoC03.setVector(theData);
	theResponses[3]->getResponse();

	theMaterial[2]->setTrialStrain(tempStrain(0), 0.0); //epslon1_bar
	theMaterial[3]->setTrialStrain(tempStrain(1), 0.0); //epslon2_bar

	sigmaOneC = theMaterial[2]->getStress();
	sigmaTwoC = theMaterial[3]->getStress();
  }
  //  else {
  //  opserr << "CSMMRCPlaneStress::(getPrincipalStressAngle) -- wrong dirStatus" << endln; // should not occur
  //}

  // set GOneTwoC = 1.0;
  if (epslonOne == epslonTwo) {
    GOneTwoC = 10000.0; // max value for GOneTwoC
  } 
  else {
    GOneTwoC = fabs((sigmaOneC - sigmaTwoC) / (epslonOne - epslonTwo));
    if (GOneTwoC > 10000.0) // if larger than max value
      GOneTwoC = 10000.0;
  }
  G12 = GOneTwoC; // record the value for output in screen

  // DC_bar(0,1) = theMaterial[2]->getPD();
  // DC_bar(1,0) = theMaterial[3]->getPD();
  theResponses[4]->getResponse();
  Information &theInfoC1 = theResponses[4]->getInformation();

  theResponses[5]->getResponse();
  Information &theInfoC2 = theResponses[5]->getInformation();

  DC_bar(0, 0) = theMaterial[2]->getTangent();
  DC_bar(0, 1) = theInfoC1.theDouble;
  DC_bar(0, 2) = 0.0;

  DC_bar(1, 0) = theInfoC2.theDouble;
  DC_bar(1, 1) = theMaterial[3]->getTangent();
  DC_bar(1, 2) = 0.0;

  DC_bar(2, 0) = 0.0;
  DC_bar(2, 1) = 0.0;
  DC_bar(2, 2) = GOneTwoC;

  //before [DC]=[v]*[TOne], now update [DC]=[Dc_bar]*[V]*[TOne]
  tempD.addMatrixProduct(0.0, DC_bar, DC, 1.0); // before here, DC = V * TOne;
  DC = tempD;

  //update [DC]=[TMOne]*[Dc_bar]*[V]*[TOne]
  tempD.addMatrixProduct(0.0, TMOne, DC, 1.0); // before here, DC = DC_Bar * V * TOne;
  DC = tempD;

  //***************** get [DSL] ******************
  //get [DSL]=[V][TOne]
  DSL.addMatrixProduct(0.0, V, TOne, 1.0);

  //get [DSL]=[TLMOne][V][TOne]
  tempD.addMatrixProduct(0.0, TSL_One, DSL, 1.0);
  DSL = tempD;

  double strainSL_b; //uniaxial strain of steel in L direction
  strainSL_b = DSL(0, 0)*Tstrain(0) + DSL(0, 1)*Tstrain(1) + DSL(0, 2)*Tstrain(2);

  //get stiffness
  double tangentSL;
  theMaterial[0]->setTrialStrain(strainSL_b);
  tangentSL = theMaterial[0]->getTangent();
  stressSL = theMaterial[0]->getStress();

  for (j = 0; j < 3; j++) {
    DSL(0, j) *= (rou1*tangentSL);
    DSL(1, j) = 0.0;
    DSL(2, j) = 0.0;
  }

  //get [DSL]=[TML][Dsl][TLMOne][V][TOne]
  tempD.addMatrixProduct(0.0, TMSL, DSL, 1.0);
  DSL = tempD;

  //**************** get [DST] ****************     
  //get [DST]=[V][TOne]
  DST.addMatrixProduct(0.0, V, TOne, 1.0);

  //get [DST]=[TST_One][V][TOne]
  tempD.addMatrixProduct(0.0, TST_One, DST, 1.0);
  DST = tempD;

  double strainST_b; //uniaxial strain of steel in T direction
  strainST_b = DST(0, 0)*Tstrain(0) + DST(0, 1)*Tstrain(1) + DST(0, 2)*Tstrain(2);

  double tangentST;
  theMaterial[1]->setTrialStrain(strainST_b);
  tangentST = theMaterial[1]->getTangent();
  stressST = theMaterial[1]->getStress();

  for (j = 0; j < 3; j++) {
    DST(0, j) *= (rou2*tangentST);
    DST(1, j) = 0.0;
    DST(2, j) = 0.0;
  }

  //get [DST]=[TMT][Dst][TTMOne][V][TOne]
  tempD.addMatrixProduct(0.0, TMST, DST, 1.0);
  DST = tempD;

  //****************** get tangent_matrix  ****************    
  // Get tangent_matrix
  tangent_matrix = DC + DSL + DST;

  tangent_matrix(0, 2) *= 0.5;
  tangent_matrix(1, 2) *= 0.5;
  tangent_matrix(2, 2) *= 0.5;

  //**************** get Tstress and stress_vec ****************
  Tstress(0) = pow(cos(citaIn), 2)*sigmaOneC + pow(sin(citaIn), 2)*sigmaTwoC
    - 2 * sin(citaIn)*cos(citaIn)*halfGammaOneTwo*GOneTwoC
    + pow(cos(citaL), 2)*rou1*stressSL + pow(cos(citaT), 2)*rou2*stressST;

  Tstress(1) = pow(sin(citaIn), 2)*sigmaOneC + pow(cos(citaIn), 2)*sigmaTwoC
    + 2 * sin(citaIn)*cos(citaIn)*halfGammaOneTwo*GOneTwoC
    + pow(sin(citaL), 2)*rou1*stressSL + pow(sin(citaT), 2)*rou2*stressST;

  Tstress(2) = cos(citaIn)*sin(citaIn)*sigmaOneC - cos(citaIn)*sin(citaIn)*sigmaTwoC
    + (pow(cos(citaIn), 2) - pow(sin(citaIn), 2))*halfGammaOneTwo*GOneTwoC
    + cos(citaL)*sin(citaL)*rou1*stressSL + cos(citaT)*sin(citaT)*rou2*stressST;

  stress_vec = Tstress;

  // get calculated principal stress direction citaOut
  double temp_citaOut;
  if ( fabs(Tstress(0)-Tstress(1)) < 1e-7 )                 citaOut = 0.25*PI;
  else { // Tstrain(0) != Tstrain(1)
	temp_citaOut = 0.5 * atan(fabs(2.0e6*Tstress[2] / (1.0e6*Tstress(0) - 1.0e6*Tstress(1))));
  	     if ( fabs(Tstress(2)) < 1e-7 )                     citaOut = 0;
	else if ((Tstress(0) > Tstress(1)) && (Tstress(2) > 0))	citaOut = temp_citaOut;
  	else if ((Tstress(0) > Tstress(1)) && (Tstress(2) < 0)) citaOut = PI - temp_citaOut;
  	else if ((Tstress(0) < Tstress(1)) && (Tstress(2) > 0)) citaOut = 0.5*PI - temp_citaOut;
  	else if ((Tstress(0) < Tstress(1)) && (Tstress(2) < 0)) citaOut = 0.5*PI + temp_citaOut;
  	else {
      opserr << "CSMMRCPlaneStress::getPrincipalStressAngle: Failure to calculate principal stress direction\n";
      opserr << " Tstress(0) = " << Tstress(0) << endln;
	  opserr << " Tstress(1) = " << Tstress(1) << endln;
	  opserr << " Tstress(2) = " << Tstress(2) << endln;
  	}
  }
  
  while ((citaOut - 0.5*PI) > 1e-8) {
    citaOut = citaOut - 0.5*PI;
  }

  citaStress = citaOut;

  return citaOut;
}

void
CSMMRCPlaneStress::determineConcreteStatus(int tempStatus)
{
	bool temp = 0;
	if (tempStatus != lastDirStatus) {
		switch (lastDirStatus) {
		case 0:
			temp = isSwapped;
			break;

		case 1:
			if (tempStatus == 3 || tempStatus == 4) {
				if (isSwapped)
					temp = 0;
				else
					temp = 1;
			}
			else {
				temp = isSwapped;
			}
			break;

		case 2:
			if (tempStatus == 3) {
				if (isSwapped)
					temp = 0;
				else
					temp = 1;
			}
			else {
				temp = isSwapped;
			}
			break;

		case 3:
			if (tempStatus == 1 || tempStatus == 2 || tempStatus == 4) {
				if (isSwapped)
					temp = 0;
				else
					temp = 1;
			}
			else {
				temp = isSwapped;
			}
			break;

		case 4:
			if (tempStatus == 1 || tempStatus == 3 || tempStatus == 5) {
				if (isSwapped)
					temp = 0;
				else
					temp = 1;
			}
			else {
				temp = isSwapped;
			}
			break;

		case 5:
			if (tempStatus == 4) {
				if (isSwapped)
					temp = 0;
				else
					temp = 1;
			}
			else {
				temp = isSwapped;
			}
			break;

		default:
			opserr << "error to determine the concrete status" << endln;
			break;
		}
		isSwapped = temp;
	}
}