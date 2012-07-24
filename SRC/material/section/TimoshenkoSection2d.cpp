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
                                                                        
// $Revision: 1.2 $
// $Date: 2010/02/04 19:10:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/TimoshenkoSection2d.cpp,v $
                                                                        
// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition for 
// TimoshenkoSection2d.h. TimoshenkoSection2d provides the abstraction of a 
// rectangular section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.
// The fiber stresses are the 11, 12, and 13 components of stress, from
// which all six beam stress resultants are obtained.

#include <stdlib.h>
#include <math.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <Fiber.h>
#include <classTags.h>
#include <TimoshenkoSection2d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <NDMaterial.h>

#include <string.h>

ID TimoshenkoSection2d::code(3);

// constructors:
TimoshenkoSection2d::TimoshenkoSection2d(int tag, int num, Fiber **fibers):
  SectionForceDeformation(tag, SEC_TAG_TimoshenkoSection2d),
  numFibers(num), theMaterials(0), matData(0), yh(0.0), zh(0.0),
  yBar(0.0), zBar(0.0), e(3), eCommit(3), s(0), ks(0)
{
  if (numFibers != 0) {
  
    theMaterials = new NDMaterial*[numFibers];
    if (theMaterials == 0) {
 	   opserr << "TimoshenkoSection2d::TimoshenkoSection2d -- failed to allocate Material pointers\n";
 	   exit(-1);
    }
    
    matData = new double [numFibers*3];
    if (matData == 0) {
 	   opserr << "TimoshenkoSection2d::TimoshenkoSection2d -- failed to allocate double array for material data\n";
 	   exit(-1);
    }
    
    double Qz = 0.0;
    double Qy = 0.0;
    double a  = 0.0;
    double yHmin=0, zHmin=0, yHmax=0, zHmax=0;
	NDMaterial *theMat;
    for (int i = 0; i < numFibers; i++) {
	  double yLoc, zLoc, Area;
	  fibers[i]->getFiberLocation(yLoc, zLoc);
	  Area = fibers[i]->getArea();
	  theMat = fibers[i]->getNDMaterial();
	  if (theMat == 0) {
	    opserr << "TimoshenkoSection3d::TimoshenkoSection3d -- failed to get fiber information" << endln;
	  }
      
 	  Qz += yLoc*Area;
 	  Qy += zLoc*Area;
 	  a  += Area;
      
 	  matData[i*3] = yLoc;
 	  matData[i*3+1] = zLoc;
 	  matData[i*3+2] = Area;
      
      theMaterials[i] = theMat->getCopy("BeamFiber2d") ; // theMat.getCopy("TimoshenkoFiber");
      if (theMaterials[i] == 0)
         opserr << "TimoshenkoSection2d::TimoshenkoSection2d -- failed to get copy of beam fiber" << endln;

	  if (yLoc < yHmin ) yHmin=yLoc;
	  if (zLoc < zHmin ) zHmin=zLoc;
	  if (yLoc > yHmax ) yHmax=yLoc;
	  if (zLoc > zHmax ) zHmax=zLoc;
    }
    yBar = Qz/a;
    zBar = Qy/a;
	zh   = yHmax - yHmin;
	yh   = zHmax - zHmin;
  } 
  
  s = new Vector(sData, 3);
  ks = new Matrix(kData, 3, 3);

  for (int i=0; i<9; i++)
	kData[i] = 0.0;
  for (int i=0; i<3; i++)
	sData[i] = 0.0;
  
  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_VY;

}

// constructor for blank object that recvSelf needs to be invoked upon
TimoshenkoSection2d::TimoshenkoSection2d():
  SectionForceDeformation(0, SEC_TAG_TimoshenkoSection2d),
  numFibers(0), theMaterials(0), matData(0), yh(0.0), zh(0.0),
  yBar(0.0), zBar(0.0), e(3), eCommit(3), s(0), ks(0)
{
  s = new Vector(sData, 3);
  ks = new Matrix(kData, 3, 3);
  
  for (int i=0; i<9; i++)
	kData[i] = 0.0;
  for (int i=0; i<3; i++)
	sData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_VY;
}

// destructor:
TimoshenkoSection2d::~TimoshenkoSection2d()
{
  if (theMaterials != 0) {
  	for (int i = 0; i < numFibers; i++)
  	  if (theMaterials[i] != 0)
  	    delete theMaterials[i];
  
  	delete [] theMaterials;
  }
  
  if (matData != 0)
  	delete [] matData;
  
  if (s != 0)
  	delete s;
  
  if (ks != 0)
	delete ks;
}

int
TimoshenkoSection2d::addFiber(Fiber &newFiber)
{
  // need to create larger arrays
  int newSize = numFibers+1;
  NDMaterial **newArray = new NDMaterial *[newSize]; 
  double *newMatData = new double [3 * newSize];
  if (newArray == 0 || newMatData == 0) {
    opserr <<"TimoshenkoSection2d::addFiber -- failed to allocate Fiber pointers\n";
    return -1;
  }

  // copy the old pointers and data
  int i;
  for (i = 0; i < numFibers; i++) {
    newArray[i] = theMaterials[i];
    newMatData[3*i] = matData[3*i];
    newMatData[3*i+1] = matData[3*i+1];
	newMatData[3*i+2] = matData[3*i+2];
  }

  // set the new pointers and data
  double yLoc, zLoc, Area;
  newFiber.getFiberLocation(yLoc, zLoc);
  Area = newFiber.getArea();
  newMatData[numFibers*3] = yLoc;
  newMatData[numFibers*3+1] = zLoc;
  newMatData[numFibers*3+2] = Area;
  NDMaterial *theMat = newFiber.getNDMaterial();
  newArray[numFibers] = theMat->getCopy("BeamFiber2d");

  if (newArray[numFibers] == 0) {
    opserr <<"TimoshenkoSection2d::addFiber -- failed to get copy of a Material\n";
    delete [] newMatData;
    return -1;
  }

  numFibers++;

  if (theMaterials != 0) {
    delete [] theMaterials;
    delete [] matData;
  }

  theMaterials = newArray;
  matData = newMatData;

  double Qz = 0.0;
  double Qy = 0.0;
  double A  = 0.0;

  // Recompute centroid
  for (i = 0; i < numFibers; i++) {
    yLoc = matData[3*i];
	zLoc = matData[3*i+1];
	Area = matData[3*i+2];
    A  += Area;
    Qz += yLoc*Area;
	Qy += zLoc*Area;
  }

  yBar = Qz/A;
  zBar = Qy/A;

  return 0;
}

// Compute fiber strain, eps, from section deformation, e, using the
// linear kinematic relationship, eps = a*e, where
// eps = [eps_11, eps_12 ]'
// e   = [eps_a, kappa_z, gamma_y ]'
// a   = [1 -y*beta  0   
//        0     0    1 ]
int TimoshenkoSection2d::setTrialSectionDeformation (const Vector &deforms)
{
  int res = 0;
  e = deforms;
  
  // Material strain vector
  static Vector eps(2);
  
  double y, z;
  
  double beta = 1.0; //5.0/6.0;

  for (int i = 0; i < numFibers; i++) {
    y = matData[i*3] - yBar;
    z = matData[i*3+1] - zBar;
      
    eps(0) = e(0) - y * e(1);
    eps(1) = e(2);

	res += theMaterials[i]->setTrialStrain(eps);
  }

  return res;
}

const Vector&
TimoshenkoSection2d::getSectionDeformation(void)
{
  return e;
}

double
TimoshenkoSection2d::getZh(void)
{
  double yHmin=0, zHmin=0, yHmax=0, zHmax=0;
  double yLoc; //, zLoc
  for (int i = 0; i < numFibers; i++) {
	yLoc=matData[i*3];
 	//zLoc=matData[i*3+1];

	if (yLoc < yHmin ) yHmin=yLoc;
	//if (zLoc < zHmin ) zHmin=zLoc;
	if (yLoc > yHmax ) yHmax=yLoc;
	//if (zLoc > zHmax ) zHmax=zLoc;
  }

  zh = yHmax - yHmin;
  //yh = zHmax - zHmin;

  return zh;
}

double
TimoshenkoSection2d::getEIz(void)
{
  double G, E, K, EIz = 0.;
  double y, z, A;
  for (int i = 0; i < numFibers; i++) {
    y = matData[i*3] - yBar;
    z = matData[i*3+1] - zBar;
	A = matData[i*3+2];

	const Matrix &Dt = theMaterials[i]->getTangent();

	//G = Dt(1,1); K= Dt(0,0) - 4./3.*G;
	E = Dt(0,0); //9.*K*G/(3.*K+G);
	EIz += E * A * pow(y,2.);
  }

  return EIz;
}

double
TimoshenkoSection2d::getGAy(void)
{
  double G, GAy = 0.;
  double y, z, A;
  for (int i = 0; i < numFibers; i++) {
	//y =matData[i*3] - yBar;
 	//z =matData[i*3+1] - zBar;
	A =matData[i*3+2];

	const Matrix &Dt = theMaterials[i]->getTangent();
	G = Dt(1,1);
	GAy += G * A;
  }

  return GAy;
}

// Compute section tangent stiffness, ks, from material tangent, Dt,
// by integrating over section area
// eps = [eps_11, eps_12 ]'
// e   = [eps_a, kappa_z, gamma_y ]'
// ks = int_A a'*Dt*a dA
// a   = [1 -y  0   
//        0  0  1 ]
const Matrix&
TimoshenkoSection2d::getSectionTangent(void)
{
  for (int i=0; i<9; i++)
	kData[i] = 0.0;
 
  double y, z, w;
  double y2, z2, yz;
  double d00, d01;
  double d10, d11;
  double tmp;
  //double five6  = 5./6.;
  //double root56 = sqrt(five6);
  
  for (int i = 0; i < numFibers; i++) {
    
    y = matData[i*3] - yBar;
    z = matData[i*3+1] - zBar;
    w = matData[i*3+2];

    y2 = y*y;
    z2 = z*z;
    yz = y*z;
    
    const Matrix &Dt = theMaterials[i]->getTangent();

	d00 = Dt(0,0)*w; d01 = Dt(0,1)*w; 
	d10 = Dt(1,0)*w; d11 = Dt(1,1)*w; 
    
    kData[0] += d00;    //0,0           P
    kData[4] += y2*d00; //1,1           M
    kData[8] += d11;    //2,2 five6*    V

    tmp = -y*d00;
    kData[1] += tmp; // 0,1
    kData[3] += tmp; // 1,0
    
	// Hit tangent terms with root56
	//d01 *= root56; d10 *= root56;

    kData[2] += d01; //0,2
    kData[6] += d10; //2,0

    kData[5] -= y*d01; //1,2
    kData[7] -= y*d10; //2,1
    
  }

  return *ks;
}

// Compute section stress resultant, s, from material stress, sig,
// by integrating over section area
// s = int_A a'*sig dA
// s = [P, Mz, Vy ]'
// a   = [1 -y  0   
//        0  0  1 ]
const Vector&
TimoshenkoSection2d::getStressResultant(void)
{
  for (int i=0; i<3; i++)
	sData[i] = 0.0;
  
  double y, z, w;
  double sig0, sig1;
  
  //double root56 = sqrt(5./6.);

  for (int i = 0; i < numFibers; i++) {
    
    y = matData[i*3] - yBar;
    z = matData[i*3+1] - zBar;
    w = matData[i*3+2];
    
    const Vector &sig = theMaterials[i]->getStress();
    
    sig0 = sig(0)*w;
    sig1 = sig(1)*w;
    
    sData[0] += sig0;
    sData[1] -= y*sig0;
    sData[2] += sig1;
  }

  return *s;
}

SectionForceDeformation*
TimoshenkoSection2d::getCopy(void)
{

  TimoshenkoSection2d *theCopy = new TimoshenkoSection2d();
  theCopy->setTag(this->getTag());
  theCopy->numFibers = numFibers;

  if (numFibers != 0) {
    theCopy->theMaterials = new NDMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
      opserr << "TimoshenkoSection2d::TimoshenkoSection2d -- failed to allocate Material pointers\n";
      exit(-1);			    
    }

    theCopy->matData = new double [numFibers*3];

    if (theCopy->matData == 0) {
      opserr << "TimoshenkoSection2d::TimoshenkoSection2d -- failed to allocate double array for material data\n";
      exit(-1);
    }
			    
    
    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[i*3] = matData[i*3];
      theCopy->matData[i*3+1] = matData[i*3+1];
      theCopy->matData[i*3+2] = matData[i*3+2];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy("BeamFiber2d");

      if (theCopy->theMaterials[i] == 0) {
	opserr << "TimoshenkoSection2d::getCopy -- failed to get copy of a Material\n";
	exit(-1);
      }
    }    
  }
  theCopy->eCommit = eCommit;
  theCopy->e = e;
  theCopy->yBar = yBar;
  theCopy->zBar = zBar;
  
  for (int i=0; i<9; i++)
	theCopy->kData[i] = kData[i];
  for (int i=0; i<3; i++)
	theCopy->sData[i] = sData[i];

  return theCopy;
}

const ID&
TimoshenkoSection2d::getType(void)
{
  return code;
}

int
TimoshenkoSection2d::getOrder(void) const
{
  return 3;
}

int
TimoshenkoSection2d::commitState(void)
{
  int err = 0;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->commitState();
  
  eCommit = e;

  return err;
}

int
TimoshenkoSection2d::revertToLastCommit(void)
{
  int err = 0;

  // Last committed section deformations
  e = eCommit;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->revertToLastCommit();
  
  return err;
}

int
TimoshenkoSection2d::revertToStart(void)
{
  int err = 0;
  
  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->revertToStart();
  
  return err;
}

int
TimoshenkoSection2d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // create an id to send objects tag and numFibers, 
  //     size 3 so no conflict with matData below if just 1 fiber
  static ID data(3);
  data(0) = this->getTag();
  data(1) = numFibers;
  int dbTag = this->getDbTag();
  res += theChannel.sendID(dbTag, commitTag, data);
  if (res < 0) {
    opserr <<  "TimoshenkoSection2d::sendSelf - failed to send ID data\n";
    return res;
  }    

  if (numFibers != 0) {
    
    // create an id containingg classTag and dbTag for each material & send it
    ID materialData(2*numFibers);
    for (int i=0; i<numFibers; i++) {
      NDMaterial *theMat = theMaterials[i];
      materialData(2*i) = theMat->getClassTag();
      int matDbTag = theMat->getDbTag();
      if (matDbTag == 0) {
	matDbTag = theChannel.getDbTag();
	if (matDbTag != 0)
	  theMat->setDbTag(matDbTag);
      }
      materialData(2*i+1) = matDbTag;
    }    
    
    res += theChannel.sendID(dbTag, commitTag, materialData);
    if (res < 0) {
      opserr <<  "TimoshenkoSection2d::sendSelf - failed to send material data\n";
      return res;
    }    

    // send the fiber data, i.e. area and loc
    Vector fiberData(matData, 3*numFibers);
    res += theChannel.sendVector(dbTag, commitTag, fiberData);
    if (res < 0) {
      opserr <<  "TimoshenkoSection2d::sendSelf - failed to send material data\n";
      return res;
    }    

    // now invoke send(0 on all the materials
    for (int j=0; j<numFibers; j++)
      theMaterials[j]->sendSelf(commitTag, theChannel);

  }

  return res;
}

int
TimoshenkoSection2d::recvSelf(int commitTag, Channel &theChannel,
		     FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static ID data(3);
  
  int dbTag = this->getDbTag();
  res += theChannel.recvID(dbTag, commitTag, data);
  if (res < 0) {
    opserr <<  "TimoshenkoSection2d::recvSelf - failed to recv ID data\n";
    return res;
  }    
  this->setTag(data(0));

  // recv data about materials objects, classTag and dbTag
  if (data(1) != 0) {
    ID materialData(2*data(1));
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
      opserr <<  "TimoshenkoSection2d::recvSelf - failed to recv material data\n";
      return res;
    }    

    // if current arrays not of correct size, release old and resize
    if (theMaterials == 0 || numFibers != data(1)) {
      // delete old stuff if outa date
      if (theMaterials != 0) {
	for (int i=0; i<numFibers; i++)
	  delete theMaterials[i];
	delete [] theMaterials;
	if (matData != 0)
	  delete [] matData;
	matData = 0;
	theMaterials = 0;
      }

      // create memory to hold material pointers and fiber data
      numFibers = data(1);
      if (numFibers != 0) {
	theMaterials = new NDMaterial *[numFibers];
	
	if (theMaterials == 0) {
	  opserr <<"TimoshenkoSection2d::recvSelf -- failed to allocate Material pointers\n";
	  exit(-1);
	}
	
	for (int j=0; j<numFibers; j++)
	  theMaterials[j] = 0;

	matData = new double [numFibers*3];

	if (matData == 0) {
	  opserr <<"TimoshenkoSection2d::recvSelf  -- failed to allocate double array for material data\n";
	  exit(-1);
	}
      }
    }

    Vector fiberData(matData, 3*numFibers);
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
      opserr <<  "TimoshenkoSection2d::recvSelf - failed to recv material data\n";
      return res;
    }    

    int i;
    for (i=0; i<numFibers; i++) {
      int classTag = materialData(2*i);
      dbTag = materialData(2*i+1);

      // if material pointed to is blank or not of corrcet type, 
      // release old and create a new one
      if (theMaterials[i] == 0)
	theMaterials[i] = theBroker.getNewNDMaterial(classTag);
      else if (theMaterials[i]->getClassTag() != classTag) {
	delete theMaterials[i];
	theMaterials[i] = theBroker.getNewNDMaterial(classTag);      
      }

      if (theMaterials[i] == 0) {
	opserr <<"TimoshenkoSection2d::recvSelf -- failed to allocate double array for material data\n";
	exit(-1);
      }

      theMaterials[i]->setDbTag(dbTag);
      res += theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    double Qz = 0.0;
	double Qy = 0.0;
    double A  = 0.0;
    double yLoc, zLoc, Area;

    // Recompute centroid
    for (i = 0; i < numFibers; i++) {
      yLoc = matData[3*i];
      zLoc = matData[3*i+1];
      Area = matData[3*i+2];
      A  += Area;
      Qz += yLoc*Area;
	  Qy += zLoc*Area;
    }
    
    yBar = Qz/A;
    zBar = Qy/A;
  }    

  return res;
}

void
TimoshenkoSection2d::Print(OPS_Stream &s, int flag)
{
  s << "\nTimoshenkoSection2d, tag: " << this->getTag() << endln;
  for (int i=0;i<numFibers; i++) {
    s << "\tFiber Material, tag: " << theMaterials[i]->getTag() << endln;
    theMaterials[0]->Print(s, flag);
  }
}

Response*
TimoshenkoSection2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  const ID &type = this->getType();
  int typeSize = this->getOrder();

  Response *theResponse =0;

  output.tag("SectionOutput");
  output.attr("secType", this->getClassType());
  output.attr("secTag", this->getTag());

  // deformations
  if (strcmp(argv[0],"deformations") == 0 || strcmp(argv[0],"deformation") == 0) {
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","kappaZ");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","eps");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","gammaY");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    theResponse =  new MaterialResponse(this, 1, this->getSectionDeformation());
  
  // forces
  } else if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0) {
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","Mz");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","P");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","Vy");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    theResponse =  new MaterialResponse(this, 2, this->getStressResultant());
  
  // force and deformation
  } else if (strcmp(argv[0],"forceAndDeformation") == 0) { 
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","kappaZ");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","eps");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","gammaY");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    for (int j=0; j<typeSize; j++) {
      int code = type(j);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","Mz");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","P");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","Vy");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }

    theResponse =  new MaterialResponse(this, 4, Vector(2*this->getOrder()));
  
  } else if (strcmp(argv[0],"fiberData") == 0) {
    int numData = numFibers*5;
    for (int j = 0; j < numFibers; j++) {
      output.tag("FiberOutput");
      output.attr("yLoc", matData[3*j]);
      output.attr("zLoc", matData[3*j+1]);
      output.attr("area", matData[3*j+2]);    
      output.tag("ResponseType","yCoord");
      output.tag("ResponseType","zCoord");
      output.tag("ResponseType","area");
      output.tag("ResponseType","stress");
      output.tag("ResponseType","strain");
      output.endTag();
    }
    Vector theResponseData(numData);
    theResponse = new MaterialResponse(this, 5, theResponseData);
  }


  else {
    if (argc > 2 || strcmp(argv[0],"fiber") == 0) {
      
      int key = numFibers;
      int passarg = 2;
      
      
      if (argc <= 3)	{  // fiber number was input directly
	
	key = atoi(argv[1]);
	
      } else if (argc > 4) {         // find fiber closest to coord. with mat tag
	int matTag = atoi(argv[3]);
	double yCoord = atof(argv[1]);
	double zCoord = atof(argv[2]);
	double closestDist = 0.0;
	double ySearch, zSearch, dy, dz;
	double distance;
	int j;
	
	// Find first fiber with specified material tag
	for (j = 0; j < numFibers; j++) {
	  if (matTag == theMaterials[j]->getTag()) {
	    ySearch =  matData[3*j];
	    zSearch =  matData[3*j+1];
	    dy = ySearch-yCoord;
	    dz = zSearch-zCoord;
	    closestDist = sqrt(dy*dy + dz*dz);
	    key = j;
	    break;
	  }
	}
	
	// Search the remaining fibers
	for ( ; j < numFibers; j++) {
	  if (matTag == theMaterials[j]->getTag()) {
	    ySearch =  matData[3*j];
	    zSearch =  matData[3*j+1];
	    dy = ySearch-yCoord;
	    dz = zSearch-zCoord;
	    distance = sqrt(dy*dy + dz*dz);
	    if (distance < closestDist) {
	      closestDist = distance;
	      key = j;
	    }
	  }
	}
	passarg = 4;
      }
      
      else {                  // fiber near-to coordinate specified
	double yCoord = atof(argv[1]);
	double zCoord = atof(argv[2]);
	double closestDist;
	double ySearch, zSearch, dy, dz;
	double distance;
	ySearch =  matData[0];
	zSearch =  matData[1];
	dy = ySearch-yCoord;
	dz = zSearch-zCoord;
	closestDist = sqrt(dy*dy + dz*dz);
	key = 0;
	for (int j = 1; j < numFibers; j++) {
	  ySearch =  matData[3*j];
	  zSearch =  matData[3*j+1];
	  dy = ySearch-yCoord;
	  dz = zSearch-zCoord;
	  distance = sqrt(dy*dy + dz*dz);
	  if (distance < closestDist) {
	    closestDist = distance;
	    key = j;
	  }
	}
	passarg = 3;
      }
      
      if (key < numFibers && key >= 0) {
	output.tag("FiberOutput");
	output.attr("yLoc",matData[3*key]);
	output.attr("zLoc",matData[3*key+1]);
	output.attr("area",matData[3*key+2]);
	
	theResponse =  theMaterials[key]->setResponse(&argv[passarg], argc-passarg, output);
	
	output.endTag();
      }
    }
  }

  output.endTag();
  return theResponse;
}

int 
TimoshenkoSection2d::getResponse(int responseID, Information &sectInfo)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  if (responseID == 5) {
    int numData = 5*numFibers;
    Vector data(numData);
    int count = 0;
    for (int j = 0; j < numFibers; j++) {
      double yLoc, zLoc, A;
	  Vector stress, strain;
      yLoc =  matData[3*j];
      zLoc =  matData[3*j+1];
      A = matData[3*j+2];
      stress = theMaterials[j]->getStress();
      strain = theMaterials[j]->getStrain();
      data(count) = yLoc; data(count+1) = zLoc; data(count+2) = A;
      data(count+3) = stress(0); data(count+4) = stress(1); data(count+5) = stress(2);
	  data(count+6) = strain(0); data(count+7) = strain(1); data(count+8) = strain(2);
      count += 8;
    }
    return sectInfo.setVector(data);	
  } else {
    return SectionForceDeformation::getResponse(responseID, sectInfo);
  }
}

int
TimoshenkoSection2d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 3)
    return -1;


  int result = 0;

  // A material parameter
  if (strstr(argv[0],"material") != 0) {

    // Get the tag of the material
    int paramMatTag = atoi(argv[1]);

    // Loop over fibers to find the right material(s)
    int ok = 0;
    for (int i = 0; i < numFibers; i++)
      if (paramMatTag == theMaterials[i]->getTag()) {
	ok = theMaterials[i]->setParameter(&argv[2], argc-2, param);
	if (ok != -1)
	  result = ok;
      }
    
    return result;
  }    

  int ok = 0;
  
  // loop over every material
  for (int i = 0; i < numFibers; i++) {
    ok = theMaterials[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  return result;
}
