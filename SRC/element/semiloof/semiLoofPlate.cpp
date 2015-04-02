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

// $Revision: 1.1.1.1 $
// $Date: 2000/09/15 08:23:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/femski/semiLoofPlate.cpp,v $

//
// Description: This file contains the implementation for the semiLoofPlate class.
// semi-loof plate

#include "semiLoofPlate.h"
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <ShellMITC4.h>
#include <R3vectors.h>
#include <Renderer.h>
#include <ElementResponse.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

static int numSemiLoofPlate = 0;

void *
OPS_NewSemiLoofPlate(void)
{
  if (numSemiLoofPlate == 0) {
    //    opserr << "Using SemiLoofPlate - Developed by: Ning Li, neallee@tju.edu.cn\n";
    numSemiLoofPlate++;
  }

  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 11) {
    opserr << "Want: element SemiLoofPlate eleTag? node1? node2? .... node9? secTag?";
    return 0;
  }

  int iData[11];
  int numData = 11;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element SemiLoofPlate \n";
    return 0;
  }

  SectionForceDeformation *theSection = OPS_GetSectionForceDeformation(iData[10]);
  if (theSection == 0) {
    opserr << "ERROR:  element SemiLoofPlate " << iData[0] << "section " << iData[10] << " not found\n";
    return 0;
  }

  double rho = 0.0;
  if (numArgs > 11) {
    double dData[1];
    numData = 1;
    if (OPS_GetDouble(&numData, dData) != 0) {
      opserr << "WARNING invalid rho parameter of element SemiLoofPlate\n";
      return 0;
    }
    rho = dData[0];
  }

  theElement = new semiLoofPlate(iData[0], iData[1], iData[2], iData[3],
    iData[4], iData[5], iData[6], iData[7], iData[8], iData[9], *theSection, rho);

  return theElement;
}

//static data
Matrix  semiLoofPlate::stiff(54, 54);  //??
Vector  semiLoofPlate::resid(54);
Matrix  semiLoofPlate::mass(54, 54);

semiLoofPlate::semiLoofPlate(int tag,
  int nd1, int nd2, int nd3, int nd4, int nd5, int nd6, int nd7, int nd8, int nd9,
  SectionForceDeformation &theMaterial, double rho)
  :semiLoofElement(tag, ELE_TAG_semiLoofPlate, 2, 3, 9, 2, 2, 0, 0),
  theLoad(0), Ki(0)
{
  for (int i = 0; i < 9; i++) {
    materialPointers[i] = theMaterial.getCopy();
    if (materialPointers[i] == 0) {
      opserr << "semiLoofPlate::constructor - failed to get a material of type: PlateSection\n";
    } //end if
  } //end for i

  (*data)(0) = rho;

  (*connectedNodes)(0) = nd1;
  (*connectedNodes)(1) = nd2;
  (*connectedNodes)(2) = nd3;
  (*connectedNodes)(3) = nd4;
}

semiLoofPlate::semiLoofPlate()
  :semiLoofElement(ELE_TAG_semiLoofPlate)
{
  // does nothing
}

semiLoofPlate::~semiLoofPlate()
{
  if (theNodes != 0)
    delete[] theNodes;
  if (theLoad != 0)
    delete theLoad;
  if (Ki != 0)
    delete Ki;
}

Node **
semiLoofPlate::getNodePtrs(void)
{
  return theNodes;
}

void
semiLoofPlate::Print(OPS_Stream &s, int flag)
{
  opserr << "semiLoofBeam::Print" << endln;
}

int
semiLoofPlate::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
  return 0;
}

//get the number of external nodes
//int
//semiLoofPlate::getNumExternalNodes() const
//{
//  return 9;
//}

//return connected external nodes
//const ID&
//semiLoofPlate::getExternalNodes()
//{
//  return *connectedNodes;
//}

int
semiLoofPlate::getNumDOF(void)
{
  return 35;  //????
}

//commit state
int
semiLoofPlate::commitState()
{
  int success = 0;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "semiLoofPlate::commitState () - failed in base class";
  }

  for (int i = 0; i < 9; i++)
    success += materialPointers[i]->commitState();

  return success;
}

//revert to last commit 
int
semiLoofPlate::revertToLastCommit()
{
  int i;
  int success = 0;

  for (i = 0; i < 9; i++)
    success += materialPointers[i]->revertToLastCommit();

  return success;
}

//revert to start 
int
semiLoofPlate::revertToStart()
{
  int i;
  int success = 0;

  for (i = 0; i < 9; i++)
    success += materialPointers[i]->revertToStart();

  return success;
}

int
semiLoofPlate::update()
{
  return 0;
}

//return secant matrix 
const Matrix&
semiLoofPlate::getInitialStiff()
{
  if (Ki != 0)
    return *Ki;

  static const int ndf = 6; //two membrane plus three bending plus one drill

  static const int nstress = 8; //three membrane, three moment, two shear

  static const int ngauss = 9;

  static const int numnodes = 9;

  int i, j, k, p, q;
  int jj, kk;
  int node;

  double volume = 0.0;

  static double xsj;  // determinant jacaobian matrix 

  static double dvol[ngauss]; //volume element

  static double shp[3][numnodes];  //shape functions at a gauss point

  static Matrix stiffJK(ndf, ndf); //nodeJK stiffness 

  static Matrix dd(nstress, nstress);  //material tangent

  //static Matrix J0(2,2) ;  //Jacobian at center

  //static Matrix J0inv(2,2) ; //inverse of Jacobian at center

  //---------B-matrices------------------------------------

  static Matrix BJ(nstress, ndf);      // B matrix node J

  static Matrix BJtran(ndf, nstress);

  static Matrix BK(nstress, ndf);      // B matrix node k

  static Matrix BJtranD(ndf, nstress);


  static Matrix Bbend(3, 3);  // bending B matrix

  static Matrix Bshear(2, 3); // shear B matrix

  static Matrix Bmembrane(3, 2); // membrane B matrix


  static double BdrillJ[ndf]; //drill B matrix

  static double BdrillK[ndf];

  double *drillPointer;

  static double saveB[nstress][ndf][numnodes];

  //-------------------------------------------------------

  //stiff.Zero();


  //compute Jacobian and inverse at center
  double L1 = 0.0;
  double L2 = 0.0;
  //computeJacobian( L1, L2, xl, J0, J0inv ) ; 

  //gauss loop 
  for (i = 0; i < ngauss; i++) {

    //get shape functions
//    shape2d(sg[i], tg[i], xl, shp, xsj);

    //volume element to also be saved
//    dvol[i] = wg[i] * xsj;

    volume += dvol[i];

    // j-node loop to compute strain 
    for (j = 0; j < numnodes; j++)  {

      //compute B matrix 

//      Bmembrane = computeBmembrane(j, shp);

//      Bbend = computeBbend(j, shp);

 //     Bshear = computeBshear(j, shp);

//      BJ = assembleB(Bmembrane, Bbend, Bshear);

      //save the B-matrix
      for (p = 0; p < nstress; p++) {
        for (q = 0; q < ndf; q++)
          saveB[p][q][j] = BJ(p, q);
      }//end for p

      //drilling B matrix
//      drillPointer = computeBdrill(j, shp);
      for (p = 0; p < ndf; p++) {
        //BdrillJ[p] = *drillPointer++ ;
        BdrillJ[p] = *drillPointer; //set p-th component
        drillPointer++;             //pointer arithmetic
      }//end for p
    } // end for j


    dd = materialPointers[i]->getInitialTangent();
    dd *= dvol[i];

    //residual and tangent calculations node loops

    jj = 0;
    for (j = 0; j < numnodes; j++) {

      //extract BJ
      for (p = 0; p < nstress; p++) {
        for (q = 0; q < ndf; q++)
          BJ(p, q) = saveB[p][q][j];
      }//end for p

      //multiply bending terms by (-1.0) for correct statement
      // of equilibrium  
      for (p = 3; p < 6; p++) {
        for (q = 3; q < 6; q++)
          BJ(p, q) *= (-1.0);
      } //end for p


      //transpose 
      //BJtran = transpose( 8, ndf, BJ ) ;
      for (p = 0; p < ndf; p++) {
        for (q = 0; q < nstress; q++)
          BJtran(p, q) = BJ(q, p);
      }//end for p

      //drilling B matrix
//      drillPointer = computeBdrill(j, shp);
      for (p = 0; p < ndf; p++) {
        BdrillJ[p] = *drillPointer;
        drillPointer++;
      }//end for p

      //BJtranD = BJtran * dd ;
      BJtranD.addMatrixProduct(0.0, BJtran, dd, 1.0);


      for (p = 0; p < ndf; p++)
//        BdrillJ[p] *= (Ktt*dvol[i]);


      kk = 0;
      for (k = 0; k < numnodes; k++) {

        //extract BK
        for (p = 0; p < nstress; p++) {
          for (q = 0; q < ndf; q++)
            BK(p, q) = saveB[p][q][k];
        }//end for p

        //drilling B matrix
 //       drillPointer = computeBdrill(k, shp);
        for (p = 0; p < ndf; p++) {
          BdrillK[p] = *drillPointer;
          drillPointer++;
        }//end for p

        //stiffJK = BJtranD * BK  ;
        // +  transpose( 1,ndf,BdrillJ ) * BdrillK ; 
        stiffJK.addMatrixProduct(0.0, BJtranD, BK, 1.0);

        for (p = 0; p < ndf; p++)  {
          for (q = 0; q < ndf; q++) {
//            stiff(jj + p, kk + q) += stiffJK(p, q)
//              + (BdrillJ[p] * BdrillK[q]);
          }//end for q
        }//end for p
        kk += ndf;
      } // end for k loop
      jj += ndf;
    } // end for j loop
  } //end for i gauss loop 
//  Ki = new Matrix(stiff);
//  return stiff;
}

const Matrix &
semiLoofPlate::getTangentStiff(void)
{
  // check for quick return
  //if (nen == 0)
  //  return (*semiLoofElementM[0]);

  // get the current load factor
  Domain *theDomain = this->getDomain();
  double dm = theDomain->getCurrentTime();

  // set ctan, ior and iow
  double ctan[3];
  ctan[0] = 1.0; ctan[1] = 0.0; ctan[2] = 0.0;
  int ior = 0; int iow = 0;

  // call the ready routine to set ul, xl, tl and ix, NH1, NH2 and NH3
  //int nstR = this->readyfRoutine(false);

  // zero the matrix
  //semiLoofElementM[nstR]->Zero();

  // invoke the Fortran subroutine
  int isw = 3;
  //int nstI = this->invokefRoutine(ior, iow, ctan, isw);

  // check nst is as determined in readyfRoutine()
  //if (nstI != nstR) {
  //  opserr << "FATAL semiLoofElement::getTangentStiff() problems with incompatable nst";
  //  opserr << " ready: " << nstR << " invoke: " << nstI << endln;
  //  exit(-1);
  //}

  // return the matrix

  //return *(semiLoofElementM[nstR]);
  Matrix temp(0, 0);
  return temp;
}

const Matrix &
semiLoofPlate::getDamp(void)
{
  // check for quick return
  //if (nen == 0)
  //  return (*semiLoofElementM[0]);

  // get the current load factor
  Domain *theDomain = this->getDomain();
  double dm = theDomain->getCurrentTime();

  // set ctan, ior and iow
  double ctan[3];
  ctan[0] = 0.0; ctan[1] = 1.0; ctan[2] = 0.0;
  int ior = 0; int iow = 0;

  // call the ready routine to set ul, xl, tl and ix, NH1, NH2 and NH3
  int NH1, NH2, NH3;
  //int nstR = this->readyfRoutine(true);

  // zero the matrix
  //semiLoofElementM[nstR]->Zero();

  // invoke the Fortran subroutine
  //int isw = 3; int nst = nen*ndf; int n = this->getTag();

  //int nstI = this->invokefRoutine(ior, iow, ctan, isw);

  // check nst is as determined in readyfRoutine()
//   if (nstI != nstR) {
//     opserr << "FATAL semiLoofElement::getTangentStiff() problems with incompatable nst";
//     opserr << " ready: " << nstR << " invoke: " << nstI << endln;
//     exit(-1);
//   }

  // return the matrix
  //return *(semiLoofElementM[nstR]);
  Matrix temp(0, 0);
  return temp;
}

const Matrix &
semiLoofPlate::getMass(void)
{
  // check for quick return
  //if (nen == 0)
  //  return (*semiLoofElementM[0]);

  // get the current load factor
  Domain *theDomain = this->getDomain();
  double dm = theDomain->getCurrentTime();

  // set ctan, ior and iow
  double ctan[3];
  ctan[0] = 0.0; ctan[1] = 0.0; ctan[2] = 1.0;
  int ior = 0; int iow = 0;

  // call the ready routine to set ul, xl, tl and ix, NH1, NH2 and NH3
  int NH1, NH2, NH3;
  //int nstR = this->readyfRoutine(true);

  // zero the matrix and vector (consistent and lumped)
  //semiLoofElementM[nstR]->Zero();
  //semiLoofElementV[nstR]->Zero();

  // invoke the Fortran subroutine
  //int isw = 5; int nst = nen*ndf; int n = this->getTag();
  //int nstI = this->invokefRoutine(ior, iow, ctan, isw);

  // check nst is as determined in readyfRoutine()
//   if (nstI != nstR) {
//     opserr << "FATAL semiLoofElement::getTangentStiff() problems with incompatable nst";
//     opserr << " ready: " << nstR << " invoke: " << nstI << endln;
//     exit(-1);
//   }

  // return the matrix
  //return *(semiLoofElementM[nstR]);
  Matrix temp(0, 0);
  return temp;
}

int
semiLoofPlate::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "semiLoofShell::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
  return -1;
}

//get residual with inertia terms
const Vector&
semiLoofPlate::getResistingForceIncInertia()
{
  static Vector res(54);
  int tang_flag = 0; //don't get the tangent

  //do tangent and residual here 
//   formResidAndTangent(tang_flag);

//   formInertiaTerms(tang_flag);

  res = resid;
  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    res += this->getRayleighDampingForces();

  // subtract external loads 
  if (theLoad != 0)
    res -= *theLoad;

  return res;
}

int
semiLoofPlate::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  // Now quad sends the ids of its materials
  int matDbTag;
  static ID idData(27);
  int i;

  for (i = 0; i < 9; i++) {
    idData(i) = materialPointers[i]->getClassTag();
    matDbTag = materialPointers[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0)
        materialPointers[i]->setDbTag(matDbTag);
    }
    idData(i + 9) = matDbTag;
  }

  idData(18) = this->getTag();
  idData(19) = (*connectedNodes)(0);
  idData(20) = (*connectedNodes)(1);
  idData(21) = (*connectedNodes)(2);
  idData(22) = (*connectedNodes)(3);
  idData(23) = (*connectedNodes)(4);
  idData(24) = (*connectedNodes)(5);
  idData(25) = (*connectedNodes)(6);
  idData(26) = (*connectedNodes)(7);
  idData(27) = (*connectedNodes)(8);
  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING semiLoofPlate::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  static Vector vectData(4);
  vectData(0) = alphaM;
  vectData(1) = betaK;
  vectData(2) = betaK0;
  vectData(3) = betaKc;

  res += theChannel.sendVector(dataTag, commitTag, vectData);
  if (res < 0) {
    opserr << "WARNING semiLoofPlate::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 9; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING semiLoofPlate::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  return res;
}

int
semiLoofPlate::recvSelf(int commitTag, Channel &theChannel,
FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dataTag = this->getDbTag();
  static ID idData(27);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING semiLoofPlate::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(18));
  (*connectedNodes)(0) = idData(19);
  (*connectedNodes)(1) = idData(20);
  (*connectedNodes)(2) = idData(21);
  (*connectedNodes)(3) = idData(22);
  (*connectedNodes)(4) = idData(23);
  (*connectedNodes)(5) = idData(24);
  (*connectedNodes)(6) = idData(25);
  (*connectedNodes)(7) = idData(26);
  (*connectedNodes)(8) = idData(27);
  static Vector vectData(4);
  res += theChannel.recvVector(dataTag, commitTag, vectData);
  if (res < 0) {
    opserr << "WARNING ShellNL::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  alphaM = vectData(0);
  betaK = vectData(1);
  betaK0 = vectData(2);
  betaKc = vectData(3);

  int i;

  if (materialPointers[0] == 0) {
    for (i = 0; i < 9; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i + 9);
      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewSection(matClassTag);
      if (materialPointers[i] == 0) {
        opserr << "semiLoofPlate::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;;
        return -1;
      }
      // Now receive materials into the newly allocated space
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "semiLoofPlate::recvSelf() - material " << i << "failed to recv itself\n";
        return res;
      }
    }
  }
  // Number of materials is the same, receive materials into current space
  else {
    for (i = 0; i < 9; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i + 9);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (materialPointers[i]->getClassTag() != matClassTag) {
        delete materialPointers[i];
        materialPointers[i] = theBroker.getNewSection(matClassTag);
        if (materialPointers[i] == 0) {
          opserr << "semiLoofPlate::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;
          exit(-1);
        }
      }
      // Receive the material
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "semiLoofPlate::recvSelf() - material " << i << "failed to recv itself\n";
        return res;
      }
    }
  }
  return res;
}

void
semiLoofPlate::setDomain(Domain *theDomain)
{
  opserr << "semiLoofPlate::setDomain" << endln;
}

void
semiLoofPlate::zeroLoad()
{
  if (theLoad != 0)
    theLoad->Zero();

  return;
}

int
semiLoofPlate::addInertiaLoadToUnbalance(const Vector &accel)
{
  opserr << "semiLoofPlate::addInertiaLoadToUnbalance" << endln;
  return 0;
}

const Vector &
semiLoofPlate::getResistingForce()
{
  opserr << "semiLoofPlate::addInertiaLoadToUnbalance" << endln;
  return 0;
}

