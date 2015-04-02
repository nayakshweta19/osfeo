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

// $Revision: 1.7 $
// $Date: 2003/04/02 22:02:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/femski/semiLoofElement.cpp,v $
//
// Description: This file contains the implementation for the semiLoofElement class.
//

#include "semiLoofElement.h"
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>

// initialize all class wise pointers to 0 and numSemiloofElements to 0
Matrix **semiLoofElement::semiLoofElementM;
Vector **semiLoofElement::semiLoofElementV;
double *semiLoofElement::r;
double *semiLoofElement::s;
double *semiLoofElement::ul;
double *semiLoofElement::xl;
double *semiLoofElement::tl;
int    *semiLoofElement::ix;
int    semiLoofElement::numSemiloofElements(0);

static double *work = 0;
static int sizeWork = 0;

#define MAX_NST 64

// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the semiLoofElement end nodes.
semiLoofElement::semiLoofElement(int tag,
  int classTag,
  int EleType,
  int sizeD, int NEN,
  int NDM, int NDF,
  int numNh1, int numNh3)
  :Element(tag, classTag), nh1(numNh1), nh3(numNh3), h(0), eleType(EleType),
  u(0), nen(NEN), ndf(NDF), ndm(NDM), d(0), data(0),
  connectedNodes(0), nrCount(0)
{
  // allocate space for h array
  if (nh1 < 0) nh1 = 0;
  if (nh3 < 0) nh3 = 0;
  if (nh1 != 0 || nh3 != 0) {
    int sizeH = 2 * nh1 + nh3;
    h = new double[sizeH];
    if (sizeWork < sizeH) {
      if (work != 0)
        delete[] work;
      work = new double[sizeH];
      if (work == 0) {
        opserr << "FATAL: semiLoofElement::semiLoofElement() - eleTag: " << tag;
        opserr << " ran out of memory creating h of size " << 2 * nh1 + nh3 << endln;
        exit(-1);
      }
      sizeWork = sizeH;
    }
    if (h == 0 || work == 0) {
      opserr << "FATAL: semiLoofElement::semiLoofElement() - eleTag: " << tag;
      opserr << " ran out of memory creating h of size " << 2 * nh1 + nh3 << endln;
      exit(-1);
    }

    for (int i = 0; i < sizeH; i++)
      h[i] = 0.0;
  }

  connectedNodes = new ID(NEN);
  d = new double[sizeD];
  for (int i = 0; i < sizeD; i++) d[i] = 0.0;
  data = new Vector(d, sizeD);

  // allocate space for static variables on creation of first instance
//  if (numSemiloofElements == 0) {
//    semiLoofElementM = new Matrix *[MAX_NST + 1];
//    semiLoofElementV = new Vector *[MAX_NST + 1];
//    s = new double[(MAX_NST + 1)*(MAX_NST + 1)];
//    r = new double[MAX_NST + 1];
//    ul = new double[(MAX_NST + 1) * 6];
//    xl = new double[MAX_NST + 1];
//    tl = new double[MAX_NST + 1];
//    ix = new int[MAX_NST + 1];
//
//    // check space was available -- otherwise exit
//    if (semiLoofElementM == 0 || semiLoofElementV == 0 || ix == 0 ||
//      r == 0 || s == 0 || ul == 0 || xl == 0 || tl == 0) {
//
//      opserr << "FATAL: semiLoofElement::semiLoofElement() - eleTag: " << tag;
//      opserr << " ran out of memory initializing static stuff\n";
//      exit(-1);
//    }
//
//    for (int i = 0; i < MAX_NST + 1; i++) {
//      semiLoofElementM[i] = 0;
//      semiLoofElementV[i] = 0;
//    }
//    semiLoofElementM[0] = new Matrix(1, 1); // dummy for error
//    semiLoofElementV[0] = new Vector(1);
//  }

  // increment number of elements
  numSemiloofElements++;
}

//semiLoofElement::semiLoofElement(int tag,
//  int classTag,
//  int EleType,
//  int sizeD, int NEN,
//  int NDM, int NDF, int iow)
//  :Element(tag, classTag), nh1(0), nh3(0), h(0), eleType(EleType),
//  u(0), nen(NEN), ndf(NDF), ndm(NDM), d(0), data(0),
//  connectedNodes(0), nrCount(0)
//{
//  connectedNodes = new ID(NEN);
//  d = new double[sizeD];
//  data = new Vector(d, sizeD);
//  if (d == 0 || data == 0) {
//    opserr << "FATAL: semiLoofElement::semiLoofElement() - eleTag: " << tag;
//    opserr << " ran out of memory creating d of size " << sizeD << endln;
//    exit(-1);
//  }
//  for (int i = 0; i < sizeD; i++) d[i] = 0.0;
//
//  // invoke the elmt() routine with isw == 1 to read in the element data
//  // and create room for the h array stuff needed by the element
////   this->invokefInit(1, iow);
//
//  // allocate space for h array
//  if (nh1 < 0) nh1 = 0;
//  if (nh3 < 0) nh3 = 0;
//  if (nh1 != 0 || nh3 != 0) {
//    int sizeH = 2 * nh1 + nh3;
//    h = new double[sizeH];
//    if (sizeWork < sizeH) {
//      if (work != 0)
//        delete[] work;
//      work = new double[sizeH];
//      if (work == 0) {
//        opserr << "FATAL: semiLoofElement::semiLoofElement() - eleTag: " << tag;
//        opserr << " ran out of memory creating h of size " << 2 * nh1 + nh3 << endln;
//        exit(-1);
//      }
//      sizeWork = sizeH;
//    }
//    if (h == 0 || work == 0) {
//      opserr << "FATAL: semiLoofElement::semiLoofElement() - eleTag: " << this->getTag();
//      opserr << " ran out of memory creating h of size " << sizeH << endln;
//      exit(-1);
//    }
//    else
//      for (int i = 0; i < sizeH; i++) h[i] = 0.0;
//  }
//
//  // allocate space for static variables on creation of first instance
//  if (numSemiloofElements == 0) {
//    semiLoofElementM = new Matrix *[MAX_NST + 1];
//    semiLoofElementV = new Vector *[MAX_NST + 1];
//    s = new double[(MAX_NST + 1)*(MAX_NST + 1)];
//    r = new double[MAX_NST + 1];
//    ul = new double[(MAX_NST + 1) * 6];
//    xl = new double[MAX_NST + 1];
//    tl = new double[MAX_NST + 1];
//    ix = new int[MAX_NST + 1];
//
//    // check space was available -- otherwise exit
//    if (semiLoofElementM == 0 || semiLoofElementV == 0 || ix == 0 ||
//      r == 0 || s == 0 || ul == 0 || xl == 0 || tl == 0) {
//
//      opserr << "FATAL: semiLoofElement::semiLoofElement() - eleTag: " << tag;
//      opserr << " ran out of memory initializing static stuff\n";
//      exit(-1);
//    }
//
//    for (int i = 0; i < MAX_NST + 1; i++) {
//      semiLoofElementM[i] = 0;
//      semiLoofElementV[i] = 0;
//    }
//    semiLoofElementM[0] = new Matrix(1, 1); // dummy for error
//    semiLoofElementV[0] = new Vector(1);
//  }
//
//  // increment number of elements
//  numSemiloofElements++;
//}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
semiLoofElement::semiLoofElement(int classTag)
  :Element(0, classTag), nh1(0), nh3(0), h(0),
  u(0), nen(0), ndf(0), ndm(0), d(0), data(0), connectedNodes(0)
{
  // does nothing
}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the material object.
semiLoofElement::~semiLoofElement()
{
  // clear up any space allocated for individual object

  if (h != 0)
    delete[] h;
  if (u != 0)
    delete[] u;
  if (data != 0)
    delete  data;
  if (connectedNodes != 0)
    delete  connectedNodes;
  if (d != 0)
    delete[] d;

  // if last element - clear up space allocated

  numSemiloofElements--;
  if (numSemiloofElements == 0) {
    for (int i = 0; i < MAX_NST + 1; i++) {
      if (semiLoofElementM[i] != 0) delete semiLoofElementM[i];
      if (semiLoofElementV[i] != 0) delete semiLoofElementV[i];
    }
    delete[] semiLoofElementM;
    delete[] semiLoofElementV;
    delete[] s;
    delete[] r;
    delete[] ul;
    delete[] xl;
    delete[] tl;
    delete[] ix;
  }
}

int
semiLoofElement::getNumExternalNodes(void) const
{
  return connectedNodes->Size();
}

const ID &
semiLoofElement::getExternalNodes(void)
{
  return *connectedNodes;
}

int
semiLoofElement::getNumDOF(void)
{
  return ndf*nen;
}

int
semiLoofElement::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "semiLoofElement::commitState () - failed in base class";
  }

  if (nh1 != 0)
    for (int i = 0; i < nh1; i++)
      h[i] = h[i + nh1];

  nrCount = 0;
  return retVal;
}

int
semiLoofElement::revertToLastCommit()
{
  if (nh1 != 0)
    for (int i = 0; i < nh1; i++)
      h[i + nh1] = h[i];

  nrCount = 0;
  return 0;
}

int
semiLoofElement::revertToStart()
{
  // does nothing
  return 0;
}

int
semiLoofElement::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "semiLoofElement::addLoad - load type unknown for truss with tag: " << this->getTag() << endln;
  return -1;
}

const Vector &
semiLoofElement::getResistingForceIncInertia()
{
  opserr << "semiLoofElement::getResistingForceIncInertia() - subclass responds this\n";
  return -1;
}

int
semiLoofElement::sendSelf(int commitTag, Channel &theChannel)
{
  opserr << "semiLoofElement::sendSelf() - not yet implemented\n";
  return -1;
}

int
semiLoofElement::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  opserr << "semiLoofElement::recvSelf() - not yet implemented\n";
  return -1;
}

//int
//semiLoofElement::displaySelf(Renderer &theViewer, int displayMode, float fact)
//{
//  return 0;
//}
//
//void
//semiLoofElement::Print(OPS_Stream &s, int flag)
//{
//  int ior = 0; int iow = 1;
//  ior = -1;
//
//  s << "semiLoofElement::Print() - can only print to cerr at present\n";
//  ior = -1;
//
//  // get the current load factor
//  Domain *theDomain = this->getDomain();
//  double dm = theDomain->getCurrentTime();
//
//  // set ctan, ior and iow
//  double ctan[3];
//  ctan[0] = 0.0; ctan[1] = 0.0; ctan[2] = 0.0;
//
//
//  // call the ready routine to set ul, xl, tl and ix, NH1, NH2 and NH3
//  int NH1, NH2, NH3;
//  int nstR = this->readyfRoutine(false);
//
//  // invoke the Fortran subroutine
//  int isw = 4; int nst = nen*ndf; int n = this->getTag();
//  int nstI = this->invokefRoutine(ior, iow, ctan, isw);
//}

#ifdef _WIN32

extern "C" int GETCOMMON(int *mynh1, int *mynh3, int *sizeH,
  double *myh);


extern "C" int FILLCOMMON(int *mynen, double *mydm, int *myn,
  int *myior, int *myiow, int *mynh1,
  int *mynh2, int *mynh3, int *sumnh,
  double *myh, double *myctan,
  int *nrCount);

extern "C" int ELMT01(double *d, double *ul, double *xl, int *ix,
  double *tl, double *s, double *r, int *ndf,
  int *ndm, int *nst, int *isw);

extern "C" int ELMT02(double *d, double *ul, double *xl, int *ix,
  double *tl, double *s, double *r, int *ndf,
  int *ndm, int *nst, int *isw);

extern "C" int ELMT03(double *d, double *ul, double *xl, int *ix,
  double *tl, double *s, double *r, int *ndf,
  int *ndm, int *nst, int *isw);

extern "C" int ELMT04(double *d, double *ul, double *xl, int *ix,
  double *tl, double *s, double *r, int *ndf,
  int *ndm, int *nst, int *isw);

extern "C" int ELMT05(double *d, double *ul, double *xl, int *ix,
  double *tl, double *s, double *r, int *ndf,
  int *ndm, int *nst, int *isw);

#define getcommon_ 	GETCOMMON
#define fillcommon_	FILLCOMMON
#define elmt01_		ELMT01
#define elmt02_		ELMT02
#define elmt03_		ELMT03
#define elmt04_		ELMT04
#define elmt05_		ELMT05

#else
extern "C" int getcommon_(int *mynh1, int *mynh3, int *sizeH, double *myh);


extern "C" int fillcommon_(int *mynen, double *mydm, int *myn, int *myior,
  int *myiow, int *mynh1, int *mynh2, int *mynh3,
  int *sumnh, double *myh, double *myctan, int *nrCount);

extern "C" int elmt01_(double *d, double *ul, double *xl, int *ix, double *tl,
  double *s, double *r, int *ndf, int *ndm, int *nst, int *isw);

extern "C" int elmt02_(double *d, double *ul, double *xl, int *ix, double *tl,
  double *s, double *r, int *ndf, int *ndm, int *nst, int *isw);

extern "C" int elmt03_(double *d, double *ul, double *xl, int *ix, double *tl,
  double *s, double *r, int *ndf, int *ndm, int *nst, int *isw);

extern "C" int elmt04_(double *d, double *ul, double *xl, int *ix, double *tl,
  double *s, double *r, int *ndf, int *ndm, int *nst, int *isw);

extern "C" int elmt05_(double *d, double *ul, double *xl, int *ix, double *tl,
  double *s, double *r, int *ndf, int *ndm, int *nst, int *isw);
#endif	       

/*int
semiLoofElement::invokefRoutine(int ior, int iow, double *ctan, int isw)
{
  // fill the common blocks
  // determine position in h of nh1, nh2 and nh3 - remember Fortran indexing
  int NH1, NH2, NH3;
  if (nh1 != 0) {
    NH1 = 1;
    NH2 = nh1 + NH1;
    NH3 = nh1 + NH2;
  }
  else {
    NH1 = 1;
    NH2 = 1;
    NH3 = 1;
  }

  int NDM = ndm;
  int NDF = ndf;

  int n = this->getTag();
  int sum = 2 * nh1 + nh3;
  int count = nrCount;

  double dm = 0.0; // load factor

  fillcommon_(&nen, &dm, &n, &ior, &iow, &NH1, &NH2, &NH3, &sum,
    h, ctan, &count);

  // invoke the Fortran subroutine

  int nst = nen*ndf;
  if (nst != 0) {
    if (eleType == 1)
      elmt01_(d, ul, xl, ix, tl, s, r, &NDF, &NDM, &nst, &isw);
    else if (eleType == 2)
      elmt02_(d, ul, xl, ix, tl, s, r, &NDF, &NDM, &nst, &isw);
    else if (eleType == 3)
      elmt03_(d, ul, xl, ix, tl, s, r, &NDF, &NDM, &nst, &isw);
    else if (eleType == 4)
      elmt04_(d, ul, xl, ix, tl, s, r, &NDF, &NDM, &nst, &isw);
    else if (eleType == 5)
      elmt05_(d, ul, xl, ix, tl, s, r, &ndf, &NDM, &nst, &isw);
    else {
      opserr << "semiLoofElement::invokefRoutine() unknown element type ";
      opserr << eleType << endln;
    }

    // now copy the stuff from common block to h array
    getcommon_(&NH1, &NH3, &sum, h);
  }

  return nst;
}


int
semiLoofElement::invokefInit(int isw, int iow)
{
  // fill the common blocks
  // determine position in h of nh1, nh2 and nh3 - remember Fortarn indexing
  int NH1 = 0;
  int NH2 = 0;
  int NH3 = 0;

  int NDM = ndm;
  int NDF = ndf;
  double ctan[3];

  int n = this->getTag();
  int sum = 0;
  int ior = 0;
  int count = nrCount;

  double dm = 0.0;

  fillcommon_(&nen, &dm, &n, &ior, &iow, &NH1, &NH2, &NH3, &sum,
    h, ctan, &count);

  // invoke the Fortran subroutine

  int nst = nen*ndf;
  if (nst != 0) {
    if (eleType == 1)
      elmt01_(d, ul, xl, ix, tl, s, r, &NDF, &NDM, &nst, &isw);
    else if (eleType == 2)
      elmt02_(d, ul, xl, ix, tl, s, r, &NDF, &NDM, &nst, &isw);
    else if (eleType == 3)
      elmt03_(d, ul, xl, ix, tl, s, r, &NDF, &NDM, &nst, &isw);
    else if (eleType == 4)
      elmt04_(d, ul, xl, ix, tl, s, r, &NDF, &NDM, &nst, &isw);
    else if (eleType == 5)
      elmt05_(d, ul, xl, ix, tl, s, r, &NDF, &NDM, &nst, &isw);
    else {
      opserr << "semiLoofElement::invokefRoutine() unknown element type ";
      opserr << eleType << endln;
    }

    if (nst < 0) {
      opserr << "FATAL: semiLoofElement::semiLoofElement() - eleTag: " << this->getTag();
      opserr << " ran out of memory creating h of size " << nst << endln;
      exit(-1);
    }
  }

  // now get the size of the state info needed by the element
  sum = 0;
  getcommon_(&NH1, &NH3, &sum, h);
  nh1 = NH1; nh3 = NH3;
  return 0;
}*/

//int
//semiLoofElement::update()
//{
//  // determine nst 
//  int nst = ndf*nen;
//
//  // increment the newton-raphson count -- needed for Prof. Fillipou's element
//  nrCount++;
//
//  return 0;
//}

