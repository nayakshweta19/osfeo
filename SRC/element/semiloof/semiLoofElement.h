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

// $Revision: 1.4 $
// $Date: 2003/02/14 23:01:08 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/femski/semiLoofElement.h,v $


#ifndef semiLoofElement_h
#define semiLoofElement_h

//
// Description: This file contains the class definition for semiLoofElement. semiLoofElement
// is a wrapper used to call Fortran element subroutine. It is an abstract class.
//
// What: "@(#) semiLoofElement.h, revA"

#include <Element.h>
#include <Matrix.h>

class Node;
class Channel;
class UniaxialMaterial;

class semiLoofElement : public Element
{
public:
  // constructors
  semiLoofElement(int tag,
    int classTag,
    int eleType,
    int sizeD, int nen,
    int ndm, int ndf,
    int nh1, int nh3);

//  semiLoofElement(int tag,
//    int classTag,
//    int eleType,
//    int sizeD, int nen,
//    int ndm, int ndf, int iow);

  semiLoofElement(int classTag);

  // destructor
  virtual ~semiLoofElement();

  // public methods for element operations
  virtual int getNumExternalNodes(void) const;
  virtual const ID &getExternalNodes(void);
//  Node **getNodePtrs(void);

  virtual int getNumDOF(void);
//  virtual void setDomain(Domain *theDomain);

  virtual int commitState(void);
  virtual int revertToLastCommit(void);
  virtual int revertToStart(void);
//  virtual int update(void);

//  virtual const Matrix &getTangentStiff(void);
//  virtual const Matrix &getInitialStiff(void);
//  virtual const Matrix &getDamp(void);
//  virtual const Matrix &getMass(void);

  virtual void zeroLoad(void){ return; };
  virtual int addLoad(ElementalLoad *theLoad, double loadFactor);
  virtual int addInertiaLoadToUnbalance(const Vector &accel){ return 0; };

  virtual const Vector &getResistingForce(void) { return Vector(0); };
  virtual const Vector &getResistingForceIncInertia(void);

  // public methods for output
  virtual int sendSelf(int commitTag, Channel &theChannel);
  virtual int recvSelf(int commitTag, Channel &theChannel,
    FEM_ObjectBroker &theBroker);
//  virtual int displaySelf(Renderer &theViewer, int displayMode, float fact);
//  virtual void Print(OPS_Stream &s, int flag = 0);

protected:
  // protected methods 
  //virtual int invokefRoutine(int ior, int iow, double *ctan, int isw);
  //virtual int readyfRoutine(bool incInertia);
  //virtual int invokefInit(int isw, int iow);

  // protected data
  Vector *data;
  ID *connectedNodes;

private:
  // private attributes - a copy for each object of the class      
  double *h;            // the elements local h array
  double *u;            // to hold u^(k-1)_(n+1) as nodes don't save this
  double *d; 	          // data stored for each element
  int eleType;          // defines which elmtnn subroutine to invoke
  int ndf;		  // number of dof at nodes  nst = nen x ndf
  int nen;              // number of element nodes
  int ndm;	          // dimension of mesh
  int nh1, nh3;         // the size of nh1(nh2) and nh3 data in h
  int nrCount;          // needed for Prof. Fillipou's Elmt05

  // static data - single copy for all objects of the class	
  static Matrix **semiLoofElementM;   // class wide matrices - use s array
  static Vector **semiLoofElementV;   // class wide vectors - use r array
  static double *s;  // element tangent (nst x nst)
  static double *r;  // element residual (ndf x nen) == (nst)
  static double *ul; // nodal responses (ndf x nen x 5) == (nst X 5)
  static double *xl; // nodal coordinates (ndf x nen) == (nst)
  static double *tl; // nodal temp (nen)
  static int    *ix; // nodal tags (nen)
  static int numSemiloofElements;
};

#endif




