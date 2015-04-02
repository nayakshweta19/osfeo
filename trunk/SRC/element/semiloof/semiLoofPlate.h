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
// $Source: /usr/local/cvs/OpenSees/SRC/element/femski/semiLoofPlate.h,v $


#ifndef semiLoofPlate_h
#define semiLoofPlate_h

// File: ~/element/fortran/semiLoofPlate.h

// Description: This file contains the class definition for semiLoofPlate. semiLoofPlate
// is a wrapper used to call Fortran element subroutine semiLoofPlate. semiLoofPlate is a 
// semi-loof plate element


#include <semiLoofElement.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <R3vectors.h>

class semiLoofPlate : public semiLoofElement
{
public:
  // constructors
  semiLoofPlate(int tag,
    int Nd1, int Nd2, int Nd3, int Nd4, int Nd5, int Nd6, int Nd7, int Nd8, int Nd9,
    SectionForceDeformation &theMaterial, double rho = 0.0);

  semiLoofPlate();

  // destructor
  ~semiLoofPlate();

  // public methods for element operations
//  int getNumExternalNodes(void) const;
//  const ID &getExternalNodes(void);
  Node **getNodePtrs(void);

  int getNumDOF(void);
  void setDomain(Domain *theDomain);

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  int update(void);

  const Matrix &getTangentStiff(void);
  const Matrix &getInitialStiff(void);
  const Matrix &getDamp(void);
  const Matrix &getMass(void);

  void zeroLoad(void);
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);

  const Vector &getResistingForce(void);
  const Vector &getResistingForceIncInertia(void);

  // public methods for output
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel,
    FEM_ObjectBroker &theBroker);
  int displaySelf(Renderer &theViewer, int displayMode, float fact);
  void Print(OPS_Stream &s, int flag = 0);

protected:
  // protected methods 
  //int invokefRoutine(int ior, int iow, double *ctan, int isw);
  //int readyfRoutine(bool incInertia);
  //int invokefInit(int isw, int iow);

private:
  Node *theNodes[9];       // pointer to the elements nodes
  Vector *theLoad;    // vector to hold the applied load P
  Matrix *Ki;

  //static data
  static Matrix stiff;
  static Vector resid;
  static Matrix mass;
  static Matrix damping;

  //material information: pointers to four materials
  SectionForceDeformation *materialPointers[9];
};

#endif



