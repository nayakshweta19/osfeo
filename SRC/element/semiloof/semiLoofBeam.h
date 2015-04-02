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
// $Source: /usr/local/cvs/OpenSees/SRC/element/femski/semiLoofBeam.h,v $


#ifndef semiLoofBeam_h
#define semiLoofBeam_h

// Description: This file contains the class definition for semiLoofBeam. semiLoofBeam
// is a wrapper used to call Fortran element subroutine elmt02. elmt02 is a 
// linear elastic 2d element
//
// What: "@(#) semiLoofBeam.h, revA"

#include <Element.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <CrdTransf.h>
#include <SectionForceDeformation.h>
#include <R3vectors.h>

class semiLoofBeam : public Element
{
public:
  // constructors
  semiLoofBeam(int tag,
    int Nd1, int Nd2, int Nd3,
    SectionForceDeformation **s,
    CrdTransf &coordTransf,
    double rho = 0.0);

  semiLoofBeam();

  // destructor
  ~semiLoofBeam();

  // public methods for element operations
  int getNumExternalNodes(void) const;
  const ID &getExternalNodes(void);
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
  Matrix getNL(int sec, double L);

protected:
  // protected methods 
  //int invokefRoutine(int ior, int iow, double *ctan, int isw);
  //int readyfRoutine(bool incInertia);
  //int invokefInit(int isw, int iow);

private:
  Node *theNodes[3];       // pointer to the elements nodes
  SectionForceDeformation *theSections[3];  //material information: pointers to four materials
  CrdTransf *crdTransf;
  ID connectedExternalNodes;
  
  Vector *theLoad;    // vector to hold the applied load P
  Matrix *Ki;
  double rho;    // Mass density per unit length
  int parameterID;

  Vector Q;		// Applied nodal loads
  Vector q;		// Basic force
  double q0[5];   // Fixed end forces in basic system (no torsion)
  double p0[5];   // Reactions in basic system (no torsion)

  //static data
  static Matrix stiff;
  static Vector resid;
  static Matrix mass;
  static Matrix damping;

//  static double workArea[];
  static double *wbeam;
  static double *frlof; //Matrix of local unit vector at loof nods
  int lnods;    //Matrix of global number DOF the element, 3nodes
  static double *elxyz; //Matrix of nodal coordinate
  double side;          //Norm of vector X, before normalization(used for numerical integration)
  int noelz;            //Max number of lnodz for act structure
  int lvabz;            //total number of dof of the actural element
  int ielem;
  int nelem;            //total num of elements
  int nfirst;
  int lnomax;
  static double *point; //global coordinates of the sampling point
  double ksi;
  int kount;           //IGAUS Number of the actual integrating point
  static double *atnod;//auxiliary point used to define local X,Y reference plane
  double excen;        //Norm of the eccentricity vector of DEXCE
  static double *dexce;//unit vector defining the eccentricity direction for beam element
  static double *shear;//constraint matrix
};

#endif



