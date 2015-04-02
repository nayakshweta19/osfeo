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
// $Source: /usr/local/cvs/OpenSees/SRC/element/femski/semiLoofShell.h,v $


#ifndef semiLoofShell_h
#define semiLoofShell_h
//
// Description: This file contains the class definition for semiLoofShell. semiLoofShell
// is a wrapper used to call Fortran element subroutine SHELL. SHELL is a 
// semi-loof shell 3d element
//

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

class semiLoofShell : public semiLoofElement
{
public:
  // constructors
  semiLoofShell(int tag,
    int nd1, int nd2, int nd3, int nd4, int nd5, int nd6, int nd7, int nd8, int nd9,
    SectionForceDeformation &theMaterial, double rho = 0.0);

  //semiLoofShell(int tag,
  //  int nd1, int nd2, int nd3, int nd4, int nd5, int nd6, int nd7, int nd8, int nd9,
  //  SectionForceDeformation &theMaterial);

  semiLoofShell();

  // destructor
  ~semiLoofShell();

  //get the number of external nodes
//  int getNumExternalNodes() const;
  //return connected external nodes
//  const ID &getExternalNodes();
  Node **getNodePtrs();
  //return number of dofs
  int getNumDOF();

  //set domain 
  void setDomain(Domain *theDomain);

  int commitState();
  int revertToLastCommit();
  int revertToStart();
  int update(void);

  //print out element data
  void Print(OPS_Stream &s, int flag);

  //return stiffness matrix 
  const Matrix &getTangentStiff();
  const Matrix &getInitialStiff();
  const Matrix &getMass();
  const Matrix &getDamp(void);
  // methods for applying loads
  void zeroLoad(void);
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);

  //get residual
  const Vector &getResistingForce();

  //get residual with inertia terms
  const Vector &getResistingForceIncInertia();

  // public methods for element output
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker
    &theBroker);


  Response* setResponse(const char **argv, int argc, OPS_Stream &output);
  int getResponse(int responseID, Information &eleInfo);

  //plotting 
  int displaySelf(Renderer &theViewer, int displayMode, float fact);

protected:

private:

  Node *theNodes[9];       // pointer to the elements nodes
  Vector *theLoad;    // vector to hold the applied load P
  Matrix *Ki;

  //static data
  static Matrix stiff;
  static Vector resid;
  static Matrix mass;
  static Matrix damping;

  //quadrature data
  static const double root3;
  static const double root3_over_root5;
  static double sg[9];
  static double tg[9];
  static double wg[9];

  //node information
  //ID connectedExternalNodes;  //nine node numbers
  //pointers to nine nodes
  //Node *nodePointers[9];

  //drilling stiffness
  double Ktt;

  //material information: pointers to four materials
  SectionForceDeformation *materialPointers[9];

  //local nodal coordinates, two coordinates for each of nine nodes
  //static double xl[][9] ; 
  double xl[2][9];

  //shell basis vectors
  double g1[3];
  double g2[3];
  double g3[3];

  //compute local coordinates and basis
  void computeBasis();

  //inertia terms
  void formInertiaTerms(int tangFlag);

  //form residual and tangent					  
  void formResidAndTangent(int tang_flag);

  //compute Jacobian matrix and inverse at point {L1,L2}
  //void  computeJacobian( double L1, double L2,const double x[2][9], 
  //                       Matrix &JJ,Matrix &JJinv ) ;

  //compute Bdrill matrix
  double* computeBdrill(int node, const double shp[3][9]);

  //assemble a B matrix 
  const Matrix& assembleB(const Matrix &Bmembrane,
    const Matrix &Bbend,
    const Matrix &Bshear);

  //compute Bmembrane matrix
  const Matrix& computeBmembrane(int node, const double shp[3][9]);

  //compute Bbend matrix
  const Matrix& computeBbend(int node, const double shp[3][9]);

  //compute standard Bshear matrix
  const Matrix&  computeBshear(int node, const double shp[3][9]);

  //Matrix transpose
  Matrix transpose(int dim1, int dim2, const Matrix &M);

  //shape function routine for four node quads
  void shape2d(double ss, double tt,
    const double x[2][9],
    double shp[3][9],
    double &xsj);

};

#endif



