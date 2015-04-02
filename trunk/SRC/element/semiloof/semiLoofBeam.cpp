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
// $Source: /usr/local/cvs/OpenSees/SRC/element/femski/semiLoofBeam.cpp,v $

// Description: This file contains the implementation for the semiLoofBeam class.
// semi-loof beam

#include "semiLoofBeam.h"
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <CrdTransf.h>
#include <LobattoBeamIntegration.h>
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

static int numSemiLoofBeam = 0;

double *semiLoofBeam::wbeam;
double *semiLoofBeam::frlof;
double *semiLoofBeam::elxyz;
double *semiLoofBeam::point;
double *semiLoofBeam::atnod;
double *semiLoofBeam::dexce;
double *semiLoofBeam::shear;

void *
OPS_NewSemiLoofBeam(void)
{
  if (numSemiLoofBeam == 0) {
    //    opserr << "Using SemiLoofBeam - Developed by: Ning Li, neallee@tju.edu.cn\n";
    numSemiLoofBeam++;
  }

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 6) {
    opserr << "Want: element SemiLoofBeam eleTag? node1? node2? node3? secTag? transfTag?";
    return 0;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element SemiLoofBeam \n";
    return 0; //return OPS_Error;
  }

  SectionForceDeformation *theSection = OPS_GetSectionForceDeformation(iData[4]);
  if (theSection == 0) {
    opserr << "ERROR:  element SemiLoofBeam " << iData[0] << "section " << iData[4] << " not found\n";
    return 0;
  }

  SectionForceDeformation **sections = new SectionForceDeformation*[3];

  if (!sections) {
    opserr << "WARNING OPS_NewSemiLoofBeam - Insufficient memory to create sections\n";
    return 0;
  }

  for (int j = 0; j<3; j++) {
    sections[j] = theSection;
  }
  CrdTransf *crdTransf = OPS_GetCrdTransf(iData[5]);
  if (crdTransf == 0) {
    opserr << "ERROR: transformation not found\n";
    opserr << "transformation: " << iData[5];
    opserr << " element: " << iData[0] << endln;
    return 0;
  }

  double rho = 0.0;
  if (numArgs > 6) {
    double dData[1];
    numData = 1;
    if (OPS_GetDouble(&numData, dData) != 0) {
      opserr << "WARNING invalid rho parameter of element SemiLoofBeam\n";
      return 0;
    }
    rho = dData[0];
  }

  Element *theElement = new semiLoofBeam(iData[0], iData[1], iData[2], iData[3],
    sections, *crdTransf, rho);

  return theElement;
}

//static data
Matrix  semiLoofBeam::stiff(17, 17); //6(A)+5(B with 2 loof)+6(C)
Vector  semiLoofBeam::resid(17);
Matrix  semiLoofBeam::mass(17, 17);

semiLoofBeam::semiLoofBeam(int tag, 
  int nd1, int nd2, int nd3, 
  SectionForceDeformation **s, 
  CrdTransf &coordTransf, 
  double rho)
  :Element(tag, ELE_TAG_semiLoofBeam),
  connectedExternalNodes(3), crdTransf(0), Q(18), q(17),
  theLoad(0), Ki(0), rho(rho), parameterID(0)
{
  connectedExternalNodes(0) = nd1;
  connectedExternalNodes(1) = nd2;
  connectedExternalNodes(2) = nd3;

  crdTransf = coordTransf.getCopy3d();
  if (crdTransf == 0) {
    opserr << "semiLoofBeam::semiLoofBeam - failed to copy coordinate transformation\n";
    exit(-1);
  }

  if (theSections == 0) {
    opserr << "semiLoofBeam::semiLoofBeam - failed to allocate section model pointer\n";
    exit(-1);
  }
  for (int i = 0; i< 3; i++){
    // Get copies of the material model for each integration point
    theSections[i] = s[i]->getCopy();

    // Check allocation
    if (theSections[i] == 0) {
      opserr << "semiLoofBeam::semiLoofBeam() --failed to get a copy of section model " << endln;
      exit(-1);
    }
  }
  theNodes[0] = 0;
  theNodes[1] = 0;
  theNodes[2] = 0;

  frlof = new double[3 * 6];
  elxyz = new double[9 * 4]; // 3*3
  side = 0.0;
  noelz = 5; //3
  lnomax = 5;
  lvabz = 17; //21;
  ielem = 1;
  nelem = 1;
  nfirst = 1;
  lnods = 3;//  lnods = new int[ielem * nelem];
  point = new double[3];
  kount = 1;
  atnod = new double[3];
  excen = 0;
  dexce = new double[3];
  shear = new double[4 * 21];
}

semiLoofBeam::semiLoofBeam()
  :Element(0, ELE_TAG_semiLoofBeam)
{
  theSections[0] = 0;
  theSections[1] = 0;
  theSections[2] = 0;

  theNodes[0] = 0;
  theNodes[1] = 0;
  theNodes[2] = 0;
}

semiLoofBeam::~semiLoofBeam()
{
  for (int i = 0; i < 3; i++) {
    if (theSections[i])
      delete theSections[i];
    if (theNodes[i])
      delete theNodes[i];
  }
  
  //if (theSections)
  //  delete [] theSections;

  //if (theNodes)
  //  delete [] theNodes;

  if (theLoad)
    delete theLoad;

  if (Ki)
    delete Ki;
  
//  delete [] wbeam;
  delete [] frlof;
  delete [] elxyz;
  delete [] point;
  delete [] atnod;
  delete [] dexce;
  delete [] shear;

}

#ifdef _WIN32

extern "C" int BWBEAM(int *nnode, int ielem, int nelem, double *lnods, double *vanpo,
  double *point, double *xlocal, double *elxyz, double *frlof,double *eside, double *shear,
  double *wbeam, double excen, double *dexcen);

//extern  {
#define cppbeams BEAMS
extern "C" struct { double side; double wbeam[12 * 27]; double wcorn[3 * 2]; double wloof[4 * 2]; } cppbeams;
//}

extern "C" int BLOFBEAM(double *frlof, int *lnods, double *elxyz, double *side,
  double *wbeam, int *noelz, int *lvabz, int *ielem, int *nelem, double *point,
  double *ksi, int *kount, double *atnod, double *excen, double *dexce, double *shear);

extern "C" int LOFBEM(double *elxyz, double *frlof, int *lnods, 
  int *noelz, int *lnomax, int *lvabz, int *ielem, int *nelem, 
  int *nfirst, double *point, double *ksi);

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

#define bwbeam_ 	BWBEAM
#define blofbeam_	BLOFBEAM
#define lofbem_		LOFBEM
#define elmt02_		ELMT02
#define elmt03_		ELMT03
#define elmt04_		ELMT04
#define elmt05_		ELMT05

#else

//extern "C" {
  extern "C" struct cppbeams {double side; double wbeam[12 * 27]; double wcorn[3 * 2]; double wloof[4 * 2];};
//}

extern "C" int bwbeam_(int *nnode, int *ielem, int *nelem, double *lnods, double *vanpo,
  double *point, double *xlocal, double *elxyz, double *frlof,double *eside, double *shear,
  double *wbeam, double *excen, double *dexcen);

extern "C" int blofbeam_(double *frlof, int *lnods, double *elxyz, double *side,
  double *wbeam, int *noelz, int *lvabz, int *ielem, int *nelem, double *point,
  double *ksi, int *kount, double *atnod, double *excen, double *dexce, double *shear);

extern "C" int lofbem_(double *elxyz, double *frlof, int *lnods, 
  int *noelz, int *lnomax, int *lvabz, int *ielem, int *nelem, 
  int *nfirst, double *point, double *ksi);

extern "C" int elmt02_(double *d, double *ul, double *xl, int *ix, double *tl,
  double *s, double *r, int *ndf, int *ndm, int *nst, int *isw);

extern "C" int elmt03_(double *d, double *ul, double *xl, int *ix, double *tl,
  double *s, double *r, int *ndf, int *ndm, int *nst, int *isw);

extern "C" int elmt04_(double *d, double *ul, double *xl, int *ix, double *tl,
  double *s, double *r, int *ndf, int *ndm, int *nst, int *isw);

extern "C" int elmt05_(double *d, double *ul, double *xl, int *ix, double *tl,
  double *s, double *r, int *ndf, int *ndm, int *nst, int *isw);
#endif	       


void
semiLoofBeam::setDomain(Domain *theDomain)
{
  // Check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    theNodes[2] = 0;
    return;
  }

  int Nd1 = connectedExternalNodes(0);
  int Nd2 = connectedExternalNodes(1);
  int Nd3 = connectedExternalNodes(2);

  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);
  theNodes[2] = theDomain->getNode(Nd3);

  for (int i = 0; i < 3; i++){
    Vector crd = theNodes[i]->getCrds();
    elxyz[i] = crd(0); // Fortran and C++ different matrix defination
    elxyz[9 + i] = crd(1);
    elxyz[18 + i] = crd(2);
  }
  atnod[0] = 0.;
  atnod[1] = 0.;
  atnod[2] = 0.;
  dexce[0] = 0.;
  dexce[1] = 0.;
  dexce[2] = 0.;

  if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0) {
    opserr << "FATAL ERROR semiLoofBeam (tag: %d), node not found in domain",
    	this->getTag();
    return;
  }

  if (crdTransf->initialize(theNodes[0], theNodes[2])) {
    // Add some error check
  }

  double L = crdTransf->getInitialLength();

  if (L == 0.0) {
    // Add some error check
  }
  this->DomainComponent::setDomain(theDomain);

  this->update();
}

void
semiLoofBeam::zeroLoad(void)
{
  if (theLoad != 0)
    theLoad->Zero();
}

int
semiLoofBeam::addInertiaLoadToUnbalance(const Vector &accel)
{
  // Check for a quick return
  if (rho == 0.0)
    return 0;

  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);

  if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
    opserr << "Timoshenko3d04::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
    return -1;
  }

  double L = crdTransf->getInitialLength();
  double m = 0.5*rho*L;

  // Want to add ( - fact * M R * accel ) to unbalance
  // Take advantage of lumped mass matrix
  Q(0) -= m*Raccel1(0);
  Q(1) -= m*Raccel1(1);
  Q(2) -= m*Raccel1(2);
  Q(6) -= m*Raccel2(0);
  Q(7) -= m*Raccel2(1);
  Q(8) -= m*Raccel2(2);
  Q(12) -= m*Raccel2(0);
  Q(13) -= m*Raccel2(1);
  Q(14) -= m*Raccel2(2);

  return 0;
}

const Vector &
semiLoofBeam::getResistingForce()
{
  opserr << "semiLoofBeam::addInertiaLoadToUnbalance()" << endln;
  return 0;
}

Node **
semiLoofBeam::getNodePtrs(void)
{
  return theNodes;
}

//get the number of external nodes
int
semiLoofBeam::getNumExternalNodes() const
{
  return 3;
}

//return connected external nodes
const ID&
semiLoofBeam::getExternalNodes()
{
  return connectedExternalNodes;
}

int
semiLoofBeam::getNumDOF(void)
{
  return 18; //17;   //6(A)+5(B with 2 loof)+6(C)
}

//commit state
int
semiLoofBeam::commitState()
{
  int err = 0; // error flag
  int i = 0; // integer for loops
  // call element commitState to do any base class stuff
  if ((err = this->Element::commitState()) != 0) {
    opserr << "semiLoofBeam::commitState () - failed in base class";
    return err;
  }

  // commit the sections
  do {
    err = theSections[i++]->commitState();
  } while (err == 0 && i < 3);

  if (err)
    return err;

  // commit the transformation between coord. systems
  if ((err = crdTransf->commitState()) != 0)
    return err;

  return err;
}

//revert to last commit 
int
semiLoofBeam::revertToLastCommit()
{
  int err;
  int i = 0;

  do {
    err = theSections[i]->revertToLastCommit();
    i++;
  } while (err == 0 && i < 3);

  if (err)
    return err;

  // revert the transformation to last commit
  if ((err = crdTransf->revertToLastCommit()) != 0)
    return err;

  return err;
}

//revert to start 
int
semiLoofBeam::revertToStart()
{
  int err;
  int i, j, k; // for loops
  i = 0;

  // revert the sections state to start
  do {
    err = theSections[i++]->revertToStart();
  } while (err == 0 && i < 3);

  if (err)
    return err;

  // revert the transformation to start
  if ((err = crdTransf->revertToStart()) != 0)
    return err;

  // revert the element state to start

  // Set initial length
  double L = crdTransf->getInitialLength();

  return err;
}

Matrix
semiLoofBeam::getNL(int sec, double L)
{
  return Matrix(1,1);
}

int
semiLoofBeam::update()
{
  int err = 0;
  // Update the transformation
  crdTransf->update();

  // Compute the natural displacements
  // Get basic deformations
  Vector def(18);
  for (int i = 1; i < 3; i++){
    const Vector &v = theNodes[i]->getTrialDisp();
    //def.Assemble(v);
  }
  
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  // Get the numerical integration weights
  double xi[5],weight[5];
  xi[0] = -.577350269; weight[0] = 1.;
  xi[1] = .577350269; weight[1] = 1.;
  xi[2] = -.7745966692; weight[2] = .555555556;
  xi[3] = .7745966692; weight[3] = .555555556;
  xi[4] = .0; weight[4] = .888888889;
  
  // Compute shape functions and their transposes
  // Loop over the integration points
  for (int i = 0; i < 5; i++){
    ksi = xi[i];
    // calculate the WBEAM matrix
//    blofbeam_(frlof, &lnods, elxyz, &side, wbeam, &noelz, &lvabz, &ielem, &nelem, point,
//      &ksi, &kount, atnod, &excen, dexce, shear);

    lofbem_(elxyz, frlof, &lnods, &noelz, &lnomax, &lvabz, 
      &ielem, &nelem, &nfirst, point, &ksi);

    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
    // determine strain vector
    Vector e = Vector(6);
    //Matrix WBEAM = Matrix(wbeam, 12, 27);
    Matrix B = Matrix(6, 17);
    for (int k = 0; k < 17; k++) {
      //cppbeam temp = cppbeam.wbeam[i * 27 + 4 - 1];
      B(1, k) = cppbeams.wbeam[k * 27 + 4 - 1];
      B(2, k) = cppbeams.wbeam[k * 27 + 11 - 1];
      B(3, k) = cppbeams.wbeam[k * 27 + 12 - 1];
      B(4, k) = cppbeams.wbeam[k * 27 + 8 - 1];
      B(5, k) = cppbeams.wbeam[k * 27 + 9 - 1];
      B(6, k) = cppbeams.wbeam[k * 27 + 10 - 1];
    }

    //P, Mz, My, Vy, Vz, T
    for (int j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P:     // axial strain
//        e(j) = oneOverL * v(0); break;
      case SECTION_RESPONSE_MZ:    // curvature
//        e(j) = N1z * v(1) + N2z * v(2); break;
      case SECTION_RESPONSE_MY:    // curvature
//        e(j) = N1y * v(3) + N2y * v(4); break;
      case SECTION_RESPONSE_VY:    // shear strain
//        e(j) = N3y * (v(1) + v(2)); break;
      case SECTION_RESPONSE_VZ:    // shear strain
//        e(j) = N3z * (v(3) + v(4)); break;
      case SECTION_RESPONSE_T:     // Torsion
//        e(j) = oneOverL*v(5); break;
      default:
        break;
      }
    }
    // Set the section deformations
    err += theSections[i]->setTrialSectionDeformation(e);

  }

  if (err != 0) {
    opserr << "semiLoofBeam::update() - failed setTrialSectionDeformations()\n";
    return err;
  }
  return 0;
}

//return secant matrix 
const Matrix&
semiLoofBeam::getInitialStiff()
{
  if (Ki != 0)
    return *Ki;

  static const int ndf = 6; //two membrane plus three bending plus one drill

  static const int nstress = 8; //three membrane, three moment, two shear

  static const int ngauss = 9;

  static const int numnodes = 9;

  int i, j, k, p, q;
  int jj, kk;

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

//   stiff.Zero();


  //compute Jacobian and inverse at center
  double L1 = 0.0;
  double L2 = 0.0;
  //computeJacobian( L1, L2, xl, J0, J0inv ) ; 

  //gauss loop 
  for (i = 0; i < ngauss; i++) {

    //get shape functions
//     shape2d(sg[i], tg[i], xl, shp, xsj);

    //volume element to also be saved
//     dvol[i] = wg[i] * xsj;

    volume += dvol[i];

    // j-node loop to compute strain 
    for (j = 0; j < numnodes; j++)  {

      //compute B matrix 

//       Bmembrane = computeBmembrane(j, shp);

//       Bbend = computeBbend(j, shp);

//       Bshear = computeBshear(j, shp);

//       BJ = assembleB(Bmembrane, Bbend, Bshear);

      //save the B-matrix
      for (p = 0; p < nstress; p++) {
        for (q = 0; q < ndf; q++)
          saveB[p][q][j] = BJ(p, q);
      }//end for p

      //drilling B matrix
//       drillPointer = computeBdrill(j, shp);
      for (p = 0; p < ndf; p++) {
        //BdrillJ[p] = *drillPointer++ ;
        BdrillJ[p] = *drillPointer; //set p-th component
        drillPointer++;             //pointer arithmetic
      }//end for p
    } // end for j


    dd = theSections[i]->getInitialTangent();
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
//       drillPointer = computeBdrill(j, shp);
      for (p = 0; p < ndf; p++) {
        BdrillJ[p] = *drillPointer;
        drillPointer++;
      }//end for p

      //BJtranD = BJtran * dd ;
      BJtranD.addMatrixProduct(0.0, BJtran, dd, 1.0);


      for (p = 0; p < ndf; p++)
//         BdrillJ[p] *= (Ktt*dvol[i]);


      kk = 0;
      for (k = 0; k < numnodes; k++) {

        //extract BK
        for (p = 0; p < nstress; p++) {
          for (q = 0; q < ndf; q++)
            BK(p, q) = saveB[p][q][k];
        }//end for p

        //drilling B matrix
//         drillPointer = computeBdrill(k, shp);
        for (p = 0; p < ndf; p++) {
          BdrillK[p] = *drillPointer;
          drillPointer++;
        }//end for p

        //stiffJK = BJtranD * BK  ;
        // +  transpose( 1,ndf,BdrillJ ) * BdrillK ; 
        stiffJK.addMatrixProduct(0.0, BJtranD, BK, 1.0);

        for (p = 0; p < ndf; p++)  {
          for (q = 0; q < ndf; q++) {
//             stiff(jj + p, kk + q) += stiffJK(p, q)
//               + (BdrillJ[p] * BdrillK[q]);
          }//end for q
        }//end for p
        kk += ndf;
      } // end for k loop
      jj += ndf;
    } // end for j loop
  } //end for i gauss loop 
  Ki = new Matrix(stiff);
  return stiff;
}

const Matrix &
semiLoofBeam::getTangentStiff(void)
{
  // check for quick return
//   if (nen == 0)
//     return (*semiLoofElementM[0]);

  // get the current load factor
  Domain *theDomain = this->getDomain();
  double dm = theDomain->getCurrentTime();

  // set ctan, ior and iow
  double ctan[3];
  ctan[0] = 1.0; ctan[1] = 0.0; ctan[2] = 0.0;
  int ior = 0; int iow = 0;

  // call the ready routine to set ul, xl, tl and ix, NH1, NH2 and NH3
//   int nstR = this->readyfRoutine(false);

  // zero the matrix
//   semiLoofElementM[nstR]->Zero();

  // invoke the Fortran subroutine
  int isw = 3;
//   int nstI = this->invokefRoutine(ior, iow, ctan, isw);

  // check nst is as determined in readyfRoutine()
//   if (nstI != nstR) {
//     opserr << "FATAL semiLoofElement::getTangentStiff() problems with incompatable nst";
//     opserr << " ready: " << nstR << " invoke: " << nstI << endln;
//     exit(-1);
//   }

  // return the matrix

//   return *(semiLoofElementM[nstR]);
  Matrix temp(0, 0);
  return temp;
}

const Matrix &
semiLoofBeam::getDamp(void)
{
  // check for quick return
//   if (nen == 0)
//     return (*semiLoofElementM[0]);

  // get the current load factor
  Domain *theDomain = this->getDomain();
  double dm = theDomain->getCurrentTime();

  // set ctan, ior and iow
  double ctan[3];
  ctan[0] = 0.0; ctan[1] = 1.0; ctan[2] = 0.0;
  int ior = 0; int iow = 0;

  // call the ready routine to set ul, xl, tl and ix, NH1, NH2 and NH3
  int NH1, NH2, NH3;
//   int nstR = this->readyfRoutine(true);

  // zero the matrix
//   semiLoofElementM[nstR]->Zero();

  // invoke the Fortran subroutine
//   int isw = 3; int nst = nen*ndf; int n = this->getTag();

//   int nstI = this->invokefRoutine(ior, iow, ctan, isw);

  // check nst is as determined in readyfRoutine()
//   if (nstI != nstR) {
//     opserr << "FATAL semiLoofElement::getTangentStiff() problems with incompatable nst";
//     opserr << " ready: " << nstR << " invoke: " << nstI << endln;
//     exit(-1);
//   }

  // return the matrix
//   return *(semiLoofElementM[nstR]);
  Matrix temp(0, 0);
  return temp;
}

const Matrix &
semiLoofBeam::getMass(void)
{
  // check for quick return
//   if (nen == 0)
//     return (*semiLoofElementM[0]);

  // get the current load factor
  Domain *theDomain = this->getDomain();
  double dm = theDomain->getCurrentTime();

  // set ctan, ior and iow
  double ctan[3];
  ctan[0] = 0.0; ctan[1] = 0.0; ctan[2] = 1.0;
  int ior = 0; int iow = 0;

  // call the ready routine to set ul, xl, tl and ix, NH1, NH2 and NH3
  int NH1, NH2, NH3;
//   int nstR = this->readyfRoutine(true);

  // zero the matrix and vector (consistent and lumped)
//   semiLoofElementM[nstR]->Zero();
//   semiLoofElementV[nstR]->Zero();

  // invoke the Fortran subroutine
//   int isw = 5; int nst = nen*ndf; int n = this->getTag();
//   int nstI = this->invokefRoutine(ior, iow, ctan, isw);

  // check nst is as determined in readyfRoutine()
//   if (nstI != nstR) {
//     opserr << "FATAL semiLoofElement::getTangentStiff() problems with incompatable nst";
//     opserr << " ready: " << nstR << " invoke: " << nstI << endln;
//     exit(-1);
//   }

  // return the matrix
//   return *(semiLoofElementM[nstR

  return Matrix(0, 0);
}

int
semiLoofBeam::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "semiLoofShell::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
  return -1;
}

//get residual with inertia terms
const Vector&
semiLoofBeam::getResistingForceIncInertia()
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
semiLoofBeam::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  // Now quad sends the ids of its materials
  int matDbTag;
  static ID idData(9);
  int i;

  for (i = 0; i < 3; i++) {
    idData(i) = theSections[i]->getClassTag();
    matDbTag = theSections[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0)
        theSections[i]->setDbTag(matDbTag);
    }
    idData(i + 3) = matDbTag;
  }

  idData(6) = this->getTag();
  idData(7) = connectedExternalNodes(0);
  idData(8) = connectedExternalNodes(1);
  idData(9) = connectedExternalNodes(2);
  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING semiLoofBeam::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  static Vector vectData(4);
  vectData(0) = alphaM;
  vectData(1) = betaK;
  vectData(2) = betaK0;
  vectData(3) = betaKc;

  res += theChannel.sendVector(dataTag, commitTag, vectData);
  if (res < 0) {
    opserr << "WARNING semiLoofBeam::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 3; i++) {
    res += theSections[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING semiLoofBeam::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  return res;
}

int
semiLoofBeam::recvSelf(int commitTag, Channel &theChannel,
FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dataTag = this->getDbTag();
  static ID idData(9);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING semiLoofBeam::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(6));
  connectedExternalNodes(0) = idData(7);
  connectedExternalNodes(1) = idData(8);
  connectedExternalNodes(2) = idData(9);
  static Vector vectData(4);
  res += theChannel.recvVector(dataTag, commitTag, vectData);
  if (res < 0) {
    opserr << "WARNING semiLoofBeam::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  alphaM = vectData(0);
  betaK = vectData(1);
  betaK0 = vectData(2);
  betaKc = vectData(3);

  int i;

  if (theSections[0] == 0) {
    for (i = 0; i < 3; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i + 3);
      // Allocate new material with the sent class tag
      theSections[i] = theBroker.getNewSection(matClassTag);
      if (theSections[i] == 0) {
        opserr << "semiLoofBeam::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;;
        return -1;
      }
      // Now receive materials into the newly allocated space
      theSections[i]->setDbTag(matDbTag);
      res += theSections[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "semiLoofBeam::recvSelf() - material " << i << "failed to recv itself\n";
        return res;
      }
    }
  }
  // Number of materials is the same, receive materials into current space
  else {
    for (i = 0; i < 3; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i + 9);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (theSections[i]->getClassTag() != matClassTag) {
        delete theSections[i];
        theSections[i] = theBroker.getNewSection(matClassTag);
        if (theSections[i] == 0) {
          opserr << "semiLoofBeam::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;
          exit(-1);
        }
      }
      // Receive the material
      theSections[i]->setDbTag(matDbTag);
      res += theSections[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "semiLoofBeam::recvSelf() - material " << i << "failed to recv itself\n";
        return res;
      }
    }
  }
  return res;
}

void
semiLoofBeam::Print(OPS_Stream &s, int flag)
{
  if (flag == 1) {
    s << "\nElement: " << this->getTag() << " Type: semiLoofBeam ";
    s << "\tConnected Nodes: " << connectedExternalNodes;
    s << "\tNumber of Sections: " << 3;
    s << "\tMass density: " << rho;
    for (int i = 0; i < 3; i++)
      s << "\nSection " << i << " :" << *theSections[i];

  }
  else if (flag == 33) {
    s << "\nElement: " << this->getTag() << " Type: semiLoofBeam ";
    double xi[3]; // location of sections or gauss points or integration points
    double wt[3]; // weights of sections or gauss points of integration points
    double L = crdTransf->getInitialLength();
    BeamIntegration *beamInt = new LobattoBeamIntegration();
    beamInt->getSectionLocations(3, L, xi);
    beamInt->getSectionWeights(3, L, wt);
    s << "\n section xi wt";
    for (int i = 0; i < 3; i++)
      s << "\n" << i << " " << xi[i] << " " << wt[i];

  }
  else {
    s << "\nElement: " << this->getTag() << " Type: mixedBeamColumn3d ";
    s << "\tConnected Nodes: " << connectedExternalNodes;
    s << "\tNumber of Sections: " << 3;
    s << "\tMass density: " << rho << endln;
  }
}

int
semiLoofBeam::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
  // first determine the end points of the quad based on
  // the display factor (a measure of the distorted image)
  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();

  static Vector v1(3);
  static Vector v2(3);

  if (displayMode >= 0) {
    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();

    for (int i = 0; i < 3; i++) {
      v1(i) = end1Crd(i) + end1Disp(i)*fact;
      v2(i) = end2Crd(i) + end2Disp(i)*fact;
    }
  }
  else {
    int mode = displayMode * -1;
    const Matrix &eigen1 = theNodes[0]->getEigenvectors();
    const Matrix &eigen2 = theNodes[1]->getEigenvectors();
    if (eigen1.noCols() >= mode) {
      for (int i = 0; i < 3; i++) {
        v1(i) = end1Crd(i) + eigen1(i, mode - 1)*fact;
        v2(i) = end2Crd(i) + eigen2(i, mode - 1)*fact;
      }

    }
    else {
      for (int i = 0; i < 3; i++) {
        v1(i) = end1Crd(i);
        v2(i) = end2Crd(i);
      }
    }
  }
  return theViewer.drawLine(v1, v2, 1.0, 1.0);
}
