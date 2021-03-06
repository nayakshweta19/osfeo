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

// $Revision: 1.1 $
// $Date: 2010-05-04 17:14:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/CompositePackages/mixedBeamColumn2dS.cpp,v $

#include "mixedBeamColumn2dS.h"
#include <elementAPI.h>
#include <OPS_Globals.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <iomanip>

#include <Information.h>
#include <MatrixUtil.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <math.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <iostream>
#include <fstream>
#include <Node.h>
#include <Message.h>

#include <LobattoBeamIntegration.h>
#include <LegendreBeamIntegration.h>
#include <RadauBeamIntegration.h>
#include <NewtonCotesBeamIntegration.h>
#include <TrapezoidalBeamIntegration.h>
#include <RegularizedHingeIntegration.h>

// Constants that define the dimensionality
#define  NDM   2                      // dimension of the problem (2d)
#define  NND   3                      // number of nodal dof's
#define  NEGD  6                      // number of element global dof's
#define  NDM_SECTION  3               // number of section dof's without torsio
#define  NDM_NATURAL  3               // number of element dof's in the basic system without torsion
#define  MAX_NUM_SECTIONS  10         // maximum number of sections allowed

using namespace std;

Matrix mixedBeamColumn2dS::theMatrix(NEGD,NEGD);
Vector mixedBeamColumn2dS::theVector(NEGD);
double mixedBeamColumn2dS::workArea[400];

Vector *mixedBeamColumn2dS::sectionDefShapeFcn = 0;
Matrix *mixedBeamColumn2dS::nldhat = 0;
Matrix *mixedBeamColumn2dS::nd1 = 0;
Matrix *mixedBeamColumn2dS::nd2 = 0;
Matrix *mixedBeamColumn2dS::nd1T = 0;
Matrix *mixedBeamColumn2dS::nd2T = 0;

#ifdef _USRDLL
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
//#define OPS_Export extern "C"
#define OPS_Export
#endif

// Documentation: Two Dimensional Mixed Beam Column Element including Shear Deformation
// element mixedBeamColumn2dS $tag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens>
//   <-integration $intType> <-doRayleigh $rFlag> <-geomLinear>
//
// Required Input Parameters:
//   $tag					integer tag identifying the element
//   $iNode, $jNode         end nodes
//   $numIntgrPts 			number of integration points along the element length
//   $secTag 				identifier for previously-defined section object
//   $transfTag   			identifier for previously-defined coordinate-transformation (CrdTransf) object
//
//
// Optional Input:
//   -mass $massDens
//       $massDens          element mass density (per unit length), from which a lumped-mass matrix is formed (optional, default=0.0)
//   -integration $intType
//       $intType           numerical integration type, options are Lobotto, Legendre, Radau, NewtonCotes, Trapezoidal (optional, default= Lobotto)
//   -doRayleigh $rFlag
//       $rFlag             optional, default = 1
//                              rFlag = 0 no rayleigh damping
//                              rFlag = 1 include rayleigh damping (default)
//   -geomLinear            perform analysis without internal geometric nonlinearity
//
//
// References:
//   1. Bulent N. Alemdar and Donald W. White, “Displacement, Flexibility, and Mixed Beam-Column Finite
//      Element Formulations for Distributed Plasticity Analysis,” Journal of Structural Engineering 131,
//      no. 12 (December 2005): 1811-1819.
//   2. Cenk Tort and Jerome F. Hajjar, “Mixed Finite Element for Three-Dimensional Nonlinear Dynamic
//      Analysis of Rectangular Concrete-Filled Steel Tube Beam-Columns,” Journal of Engineering Mechanics
//      136, no. 11 (November 0, 2010): 1329-1339.
//   3. Denavit, M. D. and Hajjar, J. F. (2010). "Nonlinear Seismic Analysis of Circular Concrete-Filled
//      Steel Tube Members and Frames," Report No. NSEL-023, Newmark Structural Laboratory Report Series
//      (ISSN 1940-9826), Department of Civil and Environmental Engineering, University of Illinois at
//      Urbana-Champaign, Urbana, Illinois, March.


OPS_Export void * OPS_mixedBeamColumn2dS() {
  // Variables to retrieve input
  int iData[10];
  double dData[10];
  int sDataLength = 40;
  char *sData  = new char[sDataLength];
  char *sData2 = new char[sDataLength];
  int numData;

  // Check the number of dimensions
  if (OPS_GetNDM() != NDM) {
   opserr << "ERROR: mixedBeamColumn2dS: invalid number of dimensions\n";
   return 0;
  }

  // Check the number of degrees of freedom
  if (OPS_GetNDF() != NND) {
   opserr << "ERROR: mixedBeamColumn2dS: invalid number of degrees of freedom\n";
   return 0;
  }

  // Check for minimum number of arguments
  if (OPS_GetNumRemainingInputArgs() < 6) {
    opserr << "ERROR: mixedBeamColumn2dS: too few arguments\n";
    return 0;
  }

  numData = 6;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data - mixedBeamColumn2dS\n";
    return 0;
  }

  int eleTag = iData[0];
  int nodeI = iData[1];
  int nodeJ = iData[2];
  int numIntgrPts = iData[3];
  int secTag = iData[4];
  int transfTag = iData[5];

  // Get the section
  SectionForceDeformation *theSection = OPS_GetSectionForceDeformation(secTag);
  if (theSection == 0) {
    opserr << "WARNING section with tag " << secTag << "not found for element " << eleTag << endln;
    return 0;
  }

  SectionForceDeformation **sections = new SectionForceDeformation *[numIntgrPts];
  for (int i = 0; i < numIntgrPts; i++)
    sections[i] = theSection;

  // Get the coordinate transformation
  CrdTransf *theTransf = OPS_GetCrdTransfPtr(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING geometric transformation with tag " << transfTag << "not found for element " << eleTag << endln;
    return 0;
  }


  // Set Default Values for Optional Input
  int doRayleigh = 1;
  double massDens = 0.0;
  bool geomLinear = false;
  BeamIntegration *beamIntegr = 0;

  // Loop through remaining arguments to get optional input
  while ( OPS_GetNumRemainingInputArgs() > 0 ) {
    if ( OPS_GetString(sData, sDataLength) != 0 ) {
      opserr << "WARNING invalid input";
      return 0;
    }

    if ( strcmp(sData,"-mass") == 0 ) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "WARNING invalid input, want: -mass $massDens \n";
      return 0;
      }
      massDens = dData[0];

    } else if ( strcmp(sData,"-integration") == 0 ) {
      if ( OPS_GetString(sData2, sDataLength) != 0 ) {
        opserr << "WARNING invalid input, want: -integration $intType";
        return 0;
      }

      if (strcmp(sData2,"Lobatto") == 0) {
        beamIntegr = new LobattoBeamIntegration();
      } else if (strcmp(sData2,"Legendre") == 0) {
        beamIntegr = new LegendreBeamIntegration();
      } else if (strcmp(sData2,"Radau") == 0) {
        beamIntegr = new RadauBeamIntegration();
      } else if (strcmp(sData2,"NewtonCotes") == 0) {
        beamIntegr = new NewtonCotesBeamIntegration();
      } else if (strcmp(sData2,"Trapezoidal") == 0) {
        beamIntegr = new TrapezoidalBeamIntegration();
      } else if (strcmp(sData2,"RegularizedLobatto") == 0 || strcmp(sData2,"RegLobatto") == 0) {
      numData = 4;
      if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid input, want: -integration RegularizedLobatto $lpI $lpJ $zetaI $zetaJ \n";
        return 0;
      }
      BeamIntegration *otherBeamInt = 0;
      otherBeamInt = new LobattoBeamIntegration();
      beamIntegr = new RegularizedHingeIntegration(*otherBeamInt, dData[0], dData[1], dData[2], dData[3]);
        if (otherBeamInt != 0)
          delete otherBeamInt;
      } else {
      opserr << "WARNING invalid integration type, element: " << eleTag;
      return 0;
      }

    } else if ( strcmp(sData,"-doRayleigh") == 0 ) {
        numData = 1;
        if (OPS_GetInt(&numData, &doRayleigh) != 0) {
          opserr << "WARNING: Invalid doRayleigh in element mixedBeamColumn2dS " << eleTag;
          return 0;
        }

    } else if ( strcmp(sData,"-geomLinear") == 0 ) {
      geomLinear = true;

    } else {
      opserr << "WARNING unknown option " << sData << "\n";
    }
  }

  // Set the beam integration object if not in options
  if (beamIntegr == 0) {
    beamIntegr = new LobattoBeamIntegration();
  }

  // now create the element and add it to the Domain
  Element *theElement = new mixedBeamColumn2dS(eleTag, nodeI, nodeJ, numIntgrPts, sections, *beamIntegr, *theTransf, massDens, doRayleigh, geomLinear);

  if (theElement == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
    return 0;
  }

  return theElement;
}

// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points.
// allocates the necessary space needed by each object
mixedBeamColumn2dS::mixedBeamColumn2dS (int tag, int nodeI, int nodeJ, int numSec,
                                      SectionForceDeformation **sec,
                                      BeamIntegration &bi,
                                      CrdTransf &coordTransf,
                                      double massDensPerUnitLength,
                                      int damp, bool geomLin):
  Element(tag,ELE_TAG_mixedBeamColumn2dS),
  connectedExternalNodes(2), beamIntegr(0), numSections(0), sections(0),
  crdTransf(0), doRayleigh(damp), geomLinear(geomLin),
  rho(massDensPerUnitLength), initialLength(0.0),
  initialFlag(0), itr(0),
  V(NDM_NATURAL), committedV(NDM_NATURAL),
  naturalForce(NDM_NATURAL), commitedNaturalForce(NDM_NATURAL),
  internalForce(NDM_NATURAL), commitedInternalForce(NDM_NATURAL),
  lastNaturalDisp(NDM_NATURAL), commitedLastNaturalDisp(NDM_NATURAL),
  sp(0),
  Hinv(NDM_NATURAL,NDM_NATURAL), commitedHinv(NDM_NATURAL,NDM_NATURAL),
  GMH(NDM_NATURAL,NDM_NATURAL), commitedGMH(NDM_NATURAL,NDM_NATURAL),
  kv(NDM_NATURAL,NDM_NATURAL),
  kvcommit(NDM_NATURAL,NDM_NATURAL),
  Ki(0),
  sectionForceFibers(0), commitedSectionForceFibers(0),
  sectionDefFibers(0), commitedSectionDefFibers(0),
  sectionFlexibility(0), commitedSectionFlexibility(0)
{
   theNodes[0] = 0;
   theNodes[1] = 0;

   connectedExternalNodes(0) = nodeI;
   connectedExternalNodes(1) = nodeJ;

   // get copy of the beam integration object
   beamIntegr = bi.getCopy();
   if(beamIntegr == 0) {
     opserr<<"Error: mixedBeamColumn2dS::mixedBeamColumn2dS: could not create copy of beam integration object" << endln;
     exit(-1);
   }

   // get copy of the transformation object
   crdTransf = coordTransf.getCopy2d();
   if (crdTransf == 0) {
      opserr << "Error: mixedBeamColumn2dS::mixedBeamColumn2dS: could not create copy of coordinate transformation object" << endln;
      exit(-1);
   }


   //this->setSectionPointers(numSec,sec);
   if (numSec > MAX_NUM_SECTIONS) {
     opserr << "Error: mixedBeamColumn2dS::setSectionPointers -- max number of sections exceeded";
   }

   numSections = numSec;

   if (sec == 0) {
     opserr << "Error: mixedBeamColumn2dS::setSectionPointers -- invalid section pointer";
   }

   sections = new SectionForceDeformation *[numSections];
   if (sections == 0) {
     opserr << "Error: mixedBeamColumn2dS::setSectionPointers -- could not allocate section pointers";
   }

   for (int i = 0; i < numSections; i++) {

     if (sec[i] == 0) {
       opserr << "Error: mixedBeamColumn2dS::setSectionPointers -- null section pointer " << i << endln;
     }

     sections[i] = (SectionForceDeformation*) sec[i]->getCopy();

     if (sections[i] == 0) {
       opserr << "Error: mixedBeamColumn2dS::setSectionPointers -- could not create copy of section " << i << endln;
     }

   }

   p0[0] = 0.0;
   p0[1] = 0.0;
   p0[2] = 0.0;

   q0[0] = 0.0;
   q0[1] = 0.0;
   q0[2] = 0.0;

   // Element vectors and matrices
   sectionForceFibers = new Vector [numSections];
   commitedSectionForceFibers = new Vector [numSections];
   sectionDefFibers = new Vector [numSections];
   commitedSectionDefFibers = new Vector [numSections];
   sectionFlexibility = new Matrix [numSections];
   commitedSectionFlexibility = new Matrix [numSections];


   for (int i = 0; i < numSections; i++){
     sectionForceFibers[i] = Vector(NDM_SECTION);
     sectionForceFibers[i].Zero();
     commitedSectionForceFibers[i] = Vector(NDM_SECTION);
     commitedSectionDefFibers[i].Zero();
     sectionDefFibers[i] = Vector(NDM_SECTION);
     sectionDefFibers[i].Zero();
     commitedSectionDefFibers[i] = Vector(NDM_SECTION);
     commitedSectionForceFibers[i].Zero();
     sectionFlexibility[i] = Matrix(NDM_SECTION,NDM_SECTION);
     sectionFlexibility[i].Zero();
     commitedSectionFlexibility[i] = Matrix(NDM_SECTION,NDM_SECTION);
     commitedSectionFlexibility[i].Zero();
   }

   V.Zero();
   naturalForce.Zero();
   internalForce.Zero();
   lastNaturalDisp.Zero();
   Hinv.Zero();
   GMH.Zero();
   kv.Zero();

   committedV.Zero();
   commitedNaturalForce.Zero();
   commitedInternalForce.Zero();
   commitedLastNaturalDisp.Zero();
   commitedHinv.Zero();
   commitedGMH.Zero();
   kvcommit.Zero();

   if (sectionDefShapeFcn == 0)
     sectionDefShapeFcn  = new Vector [MAX_NUM_SECTIONS];
   if (nldhat == 0)
     nldhat  = new Matrix [MAX_NUM_SECTIONS];
   if (nd1 == 0)
     nd1  = new Matrix [MAX_NUM_SECTIONS];
   if (nd2 == 0)
     nd2  = new Matrix [MAX_NUM_SECTIONS];
   if (nd1T == 0)
     nd1T  = new Matrix [MAX_NUM_SECTIONS];
   if (nd2T == 0)
     nd2T  = new Matrix [MAX_NUM_SECTIONS];
   if (!sectionDefShapeFcn || !nldhat || !nd1 || !nd2 || !nd1T || !nd2T ) {
     opserr << "mixedBeamColumn2dS::mixedBeamColumn2dS() -- failed to allocate static section arrays";
     exit(-1);
   }

   int i;
   for ( i=0; i<MAX_NUM_SECTIONS; i++ ){
     nd1T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
     nd2T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
   }
}

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
// CONSTRUCTOR FOR PARALLEL PROCESSING
mixedBeamColumn2dS::mixedBeamColumn2dS():
  Element(0,ELE_TAG_mixedBeamColumn2dS),
  connectedExternalNodes(2), beamIntegr(0), numSections(0), sections(0), crdTransf(0), doRayleigh(0), geomLinear(false),
  rho(0.0), initialLength(0.0),
  initialFlag(0), itr(0),
  V(NDM_NATURAL), committedV(NDM_NATURAL),
  naturalForce(NDM_NATURAL), commitedNaturalForce(NDM_NATURAL),
  internalForce(NDM_NATURAL), commitedInternalForce(NDM_NATURAL),
  lastNaturalDisp(NDM_NATURAL), commitedLastNaturalDisp(NDM_NATURAL),
  sp(0),
  Hinv(NDM_NATURAL,NDM_NATURAL), commitedHinv(NDM_NATURAL,NDM_NATURAL),
  GMH(NDM_NATURAL,NDM_NATURAL), commitedGMH(NDM_NATURAL,NDM_NATURAL),
  kv(NDM_NATURAL,NDM_NATURAL), kvcommit(NDM_NATURAL,NDM_NATURAL),
  Ki(0),
  sectionForceFibers(0), commitedSectionForceFibers(0), sectionDefFibers(0), commitedSectionDefFibers(0),
  sectionFlexibility(0), commitedSectionFlexibility(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  
  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;

  // Element vectors and matrices
  sectionForceFibers = new Vector [numSections];
  commitedSectionForceFibers = new Vector [numSections];
  sectionDefFibers = new Vector [numSections];
  commitedSectionDefFibers = new Vector [numSections];
  sectionFlexibility = new Matrix [numSections];
  commitedSectionFlexibility = new Matrix [numSections];


  for (int i = 0; i < numSections; i++){
    sectionForceFibers[i] = Vector(NDM_SECTION);
    sectionForceFibers[i].Zero();
    commitedSectionForceFibers[i] = Vector(NDM_SECTION);
    commitedSectionForceFibers[i].Zero();
    sectionDefFibers[i] = Vector(NDM_SECTION);
    sectionDefFibers[i].Zero();
    commitedSectionDefFibers[i] = Vector(NDM_SECTION);
    commitedSectionDefFibers[i].Zero();
    sectionFlexibility[i] = Matrix(NDM_SECTION,NDM_SECTION);
    sectionFlexibility[i].Zero();
    commitedSectionFlexibility[i] = Matrix(NDM_SECTION,NDM_SECTION);
    commitedSectionFlexibility[i].Zero();
  }

  V.Zero();
  naturalForce.Zero();
  internalForce.Zero();
  lastNaturalDisp.Zero();
  Hinv.Zero();
  GMH.Zero();
  kv.Zero();

  committedV.Zero();
  commitedNaturalForce.Zero();
  commitedInternalForce.Zero();
  commitedLastNaturalDisp.Zero();
  commitedHinv.Zero();
  commitedGMH.Zero();
  kvcommit.Zero();

  if (sectionDefShapeFcn == 0)
     sectionDefShapeFcn  = new Vector [MAX_NUM_SECTIONS];
  if (nldhat == 0)
     nldhat  = new Matrix [MAX_NUM_SECTIONS];
  if (nd1 == 0)
     nd1  = new Matrix [MAX_NUM_SECTIONS];
  if (nd2 == 0)
     nd2  = new Matrix [MAX_NUM_SECTIONS];
  if (nd1T == 0)
     nd1T  = new Matrix [MAX_NUM_SECTIONS];
  if (nd2T == 0)
     nd2T  = new Matrix [MAX_NUM_SECTIONS];
  if (!sectionDefShapeFcn || !nldhat || !nd1 || !nd2 || !nd1T || !nd2T ) {
    opserr << "mixedBeamColumn2dS::mixedBeamColumn2dS() -- failed to allocate static section arrays";
    exit(-1);
  }

  int i;
  for ( i=0; i<MAX_NUM_SECTIONS; i++ ){
    nd1T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
    nd2T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
  }

}

mixedBeamColumn2dS::~mixedBeamColumn2dS() {

  if (sections) {
    for (int i=0; i < numSections; i++) {
      if (sections[i]) {
        delete sections[i];
      }
    }
    delete [] sections;
  }

  if (crdTransf)
   delete crdTransf;

  if(beamIntegr != 0)
   delete beamIntegr;

  if (sp != 0)
    delete sp;

  if (Ki != 0)
   delete Ki;

  if(sectionForceFibers != 0)
   delete [] sectionForceFibers;

  if(commitedSectionForceFibers != 0)
   delete [] commitedSectionForceFibers;

  if(sectionDefFibers != 0)
   delete [] sectionDefFibers;

  if(commitedSectionDefFibers != 0)
   delete [] commitedSectionDefFibers;

  if(sectionFlexibility != 0)
   delete [] sectionFlexibility;

  if(commitedSectionFlexibility != 0)
   delete [] commitedSectionFlexibility;
}

int 
mixedBeamColumn2dS::getNumExternalNodes(void) const {
   return 2;
}

const ID &
mixedBeamColumn2dS::getExternalNodes(void) {
   return connectedExternalNodes;
}

Node **
mixedBeamColumn2dS::getNodePtrs(void) {
   return theNodes;
}

int 
mixedBeamColumn2dS::getNumDOF(void) {
   return NEGD;
}

void 
mixedBeamColumn2dS::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;

    opserr << "mixedBeamColumn2dS::setDomain:  theDomain = 0 ";
    exit(0);
  }

  // get pointers to the nodes
  int Nd1 = connectedExternalNodes(0);
  int Nd2 = connectedExternalNodes(1);

  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);

  if (theNodes[0] == 0) {
    opserr << "mixedBeamColumn2dS::setDomain: Nd1: ";
    opserr << Nd1 << "does not exist in model\n";
    exit(0);
  }

  if (theNodes[1] == 0) {
    opserr << "mixedBeamColumn2dS::setDomain: Nd2: ";
    opserr << Nd2 << "does not exist in model\n";
    exit(0);
  }

  // call the DomainComponent class method
  this->DomainComponent::setDomain(theDomain);

  // ensure connected nodes have correct number of dof's
  int dofNode1 = theNodes[0]->getNumberDOF();
  int dofNode2 = theNodes[1]->getNumberDOF();

  if ((dofNode1 != NND) || (dofNode2 != NND)) {
    opserr << "mixedBeamColumn2dS::setDomain(): Nd2 or Nd1 incorrect dof ";
    exit(0);
  }

  // initialize the transformation
  if (crdTransf->initialize(theNodes[0], theNodes[1])) {
    opserr << "mixedBeamColumn2dS::setDomain(): Error initializing coordinate transformation";
    exit(0);
  }

  // Check element length
  if (crdTransf->getInitialLength() == 0.0) {
    opserr << "mixedBeamColumn2dS::setDomain(): Zero element length:" << this->getTag();
    exit(0);
  }
}

int 
mixedBeamColumn2dS::commitState() {
  int err = 0; // error flag
  int i = 0; // integers for loops

  // call element commitState to do any base class stuff
  if ((err = this->Element::commitState()) != 0) {
   opserr << "mixedBeamColumn2dS::commitState () - failed in base class";
   return err;
  }

  // commit the sections
  do {
    err = sections[i++]->commitState();
  } while (err == 0 && i < numSections);

  if (err)
    return err;

  // commit the transformation between coord. systems
  if ((err = crdTransf->commitState()) != 0)
    return err;

  // commit the element variables state
  committedV = V;
  commitedNaturalForce = naturalForce;
  commitedInternalForce = internalForce;
  commitedLastNaturalDisp = lastNaturalDisp;
  commitedHinv = Hinv;
  commitedGMH = GMH;
  kvcommit = kv;
  for( i = 0; i < numSections; i++){
    commitedSectionForceFibers[i] = sectionForceFibers[i];
    commitedSectionDefFibers[i] = sectionDefFibers[i];
    commitedSectionFlexibility[i] = sectionFlexibility[i];
  }

  // Reset iteration counter
  itr = 0;

  return err;
}

int 
mixedBeamColumn2dS::revertToLastCommit() {
  int err;
  int i = 0;

  do {
    err = sections[i]->revertToLastCommit();
    i++;
  } while (err == 0 && i < numSections);

  if (err)
    return err;

  // revert the transformation to last commit
  if ((err = crdTransf->revertToLastCommit()) != 0)
    return err;

  // revert the element state to last commit
  V = committedV;
  internalForce = commitedInternalForce;
  naturalForce = commitedNaturalForce;
  lastNaturalDisp = commitedLastNaturalDisp;
  Hinv = commitedHinv;
  GMH = commitedGMH;
  kv   = kvcommit;
  for( i = 0; i < numSections; i++){
    sectionForceFibers[i] = commitedSectionForceFibers[i];
    sectionDefFibers[i] = commitedSectionDefFibers[i];
    sectionFlexibility[i] = commitedSectionFlexibility[i];
  }

  // Reset iteration counter
  itr = 0;

  return err;
}

int 
mixedBeamColumn2dS::revertToStart() {
  int err;
  int i,j,k; // for loops
  i = 0;

  // revert the sections state to start
  do {
     err = sections[i++]->revertToStart();

  }while (err == 0 && i < numSections);

  if (err)
    return err;

  // revert the transformation to start
  if ((err = crdTransf->revertToStart()) != 0)
    return err;

  // revert the element state to start

  // Set initial length
  initialLength = crdTransf->getInitialLength();

  // Get the numerical integration weights
  double wt[MAX_NUM_SECTIONS]; // weights of sections or gauss points of integration points
  beamIntegr->getSectionWeights(numSections, initialLength, wt);

  // Vector of zeros to use at initial natural displacements
  Vector myZeros(NDM_NATURAL);
  myZeros.Zero();

  // Set initial shape functions
  for ( i = 0; i < numSections; i++ ){
    nldhat[i] = this->getNld_hat(i, myZeros, initialLength, geomLinear);
    nd1[i] = this->getNd1(i, myZeros, initialLength, geomLinear);
    nd2[i] = this->getNd2(i, 0, initialLength);

    for( j = 0; j < NDM_SECTION; j++ ){
      for( k = 0; k < NDM_NATURAL; k++ ){
        nd1T[i](k,j) = nd1[i](j,k);
        nd2T[i](k,j) = nd2[i](j,k);
      }
    }
  }

  // Set initial and committed section flexibility and GJ
  Matrix ks(NDM_SECTION,NDM_SECTION);
  for ( i = 0; i < numSections; i++ ){
    getSectionTangent(i,2,ks);
    invertMatrix(NDM_SECTION,ks,sectionFlexibility[i]);
    commitedSectionFlexibility[i] = sectionFlexibility[i];
  }

  // Set initial and committed section forces and deformations
  for ( i = 0; i < numSections; i++ ){
    sectionForceFibers[i].Zero();
    commitedSectionForceFibers[i].Zero();
    sectionDefFibers[i].Zero();
    commitedSectionDefFibers[i].Zero();
  }

  // Compute the following matrices: G, G2, H, H12, H22, Md, Kg
  Matrix G(NDM_NATURAL,NDM_NATURAL);
  Matrix G2(NDM_NATURAL,NDM_NATURAL);
  Matrix H(NDM_NATURAL,NDM_NATURAL);
  Matrix H12(NDM_NATURAL,NDM_NATURAL);
  Matrix H22(NDM_NATURAL,NDM_NATURAL);
  Matrix Md(NDM_NATURAL,NDM_NATURAL);
  Matrix Kg(NDM_NATURAL,NDM_NATURAL);

  G.Zero();
  G2.Zero();
  H.Zero();
  H12.Zero();
  H22.Zero();
  Md.Zero();
  Kg.Zero();
  for( i = 0; i < numSections; i++ ){
    G   = G   + initialLength * wt[i] * nd1T[i] * nldhat[i];
    G2  = G2  + initialLength * wt[i] * nd2T[i] * nldhat[i];
    H   = H   + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd1[i];
    H12 = H12 + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd2[i];
    H22 = H22 + initialLength * wt[i] * nd2T[i] * sectionFlexibility[i] * nd2[i];
    // Md is zero since deformations are zero
    Kg  = Kg  + initialLength * wt[i] * this->getKg(i, 0.0, initialLength);
  }

  // Compute the inverse of the H matrix
  invertMatrix(NDM_NATURAL, H, Hinv);
  commitedHinv = Hinv;

  // Compute the GMH matrix ( G + Md - H12 ) and its transpose
  GMH = G + Md - H12;
  //GMH = G; // Omit P-small delta
  commitedGMH = GMH;

  // Compute the transposes of the following matrices: G2, GMH
  Matrix G2T(NDM_NATURAL,NDM_NATURAL);
  Matrix GMHT(NDM_NATURAL,NDM_NATURAL);
  for( i = 0; i < NDM_NATURAL; i++ ){
    for( j = 0; j < NDM_NATURAL; j++ ){
      G2T(i,j) = G2(j,i);
      GMHT(i,j) = GMH(j,i);
    }
  }

  // Compute the stiffness matrix
  kv.Zero();
  kv = ( Kg + G2 + G2T - H22 ) + GMHT * Hinv * GMH;
  kvcommit = kv;
  Ki = new Matrix(crdTransf->getInitialGlobalStiffMatrix(kv));

  // Vector V is zero at initial state
  V.Zero();
  committedV.Zero();

  // Internal force is zero at initial state
  internalForce.Zero();
  commitedInternalForce.Zero();
  naturalForce.Zero();
  commitedNaturalForce.Zero();

  // Last natural displacement is zero at initial state
  lastNaturalDisp.Zero();
  commitedLastNaturalDisp.Zero();

  // Reset iteration counter
  itr = 0;

  // Set initialFlag to 1 so update doesn't call again
  initialFlag = 1;

  return err;
}

const Matrix & 
mixedBeamColumn2dS::getInitialStiff(void) {
  // If things haven't be initialized, then do so
  if (initialFlag == 0) {
    this->revertToStart();
  }
  return *Ki;
}

const Matrix & 
mixedBeamColumn2dS::getTangentStiff(void) {
  // If things haven't be initialized, then do so
  if (initialFlag == 0) {
    this->revertToStart();
  }
  crdTransf->update();  // Will remove once we clean up the corotational 3d transformation -- MHS
  return crdTransf->getGlobalStiffMatrix(kv,internalForce);
}

const Vector & 
mixedBeamColumn2dS::getResistingForce(void) {
  crdTransf->update();  // Will remove once we clean up the corotational 3d transformation -- MHS
  Vector p0Vec(p0, NDM_NATURAL);
  return crdTransf->getGlobalResistingForce(internalForce, p0Vec);
}

int 
mixedBeamColumn2dS::update() {

  // If things haven't be initialized, then do so
  if (initialFlag == 0) {
    this->revertToStart();
  }

  int i,j,k; // integers for loops

  // Update iteration counter
  // says how many times update has been called since the last commit state
  itr++;

  // Update Coordinate Transformation
  crdTransf->update();

  // Current Length
  double currentLength;
  if (geomLinear) {
    currentLength = initialLength;
  } else {
    currentLength = crdTransf->getDeformedLength();
  }

  // Compute the natural displacements
  Vector naturalDisp = crdTransf->getBasicTrialDisp();


  Vector naturalIncrDeltaDisp(NDM_NATURAL);
  naturalIncrDeltaDisp = naturalDisp - lastNaturalDisp;
  lastNaturalDisp = naturalDisp;

  // Get the numerical integration weights
  double wt[MAX_NUM_SECTIONS]; // weights of sections or gauss points of integration points
  beamIntegr->getSectionWeights(numSections, initialLength, wt);

  // Define Variables
  Vector *sectionForceShapeFcn;
  sectionForceShapeFcn = new Vector [numSections];
  for ( i = 0; i < numSections; i++ ) {
    sectionForceShapeFcn[i] = Vector(NDM_SECTION);
  }

  // Compute shape functions and their transposes
  for ( i = 0; i < numSections; i++ ){
    nldhat[i] = this->getNld_hat(i, naturalDisp, currentLength, geomLinear);
    sectionDefShapeFcn[i] = this->getd_hat(i, naturalDisp, currentLength, geomLinear);
    nd1[i] = this->getNd1(i, naturalDisp, currentLength, geomLinear);
    if (geomLinear) {
      nd2[i].Zero();
    } else {
      nd2[i] = this->getNd2(i, internalForce(0), currentLength);
    }

    for( j = 0; j < NDM_SECTION; j++ ){
      for( k = 0; k < NDM_NATURAL; k++ ){
        nd1T[i](k,j) = nd1[i](j,k);
        nd2T[i](k,j) = nd2[i](j,k);
      }
    }
  }

  // Update natural force
  if (geomLinear) {
    naturalForce = naturalForce + Hinv * ( GMH * naturalIncrDeltaDisp + V );
  } else {
    naturalForce = naturalForce + Hinv * ( GMH * naturalIncrDeltaDisp + V );
  }

  // Update sections
  for ( i = 0; i < numSections; i++){
    // Compute section deformations
    sectionForceShapeFcn[i] = nd1[i] * naturalForce;
    if (sp != 0) {
      const Matrix &s_p = *sp;
      for ( j = 0; j < NDM_SECTION; j++ ) {
        sectionForceShapeFcn[i](j) += s_p(j,i);
      }
    }
    sectionDefFibers[i] = sectionDefFibers[i] + sectionFlexibility[i] * ( sectionForceShapeFcn[i] - sectionForceFibers[i] );

    // Send section deformation to section object
    setSectionDeformation(i,sectionDefFibers[i]);

    // Get section force vector
    getSectionStress(i,sectionForceFibers[i]);

    // Get section tangent matrix
    Matrix ks(NDM_SECTION,NDM_SECTION);
    getSectionTangent(i,1,ks);

    // Compute section flexibility matrix
    invertMatrix(NDM_SECTION,ks,sectionFlexibility[i]);

  }

  // Compute the following matrices: V, V2, G, G2, H, H12, H22, Md, Kg
  Vector V2(NDM_NATURAL);
  Matrix G(NDM_NATURAL,NDM_NATURAL);
  Matrix G2(NDM_NATURAL,NDM_NATURAL);
  Matrix H(NDM_NATURAL,NDM_NATURAL);
  Matrix H12(NDM_NATURAL,NDM_NATURAL);
  Matrix H22(NDM_NATURAL,NDM_NATURAL);
  Matrix Md(NDM_NATURAL,NDM_NATURAL);
  Matrix Kg(NDM_NATURAL,NDM_NATURAL);

  V.Zero();
  V2.Zero();
  G.Zero();
  G2.Zero();
  H.Zero();
  H12.Zero();
  H22.Zero();
  Md.Zero();
  Kg.Zero();

  for( i = 0; i < numSections; i++ ){
    V   = V   + initialLength * wt[i] * nd1T[i] * (sectionDefShapeFcn[i] - sectionDefFibers[i] - sectionFlexibility[i] * ( sectionForceShapeFcn[i] - sectionForceFibers[i] ) );
    V2  = V2  + initialLength * wt[i] * nd2T[i] * (sectionDefShapeFcn[i] - sectionDefFibers[i]);
    G   = G   + initialLength * wt[i] * nd1T[i] * nldhat[i];
    G2  = G2  + initialLength * wt[i] * nd2T[i] * nldhat[i];
    H   = H   + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd1[i];
    H12 = H12 + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd2[i];
    H22 = H22 + initialLength * wt[i] * nd2T[i] * sectionFlexibility[i] * nd2[i];
    if (!geomLinear) {
      Kg = Kg  + initialLength * wt[i] * this->getKg(i, sectionForceFibers[i](0), currentLength);
        // sectionForceFibers[i](0) is the axial load, P
      Md = Md  + initialLength * wt[i] * this->getMd(i, sectionDefShapeFcn[i], sectionDefFibers[i], currentLength);
    }
  }

  // Compute the inverse of the H matrix
  invertMatrix(NDM_NATURAL, H, Hinv);

  // Compute the GMH matrix ( G + Md - H12 ) and its transpose
  GMH = G + Md - H12;
  //GMH = G; // Omit P-small delta

  // Compute the transposes of the following matrices: G, G2, GMH
  Matrix GT(NDM_NATURAL,NDM_NATURAL);
  Matrix G2T(NDM_NATURAL,NDM_NATURAL);
  Matrix GMHT(NDM_NATURAL,NDM_NATURAL);
  for( i = 0; i < NDM_NATURAL; i++ ) {
    for( j = 0; j < NDM_NATURAL; j++ ) {
      GT(i,j) = G(j,i);
      G2T(i,j) = G2(j,i);
      GMHT(i,j) = GMH(j,i);
    }
  }


  // Define the internal force
  if (geomLinear) {
    internalForce = GT * naturalForce + V2 + GMHT * Hinv * V;
  } else {
    internalForce = GT * naturalForce + V2 + GMHT * Hinv * V;
  }


  // Compute the stiffness matrix without the torsion term
  if (geomLinear) {
    kv = ( Kg + G2 + G2T - H22 ) + GMHT * Hinv * GMH;
  } else {
    kv = ( Kg + G2 + G2T - H22 ) + GMHT * Hinv * GMH;
  }

  return 0;
}

const Matrix &
mixedBeamColumn2dS::getMass(void) {
  theMatrix.Zero();

  if (rho != 0.0) {
    theMatrix(0,0) = theMatrix(1,1) =
    theMatrix(3,3) = theMatrix(4,4) = 0.5*initialLength*rho;
  }

  return theMatrix;
}

const Matrix & 
mixedBeamColumn2dS::getDamp(void) {
  theMatrix.Zero();

  // Add the damping forces
  if ( doRayleigh == 1 && (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) ) {
    theMatrix = this->Element::getDamp();
  }

  return theMatrix;
}

void 
mixedBeamColumn2dS::zeroLoad(void) {
  if (sp != 0)
    sp->Zero();

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
}

int 
mixedBeamColumn2dS::addLoad(ElementalLoad *theLoad, double loadFactor) {

  int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  if (sp == 0) {
    sp = new Matrix(NDM_SECTION,numSections);
    if (sp == 0) {
      opserr << "mixedBeamColumn2dS::addLoad -- out of memory\n";
      exit(-1);
    }
  }

  double L = crdTransf->getInitialLength();

  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  if (type == LOAD_TAG_Beam2dUniformLoad) {
    double wy = data(0)*loadFactor;  // Transverse
    double wx = data(1)*loadFactor;  // Axial            // need revision

    Matrix &s_p = *sp;

    // Accumulate applied section forces due to element loads
    for (int i = 0; i < numSections; i++) {
      double x = xi[i]*L;
      // Axial
      s_p(0,i) += wx*(L-x);
      // Moment
      s_p(1,i) += wy*0.5*x*(x-L);
    }

	double V = 0.5*wy*L;
    double M = V*L/6.0; // wt*L*L/12
    double P = wx*L;
    // Accumulate reactions in basic system
    p0[0] -= P;
    p0[1] -= V;
    p0[2] -= V;

	// Fixed end forces in basic system
    q0[0] -= 0.5*P;
    q0[1] -= M;
    q0[2] += M;

  } else if (type == LOAD_TAG_Beam2dPointLoad) {
    double Py = data(0)*loadFactor;
    double N  = data(1)*loadFactor;
    double aOverL = data(2);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL*L;
	double b = L-a;

    double Vy2 = Py*aOverL;
    double Vy1 = Py-Vy2;

    Matrix &s_p = *sp;

    // Accumulate applied section forces due to element loads
    for (int i = 0; i < numSections; i++) {
      double x = xi[i]*L;
      if (x <= a) {
        s_p(0,i) += N;
        s_p(1,i) -= x*Vy1;
      }
      else {
        s_p(1,i) -= (L-x)*Vy2;
      }
    }

    // Accumulate reactions in basic system
    p0[0] -= N;
    p0[1] -= Vy1;
    p0[2] -= Vy2;
	
	double L2 = 1.0/(L*L);
    double a2 = a*a;
    double b2 = b*b;

    // Fixed end forces in basic system
    q0[0] -= N*aOverL;
    double M1 = -a * b2 * Py * L2;
    double M2 = a2 * b * Py * L2;
    q0[1] += M1;
    q0[2] += M2;

  } else {
    opserr << "mixedBeamColumn2dS::addLoad() -- load type unknown for element with tag: " <<
        this->getTag() << endln;

    return -1;
  }

  return 0;
}

const Vector & 
mixedBeamColumn2dS::getResistingForceIncInertia() {

  // Compute the current resisting force
  theVector = this->getResistingForce();

  // Add the inertial forces
  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();

    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;

    theVector(0) += m*accel1(0);
    theVector(1) += m*accel1(1);
    theVector(4) += m*accel2(0);
    theVector(5) += m*accel2(1);
  }

  // Add the damping forces
  if ( doRayleigh == 1 && (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) ) {
    theVector += this->getRayleighDampingForces();
  }

  return theVector;
}

void 
mixedBeamColumn2dS::Print(OPS_Stream &s, int flag)
{

  if (flag == 1) {
    s << "\nElement: " << this->getTag() << " Type: mixedBeamColumn2dS ";
    s << "\tConnected Nodes: " << connectedExternalNodes ;
    s << "\tNumber of Sections: " << numSections;
    s << "\tMass density: " << rho;
    for (int i = 0; i < numSections; i++)
       s << "\nSection "<<i<<" :" << *sections[i];

  } else if (flag == 33) {
    s << "\nElement: " << this->getTag() << " Type: mixedBeamColumn2dS ";
    double xi[MAX_NUM_SECTIONS]; // location of sections or gauss points or integration points
    beamIntegr->getSectionLocations(numSections, initialLength, xi);
    double wt[MAX_NUM_SECTIONS]; // weights of sections or gauss points of integration points
    beamIntegr->getSectionWeights(numSections, initialLength, wt);
    s << "\n section xi wt";
    for (int i = 0; i < numSections; i++)
      s << "\n"<<i<<" "<<xi[i]<<" "<<wt[i];

  } else {
    s << "\nElement: " << this->getTag() << " Type: mixedBeamColumn2dS " << endln;
    s << "\tConnected Nodes: " << connectedExternalNodes << endln;
    s << "\tNumber of Sections: " << numSections << endln;
    s << "\tMass density: " << rho << endln;

	double L = crdTransf->getInitialLength();
	// Axial
	double N = internalForce(0);
	theVector(3) =  N;
	theVector(0) = -N+p0[0];

	// Moments about z and shears along y
	double M1 = internalForce(1);
	double M2 = internalForce(2);
	theVector(2)  = M1;
	theVector(5) = M2;
	double V = (M1+M2)/L;
	theVector(1) =  V+p0[1];
	theVector(4) = -V+p0[2];

	s << "\tEnd 1 Forces (P V M): " << -theVector(0)
		<< " " << theVector(1) << " " << theVector(2) << endln;
	s << "\tEnd 2 Forces (P V M): " << theVector(3)
		<< " " << -theVector(4) << " " << theVector(5) << endln;

	beamIntegr->Print(s, flag);

	for (int i = 0; i < numSections; i++){
	  s << "\t Section \t" << i << endln; 
	  sections[i]->Print(s,flag);
	}
  }
}

OPS_Stream &
operator<<(OPS_Stream &s, mixedBeamColumn2dS &E)
{
    E.Print(s);
    return s;
}

int
mixedBeamColumn2dS::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
  // first determine the end points of the beam based on
  // the display factor (a measure of the distorted image)
  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();	

  static Vector v1(3);
  static Vector v2(3);

  if (displayMode >= 0) {
    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();
    
    for (int i = 0; i < 2; i++) {
      v1(i) = end1Crd(i) + end1Disp(i)*fact;
      v2(i) = end2Crd(i) + end2Disp(i)*fact;    
    }
  } else {
    int mode = displayMode  *  -1;
    const Matrix &eigen1 = theNodes[0]->getEigenvectors();
    const Matrix &eigen2 = theNodes[1]->getEigenvectors();
    if (eigen1.noCols() >= mode) {
      for (int i = 0; i < 2; i++) {
	v1(i) = end1Crd(i) + eigen1(i,mode-1)*fact;
	v2(i) = end2Crd(i) + eigen2(i,mode-1)*fact;    
      }    
    } else {
      for (int i = 0; i < 2; i++) {
	v1(i) = end1Crd(i);
	v2(i) = end2Crd(i);
      }    
    }
  }
  
  return theViewer.drawLine (v1, v2, 1.0, 1.0);
}

Response* 
mixedBeamColumn2dS::setResponse(const char **argv, int argc,
                                         OPS_Stream &output) {

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","mixedBeamColumn2dS");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);

  //
  // we compare argv[0] for known response types
  //

  // global force -
  if (strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"force") == 0 ||
      strcmp(argv[0],"globalForce") == 0 ||
      strcmp(argv[0],"globalForces") == 0) {

    output.tag("ResponseType","Px_1");
    output.tag("ResponseType","Py_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Px_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Mz_2");

    theResponse = new ElementResponse(this, 1, theVector);

  // local force
  }  else if (strcmp(argv[0],"localForce") == 0 ||
              strcmp(argv[0],"localForces") == 0) {

    output.tag("ResponseType","N_1");
    output.tag("ResponseType","V_1");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","N_2");
    output.tag("ResponseType","V_2");
    output.tag("ResponseType","M_2");

    theResponse = new ElementResponse(this, 2, theVector);

  // basic or natural forces
  } else if (strcmp(argv[0],"basicForce") == 0 ||
             strcmp(argv[0],"basicForces") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","M_2");

    theResponse = new ElementResponse(this, 3, Vector(3));

  } else if (strcmp(argv[0],"section") ==0) {
    if (argc > 2) {

      int sectionNum = atoi(argv[1]);
      if (sectionNum > 0 && sectionNum <= numSections) {

        double xi[MAX_NUM_SECTIONS];
        double L = crdTransf->getInitialLength();
        beamIntegr->getSectionLocations(numSections, L, xi);

        output.tag("GaussPointOutput");
        output.attr("number",sectionNum);
        output.attr("eta",xi[sectionNum-1]*L);

        theResponse =  sections[sectionNum-1]->setResponse(&argv[2], argc-2, output);

        output.endTag();
      }
    }
  }

  output.endTag();
  return theResponse;
}

int 
mixedBeamColumn2dS::getResponse(int responseID, Information &eleInfo) {
  if (responseID == 1) { // global forces
    return eleInfo.setVector(this->getResistingForce());

  } else if (responseID == 2) { // local forces
    // Axial
    double N = internalForce(0);
    theVector(3) =  N;
    theVector(0) = -N+p0[0];

    // Moments about z and shears along y
    double M1 = internalForce(1);
    double M2 = internalForce(2);
    theVector(2)  = M1;
    theVector(5) = M2;
    double L = crdTransf->getInitialLength();
    double V = (M1+M2)/L;
    theVector(1) =  V+p0[1];
    theVector(4) = -V+p0[2];

    return eleInfo.setVector(theVector);

  } else if (responseID == 3) { // basic forces
    return eleInfo.setVector(internalForce);

  } else {
    return -1;

  }
}

Vector
mixedBeamColumn2dS::getd_hat(int sec, const Vector &v, double L, bool geomLinear){
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double x, C, E, F;
  Vector D_hat(NDM_SECTION);
  D_hat.Zero();

  x = L*xi[sec];
  C = 1/L;
  E = -4/L + 6*x/(L*L);
  F = -2/L + 6*x/(L*L);

  if (geomLinear) {

    D_hat(0) =  C*v(0);
    D_hat(1) =  E*v(1) + F*v(2);
	D_hat(2) =  E*v(1) + F*v(2); // need revised

  } else {

    double A,B;
    A = 1 - 4*(x/L) + 3*pow(x/L,2);
    B =   - 2*(x/L) + 3*pow(x/L,2);

    D_hat(0) =  C*v(0) +
                0.5*A*A*v(1)*v(1) +
                    A*B*v(1)*v(2) +
                0.5*B*B*v(2)*v(2);
    D_hat(1) =  E*v(1) + F*v(2);
	D_hat(2) =  E*v(1) + F*v(2); // need revised

  }

  return D_hat;
}

Matrix
mixedBeamColumn2dS::getKg(int sec, double P, double L) {
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double x, A, B;

  x = L*xi[sec];
  A = 1 - 4*(x/L) + 3*pow(x/L,2);
  B =   - 2*(x/L) + 3*pow(x/L,2);

  Matrix kg(NDM_NATURAL,NDM_NATURAL);
  kg.Zero();

  kg(0,0) = P / ( L * L );
  kg(1,1) = P*A*A;
  kg(2,2) = P*B*B;
  kg(1,2) = P*A*B;
  kg(2,1) = P*A*B;

  return kg;
}

Matrix 
mixedBeamColumn2dS::getMd(int sec, Vector dShapeFcn, Vector dFibers, double L) {
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double x, A, B;

  Matrix md(NDM_NATURAL,NDM_NATURAL);
  md.Zero();

  x = L*xi[sec];
  A =  ( x/L - 2*pow(x/L,2) + pow(x/L,3) )*L;
  B =          (-pow(x/L,2) + pow(x/L,3) )*L;

  md(0,1) = A * ( dShapeFcn(1) - dFibers(1) );
  md(0,2) = B * ( dShapeFcn(1) - dFibers(1) );

  return md;
}

Matrix 
mixedBeamColumn2dS::getNld_hat(int sec, const Vector &v, double L, bool geomLinear) {
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double x, C, E, F;
  Matrix Nld_hat(NDM_SECTION,NDM_NATURAL);
  Nld_hat.Zero();

  x = L*xi[sec];

  C =  1/L;
  E = -4/L + 6*x/(L*L);
  F = -2/L + 6*x/(L*L);

  if (geomLinear) {

    Nld_hat(0,0) = C;
    Nld_hat(1,1) = E;
    Nld_hat(1,2) = F;
	Nld_hat(2,1) = C; // need revised
	Nld_hat(2,2) = C; // need revised

  } else {

    double A,B;
    A = 1 - 4 * ( x / L ) + 3 * pow ( ( x / L ) , 2 );
    B = - 2 * ( x / L ) + 3 * pow ( ( x / L ) , 2 );
    Nld_hat(0,0) = C + C*C*v(0);
    Nld_hat(0,1) = A*A*v(1) + A*B*v(2);
    Nld_hat(0,2) = A*B*v(1) + B*B*v(2);
    Nld_hat(1,1) = E;
    Nld_hat(1,2) = F;
	Nld_hat(2,1) = C; // need revised
	Nld_hat(2,2) = C; // need revised
  }

  return Nld_hat;
}

Matrix
mixedBeamColumn2dS::getNd2(int sec, double P, double L){
   double xi[MAX_NUM_SECTIONS];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double x, A, B;

   x = L * xi[sec];

   Matrix Nd2(NDM_SECTION,NDM_NATURAL);
   Nd2.Zero();

   A = L * ( x / L - 2 * pow( x / L, 2 ) + pow( x / L, 3 ) );
   B = L * ( -pow( x / L, 2 ) + pow( x / L, 3 ) );

   Nd2(1,1) = P * A;
   Nd2(1,2) = P * B;
   Nd2(2,1) = P * A; // need revised
   Nd2(2,2) = P * B; // need revised

   return Nd2;
}

Matrix
mixedBeamColumn2dS::getNd1(int sec, const Vector &v, double L, bool geomLinear){
   double xi[MAX_NUM_SECTIONS];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double x = L * xi[sec];

   Matrix Nd1(NDM_SECTION,NDM_NATURAL);
   Nd1.Zero();

   if (geomLinear) {

     Nd1(0,0) = 1.0;
     Nd1(1,1) = -x/L + 1.0;
     Nd1(1,2) =  x/L;
	 Nd1(2,1) = 1.0/L; // need revised
	 Nd1(2,2) = 1.0/L; // need revised

   } else {

     double A;
     A = L * ( x/L - 2*pow(x/L,2) + pow(x/L,3) ) * v[1]
              + L * ( -pow(x/L,2) + pow(x/L,3) ) * v[2];

     Nd1(0,0) = 1.0;
     Nd1(1,0) = A;
     Nd1(1,1) = -x/L + 1.0;
     Nd1(1,2) =  x/L;
	 Nd1(2,1) = 1.0/L; // need revised
	 Nd1(2,2) = 1.0/L; // need revised

   }

   return Nd1;
}

void 
mixedBeamColumn2dS::getSectionTangent(int sec,int type,Matrix &kSection) {
  int order = sections[sec]->getOrder();
  const ID &code = sections[sec]->getType();

  // Initialize formulation friendly variables
  kSection.Zero();

  // Get the stress resultant from section
  Matrix sectionTangent(order,order);
  if ( type == 1 ) {
    sectionTangent = sections[sec]->getSectionTangent();
  } else if ( type == 2 ) {
    sectionTangent = sections[sec]->getInitialTangent();
  } else {
    sectionTangent.Zero();
  }

  // Set Components of Section Tangent
  int i,j;
  for (i = 0; i < order; i++) {
    for (j = 0; j < order; j++) {
      switch(code(i)) {
        case SECTION_RESPONSE_P:
          switch(code(j)) {
            case SECTION_RESPONSE_P:
              kSection(0,0) = sectionTangent(i,j);
              break;
            case SECTION_RESPONSE_MZ:
              kSection(0,1) = sectionTangent(i,j);
              break;
			case SECTION_RESPONSE_VY:
              kSection(0,2) = sectionTangent(i,j); // need revised
              break;
            default:
              break;
          }
          break;
        case SECTION_RESPONSE_MZ:
          switch(code(j)) {
            case SECTION_RESPONSE_P:
              kSection(1,0) = sectionTangent(i,j);
              break;
            case SECTION_RESPONSE_MZ:
              kSection(1,1) = sectionTangent(i,j);
              break;
			case SECTION_RESPONSE_VY:
              kSection(1,2) = sectionTangent(i,j); // need revised
              break;
            default:
              break;
          }
          break;
		case SECTION_RESPONSE_VY:
          switch(code(j)) {
            case SECTION_RESPONSE_P:
              kSection(2,0) = sectionTangent(i,j); // need revised
              break;
            case SECTION_RESPONSE_MZ:
              kSection(2,1) = sectionTangent(i,j); // need revised
              break;
			case SECTION_RESPONSE_VY:
              kSection(2,2) = sectionTangent(i,j); // need revised
              break;
            default:
              break;
          }
          break;
        default:
          break;
      }
    }
  }
}

void 
mixedBeamColumn2dS::getSectionStress(int sec,Vector &fSection) {
  int order = sections[sec]->getOrder();
  const ID &code = sections[sec]->getType();

  // Get the stress resultant from section
  Vector stressResultant = sections[sec]->getStressResultant();

  // Initialize formulation friendly variables
  fSection.Zero();

  // Set Components of Section Stress Resultant
  int j;
  for (j = 0; j < order; j++) {
    switch(code(j)) {
      case SECTION_RESPONSE_P:
        fSection(0) = stressResultant(j);
        break;
      case SECTION_RESPONSE_MZ:
        fSection(1) = stressResultant(j);
        break;
	  case SECTION_RESPONSE_VY:
        fSection(2) = stressResultant(j); // need revised
        break; 
      default:
        break;
    }
  }
}

void
mixedBeamColumn2dS::setSectionDeformation(int sec,Vector &defSection) {
  int order = sections[sec]->getOrder();
  const ID &code = sections[sec]->getType();

  // Initialize Section Deformation Vector
  Vector sectionDeformation(order);
  sectionDeformation.Zero();

  // Set Components of Section Deformations
  int j;
  for (j = 0; j < order; j++) {
    switch(code(j)) {
      case SECTION_RESPONSE_P:
        sectionDeformation(j) = defSection(0);
        break;
      case SECTION_RESPONSE_MZ:
        sectionDeformation(j) = defSection(1);
        break;
	  case SECTION_RESPONSE_VY:
        sectionDeformation(j) = defSection(2);  // need revised
        break;
      default:
        break;
    }
  }

  // Set the section deformations
  int res = sections[sec]->setTrialSectionDeformation(sectionDeformation);
}

int 
mixedBeamColumn2dS::sendSelf(int commitTag, Channel &theChannel){
  // @todo write mixedBeamColumn2dS::sendSelf
  opserr << "Error: mixedBeamColumn2dS::sendSelf -- not yet implemented for mixedBeamColumn2dS element";
  return -1;
}

int 
mixedBeamColumn2dS::recvSelf(int commitTag, Channel &theChannel,
                                FEM_ObjectBroker &theBroker){
  // @todo write mixedBeamColumn2dS::recvSelf
  opserr << "Error: mixedBeamColumn2dS::sendSelf -- not yet implemented for mixedBeamColumn2dS element";
  return -1;
}
