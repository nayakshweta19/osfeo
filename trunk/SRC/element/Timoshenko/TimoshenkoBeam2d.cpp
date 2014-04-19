// $Source: /usr/local/cvs/OpenSees/SRC/element/Timoshenko/TimoshenkoBeam2d.cpp,v $
// $Revision: 1.1 $
// $Date: 2009/01/10 21:22:20 $

// Created: 09/09
// Created by: Li Ning (neallee@tju.edu.cn)
// Description: This file contains the class implementation of Timoshenko2d.
//              Make use of Neddy(1997) Interdependent Integration Element 
//              procecess and fiber section model.

// Reference: LI Ning, LI Zhong-Xian, XIE Li-Li. A Fiber-Section Model Based
//            Timoshenko Beam Element Using Shear-Bending Interdependent Shape 
//            Function. Earthquake Engineering & Engineering Vibration. 2013, 
//            12(3): 421-432.

#include <TimoshenkoBeam2d.h>
#include <Node.h>
#include <FiberSection2d.h>
#include <TimoshenkoSection2d.h>
#include <CrdTransf.h>
#include <TimoshenkoLinearCrdTransf2d.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <BeamIntegration.h>

Matrix TimoshenkoBeam2d::K(6,6);
Vector TimoshenkoBeam2d::P(6);
double TimoshenkoBeam2d::workArea[100];

TimoshenkoBeam2d::TimoshenkoBeam2d(int tag, 
					 int nd1, 
					 int nd2,	
					 int numSec, 
					 SectionForceDeformation **s,
					 CrdTransf &coordTransf, 
					 BeamIntegration& bi,
					 double C, double r)
  :Element (tag, ELE_TAG_TimoshenkoBeam2d), 
  numSections(numSec), theSections(0), crdTransf(0), beamInt(0),
  connectedExternalNodes(2), C1(C),
  Q(6), q(6), rho(r)
{
  // Allocate arrays of pointers to SectionForceDeformations
  theSections = new SectionForceDeformation *[numSections];
    
  if (theSections == 0) {
    opserr << "TimoshenkoBeam2d::TimoshenkoBeam2d - failed to allocate section model pointer\n";
    exit(-1);
  }

  for (int i = 0; i < numSections; i++) {
    
    // Get copies of the material model for each integration point
    SectionForceDeformation *theSection = s[i]->getCopy();
    switch (s[i]->getClassTag()) {
	case SEC_TAG_TimoshenkoSection2d:
		theSections[i] = (TimoshenkoSection2d *)theSection;
		break;
	default:
		opserr << "TimoshenkoBeam2d::TimoshenkoBeam2d() --default secTag at sec " << i+1 << endln;
		theSections[i] = theSection;
		break;
    }
  }

  CrdTransf *theCoord = coordTransf.getCopy2d();
  if (theCoord == 0 || theCoord->getClassTag() != CRDTR_TAG_TimoshenkoLinearCrdTransf2d) {
    opserr << "TimoshenkoBeam2d::TimoshenkoBeam2d -- failed to get a copy of coordinate transformation\n";
    if (theCoord == 0)
      opserr << "COPY ZERO\n";
    else
      opserr << "COPY NON _ZERO CLASTAG " << theCoord->getClassTag() << endln;
    exit(-1);
  } 
  
  crdTransf = (TimoshenkoLinearCrdTransf2d *)theCoord;
  
  if (crdTransf == 0) {
    opserr << "TimoshenkoBeam2d::TimoshenkoBeam2d - failed to copy coordinate transformation\n";
    exit(-1);
  }
  
  beamInt = bi.getCopy();

  if (beamInt == 0) {
	  opserr << "DispBeamColumn2d::DispBeamColumn2d - failed to copy beam integration\n";
	  exit(-1);
  }

  // Set connected external node IDs
  connectedExternalNodes(0) = nd1;
  connectedExternalNodes(1) = nd2;

  theNodes[0] = 0;
  theNodes[1] = 0;
  
  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;
  q0[5] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;
  p0[5] = 0.0;

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
// AddingSensitivity:END //////////////////////////////////////
}

TimoshenkoBeam2d::TimoshenkoBeam2d()
:Element (0, ELE_TAG_TimoshenkoBeam2d),
 numSections(0), theSections(0), crdTransf(0), beamInt(0),
 connectedExternalNodes(2),
  Q(6), q(6), rho(0.0)
{
    q0[0] = 0.0;
    q0[1] = 0.0;
    q0[2] = 0.0;
    q0[3] = 0.0;
    q0[4] = 0.0;
    q0[5] = 0.0;

	p0[0] = 0.0;
	p0[1] = 0.0;
	p0[2] = 0.0;
	p0[3] = 0.0;
	p0[4] = 0.0;
	p0[5] = 0.0;

    theNodes[0] = 0;
    theNodes[1] = 0;

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
// AddingSensitivity:END //////////////////////////////////////
}

TimoshenkoBeam2d::~TimoshenkoBeam2d()
{    
    for (int i = 0; i < numSections; i++) {
		if (theSections[i])
			delete theSections[i];
	}

    // Delete the array of pointers to SectionForceDeformation pointer arrays
    if (theSections)
		delete [] theSections;

	if (crdTransf)
		delete crdTransf;

	if (beamInt != 0)
		delete beamInt;
}

int
TimoshenkoBeam2d::getNumExternalNodes() const
{
    return 2;
}

const ID&
TimoshenkoBeam2d::getExternalNodes()
{
    return connectedExternalNodes;
}

Node **
TimoshenkoBeam2d::getNodePtrs()
{
    return theNodes;
}

int
TimoshenkoBeam2d::getNumDOF()
{
    return 6;
}

void
TimoshenkoBeam2d::setDomain(Domain *theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	theNodes[0] = 0;
	theNodes[1] = 0;
	return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);

    if (theNodes[0] == 0 || theNodes[1] == 0) {

	return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    
    if (dofNd1 != 3 || dofNd2 != 3) {
		//opserr << "FATAL ERROR TimoshenkoBeam2d (tag: %d), has differing number of DOFs at its nodes",
		//	this->getTag());
		return;
    }

	if (crdTransf->initialize(theNodes[0], theNodes[1])) {
		// Add some error check
	}

	double L = crdTransf->getInitialLength();

	if (L == 0.0) {
		// Add some error check
	}

    this->DomainComponent::setDomain(theDomain);

	this->update();
}

int
TimoshenkoBeam2d::commitState()
{
    int retVal = 0;


    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "TimoshenkoBeam2d::commitState () - failed in base class";
    }    

    // Loop over the integration points and commit the material states
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->commitState();

    retVal += crdTransf->commitState();					

    return retVal;
}

int
TimoshenkoBeam2d::revertToLastCommit()
{
    int retVal = 0;

	double L = crdTransf->getInitialLength();

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToLastCommit(L);

    retVal += crdTransf->revertToLastCommit();

    return retVal;
}

int
TimoshenkoBeam2d::revertToStart()
{
    int retVal = 0;
	
    //crdTransf->getInitialLength();

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToStart();

    retVal += crdTransf->revertToStart();

    return retVal;
}

int
TimoshenkoBeam2d::update(void)
{
  int err = 0;

  // Update the transformation
  crdTransf->update();
  
  // Get basic deformations
  const Vector &v = crdTransf->getBasicTrialDispInt();		
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;

  double pts[maxNumSections];
  beamInt->getSectionLocations(numSections, L, pts);

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
														
    Vector e(workArea, order);
    
	double xi = 2.0*pts[i]-1.0;

    for (int j = 0; j < order; j++) {
		
      switch(code(j)) {
      case SECTION_RESPONSE_P:   // axial strain
	e(j) = oneOverL*(-v(0)+v(3)); break;
      case SECTION_RESPONSE_MZ:	// curvature
	e(j) = oneOverL*(-1.0+3.0*(1.0-2.0*C1)*xi)*(v(2)-v(5)); break;
	  case SECTION_RESPONSE_VY:	// shear strain
	e(j) = oneOverL*(-v(1)+v(4))-C1*v(2)+(C1-1.0)*v(5); break;
	  default:
	break;
      }
	}

    // Set the section deformations
    err += theSections[i]->setTrialSectionDeformation(e,L);

  }
  if (err != 0) {
    opserr << "TimoshenkoBeam2d::update() - failed setTrialSectionDeformations(e,L)\n";
	return err;
  }

  return 0;
}

const Matrix&
TimoshenkoBeam2d::getTangentStiff(void)
{
  static Matrix kb(6,6);

  // Zero for integral
  kb.Zero();
  q.Zero();
  
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
  double pts[maxNumSections];
  beamInt->getSectionLocations(numSections, L, pts);
  double wts[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wts);
  
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

	//Matrix ka(workArea, order, 3);
	//ka.Zero();

    double xi = 2.0*pts[i]-1.0;

    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getSectionTangent();			
    const Vector &s = theSections[i]->getStressResultant();			
    
    // Perform numerical integration
	//kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
    double wti = 2.*oneOverL; //wts[i]*oneOverL;
	double d11, d12, d13, d21, d22, d23, d31, d32, d33;

	d11 = ks(0,0);	d12 = ks(0,1);	d13 = ks(0,2);	//P
	d21 = ks(1,0);	d22 = ks(1,1);	d23 = ks(1,2);	//M
	d31 = ks(2,0);	d32 = ks(2,1);	d33 = ks(2,2);	//V

	kb(0,0) += wti*(d11);
	kb(0,1) += wti*(d13);
	kb(0,2) += wti*(d21 + C1*d13*L - 3.0*d21*xi + 6.0*C1*d21*xi);
	kb(0,3) += wti*(-d11);
	kb(0,4) += wti*(-d13);
	kb(0,5) += wti*(-((-1 + C1)*d13*L) + d21*(-1.0 + (3.0 - 6.0*C1)*xi));

	kb(1,0) += wti*(d31);
	kb(1,1) += wti*(d33);
	kb(1,2) += wti*(d32 + C1*d33*L - 3.0*d32*xi + 6.0*C1*d32*xi);
	kb(1,3) += wti*(-d31);
	kb(1,4) += wti*(-d33);
	kb(1,5) += wti*(-((-1.0 + C1)*d33*L) + d32*(-1.0 + (3.0 - 6.0*C1)*xi));

	kb(2,0) +=wti*(d21 + C1*d31*L - 3.0*d21*xi + 6.0*C1*d21*xi);
	kb(2,1) +=wti*(d23 + C1*d33*L - 3.0*d23*xi + 6.0*C1*d23*xi);
	kb(2,2) +=wti*(d22*(1.0 + (-3.0 + 6.0*C1)*xi)*(1.0 + (-3.0 + 6.0*C1)*xi) + C1*L*(d23 + d32 + C1*d33*L - 3.0*d23*xi + 6.0*C1*d23*xi - 3.0*d32*xi + 6.0*C1*d32*xi));
	kb(2,3) +=wti*(-d21 - C1*d31*L + 3.0*d21*xi - 6.0*C1*d21*xi);
	kb(2,4) +=wti*(-d23 - C1*d33*L + 3.0*d23*xi - 6.0*C1*d23*xi);
	kb(2,5) +=wti*(-(d22*(1.0 + (-3.0 + 6.0*C1)*xi)*(1.0 + (-3.0 + 6.0*C1)*xi)) - L*((-1.0 + C1)*d23*(1.0 + (-3.0 + 6.0*C1)*xi) + C1*((-1.0 + C1)*d33*L + d32*(1.0 - 3.0*xi + 6.0*C1*xi))));

	kb(3,0) +=wti*(-d11);
	kb(3,1) +=wti*(-d13);
	kb(3,2) +=wti*(-d21 - C1*d13*L + 3.0*d21*xi - 6.0*C1*d21*xi);
	kb(3,3) +=wti*(d11);
	kb(3,4) +=wti*(d13);
	kb(3,5) +=wti*((-1.0 + C1)*d13*L + d21*(1.0 + (-3.0 + 6.0*C1)*xi));

	kb(4,0) +=wti*(-d31);
	kb(4,1) +=wti*(-d33);
	kb(4,2) +=wti*(-d32 - C1*d33*L + 3.0*d32*xi - 6.0*C1*d32*xi);
	kb(4,3) +=wti*(d31);
	kb(4,4) +=wti*(d33);
	kb(4,5) +=wti*((-1.0 + C1)*d33*L + d32*(1.0 + (-3.0 + 6.0*C1)*xi));

	kb(5,0) +=wti*(-((-1.0 + C1)*d31*L) + d21*(-1.0 + (3.0 - 6.0*C1)*xi));
	kb(5,1) +=wti*(-((-1.0 + C1)*d33*L) + d23*(-1.0 + (3.0 - 6.0*C1)*xi));
	kb(5,2) +=wti*(-(d22*(1.0 + (-3.0 + 6.0*C1)*xi)*(1.0 + (-3.0 + 6.0*C1)*xi)) - L*(d32*(-1.0 + 3.0*xi) + C1*(d23 + d32 - d33*L - 3.0*d23*xi - 9.0*d32*xi) + C1*C1*(d33*L + 6.0*(d23 + d32)*xi)));
	kb(5,3) +=wti*((-1.0 + C1)*d31*L + d21*(1.0 + (-3.0 + 6.0*C1)*xi));
	kb(5,4) +=wti*((-1.0 + C1)*d33*L + d23*(1.0 + (-3.0 + 6.0*C1)*xi));
	kb(5,5) +=wti*(d22*(1.0 + (-3.0 + 6.0*C1)*xi)*(1.0 + (-3.0 + 6.0*C1)*xi) + (-1.0 + C1)*L*((-1.0 + C1)*d33*L + d32*(1.0 - 3.0*xi + 6.0*C1*xi) + d23*(1.0 + (-3.0 + 6.0*C1)*xi)));
   
    double s1, s2, s3;
    double wto = wts[i];

	s1=s(0);	s2=s(1);	s3=s(2);

    q(0)+= wto*(-s1);
    q(1)+= wto*(-s3);
    q(2)+= wto*(-s2 - C1*L*s3 + 3.0*s2*xi - 6.0*C1*s2*xi);
    q(3)+= wto*(s1);
    q(4)+= wto*(s3);
    q(5)+= wto*((-1.0 + C1)*L*s3 + s2*(1.0 + (-3.0 + 6.0*C1)*xi));
  }
  
  // Add effects of element loads, q = q(v) + q0		
  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];
  q(5) += q0[5];

  // Transform to global stiffness
  K = crdTransf->getGlobalStiffMatrixInt(kb, q);

  return K;
}

const Matrix&
TimoshenkoBeam2d::getInitialBasicStiff()
{
  static Matrix kb(6,6);

  // Zero for integral
  kb.Zero();
  
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
  //double pts[maxNumSections];
  //beamInt->getSectionLocations(numSections, L, pts);
  //double wts[maxNumSections];
  //beamInt->getSectionWeights(numSections, L, wts);
  
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    //int order = theSections[i]->getOrder();
    //const ID &code = theSections[i]->getType();
  
	//Matrix ka(workArea, order, 6);
	//ka.Zero();

    double xi = 0.; //2.0*pts[i]-1.0;

    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getInitialTangent();
    
    // Perform numerical integration
    //kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);

    double wti = 2.*oneOverL; // wts[i]*oneOverL;
    double C1 = 0.5;

	double d11, d12, d13, d21, d22, d23, d31, d32, d33;
	
	d11 = ks(0,0);	d12 = ks(0,1);	d13 = ks(0,2);	
	d21 = ks(1,0);	d22 = ks(1,1);	d23 = ks(1,2);	
	d31 = ks(2,0);	d32 = ks(2,1);	d33 = ks(2,2);	

	kb(0,0) += wti*(d11);
	kb(0,1) += wti*(d13);
	kb(0,2) += wti*(d21 + C1*d13*L - 3.0*d21*xi + 6.0*C1*d21*xi);
	kb(0,3) += wti*(-d11);
	kb(0,4) += wti*(-d13);
	kb(0,5) += wti*(-((-1 + C1)*d13*L) + d21*(-1.0 + (3.0 - 6.0*C1)*xi));

	kb(1,0) += wti*(d31);
	kb(1,1) += wti*(d33);
	kb(1,2) += wti*(d32 + C1*d33*L - 3.0*d32*xi + 6.0*C1*d32*xi);
	kb(1,3) += wti*(-d31);
	kb(1,4) += wti*(-d33);
	kb(1,5) += wti*(-((-1.0 + C1)*d33*L) + d32*(-1.0 + (3.0 - 6.0*C1)*xi));

	kb(2,0) +=wti*(d21 + C1*d31*L - 3.0*d21*xi + 6.0*C1*d21*xi);
	kb(2,1) +=wti*(d23 + C1*d33*L - 3.0*d23*xi + 6.0*C1*d23*xi);
	kb(2,2) +=wti*(d22*(1.0 + (-3.0 + 6.0*C1)*xi)*(1.0 + (-3.0 + 6.0*C1)*xi) + C1*L*(d23 + d32 + C1*d33*L - 3.0*d23*xi + 6.0*C1*d23*xi - 3.0*d32*xi + 6.0*C1*d32*xi));
	kb(2,3) +=wti*(-d21 - C1*d31*L + 3.0*d21*xi - 6.0*C1*d21*xi);
	kb(2,4) +=wti*(-d23 - C1*d33*L + 3.0*d23*xi - 6.0*C1*d23*xi);
	kb(2,5) +=wti*(-(d22*(1.0 + (-3.0 + 6.0*C1)*xi)*(1.0 + (-3.0 + 6.0*C1)*xi)) - L*((-1.0 + C1)*d23*(1.0 + (-3.0 + 6.0*C1)*xi) + C1*((-1.0 + C1)*d33*L + d32*(1.0 - 3.0*xi + 6.0*C1*xi))));

	kb(3,0) +=wti*(-d11);
	kb(3,1) +=wti*(-d13);
	kb(3,2) +=wti*(-d21 - C1*d13*L + 3.0*d21*xi - 6.0*C1*d21*xi);
	kb(3,3) +=wti*(d11);
	kb(3,4) +=wti*(d13);
	kb(3,5) +=wti*((-1.0 + C1)*d13*L + d21*(1.0 + (-3.0 + 6.0*C1)*xi));

	kb(4,0) +=wti*(-d31);
	kb(4,1) +=wti*(-d33);
	kb(4,2) +=wti*(-d32 - C1*d33*L + 3.0*d32*xi - 6.0*C1*d32*xi);
	kb(4,3) +=wti*(d31);
	kb(4,4) +=wti*(d33);
	kb(4,5) +=wti*((-1.0 + C1)*d33*L + d32*(1.0 + (-3.0 + 6.0*C1)*xi));

	kb(5,0) +=wti*(-((-1.0 + C1)*d31*L) + d21*(-1.0 + (3.0 - 6.0*C1)*xi));
	kb(5,1) +=wti*(-((-1.0 + C1)*d33*L) + d23*(-1.0 + (3.0 - 6.0*C1)*xi));
	kb(5,2) +=wti*(-(d22*(1.0 + (-3.0 + 6.0*C1)*xi)*(1.0 + (-3.0 + 6.0*C1)*xi)) - L*(d32*(-1.0 + 3.0*xi) + C1*(d23 + d32 - d33*L - 3.0*d23*xi - 9.0*d32*xi) + C1*C1*(d33*L + 6.0*(d23 + d32)*xi)));
	kb(5,3) +=wti*((-1.0 + C1)*d31*L + d21*(1.0 + (-3.0 + 6.0*C1)*xi));
	kb(5,4) +=wti*((-1.0 + C1)*d33*L + d23*(1.0 + (-3.0 + 6.0*C1)*xi));
	kb(5,5) +=wti*(d22*(1.0 + (-3.0 + 6.0*C1)*xi)*(1.0 + (-3.0 + 6.0*C1)*xi) + (-1.0 + C1)*L*((-1.0 + C1)*d33*L + d32*(1.0 - 3.0*xi + 6.0*C1*xi) + d23*(1.0 + (-3.0 + 6.0*C1)*xi)));

  }

  return kb;
}

const Matrix&
TimoshenkoBeam2d::getInitialStiff()
{
  const Matrix &kb = this->getInitialBasicStiff();

  // Transform to global stiffness
  K = crdTransf->getInitialGlobalStiffMatrixInt(kb);

  return K;
}

const Matrix&
TimoshenkoBeam2d::getMass()
{
  K.Zero();

  if (rho == 0.0)
    return K;
  
  double L = crdTransf->getInitialLength();
  double m = 0.5*rho*L;
  
  K(0,0) = K(1,1) = K(3,3) = K(4,4) = m;
  
  return K;
}

void
TimoshenkoBeam2d::zeroLoad(void)
{
  Q.Zero();

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;
  q0[5] = 0.0;
  
  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;
  p0[5] = 0.0;

  return;
}

int 
TimoshenkoBeam2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L = crdTransf->getInitialLength();
  
  if (type == LOAD_TAG_Beam2dUniformLoad) {
    double wt = data(0)*loadFactor;  // Transverse (+ve upward)
    double wa = data(1)*loadFactor;  // Axial (+ve from node I to J)

	double V = 0.5*wt*L;
	double M = V*L/6.0; // wt*L*L/12
	double P = wa*L;

	// Reactions in basic system
	p0[0] -= P;
	p0[1] -= V;
	p0[2] -= M;
	p0[3] -= P;
	p0[4] -= V;
	p0[5] -= -M; //add and modified

    // Fixed end forces in basic system
    q0[0] += wa*L*0.5;
    q0[1] += wt*L*0.5;
    q0[2] += wt*L*L/12.0;
    q0[3] += wa*L*0.5;
    q0[4] += wt*L*0.5;
    q0[5] += -wt*L*L/12.0;

  }
  else if (type == LOAD_TAG_Beam2dPointLoad) {
    double V = data(0)*loadFactor;
    double N = data(1)*loadFactor;
    double aOverL = data(2);					

	if (aOverL < 0.0 || aOverL > 1.0)
		return 0;

	double a = aOverL*L;
	double b = L-a;

	// Reactions in basic system
	p0[0] -= N;
	double V1 = V*(1.0-aOverL);
	double V2 = V*aOverL;
	p0[1] -= V1;
	p0[2] -= V2;

    // Fixed end forces in basic system
   	
	double M1 = L*V*aOverL*(1.0-aOverL)*(1.0-aOverL+aOverL*2.0);
    q0[0] += N*(1.0-aOverL);
    q0[1] += V*(1.0-aOverL);
    q0[2] += M1;
    q0[3] += N*aOverL;
    q0[4] += V*aOverL;
    q0[5] += -M1;

  }
  else {
    opserr << "TimoshenkoBeam2d::TimoshenkoBeam2d -- load type unknown for element with tag: "
	   << this->getTag() << "TimoshenkoBeam2d::addLoad()\n"; 
			    
    return -1;
  }

  return 0;
}

int 
TimoshenkoBeam2d::addInertiaLoadToUnbalance(const Vector &accel)
{
	// Check for a quick return
	if (rho == 0.0) 
		return 0;

	// Get R * accel from the nodes
	const Vector &Raccel1 = theNodes[0]->getRV(accel);
	const Vector &Raccel2 = theNodes[1]->getRV(accel);

    if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
      opserr << "TimoshenkoBeam2d::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
      return -1;
    }

	double L = crdTransf->getInitialLength();
	double m = 0.5*rho*L;

    // Want to add ( - fact * M R * accel ) to unbalance
	// Take advantage of lumped mass matrix
	Q(0) -= m*Raccel1(0);
	Q(1) -= m*Raccel1(1);
	Q(3) -= m*Raccel2(0);
	Q(4) -= m*Raccel2(1);

    return 0;
}

const Vector&
TimoshenkoBeam2d::getResistingForce()
{
  double L = crdTransf->getInitialLength();

  double pts[maxNumSections];
  beamInt->getSectionLocations(numSections, L, pts);
  double wts[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wts);
  
  // Zero for integration
  q.Zero();
  
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    //int order = theSections[i]->getOrder();
    //const ID &code = theSections[i]->getType();
  
    //double xi = 2.0*pts[i]-1.0;
	double xi = pts[i];

    // Get section stress resultant
    const Vector &s = theSections[i]->getStressResultant();
    
    // Perform numerical integration on internal force
       
	//q.addMatrixTransposeVector(1.0, *B, s, wts(i));

	double s1, s2, s3;
	double wto = wts[i];
	  
	s1=s(0);
	s2=s(1);
	s3=s(2);

	q(0)+= wto*(-s1);
	q(1)+= wto*(-s3);
	q(2)+= wto*(-s2 - 0.5*L*s3);
	q(3)+= wto*(s1);
	q(4)+= wto*(s3);
	q(5)+= wto*(s2 -0.5*L*s3);

  }
  
  // Add effects of element loads, q = q(v) + q0		
  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];
  q(5) += q0[5];


  // Vector for reactions in basic system
  Vector p0Vec(p0, 6);
  P = crdTransf->getGlobalResistingForceInt(q, p0Vec);	
  
  // Subtract other external nodal loads ... P_res = P_int - P_ext
  P(0) -= Q(0);
  P(1) -= Q(1);
  P(2) -= Q(2);
  P(3) -= Q(3);
  P(4) -= Q(4);
  P(5) -= Q(5);
  
  return P;
}

const Vector&
TimoshenkoBeam2d::getResistingForceIncInertia()
{

  this->getResistingForce();
  
  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();
    
    // Compute the current resisting force
    this->getResistingForce();
    
    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;
    
    P(0) += m*accel1(0);
    P(1) += m*accel1(1);
    P(3) += m*accel2(0);
    P(4) += m*accel2(1);
    
    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P += this->getRayleighDampingForces();

  } else {
    
    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P += this->getRayleighDampingForces();
  }

  return P;
}

int
TimoshenkoBeam2d::sendSelf(int commitTag, Channel &theChannel)
{
  // place the integer data into an ID

  int dbTag = this->getDbTag();
  int i, j;
  int loc = 0;
  
  static ID idData(7);  // one bigger than needed so no clash later
  idData(0) = this->getTag();
  idData(1) = connectedExternalNodes(0);
  idData(2) = connectedExternalNodes(1);
  idData(3) = numSections;
  idData(4) = crdTransf->getClassTag();
  int crdTransfDbTag  = crdTransf->getDbTag();
  if (crdTransfDbTag  == 0) {
    crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag  != 0) 
      crdTransf->setDbTag(crdTransfDbTag);
  }
  idData(5) = crdTransfDbTag;
  
  if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
    opserr << "TimoshenkoBeam2d::sendSelf() - failed to send ID data\n";
     return -1;
  }    

  // send the coordinate transformation
  
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
     opserr << "TimoshenkoBeam2d::sendSelf() - failed to send crdTranf\n";
     return -1;
  }      

  
  //
  // send an ID for the sections containing each sections dbTag and classTag
  // if section ha no dbTag get one and assign it
  //

  ID idSections(2*numSections);
  loc = 0;
  for (i = 0; i<numSections; i++) {
    int sectClassTag = theSections[i]->getClassTag();
    int sectDbTag = theSections[i]->getDbTag();
    if (sectDbTag == 0) {
      sectDbTag = theChannel.getDbTag();
      theSections[i]->setDbTag(sectDbTag);
    }

    idSections(loc) = sectClassTag;
    idSections(loc+1) = sectDbTag;
    loc += 2;
  }

  if (theChannel.sendID(dbTag, commitTag, idSections) < 0)  {
    opserr << "DispBeamColumn2d::sendSelf() - failed to send ID data\n";
    return -1;
  }    

  // send the sections
  
  for (j = 0; j<numSections; j++) {
    if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "TimoshenkoBeam2d::sendSelf() - section " << 
	j << "failed to send itself\n";
      return -1;
    }
  }

  return 0;
}

int
TimoshenkoBeam2d::recvSelf(int commitTag, Channel &theChannel,
			   FEM_ObjectBroker &theBroker)
{
  // receive the integer data containing tag, numSections and coord transformation info
 
	int dbTag = this->getDbTag();
  int i;
  
  static ID idData(7); // one bigger than needed so no clash with section ID

  if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
    opserr << "TimoshenkoBeam2d::recvSelf() - failed to recv ID data\n";
    return -1;
  }    

  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);
  
  int crdTransfClassTag = idData(4);
  int crdTransfDbTag = idData(5);

  // create a new crdTransf object if one needed
  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
      if (crdTransf != 0)
	  delete crdTransf;

      crdTransf = new TimoshenkoLinearCrdTransf2d();

      if (crdTransf == 0) {
	opserr << "DispBeamColumn2d::recvSelf() - failed to obtain a CrdTrans object with classTag " <<
	  crdTransfClassTag << endln;
	  return -2;	  
      }
  }

  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf object
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "TimoshenkoBeam2d::sendSelf() - failed to recv crdTranf\n";
    return -3;
  }      
  
  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2*idData(3));
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "TimoshenkoBeam2d::recvSelf() - failed to recv ID data\n";
    return -1;
  }    

    // now receive the sections
   
  if (numSections != idData(3)) {

    //
    // we do not have correct number of sections, must delete the old and create
    // new ones before can recvSelf on the sections
    //

    // delete the old
    if (numSections != 0) {
      for (i=0; i<numSections; i++)
	delete theSections[i];
      delete [] theSections;
    }

    // create a new array to hold pointers
    theSections = new SectionForceDeformation *[idData(3)];
    if (theSections == 0) {
opserr << "TimoshenkoBeam2d::recvSelf() - out of memory creating sections array of size " <<
  idData(3) << endln;
      return -1;
    }    

    // create a section and recvSelf on it
    numSections = idData(3);
    loc = 0;
    
    for (i=0; i<numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;
      theSections[i] = new TimoshenkoSection2d();
      if (theSections[i] == 0) {
	opserr << "TimoshenkoBeam2d::recvSelf() - Broker could not create Section of class type " <<
	  sectClassTag << endln;
	exit(-1);
      }
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "TimoshenkoBeam2d::recvSelf() - section " << i << " failed to recv itself\n";
	return -1;
      }     
    }

  } else {

    // for each existing section, check it is of correct type
    // (if not delete old & create a new one) then recvSelf on it
      
    loc = 0;
    for (i=0; i<numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;

      // check of correct type
      if (theSections[i]->getClassTag() !=  sectClassTag) {
	// delete the old section[i] and create a new one
	delete theSections[i];
	theSections[i] = new TimoshenkoSection2d();
	if (theSections[i] == 0) {
	opserr << "TimoshenkoBeam2d::recvSelf() - Broker could not create Section of class type " <<
	  sectClassTag << endln;
	exit(-1);
	}
      }

      // recvSelf on it
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "TimoshenkoBeam2d::recvSelf() - section " << i << " failed to recv itself\n";
	return -1;
      }     
    }
  }

  return 0;
}

void
TimoshenkoBeam2d::Print(OPS_Stream &s, int flag)		
{
  s << "\nTimoshenkoBeam2d, element id:  " << this->getTag() << endln;
  s << "\tConnected external nodes:  " << connectedExternalNodes;
  s << "\tCoordTransf: " << crdTransf->getTag() << endln;
  s << "\tmass density:  " << rho << endln;

  double L = crdTransf->getInitialLength();
  double P  = q(0);
  double M1 = q(1);
  double M2 = q(2);
  double V = (M1+M2)/L;

  s << "\tEnd 1 Forces (P V M): " << -q(0)
    << " " << q(1) << " " << q(2) << endln;
  s << "\tEnd 2 Forces (P V M): " << q(3)
    << " " << -q(4) << " " << q(5) << endln;

  beamInt->Print(s, flag);

  for (int i = 0; i < numSections; i++)
	  theSections[i]->Print(s,flag);
}


int
TimoshenkoBeam2d::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	

    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();

	static Vector v1(3);
	static Vector v2(3);

	for (int i = 0; i < 2; i++) {
		v1(i) = end1Crd(i) + end1Disp(i)*fact;
		v2(i) = end2Crd(i) + end2Disp(i)*fact;    
	}
	
	return theViewer.drawLine (v1, v2, 1.0, 1.0);
}

Response*
TimoshenkoBeam2d::setResponse(const char **argv, int argc, OPS_Stream &s)
{
    // global force - 
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
		|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
		return new ElementResponse(this, 1, P);

    // local force -
    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)
		return new ElementResponse(this, 2, P);

    // chord rotation -
    else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0
	     || strcmp(argv[0],"basicDeformation") == 0)
      return new ElementResponse(this, 3, Vector(3));
    
    // plastic rotation -
    else if (strcmp(argv[0],"plasticRotation") == 0 || strcmp(argv[0],"plasticDeformation") == 0)
      return new ElementResponse(this, 4, Vector(3));
    
    // section response -
    else if (strcmp(argv[0],"section") == 0 || strcmp(argv[0],"-section") == 0) {
      if (argc <= 2)
	return 0;
      
      int sectionNum = atoi(argv[1]);
      if (sectionNum > 0 && sectionNum <= numSections)
	return theSections[sectionNum-1]->setResponse(&argv[2], argc-2, s);
      else
	return 0;
    }
    
    else
      return 0;
}

int 
TimoshenkoBeam2d::getResponse(int responseID, Information &eleInfo)		//L modify
{
  //double V;
  double L = crdTransf->getInitialLength();

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 2) {
      P(3) =  q(3);
      P(0) = -q(0);
      P(2) = q(2);
      P(5) = q(5);
      P(1) =  q(1);
      P(4) = -q(4);
      return eleInfo.setVector(P);
  }

  // Chord rotation
  else if (responseID == 3)
    return eleInfo.setVector(crdTransf->getBasicTrialDispInt());

  // Plastic rotation
  else if (responseID == 4) {
    static Vector vp(6);
    static Vector ve(6);
    const Matrix &kb = this->getInitialBasicStiff();
    kb.Solve(q, ve);
    vp = crdTransf->getBasicTrialDispInt();
    vp -= ve;
    return eleInfo.setVector(vp);
  }

  else
    return -1;
}


// AddingSensitivity:BEGIN ///////////////////////////////////

const Matrix &
TimoshenkoBeam2d::getKiSensitivity(int gradNumber)
{
	K.Zero();
	return K;
}

const Matrix &
TimoshenkoBeam2d::getMassSensitivity(int gradNumber)
{
	K.Zero();
	return K;
}


const Vector &
TimoshenkoBeam2d::getResistingForceSensitivity(int gradNumber)
{
	static Vector dummy(3);		// No distributed loads
	return dummy;
}


// NEW METHOD
int
TimoshenkoBeam2d::commitSensitivity(int gradNumber, int numGrads)
{
 	return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////


