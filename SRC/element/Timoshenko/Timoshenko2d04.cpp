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
                                                                        
// $Source: /usr/local/cvs/OpenSees/SRC/element/Timoshenko/Timoshenko2d04.cpp,v $
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

#include <Timoshenko2d04.h>
#include <Node.h>
#include <FiberSection2d.h>
#include <SectionForceDeformation.h>
#include <TimoshenkoSection2d.h>
#include <CrdTransf.h>
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
#include <CompositeResponse.h>
#include <ElementalLoad.h>
#include <BeamIntegration.h>
#include <Parameter.h>
#include <math.h>

Matrix Timoshenko2d04::K(6,6);
Vector Timoshenko2d04::P(6);
double Timoshenko2d04::workArea[100];

Matrix *Timoshenko2d04::bd = 0;
Matrix *Timoshenko2d04::nd = 0;

#define MAXITER 40

Timoshenko2d04::Timoshenko2d04(int tag, 
					 int nd1, 
					 int nd2,	
					 int numSec,
					 SectionForceDeformation **s,
					 CrdTransf &coordTransf, 
					 BeamIntegration& bi,
					 double r, double SCF)
    :Element (tag, ELE_TAG_Timoshenko2d04),
    numSections(numSec), theSections(0), crdTransf(0), beamInt(0),
    connectedExternalNodes(2),
    Q(6), q(3), rho(r), shearCF(SCF), Omega(0.0)
{
  // Allocate arrays of pointers to SectionForceDeformations
  theSections = new SectionForceDeformation *[numSections];
    
  if (theSections == 0) {
    opserr << "Timoshenko2d04::Timoshenko2d04 - failed to allocate section model pointer\n";
    exit(-1);
  }

  for (int i = 0; i< numSections; i++){
    // Get copies of the material model for each integration point
    //SectionForceDeformation *theSection = s[i]->getCopy();
    //switch (s[i]->getClassTag()) {
	//case SEC_TAG_TimoshenkoSection2d:
		theSections[i] =  s[i]->getCopy(); //(TimoshenkoSection2d *)theSection;
	//	break;
	//default:
	//	opserr << "Timoshenko2d04::Timoshenko2d04() --default secTag at sec " << i+1 << endln;
	//	theSections[i] = theSection;
	//	break;
    //}
	// Check allocation
    if (theSections[i] == 0 ) { //|| s[i]->getClassTag() != SEC_TAG_TimoshenkoSection2d
      opserr << "Timoshenko2d04::Timoshenko2d04() --failed to get a copy of section model " << endln;
      exit(-1);
    }
  }

  crdTransf = coordTransf.getCopy2d();
  
  if (crdTransf == 0) {
    opserr << "Timoshenko2d04::Timoshenko2d04 - failed to copy coordinate transformation\n";
    exit(-1);
  }
  
  beamInt = bi.getCopy();

  if (beamInt == 0) {
	  opserr << "Timoshenko2d04::Timoshenko2d04() - failed to copy beam integration\n";
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

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  if (nd == 0) 	nd = new Matrix [maxNumSections];
  if (bd == 0)	bd = new Matrix [maxNumSections];

  if (!nd || !bd ) {
    opserr << "Timoshenko2d04::Timoshenko2d04() -- failed to allocate static section arrays";
    exit(-1);
  }

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
// AddingSensitivity:END //////////////////////////////////////
}

Timoshenko2d04::Timoshenko2d04()
	:Element (0, ELE_TAG_Timoshenko2d04),
	numSections(0), theSections(0), crdTransf(0), beamInt(0),
	connectedExternalNodes(2),
	Q(6), q(3), rho(0.0), Omega(0.0)
{
  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  theNodes[0] = 0;
  theNodes[1] = 0;

  if (nd == 0) 	nd  = new Matrix [maxNumSections];
  if (bd == 0)	bd  = new Matrix [maxNumSections];

  if (!nd || !bd ) {
    opserr << "Timoshenko2d04::Timoshenko2d04() -- failed to allocate static section arrays";
    exit(-1);
  }

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
// AddingSensitivity:END //////////////////////////////////////
}

Timoshenko2d04::~Timoshenko2d04()
{    
    for (int i = 0; i<numSections; i++)
	  if (theSections[i])
		delete theSections[i];
	
    // Delete the array of pointers to SectionForceDeformation pointer arrays
    if (theSections)
		delete [] theSections;

	if (crdTransf)
		delete crdTransf;

	if (beamInt != 0)
		delete beamInt;
}

int
Timoshenko2d04::getNumExternalNodes() const
{
    return 2;
}

const ID&
Timoshenko2d04::getExternalNodes()
{
    return connectedExternalNodes;
}

Node **
Timoshenko2d04::getNodePtrs()
{
    return theNodes;
}

int
Timoshenko2d04::getNumDOF()
{
    return 6;
}

void
Timoshenko2d04::setDomain(Domain *theDomain)
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
	  opserr << "FATAL ERROR Timoshenko2d04 (tag: %d), has differing number of DOFs at its nodes",
	  this->getTag();
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
Timoshenko2d04::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "Timoshenko2d04::commitState () - failed in base class";
    }    

    // Loop over the integration points and commit the material states
    for (int i = 0; i<numSections; i++)
	  retVal += theSections[i]->commitState();

    retVal += crdTransf->commitState();					

    return retVal;
}

int
Timoshenko2d04::revertToLastCommit()
{
    int retVal = 0;

	double L = crdTransf->getInitialLength();

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i<numSections; i++) 
	  retVal += theSections[i]->revertToLastCommit(L);

    retVal += crdTransf->revertToLastCommit();

    return retVal;
}

int
Timoshenko2d04::revertToStart()
{
    int retVal = 0;
	
    //crdTransf->getInitialLength();

    // Loop over the integration points and revert states to start
    for (int i = 0; i<numSections; i++)
	  retVal += theSections[i]->revertToStart();

    retVal += crdTransf->revertToStart();

    return retVal;
}

int
Timoshenko2d04::update(void)
{
  int err = 0;

  // Update the transformation
  crdTransf->update();
  
  // Get basic deformations
  const Vector &v = crdTransf->getBasicTrialDisp();		
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;

  double pts[maxNumSections];
  beamInt->getSectionLocations(numSections, L, pts);
  double wts[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wts);
  int order = theSections[0]->getOrder();     // Section 0 for all
  const ID &code = theSections[0]->getType(); // Section 0 for all

  double mu, x, N1, N2, N3, error = 1.0, temp=0., tol = 1.e-4; //N4,
  double EI, GA, Omega1, Omega2;
  //Vector Rslt(3);
  Vector e(workArea, order);

  int iterNum=0;
  do {
	//opserr << "\tTrial Omega: "<< Omega << endln;
	temp =0.0;
    // Loop over the integration points
    for (int i = 0; i<numSections; i++) {
	  //opserr << "\t\tSection point: " << i << endln;
      mu    = 1./(1.+12.*Omega);
      x     = L * pts[i];
      //phi1  =  mu*x*(L-x)*(L-x+6.*L*Omega)                      /L/L;
      //phi1p =  mu*(3.*x*x+L*L*(1.+6.*Omega)-4.*L*(x+3.*x*Omega))/L/L;
      //phi2  = -mu*x*(L-x)*(x + 6.*L*Omega)                      /L/L;
      //phi2p =  mu*(3.*x*x-L*L* 6. *Omega +2.*L*x*(6.*Omega-1.))/L/L;
      //phi3  =  mu*(L-x)*(L-3.*x+12.*L*Omega)                    /L/L;
      //phi3p =  mu*(6.*x - 4.*L*(1.+3.*Omega))                   /L/L;
      //phi4  =  mu*x*(  3.*x+2.*L*(6.*Omega-1.))                 /L/L;
      //phi4p =  mu*2.*(3.*x+L*(6.*Omega-1.))                     /L/L;
      N1 = mu*(6.*x-4.*L*(1.+3.*Omega))/L/L;
	  N2 = 2.*mu*(3.*x+L*(6.*Omega-1.))/L/L;
	  N3 = 6.*mu*Omega; //*2.*(1-x/L)
	  //N4 = 6.*mu*Omega;   //*2.*(x/L)
      for (int j = 0; j < order; j++) {
        switch(code(j)) {
        case SECTION_RESPONSE_P:     // axial strain
      e(j) = oneOverL*v(0); break;
        case SECTION_RESPONSE_MZ:    // curvature
      e(j) = N1 * v(1) + N2* v(2); break;
        case SECTION_RESPONSE_VY:    // shear strain
      e(j) = N3 * (v(1) + v(2)); break;
        default:
      break;
        }
      }
      //opserr << "\t\tSection("<<i<<") strain vector e=" << e(0) <<" " << e(1) << " " << e(2) << endln; 
      // Set the section deformations
      err += theSections[i]->setTrialSectionDeformation(e);
	  //Rslt = theSections[i]->getStressResultant();
	  //opserr << "\t\tSection("<<i<<") stress vector e="<<Rslt(0)<<" "<<Rslt(1)<<" "<<Rslt(2)<< endln; 
	  //if (e(1) != 0 && e(2) != 0 && Rslt(1) != 0 && Rslt(2) != 0) {
	  //  EI = Rslt(1)/e(1);
	  //  GA = Rslt(2)/e(2);
	  //  opserr << "\t\t Omega["<< i <<"]="<< "M/phi*gamma/V/L^2=" << Rslt(1) << "/"<< e(1) << "*" << e(2) << "/" << Rslt(2) << "=" << EI/GA/L/L << endln;
	  //} else {
		const Matrix &ks = theSections[i]->getSectionTangent();
        EI = ks(1,1);
		GA = ks(2,2);
		//opserr << "\t\t Omega["<< i <<"]="<< "EI/GA/L^2=" << EI << "/"<< GA << "=" << EI/GA/L/L << endln;
	  //}
	  temp += EI/GA/shearCF/L/L /numSections; //(1.-wts[i])/(numSections-1)wts[i]
	}

	error = abs(temp - Omega);
	//opserr << "\tmean Omega=" << temp << endln;
	Omega1 = Omega;
	Omega = temp;
	Omega2 = temp;

	iterNum ++; 

	if (iterNum > MAXITER) {
	  opserr<<"Timoshenko2d04::update()-iterNum > "<<MAXITER<< endln;
	  break;
	}

  } while (error > tol);

  // If Omega iterative doesnot converge, bsearch
  double temp1 = 0.0, temp2 = 0.0;
  if (iterNum > MAXITER) {
    if (Omega1 > Omega2) temp = Omega2; Omega2=Omega1; Omega1 = temp;
	while (Omega2-Omega1 > tol) {
	  temp1 = 0.; temp2 = 0.;
	  // for Omega1
	  for (int i = 0; i<numSections; i++) {
        mu    = 1./(1.+12.*Omega1);
        x     = L * pts[i];
        N1 = mu*(6.*x-4.*L*(1.+3.*Omega1))/L/L;
        N2 = 2.*mu*(3.*x+L*(6.*Omega1-1.))/L/L;
        N3 = 6.*mu*Omega1;
        for (int j = 0; j < order; j++) {
          switch(code(j)) {
          case SECTION_RESPONSE_P:     // axial strain
        e(j) = oneOverL*v(0); break;
          case SECTION_RESPONSE_MZ:    // curvature
        e(j) = N1 * v(1) + N2* v(2); break;
          case SECTION_RESPONSE_VY:    // shear strain
        e(j) = N3 * (v(1) + v(2)); break;
          default:
        break;
          }
        }
        err += theSections[i]->setTrialSectionDeformation(e);
        const Matrix &ks = theSections[i]->getSectionTangent();
        EI = ks(1,1);
        GA = ks(2,2);
  	    temp1 += EI/GA/shearCF/L/L /numSections; //(1.-wts[i])/(numSections-1)wts[i]
      }
	  // for x
	  double midOmega = (Omega1+Omega2)/2;
	  for (int i = 0; i<numSections; i++) {
        mu    = 1./(1.+12.*midOmega);
        x     = L * pts[i];
        N1 = mu*(6.*x-4.*L*(1.+3.*midOmega))/L/L;
        N2 = 2.*mu*(3.*x+L*(6.*midOmega-1.))/L/L;
        N3 = 6.*mu*midOmega;
        for (int j = 0; j < order; j++) {
          switch(code(j)) {
          case SECTION_RESPONSE_P:     // axial strain
        e(j) = oneOverL*v(0); break;
          case SECTION_RESPONSE_MZ:    // curvature
        e(j) = N1 * v(1) + N2* v(2); break;
          case SECTION_RESPONSE_VY:    // shear strain
        e(j) = N3 * (v(1) + v(2)); break;
          default:
        break;
          }
        }
        err += theSections[i]->setTrialSectionDeformation(e);
        const Matrix &ks = theSections[i]->getSectionTangent();
        EI = ks(1,1);
        GA = ks(2,2);
  	    temp2 += EI/GA/shearCF/L/L /numSections; //(1.-wts[i])/(numSections-1)wts[i]
      }
      if((temp1-Omega)*(temp2-Omega)<0) {
		Omega2= midOmega;
	    Omega = temp2;
	  } else {
		Omega1= midOmega;
		Omega = temp1;
	  }
    }
  }
  //opserr << "Final Omega: "<< Omega << endln;

  if (err != 0) {
    opserr << "Timoshenko2d04::update() - failed setTrialSectionDeformations(e)\n";
	return err;
  }

  return 0;
}

const Matrix&
Timoshenko2d04::getTangentStiff(void)
{
  static Matrix kb(3,3);

  // Zero for integral
  kb.Zero();
  q.Zero();
  
  double L = crdTransf->getInitialLength();
  const Vector &v = crdTransf->getBasicTrialDisp();
  double oneOverL = 1.0/L;
  
  double pts[maxNumSections];
  beamInt->getSectionLocations(numSections, L, pts);
  double wts[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wts);

  // Loop over the integration points
  for (int i = 0; i<numSections; i++) {
    int order = theSections[i]->getOrder(); // P M V
    const ID &code = theSections[i]->getType();

	// Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getSectionTangent();			
    const Vector &s = theSections[i]->getStressResultant();			
    
    // Perform numerical integration
	bd[i] = this->getBd(i, v, L);
	kb.addMatrixTripleProduct(1.0, bd[i], ks, L*wts[i]);
    q.addMatrixTransposeVector(1.0, bd[i], s, L*wts[i]);
  }

  // Add effects of element loads, q = q(v) + q0		
  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];

  // Transform to global stiffness
  K = crdTransf->getGlobalStiffMatrix(kb, q);

  return K;
}

const Matrix&
Timoshenko2d04::getInitialBasicStiff()
{
  static Matrix kb(3,3);

  // Zero for integral
  kb.Zero();
  
  double L = crdTransf->getInitialLength();
  const Vector &v = crdTransf->getBasicTrialDisp();
  double oneOverL = 1.0/L;
  
  double pts[maxNumSections];
  beamInt->getSectionLocations(numSections, L, pts);
  double wts[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wts);

  // Loop over the integration points
  for (int i = 0; i<numSections; i++) {
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
  
    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getInitialTangent();
    
	// Perform numerical integration
	bd[i] = this->getBd(i, v, L);
    kb.addMatrixTripleProduct(1.0, bd[i], ks, L*wts[i]);
  }
  return kb;
}

const Matrix&
Timoshenko2d04::getInitialStiff()
{
  const Matrix &kb = this->getInitialBasicStiff();

  // Transform to global stiffness
  K = crdTransf->getInitialGlobalStiffMatrix(kb);

  return K;
}

const Matrix&
Timoshenko2d04::getMass()
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
Timoshenko2d04::zeroLoad(void)
{
  Q.Zero();

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  return;
}

int 
Timoshenko2d04::addLoad(ElementalLoad *theLoad, double loadFactor)
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
    p0[2] -= V;

    // Fixed end forces in basic system
    q0[0] -= 0.5*P;
    q0[1] -= M;
    q0[2] += M;

  } else if (type == LOAD_TAG_Beam2dPointLoad) {
    double P = data(0)*loadFactor;
    double N = data(1)*loadFactor;
    double aOverL = data(2);					

	if (aOverL < 0.0 || aOverL > 1.0)
		return 0;

	double a = aOverL*L;
	double b = L-a;

	// Reactions in basic system
    p0[0] -= N;
    double V1 = P*(1.0-aOverL);
    double V2 = P*aOverL;
    p0[1] -= V1;
    p0[2] -= V2;

    double L2 = 1.0/(L*L);
    double a2 = a*a;
    double b2 = b*b;

    // Fixed end forces in basic system
    q0[0] -= N*aOverL;
    double M1 = -a * b2 * P * L2;
    double M2 = a2 * b * P * L2;
    q0[1] += M1;
    q0[2] += M2;

  } else {
    opserr << "Timoshenko2d04::Timoshenko2d04 -- load type unknown for element with tag: "
	   << this->getTag() << "Timoshenko2d04::addLoad()\n"; 
			    
    return -1;
  }

  return 0;
}

int 
Timoshenko2d04::addInertiaLoadToUnbalance(const Vector &accel)
{
	// Check for a quick return
	if (rho == 0.0) 
		return 0;

	// Get R * accel from the nodes
	const Vector &Raccel1 = theNodes[0]->getRV(accel);
	const Vector &Raccel2 = theNodes[1]->getRV(accel);

    if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
      opserr << "Timoshenko2d04::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
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
Timoshenko2d04::getResistingForce()
{
  crdTransf->update();  // Will remove once we clean up the corotational 3d transformation
  double L = crdTransf->getInitialLength();

  double pts[maxNumSections];
  beamInt->getSectionLocations(numSections, L, pts);
  double wts[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wts);
  
  // Zero for integration
  q.Zero();
  
  // Loop over the integration points
  for (int i = 0; i<numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
  
    double x = L * pts[i]; //xi;
	const Vector &v = crdTransf->getBasicTrialDisp();

    // Get section stress resultant
    const Vector &s = theSections[i]->getStressResultant();
    
    // Perform numerical integration on internal force
	bd[i] = this->getBd(i, v, L);
	q.addMatrixTransposeVector(1.0, bd[i], s, L*wts[i]);
  }

  // Add effects of element loads, q = q(v) + q0		
  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  
  // Vector for reactions in basic system
  Vector p0Vec(p0, 3);
  P = crdTransf->getGlobalResistingForce(q, p0Vec);	
  
  // Subtract other external nodal loads ... P_res = P_int - P_ext
  //P.addVector(1.0, Q, -1.0);
  P(0) -= Q(0);
  P(1) -= Q(1);
  P(2) -= Q(2);
  P(3) -= Q(3);
  P(4) -= Q(4);
  P(5) -= Q(5);

  return P;
}

const Vector&
Timoshenko2d04::getResistingForceIncInertia()
{
  // Add the inertial forces
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
    
  }
  
  // Add the damping forces
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) {
    P += this->getRayleighDampingForces();
  }

  return P;
}

int
Timoshenko2d04::sendSelf(int commitTag, Channel &theChannel)
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
    opserr << "Timoshenko2d04::sendSelf() - failed to send ID data\n";
     return -1;
  }    

  // send the coordinate transformation
  
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
     opserr << "Timoshenko2d04::sendSelf() - failed to send crdTranf\n";
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
    opserr << "Timoshenko2d04::sendSelf() - failed to send ID data\n";
    return -1;
  }    

  // send the sections
  
  for (j = 0; j<numSections; j++) {
    if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "Timoshenko2d04::sendSelf() - section " << 
	j << "failed to send itself\n";
      return -1;
    }
  }

  return 0;
}

int
Timoshenko2d04::recvSelf(int commitTag, Channel &theChannel,
			   FEM_ObjectBroker &theBroker)
{
  // receive the integer data containing tag, numSections and coord transformation info
 
	int dbTag = this->getDbTag();
  int i;
  
  static ID idData(7); // one bigger than needed so no clash with section ID

  if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
    opserr << "Timoshenko2d04::recvSelf() - failed to recv ID data\n";
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

      crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);

      if (crdTransf == 0) {
	opserr << "Timoshenko2d04::recvSelf() - failed to obtain a CrdTrans object with classTag " <<
	  crdTransfClassTag << endln;
	  return -2;	  
      }
  }

  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf object
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "Timoshenko2d04::sendSelf() - failed to recv crdTranf\n";
    return -3;
  }      
  
  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2*idData(3));
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "Timoshenko2d04::recvSelf() - failed to recv ID data\n";
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
      opserr << "Timoshenko2d04::recvSelf() - out of memory creating sections array of size " <<
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
	  opserr << "Timoshenko2d04::recvSelf() --default secTag at sec " << i+1 << endln;
	  theSections[i] = new FiberSection2d();

	  if (theSections[i] == 0) {
	opserr << "Timoshenko2d04::recvSelf() - Broker could not create Section of class type " <<
	  sectClassTag << endln;
	exit(-1);
      }
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "Timoshenko2d04::recvSelf() - section " << i << " failed to recv itself\n";
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
	opserr << "Timoshenko2d04::recvSelf() --default secTag at sec " << i+1 << endln;
	theSections[i] = new FiberSection2d();
	
	if (theSections[i] == 0) {
	opserr << "Timoshenko2d04::recvSelf() - Broker could not create Section of class type " <<
	  sectClassTag << endln;
	exit(-1);
	}
      }

      // recvSelf on it
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "Timoshenko2d04::recvSelf() - section " << i << " failed to recv itself\n";
	return -1;
      }     
    }
  }

  return 0;
}

void
Timoshenko2d04::Print(OPS_Stream &s, int flag)		
{
  s << "\nTimoshenko2d03, element id:  " << this->getTag() << endln;
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
Timoshenko2d04::displaySelf(Renderer &theViewer, int displayMode, float fact)
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
Timoshenko2d04::setResponse(const char **argv, int argc, OPS_Stream &s)
{
	Response *theResponse = 0;

    s.tag("ElementOutput");
    s.attr("eleType","Timoshenko2D");
    s.attr("eleTag",this->getTag());
    s.attr("node1",connectedExternalNodes[0]);
    s.attr("node2",connectedExternalNodes[1]);

    // global force - 
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
		|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {
		
	    s.tag("ResponseType","Px_1");
        s.tag("ResponseType","Py_1");
        s.tag("ResponseType","Mz_1");
        s.tag("ResponseType","Px_2");
        s.tag("ResponseType","Py_2");
        s.tag("ResponseType","Mz_2");

		theResponse = new ElementResponse(this, 1, P);

    // local force -
    } else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {
		
        s.tag("ResponseType","N_1");
        s.tag("ResponseType","V_1");
        s.tag("ResponseType","M_1");
        s.tag("ResponseType","N_2");
        s.tag("ResponseType","V_2");
        s.tag("ResponseType","M_2");

		theResponse = new ElementResponse(this, 2, P);

    // chord rotation -
	} else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0
		|| strcmp(argv[0],"basicDeformation") == 0) {
      
		 s.tag("ResponseType","eps");
         s.tag("ResponseType","theta_1");
         s.tag("ResponseType","theta_2");

		 theResponse = new ElementResponse(this, 3, Vector(3));
    
    // plastic rotation -
	} else if (strcmp(argv[0],"plasticRotation") == 0 || strcmp(argv[0],"plasticDeformation") == 0) {

        s.tag("ResponseType","epsP");
        s.tag("ResponseType","thetaP_1");
        s.tag("ResponseType","thetaP_2");

		theResponse = new ElementResponse(this, 4, Vector(3));
    

	} else if (strcmp(argv[0],"integrationPoints") == 0)
		theResponse = new ElementResponse(this, 10, Vector(numSections));

	else if (strcmp(argv[0],"integrationWeights") == 0)
		theResponse = new ElementResponse(this, 11, Vector(numSections));

    // section response -
	else if (strcmp(argv[0],"section") == 0 || strcmp(argv[0],"-section") == 0) {

      int sectionNum = atoi(argv[1]);
      if (sectionNum > 0 && sectionNum <= numSections)
	return theSections[sectionNum-1]->setResponse(&argv[2], argc-2, s);
      else if (sectionNum == 0) { // argv[1] was not an int, we want all sections, 
	CompositeResponse *theCResponse = new CompositeResponse();
	int numResponse = 0;
	double xi[maxNumSections];
	double L = crdTransf->getInitialLength();
	beamInt->getSectionLocations(numSections, L, xi);

	for (int i=0; i<numSections; i++) {

	  s.tag("GaussPoint");
	  s.attr("number",i+1);
	  s.attr("eta",xi[i]*L);

	  Response *theSectionResponse = theSections[i]->setResponse(&argv[1], argc-1, s);

	  if (theSectionResponse != 0) {
	    numResponse = theCResponse->addResponse(theSectionResponse);
	  }

	  s.endTag();
	}

	if (numResponse == 0) // no valid responses found
	  delete theCResponse;
	else
	  theResponse = theCResponse;

	  }
    }
    
    s.endTag(); // ElementOutput
    
    return theResponse;
}

int 
Timoshenko2d04::getResponse(int responseID, Information &eleInfo)		//LN modify
{

  if (responseID == 1) // global forces
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 2) { // local forces
      // Axial
	  P(3) =  q(0);
      P(0) = -q(0)+p0[0];
      // Moments about z and shears along y
	  P(2) = q(1);
      P(5) = q(2);
	  double L = crdTransf->getInitialLength();
      double V = (q(1)+q(2))/L;
      P(1) =  V + p0[1];
      P(4) = -V + p0[2];
      return eleInfo.setVector(P);
  }

  // Chord rotation
  else if (responseID == 3)
    return eleInfo.setVector(crdTransf->getBasicTrialDisp());

  // Plastic rotation
  else if (responseID == 4) {
    static Vector vp(3);
    static Vector ve(3);
    const Matrix &kb = this->getInitialBasicStiff();
    kb.Solve(q, ve);
    vp = crdTransf->getBasicTrialDisp();
    vp -= ve;
    return eleInfo.setVector(vp);
  }
  
  // integrationPoints
  else if (responseID == 10) {
    double L = crdTransf->getInitialLength();
    double pts[maxNumSections];
    beamInt->getSectionLocations(numSections, L, pts);
    Vector locs(numSections);
    for (int i = 0; i < numSections; i++)
      locs(i) = pts[i]*L;
    return eleInfo.setVector(locs);
  }

  //integrationWeight
  else if (responseID == 11) {
    double L = crdTransf->getInitialLength();
    double wts[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wts);
    Vector weights(numSections);
    for (int i = 0; i < numSections; i++)
      weights(i) = wts[i]*L;
    return eleInfo.setVector(weights);
  }
  else
    return -1;
}

Matrix
Timoshenko2d04::getNd(int sec, const Vector &v, double L)
{
  double pts[maxNumSections];
  beamInt->getSectionLocations(numSections, L, pts);

  double mu    = 1./(1.+12.*Omega);
  double x     = L * pts[sec];
  double phi1  =  mu*x*(L-x)*(L-x+6.*L*Omega)                      /L/L;
  //double phi1p =  mu*(3.*x*x+L*L*(1.+6.*Omega)-4.*L*(x+3.*x*Omega))/L/L;
  double phi2  = -mu*x*(L-x)*(x + 6.*L*Omega)                      /L/L;
  //double phi2p =  mu*(3.*x*x-L*L* 6. *Omega  +2.*L*x*(6.*Omega-1.))/L/L;
  double phi3  =  mu*(L-x)*(L-3.*x+12.*L*Omega)                    /L/L;
  //double phi3p =  mu*(6.*x - 4.*L*(1.+3.*Omega))                   /L/L;
  double phi4  =  mu*x*(  3.*x+2.*L*(6.*Omega-1.))                 /L/L;
  //double phi4p =  mu*2.*(3.*x+L*(6.*Omega-1.))                     /L/L;

  Matrix Nd(3,3);
  Nd.Zero();
  
  Nd(0,0) = 1./L;
  Nd(1,1) = phi3; // w, transverse displacement
  Nd(1,2) = phi4; // 
  Nd(2,1) = phi1; // theta, section rotation
  Nd(2,2) = phi2; // 
  
  return Nd;
}

Matrix
Timoshenko2d04::getBd(int sec, const Vector &v, double L)
{
  double pts[maxNumSections];
  beamInt->getSectionLocations(numSections, L, pts);
  
  double mu    = 1./(1.+12.*Omega);
  double x     = L * pts[sec];
  //double phi1  =  mu*x*(L-x)*(L-x+6.*L*Omega)                      /L/L;
  //double phi1p =  mu*(3.*x*x+L*L*(1.+6.*Omega)-4.*L*(x+3.*x*Omega))/L/L;
  //double phi2  = -mu*x*(L-x)*(x + 6.*L*Omega)                      /L/L;
  //double phi2p =  mu*(3.*x*x-L*L* 6. *Omega  +2.*L*x*(6.*Omega-1.))/L/L;
  //double phi3  =  mu*(L-x)*(L-3.*x+12.*L*Omega)                    /L/L;
  //double phi3p =  mu*(6.*x - 4.*L*(1.+3.*Omega))                   /L/L;
  //double phi4  =  mu*x*(  3.*x+2.*L*(6.*Omega-1.))                 /L/L;
  //double phi4p =  mu*2.*(3.*x+L*(6.*Omega-1.))                     /L/L;

  Matrix Bd(3,3);
  Bd.Zero();
  
  Bd(0,0) = 1./L;
  Bd(1,1) = mu*(6.*x-4.*L*(1.+3.*Omega))/L/L;      // sectional curvature
  Bd(1,2) = 2.*mu*(3.*x+L*(6.*Omega-1.))/L/L;      //
  Bd(2,1) = 6.*mu*Omega; //*2.0*(x/L) shear components 
  Bd(2,2) = 6.*mu*Omega; //*2.0*(1-x/L)

  return Bd;
}

// AddingSensitivity:BEGIN ///////////////////////////////////

const Matrix &
Timoshenko2d04::getKiSensitivity(int gradNumber)
{
	K.Zero();
	return K;
}

const Matrix &
Timoshenko2d04::getMassSensitivity(int gradNumber)
{
	K.Zero();
	return K;
}

const Vector &
Timoshenko2d04::getResistingForceSensitivity(int gradNumber)
{
	static Vector dummy(3);		// No distributed loads
	return dummy;
}

// NEW METHOD
int
Timoshenko2d04::commitSensitivity(int gradNumber, int numGrads)
{
 	return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////
