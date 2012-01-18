/* ****************************************************************** **
** OpenSees - Open System for Earthquake Engineering Simulation       **
** Pacific Earthquake Engineering Research Center                     **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited. See    **
** file 'COPYRIGHT' in main directory for information on usage and    **
** redistribution, and for a DISCLAIMER OF ALL WARRANTIES.            **
**                                                                    **
** Developed by:                                                      **
** Frank McKenna (fmckenna@ce.berkeley.edu)                           **
** Gregory L. Fenves (fenves@ce.berkeley.edu)                         **
** Filip C. Filippou (filippou@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.00 $
// $Date: 2008/07/18 18:05:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/Pipe/Pipelin2.cpp,v $


// Written: Antonios Vytiniotis
// Created: 07/08
// Revision: A
//
// Description: This file contains the implementation for the Pipelin2 class.
//

#include "Pipelin2.h"
#include <Information.h>
#include <Parameter.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <ElementResponse.h>

#include <Matrix.h>
#include <Vector.h>

#include <ElasticMaterial.h>

// initial the class wide variables
Matrix Pipelin2::trussK(6,6);
Matrix Pipelin2::trussM(6,6);
Matrix Pipelin2::trussD(6,6);
Vector Pipelin2::trussR(6);

Pipelin2::Pipelin2(int tag,
	int Nd1, int Nd2,
	UniaxialMaterial &theMat,
	double a, double c3, double g)
	:Element(tag,ELE_TAG_Pipelin2),
	theMaterial(0),
	externalNodes(2),
	trans(1,4),L(0.0), A(a), C_3(c3), Gamma(g),
	end1Ptr(0), end2Ptr(0), eletag(tag),
	d_y_class(0.0)
{
	//create a copy of the material object
	theMaterial=theMat.getCopy();

	//fill in the ID containing external node info with node id's
	externalNodes(0)=Nd1;
	externalNodes(1)=Nd2;
	for (int i=0; i<2; i++)
		theNodes[i] = 0;
	trussR.Zero();
}

//constructor which should be invoked by an FE_ObjectBroker only
Pipelin2::Pipelin2()
	:Element(0,ELE_TAG_Pipelin2),
	theMaterial(0),
	externalNodes(2),
	trans(1,4), L(0.0), A(0.0), C_3(0.0), Gamma(0.0),
	end1Ptr(0), end2Ptr(0),
	d_y_class(0.0)
{
	for (int i=0; i<2; i++)
		theNodes[i] = 0;
}

Pipelin2::~Pipelin2()
{
	if (theMaterial !=0)
		delete theMaterial;
}

int Pipelin2::getNumExternalNodes(void) const
{
	return 2;
}

const ID &
	Pipelin2::getExternalNodes(void)
{
	return externalNodes;
}

int
	Pipelin2::getNumDOF(void){
		return 6;
}

Node **
	Pipelin2::getNodePtrs(void)
{
	return theNodes;
}

void
	Pipelin2::setDomain(Domain *theDomain)
{
	//first ensure nodes exist in Domain and set the node pointers
	int Nd1 =externalNodes(0);
	int Nd2 =externalNodes(1);
	end1Ptr =theDomain->getNode(Nd1);
	end2Ptr =theDomain->getNode(Nd2);
	theNodes[0] = theDomain->getNode(Nd1);
	theNodes[1] = theDomain->getNode(Nd2);

	if (theNodes[0]==0)
		return;
	if (theNodes[1]==0)
		return;

	// call the DomainComponent class method
	this->DomainComponent::setDomain(theDomain);

	//ensure connected nodes have corrent number of dof's
	int dofNd1=theNodes[0]->getNumberDOF();
	int dofNd2=theNodes[1]->getNumberDOF();
	if ((dofNd1 !=3) || (dofNd2 !=3))
		return; //don't go any further otherwise segmentation fault

	//now determine the length & transformation matrix
	const Vector &end1Crd=theNodes[0]->getCrds();
	const Vector &end2Crd=theNodes[1]->getCrds();

	double dx= end2Crd(0)-end1Crd(0);
	double dy= end2Crd(1)-end1Crd(1);

	d_y_class=dy;
	L=sqrt(dx*dx+dy*dy);

	if (L==0.0)
		return;

	double cs=dx/L;
	double sn=dy/L;
	trans(0,0)=-cs;
	trans(0,1)=-sn;
	trans(0,2)= cs;
	trans(0,3)= sn;
//	// determine the nodal mass for lumped mass approach
//	M=M*A*L/2; //M was set to rho by the constructor
}

int
	Pipelin2::commitState()
{
	return theMaterial->commitState();
}

int
	Pipelin2::revertToLastCommit()
{
	return theMaterial->revertToLastCommit();
}

int
	Pipelin2::revertToStart()
{
	return theMaterial->revertToStart();
}

int
	Pipelin2::update()
{
	//determine the current strain given trial displacements at nodes
	double strain=this->computeCurrentStrain();

	//set the strain in the materials
	theMaterial->setTrialStrain(strain);

	return 0;
}

const Matrix &
	Pipelin2::getTangentStiff(void)
{
	if (L==0) {//if length ==zero - we zero and return
		trussK.Zero();
		return trussK;
	}

	//get the current E from the material for the strain that was set
	// at the material when the update() method was invoked
	
	double E = theMaterial->getTangent();
	
	//form the tangent stiffness matrix
	Matrix K_temp(4,4);
	K_temp =trans^trans; //This is a temporary matrix containing the truss stiffness parameters
	K_temp *=A*E/L;

//	trussK.Zero();
	// Truss stiffness components:
	trussK(0,0)=K_temp(0,0);
	trussK(1,0)=K_temp(1,0);
	trussK(0,1)=K_temp(0,1);
	trussK(1,1)=K_temp(1,1);
	trussK(3,3)=K_temp(2,2);
	trussK(4,3)=K_temp(3,2);
	trussK(3,4)=K_temp(2,3);
	trussK(4,4)=K_temp(3,3);
	trussK(2,2)=0.;
	trussK(5,5)=0.;
	trussK(5,2)=0.;
	trussK(2,5)=0.;

	return trussK;
}

const Matrix &
	Pipelin2::getInitialStiff(void)
{
	if (L==0) {
		trussK.Zero();
		return trussK;
	}

	//get the current strain from the material
	double strain = theMaterial->getStrain();

	//get the current stress from the material
	double stress = theMaterial->getStress();

	//compute the tangent
	double E=stress/strain;

	//form the tangent stiffness matrix
	Matrix K_temp(4,4);
	K_temp =trans^trans; //This is a temporary matrix containing the truss stiffness parameters
	K_temp *=A*E/L;

//	trussK.Zero();
	// Truss stiffness components:
	trussK(0,0)=K_temp(0,0);
	trussK(1,0)=K_temp(1,0);
	trussK(0,1)=K_temp(0,1);
	trussK(1,1)=K_temp(1,1);
	trussK(3,3)=K_temp(2,2);
	trussK(4,3)=K_temp(3,2);
	trussK(3,4)=K_temp(2,3);
	trussK(4,4)=K_temp(3,3);
	trussK(2,2)=0.;
	trussK(5,5)=0.;
	trussK(5,2)=0.;
	trussK(2,5)=0.;

//	opserr << "Componenents of the matrix K11" << trussK(3,3) << endln;
	return trussK;
}

const Matrix &
	Pipelin2::getDamp(void)
{
	//No damping associated with this type of element
	trussD.Zero();
	trussD(2,2)=-C_3/L;
	trussD(5,5)=-C_3/L;
	trussD(5,2)=C_3/L;
	trussD(2,5)=C_3/L;
	double deleteme=C_3/L;
	double deleteme2=C_3/L;
	return trussD;
}

const Matrix &
	Pipelin2::getMass(void)
{
	if (L==0){
		trussM.Zero();
		return trussM;
	}
//	At this point we have zero lumped mass
	trussM.Zero();
	return trussM;
}

void
	Pipelin2::zeroLoad(void)
{
	//does nothing - no element load associated with this object
}

int
	Pipelin2::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	opserr <<"MyTruss::addLoad - load type unknown for truss with tag: " << this->getTag() << endln;
	 return -1;
}

int
	Pipelin2::addInertiaLoadToUnbalance(const Vector &accel)
{
	return 0;
}

const Vector &
	Pipelin2::getResistingForce()
{
	if (L==0) {//if length ==zero - zero and return
		trussR.Zero();
		return trussR;
	}
	// R=Ku-Pext
	//force =F*transformation
	double force = A* theMaterial->getStress();
	trussR(0)= trans(0,0)*force;
	trussR(1)= trans(0,1)*force;
	trussR(3)= trans(0,2)*force;
	trussR(4)= trans(0,3)*force;

	const Vector &vel1 = theNodes[0]->getTrialVel();
	const Vector &vel2 = theNodes[1]->getTrialVel();

	// This is the linear element with total disp (no need for state params)
	// Domain::update(double a, double b);
	// double dt=;
	// trussR(2)=trussD(2,2)*vel1(2)+trussD(2,5)*vel2(2)+trussD(2,2)*d_y_class*Gamma;
	// trussR(5)=trussD(5,2)*vel1(2)+trussD(5,5)*vel2(2)-trussD(5,5)*d_y_class*Gamma;
	trussR(2)=-C_3/L*vel1(2)+C_3/L*vel2(2)-C_3/L*d_y_class*Gamma;
	trussR(5)=C_3/L*vel1(2)-	C_3/L*vel2(2)+C_3/L*d_y_class*Gamma;
	return trussR;
}

const Vector &
	Pipelin2::getResistingForceIncInertia()
{
	this->getResistingForce();

	//No inertia is included in the in this element formulation
	return trussR;
}

int
	Pipelin2::sendSelf (int commitTag, Channel &theChannel)
{
	int dataTag=this->getDbTag();
	// Pipelin2 packs it's data into a Vector and sends this to theChannel
	//along with it's dbTag and the commitTag passed in the arguments
	Vector data(6);
	data(0)= this->getTag();
	data(1)=A;
	data(4)=C_3;
	data(5)=Gamma;
	data(2)=theMaterial->getClassTag();
	int matDbTag=theMaterial->getDbTag();
	if (matDbTag==0) {
		matDbTag =theChannel.getDbTag();
		if (matDbTag !=0)
			theMaterial->setDbTag(matDbTag);
	}
	data(3)=matDbTag;

	theChannel.sendVector (dataTag, commitTag, data);

	theChannel.sendID(dataTag, commitTag, externalNodes);

	theMaterial->sendSelf(commitTag, theChannel);

	return 0;
}

int
	Pipelin2::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	int dataTag= this->getDbTag();
	Vector data(6);
	theChannel.recvVector(dataTag, commitTag, data);

	this->setTag((int)data(0));
	A=data(1);
	C_3=data(4);
	Gamma=data(5);
	theChannel.recvID(dataTag, commitTag, externalNodes);

	int matClass=data(2);
	int matDb = data(3);
	theMaterial= theBroker.getNewUniaxialMaterial(matClass);

	theMaterial->setDbTag(matDb);
	theMaterial->recvSelf(commitTag, theChannel, theBroker);

	return 0;
}

int
	Pipelin2::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	const Vector &end1Crd= end1Ptr->getCrds();
	const Vector &end2Crd= end2Ptr->getCrds();
	const Vector &end1Disp=end1Ptr->getDisp();
	const Vector &end2Disp=end2Ptr->getDisp();

	Vector v1(3);
	Vector v2(3);
	for (int i=0; i<2;i++) {
		v1(i)=end1Crd(i)+end1Disp(i)*fact;
		v2(i)=end2Crd(i)+end2Disp(i)*fact;
	}
	if (displayMode==3) {
		//use the strain as the drawing measure
		double strain = theMaterial->getStrain();
		return theViewer.drawLine(v1, v2, strain, strain);
	}
	else if (displayMode==2){
		//otherwise use the material stress
		double stress =A*theMaterial->getStress();
		return theViewer.drawLine(v1,v2,stress, stress);
	}
	else{
		//use the axial force
		double force = A*theMaterial->getStress();
		return theViewer.drawLine(v1,v2,force,force);
	}
}

void
	Pipelin2::Print(OPS_Stream &s, int flag)
{
	//compute the strain and axial force in the member
	double strain, force;
	if (L==0.0) {
		strain=0;
		force=0.0;
	}
	else{
		strain = theMaterial->getStrain();
		force=A*theMaterial->getStress();
	}
	trussR(0)= trans(0,0)*force;
	trussR(1)= trans(0,1)*force;
	trussR(3)= trans(0,2)*force;
	trussR(4)= trans(0,3)*force;

	const Vector &vel1 = theNodes[0]->getVel();
	const Vector &vel2 = theNodes[1]->getVel();
	trussR(2)=-C_3/L*vel1(2)+C_3/L*vel2(2)-C_3/L*d_y_class*Gamma;
	trussR(5)=C_3/L*vel1(2)-C_3/L*vel2(2)+C_3/L*d_y_class*Gamma;

	if (flag==0) {//print everything
		s<< "Element: " <<this->getTag();
		s<< " type: My Truss iNode: "<< externalNodes(0);
		s<< " jNode: "<<externalNodes(1);
		s<< " Area: "<< A;
		if (Gamma!=0) s << "Gamma: "<<Gamma;

		s<< " \n\t strain: " <<strain;
		s<< " axial load: " <<force;
		s<< " \n\t unbalanced load: " <<trussR;
		s<< " \t Material: " << *theMaterial;
		s<< endln;
	} else if (flag==1) {//just print ele id, strain and force
		s<< this->getTag() << " " <<strain << " " << force <<endln;
	}
}

Response *
	Pipelin2::setResponse(const char **argv, int argc , OPS_Stream &s)
{
	// we compare arg(0) for known response types for the Truss
	// axial force
	if(strcmp(argv[0], "axialForce")==0)
		return new ElementResponse(this, 1, 0.0);

	//a material quantity
	else if (strcmp(argv[0], "material")==0)
		return theMaterial->setResponse(&argv[1], argc-1, s);
	else
		return 0;
}

int
	Pipelin2::getResponse(int responseID, Information &eleInformation)
{
	switch (responseID){
	case -1:
		return -1;

	case 1:
		return eleInformation.setDouble (A*theMaterial->getStress());

	default:
		return 0;
	}
}

double
	Pipelin2::computeCurrentStrain(void) const
{
	//determine the strain
	const Vector &disp1=end1Ptr->getTrialDisp();
	const Vector &disp2=end2Ptr->getTrialDisp();

	double dLength=0.0;
	for (int i=0;i<2;i++)
		dLength -= (disp2(i)-disp1(i))*trans(0,i);

	double strain =dLength/L;

	return strain;
}

