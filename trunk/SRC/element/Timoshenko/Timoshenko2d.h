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

// $Source: /usr/local/cvs/OpenSees/SRC/element/Timoshenko/Timoshenko2d.cpp,v $
// $Revision: 1.2 $
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
  
#ifndef Timoshenko2d_h
#define Timoshenko2d_h

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <BeamIntegration.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>

class Node;
class SectionForceDeformation;
class CrdTransf;
class TimoshenkoLinearCrdTransf2d;
class Response;

class Timoshenko2d : public Element
{
  public:
    Timoshenko2d(int tag, 
			int nd1, 
			int nd2,
			int numSec, 
			SectionForceDeformation **s,
			CrdTransf &coordTransf, 
			BeamIntegration &bi,
			double rho = 0.0, double shearCF = 1.0, int noIter = 0);

    Timoshenko2d();
    virtual ~Timoshenko2d();

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    int update(void);
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInfo);


    // AddingSensitivity:BEGIN //////////////////////////////////////////
    const Vector & getResistingForceSensitivity(int gradNumber);
    const Matrix & getKiSensitivity(int gradNumber);
    const Matrix & getMassSensitivity(int gradNumber);
    int            commitSensitivity(int gradNumber, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    const Matrix &getInitialBasicStiff(void);
	Matrix getNd(int sec, const Vector &v, double L);
	Matrix getBd(int sec, const Vector &v, double L);

    int numSections;

    SectionForceDeformation **theSections; // pointer to the ND material objects
    CrdTransf *crdTransf;          // pointer to coordinate transformation object 
    BeamIntegration *beamInt;
	
	ID connectedExternalNodes; // Tags of quad nodes
    
	Node *theNodes[2];
    
	static Matrix K;		// Element stiffness, damping, and mass Matrix
    static Vector P;		// Element resisting force vector

    Vector Q;		// Applied nodal loads
    Vector q;		// Basic force
    double q0[3];   // Fixed end forces in basic system
    double p0[3];   // Reactions in the basic system due to element loads

    double rho;	    // Mass density per unit length
	double shearCF; // shear corrector factor
	int    iterSwitch; // Iteration switch of the IIE shape function
	//double Omega;   // shear contribution factor
	Vector Rslt;
	Vector Defo;

	static Matrix *bd;
	static Matrix *nd;

	enum {maxNumSections = 20};
	double Omega[maxNumSections];
    static double workArea[];
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    // AddingSensitivity:END ///////////////////////////////////////////
};

#endif


