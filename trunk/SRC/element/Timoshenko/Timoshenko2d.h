// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumnInt/Timoshenko2d.h,v $
// $Revision: 1.1 $
// $Date: 2009/01/10 21:22:20 $

// Created: 09/09
// Modified by: Li Ning 
// Description: This file contains the class implementation of TimoshenkoBeam.Based on Timoshenko2d.cpp.

  
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
			double rho = 0.0, double shearCF = 1.0);

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
	double Omega;   // shear contribution factor

	static Matrix *bd;
	static Matrix *nd;

	enum {maxNumSections = 20};
    static double workArea[];
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    // AddingSensitivity:END ///////////////////////////////////////////
};

#endif


