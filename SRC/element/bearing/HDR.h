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

// Written: Manish Kumar (mkumar2@buffalo.edu)
// The high damping rubber bearing model of Grant et. al (2004) is implemented here
// Created: 07/14/2013
// Revision: A

#ifndef HDR_h
#define HDR_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class UniaxialMaterial;
class Response;

class HDR : public Element
{
public:
    // Constructor
    HDR(int tag, int Nd1, int Nd2, double qRubber, double uh, double Gr, double Kbulk, double D1, double D2, double ts, double tr, int n, 
		double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3, double c4, 
		const Vector y, const Vector x=0, double kc=10, double PhiM=0.75, double ac=1.0, double sDratio=0.5, double m=0.0, double cd=128000, double tc=0.0);
        
    HDR();
       
     // Destructor
    ~HDR();
       
    // Method to get class type
    const char *getClassType() const {return "HDR";};
   
    // Public methods to obtain information about dof & connectivity    
    int getNumExternalNodes() const;
    const ID &getExternalNodes();
    Node **getNodePtrs();
    int getNumDOF();
    void setDomain(Domain *theDomain);
       
    // Public methods to set the state of the element    
    int commitState();
    int revertToLastCommit();        
    int revertToStart();        
    int update();

    // Public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    const Matrix &getMass();
       
    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);
   
    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();
   
    // Public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact);    
    void Print(OPS_Stream &s, int flag = 0);    
       
    // Public methods for element recorder
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInfo);
   
protected:

private:
    // Private methods
    void setUp();
    double sgn(double x);
   
    // Private attributes - a copy for each object of the class
    ID connectedExternalNodes;					// contains the tags of the end nodes
    Node *theNodes[2];							// array of nodes
        UniaxialMaterial *theMaterials[3];		// array of uniaxial materials
   
    // PARAMETERS
	// Horizontal direction
	double a1, a2, a3, b1, b2, b3, c1, c2, c3, c4;
    double k0;									// Initial stiffness of hysteretic component(due to lead)
	double qYield;								// Current yield stress
    double ke;									// Stiffness of elastic component(due to rubber)
	double cd;									// Viscous damping parameter

	// Vertical direction
	double Ec;
	double Kv0;									// Stiffness at zero horizontal displacement
	double Kv;									// Elastic stiffness in compression and tension
	double kc;									// Quality index of elastomer
	double PhiM;								// Maximum reduction in the cavitation strength of elastomer
	double ac;									// Strength degradation parameter
	double Fcr;									// Critical buckling load at zero lateral deformation
	double ucr;									// Critical buckling deformation at zero lateral deformation
	double Fc;									// Initial cavitation strength
	double uc;									// Initial cavitation deformation

	double Kt, Kr;								// Torsional and rotational stiffness of bearing


	// Others
	double G;									// Shear modulus of rubber
	Vector x;									// local x direction
    Vector y;									// local y direction
    double shearDistI;							// shear distance from node I as fraction of length
    double mass;								// mass of element
	double Tr;									// height of rubber in the bearing
	double D1, D2;								// Inner and outer diameter of the bearing
    double L;									// element length
	double h;									// height of rubber + shims
	double rg;									// radius of gyration of bearing
	double A;									// Bonded rubber area of bearing
	double Ar;									// reduced bonded area due to shear displacement
	double n, ts;								// number of layers and shim thickness

    // State variables
	double Fcrn;								// Current critical buckling force at particular shear deformation
	double Fcn;									// Current cavitation strength and deformation
	double umax;								// Maximum force and deformation ever experienced by the elastomer
	double DSplus, DSminus, DS;
	double DM, Delta;
	Vector F2;


    Vector ub;									// Displacements in basic system
	Vector ubdot;								// Velocities in basic system
    Vector qb;									// Forces in basic system
    Matrix kb;									// Stiffness matrix in basic system
    Vector ul;									// Displacements in local system
    Matrix Tgl;									// Transformation matrix from global to local system
    Matrix Tlb;									// Transformation matrix from local to basic system
   
    // Committed history variables
	Vector ubC;									// Displacements in basic system
	double DSplusC, DSminusC, DSC;
	double DMC;
	Vector F2C;
    // Initial stiffness matrix in basic system
    Matrix kbInit;
   
	

    static Matrix theMatrix;
    static Vector theVector;
    static Vector theLoad;
};

#endif