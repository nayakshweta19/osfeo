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
**   Cenk Tort (tort0008@umn.edu)				      **
**   Jerome F. Hajjar (hajjar@struc.ce.umn.edu)                       **
**								      **
**   University of Minnesota - Twin Cities			      **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.1.1.1 $
// $Date: 2007-09-21 15:28:27 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/element/RCFTSTLBeamColumn3D.h,v $


// File: ~/model/element/RCFTBeamColumn3D.h

#ifndef RCFTSTLBeamColumn3D_h
#define RCFTSTLBeamColumn3D_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <SectionForceDeformation.h>
#include <RCFTSTLFiberSection3D.h>
#include <CrdTransf.h>
#include <RCFTSTLCrdTransf3D.h>
#include <BeamIntegration.h>

class Response;

class RCFTSTLBeamColumn3D: public Element
{
  public:
    RCFTSTLBeamColumn3D ();
    RCFTSTLBeamColumn3D (int tag, int nodeI, int nodeJ,
		    int numSections, RCFTSTLFiberSection3D *sectionPtrs[], BeamIntegration &bi,
		    RCFTSTLCrdTransf3D &coordTransf, double massDensPerUnitLength = 0.0);

    ~RCFTSTLBeamColumn3D();

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
    const Matrix &getMass(void);

    void zeroLoad(void);
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);

    bool isSubdomain(void);

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact);

    friend OPS_Stream &operator<<(OPS_Stream &s, RCFTSTLBeamColumn3D &E);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(char **argv, int argc, Information &eleInformation);
    int getResponse(int responseID, Information &eleInformation);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    const Vector &getBasicIncrDisp(void);
  protected:
    void setSectionPointers(int numSections, RCFTSTLFiberSection3D **secPtrs);
    int getInitialFlexibility(Matrix &fe);
  private:
    void initializeSectionHistoryVariables (void);
    Vector getLocalIncrDeltaDisp(void);
    Vector getBasicIncrDeltaDisp(void);
    double getDeformedLength(void);
    void calcDeformedLength(void);
    void calcResistingForce(void);
    
    // internal data
    ID     connectedExternalNodes; // tags of the end nodes
    int Tagg;

    BeamIntegration *beamIntegr;
    int numSections;
    RCFTSTLFiberSection3D **sections;          // array of pointers to sections
    //CrdTransf3d *crdTransf;        // pointer to coordinate tranformation object
    RCFTSTLCrdTransf3D *crdTransf;
    // (performs the transformation between the global and basic system)
    double rho;                    // mass density per unit length
    int    maxIters;               // maximum number of local iterations
    double tol;	                   // tolerance for relative energy norm for local iterations
    double deflength;

    int    initialFlag;            // indicates if the element has been initialized

    Node *theNodes[2];   // pointers to the nodes

    Matrix kv;                     // stiffness matrix in the basic system
    Vector Se;                     // element resisting forces in the basic system
    Vector Sg;                     // element resisting forces in the global system
    Vector Sglobal;

    Matrix kvcommit;               // commited stiffness matrix in the basic system
    Vector Secommit;               // commited element end forces in the basic system

    Matrix *ks;                    // array of section stiffness matrices
    Matrix *ksa;                    // array of section stiffness matrices

    Vector *dhat;
    Vector *DSQ;
    Vector *DSQa;
    Vector *f3;
    Vector *d3;
    Matrix *str_f3;
    Matrix *str_f3inv;

    Matrix *Ki;

    static Matrix theMatrix;
    static Vector theVector;
    static double workArea[];

    Vector XAxis;
    Vector YAxis;
    Vector ZAxis;

    double R[3][3];

    Matrix sr;
    Matrix ss;

    Vector df_i;
    Vector fk;
    Vector fk_incr;
    Vector ub;

    enum {maxNumSections = 10};
    // following are added for subdivision of displacement increment
    int    maxSubdivisions;       // maximum number of subdivisons of dv for local iterations
    //static int maxNumSections;
};

#endif


