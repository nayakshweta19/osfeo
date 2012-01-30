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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/element/RCFTSCHBeamColumn3D.h,v $


// File: ~/model/element/RCFTBeamColumn3D.h

#ifndef RCFTSCHBeamColumn3D_h
#define RCFTSCHBeamColumn3D_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <SectionForceDeformation.h>
#include <RCFTFiberSection3D.h>
#include <RCFTAggregator.h>
#include <CrdTransf.h>
#include <RCFTCrdTransf3D.h>
#include <BeamIntegration.h>

class Response;

class RCFTSCHBeamColumn3D: public Element
{
  public:
    RCFTSCHBeamColumn3D ();
    RCFTSCHBeamColumn3D (int tag, int nodeI, int nodeJ,
		    int numSections, RCFTAggregator *sectionPtrs[], BeamIntegration &bi,
		    RCFTCrdTransf3D &coordTransf, double massDensPerUnitLength = 0.0);

    ~RCFTSCHBeamColumn3D();

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

    friend OPS_Stream &operator<<(OPS_Stream &s, RCFTSCHBeamColumn3D &E);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(char **argv, int argc, Information &eleInformation);
    int getResponse(int responseID, Information &eleInformation);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    const Vector &getBasicIncrDisp(void);
  protected:
    void setSectionPointers(int numSections, RCFTAggregator **secPtrs);
    int getInitialFlexibility(Matrix &fe);
  private:
    void initializeSectionHistoryVariables (void);
    Vector getLocalIncrDeltaDisp(void);
    double getDeformedLength(void);
    void calcDeformedLength(void);
    void calcResistingForce(void);    

    // internal data
    ID     connectedExternalNodes; // tags of the end nodes
    int itr;

    BeamIntegration *beamIntegr;
    int numSections;
    RCFTAggregator **sections;          // array of pointers to sections
    //CrdTransf3d *crdTransf;        // pointer to coordinate tranformation object
    RCFTCrdTransf3D *crdTransf;
    // (performs the transformation between the global and basic system)
    double rho;                    // mass density per unit length
    int    maxIters;               // maximum number of local iterations
    double tol;	                   // tolerance for relative energy norm for local iterations
    double deflength;
    double Li;

    int    initialFlag;            // indicates if the element has been initialized

    Node *theNodes[2];   // pointers to the nodes

    Matrix kv;                     // stiffness matrix in the basic system
    Vector Sg;                     // element resisting forces in the global system
    Vector Sglobal;

    Matrix kvcommit;               // commited stiffness matrix in the basic system

    Matrix *ksa;                    // array of section stiffness matrices

    Vector *dhat;
    Vector *DSQa;
    Vector *CDSQa;
    Vector *DSQ;
    Vector *f4;
    Vector *d4;
    Matrix *str_f4;
    Matrix *str_f4inv;
    double slp_strn;

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
    Vector f_i;
    Vector fk;
    Vector ub;

    double p_si;
    double p_ci;
    double my_si;
    double mz_si;
    double my_ci;
    double mz_ci;

    double p_sj;
    double p_cj;
    double my_sj;
    double mz_sj;
    double my_cj;
    double mz_cj;

    int Tagg;

    enum {maxNumSections = 10};
    // following are added for subdivision of displacement increment
    int    maxSubdivisions;       // maximum number of subdivisons of dv for local iterations
    //static int maxNumSections;
};

#endif


