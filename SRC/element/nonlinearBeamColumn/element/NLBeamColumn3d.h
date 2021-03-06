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
// $Date: 2007-09-21 15:28:32 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/element/nonlinearBeamColumn/element/NLBeamColumn3d.h,v $
                                                                        
                                                                        
// File: ~/model/element/NLBeamColumn3d.h
//
// Written: Remo Magalhaes de Souza on 03/99 
// Revised: rms 06/99 (mass matrix)
//          rms 07/99 (using setDomain)
//          rms 08/99 (included P-Delta effect)
//	    fmk 10/99 setResponse() & getResponse()
//          rms 04/00 (using CrdTransf class)
//          mhs 06/00 (using new section class)
//          mhs 06/00 (using new section class w/ variable dimensions)
//          rms 06/00 (torsional stiffness considered at the section level)
//          rms 06/00 (making copy of the sections)
//          rms 06/00 (storing section history variables at the element level)
//            
// Purpose: This file contains the class definition for NLBeamColumn3d.
// NLBeamColumn3d is a materially nonlinear flexibility based frame element.

#ifndef NLBeamColumn3d_h
#define NLBeamColumn3d_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <SectionForceDeformation.h>
#include <CrdTransf3d.h>
#include <GaussLobattoQuadRule1d01.h>

class Response;

class NLBeamColumn3d: public Element
{
  public:
    NLBeamColumn3d ();
    NLBeamColumn3d (int tag, int nodeI, int nodeJ,  
                    int numSections, SectionForceDeformation *sectionPtrs[], 
                    CrdTransf3d &coordTransf, double massDensPerUnitLength = 0.0, 
		    int maxNumIters = 10, double tolerance = 1e-12);
    
    ~NLBeamColumn3d();

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
    
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact);

    
    friend OPS_Stream &operator<<(OPS_Stream &s, NLBeamColumn3d &E);        
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, Information &eleInformation);
    int getResponse(int responseID, Information &eleInformation);
    
    int setParameter(const char **argv, int argc, Information &info);
    int updateParameter(int parameterID, Information &info);

  private:
    void getGlobalDispls(Vector &dg) const;
    void getGlobalAccels(Vector &ag) const;             
    void getForceInterpolatMatrix(double xi, Matrix &b, const ID &code);
    void getDistrLoadInterpolatMatrix(double xi, Matrix &bp, const ID &code);
    void compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const;
    void initializeSectionHistoryVariables (void);
    
    
    // internal data
	
    ID     connectedExternalNodes; // tags of the end nodes
    int    nSections;              // number of sections (integration 
                                   // points) along the element
    SectionForceDeformation **sections;          // array of pointers to sections
    CrdTransf3d *crdTransf;        // pointer to coordinate tranformation object 
	                           // (performs the transformation between the global and basic system)
    Node *theNodes[2];

    double rho;                    // mass density per unit length
    int    maxIters;               // maximum number of local iterations
    double tol;	                   // tolerance for relative energy norm for local iterations

    int    initialFlag;            // indicates if the element has been initialized
    bool isTorsion;
	
    Vector load;                   // equivalent nodal loads ????

    Matrix kv;                     // stiffness matrix in the basic system 
    Vector Se;                     // element resisting forces in the basic system

    Matrix kvcommit;               // commited stiffness matrix in the basic system
    Vector Secommit;               // commited element end forces in the basic system

    Matrix *fs;                    // array of section flexibility matrices
    Vector *vs;                    // array of section deformation vectors
    Vector *Ssr;                   // array of section resisting force vectors
 
    Vector *vscommit;              // array of commited section deformation vectors

    Matrix *sp;  // Applied section forces due to element loads, 5 x nSections
    double p0[5]; // Reactions in the basic system due to element loads
    Matrix *Ki;

    static Matrix theMatrix;
    static Vector theVector;
    static GaussLobattoQuadRule1d01 quadRule;
    static double workArea[];
};

#endif

