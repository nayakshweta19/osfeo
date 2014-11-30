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
// $Date: 2007-09-21 15:28:09 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/coordTransformation/RCFTMCrdLinTransf3D.h,v $
                                                                        
                                                                        
// File: ~/crdTransf/RCFTCrdTransf3d.h
//
//
// Description: This file contains the class definition for
// RCFTMCrdLinTransf3d.h. RCFTMCrdLinTransf3d provides the
// abstraction of a linear transformation for a spatial frame
// between the global and basic coordinate systems

// What: "@(#) RCFTMCrdLinTransf3d.h, revA"

#ifndef RCFTMCrdLinTransf3D_h
#define RCFTMCrdLinTransf3D_h

#include <CrdTransf.h>
#include <Vector.h>
#include <Matrix.h>

class RCFTMCrdLinTransf3D: public CrdTransf
{
  public:
    RCFTMCrdLinTransf3D (int tag, const Vector &vecInLocXZPlane);
    //RCFTCrdTransf3D (int tag, const Vector &vecInLocXZPlane);
    
    RCFTMCrdLinTransf3D();
    ~RCFTMCrdLinTransf3D();

    int    initialize(Node *node1Pointer, Node *node2Pointer);
    int    update(void);
    double getInitialLength(void);
    double getDeformedLength(void);

    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);
    
	const Vector &getBasicTrialDisp     (void);
	const Vector &getBasicIncrDisp      (void);
	const Vector &getBasicIncrDeltaDisp (void);
	const Vector &getBasicTrialVel      (void);
	const Vector &getBasicTrialAccel    (void);

    const Vector &getGlobalResistingForce (const Vector &basicForce, const Vector &ulocal);
    const Matrix &getLocalStiffMatrix     (const Matrix &basicStiff, const Vector &basicForce);
    const Matrix &getGlobalStiffMatrix    (const Matrix &basicStiff, const Vector &basicForce);
    const Matrix &getInitialGlobalStiffMatrix(const Matrix &kb);

    CrdTransf *getCopy3d(void);
    
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);
  
    // method used to rotate consistent mass matrix
    const Matrix &getGlobalMatrixFromLocal(const Matrix &local);

    // methods used in post-processing only
    const Vector &getPointGlobalCoordFromLocal(const Vector &localCoords);
    const Vector &getPointGlobalDisplFromBasic(double xi, const Vector &basicDisps);

    int  getLocalAxes(Vector &xAxis, Vector &yAxis, Vector &zAxis);
  
  private:
    int  computeElemtLengthAndOrient(void);
    void compTransfMatrixLocalGlobal(Matrix &Tlg);
    int  getLocalAxes (void);

    Matrix kl;
    Vector pg;

    // internal data
    Node *nodeIPtr, *nodeJPtr;          // pointers to the element two endnodes

    double R[3][3];	// Transformation matrix
    double L;        // undeformed element length

    static Matrix Tlg;  // matrix that transforms from global to local coordinates
    static Matrix kg;   // global stiffness matrix

};

#endif

