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
// $Date: 2007-09-21 15:28:38 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/RCFTAggregator.h,v $
                                                                        
                                                                        
// File: ~/section/SectionAggregator.h
//
// Written: MHS
// Created: Jun 2000
// Revision: A
//
// Description: This file contains the class definition for 
// SectionAggregator.  SectionAggregator decorates an MP
// section (couple bending and axial) with an uncoupled shear
// relation.
//
// What: "@(#) SectionAggregator.h, revA"

#ifndef RCFTAggregator_h
#define RCFTAggregator_h

#include <SectionForceDeformation.h>
#include <RCFTFiberSection3D.h>
#include <UniaxialMaterial.h>
#include <RCFTSlipMaterial.h>

#include <Vector.h>
#include <Matrix.h>

class RCFTAggregator : public SectionForceDeformation
{
  public:
    RCFTAggregator(); 

    RCFTAggregator(int tag, SectionForceDeformation &theSection,
		      int numAdditions, UniaxialMaterial **theAdditions,
		      const ID &code); 
    RCFTAggregator(int tag, int numAdditions,
		      UniaxialMaterial **theAdditions, const ID &code); 
    RCFTAggregator(int tag, SectionForceDeformation &thesection,
		      UniaxialMaterial &theAddition, int c);

    ~RCFTAggregator();

    int   setTrialSectionDeformation(const Vector &deforms); 
    const Vector &getSectionDeformation(void);

    const Vector &getStressResultant(void);
    const Matrix &getSectionTangent(void);
    const Matrix &getSectionCommitTangent(void);
    const Matrix &getInitialTangent(void);
    const Matrix &getSectionFlexibility(void);
    const Matrix &getInitialFlexibility(void);

    int   commitState(void);
    int   revertToLastCommit(void);    
    int   revertToStart(void);
 
    SectionForceDeformation *getCopy(void);
    const ID &getType (void);
    int getOrder (void) const;

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
 
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &output);
    int getResponse(int responseID, Information &info);

    int setVariable(const char *argv);
    int getVariable(int variableID, double &info);

  protected:
    
  private:
    
    SectionForceDeformation *theSection;
    UniaxialMaterial **theAdditions;

    ID *matCodes;
    int numMats;
    
    Vector *e;        // Storage for section deformations
    Vector *s;        // Storage for stress resultants
    Matrix *ks;       // Storage for section stiffness
    Matrix *ksCommit; // Storage for committed section stiffness
    Matrix *fs;       // Storage for section flexibility
    Matrix *fsCommit; // Storage for section commitflexibility
    ID     *theCode;  // Storage for section type information
   
    int otherDbTag;

    static double workArea[];
    static int codeArea[];
};

#endif
