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
// $Date: 2007-09-21 15:28:43 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/uniaxial/RCFTSlipMaterial.h,v $
                                                                        
#ifndef RCFTSlipMaterial_h
#define RCFTSlipMaterial_h

// File: ~/material/ElasticPPMaterial.h
//
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for 
// ElasticPPMaterial. ElasticPPMaterial provides the abstraction
// of an elastic perfectly plastic uniaxial material, 
//
// What: "@(#) ElasticPPMaterial.h, revA"

#include <UniaxialMaterial.h>

class RCFTSlipMaterial : public UniaxialMaterial
{
  public:
    RCFTSlipMaterial(int tag, double fy, double e, double esh);    
    RCFTSlipMaterial();    

    ~RCFTSlipMaterial();

    int setTrialStrain(double strain, double strainRate); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
    double getCommitTangent(void);

    double getInitialTangent(void) {return E;};

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);    

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
  protected:
    
  private:
    double fyp, fyn;	// positive and negative yield stress
    double E,Esh,Fy;		// elastic modulus
    double trialStrain;	// trial strain
    double commitStrain;
    double commitTangent; 
    double commitStress;
    double commitX, trialX;

    double trialStress;      // current trial stress
    double trialTangent;     // current trial tangent
};


#endif



