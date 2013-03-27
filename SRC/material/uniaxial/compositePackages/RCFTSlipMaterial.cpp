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
                                                                        
// $Revision: 1.2 $
// $Date: 2008-07-03 18:08:38 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/uniaxial/RCFTSlipMaterial.cpp,v $
                                                                        
                                                                        
// File: ~/material/RCFTSlipMaterial.cpp
//
// Written: Cenk Tort
// Created: 03/05
// Revision: A
//
// Description: This file contains the class implementation for 
// RCFTSlipMaterial. 
//
// What: "@(#) RCFTSlipMaterial.C, revA"


#include "RCFTSlipMaterial.h"
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>
#include <fstream>

//using std::ofstream;
//using std::endl;
//using std::ios;
using namespace std;


RCFTSlipMaterial::RCFTSlipMaterial(int tag, double fy, double e, double esh)
:UniaxialMaterial(tag,MAT_TAG_RCFTSlipMaterial),
 trialX(0.0), commitX(0.0), E(e), Esh(esh), fyp(fy), fyn(-fy), trialStrain(0.0), commitStrain(0.0),
 trialStress(0.0), commitStress(0.0), trialTangent(E), commitTangent(E)
{
}

RCFTSlipMaterial::RCFTSlipMaterial()
:UniaxialMaterial(0,MAT_TAG_RCFTSlipMaterial),
 trialX(0.0), commitX(0.0), E(0.0), Esh(0.0), fyp(0.0), fyn(0.0), trialStrain(0.0), commitStrain(0.0),
 trialStress(0.0), commitStress(0.0), trialTangent(0.0), commitTangent(0.0)
{
}

RCFTSlipMaterial::~RCFTSlipMaterial()
{
  // does nothing
}

int 
RCFTSlipMaterial::setTrialStrain(double strain, double strainRate)
{
#ifdef COMPOSITE_DEBUG
    ofstream output;
    output.open("setslip.dat",ios::app);
#endif

    trialStrain = commitStrain + strain;

    double ey_t = ( commitX + fyp - commitStress ) / E;

    double ey_c = ( commitX + fyn - commitStress ) / E;

    if ( strain > ey_t )
    {
       trialX = commitX + Esh * ( strain - ey_t );
       trialStress = trialX + fyp;
       trialTangent = Esh;
    }
    
    else if ( strain < ey_c )
    {
       trialX = commitX + Esh * ( strain - ey_c );
       trialStress = trialX + fyn;
       trialTangent = Esh;
    }
    else
    {
       trialX = commitX;
       trialStress = commitStress + E*strain;
       trialTangent = E;
    }
    
    return 0;
}

double 
RCFTSlipMaterial::getStrain(void)
{
  return trialStrain;
}

double 
RCFTSlipMaterial::getStress(void)
{
  return trialStress;
}


double 
RCFTSlipMaterial::getTangent(void)
{
  return trialTangent;
}

double
RCFTSlipMaterial::getCommitTangent(void)
{
  return commitTangent;
}

int 
RCFTSlipMaterial::commitState(void)
{
#ifdef COMPOSITE_DEBUG
	ofstream RCFTslip;
    RCFTslip.open("RCFTslip.dat",ios::app);
    RCFTslip<<trialStrain<<"  "<<trialStress<<"  "<<trialTangent<<"  "<<commitStrain<<"  "<<commitStress<<"  "<<commitTangent<<"  "<<E<<endl;
#endif
	commitTangent = trialTangent; 
    commitStrain  = trialStrain;
    commitStress  = trialStress;
    commitX = trialX;
    return 0;
}	


int 
RCFTSlipMaterial::revertToLastCommit(void)
{
    trialTangent = commitTangent;
    trialStrain  = commitStrain;
    trialStress  = commitStress;
    trialX = commitX;
    return 0;
}


int 
RCFTSlipMaterial::revertToStart(void)
{
    return 0;
}


UniaxialMaterial *
RCFTSlipMaterial::getCopy(void)
{
  RCFTSlipMaterial *theCopy =
    new RCFTSlipMaterial(this->getTag(),Fy,E,Esh);
  theCopy->fyp = this->fyp;
  theCopy->fyn = this->fyn;
  theCopy->E = this->E;
  theCopy->Esh = this->Esh;
  theCopy->Fy = this->Fy;
  theCopy->trialStrain = this->trialStrain;
  theCopy->commitStrain = this->commitStrain;
  theCopy->commitTangent = this->commitTangent;
  theCopy->commitStress = this->commitStress;
  theCopy->commitX = this->commitX;
  theCopy->trialX = this->trialX; 
  theCopy->trialStress = this->trialStress;
  theCopy->trialTangent = this->trialTangent;    
  
  return theCopy;
}


int 
RCFTSlipMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  return res;
}

int 
RCFTSlipMaterial::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  return res;
}

void 
RCFTSlipMaterial::Print(OPS_Stream &s, int flag)
{
    s << "RCFTSlipMaterial tag: " << this->getTag() << endln;
    s << "  E: " << E << endln;
    s << "  Otress: " << trialStress << " tangent: " << trialTangent << endln;
}


