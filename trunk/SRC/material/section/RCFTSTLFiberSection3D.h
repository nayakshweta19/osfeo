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
// $Date: 2007-09-21 15:28:40 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/RCFTSTLFiberSection3D.h,v $
                                                                        
// Written: fmk
// Created: 04/01
//
// Description: This file contains the class definition for 
// RCFTSTLFiberSection3D.h. RCFTSTLFiberSection3D provides the abstraction of a 
// 3d beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.

#ifndef RCFTSTLFiberSection3D_h
#define RCFTSTLFiberSection3D_h

#include <SectionForceDeformation.h>
#include <Vector.h>
#include <Matrix.h>

class UniaxialMaterial;
class Fiber;
class Response;

class RCFTSTLFiberSection3D : public SectionForceDeformation
{
 public:
  RCFTSTLFiberSection3D(); 
  RCFTSTLFiberSection3D(int tag, int numFibers, Fiber **fibers, double GJ = 1.0e10); 
  ~RCFTSTLFiberSection3D();
  
  int   setTrialSectionDeformation(const Vector &deforms); 
  const Vector &getSectionDeformation(void);
  
  const Vector &getStressResultant(void);
  const Matrix &getSectionTangent(void);
  const Matrix &getInitialTangent(void);
  
  int   commitState(void);
  int   revertToLastCommit(void);    
  int   revertToStart(void);
  
  SectionForceDeformation *getCopy(void);
  const ID &getType (void);
  int getOrder (void) const;
  
  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag = 0);
  
  Response *setResponse(const char **argv, int argc, OPS_Stream &output);
  int getResponse(int responseID, Information &info);
  
  int addFiber(Fiber &theFiber);
  
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  
 protected:
  
 private:
  int numFibers;                   // number of fibers in the section
  UniaxialMaterial **theMaterials; // array of pointers to materials
  double *matData;               // data for the materials [yloc and area]
  double kData[16];               // data for ks matrix 
  double sData[3];               // data for s vector 
  
  double yBar;       // Section centroid
  double zBar;
  
  Vector e;          // trial section deformations 
  Vector eCommit;    // committed section deformations 
  
  static ID code;
  static Vector s;         // section resisting forces
  static Matrix ks;        // section stiffness
  
  double GJ;
};

#endif
