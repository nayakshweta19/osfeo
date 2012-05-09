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
                                                                        
// $Revision: 1.1 $
// $Date: 2008/12/03 23:46:59 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/TimoshenkoSection2d.h,v $
                                                                        
// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition for 
// TimoshenkoSection2d.h. TimoshenkoSection2d provides the abstraction of a 
// rectangular section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.
// The fiber stresses are the 11, 12, and 13 components of stress, from
// which all six beam stress resultants are obtained.

#ifndef TimoshenkoSection2d_h
#define TimoshenkoSection2d_h

#include <SectionForceDeformation.h>
#include <Vector.h>
#include <Matrix.h>

class NDMaterial;
class Fiber;
class Response;

class TimoshenkoSection2d : public SectionForceDeformation
{
 public:
  TimoshenkoSection2d(); 
  TimoshenkoSection2d(int tag, int num, Fiber **fibers);
  ~TimoshenkoSection2d();
  
  const char *getClassType(void) const {return "TimoshenkoSection2d";};

  int   setTrialSectionDeformation(const Vector &deforms); 
  const Vector &getSectionDeformation(void);
  
  const Vector &getStressResultant(void);
  const Matrix &getSectionTangent(void);
  const Matrix &getInitialTangent(void) {return this->getSectionTangent();}
  
  int   commitState(void);
  int   revertToLastCommit(void);    
  int   revertToStart(void);
  
  SectionForceDeformation *getCopy(void);
  const ID &getType(void);
  int getOrder(void) const;
  double getZh(void);
  double getEIz(void);
  double getGAy(void);

  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag = 0);
  
  Response *setResponse(const char **argv, int argc, 
			  OPS_Stream &s);
  int getResponse(int responseID, Information &info);

  int addFiber(Fiber &theFiber);
  
  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int setParameter(const char **argv, int argc, Parameter &param);

 protected:
  
 private:
  
  int numFibers;                   // number of fibers in the section
  NDMaterial **theMaterials;       // array of pointers to materials
  double   *matData;               // data for the materials [yloc and area]
  double   kData[9];               // data for ks matrix 
  double   sData[3];               // data for s vector 

  double yBar;       // Section centroid
  double zBar;

  double yh;         // Section Nominal Hight
  double zh;
  
  static ID code;

  Vector e;          // trial section deformations 
  Vector eCommit;    // committed section deformations 
  Vector *s;         // section resisting forces  (axial force, bending moment)
  Matrix *ks;        // section stiffness
};

#endif
