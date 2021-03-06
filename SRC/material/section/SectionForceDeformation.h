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
                                                                        
// $Revision: 1.14 $
// $Date: 2008/08/26 16:45:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/SectionForceDeformation.h,v $
                                                                        
                                                                        
#ifndef SectionForceDeformation_h
#define SectionForceDeformation_h

// Written: MHS
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class definition for SectionForceDeformation.
// SectionForceDeformation is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. 
//
// What: "@(#) SectionForceDeformation.h, revA"

#include <Material.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Information;
class Response;

#define MAX_SECTION_RESPONSE_ID 10000

#define SECTION_RESPONSE_MZ		1
#define SECTION_RESPONSE_P		2
#define SECTION_RESPONSE_VY		3
#define SECTION_RESPONSE_MY		4
#define SECTION_RESPONSE_VZ		5
#define SECTION_RESPONSE_T		6
#define SECTION_RESPONSE_R              7
#define SECTION_RESPONSE_Q              8

#define SECTION_RESPONSE_S              17

#define SECTION_RESPONSE_Pc             18
#define SECTION_RESPONSE_MYc            19
#define SECTION_RESPONSE_MZc            20
#define SECTION_RESPONSE_LYYc           21
#define SECTION_RESPONSE_LZZc           22
#define SECTION_RESPONSE_LYZc           23

#define SECTION_RESPONSE_Ps             24
#define SECTION_RESPONSE_MYs            25
#define SECTION_RESPONSE_MZs            26
#define SECTION_RESPONSE_LYYs           27
#define SECTION_RESPONSE_LZZs           28
#define SECTION_RESPONSE_LYZs           29

class SectionForceDeformation : public Material
{
 public:
  SectionForceDeformation (int tag, int classTag);
  SectionForceDeformation ();
  virtual ~SectionForceDeformation ();
  
  //virtual int setTrialSectionDeformation (const Vector&) = 0;
  virtual int setTrialSectionDeformation (const Vector&) ; //the default valuoe 0 is removeed byJZ ,UoE 
  virtual const Vector &getSectionDeformation (void) = 0;
  
  virtual const Vector &getStressResultant (void) = 0;
  virtual const Matrix &getSectionTangent (void) = 0;
  virtual const Matrix &getInitialTangent (void) = 0;
  virtual const Matrix &getSectionFlexibility (void);
  virtual const Matrix &getInitialFlexibility (void);
  
  virtual double getRho(void);
  virtual double getZh(void) {return 0;};
  virtual double getYh(void) {return 0;};
  virtual double getEIy(void) {return 0;};
  virtual double getEIz(void) {return 0;};
  virtual double getGAy(void) {return 0;};
  virtual double getGAz(void) {return 0;};

  virtual int commitState (void) = 0;
  virtual int revertToLastCommit (void) = 0;
  virtual int revertToStart (void) = 0;
  
  virtual SectionForceDeformation *getCopy (void) = 0;
  virtual const ID &getType (void) = 0;
  virtual int getOrder (void) const = 0;
  
  virtual Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  virtual int getResponse(int responseID, Information &info);

  virtual int getResponseSensitivity(int responseID, int gradIndex,
				     Information &info);
  
  virtual int revertToLastCommit(double) {return 0;};
  virtual int setTrialSectionDeformation(const Vector&, double L) {return 0;};

  // AddingSensitivity:BEGIN //////////////////////////////////////////
  virtual const Vector &getStressResultantSensitivity(int gradIndex,
						      bool conditional);
  virtual const Vector &getSectionDeformationSensitivity(int gradIndex);
  virtual const Matrix &getSectionTangentSensitivity(int gradIndex);
  virtual const Matrix &getSectionFlexibilitySensitivity(int gradIndex);
  virtual const Matrix &getInitialTangentSensitivity(int gradIndex);
  virtual const Matrix &getInitialFlexibilitySensitivity(int gradIndex);
  virtual double getRhoSensitivity(int gradIndex);
  virtual int commitSensitivity(const Vector& sectionDeformationGradient,
				int gradIndex, int numGrads);
  // AddingSensitivity:END ///////////////////////////////////////////
  
  //--- Adding Thermal Materials:[BEGIN]   by UoE OpenSees Group ----//  
  virtual int setTrialSectionDeformation(const Vector&, const Vector &); //JZ
  virtual const Vector &getTemperatureStress(const Vector &tData);//27 is for 'FireLoadPattern'
  virtual int setTrialSectionDeformationTemperature(const Vector&, double *); //JZ
  //--- Adding Thermal Functions:[END]   by UoE OpenSees Group ----//

 protected:
  Matrix *fDefault;	// Default flexibility matrix
  Vector *sDefault;
  
 private:

};

extern bool OPS_addSectionForceDeformation(SectionForceDeformation *newComponent);
extern SectionForceDeformation *OPS_getSectionForceDeformation(int tag);
extern void OPS_clearAllSectionForceDeformation(void);

#endif
