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
// $Date: 2011-07-18 10:11:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam3dThermalAction.h,v $


// Written: Jian Jiang of Edinburgh University
// Created: Jul. 2011


// Description: This file contains the class definition for Beam3dThermalAction.
// Beam3dThermalAction is a thermal field class created to store the temperature
// distribution through the depth of section defined by temperature and location.


#ifndef Beam3dThermalAction_h
#define Beam3dThermalAction_h


#include <ElementalLoad.h>

class Beam3dThermalAction : public ElementalLoad
{
  public:
  // Constructors based on 9, 5 or 2 temperature points 
  // t-temperature; locY-coordinate through the depth of section
  Beam3dThermalAction(int tag,
                double t1, double locY1, double t2, double locY2,
                double t3, double locY3, double t4, double locY4,
                double t5, double locY5, double t6, double locZ1,
                double t7, double locZ2, double t8, double locZ3,
                double t9, double locZ4, double t10, double locZ5,
		        int theElementTag);

  Beam3dThermalAction();    

  ~Beam3dThermalAction();
  //--Adding declaration of 'setfactors' for 'UoE-FireLoadPattern' [-BEGIN-]: by L.J&P.K(university of Edinburgh)-07-MAY-2012-///
  virtual void setfactors(int value, double loadFactor2, double loadFactor3, double loadFactor4, 
	  double loadFactor5, double loadFactor6, double loadFactor7, 
	  double loadFactor8, double loadFactor9);  //Not Uodated
  Vector &getfactors();
  //--Adding declaration of 'setfactors' for 'UoE-FireLoadPattern' [-END-]: by L.J&P.K(university of Edinburgh)-07-MAY-2012-///
  
  const Vector &getData(int &type, double loadFactor);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel,  
	       FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);       
  
 protected:
  
 private:
  double T1; //Temperature
  double LocY1; // Location through the depth of section
  double T2;
  double LocY2;
  double T3;
  double LocY3;
  double T4;
  double LocY4;
  double T5;
  double LocY5;
  double T6;
  double LocZ1;
  double T7;
  double LocZ2;
  double T8;
  double LocZ3;
  double T9;
  double LocZ4;
  double T10;
  double LocZ5;
  static Vector data; // data for temperature and locations

  //--Adding a factor vector for FireLoadPattern [-BEGIN-]: by L.J&P.K(university of Edinburgh)-07-MAY-2012-///
  static Vector factors;
  int indicator; //indicator if fireloadpattern was called
  double Factor2;
  double Factor3;
  double Factor4;
  double Factor5;
  double Factor6;
  double Factor7;
  double Factor8;
  double Factor9;
  double Factor10;
  //--Adding a factor vector for FireLoadPattern [-END-]: by L.J&P.K(university of Edinburgh)-07-MAY-2012-///
};

#endif

