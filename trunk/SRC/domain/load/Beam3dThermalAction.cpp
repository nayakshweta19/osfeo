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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam3dThermalAction.cpp,v $


// Written: Jian Jiang of Edinburgh University
// Created: Jul. 2011


// Description: This file contains the class implementation for Beam3dThermalAction.
// Beam3dThermalAction is a thermal field class created to store the temperature
// distribution through the depth of section defined by temperature and location.


#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <Beam3dThermalAction.h>
#include <Vector.h>

Vector Beam3dThermalAction::data(20);
Vector Beam3dThermalAction::factors(9);

Beam3dThermalAction::Beam3dThermalAction(int tag, 
                         double t1, double locY1, double t2, double locY2,
                         double t3, double locY3, double t4, double locY4,
                         double t5, double locY5, double t6, double locZ1,
                         double t7, double locZ2, double t8, double locZ3,
                         double t9, double locZ4, double t10,double locZ5,
			             int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam3dThermalAction, theElementTag), 
T1(t1),LocY1(locY1),T2(t2),LocY2(locY2),T3(t3),LocY3(locY3),T4(t4),LocY4(locY4),
T5(t5),LocY5(locY5),T6(t6),LocZ1(locZ1),T7(t7),LocZ2(locZ2),T8(t8),LocZ3(locZ3),
T9(t9),LocZ4(locZ4),T10(t10),LocZ5(locZ5)
{

}




Beam3dThermalAction::~Beam3dThermalAction()
{

}

const Vector &
Beam3dThermalAction::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam3dThermalAction;
  data(0) = T1;
  data(1) = LocY1;
  data(2) = T2;
  data(3) = LocY2;
  data(4) = T3;
  data(5) = LocY3;
  data(6) = T4;
  data(7) = LocY4;
  data(8) = T5;
  data(9) = LocY5;
  data(10) = T6;
  data(11) = LocZ1;
  data(12) = T7;
  data(13) = LocZ2;
  data(14) = T8;
  data(15) = LocZ3;
  data(16) = T9;
  data(17) = LocZ4;
  data(18) = T10;
  data(19) = LocZ5;
  return data;
}
///Adding Loadfactors to Load 'Beam2dThermalAction' for FireLoadPattern [-BEGIN-]:by L.J&P.K(university of Edinburgh)-07-MAY-2012-///
void
Beam3dThermalAction::setfactors(int value, double loadFactor2, double loadFactor3, double loadFactor4, 
	  double loadFactor5, double loadFactor6, double loadFactor7, 
	  double loadFactor8, double loadFactor9)
{
indicator=value; //1 indicates fireloadpattern was called
Factor2=loadFactor2;
Factor3=loadFactor3;
Factor4=loadFactor4;
Factor5=loadFactor5;
Factor6=loadFactor6;
Factor7=loadFactor7;
Factor8=loadFactor8;
Factor9=loadFactor9;
Factor10=loadFactor9;  //Not Updated, Liming,2013
//opserr << "beam2d loadfactor2\n";
//opserr << loadFactor2;

factors(0) = indicator;
factors(1) = Factor2;
factors(2) = Factor3;
factors(3) = Factor4;
factors(4) = Factor5;
factors(5) = Factor6;
factors(6) = Factor7;
factors(7) = Factor8;
factors(8) = Factor9;
factors(9) = Factor10;
//opserr << "Beam2dTemperatureLoad set factors: ";
  
}

Vector &
Beam3dThermalAction::getfactors()
{
 // opserr << "Beam2dThermalAction get factors: ";
  return factors;
  //opserr << "factors2: "
	//<< factors << "beam2dThermal\n";
}
///Adding Loadfactors to Load 'Beam2dThermalAction' for FireLoadPattern [-END-]:by L.J&P.K(university of Edinburgh)-07-MAY-2012-///

int 
Beam3dThermalAction::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
Beam3dThermalAction::recvSelf(int commitTag, Channel &theChannel,  
			 FEM_ObjectBroker &theBroker)
{
  return -1;
}

// do it later
void 
Beam3dThermalAction::Print(OPS_Stream &s, int flag)
{
  s << "Beam3dThermalAction - reference load : " << T1 <<" change  temp of bot\n";
  s <<  T2 << " change  temp at top\n";
  s << "  element acted on: " << eleTag << endln;
}

