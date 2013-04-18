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

                                                                        
#ifndef ShellThermalAction_h
#define ShellThermalAction_h

// JZ 07/10


#include <ElementalLoad.h>

class ShellThermalAction : public ElementalLoad
{
  public:
  // Constructors based on 9, 2 or 0 temperature changes given
  ShellThermalAction(int tag,
                double t1, double locY1, double t2, double locY2,
                double t3, double locY3, double t4, double locY4,
                double t5, double locY5, double t6, double locY6,
                double t7, double locY7, double t8, double locY8,
                double t9, double locY9, 
		        int theElementTag);

  ShellThermalAction(int tag,
                double t1, double locY1, double t2, double locY2,
				double t3, double locY3, double t4, double locY4,
                double t5, double locY5, int theElementTag);
    ShellThermalAction(int tag,
                double t1, double locY1, double t2, double locY2,
                int theElementTag);

  ShellThermalAction();    

  ~ShellThermalAction();
  
  const Vector &getData(int &type, double loadFactor);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel,  
	       FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);       
  
 protected:
  
 private:
  double T1;
  double LocY1;
  double T2;
  double LocY2;
  double T3;
  double LocY3;
  double T4;
  double LocY4;
  double T5;
  double LocY5;
  double T6;
  double LocY6;
  double T7;
  double LocY7;
  double T8;
  double LocY8;
  double T9;
  double LocY9;
  static Vector data; // data for temp loads and locs
};

#endif

