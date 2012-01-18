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
                                                                        
// $Revision: 1.9 $
// $Date: 2011/06/16 01:18:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/fiber/BiaxialFiber2d.h,v $
                                                                        
                                                                        
// File: ~/fiber/BiaxialFiber2d.h
//
// Written: Neallee
// Created: 2011
// Revision: 
//
// Description: This file contains the class definition for 
// BiaxialFiber2d.h. BiaxialFiber2d provides the abstraction of a
// biaxial fiber whose position is defined with only one coordinate.
// The BiaxialFiber2d is subjected to a stress state with 
// nonzero axial and shear stresses and corresponding strains.
//
// What: "@(#) BiaxialFiber2d.h, revA"

#ifndef biaxialFiber2d_h
#define biaxialFiber2d_h

#include <Fiber.h>
#include <Vector.h>
#include <Matrix.h>

class UniaxialMaterial;
class NDMaterial;
class Parameter;

class BiaxialFiber2d : public Fiber
{
  public:
    BiaxialFiber2d ();   
    BiaxialFiber2d (int tag, NDMaterial &theMat, double Area, double position);
    ~BiaxialFiber2d();

    
    int   setTrialFiberStrain(const Vector &vs);
    Vector &getFiberStressResultants (void);
    Matrix &getFiberTangentStiffContr (void);

    int   commitState(void);
    int   revertToLastCommit(void);    
    int   revertToStart(void);
    
    Fiber *getCopy(void);
    int getOrder(void);
    const ID &getType(void);

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0);
    
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &info);
	
    void getFiberLocation(double &y, double &z);
	UniaxialMaterial *getMaterial(void) {return 0;};
    NDMaterial *getNDMaterial(void) {return theMaterial;};
    double getArea(void) {return area;};

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    int activateParameter(int parameterID);

    const Vector &getFiberSensitivity(int gradNumber, bool cond);
    int commitSensitivity(const Vector &dedh, int gradNumber,
			  int numGrads);

  protected:
    
  private:
    NDMaterial *theMaterial;   // pointer to a material
    double area;                          // area of the fiber 
    double y;		// fiber location

    static Matrix ks;       // static class wide matrix object for returns
    static Vector fs;	    // static class wide vector object for returns

    static ID code;
};


#endif






