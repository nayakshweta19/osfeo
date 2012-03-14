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
// $Date: 2007/02/02 01:18:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/fiber/TriaxialFiber.h,v $
                                                                        
                                                                        
// File: ~/fiber/TriaxialFiber.h
//
// Written: Remo Magalhaes de Souza
// Created: 10/98
// Revision: 
//
// Description: This file contains the class definition for 
// TriaxialFiber.h. TriaxialFiber provides the abstraction of a
// uniaxial fiber that forms a fiber section for 3d frame elements (the 
// fiber position inside the section is defined by two coordinates)
// The TriaxialFiber is subjected to a stress state with 
// only one nonzero axial stress and corresponding axial strain.
//
// What: "@(#) TriaxialFiber.h, revA"

#ifndef TriaxialFiber_h
#define TriaxialFiber_h

#include <Matrix.h>
#include <Fiber.h>

class UniaxialMaterial;
class NDMaterial;
class Response;

class TriaxialFiber: public Fiber
{
  public:
    TriaxialFiber ();    
    TriaxialFiber (int tag, NDMaterial &theMat, double Area, 
                     const Vector &position);
 
    ~TriaxialFiber();

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
    
    Response *setResponse(const char **argv, int argc, OPS_Stream &S);
    int getResponse(int responseID, Information &info);

    void getFiberLocation(double &y, double &z);
    UniaxialMaterial *getMaterial(void) {return 0;};
	NDMaterial *getNDMaterial(void) {return theMaterial;};
    double getArea(void) {return area;};

  protected:
    
  private:
    NDMaterial *theMaterial;   // pointer to a material
    double area;                          // area of the fiber 
    double as[2];                            // matrix that transforms
	                            // section deformations into fiber strain	
    static Matrix ks;       // static class wide matrix object for returns
    static Vector fs;	    // static class wide vector object for returns					
    static ID code;
};


#endif

