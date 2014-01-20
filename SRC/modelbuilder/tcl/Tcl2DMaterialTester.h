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
                                                                        
// $Revision: 1.4 $
// $Date: 2007/06/26 20:13:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/TclnDMaterialTester.h,v $
                                                                        
// File: ~/modelbuilder/tcl/TclnDMaterialTester.h
// 
// Written: fmk 
// Created: 03/01
// Revision: A
//
// Description: This file contains the class definition for TclnDMaterialTester.
// A TclnDMaterialTester adds the commands to create and test nD materials
//
// What: "@(#) TclnDMaterialTester.h, revA"

#ifndef Tcl2DMaterialTester_h
#define Tcl2DMaterialTester_h

#include <TclModelBuilder.h>

class SectionForceDeformation;
class SectionRepres;
class UniaxialMaterial;
class NDMaterial;
class TaggedObjectStorage;

class CrdTransf2d;
class CrdTransf3d;

#include <tcl.h>

class Tcl2DMaterialTester : public TclModelBuilder
{
  public:
    Tcl2DMaterialTester(Domain &theDomain,Tcl_Interp *interp, int count=1);
    ~Tcl2DMaterialTester();    

  protected:

  private:
    Tcl_Interp *theInterp;
};

#endif







