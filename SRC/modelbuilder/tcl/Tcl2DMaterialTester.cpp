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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-02-25 23:34:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/Tcl2DMaterialTester.cpp,v $
                                                                        
// File: ~/modelbuilder/tcl/Tcl2DMaterialTester.C
// 
// Written: fmk 
// Created: 03/01
//
// Description: This file contains the implementation of the Tcl2DMaterialTester class.
//
// What: "@(#) Tcl2DMaterialTester.C, revA"

#include <stdlib.h>
#include <string.h>

#include <ArrayOfTaggedObjects.h>
#include <nDMaterial.h>
#include <Tcl2DMaterialTester.h>

#include <Matrix.h>
#include <elementAPI.h>
//
// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER
//

static Tcl2DMaterialTester *theTclBuilder = 0;
static NDMaterial *theTestingNDMaterial = 0;

// 
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//

int  Tcl2DMaterialTester_setNDMaterial(ClientData clientData, Tcl_Interp *interp, 
							int argc,   TCL_Char **argv);
				    
int  Tcl2DMaterialTester_setStrainNDMaterial(ClientData clientData, Tcl_Interp *interp,
							int argc,   TCL_Char **argv);

int  Tcl2DMaterialTester_getStressNDMaterial(ClientData clientData, Tcl_Interp *interp,
							int argc,   TCL_Char **argv);

int  Tcl2DMaterialTester_getTangNDMaterial(ClientData clientData, Tcl_Interp *interp,
							int argc,   TCL_Char **argv);

//
// CLASS CONSTRUCTOR & DESTRUCTOR
//

static int count;
static int countsTillCommit;
				    
// constructor: the constructor will add certain commands to the interpreter
Tcl2DMaterialTester::Tcl2DMaterialTester(Domain &theDomain, Tcl_Interp *interp, int cTC)
  :TclModelBuilder(theDomain, interp, 1, 1), theInterp(interp)
{
  countsTillCommit = cTC;
  Tcl_CreateCommand(interp, "planeStressTest", Tcl2DMaterialTester_setNDMaterial,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "strain2DTest", Tcl2DMaterialTester_setStrainNDMaterial,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "stress2DTest", Tcl2DMaterialTester_getStressNDMaterial,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "tang2DTest", Tcl2DMaterialTester_getTangNDMaterial,
		    (ClientData)NULL, NULL);
  
  
  // set the static pointers in this file
  theTclBuilder = this;
}

Tcl2DMaterialTester::~Tcl2DMaterialTester()
{

  theTclBuilder =0;

  Tcl_DeleteCommand(theInterp, "planeStressTest");
  Tcl_DeleteCommand(theInterp, "strain2DTest");
  Tcl_DeleteCommand(theInterp, "stress2DTest");
  Tcl_DeleteCommand(theInterp, "tang2DTest");
}

//
// THE FUNCTIONS INVOKED BY THE INTERPRETER
//

int
Tcl2DMaterialTester_setNDMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   
					      TCL_Char **argv)
{
  count = 1;
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed \n";
    return TCL_ERROR;
  }

  // check number of arguments in command line
  if (argc < 2) {
    opserr << "WARNING bad command - want: planeStressTest matID?\n";
    return TCL_ERROR;
  }

  // get the matID form command line
  int matID;
  if (Tcl_GetInt(interp, argv[1], &matID) != TCL_OK) {
    opserr << "WARNING could not read matID: planeStressTest matID?\n";
    return TCL_ERROR;
  }

  // delete the old testing material
  if (theTestingNDMaterial !=0) {
    delete theTestingNDMaterial;
    theTestingNDMaterial = 0;
  }

  // get the material from the modelbuilder with matID 
  // and set the testing material to point to a copy of it
  NDMaterial *theOrigMaterial = theTclBuilder->getNDMaterial(matID);
  if (theOrigMaterial == 0) {
    opserr << "WARNING no material found with matID\n";
    return TCL_ERROR;
  }  else {
    theTestingNDMaterial = theOrigMaterial;
  }

  return TCL_OK;
}

int  
Tcl2DMaterialTester_setStrainNDMaterial(ClientData clientData, Tcl_Interp *interp,
						    int argc,   TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";
    return TCL_ERROR;
  }

  // check number of arguments in command line
  if (argc < 4) {
	opserr << "WARNING bad command - want: strain2DTest eps_x eps_y gamma?\n";
    return TCL_ERROR;
  }    

  // get the mat strain form command line
  static Vector data(3);
  double strain[3];
  int numData = 3;
  if (Tcl_GetDouble(interp, argv[1], &strain[0]) != TCL_OK) {
	opserr << "WARNING could not read strain: strainNDTest eps_x?\n";
	return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &strain[1]) != TCL_OK) {
	  opserr << "WARNING could not read strain: strainNDTest eps_y?\n";
	  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &strain[2]) != TCL_OK) {
	  opserr << "WARNING could not read strain: strainNDTest gamma?\n";
	  return TCL_ERROR;
  }
  data = Vector(strain, numData);

  // delete the old testing material
  if (theTestingNDMaterial !=0) {
    theTestingNDMaterial->setTrialStrain(data);
    if (count == countsTillCommit) {
      theTestingNDMaterial->commitState();    
      count = 1;
    } else count++;
  }
  return TCL_OK;
}

int
Tcl2DMaterialTester_getStressNDMaterial(ClientData clientData, Tcl_Interp *interp,
							 int argc,   TCL_Char **argv)
{
  char buffer[120];
  if (theTestingNDMaterial !=0) {
	const Vector &stress = theTestingNDMaterial->getStress();
    for (int i=0; i<stress.Size(); i++) {
      sprintf(buffer,"%.10e ",stress(i));
	  Tcl_AppendResult(interp, buffer, NULL);
	}
	return TCL_OK;
  } else {
    opserr << "WARNING no active NDMaterial - use planeStressTest command\n";    
    return TCL_ERROR;
  }
}

int
Tcl2DMaterialTester_getTangNDMaterial(ClientData clientData, Tcl_Interp *interp,
						       int argc,   TCL_Char **argv)
{
  char buffer[120];
  if (theTestingNDMaterial !=0) {
	const Matrix & tangent = theTestingNDMaterial->getTangent();
	for (int i = 0; i < tangent.noRows(); i++) {
	  for (int j = 0; j < tangent.noCols(); j++) {
	    sprintf(buffer, "%.10e ", tangent(i, j));
		Tcl_AppendResult(interp, buffer, NULL);
	  }
	  Tcl_AppendResult(interp, endln, NULL);
	}	
    return TCL_OK;
  } else {
    opserr << "WARNING no active NDMaterial - use planeStressTest command\n";    
    return TCL_ERROR;
  }
}
