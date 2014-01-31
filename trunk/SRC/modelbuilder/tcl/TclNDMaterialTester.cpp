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
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/TclnDMaterialTester.cpp,v $
                                                                        
// File: ~/modelbuilder/tcl/TclnDMaterialTester.C
// 
// Written: fmk 
// Created: 03/01
//
// Description: This file contains the implementaion of the TclnDMaterialTester class.
//
// What: "@(#) TclnDMaterialTester.C, revA"

#include <stdlib.h>
#include <string.h>

#include <ArrayOfTaggedObjects.h>
#include <nDMaterial.h>
#include <TclNDMaterialTester.h>

#include <Matrix.h>
#include <elementAPI.h>
//
// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER
//

static TclNDMaterialTester *theTclBuilder =0;
static NDMaterial *theTestingNDMaterial = 0;

// 
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//

int  TclNDMaterialTester_setNDMaterial(ClientData clientData, Tcl_Interp *interp, 
						   int argc,   TCL_Char **argv);
				    
int  TclNDMaterialTester_setStrainNDMaterial(ClientData clientData, Tcl_Interp *interp,
							 int argc,   TCL_Char **argv);

int  TclNDMaterialTester_getStressNDMaterial(ClientData clientData, Tcl_Interp *interp,
							 int argc,   TCL_Char **argv);

int  TclNDMaterialTester_getTangNDMaterial(ClientData clientData, Tcl_Interp *interp,
						       int argc,   TCL_Char **argv);

//
// CLASS CONSTRUCTOR & DESTRUCTOR
//

static int count;
static int countsTillCommit;
				    
// constructor: the constructor will add certain commands to the interpreter
TclNDMaterialTester::TclNDMaterialTester(Domain &theDomain, Tcl_Interp *interp, int cTC)
  :TclModelBuilder(theDomain, interp, 1, 1), theInterp(interp)
{
  countsTillCommit = cTC;
  Tcl_CreateCommand(interp, "nDTest", TclNDMaterialTester_setNDMaterial,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "strainNDTest", TclNDMaterialTester_setStrainNDMaterial,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "stressNDTest", TclNDMaterialTester_getStressNDMaterial,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "tangNDTest", TclNDMaterialTester_getTangNDMaterial,
		    (ClientData)NULL, NULL);
  
  
  // set the static pointers in this file
  theTclBuilder = this;
}

TclNDMaterialTester::~TclNDMaterialTester()
{

  theTclBuilder =0;

  Tcl_DeleteCommand(theInterp, "nDTest");
  Tcl_DeleteCommand(theInterp, "strainNDTest");
  Tcl_DeleteCommand(theInterp, "stressNDTest");
  Tcl_DeleteCommand(theInterp, "tangNDTest");
}


//
// THE FUNCTIONS INVOKED BY THE INTERPRETER
//

int
TclNDMaterialTester_setNDMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   
					      TCL_Char **argv)
{
  count = 1;
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    Tcl_SetResult(interp, "WARNING builder has been destroyed", TCL_STATIC);
    return TCL_ERROR;
  }

  // check number of arguments in command line
  if (argc < 2) {
    Tcl_SetResult(interp, "WARNING bad command - want: nDTest matID?", TCL_STATIC);
    return TCL_ERROR;
  }    

  // get the matID form command line
  int matID;
  if (Tcl_GetInt(interp, argv[1], &matID) != TCL_OK) {
    Tcl_SetResult(interp, "WARNING could not read matID: nDTest matID?", TCL_STATIC);
    return TCL_ERROR;
  }

  // delete the old testing material
  if (theTestingNDMaterial !=0) {
    delete [] theTestingNDMaterial;
    theTestingNDMaterial = 0;
  }

  // get the material from the modelbuilder with matID 
  // and set the testing material to point to a copy of it
  NDMaterial *theOrigMaterial = OPS_GetNDMaterial(matID);
  if (theOrigMaterial == 0) {
    Tcl_SetResult(interp, "WARNING no material found with matID", TCL_STATIC);
    return TCL_ERROR;
  }  else {
    theTestingNDMaterial = theOrigMaterial->getCopy();
  }

  return TCL_OK;
}

int  
TclNDMaterialTester_setStrainNDMaterial(ClientData clientData, Tcl_Interp *interp,
						    int argc,   TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    Tcl_SetResult(interp, "WARNING builder has been destroyed", TCL_STATIC);
    return TCL_ERROR;
  }

  // check number of arguments in command line
  if (argc < 7) {
    Tcl_SetResult(interp, "WARNING bad command - want: strainNDTest eps_11 eps_22 eps_33 eps_12 eps_23 eps_13?", TCL_STATIC);
    return TCL_ERROR;
  }    

  // get the mat strain form command line
  static Vector data(6);
  double strain;
  for (int i = 1; i < argc; i++) {
	if (Tcl_GetDouble(interp, argv[i], &strain) != TCL_OK) {
	  Tcl_SetResult(interp, "WARNING could not read strain: strainNdTest strain?", TCL_STATIC);
    return TCL_ERROR;
	}
	data(i - 1) = strain;
  }
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
TclNDMaterialTester_getStressNDMaterial(ClientData clientData, Tcl_Interp *interp,
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
    Tcl_SetResult(interp, "WARNING no active NDMaterial - use nDTest command", TCL_STATIC);    
    return TCL_ERROR;
  }
}

int
TclNDMaterialTester_getTangNDMaterial(ClientData clientData, Tcl_Interp *interp,
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
    Tcl_SetResult(interp, "WARNING no active NDMaterial - use nDTest command", TCL_STATIC);    
    return TCL_ERROR;
  }
}
