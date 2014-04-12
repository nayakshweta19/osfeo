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
                                                                        
// $Revision: 1.3 $
// $Date: 2003/02/25 23:32:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/feap/TclFeapElementCommand.cpp,v $
                                                                        
                                                                        
// File: ~/element/TclElementCommands.C
// 
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the implementation of the TclElementCommands.
// The file contains the routine TclElementCommands which is invoked by the
// TclModelBuilder.
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <fElmt01.h>
#include <fElmt02.h>
#include <fElmt03.h>

#include <TclModelBuilder.h>

extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_addFeapTruss(ClientData clientData, Tcl_Interp *interp, int argc, 
			     TCL_Char **argv, Domain *theTclDomain, TclModelBuilder *theTclBuilder,
			     int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }

  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  if (ndm != 2 && ndf != 2) {
	  opserr << "WARNING - fTruss eleTag? iNode? jNode? A? E? needs ndm=2, ndf=2\n";
	  return TCL_ERROR;
  }

  // check the number of arguments is correct
  if ((argc-eleArgStart) < 5) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: element fTruss eleTag? iNode? jNode? A? E?\n";

    return TCL_ERROR;
  }    


  // get the id and end nodes 
  int trussId, iNode, jNode;
  double A,E;
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &trussId) != TCL_OK) {
    opserr << "WARNING invalid truss eleTag" << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    opserr << "truss element: " << trussId << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
     opserr << "WARNING invalid jNode\n";
     opserr << "truss element: " << trussId << endln;
     return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4+eleArgStart], &A) != TCL_OK) {
	  opserr << "WARNING invalid A\n";
	  opserr << "truss element: " << trussId << endln;
	  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[5+eleArgStart], &E) != TCL_OK) {
	  opserr << "WARNING invalid E\n";
	  opserr << "truss element: " << trussId << endln;
	  return TCL_ERROR;
  }

 
   // now create the truss and add it to the Domain
   fElmt02 *theTruss = new fElmt02(trussId,iNode,jNode,A,E);
   if (theTruss == 0) {
	  opserr << "WARNING ran out of memory creating element\n";
	  opserr << "truss element: " << trussId << endln;
	  return TCL_ERROR;
   }

   if (theTclDomain->addElement(theTruss) == false) {
	  opserr << "WARNING could not add element to the domain\n";
	  opserr << "truss element: " << trussId << endln;
	  delete theTruss;
	  return TCL_ERROR;
   }
   return TCL_OK;
}

int
TclModelBuilder_addFeapFrame2D(ClientData clientData, Tcl_Interp *interp, int argc, 
			     TCL_Char **argv, Domain *theTclDomain, TclModelBuilder *theTclBuilder,
			     int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }

  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  if (ndm != 2 && ndf != 3) {
	  opserr << "WARNING - fFrame2D eleTag? iNode? jNode? A? E? needs ndm=2, ndf=3\n";
	  return TCL_ERROR;
  }

  // check the number of arguments is correct
  if ((argc-eleArgStart) < 5) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: element fFrame2D eleTag? iNode? jNode? A? E?\n";

    return TCL_ERROR;
  }    

  // get the id and end nodes 
  int eleId, iNode, jNode;
  double A,E;
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &eleId) != TCL_OK) {
    opserr << "WARNING invalid fFrame2D eleTag" << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    opserr << "fFrame2D element: " << eleId << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
     opserr << "WARNING invalid jNode\n";
     opserr << "fFrame2D element: " << eleId << endln;
     return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4+eleArgStart], &A) != TCL_OK) {
	  opserr << "WARNING invalid A\n";
	  opserr << "fFrame2D element: " << eleId << endln;
	  return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[5+eleArgStart], &E) != TCL_OK) {
	  opserr << "WARNING invalid E\n";
	  opserr << "fFrame2D element: " << eleId << endln;
	  return TCL_ERROR;
  }

   // now create the truss and add it to the Domain
   fElmt01 *theEle = new fElmt01(eleId,iNode,jNode,A,E);
   if (theEle == 0) {
	  opserr << "WARNING ran out of memory creating element\n";
	  opserr << "fFrame2D element: " << eleId << endln;
	  return TCL_ERROR;
   }

   if (theTclDomain->addElement(theEle) == false) {
	  opserr << "WARNING could not add element to the domain\n";
	  opserr << "fFrame2D element: " << eleId << endln;
	  delete theEle;
	  return TCL_ERROR;
   }
   return TCL_OK;
}
  
int
TclModelBuilder_addFeapShell(ClientData clientData, Tcl_Interp *interp, int argc,
TCL_Char **argv, Domain *theTclDomain, TclModelBuilder *theTclBuilder,
int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";
    return TCL_ERROR;
  }

  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  if (ndf != 3) {
    opserr << "WARNING - fShell eleTag? iNode? jNode? kNode? lNode? secTag? needs ndf=3\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  if ((argc - eleArgStart) < 5) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: element fShell eleTag? iNode? jNode? kNode? lNode? secTag?\n";

    return TCL_ERROR;
  }

  // get the id and end nodes 
  int eleId, iNode, jNode, kNode, lNode, secTag;

  if (Tcl_GetInt(interp, argv[1 + eleArgStart], &eleId) != TCL_OK) {
    opserr << "WARNING invalid fShell eleTag" << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2 + eleArgStart], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    opserr << "fShell element: " << eleId << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3 + eleArgStart], &jNode) != TCL_OK) {
    opserr << "WARNING invalid jNode\n";
    opserr << "fShell element: " << eleId << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[4 + eleArgStart], &kNode) != TCL_OK) {
    opserr << "WARNING invalid kNode\n";
    opserr << "fShell element: " << eleId << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[5 + eleArgStart], &lNode) != TCL_OK) {
    opserr << "WARNING invalid lNode\n";
    opserr << "fShell element: " << eleId << endln;
    return TCL_ERROR;
  }
  //if (Tcl_GetDouble(interp, argv[5 + eleArgStart], &E) != TCL_OK) {
  //  opserr << "WARNING invalid E\n";
  //  opserr << "fShell element: " << eleId << endln;
  //  return TCL_ERROR;
  //}
  if (Tcl_GetInt(interp, argv[6 + eleArgStart], &secTag) != TCL_OK) {
    opserr << "WARNING invalid secTag\n";
    opserr << "fShell element: " << eleId << endln;
    return TCL_ERROR;
  }

  SectionForceDeformation *theSection = theTclBuilder->getSection(secTag);

  // now create the truss and add it to the Domain
  fElmt03 *theEle = new fElmt03(eleId, iNode, jNode, kNode, lNode, secTag, 0.0, 1);

  if (theEle == 0) {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "fShell element: " << eleId << endln;
    return TCL_ERROR;
  }

  if (theTclDomain->addElement(theEle) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "fShell element: " << eleId << endln;
    delete theEle;
    return TCL_ERROR;
  }
  return TCL_OK;
}