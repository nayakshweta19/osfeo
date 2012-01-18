/* ****************************************************************** **
** OpenSees - Open System for Earthquake Engineering Simulation       **
** Pacific Earthquake Engineering Research Center                     **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited. See    **
** file 'COPYRIGHT' in main directory for information on usage and    **
** redistribution, and for a DISCLAIMER OF ALL WARRANTIES.            **
**                                                                    **
** Developed by:                                                      **
** Frank McKenna (fmckenna@ce.berkeley.edu)                           **
** Gregory L. Fenves (fenves@ce.berkeley.edu)                         **
** Filip C. Filippou (filippou@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */

// $Revision: 1. $
// $Date: 2008/07/20 19:20:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/pipe/TclPipe3Command.cpp,v $


// File: ~/element/TclPipe3Command.C
//
// Written: avytin
// Created: 09/08
// Revision: A
//
// Description: This file contains the implementation of the TclModelBuilder_Pipe3()
// command.
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include "Pipe3.h"
#include <TrussSection.h>
#include <TclModelBuilder.h>
#include <CorotTruss.h>
#include <CorotTrussSection.h>
#include <UniaxialMaterial.h>

extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_Pipe3(ClientData clientData, Tcl_Interp *interp, int argc,
	TCL_Char **argv, Domain*theTclDomain,
	TclModelBuilder *theTclBuilder,
	int eleArgStart)
{
	//make sure at least one other argument to contain type of system
	if (argc!=10){
		interp->result = "WARNING bad command - Pipe3 eleId iNode jNode matID Area c_3 Gamma d_c";
		return TCL_ERROR;
	}

	//get the id, x_loc, y_loc
	int trussId, iNode, jNode, matID;
	double A, C_3, Gamma, D_C;

	if (Tcl_GetInt(interp, argv[2], &trussId)!= TCL_OK){
		interp->result = "WARNING invalid eleId - Pipe3 eleId iNode jNode matID Area c_3 Gamma d_c";
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
		interp->result = "WARNING invalid iNode - Pipe3 eleId iNode jNode matID Area c_3 Gamma d_c";
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
		interp->result = "WARNING invalid jNode - Pipe3 eleId iNode jNode matID Area c_3 Gamma d_c";
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[5], &matID) != TCL_OK) {
		interp->result = "WARNING invalid matID - Pipe3 eleId iNode jNode matID Area c_3 Gamma d_c";
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &A) != TCL_OK) {
		interp->result = "WARNING invalid Area - Pipe3 eleId iNode jNode matID Area c_3 Gamma d_c";
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[7], &C_3) != TCL_OK) {
		interp->result = "WARNING invalid C_3 - Pipe3 eleId iNode jNode matID Area c_3 Gamma d_c";
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[8], &Gamma) != TCL_OK) {
		interp->result = "WARNING invalid Gamma - Pipe3 eleId iNode jNode matID Area c_3 Gamma d_c";
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[9], &D_C) != TCL_OK) {
		interp->result = "WARNING invalid d_c - Pipe3 eleId iNode jNode matID Area c_3 Gamma d_c";
		return TCL_ERROR;
	}

	UniaxialMaterial *theMaterial = OPS_getUniaxialMaterial(matID);

	if (theMaterial ==0) {
		opserr << "WARNING TclPipe3 - Pipe3 - no Material found with tag ";
		opserr << matID << endln;
		return TCL_ERROR;
	}

	//now create the truss and add it to the domain
	// MyTruss *theTruss = new MyTruss(trussId,iNode,jNode,*theMaterial,A,M);

	Element *theTruss = 0;
	theTruss=new Pipe3(trussId,iNode,jNode,*theMaterial,A,C_3,Gamma, D_C);

	if (theTruss==0) {
		opserr << "WARNING TclPipe3 - Pipe3 - ran out of memory for node ";
		opserr << trussId << endln;
		return TCL_ERROR;
	}
	if (theTclDomain->addElement(theTruss)==false) {
		delete theTruss;
		opserr << "WARNING TclPipe3 - Pipe3 - could not add Pipe3 to the domain";
		opserr << trussId << endln;
		return TCL_ERROR;
	}

	//Everything is OK
	return TCL_OK;
}