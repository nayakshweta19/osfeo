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

// $Revision: 4967 $
// $Date: 2012-08-13 13:39:44 +0800 (星期一, 13 八月 2012) $
// $URL: svn://opensees.berkeley.edu/usr/local/svn/OpenSees/trunk/SRC/element/adapter/TclActuatorCorotCommand.cpp $

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 09/07
// Revision: A
//
// Description: This file contains the function to parse the TCL input
// for the corotational actuator element.

#include <TclModelBuilder.h>

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <ActuatorCorot.h>


extern void printCommand(int argc, TCL_Char **argv);

int TclModelBuilder_addActuatorCorot(ClientData clientData, Tcl_Interp *interp,  int argc, 
    TCL_Char **argv, Domain*theTclDomain,
    TclModelBuilder *theTclBuilder, int eleArgStart)
{
    // ensure the destructor has not been called
    if (theTclBuilder == 0) {
        opserr << "WARNING builder has been destroyed - corotActuator\n";
        return TCL_ERROR;
    }
    
    // check the number of arguments is correct
    if ((argc-eleArgStart) < 6) {
        opserr << "WARNING insufficient arguments\n";
        printCommand(argc, argv);
        opserr << "Want: element corotActuator eleTag iNode jNode EA ipPort <-doRayleigh> <-rho rho>\n";
        return TCL_ERROR;
    }
    
    Element *theElement = 0;
    int ndm = theTclBuilder->getNDM();
    
    // get the id and end nodes 
    int tag, iNode, jNode;
    double EA;
    int ipPort;
    int doRayleigh = 0;
    double rho = 0.0;
    
    if (Tcl_GetInt(interp, argv[1+eleArgStart], &tag) != TCL_OK) {
        opserr << "WARNING invalid corotActuator eleTag" << endln;
        return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
        opserr << "WARNING invalid iNode\n";
        opserr << "corotActuator element: " << tag << endln;
        return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
        opserr << "WARNING invalid jNode\n";
        opserr << "corotActuator element: " << tag << endln;
        return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4+eleArgStart], &EA) != TCL_OK) {
        opserr << "WARNING invalid EA\n";
        opserr << "corotActuator element: " << tag << endln;
        return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[5+eleArgStart], &ipPort) != TCL_OK) {
        opserr << "WARNING invalid ipPort\n";
        opserr << "corotActuator element: " << tag << endln;
        return TCL_ERROR;
    }
    for (int i = 6+eleArgStart; i < argc; i++)  {
        if (strcmp(argv[i], "-doRayleigh") == 0)
            doRayleigh = 1;
    }
    for (int i = 6+eleArgStart; i < argc; i++) {
        if (i+1 < argc && strcmp(argv[i], "-rho") == 0) {
            if (Tcl_GetDouble(interp, argv[i+1], &rho) != TCL_OK) {
                opserr << "WARNING invalid rho\n";
                opserr << "corotActuator element: " << tag << endln;
                return TCL_ERROR;
            }
        }
    }
    
    // now create the corotActuator and add it to the Domain
    theElement = new ActuatorCorot(tag, ndm, iNode, jNode, EA, ipPort,
        doRayleigh, rho);
    
    if (theElement == 0) {
        opserr << "WARNING ran out of memory creating element\n";
        opserr << "corotActuator element: " << tag << endln;
        return TCL_ERROR;
    }
    
    if (theTclDomain->addElement(theElement) == false) {
        opserr << "WARNING could not add element to the domain\n";
        opserr << "corotActuator element: " << tag << endln;
        delete theElement;
        return TCL_ERROR;
    }
    
    // if get here we have sucessfully created the corotActuator and added it to the domain
    return TCL_OK;
}
