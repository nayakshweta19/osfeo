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
                                                                        
// $Revision: 1.2 $
// $Date: 2008-07-03 18:03:49 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/element/TclRCFTLMBeamColumnCommand.cpp,v $
                                                                        
// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the implementation of the 
// TclModelBuilder_addDispBeamColumn() command. 

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <TclModelBuilder.h>

#include <RCFTLMBeamColumn3D.h>

#include <LobattoBeamIntegration.h>

#include <fstream>

//using std::ofstream;
//using std::ios;
//using std::endl;
using namespace std;

extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_addRCFTLMBeamColumn(ClientData clientData, Tcl_Interp *interp,  
				   int argc, 
				   TCL_Char **argv, 
				   Domain*theTclDomain,
				   TclModelBuilder *theTclBuilder)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }
  
  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();
  
  int ok = 0;
  if (ndm == 3 && ndf == 9)
    ok = 1;
  
  if (ok == 0) {
    opserr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
	 << " not compatible with RCFTLMBeamColumn element" << endln;
    return TCL_ERROR;
  }
  

  if (argc < 6) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: element RCFTLMBeamColumn eleTag? iNode? jNode? transfTag? ...\n";
    return TCL_ERROR;
  }

  int eleTag, iNode, jNode, transfTag;
  CrdTransf *theTransf3d = 0;
  Element *theElement = 0;


  if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
    opserr << "WARNING invalid RCFTLMBeamColumn eleTag" << endln;
    return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    opserr << "RCFTLMBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
    opserr << "WARNING invalid jNode\n";
    opserr << "RCFTLMBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }


  //
  // fmk UNDOCUMENTED FEATURE - 
  // all to take similar command to nonlinearBeamColumn & dispBeamColumn 
  // 

  if ((strcmp(argv[6],"Lobatto") != 0)) {

    int nIP, secTag;

    if (Tcl_GetInt(interp, argv[5], &nIP) != TCL_OK) {
      opserr << "WARNING invalid nIP\n";
      opserr << "RCFTLMBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[6], &secTag) != TCL_OK) {
      opserr << "WARNING invalid secTag\n";
      opserr << "RCFTLMBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }

    SectionForceDeformation *theSection = theTclBuilder->getSection(secTag);
 
    if ( theTclBuilder->getSection(secTag) == 0) {
      opserr << "WARNING section not found Cenk\n";
      opserr << "Section: " << secTag;
      opserr << "\nRCFTLMBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &transfTag) != TCL_OK) {
      opserr << "WARNING invalid transfTag\n";
      opserr << "RCFTLMBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    if (ndm == 3) {
      
      theTransf3d = OPS_GetCrdTransf(transfTag);
      
      if (theTransf3d == 0) {
	opserr << "WARNING transformation not found\n";
	opserr << "transformation: " << transfTag;
	opserr << "\nRCFTLMBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
    }

    
    //SectionForceDeformation **sections = new SectionForceDeformation *[nIP];

#ifdef COMPOSITE_DEBUG
	ofstream output;
    output.open("newton.dat",ios::app);

    output<<"\n TclRCFTLMBeamColumn \n"<<endl;
#endif

    RCFTAggregator **sections = new RCFTAggregator *[nIP];
    
    for (int i = 0; i < nIP; i++)
      sections[i] = (RCFTAggregator*) theSection;

    LobattoBeamIntegration beamIntegr;

    if (ndm == 3) {
	theElement = new RCFTLMBeamColumn3D(eleTag, iNode, jNode, nIP, sections,
					   beamIntegr, *theTransf3d);
    }

    delete [] sections;    
    if (theElement == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "RCFTLMBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    //output<<"\n TclDomain->addElement \n"<<endl;

    if (theTclDomain->addElement(theElement) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "RCFTLMBeamColumn element: " << eleTag << endln;
      delete theElement;
      return TCL_ERROR;
    }

    return TCL_OK;
  } 

  
  //
  // otherwise use correct format of command as found in current documentation
  //

  if (Tcl_GetInt(interp, argv[5], &transfTag) != TCL_OK) {
    opserr << "WARNING invalid transfTag\n";
    opserr << "RCFTLMBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  if (ndm == 3) {
    
    theTransf3d = OPS_GetCrdTransf(transfTag);
    
    if (theTransf3d == 0) {
      opserr << "WARNING transformation not found\n";
      opserr << "transformation: " << transfTag;
      opserr << "\nRCFTLMBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
  }
  
  if (strcmp(argv[6],"Lobatto") == 0) {
    int secTag, nIP;
    
    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: element RCFTLMBeamColumn eleTag? iNode? jNode? transfTag? Lobatto secTag? nIP?\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &secTag) != TCL_OK) {
      opserr << "WARNING invalid secTag\n";
      opserr << "RCFTLMBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[8], &nIP) != TCL_OK) {
      opserr << "WARNING invalid nIP\n";
      opserr << "RCFTLMBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    SectionForceDeformation *theSection = theTclBuilder->getSection(secTag);
    if (theSection == 0) {
      opserr << "WARNING section not found\n";
      opserr << "Section: " << secTag;
      opserr << "\nRCFTLMBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    RCFTAggregator **sections = new RCFTAggregator *[nIP];
    for (int i = 0; i < nIP; i++)
      sections[i] = (RCFTAggregator*) theSection;
    
    LobattoBeamIntegration beamIntegr;

    theElement = new RCFTLMBeamColumn3D(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf3d);
    delete [] sections;
  }
}
