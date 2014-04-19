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

// $Source: /usr/local/cvs/OpenSees/SRC/element/Timoshenko/Timoshenko2d.cpp,v $
// $Revision: 1.2 $
// $Date: 2009/01/10 21:22:20 $

// Created: 09/09
// Created by: Li Ning (neallee@tju.edu.cn)
// Description: This file contains the class implementation of Timoshenko2d.
//              Make use of Neddy(1997) Interdependent Integration Element 
//              procecess and fiber section model.

// Reference: LI Ning, LI Zhong-Xian, XIE Li-Li. A Fiber-Section Model Based
//            Timoshenko Beam Element Using Shear-Bending Interdependent Shape 
//            Function. Earthquake Engineering & Engineering Vibration. 2013, 
//            12(3): 421-432.

#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <elementAPI.h>
#include <Information.h>

#include "TimoshenkoBeam2d.h"
#include "Timoshenko2d01.h"
#include "Timoshenko2d02.h"
#include "Timoshenko2d03.h"
#include "Timoshenko2d04.h"
#include "Timoshenko2d.h"
#include "Timoshenko3d.h"
#include "Timoshenko3d01.h"
#include "Timoshenko3d04.h"

#include <LobattoBeamIntegration.h>
#include <LegendreBeamIntegration.h>
#include <RadauBeamIntegration.h>
#include <NewtonCotesBeamIntegration.h>
#include <UserDefinedBeamIntegration.h>

#include <HingeMidpointBeamIntegration.h>
#include <HingeEndpointBeamIntegration.h>
#include <HingeRadauBeamIntegration.h>
#include <HingeRadauTwoBeamIntegration.h>
#include <UserDefinedHingeIntegration.h>
#include <DistHingeIntegration.h>
#include <RegularizedHingeIntegration.h>

#include <TrapezoidalBeamIntegration.h>
#include <FixedLocationBeamIntegration.h>
#include <LowOrderBeamIntegration.h>
#include <MidDistanceBeamIntegration.h>

#include <CrdTransf.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, TCL_Char **argv);

// element Timoshenko2d01 eleTag? iNode? jNode? secTag? transfTag? C1?
int
TclModelBuilder_addTimoshenko2d01(ClientData clientData, Tcl_Interp *interp,  
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
	if (ndm == 2 && ndf == 3)
	  ok = 1;

	if (ok == 0) {
	  opserr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
	    << " not compatible with Timoshenko2d01 element" << endln;
	  return TCL_ERROR;
	}

	if (argc < 8) {			//8
	  opserr << "WARNING insufficient arguments\n";
	  printCommand(argc, argv);
	  opserr << "Want: element Timoshenko2d01 eleTag? iNode? jNode? secTag? transfTag? C1?\n";
	  return TCL_ERROR;
	}

	// get the id and end nodes 
	int eleTag, iNode, jNode, transfTag, secTag;
	double C1;
	int argi = 2;

	if (Tcl_GetInt(interp, argv[argi++], &eleTag) != TCL_OK) {
	  opserr << "WARNING invalid Timoshenko2d01 eleTag" << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &iNode) != TCL_OK) {
	  opserr << "WARNING invalid iNode ";
	  opserr << "Timoshenko2d01 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &jNode) != TCL_OK) {
	  opserr << "WARNING invalid jNode ";
	  opserr << "Timoshenko2d01 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &secTag) != TCL_OK) {
	  opserr << "WARNING invalid secTag - element Timoshenko2d01 eleTag? iNode? jNode? secTag? transfTag? C1?\n";
	  return TCL_ERROR;
	}

	if (argi >= argc || Tcl_GetInt(interp, argv[argi++], &transfTag) != TCL_OK) {
	  opserr << "WARNING invalid transfTag? - element Timoshenko2d01 eleTag? iNode? jNode? secTag? transfTag? C1?\n";
	  return TCL_ERROR;
	}

	if (argi >= argc || Tcl_GetDouble(interp, argv[argi++], &C1) != TCL_OK) {
	  opserr << "WARNING invalid transfTag? - element Timoshenko2d01 eleTag? iNode? jNode? secTag? transfTag? C1?\n";
	  return TCL_ERROR;
	}

	double massDens = 0.0;

	while (argi != argc) {
	  if (strcmp(argv[argi++],"-mass") == 0 && argi < argc) {
	  	if (Tcl_GetDouble(interp, argv[argi++], &massDens) != TCL_OK) {
	  	  opserr << "WARNING invalid massDens - element Timoshenko2d01 eleTag? iNode? jNode? "
			  << "secTag? transfTag? C1? -mass massDens?\n";
	  	  return TCL_ERROR;
	  	}
	  }
	}

	SectionForceDeformation **sections = new SectionForceDeformation* [1];

	if (!sections) {
	  opserr << "WARNING TclTimoshenkoBeamCommand - Insufficient memory to create sections\n";
	  return TCL_ERROR;
	}

	SectionForceDeformation *theSection = theTclBuilder->getSection(secTag);

	if (theSection == 0) {
	  opserr << "WARNING TclTimoshenkoBeamCommand - no Section found with tag ";
	  opserr << secTag << endln;
	  delete [] sections;
	  return TCL_ERROR;
	}

	sections[0] = theSection;
	

	CrdTransf *theTransf2d = 0;
	//CrdTransf *theTransf3d = 0;
	Element *theElement = 0;

	BeamIntegration *beamIntegr = 0;
	beamIntegr = new LegendreBeamIntegration();

	if (ndm == 2) {
	  theTransf2d = OPS_GetCrdTransf(transfTag);
	  
	  if (theTransf2d == 0) {
	  	opserr << "WARNING transformation not found\n";
	  	opserr << "transformation: " << transfTag;
	  	opserr << argv[1] << " element: " << eleTag << endln;
	  	return TCL_ERROR;
	  }
	  
      // now create the Timoshenko2d01 and add it to the Domain
	  theElement = new Timoshenko2d01(eleTag,iNode,jNode,sections,*theTransf2d,*beamIntegr,C1,massDens);
	  
	  delete [] sections;
	}

	if (theElement == 0) {
	  opserr << "WARNING ran out of memory creating element\n";
	  opserr << "Timoshenko2d01 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (theTclDomain->addElement(theElement) == false) {
	  opserr << "WARNING could not add element to the domain\n";
	  opserr << "Timoshenko2d01 element: " << eleTag << endln;
	  delete theElement;
	  return TCL_ERROR;
	}

	// if get here we have successfully created the element and added it to the domain
	return TCL_OK;
}

//element Timoshenko2d02 eleTag? iNode? jNode? nIP? secTag? transfTag?
int
TclModelBuilder_addTimoshenko2d02(ClientData clientData, Tcl_Interp *interp,  
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
	if (ndm == 2 && ndf == 3)
	  ok = 1;

	if (ok == 0) {
	  opserr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
	    << " not compatible with Timoshenko2d02 element" << endln;
	  return TCL_ERROR;
	}

	if (argc < 8) {			//8
	  opserr << "WARNING insufficient arguments\n";
	  printCommand(argc, argv);
	  opserr << "Want: element Timoshenko2d02 eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	  return TCL_ERROR;
	}

	// get the id and end nodes 
	int eleTag, iNode, jNode, nIP, secTag, transfTag;
	int argi = 2;

	if (Tcl_GetInt(interp, argv[argi++], &eleTag) != TCL_OK) {
	  opserr << "WARNING invalid Timoshenko2d02 eleTag" << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &iNode) != TCL_OK) {
	  opserr << "WARNING invalid iNode ";
	  opserr << "Timoshenko2d02 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &jNode) != TCL_OK) {
	  opserr << "WARNING invalid jNode ";
	  opserr << "Timoshenko2d02 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &nIP) != TCL_OK) {
		opserr << "WARNING invalid nIP ";
		opserr << "Timoshenko2d02 element: " << eleTag << endln;
		return TCL_ERROR;
	}  

	if (Tcl_GetInt(interp, argv[argi++], &secTag) != TCL_OK) {
	  opserr << "WARNING invalid secTag - element Timoshenko2d02 eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	  return TCL_ERROR;
	}

	if (argi >= argc || Tcl_GetInt(interp, argv[argi++], &transfTag) != TCL_OK) {
	  opserr << "WARNING invalid transfTag? - element Timoshenko2d02 eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	  return TCL_ERROR;
	}

	double massDens = 0.0;

	while (argi != argc) {
	  if (strcmp(argv[argi++],"-mass") == 0 && argi < argc) {
	  	if (Tcl_GetDouble(interp, argv[argi++], &massDens) != TCL_OK) {
	  	  opserr << "WARNING invalid massDens - element Timoshenko2d02 eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	  	  return TCL_ERROR;
	  	}
	  }
	}

	SectionForceDeformation **sections = new SectionForceDeformation* [nIP];

	if (!sections) {
	  opserr << "WARNING TclTimoshenkoBeamCommand - Insufficient memory to create sections\n";
	  return TCL_ERROR;
	}

	for (int j=0; j<nIP; j++) {
	  SectionForceDeformation *theSection = theTclBuilder->getSection(secTag);
	  
	  if (theSection == 0) {
	    opserr << "WARNING TclTimoshenkoBeamCommand - no Section found with tag ";
	    opserr << secTag << endln;
	    delete [] sections;
	    return TCL_ERROR;
	  }
	  
	  sections[j] = theSection;
	}

	CrdTransf *theTransf2d = 0;
	CrdTransf *theTransf3d = 0;
	Element *theElement = 0;

	BeamIntegration *beamIntegr = 0;
	beamIntegr = new LegendreBeamIntegration();

	if (ndm == 2) {
	  theTransf2d = OPS_GetCrdTransf(transfTag);
	  
	  if (theTransf2d == 0) {
	  	opserr << "WARNING transformation not found\n";
	  	opserr << "transformation: " << transfTag;
	  	opserr << argv[1] << " element: " << eleTag << endln;
	  	return TCL_ERROR;
	  }
	  
      // now create the Timoshenko2d02 and add it to the Domain
	  theElement = new Timoshenko2d02(eleTag, iNode, jNode, nIP, sections,
		                              *theTransf2d,*beamIntegr,massDens);
	  delete [] sections;
	}

	if (theElement == 0) {
	  opserr << "WARNING ran out of memory creating element\n";
	  opserr << "Timoshenko2d02 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (theTclDomain->addElement(theElement) == false) {
	  opserr << "WARNING could not add element to the domain\n";
	  opserr << "Timoshenko2d02 element: " << eleTag << endln;
	  delete theElement;
	  return TCL_ERROR;
	}

	// if get here we have successfully created the element and added it to the domain
	return TCL_OK;
}

//element Timoshenko2d03 eleTag? iNode? jNode? $beamSecTag $transfTag C?
int
TclModelBuilder_addTimoshenko2d03(ClientData clientData, Tcl_Interp *interp,  
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
	if (ndm == 2 && ndf == 3)
	  ok = 1;

	if (ok == 0) {
	  opserr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
	    << " not compatible with Timoshenko2d03 element" << endln;
	  return TCL_ERROR;
	}

	if (argc < 8) {			//8
	  opserr << "WARNING insufficient arguments\n";
	  printCommand(argc, argv);
	  opserr << "Want: element Timoshenko2d03 eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	  return TCL_ERROR;
	}

	// get the id and end nodes 
	int eleTag, iNode, jNode, nIP, transfTag, secTag;
	int argi = 2;

	if (Tcl_GetInt(interp, argv[argi++], &eleTag) != TCL_OK) {
	  opserr << "WARNING invalid Timoshenko2d03 eleTag" << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &iNode) != TCL_OK) {
	  opserr << "WARNING invalid iNode ";
	  opserr << "Timoshenko2d03 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &jNode) != TCL_OK) {
	  opserr << "WARNING invalid jNode ";
	  opserr << "Timoshenko2d03 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &nIP) != TCL_OK) {
	  opserr << "WARNING invalid nIP? ";
	  opserr << "Timoshenko2d03 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &secTag) != TCL_OK) {
	  opserr << "WARNING invalid secTag - element Timoshenko2d03 eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	  return TCL_ERROR;
	}

	if (argi >= argc || Tcl_GetInt(interp, argv[argi++], &transfTag) != TCL_OK) {
	  opserr << "WARNING invalid transfTag? - element Timoshenko2d03 eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	  return TCL_ERROR;
	}

	double massDens = 0.0;

	while (argi != argc) {
	  if (strcmp(argv[argi++],"-mass") == 0 && argi < argc) {
	  	if (Tcl_GetDouble(interp, argv[argi++], &massDens) != TCL_OK) {
	  	  opserr << "WARNING invalid massDens - element Timoshenko2d03 eleTag? iNode? jNode? "
			  << "nIP? secTag? transfTag? -mass massDens?\n";
	  	  return TCL_ERROR;
	  	}
	  }
	}

	SectionForceDeformation **sections = new SectionForceDeformation* [nIP];

	if (!sections) {
	  opserr << "WARNING TclTimoshenkoBeamCommand - Insufficient memory to create sections\n";
	  return TCL_ERROR;
	}

	SectionForceDeformation *theSection = theTclBuilder->getSection(secTag);

	if (theSection == 0) {
	  opserr << "WARNING TclTimoshenkoBeamCommand - no Section found with tag ";
	  opserr << secTag << endln;
	  delete [] sections;
	  return TCL_ERROR;
	}

	for (int i = 0; i < nIP; i++)
      sections[i] = theSection;
	
	CrdTransf *theTransf2d = 0;
	//CrdTransf *theTransf3d = 0;
	Element *theElement = 0;

	BeamIntegration *beamIntegr = 0;
	beamIntegr = new LobattoBeamIntegration();

	if (ndm == 2) {
	  theTransf2d = OPS_GetCrdTransf(transfTag);
	  
	  if (theTransf2d == 0) {
	  	opserr << "WARNING transformation not found\n";
	  	opserr << "transformation: " << transfTag;
	  	opserr << argv[1] << " element: " << eleTag << endln;
	  	return TCL_ERROR;
	  }
	  
      // now create the Timoshenko2d03 and add it to the Domain
	  theElement = new Timoshenko2d03(eleTag,iNode,jNode,nIP,sections,*theTransf2d,*beamIntegr,massDens);
	  
	  delete [] sections;
	}

	if (theElement == 0) {
	  opserr << "WARNING ran out of memory creating element\n";
	  opserr << "Timoshenko2d03 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (theTclDomain->addElement(theElement) == false) {
	  opserr << "WARNING could not add element to the domain\n";
	  opserr << "Timoshenko2d03 element: " << eleTag << endln;
	  delete theElement;
	  return TCL_ERROR;
	}

	// if get here we have successfully created the element and added it to the domain
	return TCL_OK;
}

//element Timoshenko2d04 eleTag? iNode? jNode? nIP? secTag? transfTag?
int
TclModelBuilder_addTimoshenko2d04(ClientData clientData, Tcl_Interp *interp,  
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
	if (ndm == 2 && ndf == 3)
	  ok = 1;

	if (ok == 0) {
	  opserr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
	    << " not compatible with Timoshenko2d04 element" << endln;
	  return TCL_ERROR;
	}

	if (argc < 8) {			//8
	  opserr << "WARNING insufficient arguments\n";
	  printCommand(argc, argv);
	  opserr << "Want: element Timoshenko2d04 eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	  return TCL_ERROR;
	}

	// get the id and end nodes 
	int eleTag, iNode, jNode, nIP, secTag, transfTag;
	int numData;
    // Variables to retrieve input
    int iData[10];
    double dData[10];
    int sDataLength = 40;
    char *sData  = new char[sDataLength]; 
	char *sData2  = new char[sDataLength];

    numData = 6;
    if (OPS_GetIntInput(&numData, iData) != 0) {
      opserr << "WARNING invalid element data - Timoshenko2d04\n";
      return 0;
    }

    eleTag = iData[0];
    iNode = iData[1];
    jNode = iData[2];
    nIP = iData[3];
    secTag = iData[4];
    transfTag = iData[5];

    SectionForceDeformation **sections = new SectionForceDeformation *[nIP];
	if (!sections) {
	  opserr << "WARNING TclTimoshenkoBeamCommand - Insufficient memory to create sections\n";
	  return TCL_ERROR;
	}

    for (int i = 0; i < nIP; i++) {
      // Get the section
      SectionForceDeformation *theSection = OPS_GetSectionForceDeformation(secTag);
      if (theSection == 0) {
	    opserr << "WARNING TclTimoshenkoBeamCommand - no Section found with tag ";
	    opserr << secTag << endln;
	    delete [] sections;
	    return TCL_ERROR;
      }
    
      sections[i] = theSection;
	}
    // Get the coordinate transformation
    CrdTransf *theTransf = OPS_GetCrdTransfPtr(transfTag);
    if (theTransf == 0) {
      opserr << "WARNING geometric transformation with tag " << transfTag << "not found for element " << eleTag << endln;
      return 0;
    }

	double massDens = 0.0;
	double shearCF = 1.0;
	BeamIntegration *beamIntegr = 0;

    // Loop through remaining arguments to get optional input
    while ( OPS_GetNumRemainingInputArgs() > 0 ) {
      if ( OPS_GetString(sData, sDataLength) != 0 ) {
        opserr << "WARNING invalid input";
        return 0;
      }
    
      if ( strcmp(sData,"-mass") == 0 ) {
        numData = 1;
        if (OPS_GetDoubleInput(&numData, dData) != 0) {
          opserr << "WARNING invalid input, want: -mass $massDens \n";
          return 0;
        }
        massDens = dData[0];
    
      } else if ( strcmp(sData,"-integration") == 0 ) {
        if ( OPS_GetString(sData2, sDataLength) != 0 ) {
          opserr << "WARNING invalid input, want: -integration $intType";
          return 0;
        }
    
        if (strcmp(sData2,"Lobatto") == 0) {
          beamIntegr = new LobattoBeamIntegration();
        } else if (strcmp(sData2,"Legendre") == 0) {
          beamIntegr = new LegendreBeamIntegration();
        } else if (strcmp(sData2,"Radau") == 0) {
          beamIntegr = new RadauBeamIntegration();
        } else if (strcmp(sData2,"NewtonCotes") == 0) {
          beamIntegr = new NewtonCotesBeamIntegration();
        } else if (strcmp(sData2,"Trapezoidal") == 0) {
          beamIntegr = new TrapezoidalBeamIntegration();
        } else if (strcmp(sData2,"RegularizedLobatto") == 0 || strcmp(sData2,"RegLobatto") == 0) {
          numData = 4;
          if (OPS_GetDoubleInput(&numData, dData) != 0) {
            opserr << "WARNING invalid input, want: -integration RegularizedLobatto $lpI $lpJ $zetaI $zetaJ \n";
            return 0;
          }
          BeamIntegration *otherBeamInt = 0;
          otherBeamInt = new LobattoBeamIntegration();
          beamIntegr = new RegularizedHingeIntegration(*otherBeamInt, dData[0], dData[1], dData[2], dData[3]);
          if (otherBeamInt != 0)
            delete otherBeamInt;

        } else {
        opserr << "WARNING invalid integration type, element: " << eleTag;
          return 0;
        }
    
      } else if ( strcmp(sData, "-shearCF") == 0 ) {
		numData = 1;
		if ( OPS_GetDoubleInput(&numData, dData) != 0 ) {
		  opserr << "WARNING invalid input, want: -shearCF $shearCorrectFactor";
		  return 0;
		}
		shearCF = dData[0];

	  }else {
        opserr << "WARNING unknown option " << sData << "\n";
      }
    }
    
    // Set the beam integration object if not in options
    if (beamIntegr == 0) {
      beamIntegr = new LobattoBeamIntegration();
    }

    theTransf = OPS_GetCrdTransf(transfTag);
    
    if (theTransf == 0) {
    	opserr << "WARNING transformation not found\n";
    	opserr << "transformation: " << transfTag;
    	opserr << argv[1] << " element: " << eleTag << endln;
    	return TCL_ERROR;
    }
    
    // now create the Timoshenko2d04 and add it to the Domain
    Element *theElement = new Timoshenko2d04(eleTag, iNode, jNode, nIP, sections,
                                  *theTransf,*beamIntegr,massDens,shearCF);
    
	if (theElement == 0) {
	  opserr << "WARNING ran out of memory creating element\n";
	  opserr << "Timoshenko2d04 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (theTclDomain->addElement(theElement) == false) {
	  opserr << "WARNING could not add element to the domain\n";
	  opserr << "Timoshenko2d04 element: " << eleTag << endln;
	  delete theElement;
	  return TCL_ERROR;
	}

	// if get here we have successfully created the element and added it to the domain
	return TCL_OK;
}

//element Timoshenko2d eleTag? iNode? jNode? nIP? secTag? transfTag?
int
TclModelBuilder_addTimoshenko2d(ClientData clientData, Tcl_Interp *interp,  
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
	if (ndm == 2 && ndf == 3)
	  ok = 1;

	if (ok == 0) {
	  opserr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
	    << " not compatible with Timoshenko2d04 element" << endln;
	  return TCL_ERROR;
	}

	if (argc < 8) {			//8
	  opserr << "WARNING insufficient arguments\n";
	  printCommand(argc, argv);
	  opserr << "Want: element Timoshenko2d05 eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	  return TCL_ERROR;
	}

	// get the id and end nodes 
	int eleTag, iNode, jNode, nIP, secTag, transfTag;
	int numData;
    // Variables to retrieve input
    int iData[10];
    double dData[10];
    int sDataLength = 40;
    char *sData  = new char[sDataLength]; 
	char *sData2  = new char[sDataLength];

    numData = 6;
    if (OPS_GetIntInput(&numData, iData) != 0) {
      opserr << "WARNING invalid element data - Timoshenko2d05\n";
      return 0;
    }

    eleTag = iData[0];
    iNode = iData[1];
    jNode = iData[2];
    nIP = iData[3];
    secTag = iData[4];
    transfTag = iData[5];

    // Get the section
    SectionForceDeformation *theSection = OPS_GetSectionForceDeformation(secTag);
    if (theSection == 0) {
      opserr << "WARNING section with tag " << secTag << "not found for element " << eleTag << endln;
      return 0;
    }
    
    SectionForceDeformation **sections = new SectionForceDeformation *[nIP];
    for (int i = 0; i < nIP; i++)
      sections[i] = theSection;

    // Get the coordinate transformation
    CrdTransf *theTransf = OPS_GetCrdTransfPtr(transfTag);
    if (theTransf == 0) {
      opserr << "WARNING geometric transformation with tag " << transfTag << "not found for element " << eleTag << endln;
      return 0;
    }

	double massDens = 0.0;
	double shearCF = 1.0;
	int noIter = 0;
	BeamIntegration *beamIntegr = 0;

    // Loop through remaining arguments to get optional input
    while ( OPS_GetNumRemainingInputArgs() > 0 ) {
      if ( OPS_GetString(sData, sDataLength) != 0 ) {
        opserr << "WARNING invalid input";
        return 0;
      }
    
      if ( strcmp(sData,"-mass") == 0 ) {
        numData = 1;
        if (OPS_GetDoubleInput(&numData, dData) != 0) {
          opserr << "WARNING invalid input, want: -mass $massDens \n";
          return 0;
        }
        massDens = dData[0];
    
      } else if ( strcmp(sData,"-integration") == 0 ) {
        if ( OPS_GetString(sData2, sDataLength) != 0 ) {
          opserr << "WARNING invalid input, want: -integration $intType";
          return 0;
        }
    
        if (strcmp(sData2,"Lobatto") == 0) {
          beamIntegr = new LobattoBeamIntegration();
        } else if (strcmp(sData2,"Legendre") == 0) {
          beamIntegr = new LegendreBeamIntegration();
        } else if (strcmp(sData2,"Radau") == 0) {
          beamIntegr = new RadauBeamIntegration();
        } else if (strcmp(sData2,"NewtonCotes") == 0) {
          beamIntegr = new NewtonCotesBeamIntegration();
        } else if (strcmp(sData2,"Trapezoidal") == 0) {
          beamIntegr = new TrapezoidalBeamIntegration();
        } else if (strcmp(sData2,"RegularizedLobatto") == 0 || strcmp(sData2,"RegLobatto") == 0) {
          numData = 4;
          if (OPS_GetDoubleInput(&numData, dData) != 0) {
            opserr << "WARNING invalid input, want: -integration RegularizedLobatto $lpI $lpJ $zetaI $zetaJ \n";
            return 0;
          }
          BeamIntegration *otherBeamInt = 0;
          otherBeamInt = new LobattoBeamIntegration();
          beamIntegr = new RegularizedHingeIntegration(*otherBeamInt, dData[0], dData[1], dData[2], dData[3]);
          if (otherBeamInt != 0)
            delete otherBeamInt;

        } else {
        opserr << "WARNING invalid integration type, element: " << eleTag;
          return 0;
        }
    
      } else if ( strcmp(sData, "-shearCF") == 0 ) {
		numData = 1;
		if ( OPS_GetDoubleInput(&numData, dData) != 0 ) {
		  opserr << "WARNING invalid input, want: -shearCF $shearCorrectFactor";
		  return 0;
		}
		shearCF = dData[0];

	  } else if ( strcmp(sData, "-noIterative") == 0 ) {
	    noIter = 1;
	  
	  } else {
        opserr << "WARNING unknown option " << sData << "\n";
      }
    }
    
    // Set the beam integration object if not in options
    if (beamIntegr == 0) {
      beamIntegr = new LobattoBeamIntegration();
    }

    theTransf = OPS_GetCrdTransf(transfTag);
    
    if (theTransf == 0) {
    	opserr << "WARNING transformation not found\n";
    	opserr << "transformation: " << transfTag;
    	opserr << argv[1] << " element: " << eleTag << endln;
    	return TCL_ERROR;
    }
    
    // now create the Timoshenko2d and add it to the Domain
    Element *theElement = new Timoshenko2d(eleTag, iNode, jNode, nIP, sections,
                                  *theTransf, *beamIntegr, massDens, shearCF, noIter);
    
	if (theElement == 0) {
	  opserr << "WARNING ran out of memory creating element\n";
	  opserr << "Timoshenko2d05 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (theTclDomain->addElement(theElement) == false) {
	  opserr << "WARNING could not add element to the domain\n";
	  opserr << "Timoshenko2d05 element: " << eleTag << endln;
	  delete theElement;
	  return TCL_ERROR;
	}

	// if get here we have successfully created the element and added it to the domain
	return TCL_OK;
}


//element Timoshenko3d01 eleTag? iNode? jNode? $beamSecTag $transfTag C1?
int
TclModelBuilder_addTimoshenko3d01(ClientData clientData, Tcl_Interp *interp,  
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
	if (ndm == 3 && ndf == 6)
	  ok = 1;

	if (ok == 0) {
	  opserr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
	    << " not compatible with Timoshenko3d01 element" << endln;
	  return TCL_ERROR;
	}

	if (argc < 8) {			//8
	  opserr << "WARNING insufficient arguments\n";
	  printCommand(argc, argv);
	  opserr << "Want: element Timoshenko3d01 eleTag? iNode? jNode? secTag? transfTag? C1?\n";
	  return TCL_ERROR;
	}

	// get the id and end nodes 
	int eleTag, iNode, jNode, transfTag, secTag;
	double C1;
	int argi = 2;

	if (Tcl_GetInt(interp, argv[argi++], &eleTag) != TCL_OK) {
	  opserr << "WARNING invalid Timoshenko3d01 eleTag" << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &iNode) != TCL_OK) {
	  opserr << "WARNING invalid iNode ";
	  opserr << "Timoshenko3d01 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &jNode) != TCL_OK) {
	  opserr << "WARNING invalid jNode ";
	  opserr << "Timoshenko3d01 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &secTag) != TCL_OK) {
	  opserr << "WARNING invalid secTag - element Timoshenko3d01 eleTag? iNode? jNode? secTag? transfTag? C1?\n";
	  return TCL_ERROR;
	}

	if (argi >= argc || Tcl_GetInt(interp, argv[argi++], &transfTag) != TCL_OK) {
	  opserr << "WARNING invalid transfTag? - element Timoshenko3d01 eleTag? iNode? jNode? secTag? transfTag? C1?\n";
	  return TCL_ERROR;
	}

	if (argi >= argc || Tcl_GetDouble(interp, argv[argi++], &C1) != TCL_OK) {
	  opserr << "WARNING invalid transfTag? - element Timoshenko3d01 eleTag? iNode? jNode? secTag? transfTag? C1?\n";
	  return TCL_ERROR;
	}

	double massDens = 0.0;

	while (argi != argc) {
	  if (strcmp(argv[argi++],"-mass") == 0 && argi < argc) {
	  	if (Tcl_GetDouble(interp, argv[argi++], &massDens) != TCL_OK) {
	  	  opserr << "WARNING invalid massDens - element Timoshenko3d01 eleTag? iNode? jNode? "
			  << "secTag? transfTag? C1? -mass massDens?\n";
	  	  return TCL_ERROR;
	  	}
	  }
	}

	SectionForceDeformation **sections = new SectionForceDeformation* [1];

	if (!sections) {
	  opserr << "WARNING TclTimoshenkoBeamCommand - Insufficient memory to create sections\n";
	  return TCL_ERROR;
	}

	SectionForceDeformation *theSection = theTclBuilder->getSection(secTag);

	if (theSection == 0) {
	  opserr << "WARNING TclTimoshenkoBeamCommand - no Section found with tag ";
	  opserr << secTag << endln;
	  delete [] sections;
	  return TCL_ERROR;
	}

	sections[0] = theSection;
	

	//CrdTransf *theTransf2d = 0;
	CrdTransf *theTransf3d = 0;
	Element *theElement = 0;

	BeamIntegration *beamIntegr = 0;
	beamIntegr = new LegendreBeamIntegration();

	if (ndm == 3) {
	  theTransf3d = OPS_GetCrdTransf(transfTag);
	  
	  if (theTransf3d == 0) {
	  	opserr << "WARNING transformation not found\n";
	  	opserr << "transformation: " << transfTag;
	  	opserr << argv[1] << " element: " << eleTag << endln;
	  	return TCL_ERROR;
	  }
	  
      // now create the Timoshenko3d01 and add it to the Domain
	  theElement = new Timoshenko3d01(eleTag,iNode,jNode,sections,*theTransf3d,*beamIntegr,C1,massDens);
	  
	  delete [] sections;
	}

	if (theElement == 0) {
	  opserr << "WARNING ran out of memory creating element\n";
	  opserr << "Timoshenko3d01 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (theTclDomain->addElement(theElement) == false) {
	  opserr << "WARNING could not add element to the domain\n";
	  opserr << "Timoshenko3d01 element: " << eleTag << endln;
	  delete theElement;
	  return TCL_ERROR;
	}

	// if get here we have successfully created the element and added it to the domain
	return TCL_OK;
}

//element Timoshenko3d04 eleTag? iNode? jNode? nIP? secTag? transfTag?
int
TclModelBuilder_addTimoshenko3d04(ClientData clientData, Tcl_Interp *interp,  
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
	if (ndm == 3 && ndf == 6)
	  ok = 1;

	if (ok == 0) {
	  opserr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
	    << " not compatible with Timoshenko3d04 element" << endln;
	  return TCL_ERROR;
	}

	if (argc < 8) {			//8
	  opserr << "WARNING insufficient arguments\n";
	  printCommand(argc, argv);
	  opserr << "Want: element Timoshenko3d04 eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	  return TCL_ERROR;
	}

	// get the id and end nodes 
	int eleTag, iNode, jNode, nIP, transfTag, secTag;
	int argi = 2;

	if (Tcl_GetInt(interp, argv[argi++], &eleTag) != TCL_OK) {
	  opserr << "WARNING invalid Timoshenko3d04 eleTag" << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &iNode) != TCL_OK) {
	  opserr << "WARNING invalid iNode ";
	  opserr << "Timoshenko3d04 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &jNode) != TCL_OK) {
	  opserr << "WARNING invalid jNode ";
	  opserr << "Timoshenko3d04 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &nIP) != TCL_OK) {
	  opserr << "WARNING invalid nIP ";
	  opserr << "Timoshenko3d04 element: " << eleTag << endln;
	  return TCL_ERROR;
	}  

	if (Tcl_GetInt(interp, argv[argi++], &secTag) != TCL_OK) {
	  opserr << "WARNING invalid secTag - element Timoshenko3d04 eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	  return TCL_ERROR;
	}

	if (argi >= argc || Tcl_GetInt(interp, argv[argi++], &transfTag) != TCL_OK) {
	  opserr << "WARNING invalid transfTag? - element Timoshenko3d04 eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	  return TCL_ERROR;
	}

	double massDens = 0.0;

	while (argi != argc) {
	  if (strcmp(argv[argi++],"-mass") == 0 && argi < argc) {
	  	if (Tcl_GetDouble(interp, argv[argi++], &massDens) != TCL_OK) {
	  	  opserr << "WARNING invalid massDens - element Timoshenko3d04 eleTag? iNode? jNode?"
			  << " nIP? secTag? transfTag? -mass massDens?\n";
	  	  return TCL_ERROR;
	  	}
	  }
	}

	SectionForceDeformation **sections = new SectionForceDeformation* [nIP];

	if (!sections) {
	  opserr << "WARNING TclTimoshenkoBeamCommand - Insufficient memory to create sections\n";
	  return TCL_ERROR;
	}

	for (int j=0; j<nIP; j++) {
	  SectionForceDeformation *theSection = theTclBuilder->getSection(secTag);
	  
	  if (theSection == 0) {
	    opserr << "WARNING TclTimoshenkoBeamCommand - no Section found with tag ";
	    opserr << secTag << endln;
	    delete [] sections;
	    return TCL_ERROR;
	  }
	  
	  sections[j] = theSection;
	}
	
	CrdTransf *theTransf3d = 0;
	Element *theElement = 0;

	BeamIntegration *beamIntegr = 0;
	beamIntegr = new LobattoBeamIntegration();

	if (ndm == 3) {
	  theTransf3d = OPS_GetCrdTransf(transfTag);
	  
	  if (theTransf3d == 0) {
	  	opserr << "WARNING transformation not found\n";
	  	opserr << "transformation: " << transfTag;
	  	opserr << argv[1] << " element: " << eleTag << endln;
	  	return TCL_ERROR;
	  }
	  
      // now create the Timoshenko3d04 and add it to the Domain
	  theElement = new Timoshenko3d04(eleTag,iNode,jNode,nIP, sections,*theTransf3d,*beamIntegr,massDens);
	  
	  delete [] sections;
	}

	if (theElement == 0) {
	  opserr << "WARNING ran out of memory creating element\n";
	  opserr << "Timoshenko3d01 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (theTclDomain->addElement(theElement) == false) {
	  opserr << "WARNING could not add element to the domain\n";
	  opserr << "Timoshenko3d01 element: " << eleTag << endln;
	  delete theElement;
	  return TCL_ERROR;
	}

	// if get here we have successfully created the element and added it to the domain
	return TCL_OK;
}

//element Timoshenko3d05 eleTag? iNode? jNode? nIP? secTag? transfTag?
int
TclModelBuilder_addTimoshenko3d(ClientData clientData, Tcl_Interp *interp,  
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
	if (ndm == 3 && ndf == 6)
	  ok = 1;

	if (ok == 0) {
	  opserr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
	    << " not compatible with Timoshenko3d element" << endln;
	  return TCL_ERROR;
	}

	if (argc < 8) {			//8
	  opserr << "WARNING insufficient arguments\n";
	  printCommand(argc, argv);
	  opserr << "Want: element Timoshenko3d eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	  return TCL_ERROR;
	}

	// get the id and end nodes 
	int eleTag, iNode, jNode, nIP, transfTag, secTag;
	int argi = 2;

	if (Tcl_GetInt(interp, argv[argi++], &eleTag) != TCL_OK) {
	  opserr << "WARNING invalid Timoshenko3d eleTag" << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &iNode) != TCL_OK) {
	  opserr << "WARNING invalid iNode ";
	  opserr << "Timoshenko3d element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &jNode) != TCL_OK) {
	  opserr << "WARNING invalid jNode ";
	  opserr << "Timoshenko3d element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &nIP) != TCL_OK) {
	  opserr << "WARNING invalid nIP ";
	  opserr << "Timoshenko3d element: " << eleTag << endln;
	  return TCL_ERROR;
	}  

	if (Tcl_GetInt(interp, argv[argi++], &secTag) != TCL_OK) {
	  opserr << "WARNING invalid secTag - element Timoshenko3d eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	  return TCL_ERROR;
	}

	if (argi >= argc || Tcl_GetInt(interp, argv[argi++], &transfTag) != TCL_OK) {
	  opserr << "WARNING invalid transfTag? - element Timoshenko3d eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	  return TCL_ERROR;
	}

	double massDens = 0.0;

	while (argi != argc) {
	  if (strcmp(argv[argi++],"-mass") == 0 && argi < argc) {
	  	if (Tcl_GetDouble(interp, argv[argi++], &massDens) != TCL_OK) {
	  	  opserr << "WARNING invalid massDens - element Timoshenko3d eleTag? iNode? jNode?"
			  << " nIP? secTag? transfTag? -mass massDens?\n";
	  	  return TCL_ERROR;
	  	}
	  }
	}

	SectionForceDeformation **sections = new SectionForceDeformation* [nIP];

	if (!sections) {
	  opserr << "WARNING TclTimoshenkoBeamCommand - Insufficient memory to create sections\n";
	  return TCL_ERROR;
	}

	SectionForceDeformation *theSection = theTclBuilder->getSection(secTag);

	if (theSection == 0) {
	  opserr << "WARNING TclTimoshenkoBeamCommand - no Section found with tag ";
	  opserr << secTag << endln;
	  delete [] sections;
	  return TCL_ERROR;
	}

	for (int j=0; j<nIP; j++) {
	  
	  if (theSection == 0) {
	    opserr << "WARNING TclTimoshenkoBeamCommand - no Section found with tag ";
	    opserr << secTag << endln;
	    delete [] sections;
	    return TCL_ERROR;
	  }
	  
	  sections[j] = theSection;
	}
	
	CrdTransf *theTransf3d = 0;
	Element *theElement = 0;

	BeamIntegration *beamIntegr = 0;
	beamIntegr = new LobattoBeamIntegration();

	if (ndm == 3) {
	  theTransf3d = OPS_GetCrdTransf(transfTag);
	  
	  if (theTransf3d == 0) {
	  	opserr << "WARNING transformation not found\n";
	  	opserr << "transformation: " << transfTag;
	  	opserr << argv[1] << " element: " << eleTag << endln;
	  	return TCL_ERROR;
	  }
	  
      // now create the Timoshenko3d04 and add it to the Domain
	  theElement = new Timoshenko3d(eleTag,iNode,jNode,nIP, sections,*theTransf3d,*beamIntegr,massDens);
	  
	  delete [] sections;
	}

	if (theElement == 0) {
	  opserr << "WARNING ran out of memory creating element\n";
	  opserr << "Timoshenko3d01 element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	if (theTclDomain->addElement(theElement) == false) {
	  opserr << "WARNING could not add element to the domain\n";
	  opserr << "Timoshenko3d01 element: " << eleTag << endln;
	  delete theElement;
	  return TCL_ERROR;
	}

	// if get here we have successfully created the element and added it to the domain
	return TCL_OK;
}


int
TclModelBuilder_addTimoshenkoBeamColumn(ClientData clientData, Tcl_Interp *interp,  
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
	if (ndm == 2 && ndf == 3)
		ok = 1;
	if (ndm == 3 && ndf == 6)
	    ok = 1;

	if (ok == 0) {
		opserr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
			<< " not compatible with Timoshenko Beam-Column element" << endln;
		return TCL_ERROR;
	}

	if (argc < 9) {			//8
		opserr << "WARNING insufficient arguments\n";
		printCommand(argc, argv);
		opserr << "Want: element Timoshenko eleTag? iNode? jNode? nIP? secTag? transfTag? <-mass> rho?\n";
		return TCL_ERROR;
	}

	// get the id and end nodes 
	int eleTag, iNode, jNode, nIP, transfTag;
	double C1;
	int secTag[10]; // Max size of integration rule ... can change if needed
	int argi = 2;

	if (Tcl_GetInt(interp, argv[argi++], &eleTag) != TCL_OK) {
		opserr << "WARNING invalid Timoshenko eleTag" << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &iNode) != TCL_OK) {
		opserr << "WARNING invalid iNode ";
		opserr << "Timoshenko element: " << eleTag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &jNode) != TCL_OK) {
		opserr << "WARNING invalid jNode ";
		opserr << "Timoshenko element: " << eleTag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &nIP) != TCL_OK) {
		opserr << "WARNING invalid nIP ";
		opserr << "Timoshenko element: " << eleTag << endln;
		return TCL_ERROR;
	}  

	if (strcmp(argv[argi], "-sections") == 0) {
		argi++;
		if (argi+nIP > argc) {
			opserr << "WARNING insufficient number of section tags - element Timoshenko eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
			return TCL_ERROR;
		}
		int section;
		for (int i = 0; i < nIP; i++) {
			if (Tcl_GetInt(interp, argv[argi+i], &section) != TCL_OK) {
				opserr << "WARNING invalid secTag - element Timoshenko eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
				return TCL_ERROR;
			}
			secTag[i] = section;
		}
		argi += nIP;
	}

	else {
		int section;
		if (Tcl_GetInt(interp, argv[argi++], &section) != TCL_OK) {
			opserr << "WARNING invalid secTag - element Timoshenko eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
			return TCL_ERROR;
		}
		for (int i = 0; i < nIP; i++)
			secTag[i] = section;
	}

	if (argi >= argc || Tcl_GetInt(interp, argv[argi++], &transfTag) != TCL_OK) {
		opserr << "WARNING invalid transfTag? - element Timoshenko eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[argi++], &C1) != TCL_OK) {
		opserr << "WARNING invalid dispBeamColumn C1" << endln;
		return TCL_ERROR;
	}

	double massDens = 0.0;

	while (argi != argc) {
		if (strcmp(argv[argi++],"-mass") == 0 && argi < argc) {
			if (Tcl_GetDouble(interp, argv[argi++], &massDens) != TCL_OK) {
				opserr << "WARNING invalid massDens - element Timoshenko eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
				return TCL_ERROR;
			}
		}
	}

	SectionForceDeformation **sections = new SectionForceDeformation* [nIP];

	if (!sections) {
		opserr << "WARNING TclElmtBuilder - addFrameElement - Insufficient memory to create sections\n";
		return TCL_ERROR;
	}

	for (int j=0; j<nIP; j++) {
		SectionForceDeformation *theSection = theTclBuilder->getSection(secTag[j]);

		if (theSection == 0) {
			opserr << "WARNING TclElmtBuilder - frameElement - no Section found with tag ";
			opserr << secTag[j] << endln;
			delete [] sections;
			return TCL_ERROR;
		}

		sections[j] = theSection;
	}

	CrdTransf *theTransf2d = 0;
	CrdTransf *theTransf3d = 0;
	Element *theElement = 0;

	BeamIntegration *beamIntegr = 0;
	beamIntegr = new LegendreBeamIntegration();

	if (ndm == 2) {
		
		theTransf2d = OPS_GetCrdTransf(transfTag);

		if (theTransf2d == 0) {
			opserr << "WARNING transformation not found\n";
			opserr << "transformation: " << transfTag;
			opserr << argv[1] << " element: " << eleTag << endln;
			return TCL_ERROR;
		}

        // now create the Timoshenko and add it to the Domain
		if ( nIP == 1 ) {
		  theElement = new TimoshenkoBeam2d(eleTag,iNode,jNode,nIP,sections,*theTransf2d,*beamIntegr,C1,massDens);
		} else {
		  opserr << "Timoshenko element processes only 1 section model for now!" << endln;
		}
		delete [] sections;
	}

	if (theElement == 0) {
		opserr << "WARNING ran out of memory creating element\n";
		opserr << "Timoshenko element: " << eleTag << endln;
		return TCL_ERROR;
	}

	if (theTclDomain->addElement(theElement) == false) {
		opserr << "WARNING could not add element to the domain\n";
		opserr << "Timoshenko element: " << eleTag << endln;
		delete theElement;
		return TCL_ERROR;
	}

	// if get here we have successfully created the element and added it to the domain
	return TCL_OK;
}




