// $Source: /usr/local/cvs/OpenSees/SRC/element/Timoshenko/TclTimoshenkoBeamCommand.cpp,v $
// Created: 09/09
// Modified by: Li Ning 
// Description: This file contains the class implementation of TclModelBuilder_addTimoshenkoBeam(). Based on TclModelBuilder_addDispBeamColumn().

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include "TimoshenkoBeam2d.h"
#include "Timoshenko2d01.h"
//#include "TimoshenkoBeam3d.h"

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

	if (argc < 8) {			//8
		opserr << "WARNING insufficient arguments\n";
		printCommand(argc, argv);
		opserr << "Want: element Timoshenko eleTag? iNode? jNode? nIP? secTag? transfTag? t1? NStrip1? t2? NStrip2? t3? NStrip3?\n";
		return TCL_ERROR;
	}

	// get the id and end nodes 
	int eleTag, iNode, jNode, nIP, transfTag;
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

	double massDens = 0.0;

	while (argi != argc) {
		if (strcmp(argv[argi++],"-mass") == 0 && argi < argc) {
			if (Tcl_GetDouble(interp, argv[argi++], &massDens) != TCL_OK) {
				opserr << "WARNING invalid massDens - element Timoshenko eleTag? iNode? jNode? nIP? secTag? transfTag? t? NStrip?\n";
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
	beamIntegr = new LobattoBeamIntegration();

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
		  theElement = new TimoshenkoBeam2d(eleTag,iNode,jNode,nIP,sections,*theTransf2d,*beamIntegr,massDens);
		} else {
		  theElement = new Timoshenko2d01(eleTag,iNode,jNode,nIP,sections,*theTransf2d,*beamIntegr,massDens);
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

	// if get here we have sucessfully created the element and added it to the domain
	return TCL_OK;
}



