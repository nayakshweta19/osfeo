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

// $Revision: 1.8 $
// $Date: 2007-02-15 23:03:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/TclFedeasMaterialCommand.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the implementation of the
// TclModelBuilder_addFedeasMaterial() function. 

#include <PlasticDamageMaterial.h>

#include <tcl.h>
#include <Vector.h>
#include <string.h>

static void printCommand(int argc, TCL_Char **argv)
{
    opserr << "Input command: ";
    for (int i=0; i<argc; i++)
	opserr << argv[i] << " ";
    opserr << endln;
} 

UniaxialMaterial *
TclModelBuilder_addFedeasMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	if (argc < 3) {
		opserr << "WARNING insufficient number of arguments\n";
		printCommand(argc, argv);
		return 0;
	}

	int tag;
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial tag\n";
		printCommand(argc, argv);
	    return 0;
	}

	UniaxialMaterial *theMaterial = 0;

    if (strcmp(argv[1],"ConcretePlasticDamage") == 0 || strcmp(argv[1],"PlasticDamage") == 0) {
		if (argc < 11) {
			opserr << "WARNING invalid number of arguments\n";
			opserr << "Want: uniaxialMaterial ConcretePlasticDamage tag? $Ec $Gf $Gc $ft $fcy $fc $ktcr $relax" << endln;
			return 0;
		}    

		double Ec, Ft, Fc, ft_max, fcy, fc, ktcr, relax;

		if (Tcl_GetDouble(interp, argv[3], &Ec) != TCL_OK) {
		  opserr << "WARNING uniaxialMaterial ConcretePlasticDamage tag? $Ec $Ft $Fc $ft_max $fcy $fc $ktcr $relax - invalid Ec\n";
		  return 0;	
		}
		if (Tcl_GetDouble(interp, argv[4], &Ft) != TCL_OK) {
		  opserr << "WARNING uniaxialMaterial ConcretePlasticDamage tag? $Ec $Ft $Fc $ft_max $fcy $fc $ktcr $relax - invalid Ft\n";
		  return 0;	
		}
		if (Tcl_GetDouble(interp, argv[5], &Fc) != TCL_OK) {
		  opserr << "WARNING uniaxialMaterial ConcretePlasticDamage tag? $Ec $Ft $Fc $ft_max $fcy $fc $ktcr $relax - invalid Fc\n";
		  return 0;	
		}
		if (Tcl_GetDouble(interp, argv[6], &ft_max) != TCL_OK) {
		  opserr << "WARNING uniaxialMaterial ConcretePlasticDamage tag? $Ec $Ft $Fc $ft_max $fcy $fc $ktcr $relax - invalid ft_max\n";
		  return 0;	
		}
		if (Tcl_GetDouble(interp, argv[7], &fcy) != TCL_OK) {
		  opserr << "WARNING uniaxialMaterial ConcretePlasticDamage tag? $Ec $Ft $Fc $ft_max $fcy $fc $ktcr $relax - invalid fcy\n";
		  return 0;	
		}
		if (Tcl_GetDouble(interp, argv[8], &fc) != TCL_OK) {
		  opserr << "WARNING uniaxialMaterial ConcretePlasticDamage tag? $Ec $Ft $Fc $ft_max $fcy $fc $ktcr $relax - invalid fc\n";
		  return 0;	
		}
		if (Tcl_GetDouble(interp, argv[9], &ktcr) != TCL_OK) {
		  opserr << "WARNING uniaxialMaterial ConcretePlasticDamage tag? $Ec $Ft $Fc $ft_max $fcy $fc $ktcr $relax - invalid Ktcr\n";
		  return 0;	
		}
		if (Tcl_GetDouble(interp, argv[10], &relax) != TCL_OK) {
		  opserr << "WARNING uniaxialMaterial ConcretePlasticDamage tag? $Ec $Ft $Fc $ft_max $fcy $fc $ktcr $relax - invalid relax\n";
		  return 0;	
		}

		theMaterial = new PlasticDamageMaterial(tag, Ec, Ft, Fc, ft_max, fcy, fc, ktcr, relax);
		
	}

	return theMaterial;
}
