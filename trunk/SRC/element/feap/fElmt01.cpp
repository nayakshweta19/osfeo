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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000/09/15 08:23:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/feap/fElmt01.cpp,v $
                                                                        
                                                                        
// File: ~/element/fortran/fElmt01.C
// 
// Written: fmk 
// Created: 03/99
// Revision: A
//
// Description: This file contains the implementation for the fElmt01 class.
//
// What: "@(#) fElement.C, revA"

#include "fElmt01.h"
#include <ID.h>
#include <Vector.h>

fElmt01::fElmt01(int tag, int nd1, int nd2, double A, double E, double rho)
:fElement(tag, ELE_TAG_fElmt01, 1, 3, 2, 2, 2, 0, 0)
{
    (*data)(0) = A;
    (*data)(1) = E;
    (*data)(2) = rho;
    
    (*connectedNodes)(0) = nd1; 
    (*connectedNodes)(1) = nd2;   
}

fElmt01::fElmt01(int tag, int nd1, int nd2, int iow)
:fElement(tag, ELE_TAG_fElmt01, 1, 3, 2, 2, 2, iow)
{
    (*connectedNodes)(0) = nd1; 
    (*connectedNodes)(1) = nd2;   
}
    
fElmt01::fElmt01()
:fElement(ELE_TAG_fElmt01)    
{
    // does nothing
}

fElmt01::~fElmt01()
{
    // does nothing
}

