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
// $Source: /usr/local/cvs/OpenSees/SRC/element/feap/fElmt03.cpp,v $
                                                                        
                                                                        
// File: ~/element/fortran/fElmt03.C
// 
// Written: fmk 
// Created: 03/99
// Revision: A
//
// Description: This file contains the implementation for the fElmt02 class.
//
// What: "@(#) fElement.C, revA"

#include "fElmt03.h"
#include <ID.h>
#include <Vector.h>

fElmt03::fElmt03(int tag, int nd1, int nd2, int nd3, int nd4, int secTag, double rho)
:fElement(tag, ELE_TAG_fElmt03, 2, 3, 2, 2, 2, 0, 0)
{
    (*data)(0) = rho;
    //(*data)(1) = E;
    //(*data)(2) = rho;
    
    (*connectedNodes)(0) = nd1;
    (*connectedNodes)(1) = nd2;
    (*connectedNodes)(2) = nd3;
    (*connectedNodes)(3) = nd4;
}

fElmt03::fElmt03(int tag, int nd1, int nd2, int nd3, int nd4, int secTag, double rho, int iow)
:fElement(tag, ELE_TAG_fElmt03, 2, 3, 2, 2, 2, iow)
{
    (*connectedNodes)(0) = nd1;
    (*connectedNodes)(1) = nd2;
    (*connectedNodes)(2) = nd3;
    (*connectedNodes)(3) = nd4;
}
    
fElmt03::fElmt03()
:fElement(ELE_TAG_fElmt03)    
{
    // does nothing
}

fElmt03::~fElmt03()
{
    // does nothing
}

