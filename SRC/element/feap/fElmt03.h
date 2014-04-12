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
// $Source: /usr/local/cvs/OpenSees/SRC/element/feap/fElmt03.h,v $
                                                                        
                                                                        
#ifndef fElmt03_h
#define fElmt03_h

// File: ~/element/fortran/fElmt03.h
// 
// Written: fmk 
// Created: 03/99
// Revision: A
//
// Description: 
//
// What: "@(#) fElmt03.h, revA"

#include <fElement.h>

class fElmt03 : public fElement
{
  public:
    // constructors
    fElmt03(int tag,
	    int Nd1, int Nd2, int Nd3, int Nd4, 
        int secTag, double rho = 0.0);
    
    fElmt03(int tag,
	    int Nd1, int Nd2, int Nd3, int Nd4, 
        int secTag, double rho = 0.0, int iow = 0);
    
    fElmt03();
    
    // destructor
    ~fElmt03();

  protected:
	     
  private:

};

#endif



