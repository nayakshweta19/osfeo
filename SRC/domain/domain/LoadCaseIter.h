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
// $Date: 2000/09/15 08:23:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/LoadCaseIter.h,v $
                                                                        
                                                                        
// File: ~/domain/domain/LoadCaseIter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for LoadCaseIter.
// LoadCaseIter is an abstract base class. A LoadCAseIter object is
// an iter for returning the elements of an object of class Domain.
// LoadCaseIters must be written for each subclass of Domain, hence the
// abstract base class.

#ifndef LoadCaseIter_h
#define LoadCaseIter_h

class LoadCase;

class LoadCaseIter
{
  public:
    LoadCaseIter() {};
    virtual ~LoadCaseIter() {};
    
    virtual LoadCase *operator()(void) =0;
    
  protected:

  private:

};

#endif





