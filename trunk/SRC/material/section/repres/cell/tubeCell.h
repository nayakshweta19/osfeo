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
**   Cenk Tort (tort0008@umn.edu)				      **
**   Jerome F. Hajjar (hajjar@struc.ce.umn.edu)                       **
**								      **
**	 University of Minnesota - Twin Cities			      **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2007-09-21 15:28:40 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/repres/cell/tubeCell.h,v $
                                                                        
                                                                        
// File: tubeCell.h
//
// Written by Cenk Tort
// September 2004


#ifndef tubeCell_h 
#define tubeCell_h 

#include <Cell.h>
#include <Vector.h>

class Matrix;
class Vector;


class tubeCell: public Cell
{
  public:

    tubeCell();
    tubeCell(const int &type, const Vector &centroid, const double &area );
        
    ~tubeCell();
    
    // reinforcing bar inquiring functions
    
    double getArea(void) const;
    int getType(void) const;
    const  Vector &getCentroidPosition(void);

    void Print(OPS_Stream &s, int flag = 0) const;   
    friend OPS_Stream &operator<<(OPS_Stream &s, const tubeCell &tube_Cell);    
    
  protected:
    
  private:
    int Type;
    Vector Centroid;
    double Area;
};


#endif

