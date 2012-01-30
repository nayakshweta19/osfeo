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
**   Cenk Tort (tort0008@umn.edu)			              **
**   Jerome F. Hajjar (hajjar@struc.ce.umn.edu)                       **
**								      **
**   University of Minnesota - Twin Cities                            **	
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2007-09-21 15:28:40 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/repres/cell/iCell.cpp,v $
                                                                        
                                                                        
// File: iCell.cpp
// Written by Cenk Tort
// September 2004

#include <Matrix.h>
#include <Vector.h>

#include <iCell.h>


iCell::iCell(void):Centroid(2)                 
{

}


iCell::iCell(const int &type, const Vector &centroid, const double &area): 
                   Type(type), Centroid(centroid), Area(area)
{

}


iCell::~iCell()
{

}


double iCell::getArea (void) const
{

   return Area;
                
}

int iCell::getType(void) const
{
   return Type;
}


const Vector & 
iCell::getCentroidPosition(void)
{
   return Centroid;
}



void iCell::Print(OPS_Stream &s, int flag) const
{
   s << "\nCell Type: iCell";
}


OPS_Stream &operator<<(OPS_Stream &s, const iCell &i_Cell)
{
   i_Cell.Print(s);
   return s;
}    

