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
**   Univeristy of Minnesota - Twin Cities			      **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.2 $
// $Date: 2008-07-03 18:08:38 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/repres/patch/isPatch.cpp,v $
                                                                        
                                                                        
// File: iPatch.C
// Written by Cenk Tort
// September 2004

#include <Matrix.h>
#include <Patch.h>
#include <isPatch.h>
#include <isCell.h>
#include <iostream>


isPatch::isPatch(): matID(0)
{

}


isPatch::isPatch(int materialID, double D, double B, double TF, int N, double AREA):
                       matID(materialID), d(D), b(B), tf(TF), n(N), area(AREA)
{

}  


isPatch::~isPatch()
{

}

void isPatch::setMaterialID(int materialID)
{
   matID = materialID;
}


void isPatch::setDiscretization(double D, double B, double TF, int N, double AREA)
{
   d = D;
   b = B;
   tf = TF;
   n = N;
   area = AREA;
}


int isPatch::getMaterialID(void) const
{
   return matID;
}
 
void isPatch::getDiscretization(double &D, double &B, double &TF, int &N, double &AREA) const
{
   D = d;
   B = b;
   TF = tf;
   N = n;
   AREA = area;
}


int isPatch::getNumCells (void) const
{
   return n;	      
}


Cell ** isPatch::getCells (void) const
{ 
  Vector cellcentroid(2);
  int i, j, k;
  int numCells;
  Cell **cells;
  double centroid 	= 0.0;	/* CENTROID OF FIBER		*/
  double area_temp	= 0.0;	/* FIBER AREA STORAGE		*/
  double y_temp		= 0.0;	/* Y-COORDINATE STORAGE		*/
  double z_temp		= 0.0;	/* Z-COORDINATE STORAGE		*/
  int index;
  numCells  = this->getNumCells();
  cells = new Cell*  [numCells];
  if (!cells)
   return 0;

  k = 0;

  index 		= 0;	/* FIBER NUMBER INDEX		*/

  for ( i = 1; i <= n; i++ ){
       index = i;
       cellcentroid(0) = d/2 + tf/2.0;
       cellcentroid(1) = -b/2 + i * b / (n+1);
       cells[index-1] = new isCell(matID, cellcentroid, area);
  }
  
  return cells;
}


Patch * 
isPatch::getCopy (void) const
{
   isPatch *theCopy = new isPatch (matID, d, b, tf, n, area);
   return theCopy;
}
 
void isPatch::Print(OPS_Stream &s, int flag) const
{
   s << "\nPatch Type: iPatch";
   s << "\nMaterial Id: " << matID;
   s << "\nDepth: "<< d;
   s << "\nWidth: "<< b;
   s << "\nFlange Thickness: "<< tf;
   s << "\nNumber of reinforcement " << n;
   s << "\narea of rebars " << area;
}


OPS_Stream &operator<<(OPS_Stream &s, isPatch &is_Patch)
{
   is_Patch.Print(s);
   return s;
}

