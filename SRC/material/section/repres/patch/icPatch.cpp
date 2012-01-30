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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/repres/patch/icPatch.cpp,v $
                                                                        
                                                                        
// File: iPatch.C
// Written by Cenk Tort
// September 2004

#include <Matrix.h>
#include <Patch.h>
#include <icPatch.h>
#include <icCell.h>
#include <iostream>


icPatch::icPatch(): matID(0)
{

}


icPatch::icPatch(int materialID, double D, double B, double TF, int SFL_D, int SFL_B):
                       matID(materialID), d(D), b(B), tf(TF), sfl_d(SFL_D), sfl_b(SFL_B)
{

}  


icPatch::~icPatch()
{

}

void icPatch::setMaterialID(int materialID)
{
   matID = materialID;
}


void icPatch::setDiscretization(double D, double B, double TF, int SFL_D, int SFL_B)
{
   d = D;
   b = B;
   tf = TF;
   sfl_d = SFL_D;
   sfl_b = SFL_B;
}


int icPatch::getMaterialID(void) const
{
   return matID;
}
 
void icPatch::getDiscretization(double &D, double &B, double &TF, int &SFL_D, int &SFL_B) const
{
   D = d;
   B = b;
   TF = tf;
   SFL_D = sfl_d;
   SFL_B = sfl_b;
}


int icPatch::getNumCells (void) const
{
   return sfl_d * sfl_b;	      
}


Cell ** icPatch::getCells (void) const
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
  double left, inc, top;
  numCells  = this->getNumCells();
  cells = new Cell*  [numCells];
  if (!cells)
   return 0;

  k = 0;

  index 		= 0;	/* FIBER NUMBER INDEX				*/

  /************************/
  /*                      */
  /*      FLANGES         */
  /*                      */
  /************************/

  left = - b / 2;
  inc = b / sfl_b;
  area_temp =  b * tf / ( sfl_d * sfl_b );

  for ( j = 1; j <= sfl_d; j++ ){
      y_temp = d / 2 + tf - ( ( tf / ( 2 * sfl_d ) ) * ( ( 2 * j ) - 1 ) );
      for ( k = 1; k <= sfl_b; k++ ){
           index = ( sfl_b * ( j - 1 ) + k );
           cellcentroid(0) = y_temp;
	   cellcentroid(1) = left + ( inc / 2 ) * ( 2 * k - 1 );
	   cells[index-1] = new icCell(matID, cellcentroid, area_temp);
      }

  }
  
  return cells;
}


Patch * 
icPatch::getCopy (void) const
{
   icPatch *theCopy = new icPatch (matID, d, b, tf, sfl_d, sfl_b);
   return theCopy;
}
 
void icPatch::Print(OPS_Stream &s, int flag) const
{
   s << "\nPatch Type: iPatch";
   s << "\nMaterial Id: " << matID;
   s << "\nDepth: "<< d;
   s << "\nWidth: "<< b;
   s << "\nFlange Thickness: "<< tf;
   s << "\nNumber of layers in flange in depth direction: " << sfl_d;
   s << "\nNumber of layers in flange in width direction: " << sfl_b;
}


OPS_Stream &operator<<(OPS_Stream &s, icPatch &ic_Patch)
{
   ic_Patch.Print(s);
   return s;
}

