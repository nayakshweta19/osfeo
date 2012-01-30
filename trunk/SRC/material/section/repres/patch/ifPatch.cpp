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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/repres/patch/ifPatch.cpp,v $
                                                                        
                                                                        
// File: iPatch.C
// Written by Cenk Tort
// September 2004

#include <Matrix.h>
#include <Patch.h>
#include <ifPatch.h>
#include <ifCell.h>
#include <iostream>


ifPatch::ifPatch(): matID(0)
{

}


ifPatch::ifPatch(int materialID, double D, double B, double TF, double TW, double FILLET, int SFL_D, int SFL_B,
		int SWL_D, int SWL_B):
                       matID(materialID), d(D), b(B), tf(TF), tw(TW), fillet(FILLET), sfl_d(SFL_D), sfl_b(SFL_B),
		       swl_d(SWL_D), swl_b(SWL_B)
{

}  


ifPatch::~ifPatch()
{

}

void ifPatch::setMaterialID(int materialID)
{
   matID = materialID;
}


void ifPatch::setDiscretization(double D, double B, double TF, double TW, double FILLET, int SFL_D, int SFL_B, 
		int SWL_D, int SWL_B)
{
   d = D;
   b = B;
   tf = TF;
   tw = TW;
   fillet = FILLET;
   sfl_d = SFL_D;
   sfl_b = SFL_B;
   swl_d = SWL_D;
   swl_b = SWL_B;
   fillet = FILLET;
}


int ifPatch::getMaterialID(void) const
{
   return matID;
}
 
void ifPatch::getDiscretization(double &D, double &B, double &TF, double &TW, double &FILLET, int &SFL_D, 
		int &SFL_B, int &SWL_D, int &SWL_B) const
{
   D = d;
   B = b;
   TF = tf;
   TW = tw;
   FILLET = fillet;
   SFL_D = sfl_d;
   SFL_B = sfl_b;
   SWL_D = swl_d;
   SWL_B = swl_b;
}


int ifPatch::getNumCells (void) const
{
   return 4 + 2* sfl_d * sfl_b;	      
}


Cell ** ifPatch::getCells (void) const
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
  const double pi = 3.141592654;   
  if (!cells)
   return 0;

  k = 0;

  index 		= 0;	/* FIBER NUMBER INDEX				*/

  /************************/
  /*                      */
  /*      FILLETS         */
  /*                      */
  /************************/

  area_temp= ( 1 - ( pi / 4 ) ) * pow (fillet,2);
  centroid = fillet * ( 1 - ( 1 / (6 * (1 - (pi / 4)))));
  y_temp = d/2 - tf - centroid;
  z_temp = tw/2 + centroid;
   
  cellcentroid(0) = y_temp;
  cellcentroid(1) = - z_temp;
  cells[0] = new ifCell(matID, cellcentroid, area_temp);
      
  cellcentroid(0) = y_temp;
  cellcentroid(1) = z_temp;
  cells[1] = new ifCell(matID, cellcentroid, area_temp);

  cellcentroid(0) = - y_temp;
  cellcentroid(1) = z_temp;
  cells[2] = new ifCell(matID, cellcentroid, area_temp);

  cellcentroid(0) = - y_temp;
  cellcentroid(1) = - z_temp;
  cells[3] = new ifCell(matID, cellcentroid, area_temp);

  /************************/
  /*                      */
  /*      FLANGES         */
  /*                      */
  /************************/

  left = - b / 2;
  inc = b / sfl_b;
  area_temp =  b * tf / ( sfl_d * sfl_b );

  for ( j = 1; j <= sfl_d; j++ )
  {
      y_temp = d / 2 - ( ( tf / ( 2 * sfl_d ) ) * ( ( 2 * j ) - 1 ) );

      for ( k = 1; k <= sfl_b; k++ )
      {
           index = 4  + ( 2 * sfl_b * ( j - 1 ) + k );
           cellcentroid(0) = y_temp;
	   cellcentroid(1) = left + ( inc / 2 ) * ( 2 * k - 1 );
	   cells[index-1] = new ifCell(matID, cellcentroid, area_temp);

           index = 4 + ( ( 2 * sfl_b ) * j  - k + 1 );
           cellcentroid(0) = - y_temp;
	   cellcentroid(1) = left + ( inc / 2 ) * ( 2 * k - 1 );
	   cells[index-1] = new ifCell(matID, cellcentroid, area_temp);
      }

  }
  
  return cells;
}


Patch * 
ifPatch::getCopy (void) const
{
   ifPatch *theCopy = new ifPatch (matID, d, b, tf, tw, fillet, sfl_d, sfl_b, swl_d, swl_b);
   return theCopy;
}
 
void ifPatch::Print(OPS_Stream &s, int flag) const
{
   s << "\nPatch Type: iPatch";
   s << "\nMaterial Id: " << matID;
   s << "\nDepth: "<< d;
   s << "\nWidth: "<< b;
   s << "\nFlange Thickness: "<< tf;
   s << "\nWeb Thickness: "<<tw;
   s << "\nFillet: "<<fillet; 
   s << "\nNumber of layers in flange in depth direction: " << sfl_d;
   s << "\nNumber of layers in flange in width direction: " << sfl_b;
   s << "\nNumber of layers in web in depth direction: " << swl_d;
   s << "\nNumber of layers in web in width direction: " << swl_b;
}


OPS_Stream &operator<<(OPS_Stream &s, ifPatch &if_Patch)
{
   if_Patch.Print(s);
   return s;
}

