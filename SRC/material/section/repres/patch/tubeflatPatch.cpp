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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/repres/patch/tubeflatPatch.cpp,v $
                                                                        
                                                                        
#include <Matrix.h>
#include <Patch.h>
#include <tubeflatPatch.h>
#include <tubeflatCell.h>
#include <iostream>


tubeflatPatch::tubeflatPatch(): matID(0)
{

}


tubeflatPatch::tubeflatPatch(int materialID, double D, double B, double T, int SFL_D, int SFL_B, int SWL_D, int SWL_B):
                       matID(materialID), d(D), b(B), t(T), sfl_d(SFL_D), sfl_b(SFL_B), swl_d(SWL_D), swl_b(SWL_B)
{

}  


tubeflatPatch::~tubeflatPatch()
{

}

void tubeflatPatch::setMaterialID(int materialID)
{
   matID = materialID;
}


void tubeflatPatch::setDiscretization(double D, double B, double T, int SFL_D, int SFL_B, int SWL_D, int SWL_B)
{
   d = D;
   b = B;
   t = T;
   sfl_d = SFL_D;
   sfl_b = SFL_B;
   swl_d = SWL_D;
   swl_b = SWL_B;
}


int tubeflatPatch::getMaterialID(void) const
{
   return matID;
}
 
void tubeflatPatch::getDiscretization(double &D, double &B, double &T, int &SFL_D, int &SFL_B, int &SWL_D, int &SWL_B) const
{
   D = d;
   B = b;
   T = t;
   SFL_D = sfl_d;
   SFL_B = sfl_b;
   SWL_D = swl_d;
   SWL_B = swl_b;
}


int tubeflatPatch::getNumCells (void) const
{
   return 2 * (sfl_d * sfl_b + swl_d * swl_b);
}


Cell ** tubeflatPatch::getCells (void) const
{ 
   Vector cellcentroid(2);
   double cellarea;
   int i, j, k;
   int numCells;
   Cell **cells;
   double f_layer_t	= 0.0;	/* LAYER THICKNESS IN FLANGE*/
   double w_layer_t	= 0.0;	/* LAYER THICKNESS IN WEB	*/
   double inner_r  	= 0.0;	/* FILLET RADIUS TO INSIDE	*/
   double outer_r  	= 0.0;	/* FILLET RADIUS TO OUTSIDE	*/
   double centroid 	= 0.0;	/* CENTROID OF FIBER		*/
   double area_temp	= 0.0;	/* FIBER AREA STORAGE		*/
   double y_temp	= 0.0;	/* Y-COORDINATE STORAGE		*/
   double z_temp	= 0.0;	/* Z-COORDINATE STORAGE		*/
   int index;
   double left, inc, top;

   int error = 0;

   if (d > 0  && b > 0 && t > 0 && sfl_d > 0 && sfl_b > 0 && swl_d > 0 && swl_b > 0)
   {
      numCells  = this->getNumCells();

      cells = new Cell*  [numCells];
      
      if (!cells)
         return 0;

      k = 0;

      index 		= 0;	/* FIBER NUMBER INDEX				*/
      f_layer_t	= t / sfl_d;
      /* THICKNESS OF FLANGE FIBERS			*/
      w_layer_t	= t / swl_b;
      /* THICKNESS OF WEB FIBERS			*/

     /************************/
     /*			     */
     /*	FLANGES		     */
     /*			     */
     /************************/

     left = ( 4 * t - b )/ 2.0;
     inc = ( b - 4 * t )/ sfl_b;
     area_temp = ( b - 4 * t ) * ( t ) / ( sfl_b * sfl_d );

     /* FOR EACH THROUGH-THICKNESS LAYER...				*/
     for ( j = 1; j <= sfl_d; j++ )
     {
	 y_temp = d / 2 - ( t / ( 2 * sfl_d )) * ( ( 2 * j ) - 1 );

	 /* FOR EACH FIBER IN THAT LAYER...				*/
	 for ( k = 1; k <= sfl_b; k++ )
	 {
		/********************************************************/
		/* VALUES FOR THE TOP FLANGE				*/
		/********************************************************/
		index = ( 2 * sfl_b * ( j - 1 ) + k );

		cellarea = area_temp;
		cellcentroid(0) = y_temp;
		cellcentroid(1) = left + ( inc / 2 ) * ( 2 * k - 1 );
		cells[index - 1] = new tubeflatCell(matID, cellcentroid, cellarea);

		/********************************************************/
		/* VALUES FOR THE BOT FLANGE				*/
		/********************************************************/
		index = ( ( 2 * sfl_b ) * j  - k + 1 );
			
		cellarea = area_temp;
		cellcentroid(0) = - y_temp;
		cellcentroid(1) = left + ( inc / 2 ) * ( 2 * k - 1 );
		cells[index - 1] = new tubeflatCell(matID, cellcentroid, cellarea);
	 }
      }

      /************************/
      /*		      */
      /*	WEBS   	      */
      /*		      */
      /************************/

      top = ( d -  4 * t )/ 2.0;
      inc = ( d - 4 * t )/ swl_d;
      area_temp = ( ( d - 4 * t ) * t ) / ( swl_b * swl_d );

      /* FOR EACH THROUGH-THICKNESS LAYER...					*/
      for ( j = 1; j <= swl_b; j++ )
      {
		z_temp = b / 2 - ( t / ( 2 * swl_b )) * ( ( 2 * j ) - 1 );

		y_temp = top - inc / 2.0;

		/* FOR EACH FIBER IN THAT LAYER...				*/
		for ( k = 1; k <= swl_d; k++ )
		{
			index = 2 * sfl_b * sfl_d + ( 2 * swl_d * ( j - 1 ) + k );

			cellarea = area_temp;
			cellcentroid(0) = y_temp;
			cellcentroid(1) = z_temp;
			cells[index - 1] = new tubeflatCell(matID, cellcentroid, cellarea);
			
                        index = 2 * sfl_b * sfl_d + ( 2 * swl_d * j  - k + 1 );
		
			cellarea = area_temp;
			cellcentroid(0) = y_temp;
			cellcentroid(1) = - z_temp;
			cells[index - 1] = new tubeflatCell(matID, cellcentroid, cellarea);	
			y_temp -= inc;
		}
      }
   }
   else
      return 0;

   
   return cells;
}


Patch * 
tubeflatPatch::getCopy (void) const
{
   tubeflatPatch *theCopy = new tubeflatPatch (matID, d, b, t, sfl_d, sfl_b, swl_d, swl_b);
   return theCopy;
}
 
void tubeflatPatch::Print(OPS_Stream &s, int flag) const
{
   s << "\nPatch Type: tubeflatPatch";
   s << "\nMaterial Id: " << matID;
   s << "\nDepth: "<< d;
   s << "\nWidth: "<< b;
   s << "\nThickness: "<< t;
   s << "\nNumber of layers in flange in depth direction: " << sfl_d;
   s << "\nNumber of layers in flange in width direction: " << sfl_b;
   s << "\nNumber of layers in web in depth direction: " << swl_d;
   s << "\nNumber of layers in web in width direction: " << swl_b;
}


OPS_Stream &operator<<(OPS_Stream &s, tubeflatPatch &tubeflat_Patch)
{
   tubeflat_Patch.Print(s);
   return s;
}

