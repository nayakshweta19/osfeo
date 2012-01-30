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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/repres/patch/tubecornerPatch.cpp,v $
                                                                        
                                                                        
#include <Matrix.h>
#include <Patch.h>
#include <tubecornerPatch.h>
#include <tubecornerCell.h>
#include <iostream>


tubecornerPatch::tubecornerPatch(): matID(0)
{

}


tubecornerPatch::tubecornerPatch(int materialID, double D, double B, double T, int SFL_D):
                       matID(materialID), d(D), b(B), t(T), sfl_d(SFL_D)
{

}  


tubecornerPatch::~tubecornerPatch()
{

}

void tubecornerPatch::setMaterialID(int materialID)
{
   matID = materialID;
}


void tubecornerPatch::setDiscretization(double D, double B, double T, int SFL_D)
{
   d = D;
   b = B;
   t = T;
   sfl_d = SFL_D;
}


int tubecornerPatch::getMaterialID(void) const
{
   return matID;
}
 
void tubecornerPatch::getDiscretization(double &D, double &B, double &T, int &SFL_D) const
{
   D = d;
   B = b;
   T = t;
   SFL_D = sfl_d;
}


int tubecornerPatch::getNumCells (void) const
{
   return 4 * sfl_d;
}


Cell ** tubecornerPatch::getCells (void) const
{ 
   Vector cellcentroid(2);
   double cellarea;
   int i, j, k;
   int numCells;
   Cell **cells;
   double f_layer_t	= 0.0;	/* LAYER THICKNESS IN FLANGE*/
   double inner_r  	= 0.0;	/* FILLET RADIUS TO INSIDE	*/
   double outer_r  	= 0.0;	/* FILLET RADIUS TO OUTSIDE	*/
   double centroid 	= 0.0;	/* CENTROID OF FIBER		*/
   double area_temp	= 0.0;	/* FIBER AREA STORAGE		*/
   double y_temp	= 0.0;	/* Y-COORDINATE STORAGE		*/
   double z_temp	= 0.0;	/* Z-COORDINATE STORAGE		*/
   int index;

   int error = 0;

   if (d > 0  && b > 0 && t > 0 && sfl_d > 0)
   {
      numCells  = this->getNumCells();

      //cout<<"tubePatch  "<<numCells<<endl;

      cells = new Cell*  [numCells];
      
      if (!cells)
         return 0;

      k = 0;

      index 		= 0;	/* FIBER NUMBER INDEX				*/
      f_layer_t	= t / sfl_d;
      /* THICKNESS OF FLANGE FIBERS			*/

      /***************************/
      /*			 */
      /*	CORNERS		 */
      /*			 */
      /***************************/

      for ( j = 1; j <= sfl_d; j++ )
      {

		/* CORNERS HAVE THE SAME NUMBER OF LAYERS IN THE THICKNESS DIR	*/
		/* AS THE FLANGES.  FOR EACH LAYER, THE INNER AND OUTER RADIUS	*/
		/* OF THE CORNER IS COMPUTED, AS WELL AS THE CENTROID POSITION  */
		/* AND THE AREA OF THE FIBERS.					*/

		inner_r 	= 2 * t - f_layer_t * j;
		outer_r 	= inner_r + f_layer_t;
		area_temp	= 3.14159 * ( pow (outer_r,2) - pow (inner_r,2) ) / 4;
		centroid 	= ( pow (outer_r,3) - pow (inner_r,3) ) / ( 3 * area_temp );

		y_temp = ( d / 2 - t * 2) + centroid;		
		z_temp = ( b / 2 - t * 2) + centroid;		

		/* FOR EACH OF THE 4 CORNERS, FILL UP THE FIBER DATA STRUCTURE	*/
		for ( k = 1; k <= 4; k++ )
		{
			index = index + 1;
			cellarea = area_temp;
		
			if ( k == 1 ){
				cellcentroid(0) = y_temp;
				cellcentroid(1) = - z_temp;
				cells[index - 1] = new tubecornerCell(matID, cellcentroid, cellarea);
			}

			if ( k == 2 ){
				cellcentroid(0) = y_temp;
				cellcentroid(1) = z_temp;
				cells[index - 1] = new tubecornerCell(matID, cellcentroid, cellarea);
			}

			if ( k == 3 ){
				cellcentroid(0) = - y_temp;
				cellcentroid(1) = z_temp;
				cells[index - 1] = new tubecornerCell(matID, cellcentroid, cellarea);
			}

			if ( k == 4 ){
				cellcentroid(0) = - y_temp;
				cellcentroid(1) = - z_temp;
				cells[index - 1] = new tubecornerCell(matID, cellcentroid, cellarea);
			}

		}  //for ( k = 1; k <= 4; k++ ) 
     }//for ( j = 1; j <= (cftxs[ ctr ]).sfl_d; j++ )
   }
   else
      return 0;

   return cells;
}


Patch * 
tubecornerPatch::getCopy (void) const
{
   tubecornerPatch *theCopy = new tubecornerPatch (matID, d, b, t, sfl_d);
   return theCopy;
}
 
void tubecornerPatch::Print(OPS_Stream &s, int flag) const
{
   s << "\nPatch Type: tubecornerPatch";
   s << "\nMaterial Id: " << matID;
   s << "\nDepth: "<< d;
   s << "\nWidth: "<< b;
   s << "\nThickness: "<< t;
   s << "\nNumber of layers in flange in depth direction: " << sfl_d;
}


OPS_Stream &operator<<(OPS_Stream &s, tubecornerPatch &tubecorner_Patch)
{
   tubecorner_Patch.Print(s);
   return s;
}

