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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/repres/patch/tubePatch.cpp,v $
                                                                        
                                                                        
// File: RCFT_stl_Patch.C
// Written by Cenk Tort
// September 2004

#include <Matrix.h>
#include <Patch.h>
#include <tubePatch.h>
#include <tubeCell.h>
#include <iostream>


tubePatch::tubePatch(): matID(0)
{

}


tubePatch::tubePatch(int materialID, double D, double B, double T, int SFL_D, int SFL_B, int SWL_D, int SWL_B):
                       matID(materialID), d(D), b(B), t(T), sfl_d(SFL_D), sfl_b(SFL_B), swl_d(SWL_D), swl_b(SWL_B)
{

}  


tubePatch::~tubePatch()
{

}

void tubePatch::setMaterialID(int materialID)
{
   matID = materialID;
}


void tubePatch::setDiscretization(double D, double B, double T, int SFL_D, int SFL_B, int SWL_D, int SWL_B)
{
   d = D;
   b = B;
   t = T;
   sfl_d = SFL_D;
   sfl_b = SFL_B;
   swl_d = SWL_D;
   swl_b = SWL_B;
}


int tubePatch::getMaterialID(void) const
{
   return matID;
}
 
void tubePatch::getDiscretization(double &D, double &B, double &T, int &SFL_D, int &SFL_B, int &SWL_D, int &SWL_B) const
{
   D = d;
   B = b;
   T = t;
   SFL_D = sfl_d;
   SFL_B = sfl_b;
   SWL_D = swl_d;
   SWL_B = swl_b;
}


int tubePatch::getNumCells (void) const
{
   return 4 * sfl_d + 2 * (sfl_d * sfl_b + swl_d * swl_b);
}


Cell ** tubePatch::getCells (void) const
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

      //cout<<"tubePatch  "<<numCells<<endl;

      cells = new Cell*  [numCells];
      
      if (!cells)
         return 0;

      k = 0;

      index 		= 0;	/* FIBER NUMBER INDEX				*/
      f_layer_t	= t / sfl_d;
      /* THICKNESS OF FLANGE FIBERS			*/
      w_layer_t	= t / swl_b;
      /* THICKNESS OF WEB FIBERS			*/

      //cout<<"d "<<d<<" b "<<b<<" t "<<t<<" sfl_d "<<sfl_d<<" sfl_b "<<sfl_b<<" swl_d "<<swl_d<<" swl_b "<<swl_b<<endl;

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
			//cout<<"cellarea "<<cellarea<<endl;
		
			if ( k == 1 ){
				cellcentroid(0) = y_temp;
				cellcentroid(1) = - z_temp;
				cells[index - 1] = new tubeCell(matID, cellcentroid, cellarea);
			}

			if ( k == 2 ){
				cellcentroid(0) = y_temp;
				cellcentroid(1) = z_temp;
				cells[index - 1] = new tubeCell(matID, cellcentroid, cellarea);
			}

			if ( k == 3 ){
				cellcentroid(0) = - y_temp;
				cellcentroid(1) = z_temp;
				cells[index - 1] = new tubeCell(matID, cellcentroid, cellarea);
			}

			if ( k == 4 ){
				cellcentroid(0) = - y_temp;
				cellcentroid(1) = - z_temp;
				cells[index - 1] = new tubeCell(matID, cellcentroid, cellarea);
			}

			//cout<<"tube_patch_count  "<<index-1<<endl;

		}  //for ( k = 1; k <= 4; k++ ) 
     }//for ( j = 1; j <= (cftxs[ ctr ]).sfl_d; j++ )

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
		index = 4 * sfl_d + ( 2 * sfl_b * ( j - 1 ) + k );

		cellarea = area_temp;
		//cout<<"cellarea "<<cellarea<<endl;
		cellcentroid(0) = y_temp;
		cellcentroid(1) = left + ( inc / 2 ) * ( 2 * k - 1 );
		cells[index - 1] = new tubeCell(matID, cellcentroid, cellarea);
		//cout<<"tube_patch_count  "<<index-1<<endl;

		/********************************************************/
		/* VALUES FOR THE BOT FLANGE				*/
		/********************************************************/
		index = 4 * sfl_d + ( ( 2 * sfl_b ) * j  - k + 1 );
			
		cellarea = area_temp;
		//cout<<"cellarea "<<cellarea<<endl;
		cellcentroid(0) = - y_temp;
		cellcentroid(1) = left + ( inc / 2 ) * ( 2 * k - 1 );
		cells[index - 1] = new tubeCell(matID, cellcentroid, cellarea);
		//cout<<"tube_patch_count  "<<index-1<<endl;
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
			index = 4 * sfl_d + 2 * sfl_b * sfl_d + ( 2 * swl_d * ( j - 1 ) + k );

			cellarea = area_temp;
			//cout<<"cellarea "<<cellarea<<endl;
			cellcentroid(0) = y_temp;
			cellcentroid(1) = z_temp;
			cells[index - 1] = new tubeCell(matID, cellcentroid, cellarea);
                        //cout<<"tube_patch_count  "<<index-1<<endl;
			
                        index = 4 * sfl_d + 2 * sfl_b * sfl_d + ( 2 * swl_d * j  - k + 1 );
		
			cellarea = area_temp;
			//cout<<"cellarea "<<cellarea<<endl;
			cellcentroid(0) = y_temp;
			cellcentroid(1) = - z_temp;
			cells[index - 1] = new tubeCell(matID, cellcentroid, cellarea);	
                        //cout<<"tube_patch_count  "<<index-1<<endl;
			y_temp -= inc;
		}
      }
   }
   else
      return 0;

   
   //cout<<"call getArea 3"<<cells[3]->getArea()<<endl;
   
   return cells;
}


Patch * 
tubePatch::getCopy (void) const
{
   tubePatch *theCopy = new tubePatch (matID, d, b, t, sfl_d, sfl_b, swl_d, swl_b);
   return theCopy;
}
 
void tubePatch::Print(OPS_Stream &s, int flag) const
{
   s << "\nPatch Type: tubePatch";
   s << "\nMaterial Id: " << matID;
   s << "\nDepth: "<< d;
   s << "\nWidth: "<< b;
   s << "\nThickness: "<< t;
   s << "\nNumber of layers in flange in depth direction: " << sfl_d;
   s << "\nNumber of layers in flange in width direction: " << sfl_b;
   s << "\nNumber of layers in web in depth direction: " << swl_d;
   s << "\nNumber of layers in web in width direction: " << swl_b;
}


OPS_Stream &operator<<(OPS_Stream &s, tubePatch &tube_Patch)
{
   tube_Patch.Print(s);
   return s;
}

