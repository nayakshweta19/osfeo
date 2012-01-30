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
**							              **
**   Univeristy of Minnesota - Twin Cities		              **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.2 $
// $Date: 2008-07-03 18:08:38 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/repres/patch/corePatch.cpp,v $
                                                                        
                                                                        
// File: RCFTPatch.C
// Written by Cenk Tort
// September 2004

#include <Matrix.h>
#include <Patch.h>
#include <corePatch.h>
#include <coreCell.h>
#include <iostream>


corePatch::corePatch():matID(0)
{

}


corePatch::corePatch(int materialID, double D, double B, double T, int CFL_D, int CFL_B, int CWL_D, int CWL_B, int CC_D, int CC_B):
                       matID(materialID), d(D), b(B), t(T), cfl_d(CFL_D), cfl_b(CFL_B), cwl_d(CWL_D), cwl_b(CWL_B), cc_d(CC_D), cc_b(CC_B)

{

}  


corePatch::~corePatch()
{

}

void corePatch::setMaterialID(int materialID)
{
   matID = materialID;
}


void corePatch::setDiscretization(double D, double B, double T, int CFL_D, int CFL_B, int CWL_D, int CWL_B, int CC_D,  int CC_B)
{ 
   d = D;
   b = B;
   t = T;
   cfl_d = CFL_D;
   cfl_b = CFL_B; 
   cwl_d = CWL_D;
   cwl_b = CWL_B;
   cc_d  = CC_D;
   cc_b  = CC_B;
}


int corePatch::getMaterialID(void) const
{
   return matID;
}
 
void corePatch::getDiscretization(double &D, double &B, double &T, int &CFL_D, int &CFL_B, int &CWL_D, int &CWL_B, int &CC_D, int &CC_B) const
{
   D = d;
   B = b;
   T = t;
   CFL_D = cfl_d;
   CFL_B = cfl_b;
   CWL_D = cwl_d;
   CWL_B = cwl_b;
   CC_D  = cc_d;
   CC_B  = cc_b;
}


int corePatch::getNumCells (void) const
{
   int num_conc_fib;

   num_conc_fib = 4 * cfl_d + 2 * ( cfl_d * cfl_b + cwl_d * cwl_b ) + ( cc_d * cc_b );

   return num_conc_fib;
}


Cell **
corePatch::getCells (void) const
{
   double deltaXi;
   double deltaEta; 
   int i, j, k, r, s;
   int numCells;
   int index;
   Cell **cells;
   int error = 0;

   double f_layer_t, w_layer_t, inner_r, outer_r, centroid, y_temp, z_temp, area_temp;
   double cellarea, left, inc, top, inc_y, inc_z;

   Vector cellcentroid(2);
   
  
   if (d > 0  && b > 0 && t > 0 && cfl_d > 0 && cfl_b > 0 && cwl_b > 0 && cc_d > 0 && cc_b > 0)
   {
      numCells  = this->getNumCells();
      //cout<<"corePatch  "<<" numCells  "<<numCells<<" d  "<<d<<" b  "<<b<<" t  "<<t<<" cfl_d  "<<cfl_d<<" cfl_b  "<<cfl_b<<endl;
      //cout<<" cwl_b "<<cwl_b<<" cc_d "<<cc_d<<" cc_b "<<cc_b<<endl;

      cells = new Cell*  [numCells];
      
      if (!cells)
         return 0;

      /************************************************************************/
      /*								      */
      /*	FILL UP CONCRETE FIBER ARRAYS WITH THEIR PARAMETERS 	      */
      /*							  	      */
      /*	ASSUMPTIONS:				 		      */
      /*		1) Corners are the same # of layers as flanges	      */
      /*		2) Outside corner radius is equal to the thickness    */
      /*		3) Corners numbered first, then flanges, then webs    */
      /*			, then core elements			      */
      /*		4) Outside layers have lowest numbers, then move in   */
      /*		5) Numbering starts in upper left corner, proceeds CW */
      /*								      */
      /************************************************************************/

      index = 0;
      f_layer_t=  t / cfl_d;
      w_layer_t=  t / cwl_b;

      /************************/
      /*		      */
      /*	CORNERS	      */
      /*		      */
      /************************/

      for ( j = 1; j <= cfl_d; j++ )
      {
	inner_r = t - f_layer_t * j;
	outer_r = inner_r + f_layer_t;
	area_temp = PI * ( pow (outer_r,2) - pow (inner_r,2) ) / 4;
	centroid = ( pow (outer_r,3) - pow (inner_r,3) ) / ( 3 * area_temp );
	y_temp = (d / 2 - t * 2) + centroid;		
	z_temp = (b / 2 - t * 2) + centroid;		
	//cout<<"inner_r "<<inner_r<<" outer_r "<<outer_r<<" area_temp "<<area_temp<<" centroid "<<centroid<<endl;
	//cout<<"y_temp "<<y_temp<<" z_temp "<<z_temp<<" PI "<<PI<<endl;

	for ( k = 1; k <= 4; k++ )
	{
		index = index + 1;
		cellarea = area_temp;
		//cout<<"cellarea  "<<cellarea<<" index "<<index<<endl;
		
		if ( k == 1 )
		{
			//cout<<"index-1 "<<index - 1<<"   "<<matID<<endl;
			cellcentroid(0) = y_temp;
			cellcentroid(1) = - z_temp;
			//cout<<"cellcentroid0 "<<cellcentroid(0)<<" cellcentroid1 "<<cellcentroid(1)<<endl;
			cells[index - 1] = new coreCell(matID, cellcentroid, cellarea);
			//cout<<"index-1"<<index - 1<<endl;
		}

		if ( k == 2 )
		{
			cellcentroid(0) = y_temp;
			cellcentroid(1) = z_temp;
			cells[index - 1] = new coreCell(matID, cellcentroid, cellarea);
			//cout<<"index-1"<<index  - 1<<endl;
		}

		if ( k == 3 )
		{
			cellcentroid(0) = - y_temp;
			cellcentroid(1) = z_temp;
			cells[index - 1] = new coreCell(matID, cellcentroid, cellarea);
			//cout<<"index-1"<<index - 1<<endl;
		}

		if ( k == 4 )
		{
			cellcentroid(0) = - y_temp;
			cellcentroid(1) = - z_temp;
			cells[index - 1] = new coreCell(matID, cellcentroid, cellarea);
			//cout<<"index-1"<<index - 1<<endl;
		}
		//cout<<"core_patch "<<index-1<<"  "<<cellarea<<endl;
	}
      }
      /****************************************************************/
      /*							      */
      /*	CONCRETE LAYERS WITHIN CORNER RADIUS ALONG FLANGES    */
      /*							      */
      /****************************************************************/
      
      left = ( 4 * t - b )/ 2.0;
      inc = ( b - 4 * t )/ cfl_b;
      area_temp = ( ( b - 4 * t ) * t ) / ( cfl_b * cfl_d );

      for ( j = 1; j <= cfl_d; j++ )
      {
	y_temp = d / 2 - ( t / ( 2 * cfl_d ) ) * ( ( 2 * j ) - 1 ) - t;

	for ( k = 1; k <= cfl_b; k++ )
	{
		index = 4 * cfl_d + ( 2 * cfl_b * ( j - 1 ) + k );

		
		cellcentroid(0) = y_temp;
		cellcentroid(1) = left + ( inc / 2 ) * ( 2 * k - 1 );
		cellarea = area_temp;
		cells[index - 1] = new coreCell(matID, cellcentroid, cellarea);

		index = 4 * cfl_d + ( ( 2 * cfl_b ) * j  - k + 1 );

		cellcentroid(0) = - y_temp;
		cellcentroid(1) = left + ( inc / 2 ) * ( 2 * k - 1 );
		cellarea = area_temp;
		cells[index - 1] = new coreCell(matID, cellcentroid, cellarea);
	}
      }
      /****************************************************************/
      /*							      */
      /*	CONCRETE LAYERS WITHIN CORNER RADIUS ALONG WEBS       */
      /*							      */
      /****************************************************************/

      top = ( d -  4 * t )/ 2.0;
      inc = ( d - 4 * t )/ cwl_d;
      area_temp = ( ( d - 4 * t ) * t ) / ( cwl_b * cwl_d );

      for ( j = 1; j <= cwl_b; j++ )
      {
		z_temp = b / 2 - ( t / ( 2 * cwl_b ) ) * ( ( 2 * j ) - 1 ) - t;

		for ( k = 1; k <= cwl_d; k++ )
		{
			index = 4 * cfl_d + 2 * cfl_b * cfl_d + ( 2 * cwl_d * ( j - 1 ) + k );

			cellcentroid(0) = top - ( inc / 2 ) * ( 2 * k - 1 );
			cellcentroid(1) = z_temp;
			cellarea = area_temp;
			cells[index - 1] = new coreCell(matID, cellcentroid, cellarea);

			index = 4 * cfl_d + 2 * cfl_b * cfl_d + ( ( 2 * cwl_d ) * j  - k + 1 );

			cellcentroid(0) = top - ( inc / 2 ) * ( 2 * k - 1 );
			cellcentroid(1) = - z_temp;
			cellarea = area_temp;
			cells[index - 1] = new coreCell(matID, cellcentroid, cellarea);
		}
      }
      /************************/
      /*		      */
      /*	CORE   	      */
      /*		      */
      /************************/
      index = numCells - cc_b * cc_d;
      area_temp = (b - 4 * t) *( d - 4 * t ) / ( cc_b * cc_d );
      top = ( d -  4 * t ) / 2.0;
      left = ( 4 * t - b ) / 2.0;
      inc_y = ( d - 4 * t ) / cc_d;
      inc_z = ( b - 4 * t ) / cc_b;

      for ( i = 1; i <= cc_d; i++ )
      {
		y_temp = top - ( ( 2 * i ) - 1 ) * ( inc_y / 2 );
		z_temp = left - ( inc_z / 2 );

		for ( j = 1; j <= cc_b; j++ )
		{
			index = index + 1;
			z_temp = z_temp + inc_z;

			cellcentroid(0) = y_temp;
			cellcentroid(1) = z_temp;
			cellarea = area_temp;
			cells[index - 1] = new coreCell(matID, cellcentroid, cellarea);
		}
      }

   }
   else
      return 0;

   return cells;
}


Patch * 
corePatch::getCopy (void) const
{
   corePatch *theCopy = new corePatch (matID, d, b, t, cfl_d, cfl_b, cwl_d, cwl_b, cc_d, cc_b);
   return theCopy;
}
 
void corePatch::Print(OPS_Stream &s, int flag) const
{
   s << "\nPatch Type: corePatch";
   s << "\nMaterial Id: " << matID;
   s << "\ndepth: "<<d;
   s << "\ndepth: "<<b;
   s << "\ndepth: "<<t;
   s << "\nNumber of layers in flange in depth direction: " << cfl_d;
   s << "\nNumber of layers in flange in width direction: " << cfl_b;
   s << "\nNumber of layers in web in depth direction: " << cwl_d;
   s << "\nNumber of layers in web in width direction: " << cwl_b;
   s << "\nNumber of layers in core in depth direction: " << cc_d;
   s << "\nNumber of layers in core in width direction: " << cc_b;
}


OPS_Stream &operator<<(OPS_Stream &s, corePatch &core_Patch)
{
   core_Patch.Print(s);
   return s;
}

