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
**   University of Minnesota - Twin Cities                            **
** ****************************************************************** */

// $Revision: 1.1.1.1 $
// $Date: 2007-09-21 15:28:09 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/coordTransformation/RCFTCrdLinTransf3D.cpp,v $


// File: ~/crdTransf/RCFTCrdTransf3D.cpp
//
// Purpose: This file contains the implementation for the
// RCFTCrdTransf3D class. RCFTCrdTransf3D is a linear
// transformation for a planar frame between the global
// and basic coordinate systems


#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <Channel.h>
//#include <iostream.h>
//#include <fstream.h>

#include <RCFTCrdLinTransf3D.h>

//using std::ios;
//using std::endl;
//using std::ofstream;

// constructor:
RCFTCrdLinTransf3D::RCFTCrdLinTransf3D(int tag, const Vector &vecInLocXZPlane):
  CrdTransf(tag, CRDTR_TAG_RCFTCrdLinTransf3D),
  nodeIPtr(0), nodeJPtr(0), L(0), kl(18,18), pg(18)
{
  	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			R[i][j] = 0.0;

	R[1][0] = vecInLocXZPlane(0);
	R[1][1] = vecInLocXZPlane(1);
	R[1][2] = vecInLocXZPlane(2);

}


// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
RCFTCrdLinTransf3D::RCFTCrdLinTransf3D():
  CrdTransf(0, CRDTR_TAG_RCFTCrdLinTransf3D),
  nodeIPtr(0), nodeJPtr(0), L(0), pg(18)
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			R[i][j] = 0.0;
}


// destructor:
RCFTCrdLinTransf3D::~RCFTCrdLinTransf3D()
{

}


int
RCFTCrdLinTransf3D::commitState(void)
{
   return 0;
}


int
RCFTCrdLinTransf3D::revertToLastCommit(void)
{
   return 0;
}


int
RCFTCrdLinTransf3D::revertToStart(void)
{
   return 0;
}


int
RCFTCrdLinTransf3D::initialize(Node *nodeIPointer, Node *nodeJPointer)
{
   int error;

   nodeIPtr = nodeIPointer;
   nodeJPtr = nodeJPointer;

   if ((!nodeIPtr) || (!nodeJPtr))
   {
      opserr << "\nRCFTCrdTransf3D::initialize";
      opserr << "\ninvalid pointers to the element nodes\n";
      return -1;
   }
 
   // get element length and orientation
   if ((error = this->computeElemtLengthAndOrient()))
      return error;

   static Vector XAxis(3);
   static Vector YAxis(3);
   static Vector ZAxis(3);

   const Vector &ndICoords = nodeIPtr->getCrds();
   const Vector &ndJCoords = nodeJPtr->getCrds();
   
   Vector dx(3);

   dx(0) = ndJCoords(0) - ndICoords(0);
   dx(1) = ndJCoords(1) - ndICoords(1);
   dx(2) = ndJCoords(2) - ndICoords(2);
 
   double L = sqrt( dx(0) * dx(0) + dx(1) * dx(1) + dx(2) * dx(2) );

   double cx = dx(0) / L;
   double cy = dx(1) / L;
   double cz = dx(2) / L;

   XAxis(0) = cx; XAxis(1) = cy; XAxis(2) = cz;

   YAxis(0) = R[1][0];  YAxis(1) = R[1][1];     YAxis(2) = R[1][2];

   /******************************************************/
   /* COMPUTE OUT-OF-PLANE VECTOR: {outp} = {c} x {iend} */
   /******************************************************/

   ZAxis(0)   = ( XAxis(1) * YAxis(2) ) - ( XAxis(2) * YAxis(1) );
   ZAxis(1)   = ( XAxis(2) * YAxis(0) ) - ( XAxis(0) * YAxis(2) );
   ZAxis(2)   = ( XAxis(0) * YAxis(1) ) - ( XAxis(1) * YAxis(0) );

   double znorm = ZAxis.Norm();

   if (znorm == 0) {
         opserr << "\nRCFTCrdLinTransf3D::getLocalAxes3D";
         opserr << "\nvector v that defines plane xz is parallel to x axis\n";
         return -3;
   }

   ZAxis /= znorm;

   // Fill in transformation matrix

   R[0][0]   = XAxis(0);
   R[0][1]   = XAxis(1);
   R[0][2]   = XAxis(2);

   R[2][0]   = ZAxis(0);
   R[2][1]   = ZAxis(1);
   R[2][2]   = ZAxis(2);

   return 0;
}


int
RCFTCrdLinTransf3D::update(void)
{
   return 0;
}


int
RCFTCrdLinTransf3D::computeElemtLengthAndOrient()
{
   // element projection
   static Vector dx(3);

   const Vector &ndICoords = nodeIPtr->getCrds();
   const Vector &ndJCoords = nodeJPtr->getCrds();

   dx(0) = ndJCoords(0) - ndICoords(0);
   dx(1) = ndJCoords(1) - ndICoords(1);
   dx(2) = ndJCoords(2) - ndICoords(2);

   // calculate the element length
   L = dx.Norm();

   if (L == 0.0) {
      opserr << "\nRCFTCrdLinTransf3d::computeElemtLengthAndOrien: 0 length\n";
      return -2;
   }

   // calculate the element local x axis components (direction cosines)
   // wrt to the global coordinates
   R[0][0]   = dx(0)/L;
   R[0][1]   = dx(1)/L;
   R[0][2]   = dx(2)/L;

   return 0;
}

int
RCFTCrdLinTransf3D::getLocalAxes(Vector &XAxis, Vector &YAxis, Vector &ZAxis)
{
   XAxis(0) = R[0][0]; XAxis(1) = R[0][1]; XAxis(2) = R[0][2];
   YAxis(0) = R[1][0]; YAxis(1) = R[1][1]; YAxis(2) = R[1][2];
   ZAxis(0) = R[2][0]; ZAxis(1) = R[2][1]; ZAxis(2) = R[2][2];

   return 0;
}


double
RCFTCrdLinTransf3D::getInitialLength(void)
{
   return L;
}

double
RCFTCrdLinTransf3D::getDeformedLength(void)
{
   return L;
}


const Vector &
RCFTCrdLinTransf3D::getBasicTrialDisp (void)
{
   Vector ub(12);
   ub.Zero();
   return ub;
}


const Vector &
RCFTCrdLinTransf3D::getBasicIncrDisp (void)
{
   Vector ub(12);
   ub.Zero();
   return ub;
}


const Vector &
RCFTCrdLinTransf3D::getBasicIncrDeltaDisp (void)
{
   Vector ub(12);
   ub.Zero();
   return ub;
}


const Vector &
RCFTCrdLinTransf3D::getBasicTrialVel(void)
{
   Vector ub(12);
   ub.Zero();
   return ub;
}


const Vector &
RCFTCrdLinTransf3D::getBasicTrialAccel(void)
{
   Vector ub(12);
   ub.Zero();
   return ub;
}


const Vector&
RCFTCrdLinTransf3D::getGlobalResistingForce(const Vector &sg, const Vector &p0)
{
  pg(0) = sg(0)*R[0][0] + sg(1)*R[1][0]  + sg(2)*R[2][0];
  pg(1) = sg(0)*R[0][1] + sg(1)*R[1][1]  + sg(2)*R[2][1];
  pg(2) = sg(0)*R[0][2] + sg(1)*R[1][2]  + sg(2)*R[2][2];

  pg(3) = sg(3)*R[0][0] + sg(4)*R[1][0]  + sg(5)*R[2][0];
  pg(4) = sg(3)*R[0][1] + sg(4)*R[1][1]  + sg(5)*R[2][1];
  pg(5) = sg(3)*R[0][2] + sg(4)*R[1][2]  + sg(5)*R[2][2];

  pg(6) = sg(6)*R[0][0] + sg(7)*R[1][0]  + sg(8)*R[2][0];
  pg(7) = sg(6)*R[0][1] + sg(7)*R[1][1]  + sg(8)*R[2][1];
  pg(8) = sg(6)*R[0][2] + sg(7)*R[1][2]  + sg(8)*R[2][2];

  pg(9)  = sg(9)*R[0][0] + sg(10)*R[1][0] + sg(11)*R[2][0];
  pg(10) = sg(9)*R[0][1] + sg(10)*R[1][1] + sg(11)*R[2][1];
  pg(11) = sg(9)*R[0][2] + sg(10)*R[1][2] + sg(11)*R[2][2];

  pg(12) = sg(12)*R[0][0] + sg(13)*R[1][0] + sg(14)*R[2][0];
  pg(13) = sg(12)*R[0][1] + sg(13)*R[1][1] + sg(14)*R[2][1];
  pg(14) = sg(12)*R[0][2] + sg(13)*R[1][2] + sg(14)*R[2][2];

  pg(15) = sg(15)*R[0][0] + sg(16)*R[1][0] + sg(17)*R[2][0];
  pg(16) = sg(15)*R[0][1] + sg(16)*R[1][1] + sg(17)*R[2][1];
  pg(17) = sg(15)*R[0][2] + sg(16)*R[1][2] + sg(17)*R[2][2];

  return pg;
}


const Matrix &
RCFTCrdLinTransf3D::getInitialGlobalStiffMatrix (const Matrix &KB)
{
  static Matrix kg(18,18);	// Global stiffness for return
  static double kb[12][12];
  Matrix kl(18,18);
  static double tmp[18][18];	// Temporary storage

  double oneOverL = 1.0/L;

  static double nat_to_local[18][12];   //Natural to local transformation
  static double nat_to_local_T[12][18]; //Transpose of nat_to_local
  static double temp_rt[12][18];        //Temporary storage for transformation

  int i,j,k;
  for (i = 0; i < 12; i++)
	for (j = 0; j < 12; j++)
	    kb[i][j] = KB(i,j);

  /************************************************************************/
  /* CALCULATE THE NATURAL-TO-LOCAL TRANSFORMATION MATRIX 		  */
  /************************************************************************/
  nat_to_local[0][0]     = -1.0;
  nat_to_local[0][6]     = 1.0;
  nat_to_local[1][6]     = - 1.0;
  nat_to_local[1][15]    = 1.0;
  nat_to_local[2][5]     = 1.0;
  nat_to_local[2][7]     = 1.0 / L;
  nat_to_local[2][16]    = - 1.0 / L;
  nat_to_local[3][4]     = - 1.0;
  nat_to_local[3][8]     = 1.0 / L;
  nat_to_local[3][17]    = - 1.0 / L;
  nat_to_local[4][7]     = 1.0 / L;
  nat_to_local[4][14]    = 1.0;
  nat_to_local[4][16]    = - 1.0 / L;
  nat_to_local[5][8]     = 1.0 / L;
  nat_to_local[5][13]    = - 1.0;
  nat_to_local[5][17]    = - 1.0 / L;
  nat_to_local[6][0]     = -1.0;
  nat_to_local[6][9]     = 1.0;
  nat_to_local[7][1]     = 1.0 / L;
  nat_to_local[7][5]     = 1.0;
  nat_to_local[7][10]    = - 1.0 / L;
  nat_to_local[8][2]     = 1.0 / L;
  nat_to_local[8][4]     = - 1.0;
  nat_to_local[8][11]    = - 1.0 / L;
  nat_to_local[9][1]     = 1.0 / L;
  nat_to_local[9][10]    = - 1.0 / L;
  nat_to_local[9][14]    = 1.0;
  nat_to_local[10][2]    = 1.0 / L;
  nat_to_local[10][11]   = - 1.0 / L;
  nat_to_local[10][13]   = - 1.0;
  nat_to_local[11][3]    = - 1.0;
  nat_to_local[11][12]   = 1.0;

  /************************************************************************/
  /* AND ITS TRANSPOSE [nat_to_local](T)		  		  */
  /************************************************************************/
  for ( i = 0; i < 12; i++ ){
    for ( k = 0; k < 18; k++ ){
		nat_to_local_T[i][k] = nat_to_local[k][i];
    }
  }

  /************************************************************************/
  /* MATRIX MULTIPLICATION  [k] * [nat_to_local](T)		          */
  /************************************************************************/
  for ( i = 0; i < 12; i++ ){
	for ( k = 0; k < 18; k++ ){
		for ( j = 0; j < 12; j++ ){
		temp_rt[i][k] += kb[i][j] * nat_to_local_T[j][k];
		}
	}
  }

  /*************************************************************************/
  /* MATRIX MULTIPLICATION  [nat_to_local] * [k] * [nat_to_local](T)	   */
  /*************************************************************************/
  for ( i = 0; i < 18; i++ ){
     for ( k = 0; k < 18; k++ ){
	for ( j = 0; j < 12; j++ ){
	    kl(i, k) += nat_to_local[i][j] * temp_rt[j][k];
	}
     }
  }

  // Transform local stiffness to global system
  // First compute kl*T_{lg}
  for (i = 0; i < 18; i++) {
	tmp[i][0] = kl(i,0)*R[0][0] + kl(i,1)*R[1][0]  + kl(i,2)*R[2][0];
	tmp[i][1] = kl(i,0)*R[0][1] + kl(i,1)*R[1][1]  + kl(i,2)*R[2][1];
	tmp[i][2] = kl(i,0)*R[0][2] + kl(i,1)*R[1][2]  + kl(i,2)*R[2][2];

	tmp[i][3] = kl(i,3)*R[0][0] + kl(i,4)*R[1][0]  + kl(i,5)*R[2][0];
	tmp[i][4] = kl(i,3)*R[0][1] + kl(i,4)*R[1][1]  + kl(i,5)*R[2][1];
	tmp[i][5] = kl(i,3)*R[0][2] + kl(i,4)*R[1][2]  + kl(i,5)*R[2][2];

	tmp[i][6] = kl(i,6)*R[0][0] + kl(i,7)*R[1][0]  + kl(i,8)*R[2][0];
	tmp[i][7] = kl(i,6)*R[0][1] + kl(i,7)*R[1][1]  + kl(i,8)*R[2][1];
	tmp[i][8] = kl(i,6)*R[0][2] + kl(i,7)*R[1][2]  + kl(i,8)*R[2][2];

	tmp[i][9]  = kl(i,9)*R[0][0] + kl(i,10)*R[1][0] + kl(i,11)*R[2][0];
	tmp[i][10] = kl(i,9)*R[0][1] + kl(i,10)*R[1][1] + kl(i,11)*R[2][1];
	tmp[i][11] = kl(i,9)*R[0][2] + kl(i,10)*R[1][2] + kl(i,11)*R[2][2];

	tmp[i][12] = kl(i,12)*R[0][0] + kl(i,13)*R[1][0] + kl(i,14)*R[2][0];
	tmp[i][13] = kl(i,12)*R[0][1] + kl(i,13)*R[1][1] + kl(i,14)*R[2][1];
	tmp[i][14] = kl(i,12)*R[0][2] + kl(i,13)*R[1][2] + kl(i,14)*R[2][2];

	tmp[i][15] = kl(i,15)*R[0][0] + kl(i,16)*R[1][0] + kl(i,17)*R[2][0];
	tmp[i][16] = kl(i,15)*R[0][1] + kl(i,16)*R[1][1] + kl(i,17)*R[2][1];
	tmp[i][17] = kl(i,15)*R[0][2] + kl(i,16)*R[1][2] + kl(i,17)*R[2][2];

  }

  // Now compute T'_{lg}*(kl*T_{lg})
  for (i = 0; i < 18; i++) {
	kg(0,i) = R[0][0]*tmp[0][i] + R[1][0]*tmp[1][i]  + R[2][0]*tmp[2][i];
	kg(1,i) = R[0][1]*tmp[0][i] + R[1][1]*tmp[1][i]  + R[2][1]*tmp[2][i];
	kg(2,i) = R[0][2]*tmp[0][i] + R[1][2]*tmp[1][i]  + R[2][2]*tmp[2][i];

	kg(3,i) = R[0][0]*tmp[3][i] + R[1][0]*tmp[4][i]  + R[2][0]*tmp[5][i];
	kg(4,i) = R[0][1]*tmp[3][i] + R[1][1]*tmp[4][i]  + R[2][1]*tmp[5][i];
	kg(5,i) = R[0][2]*tmp[3][i] + R[1][2]*tmp[4][i]  + R[2][2]*tmp[5][i];

	kg(6,i) = R[0][0]*tmp[6][i] + R[1][0]*tmp[7][i]  + R[2][0]*tmp[8][i];
	kg(7,i) = R[0][1]*tmp[6][i] + R[1][1]*tmp[7][i]  + R[2][1]*tmp[8][i];
	kg(8,i) = R[0][2]*tmp[6][i] + R[1][2]*tmp[7][i]  + R[2][2]*tmp[8][i];

	kg(9,i)  = R[0][0]*tmp[9][i] + R[1][0]*tmp[10][i] + R[2][0]*tmp[11][i];
	kg(10,i) = R[0][1]*tmp[9][i] + R[1][1]*tmp[10][i] + R[2][1]*tmp[11][i];
	kg(11,i) = R[0][2]*tmp[9][i] + R[1][2]*tmp[10][i] + R[2][2]*tmp[11][i];

	kg(12,i) = R[0][0]*tmp[12][i] + R[1][0]*tmp[13][i] + R[2][0]*tmp[14][i];
	kg(13,i) = R[0][1]*tmp[12][i] + R[1][1]*tmp[13][i] + R[2][1]*tmp[14][i];
	kg(14,i) = R[0][2]*tmp[12][i] + R[1][2]*tmp[13][i] + R[2][2]*tmp[14][i];

       	kg(15,i) = R[0][0]*tmp[15][i] + R[1][0]*tmp[16][i] + R[2][0]*tmp[17][i];
	kg(16,i) = R[0][1]*tmp[15][i] + R[1][1]*tmp[16][i] + R[2][1]*tmp[17][i];
	kg(17,i) = R[0][2]*tmp[15][i] + R[1][2]*tmp[16][i] + R[2][2]*tmp[17][i];

  }
  return kg;
}

const Matrix &
RCFTCrdLinTransf3D::getLocalStiffMatrix (const Matrix &KB, const Vector &fk)
{
   static Matrix kg(18,18);	// Global stiffness for return
   static double kb[12][12];
   
   Matrix kl(18,18);
   static double tmp[18][18];	// Temporary storage

   double oneOverL = 1.0/L;

   double nat_to_local[12][18];   //Natural to local transformation
   double nat_to_local_T[18][12]; //Transpose of nat_to_local
   double temp_rt[12][18];        //Temporary storage for transformation

   int i,j,k;

   for(i=0; i < 12; i++){
        for(j=0; j < 18; j++){
	    nat_to_local[i][j] = 0.0;
	    nat_to_local_T[j][i] = 0.0;
	    temp_rt[i][j] = 0.0;
	} 
   }
   
   for (i = 0; i < 12; i++)
	for (j = 0; j < 12; j++)
	    kb[i][j] = KB(i,j);

   /************************************************************************/
   /* CALCULATE THE NATURAL-TO-LOCAL TRANSFORMATION MATRIX 		   */
   /************************************************************************/
   nat_to_local[0][0]     = -1.0;
   nat_to_local[0][6]     = 1.0;
   nat_to_local[1][6]     = - 1.0;
   nat_to_local[1][15]    = 1.0;
   nat_to_local[2][5]     = 1.0;
   nat_to_local[2][7]     = 1.0 / L;
   nat_to_local[2][16]    = - 1.0 / L;
   nat_to_local[3][4]     = - 1.0;
   nat_to_local[3][8]     = 1.0 / L;
   nat_to_local[3][17]    = - 1.0 / L;
   nat_to_local[4][7]     = 1.0 / L;
   nat_to_local[4][14]    = 1.0;
   nat_to_local[4][16]    = - 1.0 / L;
   nat_to_local[5][8]     = 1.0 / L;
   nat_to_local[5][13]    = - 1.0;
   nat_to_local[5][17]    = - 1.0 / L;
   nat_to_local[6][0]     = -1.0;
   nat_to_local[6][9]     = 1.0;
   nat_to_local[7][1]     = 1.0 / L;
   nat_to_local[7][5]     = 1.0;
   nat_to_local[7][10]    = - 1.0 / L;
   nat_to_local[8][2]     = 1.0 / L;
   nat_to_local[8][4]     = - 1.0;
   nat_to_local[8][11]    = - 1.0 / L;
   nat_to_local[9][1]     = 1.0 / L;
   nat_to_local[9][10]    = - 1.0 / L;
   nat_to_local[9][14]    = 1.0;
   nat_to_local[10][2]    = 1.0 / L;
   nat_to_local[10][11]   = - 1.0 / L;
   nat_to_local[10][13]   = - 1.0;
   nat_to_local[11][3]    = - 1.0;
   nat_to_local[11][12]   = 1.0;

   /************************************************************************/
   /* AND ITS TRANSPOSE [nat_to_local](T)				   */
   /************************************************************************/
   for ( i = 0; i < 12; i++ ){
        for ( k = 0; k < 18; k++ ){
		nat_to_local_T[ k ][ i ] = nat_to_local[ i ][ k ];
	}
   }

   /************************************************************************/
   /* MATRIX MULTIPLICATION  [k] * [nat_to_local](T)		       	   */
   /************************************************************************/
   for ( i = 0; i < 12; i++ ){
	for ( k = 0; k < 18; k++ ){
		for ( j = 0; j < 12; j++ ){
		temp_rt[ i ][ k ] += kb[ i ][ j ] * nat_to_local[ j ][ k ];
		}
	}
   }

   /*************************************************************************/
   /* MATRIX MULTIPLICATION  [nat_to_local] * [k] * [nat_to_local](T)	    */
   /*************************************************************************/
   for ( i = 0; i < 18; i++ ){
	for ( k = 0; k < 18; k++ ){
		for ( j = 0; j < 12; j++ ){
		    kl(i, k) += nat_to_local_T[ i ][ j ] * temp_rt[ j ][ k ];
		}
	}
   }
   return kl;
}


const Matrix &
RCFTCrdLinTransf3D::getGlobalStiffMatrix (const Matrix &KB, const Vector &fk)
{
#ifdef COMPOSITE_DEBUG
   ofstream check4;
   check4.open("check4.dat",ios::app);

   check4>>KB;
#endif

   static Matrix kg(18,18);	// Global stiffness for return
   static double kb[12][12];
   
   Matrix kl(18,18);
   static double tmp[18][18];	// Temporary storage

   double oneOverL = 1.0/L;

   double nat_to_local[12][18];   //Natural to local transformation
   double nat_to_local_T[18][12]; //Transpose of nat_to_local
   double temp_rt[12][18];        //Temporary storage for transformation

   int i,j,k;

   for(i=0; i < 12; i++){
        for(j=0; j < 18; j++){
	    nat_to_local[i][j] = 0.0;
	    nat_to_local_T[j][i] = 0.0;
	    temp_rt[i][j] = 0.0;
	} 
   }
   
   for (i = 0; i < 12; i++)
	for (j = 0; j < 12; j++)
	    kb[i][j] = KB(i,j);

   /************************************************************************/
   /* CALCULATE THE NATURAL-TO-LOCAL TRANSFORMATION MATRIX 		   */
   /************************************************************************/
   nat_to_local[0][0]     = -1.0;
   nat_to_local[0][6]     = 1.0;
   nat_to_local[1][6]     = - 1.0;
   nat_to_local[1][15]    = 1.0;
   nat_to_local[2][5]     = 1.0;
   nat_to_local[2][7]     = 1.0 / L;
   nat_to_local[2][16]    = - 1.0 / L;
   nat_to_local[3][4]     = - 1.0;
   nat_to_local[3][8]     = 1.0 / L;
   nat_to_local[3][17]    = - 1.0 / L;
   nat_to_local[4][7]     = 1.0 / L;
   nat_to_local[4][14]    = 1.0;
   nat_to_local[4][16]    = - 1.0 / L;
   nat_to_local[5][8]     = 1.0 / L;
   nat_to_local[5][13]    = - 1.0;
   nat_to_local[5][17]    = - 1.0 / L;
   nat_to_local[6][0]     = -1.0;
   nat_to_local[6][9]     = 1.0;
   nat_to_local[7][1]     = 1.0 / L;
   nat_to_local[7][5]     = 1.0;
   nat_to_local[7][10]    = - 1.0 / L;
   nat_to_local[8][2]     = 1.0 / L;
   nat_to_local[8][4]     = - 1.0;
   nat_to_local[8][11]    = - 1.0 / L;
   nat_to_local[9][1]     = 1.0 / L;
   nat_to_local[9][10]    = - 1.0 / L;
   nat_to_local[9][14]    = 1.0;
   nat_to_local[10][2]    = 1.0 / L;
   nat_to_local[10][11]   = - 1.0 / L;
   nat_to_local[10][13]   = - 1.0;
   nat_to_local[11][3]    = - 1.0;
   nat_to_local[11][12]   = 1.0;

   /************************************************************************/
   /* AND ITS TRANSPOSE [nat_to_local](T)				   */
   /************************************************************************/
   for ( i = 0; i < 12; i++ ){
        for ( k = 0; k < 18; k++ ){
		nat_to_local_T[ k ][ i ] = nat_to_local[ i ][ k ];
	}
   }

   /************************************************************************/
   /* MATRIX MULTIPLICATION  [k] * [nat_to_local](T)		       	   */
   /************************************************************************/
   for ( i = 0; i < 12; i++ ){
	for ( k = 0; k < 18; k++ ){
		for ( j = 0; j < 12; j++ ){
		temp_rt[ i ][ k ] += kb[ i ][ j ] * nat_to_local[ j ][ k ];
		}
	}
   }

   /*************************************************************************/
   /* MATRIX MULTIPLICATION  [nat_to_local] * [k] * [nat_to_local](T)	    */
   /*************************************************************************/
   for ( i = 0; i < 18; i++ ){
	for ( k = 0; k < 18; k++ ){
		for ( j = 0; j < 12; j++ ){
		    kl(i, k) += nat_to_local_T[ i ][ j ] * temp_rt[ j ][ k ];
		}
	}
   }

   // Transform local stiffness to global system
   // First compute kl*T_{lg}
   for (i = 0; i < 18; i++) {
	tmp[i][0] = kl(i,0)*R[0][0] + kl(i,1)*R[1][0]  + kl(i,2)*R[2][0];
	tmp[i][1] = kl(i,0)*R[0][1] + kl(i,1)*R[1][1]  + kl(i,2)*R[2][1];
	tmp[i][2] = kl(i,0)*R[0][2] + kl(i,1)*R[1][2]  + kl(i,2)*R[2][2];

	tmp[i][3] = kl(i,3)*R[0][0] + kl(i,4)*R[1][0]  + kl(i,5)*R[2][0];
	tmp[i][4] = kl(i,3)*R[0][1] + kl(i,4)*R[1][1]  + kl(i,5)*R[2][1];
	tmp[i][5] = kl(i,3)*R[0][2] + kl(i,4)*R[1][2]  + kl(i,5)*R[2][2];

	tmp[i][6] = kl(i,6)*R[0][0] + kl(i,7)*R[1][0]  + kl(i,8)*R[2][0];
	tmp[i][7] = kl(i,6)*R[0][1] + kl(i,7)*R[1][1]  + kl(i,8)*R[2][1];
	tmp[i][8] = kl(i,6)*R[0][2] + kl(i,7)*R[1][2]  + kl(i,8)*R[2][2];

	tmp[i][9]  = kl(i,9)*R[0][0] + kl(i,10)*R[1][0] + kl(i,11)*R[2][0];
	tmp[i][10] = kl(i,9)*R[0][1] + kl(i,10)*R[1][1] + kl(i,11)*R[2][1];
	tmp[i][11] = kl(i,9)*R[0][2] + kl(i,10)*R[1][2] + kl(i,11)*R[2][2];

	tmp[i][12] = kl(i,12)*R[0][0] + kl(i,13)*R[1][0] + kl(i,14)*R[2][0];
	tmp[i][13] = kl(i,12)*R[0][1] + kl(i,13)*R[1][1] + kl(i,14)*R[2][1];
	tmp[i][14] = kl(i,12)*R[0][2] + kl(i,13)*R[1][2] + kl(i,14)*R[2][2];

	tmp[i][15] = kl(i,15)*R[0][0] + kl(i,16)*R[1][0] + kl(i,17)*R[2][0];
	tmp[i][16] = kl(i,15)*R[0][1] + kl(i,16)*R[1][1] + kl(i,17)*R[2][1];
	tmp[i][17] = kl(i,15)*R[0][2] + kl(i,16)*R[1][2] + kl(i,17)*R[2][2];
  }
  // Now compute T'_{lg}*(kl*T_{lg})
  for (i = 0; i < 18; i++) {
	kg(0,i) = R[0][0]*tmp[0][i] + R[1][0]*tmp[1][i] + R[2][0]*tmp[2][i];
	kg(1,i) = R[0][1]*tmp[0][i] + R[1][1]*tmp[1][i] + R[2][1]*tmp[2][i];
	kg(2,i) = R[0][2]*tmp[0][i] + R[1][2]*tmp[1][i] + R[2][2]*tmp[2][i];

	kg(3,i) = R[0][0]*tmp[3][i] + R[1][0]*tmp[4][i] + R[2][0]*tmp[5][i];
	kg(4,i) = R[0][1]*tmp[3][i] + R[1][1]*tmp[4][i] + R[2][1]*tmp[5][i];
	kg(5,i) = R[0][2]*tmp[3][i] + R[1][2]*tmp[4][i] + R[2][2]*tmp[5][i];

	kg(6,i) = R[0][0]*tmp[6][i] + R[1][0]*tmp[7][i] + R[2][0]*tmp[8][i];
	kg(7,i) = R[0][1]*tmp[6][i] + R[1][1]*tmp[7][i] + R[2][1]*tmp[8][i];
	kg(8,i) = R[0][2]*tmp[6][i] + R[1][2]*tmp[7][i] + R[2][2]*tmp[8][i];

	kg(9,i)  = R[0][0]*tmp[9][i] + R[1][0]*tmp[10][i] + R[2][0]*tmp[11][i];
	kg(10,i) = R[0][1]*tmp[9][i] + R[1][1]*tmp[10][i] + R[2][1]*tmp[11][i];
	kg(11,i) = R[0][2]*tmp[9][i] + R[1][2]*tmp[10][i] + R[2][2]*tmp[11][i];

	kg(12,i) = R[0][0]*tmp[12][i] + R[1][0]*tmp[13][i] + R[2][0]*tmp[14][i];
	kg(13,i) = R[0][1]*tmp[12][i] + R[1][1]*tmp[13][i] + R[2][1]*tmp[14][i];
	kg(14,i) = R[0][2]*tmp[12][i] + R[1][2]*tmp[13][i] + R[2][2]*tmp[14][i];

        kg(15,i) = R[0][0]*tmp[15][i] + R[1][0]*tmp[16][i] + R[2][0]*tmp[17][i];
        kg(16,i) = R[0][1]*tmp[15][i] + R[1][1]*tmp[16][i] + R[2][1]*tmp[17][i];
	kg(17,i) = R[0][2]*tmp[15][i] + R[1][2]*tmp[16][i] + R[2][2]*tmp[17][i];
   }

#ifdef COMPOSITE_DEBUG
   check4<<R[0][0]<<"   "<<R[1][0]<<"   "<<R[2][0]<<endl;
   check4<<R[0][1]<<"   "<<R[1][1]<<"   "<<R[2][1]<<endl;
   check4<<R[0][2]<<"   "<<R[1][2]<<"   "<<R[2][2]<<endl;
#endif

   return kg;
}


CrdTransf *
RCFTCrdLinTransf3D::getCopy3d(void)
{
  // create a new instance of PDeltaCrdTransf3d

  RCFTCrdLinTransf3D *theCopy;

  static Vector xz(3);
  xz(0) = R[2][0];
  xz(1) = R[2][1];
  xz(2) = R[2][2];

  theCopy = new RCFTCrdLinTransf3D(this->getTag(), xz);

  theCopy->L = L;

  for (int i = 0; i < 3; i++)
	  for (int j = 0; j < 3; j++)
		  theCopy->R[i][j] = R[i][j];


  return theCopy;
}


int
RCFTCrdLinTransf3D::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

  static Vector data(7);
  data(0) = this->getTag();
  data(1) = L;
  data(2) = 0.0;
  data(3) = 0.0;
  data(4) = R[2][0];
  data(5) = R[2][1];
  data(6) = R[2][2];

  res += theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr<<"%s - failed to send Vector - RCFTCrdTransf3D::sendSelf";
    return res;
  }

  return res;
}


int
RCFTCrdLinTransf3D::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(7);

  res += theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr<<"%s - failed to receive Vector - RCFTCrdTransf3D::recvSelf";
    return res;
  }

  this->setTag((int)data(0));
  L = data(1);
  data(0) = this->getTag();
  data(1) = L;
  R[2][0] = data(4);
  R[2][1] = data(5);
  R[2][2] = data(6);

  return res;
}


const Vector &
RCFTCrdLinTransf3D::getPointGlobalCoordFromLocal(const Vector &xl)
{
   static Vector xg(3);

   xg = nodeIPtr->getCrds();

   xg(0) += R[0][0]*xl(0) + R[1][0]*xl(1) + R[2][0]*xl(2);
   xg(1) += R[0][1]*xl(0) + R[1][1]*xl(1) + R[2][1]*xl(2);
   xg(2) += R[0][2]*xl(0) + R[1][2]*xl(1) + R[2][2]*xl(2);

   return xg;
}


const Vector &
RCFTCrdLinTransf3D::getPointGlobalDisplFromBasic (double xi, const Vector &uxb)
{
   return 0;
}



void
RCFTCrdLinTransf3D::Print(OPS_Stream &s, int flag)
{
   s << "\nCrdTransf: " << this->getTag() << " Type: RCFTLinCrdTransf3D" << endln;
}








