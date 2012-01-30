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
// $Date: 2007-09-21 15:28:10 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/coordTransformation/RCFTSTLLinCrdTransf3D.cpp,v $


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

#include <RCFTSTLLinCrdTransf3D.h>

//using std::ios;
//using std::endl;
//using std::ofstream;

// constructor:
RCFTSTLLinCrdTransf3D::RCFTSTLLinCrdTransf3D(int tag, const Vector &vecInLocXZPlane):
  CrdTransf(tag, CRDTR_TAG_RCFTSTLLinCrdTransf3D),
  nodeIPtr(0), nodeJPtr(0), L(0), kl(12,12), pg(12)
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
RCFTSTLLinCrdTransf3D::RCFTSTLLinCrdTransf3D():
  CrdTransf(0, CRDTR_TAG_RCFTSTLLinCrdTransf3D),
  nodeIPtr(0), nodeJPtr(0), L(0), pg(12)
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			R[i][j] = 0.0;
}


// destructor:
RCFTSTLLinCrdTransf3D::~RCFTSTLLinCrdTransf3D()
{

}


int
RCFTSTLLinCrdTransf3D::commitState(void)
{
   return 0;
}


int
RCFTSTLLinCrdTransf3D::revertToLastCommit(void)
{
   return 0;
}


int
RCFTSTLLinCrdTransf3D::revertToStart(void)
{
   return 0;
}


int
RCFTSTLLinCrdTransf3D::initialize(Node *nodeIPointer, Node *nodeJPointer)
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
         opserr << "\nRCFTCrdTransf3D::getLocalAxes3D";
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
RCFTSTLLinCrdTransf3D::update(void)
{
   return 0;
}


int
RCFTSTLLinCrdTransf3D::computeElemtLengthAndOrient()
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
      opserr << "\nRCFTCrdTransf3d::computeElemtLengthAndOrien: 0 length\n";
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
RCFTSTLLinCrdTransf3D::getLocalAxes(Vector &XAxis, Vector &YAxis, Vector &ZAxis)
{
   XAxis(0) = R[0][0]; XAxis(1) = R[0][1]; XAxis(2) = R[0][2];
   YAxis(0) = R[1][0]; YAxis(1) = R[1][1]; YAxis(2) = R[1][2];
   ZAxis(0) = R[2][0]; ZAxis(1) = R[2][1]; ZAxis(2) = R[2][2];

   return 0;
}


double
RCFTSTLLinCrdTransf3D::getInitialLength(void)
{
   return L;
}

double
RCFTSTLLinCrdTransf3D::getDeformedLength(void)
{
   return L;
}


const Vector &
RCFTSTLLinCrdTransf3D::getBasicTrialDisp (void)
{
   static Vector ub(7);
   ub.Zero();
   return ub;
}


const Vector &
RCFTSTLLinCrdTransf3D::getBasicIncrDisp (void)
{
   static Vector ub(7);
   ub.Zero();
   return ub;
}


const Vector &
RCFTSTLLinCrdTransf3D::getBasicIncrDeltaDisp(void)
{
   static Vector ub(7);
   ub.Zero();
   return ub;
}


const Vector &
RCFTSTLLinCrdTransf3D::getBasicTrialVel(void)
{
   static Vector ub(7);
   ub.Zero();
   return ub;
}


const Vector &
RCFTSTLLinCrdTransf3D::getBasicTrialAccel(void)
{
   static Vector ub(7);
   ub.Zero();
   return ub;
}

const Vector&
RCFTSTLLinCrdTransf3D::getGlobalResistingForce(const Vector &sg, const Vector &p0)
{
  //ofstream unbal;
  //unbal.open("unbal.dat",ios::app);
  //unbal<<"\n RCFTCrdTrnsf3D::getGlobalResistingForce"<<endl;
  //unbal>>sg;
  //unbal<<"\n R vectors"<<endl;
  //unbal<<R[0][0]<<"   "<<R[0][1]<<"   "<<R[0][2]<<endl;
  //unbal<<R[1][0]<<"   "<<R[1][1]<<"   "<<R[1][2]<<endl;
  //unbal<<R[2][0]<<"   "<<R[2][1]<<"   "<<R[2][2]<<endl;

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

  //unbal<<"\n RCFTCrdTrbsf3D::after transformation "<<endl;
  //unbal>>pg;

  return pg;
}


const Matrix &
RCFTSTLLinCrdTransf3D::getInitialGlobalStiffMatrix (const Matrix &KB)
{
  static Matrix kg(12,12);	// Global stiffness for return
  static double kb[6][6];
  Matrix kl(12,12);
  static double tmp[12][12];	// Temporary storage

  double oneOverL = 1.0/L;

  static double lambda2[12][6];   //Natural to local transformation
  static double lambda1[6][12]; //Transpose of nat_to_local
  static double temp_rt[6][12];        //Temporary storage for transformation

  int i,j,k;
  for (i = 0; i < 6; i++)
	for (j = 0; j < 6; j++)
	    kb[i][j] = KB(i,j);

  /************************************************************************/
  /* CALCULATE THE NATURAL-TO-LOCAL TRANSFORMATION MATRIX 		  */
   /************************************************************************/
  lambda1[0][0]     = -1.0;
  lambda1[0][6]     = 1.0;
  lambda1[1][1]     = 1.0 / L;
  lambda1[1][5]     = 1.0;
  lambda1[1][7]     = - 1.0 / L;
  lambda1[2][2]     = 1.0 / L;
  lambda1[2][4]     = - 1.0;
  lambda1[2][8]     = - 1.0 / L;
  lambda1[3][1]     = 1.0 / L;
  lambda1[3][7]     = - 1.0 / L;
  lambda1[3][11]    = 1.0;
  lambda1[4][2]     = 1.0 / L;
  lambda1[4][8]     = - 1.0 / L;
  lambda1[4][10]    = - 1.0;
  lambda1[5][3]     = - 1.0;
  lambda1[5][9]     = 1.0;
  
  /************************************************************************/
  /* AND ITS TRANSPOSE [nat_to_local](T)		  		  */
  /************************************************************************/
  for ( i = 0; i < 6; i++ ){
    for ( k = 0; k < 12; k++ ){
		lambda2[k][i] = lambda1[i][k];
    }
  }

  /************************************************************************/
  /* MATRIX MULTIPLICATION  [k] * [nat_to_local](T)		          */
  /************************************************************************/
  for ( i = 0; i < 6; i++ ){
	for ( k = 0; k < 12; k++ ){
		for ( j = 0; j < 6; j++ ){
		temp_rt[i][k] += kb[i][j] * lambda1[j][k];
		}
	}
  }

  /*************************************************************************/
  /* MATRIX MULTIPLICATION  [nat_to_local] * [k] * [nat_to_local](T)	   */
  /*************************************************************************/
  for ( i = 0; i < 12; i++ ){
     for ( k = 0; k < 12; k++ ){
	for ( j = 0; j < 6; j++ ){
	    kl(i, k) += lambda2[i][j] * temp_rt[j][k];
	}
     }
  }

  // Transform local stiffness to global system
  // First compute kl*T_{lg}
  for (i = 0; i < 12; i++) {
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
  }

  // Now compute T'_{lg}*(kl*T_{lg})
  for (i = 0; i < 12; i++) {
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
  }
  return kg;
}


const Matrix &
RCFTSTLLinCrdTransf3D::getLocalStiffMatrix (void)
{
  return kl;	
}


const Matrix &
RCFTSTLLinCrdTransf3D::getGlobalStiffMatrix (const Matrix &KB, const Vector &fk)
{

   //ofstream output;
   //output.open("crdlocal.dat",ios::app);

   //ofstream lstiff;
   //lstiff.open("lstiff.dat",ios::app);

   //ofstream lstiff2;
   //lstiff2.open("lstiff2.dat",ios::app);
	
   static Matrix kg(12,12);	// Global stiffness for return
   static double kb[6][6];
   
   Matrix kl(12,12);
   static double tmp[12][12];	// Temporary storage

   double oneOverL = 1.0/L;

   double lambda1[6][12]; 	//Natural to local transformation
   double lambda2[12][6]; 	//Transpose of nat_to_local
   double temp_rt[6][12]; 	//Temporary storage for transformation

   int i,j,k;

   for(i=0; i < 6; i++){
        for(j=0; j < 12; j++){
	    lambda1[i][j] = 0.0;
	    lambda2[j][i] = 0.0;
	    temp_rt[i][j] = 0.0;
	} 
   }
   
   for (i = 0; i < 6; i++)
	for (j = 0; j < 6; j++)
	    kb[i][j] = KB(i,j);

   /************************************************************************/
   /* CALCULATE THE NATURAL-TO-LOCAL TRANSFORMATION MATRIX 		   */
   /************************************************************************/
   lambda1[0][0]     = -1.0;
   lambda1[0][6]     = 1.0;
   lambda1[1][1]     = 1.0 / L;
   lambda1[1][5]     = 1.0;
   lambda1[1][7]     = - 1.0 / L;
   lambda1[2][2]     = 1.0 / L;
   lambda1[2][4]     = - 1.0;
   lambda1[2][8]     = - 1.0 / L;
   lambda1[3][1]     = 1.0 / L;
   lambda1[3][7]     = - 1.0 / L;
   lambda1[3][11]    = 1.0;
   lambda1[4][2]     = 1.0 / L;
   lambda1[4][8]     = - 1.0 / L;
   lambda1[4][10]    = - 1.0;
   lambda1[5][3]     = - 1.0;
   lambda1[5][9]     = 1.0;
				  
   /************************************************************************/
   /* AND ITS TRANSPOSE [nat_to_local](T)				   */
   /************************************************************************/
   for ( i = 0; i < 6; i++ ){
        for ( k = 0; k < 12; k++ ){
		lambda2[k][i] = lambda1[i][k];
	}
   }

   /************************************************************************/
   /* MATRIX MULTIPLICATION  [k] * [nat_to_local](T)		       	   */
   /************************************************************************/
   for ( i = 0; i < 6; i++ ){
	for ( k = 0; k < 12; k++ ){
		for ( j = 0; j < 6; j++ ){
		temp_rt[i][k] += kb[i][j] * lambda1[j][k];
		}
	}
   }

   /*************************************************************************/
   /* MATRIX MULTIPLICATION  [nat_to_local] * [k] * [nat_to_local](T)	    */
   /*************************************************************************/
   for ( i = 0; i < 12; i++ ){
	for ( k = 0; k < 12; k++ ){
		for ( j = 0; j < 6; j++ ){
		    kl(i, k) += lambda2[ i ][ j ] * temp_rt[ j ][ k ];
		}
	}
   }

   for (i = 0; i < 12; i++) {
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
  }
  // Now compute T'_{lg}*(kl*T_{lg})
  for (i = 0; i < 12; i++) {
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
   }
   return kg;
}


CrdTransf *
RCFTSTLLinCrdTransf3D::getCopy(void)
{
  // create a new instance of PDeltaCrdTransf3d

  RCFTSTLLinCrdTransf3D *theCopy;

  static Vector xz(3);
  xz(0) = R[2][0];
  xz(1) = R[2][1];
  xz(2) = R[2][2];

  theCopy = new RCFTSTLLinCrdTransf3D(this->getTag(), xz);

  theCopy->L = L;

  for (int i = 0; i < 3; i++)
	  for (int j = 0; j < 3; j++)
		  theCopy->R[i][j] = R[i][j];


  return theCopy;
}


int
RCFTSTLLinCrdTransf3D::sendSelf(int cTag, Channel &theChannel)
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
RCFTSTLLinCrdTransf3D::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
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
RCFTSTLLinCrdTransf3D::getPointGlobalCoordFromLocal(const Vector &xl)
{
   static Vector xg(3);

   //xg = nodeIPtr->getCrds() + nodeIOffset;
   xg = nodeIPtr->getCrds();

   // xg = xg + Rlj'*xl
   //xg.addMatrixTransposeVector(1.0, Rlj, xl, 1.0);
   xg(0) += R[0][0]*xl(0) + R[1][0]*xl(1) + R[2][0]*xl(2);
   xg(1) += R[0][1]*xl(0) + R[1][1]*xl(1) + R[2][1]*xl(2);
   xg(2) += R[0][2]*xl(0) + R[1][2]*xl(1) + R[2][2]*xl(2);

   return xg;
}


const Vector &
RCFTSTLLinCrdTransf3D::getPointGlobalDisplFromBasic (double xi, const Vector &uxb)
{
   return 0;
}



void
RCFTSTLLinCrdTransf3D::Print(OPS_Stream &s, int flag)
{
   s << "\nCrdTransf: " << this->getTag() << " Type: RCFTCrdTransf3D" << endln;

}








