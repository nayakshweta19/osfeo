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
**   Cenk Tort (tort0008@umn.edu)   				      **
**   Jerome F. Hajjar (hajjar@struc.ce.umn.edu)                       **
**                                                                    **
**   University of Minnesota -- Twin Cities                           **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.2 $
// $Date: 2008-07-03 18:03:49 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/element/RCFTSTLBeamColumn3D.cpp,v $



#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <iomanip>

#include <Information.h>
#include <Parameter.h>
#include <RCFTSTLBeamColumn3D.h>
#include <MatrixUtil.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <math.h>
#include <G3Globals.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <iostream>
#include <fstream>

#define  NDM   3         // dimension of the problem (3d)
#define  NL    2         // size of uniform load vector
#define  NND   6         // number of nodal dof's
#define  NEGD  12        // number of element global dof's
#define  NEBD  6        // number of element dof's in the basic system

//using std::endl;
//using std::ios;
//using std::ifstream;
//using std::ofstream;
using namespace std;

Matrix RCFTSTLBeamColumn3D::theMatrix(12,12);
Vector RCFTSTLBeamColumn3D::theVector(12);
double RCFTSTLBeamColumn3D::workArea[400];

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
// CONSTRUCTOR FOR PARALLEL PROCESSING

RCFTSTLBeamColumn3D::RCFTSTLBeamColumn3D():
Element(0,ELE_TAG_RCFTSTLBeamColumn3D), connectedExternalNodes(2),
beamIntegr(0), numSections(0), sections(0), crdTransf(0),
rho(0.0), initialFlag(0), kv(NEBD,NEBD), Se(NEBD), Sg(NEGD), df_i(NEGD), Sglobal(NEGD), kvcommit(NEBD,NEBD), Secommit(NEBD),
ks(0), Ki(0), XAxis(3), YAxis(3), ZAxis(3), sr(1,6), ss(1,1), fk(6), fk_incr(6), ub(7), Tagg(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;

  kv.Zero();
  sr.Zero();
  ss.Zero();
  fk.Zero();
  fk_incr.Zero();
  ub.Zero();
  df_i.Zero();
  deflength = 0.0;

}


// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points.
// allocates the necessary space needed by each object
RCFTSTLBeamColumn3D::RCFTSTLBeamColumn3D (int tag, int nodeI, int nodeJ,
                                      int numSec, RCFTSTLFiberSection3D **sec,
				      BeamIntegration &bi,
                                      RCFTSTLCrdTransf3D &coordTransf, double massDensPerUnitLength):
Element(tag,ELE_TAG_RCFTSTLBeamColumn3D), connectedExternalNodes(2),
beamIntegr(0), numSections(0), sections(0), crdTransf(0),
rho(massDensPerUnitLength),
initialFlag(0),
kv(NEBD,NEBD), Se(NEBD), Sg(NEGD), df_i(NEGD), Sglobal(NEGD), kvcommit(NEBD,NEBD), Secommit(NEBD),
ks(0), Ki(0), XAxis(3), YAxis(3), ZAxis(3), sr(1,6), ss(1,1), fk(6), fk_incr(6), ub(7), Tagg(tag)
{

   theNodes[0] = 0;
   theNodes[1] = 0;

   connectedExternalNodes(0) = nodeI;
   connectedExternalNodes(1) = nodeJ;

   // get copy of the sections

   beamIntegr = bi.getCopy();
   if(beamIntegr == 0) {
     opserr<<"Error: RCFTBeamColumn3D::RCFTBeamColumn3D: could not create copy of beam integration object" << endln;
     exit(-1);
   }

   // get copy of the transformation object
   
   //crdTransf = coordTransf.getCopy();

   crdTransf = (RCFTSTLCrdTransf3D*) coordTransf.getCopy();

   //deflength = crdTransf->getInitialLength();
   
   if (crdTransf == 0) {
      opserr << "Error: RCFTBeamColumn3D::RCFTBeamColumn3D: could not create copy of coordinate transformation object" << endln;
      exit(-1);
   }

   this->setSectionPointers(numSec,sec);

   //Fill in transformation matrix
   R[0][0]   = 0.0; R[0][1] = 0.0; R[0][2] = 0.0;
   R[1][0]   = 0.0; R[1][1] = 0.0; R[1][2] = 0.0;
   R[2][0]   = 0.0; R[2][1] = 0.0; R[2][2] = 0.0;

   ss.Zero();
   sr.Zero();
   
   kv.Zero();
   fk.Zero();
   df_i.Zero();
   fk_incr.Zero();
   ub.Zero();
   deflength = 0.0;
}

// ~RCFT_BeamColumn3D():
// 	destructor
//  delete must be invoked on any objects created by the object
RCFTSTLBeamColumn3D::~RCFTSTLBeamColumn3D()
{

   if (sections) {
      for (int i=0; i < numSections; i++)
         if (sections[i])
            delete sections[i];
      delete [] sections;
   }

   if(ks != 0)
     delete [] ks;

   if(ksa != 0)
     delete [] ksa;

   if(dhat != 0)
     delete [] dhat;

   if(DSQ != 0)
     delete [] DSQ;

   if(DSQa != 0)
     delete [] DSQa;
  
   if(f3 != 0)
     delete [] f3;

   if(d3 != 0)
     delete [] d3;

   if(str_f3 != 0) 
     delete [] str_f3;

   if(str_f3inv != 0)
     delete [] str_f3inv;

   if (crdTransf)
     delete crdTransf;

   if(beamIntegr != 0)
     delete beamIntegr;

   if (Ki != 0)
     delete Ki;

}


int
RCFTSTLBeamColumn3D::getNumExternalNodes(void) const
{
   return 2;
}


const ID &
RCFTSTLBeamColumn3D::getExternalNodes(void)
{
   return connectedExternalNodes;
}


Node **
RCFTSTLBeamColumn3D::getNodePtrs(void)
{
   return theNodes;
}

int
RCFTSTLBeamColumn3D::getNumDOF(void)
{
   return NEGD;
}


void
RCFTSTLBeamColumn3D::setDomain(Domain *theDomain)
{
   // check Domain is not null - invoked when object removed from a domain
 
   if (theDomain == 0)
   {
      theNodes[0] = 0;
      theNodes[1] = 0;

      opserr << "RCFTSTLBeamColumn3D::setDomain:  theDomain = 0 ";
      exit(0);
   }

   // get pointers to the nodes

   int Nd1 = connectedExternalNodes(0);
   int Nd2 = connectedExternalNodes(1);

   theNodes[0] = theDomain->getNode(Nd1);
   theNodes[1] = theDomain->getNode(Nd2);

   if (theNodes[0] == 0)
   {
      opserr << "RCFTSTLBeamColumn3D::setDomain: Nd1: ";
      opserr << Nd1 << "does not exist in model\n";
      exit(0);
   }

   if (theNodes[1] == 0)
   {
      opserr << "RCFTSTLBeamColumn3D::setDomain: Nd2: ";
      opserr << Nd2 << "does not exist in model\n";
      exit(0);
   }

   // call the DomainComponent class method
   this->DomainComponent::setDomain(theDomain);

   // ensure connected nodes have correct number of dof's
   int dofNode1 = theNodes[0]->getNumberDOF();
   int dofNode2 = theNodes[1]->getNumberDOF();

   if ((dofNode1 != NND) || (dofNode2 != NND))
   {
      opserr << "RCFTSTLBeamColumn3D::setDomain(): Nd2 or Nd1 incorrect dof ";
      exit(0);
   }

   // initialize the transformation
   if (crdTransf->initialize(theNodes[0], theNodes[1]))
   {
      opserr << "RCFTBeamColumn3D::setDomain(): Error initializing coordinate transformation";
      exit(0);
   }

   deflength = crdTransf->getInitialLength();

   crdTransf->getLocalAxes(XAxis, YAxis, ZAxis);

   R[0][0] = XAxis(0); R[0][1] = XAxis(1); R[0][2] = XAxis(2);
   R[1][0] = YAxis(0); R[1][1] = YAxis(1); R[1][2] = YAxis(2);
   R[2][0] = ZAxis(0); R[2][1] = ZAxis(1); R[2][2] = ZAxis(2);

   // get element length
   double L = crdTransf->getInitialLength();
   if (L == 0.0)
   {
      opserr << "RCFTBeamColumn3D::setDomain(): Zero element length:" << this->getTag();
      exit(0);
   }

   if (initialFlag == 0)
     this->initializeSectionHistoryVariables();

}



int
RCFTSTLBeamColumn3D::commitState()
{
   ofstream dunhat;
   dunhat.open("dunhat.dat",ios::app); 

   ofstream output;
   output.open("stlcon.dat",ios::app);
   
   int err = 0;
   int i = 0;
   int j = 0;

   // call element commitState to do any base class stuff
   if ((err = this->Element::commitState()) != 0) {
     opserr << "RCFTBeamColumn3D::commitState () - failed in base class";
     return err;
   }

   do {
      output<<"section #"<<i<<endl;	   
      err = sections[i++]->commitState();

   } while (err == 0 && i < numSections);

   if (err)
      return err;

   // commit the transformation between coord. systems
   if ((err = crdTransf->commitState()) != 0)
      return err;

   // commit the element variables state
   kvcommit = kv;
   Secommit = Se;

   ub.Zero();
   return err;
}


int RCFTSTLBeamColumn3D::revertToLastCommit()
{
   int err;
   int i = 0;

   do {
      err = sections[i]->revertToLastCommit();

      i++;
   } while (err == 0 && i < numSections);


   if (err)
      return err;

   // revert the transformation to last commit
   if ((err = crdTransf->revertToLastCommit()) != 0)
      return err;

   // revert the element state to last commit
   Se   = Secommit;
   kv   = kvcommit;

   initialFlag = 0;
   // this->update();

   return err;
}


int RCFTSTLBeamColumn3D::revertToStart()
{
   // revert the sections state to start
   int err;
   int i = 0;

   do {
       err = sections[i++]->revertToStart();

   }while (err == 0 && i < numSections);

   if (err)
      return err;

   // revert the transformation to start
   if ((err = crdTransf->revertToStart()) != 0)
      return err;

   // revert the element state to start
   Secommit.Zero();
   kvcommit.Zero();

   Se.Zero();
   kv.Zero();

   initialFlag = 0;
   // this->update();
   return err;
}

const Matrix &
RCFTSTLBeamColumn3D::getInitialStiff(void)
{
  ofstream newton; 
  newton.open("newton.dat",ios::app);

  ofstream lstiff;
  lstiff.open("lstiff.dat",ios::app);
 
  if (Ki != 0)
    return *Ki;

  int i, j, k;

  double L = crdTransf->getInitialLength();
  double oneOverL  = 1.0/L;

  /********************************************************/
  /*   GET CROSS-SECTION STIFFNESS AND FORCES FOR THE     */
  /*			  I-END AND J-END	   	  */
  /********************************************************/
  const Matrix &kcrsi = sections[0]->getInitialTangent();
  const Matrix &kcrsj = sections[1]->getInitialTangent();

  Matrix kvInit(6,6);

  kvInit.Zero();

  Matrix kt(7,7);
  Matrix temp_sr(1,6);
  Matrix temp_rr(6,6);

  /************************************************************************/
  /* TO FACILITATE SIMPLER SYNTAX IN THE EQUATIONS, DEFINE SOME LOCAL	  */
  /* VARIABLES TO REPRESENT RIGIDITY TERMS				  */
  /************************************************************************/
  /* FORCE TERMS								 */
  double	p 	= 0.0;	/* AVERAGE AXIAL FORCE 		 		 */
  double	p_diff	= 0.0;	/* AVERAGE AXIAL FORCE 		 	 	 */
  double	pa 	= 0.0;	/* AXIAL FORCE IN AT END A			 */
  double	pb 	= 0.0;	/* AXIAL FORCE IN AT END B			 */
  double	my 	= 0.0;	/* Y_AXIS MOMENT 				 */
  double	mya 	= 0.0;	/* Y_AXIS MOMENT 				 */
  double	myb 	= 0.0;	/* Y_AXIS MOMENT 				 */
  double	mz 	= 0.0;	/* Z_AXIS MOMENT 				 */
  double	mza 	= 0.0;	/* Z_AXIS MOMENT 				 */
  double	mzb 	= 0.0;	/* Z_AXIS MOMENT 				 */
  double	lzz 	= 0.0;	/* ZZ_AXIS SECOND MOMENT OF THE AREA  		 */
  double	lyy 	= 0.0;	/* YY_AXIS SECOND MOMENT OF THE AREA  		 */
  double	lyz 	= 0.0;	/* YZ_AXIS SECOND MOMENT OF THE AREA 		 */
  double	dmy 	= 0.0;	/* Y_AXIS MOMENT 				 */
  double	dmz 	= 0.0;	/* Z_AXIS MOMENT 		 		 */
  double	dlzz 	= 0.0;	/* ZZ_AXIS SECOND MOMENT OF THE AREA  		 */
  double	dlyy 	= 0.0;	/* YY_AXIS SECOND MOMENT OF THE AREA  		 */
  double	dlyz 	= 0.0;	/* YZ_AXIS SECOND MOMENT OF THE AREA E 		 */

  /* RIGIDITY TERMS										 */
  double	gj 	= 0.0;	/* TORSIONAL RIGIDITY						 */

  double	ea 	= 0.0;	/* AXIAL RIGIDITY IN CONCRETE					 */
  double	eqy	= 0.0;	/* FIRST MOMENT RIGIDITY IN CONCRETE			 	 */
  double	eqz	= 0.0;	/* FIRST MOMENT RIGIDITY IN CONCRETE			 	 */
  double	eiyy	= 0.0;	/* MOM. OF INERTIA IN CONCRETE					 */
  double	eizz	= 0.0;	/* MOM. OF INERTIA IN CONCRETE			 		 */
  double	eiyz	= 0.0;	/* MOM. OF INERTIA IN CONCRETE					 */
  double	dea 	= 0.0;	/* DIFF. BETWEEN i & j AXIAL RIGIDITY IN CONCRETE		 */
  double	deqy	= 0.0;	/* DIFF IN FIRST MOMENT RIGIDITY IN CONCRETE	 		 */
  double	deqz	= 0.0;	/* DIFF IN FIRST MOMENT RIGIDITY IN CONCRETE	 		 */
  double	deiyy	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN CONCRETE			 	 */
  double	deizz	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN CONCRETE			 	 */
  double	deiyz	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN CONCRETE			 	 */

  /************************************************************************/
  /* CALCULATE THE VALUES OF THE LOCAL VARIABLES			  */
  /************************************************************************/
  /* FORCE TERMS							  */
  /************************************************************************/
  /* NOTE THAT THESE TERMS ARE INTERNAL FORCES THAT FOLLOW THE INTERNAL	  */
  /* FORCE POSITIVE SIGN CONVENTION.  THEY DO NOT REPRESENT THE POSITIVE  */
  /* ELEMENT END FORCE SIGN CONVENTION.					  */
  /************************************************************************/
  Matrix rs(6,1);
  Matrix rr(6,6);
  
  p 	= ( fk(0) + fk(3) ) / 2.0 ;
  pa 	= fk(0);
  pb 	= fk(3);
  my 	= fk(2);
  mya = fk(2);
  myb = fk(5);
  mz 	= fk(1);
  mza = fk(1);
  mzb = fk(4);
  dmy = myb - mya;
  dmz = mzb - mza;

  gj 	 = kcrsi(3,3);

  ea 	 = kcrsi(0,0);
  eqy	 = kcrsi(0,2);
  eqz	 = kcrsi(1,0);
  eiyy   = kcrsi(2,2);
  eizz   = kcrsi(1,1);
  eiyz   = kcrsi(1,2);
  dea    = kcrsj(0,0) - kcrsi(0,0);
  deqy   = kcrsj(0,2) - kcrsi(0,2);
  deqz   = kcrsj(1,0) - kcrsi(1,0);
  deiyy  = kcrsj(2,2) - kcrsi(2,2);
  deizz  = kcrsj(1,1) - kcrsi(1,1);
  deiyz  = kcrsj(1,2) - kcrsi(1,2);

  /*****************************************************************************************/
  /* INITIALIZE THE NATURAL ELEMENT STIFFNESS MATRIX AND THE MATRICES USED IN CONDENSATION */
  /*****************************************************************************************/
  kt.Zero();
  ss.Zero();
  sr.Zero();
  temp_sr.Zero();
  rs.Zero();
  rr.Zero();
  temp_rr.Zero();

  /************************************************************/
  /* GENERATE ELASTIC TERMS IN UPPER TRIANGULAR PORTION OF kt */
  /************************************************************/

  kt(0,0)  = ( 7.0 * ea / ( 3.0 * L ) ) + ( 11.0 * dea / ( 6.0 * L ) );
  kt(0,1)  = ( eqz / L ) + ( 2.0 * deqz / ( 3.0 * L ) );
  kt(0,2)  = ( eqy / L ) + ( 2.0 * deqy / ( 3.0 * L ) );
  kt(0,3)  = ( 3.0 * eqz / L ) + ( 7.0 * deqz / ( 3.0 * L ) );
  kt(0,4)  = ( 3.0 * eqy / L ) + ( 7.0 * deqy / ( 3.0 * L ) );
  kt(0,6)  = - ( 8.0 * ea / ( 3.0 * L ) ) - ( 2.0 * dea / L );
  kt(1,1)  = ( 4.0 * eizz / L ) + ( deizz / L );
  kt(1,2)  = ( 4.0 * eiyz / L ) + ( deiyz / L );
  kt(1,3)  = ( 2.0 * eizz / L ) + ( deizz / L );
  kt(1,4)  = ( 2.0 * eiyz / L ) + ( deiyz / L );
  kt(1,6)  = - ( 4.0 * eqz / L ) - ( 4.0 * deqz / ( 3.0 * L ) );
  kt(2,2)  = ( 4.0 * eiyy / L ) + ( deiyy / L );
  kt(2,3)  = ( 2.0 * eiyz / L ) + ( deiyz / L );
  kt(2,4)  = ( 2.0 * eiyy / L ) + ( deiyy / L );
  kt(2,6)  = - ( 4.0 * eqy / L ) - ( 4.0 * deqy / ( 3.0 * L ) );
  kt(3,3)  = ( 4.0 * eizz / L ) + ( 3.0 * deizz / L );
  kt(3,4)  = ( 4.0 * eiyz / L ) + ( 3.0 * deiyz / L );
  kt(3,6)  = - ( 4.0 * eqz / L ) - ( 8.0 * deqz / ( 3.0 * L ) );
  kt(4,4)  = ( 4.0 * eiyy / L ) + ( 3.0 * deiyy / L );
  kt(4,6)  = - ( 4.0 * eqy / L ) - ( 8.0 * deqy / ( 3.0 * L ) );
  kt(5,5)  = gj / L;
  kt(6,6)  = ( 16.0 * ea / ( 3.0 * L ) ) + ( 8.0 * dea / ( 3.0 * L ) );
  kt(1,1)  += 2.0 * p * L / 15.0;
  kt(1,3)  += - p * L / 30.0;
  kt(2,2)  += 2.0 * p * L / 15.0;
  kt(2,4)  += - p * L / 30.0;
  kt(3,3)  += 2.0 * p * L / 15.0;
  kt(4,4)  += 2.0 * p * L / 15.0;
  

  kt(0,0)  += ( 7.0 * p / ( 3.0 * L ) );
  kt(0,1)  += ( mz / L ) + ( 2.0 * dmz / ( 3.0 * L ) );
  kt(0,2)  += ( my / L ) + ( 2.0 * dmy / ( 3.0 * L ) );
  kt(0,3)  += ( 3.0 * mz / L ) + ( 7.0 * dmz / ( 3.0 * L ) );
  kt(0,4)  += ( 3.0 * my / L ) + ( 7.0 * dmy / ( 3.0 * L ) );
  kt(0,6)  -= ( 8.0 * p / ( 3.0 * L ) );
  kt(1,5)  += - dmy / 12.0;
  kt(1,6)  += - ( 4.0 * mz / L ) - ( 4.0 * dmz / ( 3.0 * L ) );
  kt(2,5)  += dmz / 12.0;
  kt(2,6)  -= ( 4.0 * my / L ) + ( 4.0 * dmy / ( 3.0 * L ) );
  kt(3,5)  += dmy / 12.0;
  kt(3,6)  += - ( 4.0 * mz / L ) - ( 8.0 * dmz / ( 3.0 * L ) );
  kt(4,5)  += - dmz / 12.0;
  kt(4,6)  -= ( 4.0 * my / L ) + ( 8.0 * dmy / ( 3.0 * L ) );
  kt(6,6)  += ( 16.0 * p / ( 3.0 * L ) );
  
  /****************************************************/
  /* GENERATE TERMS IN LOWER TRIANGULAR PORTION OF kt */
  /****************************************************/
  for ( i = 0; i < 6; i++ ){
	for ( j = ( i + 1 ); j < 7; j++ ){
		kt(j,i) = kt(i,j);
	}
  }

  /************************************************************************/
  /*									  */
  /* STATIC CONDENSATION PROCEDURE - SEE DON WHITES THESIS		  */
  /*									  */
  /*	| |rr| {rs} |   { D  }   	{F}				  */
  /*	|	    | 		  = 		  		  	  */
  /*	| {sr} |ss| |   { Dm }   	{0}				  */
  /*									  */
  /* so condensed set of equations can be written as:			  */
  /*									  */
  /* |k_condensed| { D } = { F }					  */
  /*									  */
  /* where:								  */
  /*									  */
  /* |k_condensed| = |rr| - {rs} inv|ss| {sr}				  */
  /*									  */
  /************************************************************************/
  double ssdet;

  /*************************************************************/
  /* GENERATE TERMS IN OTHER MATRICIES FOR CONDENSATION	       */
  /*************************************************************/
  ssdet = kt(6,6);

  ss(0,0) = 1.0 / ssdet;

  /*************************************************************/
  /* GENERATE TERMS IN OTHER MATRICIES FOR CONDENSATION        */
  /*************************************************************/

  for ( i = 0; i < 6; i++ ){
         for ( j = 0; j < 1; j++ ){
              rs(i,j) = kt(i,6);
              sr(j,i) = rs(i,j);
         }
         for ( j = 0; j < 6; j++ ){
	      rr(i,j) = kt(i,j);
         }
  }

  /*********************************************/
  /* GENERATE CONDENSED NATURAL ELEMENT MATRIX */
  /*********************************************/

  /* MATRIX MULTIPLICATION */
  for ( i = 0; i < 1; i++ )
  {
         for ( k = 0; k < 6; k++ )
         {
                for ( j = 0; j < 1; j++ )
                {
                        temp_sr(i,k) += ss(i,j) * sr(j,k);
                }
         }
  }

  /* MATRIX MULTIPLICATION */
  for ( i = 0; i < 6; i++ )
  {
         for ( k = 0; k < 6; k++ )
         {
                for ( j = 0; j < 1; j++ )
                {
                        temp_rr(i,k) += rs(i,j) * temp_sr(j,k);
                }
	  /* FINAL FORM OF THE CONDENSED NATURAL ELEMENT STIFFNESS MATRIX */
                kvInit(i,k) = kt(i,k) - temp_rr(i,k);
          }
  }

  Ki = new Matrix(crdTransf->getInitialGlobalStiffMatrix(kvInit));

  return *Ki;

}

const Matrix &
RCFTSTLBeamColumn3D::getTangentStiff(void)
{
  int i;  
  crdTransf->update();  // Will remove once we clean up the corotational 2d transformation -- MHS
  const Matrix &KV = crdTransf->getGlobalStiffMatrix(kv,fk);
  ofstream output;
  output.open("tangentstiff.dat", ios::app);
  output<<"\n number of element"<<Tagg<<endl;

  //for(i=0; i<12; i++){
  //  output<<kv(i,0)<<"  "<<kv(i,1)<<"  "<<kv(i,2)<<"  "<<kv(i,3)<<"  "<<kv(i,4)<<"   "
  //        <<kv(i,5)<<"  "<<kv(i,6)<<"  "<<kv(i,7)<<"  "<<kv(i,8)<<"  "<<kv(i,9)<<"   "
  //        <<kv(i,10)<<"  "<<kv(i,11)<<endl;
  //}
	
  return KV;
}

const Vector &
RCFTSTLBeamColumn3D::getResistingForce(void)
{
  return Sglobal;
}

void RCFTSTLBeamColumn3D::calcResistingForce(void)
{
  //ofstream intforce;
  //intforce.open("intforce.dat",ios::app);
  crdTransf->update();
  Vector p0(12);
  p0.Zero();
  Sg = crdTransf->getGlobalResistingForce(df_i, p0);
  Sglobal  = Sglobal + Sg;
  //intforce>>Sglobal;
}

void
RCFTSTLBeamColumn3D::initializeSectionHistoryVariables (void)
{
    for (int i = 0; i < numSections; i++){
	ks[i]       = Matrix(3,3);
        ksa[i]      = Matrix(4,4);
        dhat[i]     = Vector(3);
        DSQ[i]      = Vector(3);
        DSQa[i]     = Vector(4);
        f3[i]       = Vector(3);
        d3[i]       = Vector(3);
        str_f3[i]   = Matrix(3,3);
        str_f3inv[i]= Matrix(3,3); 
    }
}

/********* NEWTON , SUBDIVIDE AND INITIAL ITERATIONS *********************/
int RCFTSTLBeamColumn3D::update()
{
  //ofstream geom;
  //geom.open("geom.dat",ios::app);

  //ofstream unbal;
  //unbal.open("unbal.dat",ios::app);

  //ofstream mpls;
  //mpls.open("mpls.dat",ios::app);

  //ofstream lstiff;
  //lstiff.open("lstiff.dat",ios::app);

  int i,j,k;

  this->calcDeformedLength();

  double L = getDeformedLength();

  double oneOverL  = 1.0/L;

  if( initialFlag == 2 )
  this->revertToLastCommit();

  const Vector du = getLocalIncrDeltaDisp();

  //mpls<<"du"<<endl;

  //mpls>>du;

  const Matrix &KL = crdTransf->getLocalStiffMatrix(kv,fk);
  df_i = KL * du;

  //mpls<<"df_i"<<endl;
  //mpls>>df_i;

  /* CALCULATE THE INCREMENTAL ELONGATIONS IN THE X, Y, AND Z GLOBAL SPACE*/
  double d6_0, d7_1, d8_2;
  d6_0 = du(6) - du(0);
  d7_1 = du(7) - du(1);
  d8_2 = du(8) - du(2);

  /************************************************************************/
  /* CALCULATE THE ELONGATION OF STEEL           			  */
  /************************************************************************/

  Vector dub(7);  

  dub(0) = (  2.0 * L * d6_0 + pow( d6_0 , 2 ) + pow( d7_1 , 2 ) + pow( d8_2 , 2 ) ) / ( L + L );

  double theta_rigid_z = ( du(7) - du(1) ) / L;
  double theta_rigid_y = ( du(8) - du(2) ) / L;

  dub(1) = du(5) - theta_rigid_z;
  dub(2) = - du(4) - theta_rigid_y;
  dub(3) = du(11) - theta_rigid_z;
  dub(4) = - du(10) - theta_rigid_y;
  dub(5) = du(9) - du(3);

  /* FIRST CALCULATE THE INCREMENT IN AXIAL FORCE AND ADD THIS TO THE     */
  /* GEOMETRIC STIFFNESS MATRIX BEFORE RECOVERING MOMENTS                 */
  /* ONLY IF GEOMETRICALLY NONLINEAR                                      */
  Vector df_nat(6);

  df_nat(0)    = 0.0;

  for ( int i = 0; i < 6; i++ ){
      df_nat(0) += dub(i) * kv(0,i);
  }
  
  kv(1,1) += 2.0 * df_nat(0) * L / 15.0;
  kv(1,3) += - df_nat(0) * L / 30.0;
  kv(2,2) += 2.0 * df_nat(0) * L / 15.0;
  kv(2,4) += - df_nat(0) * L / 30.0;
  kv(3,3) += 2.0 * df_nat(0) * L / 15.0;
  kv(3,1) += - df_nat(0) * L / 30.0;
  kv(4,4) += 2.0 * df_nat(0) * L / 15.0;
  kv(4,2) += - df_nat(0) * L / 30.0;

  for ( int i = 0; i < 6; i++ ){
      df_nat(i)  = 0.0;
      for ( int j = 1; j < 6; j++ ){
           df_nat(i) += dub(j) * kv(i,j);
      }
  }

  /************************************************************************/
  /* COMPUTE THE NATURAL DEFORMATION OF THE MIDPOINT DEGREE OF FREEDOM.   */
  /* THIS IS THE ELONGATION OF THE MEMBER BETWEEN THE MIDPOINT AND THE    */
  /* i-end OF THE MEMBER.  THIS IS ACCOMPLISHED BY UNDOING THE STATIC     */
  /* CONDENSATION THAT WAS DONE TO THE NATURAL STIFFNESS MATRIX IN THE    */
  /* FUNCTION b_stl_tangent_k.c                                           */
  /************************************************************************/
  dub(6) = 0.0;

  for ( int i = 0; i < 6; i++ ){
        dub(6) -=  sr(0,i) * dub(i) * ss(0,0);
  }

  df_i(0) = - df_nat(0);
  df_i(1) = ( df_nat(1) + df_nat(3) - fk(0) * d7_1 ) / L;
  df_i(2) = ( df_nat(2) + df_nat(4) - fk(0) * d8_2 ) / L;
  df_i(4)  = - df_nat(2);
  df_i(5)  = df_nat(1);
  df_i(6)  = - df_i(0);
  df_i(7)  = - df_i(1);
  df_i(8)  = - df_i(2);
  df_i(10) = - df_nat(4);
  df_i(11) = df_nat(3);

  for ( i = 0; i < numSections; i++){
    ksa[i] = sections[i]->getSectionTangent();
    for( j = 0; j < 3; j++ ){
       for( k = 0; k < 3; k++ ){
          ks[i](i,j) = ksa[i](i,j);
       }
    }
    str_f3[i](0,0) = ks[i](0,0);
    str_f3[i](0,1) = ks[i](0,1);
    str_f3[i](0,2) = ks[i](0,2);

    str_f3[i](1,0) = ks[i](1,0);
    str_f3[i](1,1) = ks[i](1,1);
    str_f3[i](1,2) = ks[i](1,2);

    str_f3[i](2,0) = ks[i](2,0);
    str_f3[i](2,1) = ks[i](2,1);
    str_f3[i](2,2) = ks[i](2,2);

    invertMatrix(2,str_f3[i],str_f3inv[i]);

    if( i == 0 ){
    f3[i](0) =  df_i(0);
    f3[i](1) = - df_i(1);
    f3[i](2) = - df_i(2);
    }
    else if( i == 1 ){
    f3[i](0) =  df_i(0);
    f3[i](1) =  df_i(3);
    f3[i](2) =  df_i(4);
    }
    d3[i] = str_f3inv[i] * f3[i];
    dhat[i](0) = d3[i](0);
    dhat[i](1) = d3[i](2);
    dhat[i](2) = d3[i](3);

    Vector aggdhat(4);

    aggdhat(0) = dhat[i](0);
    aggdhat(1) = dhat[i](1);
    aggdhat(2) = dhat[i](2);

    if( i == 0){
      aggdhat(3) = du(3);
    }
    if( i == 1){
      aggdhat(3) = du(9);
    }
    int res = sections[i]->setTrialSectionDeformation(aggdhat);
    DSQa[i] = sections[i]->getStressResultant();
    for( j = 0; j < 3; j++ ){
      DSQ[i](j) = DSQa[i](j);
    }
  }
 
  this->calcResistingForce();

  crdTransf->update();

  crdTransf->getLocalAxes(XAxis, YAxis, ZAxis);

  R[0][0]   = XAxis(0);
  R[0][1]   = XAxis(1);
  R[0][2]   = XAxis(2);

  R[1][0]   = YAxis(0);
  R[1][1]   = YAxis(1);
  R[1][2]   = YAxis(2);

  R[2][0]   = ZAxis(0);
  R[2][1]   = ZAxis(1);
  R[2][2]   = ZAxis(2);

  fk(0) = DSQ[0](0);
  fk(1) = DSQ[0](1);
  fk(2) = DSQ[0](2);
  fk(3) = DSQ[numSections-1](0);
  fk(4) = DSQ[numSections-1](1);
  fk(5) = DSQ[numSections-1](2);

  /********************************************************/
  /*   GET CROSS-SECTION STIFFNESS AND FORCES FOR THE     */
  /*			  I-END AND J-END	   	  */
  /********************************************************/
  
  Matrix &kcrsi = ksa[0];
  Matrix &kcrsj = ksa[numSections-1];

  Matrix kt(7,7);
  Matrix temp_sr(1,6);
  Matrix temp_rr(6,6);

  /************************************************************************/
  /* TO FACILITATE SIMPLER SYNTAX IN THE EQUATIONS, DEFINE SOME LOCAL	  */
  /* VARIABLES TO REPRESENT RIGIDITY TERMS				  */
  /************************************************************************/
  /* FORCE TERMS								 */
  double	p 	= 0.0;	/* AVERAGE AXIAL FORCE 		 		 */
  double	p_diff	= 0.0;	/* AVERAGE AXIAL FORCE 		 	 	 */
  double	pa 	= 0.0;	/* AXIAL FORCE IN AT END A			 */
  double	pb 	= 0.0;	/* AXIAL FORCE IN AT END B			 */
  double	my 	= 0.0;	/* Y_AXIS MOMENT 				 */
  double	mya 	= 0.0;	/* Y_AXIS MOMENT 				 */
  double	myb 	= 0.0;	/* Y_AXIS MOMENT 				 */
  double	mz 	= 0.0;	/* Z_AXIS MOMENT 				 */
  double	mza 	= 0.0;	/* Z_AXIS MOMENT 				 */
  double	mzb 	= 0.0;	/* Z_AXIS MOMENT 				 */
  double	dmy 	= 0.0;	/* Y_AXIS MOMENT 				 */
  double	dmz 	= 0.0;	/* Z_AXIS MOMENT 		 		 */

  /* RIGIDITY TERMS										 */
  double	gj 	= 0.0;	/* TORSIONAL RIGIDITY						 */

  double	ea 	= 0.0;	/* AXIAL RIGIDITY IN CONCRETE					 */
  double	eqy	= 0.0;	/* FIRST MOMENT RIGIDITY IN CONCRETE			 	 */
  double	eqz	= 0.0;	/* FIRST MOMENT RIGIDITY IN CONCRETE			 	 */
  double	eiyy	= 0.0;	/* MOM. OF INERTIA IN CONCRETE					 */
  double	eizz	= 0.0;	/* MOM. OF INERTIA IN CONCRETE			 		 */
  double	eiyz	= 0.0;	/* MOM. OF INERTIA IN CONCRETE					 */
  double	dea 	= 0.0;	/* DIFF. BETWEEN i & j AXIAL RIGIDITY IN CONCRETE		 */
  double	deqy	= 0.0;	/* DIFF IN FIRST MOMENT RIGIDITY IN CONCRETE	 		 */
  double	deqz	= 0.0;	/* DIFF IN FIRST MOMENT RIGIDITY IN CONCRETE	 		 */
  double	deiyy	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN CONCRETE			 	 */
  double	deizz	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN CONCRETE			 	 */
  double	deiyz	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN CONCRETE			 	 */

  /************************************************************************/
  /* CALCULATE THE VALUES OF THE LOCAL VARIABLES			  */
  /************************************************************************/
  /* FORCE TERMS							  */
  /************************************************************************/
  /* NOTE THAT THESE TERMS ARE INTERNAL FORCES THAT FOLLOW THE INTERNAL	  */
  /* FORCE POSITIVE SIGN CONVENTION.  THEY DO NOT REPRESENT THE POSITIVE  */
  /* ELEMENT END FORCE SIGN CONVENTION.					  */
  /************************************************************************/
  Matrix rs(6,1);
  Matrix rr(6,6);
  
  p 	= ( fk(0) + fk(3) ) / 2.0 ;
  pa 	= fk(0);
  pb 	= fk(3);
  my 	= fk(2);
  mya = fk(2);
  myb = fk(5);
  mz 	= fk(1);
  mza = fk(1);
  mzb = fk(4);
  dmy = myb - mya;
  dmz = mzb - mza;

  gj 	 = kcrsi(3,3);

  ea 	 = kcrsi(0,0);
  eqy	 = kcrsi(0,2);
  eqz	 = kcrsi(1,0);
  eiyy   = kcrsi(2,2);
  eizz   = kcrsi(1,1);
  eiyz   = kcrsi(1,2);
  dea    = kcrsj(0,0) - kcrsi(0,0);
  deqy   = kcrsj(0,2) - kcrsi(0,2);
  deqz   = kcrsj(1,0) - kcrsi(1,0);
  deiyy  = kcrsj(2,2) - kcrsi(2,2);
  deizz  = kcrsj(1,1) - kcrsi(1,1);
  deiyz  = kcrsj(1,2) - kcrsi(1,2);

  /*****************************************************************************************/
  /* INITIALIZE THE NATURAL ELEMENT STIFFNESS MATRIX AND THE MATRICES USED IN CONDENSATION */
  /*****************************************************************************************/
  kt.Zero();
  ss.Zero();
  sr.Zero();
  temp_sr.Zero();
  rs.Zero();
  rr.Zero();
  temp_rr.Zero();

  /************************************************************/
  /* GENERATE ELASTIC TERMS IN UPPER TRIANGULAR PORTION OF kt */
  /************************************************************/

  kt(0,0)  = ( 7.0 * ea / ( 3.0 * L ) ) + ( 11.0 * dea / ( 6.0 * L ) );
  kt(0,1)  = ( eqz / L ) + ( 2.0 * deqz / ( 3.0 * L ) );
  kt(0,2)  = ( eqy / L ) + ( 2.0 * deqy / ( 3.0 * L ) );
  kt(0,3)  = ( 3.0 * eqz / L ) + ( 7.0 * deqz / ( 3.0 * L ) );
  kt(0,4)  = ( 3.0 * eqy / L ) + ( 7.0 * deqy / ( 3.0 * L ) );
  kt(0,6)  = - ( 8.0 * ea / ( 3.0 * L ) ) - ( 2.0 * dea / L );
  kt(1,1)  = ( 4.0 * eizz / L ) + ( deizz / L );
  kt(1,2)  = ( 4.0 * eiyz / L ) + ( deiyz / L );
  kt(1,3)  = ( 2.0 * eizz / L ) + ( deizz / L );
  kt(1,4)  = ( 2.0 * eiyz / L ) + ( deiyz / L );
  kt(1,6)  = - ( 4.0 * eqz / L ) - ( 4.0 * deqz / ( 3.0 * L ) );
  kt(2,2)  = ( 4.0 * eiyy / L ) + ( deiyy / L );
  kt(2,3)  = ( 2.0 * eiyz / L ) + ( deiyz / L );
  kt(2,4)  = ( 2.0 * eiyy / L ) + ( deiyy / L );
  kt(2,6)  = - ( 4.0 * eqy / L ) - ( 4.0 * deqy / ( 3.0 * L ) );
  kt(3,3)  = ( 4.0 * eizz / L ) + ( 3.0 * deizz / L );
  kt(3,4)  = ( 4.0 * eiyz / L ) + ( 3.0 * deiyz / L );
  kt(3,6)  = - ( 4.0 * eqz / L ) - ( 8.0 * deqz / ( 3.0 * L ) );
  kt(4,4)  = ( 4.0 * eiyy / L ) + ( 3.0 * deiyy / L );
  kt(4,6)  = - ( 4.0 * eqy / L ) - ( 8.0 * deqy / ( 3.0 * L ) );
  kt(5,5)  = gj / L;
  
  kt(6,6)  = ( 16.0 * ea / ( 3.0 * L ) ) + ( 8.0 * dea / ( 3.0 * L ) );

  kt(1,1)  += 2.0 * p * L / 15.0;
  kt(1,3)  += - p * L / 30.0;
  kt(2,2)  += 2.0 * p * L / 15.0;
  kt(2,4)  += - p * L / 30.0;
  kt(3,3)  += 2.0 * p * L / 15.0;
  kt(4,4)  += 2.0 * p * L / 15.0;

  kt(0,0)  += ( 7.0 * p / ( 3.0 * L ) );
  kt(0,1)  += ( mz / L ) + ( 2.0 * dmz / ( 3.0 * L ) );
  kt(0,2)  += ( my / L ) + ( 2.0 * dmy / ( 3.0 * L ) );
  kt(0,3)  += ( 3.0 * mz / L ) + ( 7.0 * dmz / ( 3.0 * L ) );
  kt(0,4)  += ( 3.0 * my / L ) + ( 7.0 * dmy / ( 3.0 * L ) );
  kt(0,6)  -= ( 8.0 * p / ( 3.0 * L ) );
  kt(1,5)  += - dmy / 12.0;
  kt(1,6)  += - ( 4.0 * mz / L ) - ( 4.0 * dmz / ( 3.0 * L ) );
  kt(2,5)  += dmz / 12.0;
  kt(2,6)  -= ( 4.0 * my / L ) + ( 4.0 * dmy / ( 3.0 * L ) );
  kt(3,5)  += dmy / 12.0;
  kt(3,6)  += - ( 4.0 * mz / L ) - ( 8.0 * dmz / ( 3.0 * L ) );
  kt(4,5)  += - dmz / 12.0;
  kt(4,6)  -= ( 4.0 * my / L ) + ( 8.0 * dmy / ( 3.0 * L ) );
  kt(6,6)  += ( 16.0 * p / ( 3.0 * L ) );
  
  /****************************************************/
  /* GENERATE TERMS IN LOWER TRIANGULAR PORTION OF kt */
  /****************************************************/
  for ( i = 0; i < 6; i++ ){
	for ( j = ( i + 1 ); j < 7; j++ ){
		kt(j,i) = kt(i,j);
	}
  }

  /************************************************************************/
  /*									  */
  /* STATIC CONDENSATION PROCEDURE - SEE DON WHITES THESIS		  */
  /*									  */
  /*	| |rr| {rs} |   { D  }   	{F}				  */
  /*	|	    | 		  = 		  		  	  */
  /*	| {sr} |ss| |   { Dm }   	{0}				  */
  /*									  */
  /* so condensed set of equations can be written as:			  */
  /*									  */
  /* |k_condensed| { D } = { F }					  */
  /*									  */
  /* where:								  */
  /*									  */
  /* |k_condensed| = |rr| - {rs} inv|ss| {sr}				  */
  /*									  */
  /************************************************************************/
  double ssdet;

  /*************************************************************/
  /* GENERATE TERMS IN OTHER MATRICIES FOR CONDENSATION	       */
  /*************************************************************/
  ssdet = kt(6,6);

  ss(0,0) = 1.0 / ssdet;

  /*************************************************************/
  /* GENERATE TERMS IN OTHER MATRICIES FOR CONDENSATION        */
  /*************************************************************/

  for ( i = 0; i < 6; i++ ){
         for ( j = 0; j < 1; j++ ){
              rs(i,j) = kt(i,6);
              sr(j,i) = rs(i,j);
         }
         for ( j = 0; j < 6; j++ ){
	      rr(i,j) = kt(i,j);
         }
  }

  /*********************************************/
  /* GENERATE CONDENSED NATURAL ELEMENT MATRIX */
  /*********************************************/

  /* MATRIX MULTIPLICATION */
  for ( i = 0; i < 1; i++ ){
         for ( k = 0; k < 6; k++ ){
                for ( j = 0; j < 1; j++ ){
                        temp_sr(i,k) += ss(i,j) * sr(j,k);
                }
         }
  }

  /* MATRIX MULTIPLICATION */
  for ( i = 0; i < 6; i++ ){
         for ( k = 0; k < 6; k++ ){
                for ( j = 0; j < 1; j++ ){
                        temp_rr(i,k) += rs(i,j) * temp_sr(j,k);
                }
  /* FINAL FORM OF THE CONDENSED NATURAL ELEMENT STIFFNESS MATRIX */
                kv(i,k) = kt(i,k) - temp_rr(i,k);
          }
  }
  
  return 0;
}


const Matrix &
RCFTSTLBeamColumn3D::getMass(void)
{
  theMatrix.Zero();

  return theMatrix;
}



void
RCFTSTLBeamColumn3D::zeroLoad(void)
{
	
}

int
RCFTSTLBeamColumn3D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  return 0;
}


int
RCFTSTLBeamColumn3D::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}

const Vector &
RCFTSTLBeamColumn3D::getResistingForceIncInertia()
{
  // Compute the current resisting force
  theVector = this->getResistingForce();

  // add the damping forces if rayleigh damping
  if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
     theVector += this->getRayleighDampingForces();

  return theVector;
}



bool
RCFTSTLBeamColumn3D::isSubdomain(void)
{
    return false;
}


void
RCFTSTLBeamColumn3D::Print(OPS_Stream &s, int flag)
{
   if (flag == 1)
   {
      s << "\nElement: " << this->getTag() << " Type: RCFTSTLBeamColumn3D ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << numSections;
      s << "\tMass density: " << rho;
      for (int i = 0; i < numSections; i++)
         s << "\nSection "<<i<<" :" << *sections[i];
    }
   else
   {
      s << "\nElement: " << this->getTag() << " Type: RCFTSTLBeamColumn3D ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << numSections;
      s << "\tMass density: " << rho << endln;
   }
}


OPS_Stream &operator<<(OPS_Stream &s, RCFTSTLBeamColumn3D &E)
{
    E.Print(s);
    return s;
}


int
RCFTSTLBeamColumn3D::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // first determine the end points of the beam based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();

    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();

    static Vector v1(3);
    static Vector v2(3);

    for (int i = 0; i < 2; i++) {
	v1(i) = end1Crd(i) + end1Disp(i)*fact;
	v2(i) = end2Crd(i) + end2Disp(i)*fact;
    }

    return theViewer.drawLine (v1, v2, 1.0, 1.0);

}

Response*
RCFTSTLBeamColumn3D::setResponse(char **argv, int argc, Information &eleInformation)
{
    // global force -
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
	|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
	return new ElementResponse(this, 1, theVector);
    else
      return 0;
}

int
RCFTSTLBeamColumn3D::getResponse(int responseID, Information &eleInfo)
{
    switch (responseID) {
      case 1:  // global forces
         return eleInfo.setVector(this->getResistingForce());
      default:
         return -1;
  }
}

int
RCFTSTLBeamColumn3D::setParameter (const char **argv, int argc, Parameter &param)
{
     // From the parameterID value it should be possible to extract
     // information about:
     //  1) Which parameter is in question. The parameter could
     //     be at element, section, or material level.
     //  2) Which section and material number (tag) it belongs to.
     //
     // To accomplish this the parameterID is given the following value:
     //     parameterID = type + 1000*matrTag + 100000*sectionTag
     // ...where 'type' is an integer in the range (1-99) and added 100
     // for each level (from material to section to element).
     //
     // Example:
     //    If 'E0' (case 2) is random in material #3 of section #5
     //    the value of the parameterID at this (element) level would be:
     //    parameterID = 2 + 1000*3 + 100000*5 = 503002
     //    As seen, all given information can be extracted from this number.
     //

     // Initial declarations
     int parameterID = 0;

     // If the parameter belongs to the element itself
     if (strcmp(argv[0],"rho") == 0) 
		 return param.addObject(1,this);

     // If the parameter is belonging to a section or lower
     else if (strcmp(argv[0],"section") == 0) {

	// For now, no parameters of the section itself:
	if (argc<5) {
		cerr << "For now: cannot handle parameters of the section itself." << endl;
		return -1;
	}

	// Get section and material tag numbers from user input
	int paramSectionTag = atoi(argv[1]);

	// Find the right section and call its setParameter method
	for (int i=0; i<numSections; i++) {
		if (paramSectionTag == sections[i]->getTag()) {
			parameterID = sections[i]->setParameter(&argv[2], argc-2, param);
		}
	}

	// Check if the parameterID is valid
	if (parameterID < 0) {
		cerr << "RCFTBeamColumn3D::setParameter() - could not set parameter. " << endl;
		return -1;
	}
	else {
		// Return the parameterID value (according to the above comments)
		return parameterID;
	}
     }

     // Otherwise parameter is unknown for this class
     else {
	return -1;
     }
}

int
RCFTSTLBeamColumn3D::updateParameter (int parameterID, Information &info)
{
     // If the parameterID value is not equal to 1 it belongs
     // to section or material further down in the hierarchy.

     if (parameterID == 1) {
	this->rho = info.theDouble;
	return 0;
     }
     else if (parameterID > 0 ) {
	// Extract the section number
	int sectionNumber = (int)( floor((double)parameterID) / (100000) );

	int ok = -1;
	for (int i=0; i<numSections; i++) {
		if (sectionNumber == sections[i]->getTag()) {
			ok = sections[i]->updateParameter(parameterID, info);
		}
	}

	if (ok < 0) {
		cerr << "RCFTBeamColumn3D::updateParameter() - could not update parameter. " << endl;
		return ok;
	}
	else {
		return ok;
	}
     }
     else {
	cerr << "RCFTBeamColumn3D::updateParameter() - could not update parameter. " << endl;
	return -1;
     }
}

void
RCFTSTLBeamColumn3D::setSectionPointers(int numSec, RCFTSTLFiberSection3D **secPtrs)
{
  if (numSec > maxNumSections) {
    opserr << "Error: RCFTBeamColumn3D::setSectionPointers -- max number of sections exceeded";
  }

  numSections = numSec;

  if (secPtrs == 0) {
    opserr << "Error: RCFTBeamColumn3D::setSectionPointers -- invalid section pointer";
  }

  sections = new RCFTSTLFiberSection3D *[numSections];
  if (sections == 0) {
    opserr << "Error: RCFTBeamColumn3D::setSectionPointers -- could not allocate section pointers";
  }

  for (int i = 0; i < numSections; i++) {

    if (secPtrs[i] == 0) {
      opserr << "Error: RCFTBeamColumn3D::setSectionPointers -- null section pointer " << i << endln;
    }

    sections[i] = (RCFTSTLFiberSection3D*) secPtrs[i]->getCopy();

    if (sections[i] == 0) {
      opserr << "Error: RCFTBeamColumn3D::setSectionPointers -- could not create copy of section " << i << endln;
    }

  }

  ks  = new Matrix [numSections];
  if (ks == 0) {
    opserr << "RCFTBeamColumn3D::setSectionPointers -- failed to allocate ks array";
  }
  
  ksa  = new Matrix [numSections];
  if (ksa == 0) {
    opserr << "RCFTBeamColumn3D::setSectionPointers -- failed to allocate ksa array";
  }

  dhat  = new Vector [numSections];
  if (dhat == 0) {
    opserr << "RCFTBeamColumn3D::setSectionPointers -- failed to allocate dhat array";
  }

  DSQ  = new Vector [numSections];
  if (DSQ == 0) {
    opserr << "RCFTBeamColumn3D::setSectionPointers -- failed to allocate DSQ array";
  }

  DSQa  = new Vector [numSections];
  if (DSQa == 0) {
    opserr << "RCFTBeamColumn3D::setSectionPointers -- failed to allocate DSQa array";
  }
 
  f3  = new Vector [numSections];
  if (f3 == 0) {
    opserr << "RCFTBeamColumn3D::setSectionPointers -- failed to allocate f3 array";
  }

  d3  = new Vector [numSections];
  if (d3 == 0) {
    opserr << "RCFTBeamColumn3D::setSectionPointers -- failed to allocate d3 array";
  }

  str_f3  = new Matrix [numSections];
  if (str_f3 == 0) {
    opserr << "RCFTBeamColumn3D::setSectionPointers -- failed to allocate str_f3 array";
  }

  str_f3inv  = new Matrix [numSections];
  if (str_f3inv == 0) {
    opserr << "RCFTBeamColumn3D::setSectionPointers -- failed to allocate str_f3inv array";
  }
}

Vector
RCFTSTLBeamColumn3D::getLocalIncrDeltaDisp(void)
{

    //ofstream newton;
    //newton.open("newton.dat",ios::app);
	
    const Vector &disp1 = theNodes[0]->getIncrDeltaDisp();
    const Vector &disp2 = theNodes[1]->getIncrDeltaDisp();

    double ug[12];
    for (int i = 0; i < 6; i++) {
        ug[i]   = disp1(i);
        ug[i+6] = disp2(i);
    }

    //newton<<"\n global displ. \n"<<endl;
    
   //for(int i = 0; i < 12; i++) {
   //    newton<<ug[i]<<endl;
   //}
     
    //double L = crdTransf->getInitialLength();
  
    double L = getDeformedLength();
    
    double oneOverL = 1.0/L;

    Vector ul(12);

    //newton<<"\n rotation matrix \n"<<endl;
    //
    //for(int i = 0; i < 3; i++) {
	// newton<<R[i][0]<<"   "<<R[i][1]<<"   "<<R[i][2]<<endl;
    //}

    ul(0)  = R[0][0]*ug[0] + R[0][1]*ug[1] + R[0][2]*ug[2];
    ul(1)  = R[1][0]*ug[0] + R[1][1]*ug[1] + R[1][2]*ug[2];
    ul(2)  = R[2][0]*ug[0] + R[2][1]*ug[1] + R[2][2]*ug[2];

    ul(3)  = R[0][0]*ug[3] + R[0][1]*ug[4] + R[0][2]*ug[5];
    ul(4)  = R[1][0]*ug[3] + R[1][1]*ug[4] + R[1][2]*ug[5];
    ul(5)  = R[2][0]*ug[3] + R[2][1]*ug[4] + R[2][2]*ug[5];

    ul(6)  = R[0][0]*ug[6] + R[0][1]*ug[7] + R[0][2]*ug[8];
    ul(7)  = R[1][0]*ug[6] + R[1][1]*ug[7] + R[1][2]*ug[8];
    ul(8)  = R[2][0]*ug[6] + R[2][1]*ug[7] + R[2][2]*ug[8];

    ul(9)  = R[0][0]*ug[9] + R[0][1]*ug[10] + R[0][2]*ug[11];
    ul(10) = R[1][0]*ug[9] + R[1][1]*ug[10] + R[1][2]*ug[11];
    ul(11) = R[2][0]*ug[9] + R[2][1]*ug[10] + R[2][2]*ug[11];

    //newton<<"\n local displacement \n"<<endl;

    //newton>>ul;

    return ul;
}

const Vector& 
RCFTSTLBeamColumn3D::getBasicIncrDisp(void)
{
    return ub;									    
}	


Vector
RCFTSTLBeamColumn3D::getBasicIncrDeltaDisp(void)
{
   //ofstream output;
   //output.open("localaxes.dat",ios::app);

   //ofstream newton;
   //newton.open("newton101.dat",ios::app);

   //ofstream basic;
   //basic.open("basic.dat",ios::app);

   //ofstream mpls;
   //mpls.open("mpls.dat",ios::app);
	
   const Vector &disp1 = theNodes[0]->getIncrDeltaDisp();
   const Vector &disp2 = theNodes[1]->getIncrDeltaDisp();

   //output<<"\nnode 1 incredeltadisp\n";
   //output>>disp1;

   //output<<"\nnode 2 incredeltadisp\n";
   //output>>disp2;

   const Vector &disp3 = theNodes[0]->getIncrDisp();
   const Vector &disp4 = theNodes[1]->getIncrDisp();

   //output<<"\nnode 1 incrdisp\n";
   //output>>disp3;

   //output<<"\nnode 2 incrdisp\n";
   //output>>disp4;

   const Vector &disp5 = theNodes[0]->getTrialDisp();
   const Vector &disp6 = theNodes[1]->getTrialDisp();

   //output<<"\nnode 1 trialdisp\n";
   //output>>disp5;

   //output<<"\nnode 2 trialdisp\n";
   //output>>disp6;

   double ug[12];
   for (int i = 0; i < 6; i++) {
       ug[i]   = disp1(i);
       ug[i+6] = disp2(i);
   }

   //mpls<<"\n global disp"<<endl;
   //mpls<<ug[0]<<"  "<<ug[1]<<"   "<<ug[2]<<"   "<<ug[3]<<"   "<<ug[4]<<"   "<<ug[5]<<"   "<<ug[6]
   //    <<"  "<<ug[7]<<"   "<<ug[8]<<"   "<<ug[9]<<"   "<<ug[10]<<"   "<<ug[11]<<endl;

   double L = getDeformedLength();

   double oneOverL = 1.0/L;

   Vector dub(7);

   double ul[12];

   //basic<<R[0][0]<<"   "<<R[0][1]<<"   "<<R[0][2]<<endl;
   //basic<<R[1][0]<<"   "<<R[1][1]<<"   "<<R[1][2]<<endl;
   //basic<<R[2][0]<<"   "<<R[2][1]<<"   "<<R[2][2]<<endl;

   ul[0]  = R[0][0]*ug[0] + R[0][1]*ug[1] + R[0][2]*ug[2];
   ul[1]  = R[1][0]*ug[0] + R[1][1]*ug[1] + R[1][2]*ug[2];
   ul[2]  = R[2][0]*ug[0] + R[2][1]*ug[1] + R[2][2]*ug[2];

   ul[3]  = R[0][0]*ug[3] + R[0][1]*ug[4] + R[0][2]*ug[5];
   ul[4]  = R[1][0]*ug[3] + R[1][1]*ug[4] + R[1][2]*ug[5];
   ul[5]  = R[2][0]*ug[3] + R[2][1]*ug[4] + R[2][2]*ug[5];

   ul[6]  = R[0][0]*ug[6] + R[0][1]*ug[7] + R[0][2]*ug[8];
   ul[7]  = R[1][0]*ug[6] + R[1][1]*ug[7] + R[1][2]*ug[8];
   ul[8]  = R[2][0]*ug[6] + R[2][1]*ug[7] + R[2][2]*ug[8];

   ul[9]  = R[0][0]*ug[9] + R[0][1]*ug[10] + R[0][2]*ug[11];
   ul[10] = R[1][0]*ug[9] + R[1][1]*ug[10] + R[1][2]*ug[11];
   ul[11] = R[2][0]*ug[9] + R[2][1]*ug[10] + R[2][2]*ug[11];

   //mpls<<"\n local disp"<<endl;
   //mpls<<ul[0]<<"  "<<ul[1]<<"   "<<ul[2]<<"   "<<ul[3]<<"   "<<ul[4]<<"   "<<ul[5]<<"   "<<ul[6]
   //    <<"  "<<ul[7]<<"   "<<ul[8]<<"   "<<ul[9]<<"   "<<ul[10]<<"   "<<ul[11]<<endl;

   /************************************************************************/
   /* THIS SECTION COMPUTES THE ELONGATION OF THE STEEL (NATURAL DOF 7)    */
   /* THEN THE ELONGATION OF THE CONCRETE (NATURAL DOF 2), THEN THE DIFF   */
   /* IN THE STEEL AND CONCRETE AT END i (NATURAL DOF 1), THEN THE NATURAL */
   /* ROTATION OF EACH END OF THE ELEMENT ENDS.  THESE CALCULATIONS FOLLOW */
   /* MORALES, p. 50, Yang and Kuo, p. 179, Cook, p. 353, etc.             */
   /************************************************************************/

   /************************************************************************/
   /* CALCULATE THE INCREMENTAL ELONGATIONS IN THE X, Y, AND Z GLOBAL SPACE*/
   /************************************************************************/

   double d6_0, d7_1, d8_2;

   d6_0 = ul[6] - ul[0];
   d7_1 = ul[7] - ul[1];
   d8_2 = ul[8] - ul[2];

   /************************************************************************/
   /* CALCULATE THE ELONGATION OF THE STEEL AND THE CONCRETE               */
   /* STEEL                                                                */
   /************************************************************************/
   dub(0) = (  2.0 * L * d6_0 + pow( d6_0 , 2 ) + pow( d7_1 , 2 ) + pow( d8_2 , 2 ) ) /( L + L );
    
   /************************************************************************/
   /* CALCULATE THE RIGID BODY ROTATION OF THE ELEMENT, ASSUMING THAT THE  */
   /* MOVEMENT OF THE STEEL DESCRIBES THE CONCRETE ALSO.                   */
   /************************************************************************/
   double theta_rigid_z, theta_rigid_y;

   theta_rigid_z = ( ul[7] - ul[1] ) / L;

   theta_rigid_y = ( ul[8] - ul[2] ) / L;

   /************************************************************************/
   /* CALCULATE THE NATURAL ROTATIONS OF THE ELEMENT                       */
   /* STEEL ROTATIONS                                                      */
   /************************************************************************/
   dub(1) = ul[5] - theta_rigid_z;
   dub(2) = - ul[4] - theta_rigid_y;
   dub(3) = ul[11] - theta_rigid_z;
   dub(4) = - ul[10] - theta_rigid_y;
   dub(5) = ul[9] - ul[3];

   dub(6) = 0.0;

   //mpls<<"\n $$$$$$$$ sr $$$$$$$$"<<endl;

   //mpls>>sr;

   //mpls>>ss;
   
   for ( int i = 0; i < 6; i++ ){
       dub(6) -= sr(0,i) * dub(i) * ss(0,0);
   }

   return dub;									    
}	

void RCFTSTLBeamColumn3D::calcDeformedLength(void)
{

   //ofstream uiuc;
   //uiuc.open("uiuc.dat",ios::app);

   const Vector &dispi = theNodes[0]->getTrialDisp();
   const Vector &dispj = theNodes[1]->getTrialDisp();

   const Vector &crdi = theNodes[0]->getCrds();
   const Vector &crdj = theNodes[1]->getCrds();

   double ix = crdi(0) + dispi(0);
   double iy = crdi(1) + dispi(1);
   double iz = crdi(2) + dispi(2);
   double jx = crdj(0) + dispj(0);
   double jy = crdj(1) + dispj(1);
   double jz = crdj(2) + dispj(2);
 
   //uiuc<<"\n inside deformed length "<<endl;
   //
   //uiuc>>dispi;
   //uiuc>>dispj;
   //uiuc>>crdi;
   //uiuc>>crdj;

   deflength = sqrt((ix-jx)*(ix-jx)+(iy-jy)*(iy-jy)+(iz-jz)*(iz-jz));
}


double RCFTSTLBeamColumn3D::getDeformedLength(void)
{
   return deflength; 
}

int
RCFTSTLBeamColumn3D::sendSelf(int commitTag, Channel &theChannel)
{  
  return 0;
}
    
int
RCFTSTLBeamColumn3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}






    
									    

