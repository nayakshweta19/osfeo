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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/element/RCFTSTLDBeamColumn3D.cpp,v $



#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <iomanip>

#include <Information.h>
#include <Parameter.h>
#include <RCFTSTLDBeamColumn3D.h>
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

Matrix RCFTSTLDBeamColumn3D::theMatrix(12,12);
Vector RCFTSTLDBeamColumn3D::theVector(12);
double RCFTSTLDBeamColumn3D::workArea[400];

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
// CONSTRUCTOR FOR PARALLEL PROCESSING

RCFTSTLDBeamColumn3D::RCFTSTLDBeamColumn3D():
Element(0,ELE_TAG_RCFTSTLDBeamColumn3D), connectedExternalNodes(2),
beamIntegr(0), numSections(0), sections(0), crdTransf(0),
rho(0.0), initialFlag(0), kv(NEBD,NEBD), Se(NEBD), Sg(NEGD), df_i(NEGD), Sglobal(NEGD), kvcommit(NEBD,NEBD), Secommit(NEBD),
ks(0), Ki(0), XAxis(3), YAxis(3), ZAxis(3), sr(1,6), ss(1,1), T(5), G(5,6), H(5,5), Hinv(5,5), fint2(6), fnat2(5),fk(12), fk_incr(12), ub(6), Tagg(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;

  T.Zero();
  G.Zero();
  H.Zero();
  Hinv.Zero();
  fint2.Zero();
  fnat2.Zero();
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
RCFTSTLDBeamColumn3D::RCFTSTLDBeamColumn3D (int tag, int nodeI, int nodeJ,
                                      int numSec, RCFTSTLFiberSection3D **sec,
                                      BeamIntegration &bi,
                                      CrdTransf &coordTransf, double massDensPerUnitLength):
Element(tag,ELE_TAG_RCFTSTLDBeamColumn3D), connectedExternalNodes(2),
beamIntegr(0), numSections(0), sections(0), crdTransf(0),
rho(massDensPerUnitLength),
initialFlag(0),
kv(NEBD,NEBD), Se(NEBD), Sg(NEGD), df_i(NEGD), Sglobal(NEGD), kvcommit(NEBD,NEBD), Secommit(NEBD),
ks(0), Ki(0), XAxis(3), YAxis(3), ZAxis(3), sr(1,6), ss(1,1), T(5), G(5,6), H(5,5), Hinv(5,5), fint2(6), fnat2(5), fk(12), fk_incr(12), ub(6), Tagg(tag)
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
   
   crdTransf = coordTransf.getCopy3d();

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

   T.Zero();
   G.Zero();
   H.Zero();
   Hinv.Zero();
   fint2.Zero();
   fnat2.Zero();

   kv.Zero();
   fk.Zero();
   df_i.Zero();
   fk_incr.Zero();
   ub.Zero();
   deflength = 0.0;
}

RCFTSTLDBeamColumn3D::~RCFTSTLDBeamColumn3D()
{

   if (sections) {
      for (int i=0; i < numSections; i++)
         if (sections[i])
            delete sections[i];
      delete [] sections;
   }

   if(ks != 0)
     delete [] ks;

   if(fs != 0)
     delete [] fs;

   if(ksa != 0)
     delete [] ksa;

   if(fsa != 0)
     delete [] fsa;

   if(nldhat != 0)
     delete [] nldhat;

   if(dhat != 0)
     delete [] dhat;

   if(DQ != 0)
     delete [] DQ;

   if(nd1 != 0)
     delete [] nd1;

   if(nd1T != 0)
     delete [] nd1T;

   if(nd1Tf != 0)
     delete [] nd1Tf;

   if(nd1Tfnd1 != 0)
     delete [] nd1Tfnd1;

   if (crdTransf)
     delete crdTransf;

   if(beamIntegr != 0)
     delete beamIntegr;

   if (Ki != 0)
     delete Ki;

}


int
RCFTSTLDBeamColumn3D::getNumExternalNodes(void) const
{
   return 2;
}


const ID &
RCFTSTLDBeamColumn3D::getExternalNodes(void)
{
   return connectedExternalNodes;
}


Node **
RCFTSTLDBeamColumn3D::getNodePtrs(void)
{
   return theNodes;
}

int
RCFTSTLDBeamColumn3D::getNumDOF(void)
{
   return NEGD;
}


void
RCFTSTLDBeamColumn3D::setDomain(Domain *theDomain)
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
RCFTSTLDBeamColumn3D::commitState()
{
#ifdef COMPOSITE_DEBUG
   ofstream dunhat;
   dunhat.open("dunhat.dat",ios::app); 

   ofstream output;
   output.open("stlcon.dat",ios::app);
#endif

   int err = 0;
   int i = 0;
   int j = 0;

   // call element commitState to do any base class stuff
   if ((err = this->Element::commitState()) != 0) {
     opserr << "RCFTBeamColumn3D::commitState () - failed in base class";
     return err;
   }

   do {
#ifdef COMPOSITE_DEBUG
      output<<"section #"<<i<<endl;	   
#endif
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
   fnat2.Zero();
   fint2.Zero();
   for( i = 0; i < numSections; i++){
        dhat[i].Zero();
   }
   return err;
}


int RCFTSTLDBeamColumn3D::revertToLastCommit()
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


int RCFTSTLDBeamColumn3D::revertToStart()
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
RCFTSTLDBeamColumn3D::getInitialStiff(void)
{
#ifdef COMPOSITE_DEBUG
  ofstream newton; 
  newton.open("newton.dat",ios::app);

  ofstream lstiff;
  lstiff.open("lstiff.dat",ios::app);
 
  lstiff<<"\n inside getInitialStiff"<<endl;
#endif

  if (Ki != 0)
    return *Ki;

  int i, j, k;

  double L = crdTransf->getInitialLength();
  double oneOverL  = 1.0/L;

#ifdef COMPOSITE_DEBUG
  lstiff<<L<<endl;
#endif

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
  
  p 	= ( fk(0) + fk(6) ) / 2.0 ;
  pa 	= fk(0);
  pb 	= fk(6);
  my 	= fk(1);
  mya = fk(1);
  myb = fk(7);
  mz 	= fk(2);
  mza = fk(2);
  mzb = fk(8);
  lyy = fk(3);
  lzz = fk(4);
  lyz = fk(5);
  dmy = myb - mya;
  dmz = mzb - mza;
  dlyy= fk(9) - fk(3);
  dlzz= fk(10) - fk(4);
  dlyz= fk(11) - fk(5);

  gj 	 = kcrsi(3,3);

  ea 	 = kcrsi(0,0);
  eqy	 = -kcrsi(0,2);
  eqz	 = -kcrsi(1,0);
  eiyy   = kcrsi(2,2);
  eizz   = kcrsi(1,1);
  eiyz   = kcrsi(1,2);
  dea    = kcrsj(0,0) - kcrsi(0,0);
  deqy   = kcrsj(0,2) - kcrsi(0,2);
  deqz   = kcrsj(1,0) - kcrsi(1,0);
  deiyy  = kcrsj(2,2) - kcrsi(2,2);
  deizz  = kcrsj(1,1) - kcrsi(1,1);
  deiyz  = kcrsj(1,2) - kcrsi(1,2);

#ifdef COMPOSITE_DEBUG
  lstiff<<"\nea "<<ea<<" eqy "<<eqy<<" eqz "<<eqz<<" eiyy "<<eiyy<<" eizz "
        <<eizz<<" eiyz "<<eiyz<<endl;

  lstiff<<"\np "<<p<<" pa "<<pa<<" pb "<<pb<<endl;
#endif

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
  kt(1,1)  += ( 4.0 * lzz / L ) + ( dlzz / L );
  kt(1,2)  += ( 4.0 * lyz / L ) + ( dlyz / L );
  kt(1,3)  += ( 2.0 * lzz / L ) + ( dlzz / L );
  kt(1,4)  += ( 2.0 * lyz / L ) + ( dlyz / L );
  kt(1,5)  += - dmy / 12.0;
  kt(1,6)  += - ( 4.0 * mz / L ) - ( 4.0 * dmz / ( 3.0 * L ) );
  kt(2,2)  += ( 4.0 * lyy / L ) + ( dlyy / L );
  kt(2,3)  += ( 2.0 * lyz / L ) + ( dlyz / L );
  kt(2,4)  += ( 2.0 * lyy / L ) + ( dlyy / L );
  kt(2,5)  += dmz / 12.0;
  kt(2,6)  -= ( 4.0 * my / L ) + ( 4.0 * dmy / ( 3.0 * L ) );
  kt(3,3)  += ( 4.0 * lzz / L ) + ( 3.0 * dlzz / L );
  kt(3,4)  += ( 4.0 * lyz / L ) + ( 3.0 * dlyz / L );
  kt(3,5)  += dmy / 12.0;
  kt(3,6)  += - ( 4.0 * mz / L ) - ( 8.0 * dmz / ( 3.0 * L ) );
  kt(4,4)  += ( 4.0 * lyy / L ) + ( 3.0 * dlyy / L );
  kt(4,5)  += - dmz / 12.0;
  kt(4,6)  -= ( 4.0 * my / L ) + ( 8.0 * dmy / ( 3.0 * L ) );
  kt(5,5)  += ( ( lzz + lyy ) / L ) + ( ( dlzz + dlyy ) / ( 2.0 * L ) );
  kt(6,6)  += ( 16.0 * p / ( 3.0 * L ) );
  
  /****************************************************/
  /* GENERATE TERMS IN LOWER TRIANGULAR PORTION OF kt */
  /****************************************************/
  for ( i = 0; i < 6; i++ ){
	for ( j = ( i + 1 ); j < 7; j++ ){
		kt(j,i) = kt(i,j);
	}
  }

#ifdef COMPOSITE_DEBUG
  lstiff<<"\n natural stiffness"<<endl;

  lstiff>>kv;
#endif

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

#ifdef COMPOSITE_DEBUG
  lstiff<<"\nreduced stiffness matrix"<<endl;
  
  lstiff>>kvInit;
  
  Ki = new Matrix(crdTransf->getInitialGlobalStiffMatrix(kvInit));

  lstiff<<"\n INITIAL GLOBAL STIFFNESS MATRIX"<<Tagg<<endl;
  
  lstiff>>(*Ki);

  lstiff<<"\n THE END \n"<<endl;
#endif

  return *Ki;

}

const Matrix &
RCFTSTLDBeamColumn3D::getTangentStiff(void)
{
//  int i;  
  crdTransf->update();  // Will remove once we clean up the corotational 2d transformation -- MHS
  const Matrix &KV = crdTransf->getGlobalStiffMatrix(kv,fk);
#ifdef COMPOSITE_DEBUG
  ofstream output;
  output.open("tangentstiff.dat", ios::app);
  output<<"\n number of element"<<Tagg<<endl;

  for(i=0; i<12; i++){
    output<<kv(i,0)<<"  "<<kv(i,1)<<"  "<<kv(i,2)<<"  "<<kv(i,3)<<"  "<<kv(i,4)<<"   "
          <<kv(i,5)<<"  "<<kv(i,6)<<"  "<<kv(i,7)<<"  "<<kv(i,8)<<"  "<<kv(i,9)<<"   "
          <<kv(i,10)<<"  "<<kv(i,11)<<endl;
  }
#endif

  return KV;
}

const Vector &
RCFTSTLDBeamColumn3D::getResistingForce(void)
{
  return Sglobal;
}

void RCFTSTLDBeamColumn3D::calcResistingForce(void)
{
#ifdef COMPOSITE_DEBUG
  ofstream intforce;
  intforce.open("intforce.dat",ios::app);
#endif
  crdTransf->update();
  Vector p0(12);
  p0.Zero();
  Sg = crdTransf->getGlobalResistingForce(df_i, p0);
  Sglobal  = Sglobal + Sg;
#ifdef COMPOSITE_DEBUG
  intforce>>Sglobal;
#endif
}

void
RCFTSTLDBeamColumn3D::initializeSectionHistoryVariables (void)
{
    for (int i = 0; i < numSections; i++){
        ksa[i]       = Matrix(4,4);
        fsa[i]       = Matrix(4,4);
        ks[i]        = Matrix(3,3);
        fs[i]        = Matrix(3,3);
        dhat[i]      = Vector(3);
        DQ[i]        = Vector(3);
        nldhat[i]    = Matrix(3,6);
        nd1[i]       = Matrix(3,5);
        nd1T[i]      = Matrix(5,3);
        nd1Tf[i]     = Matrix(5,3);
        nd1Tfnd1[i]  = Matrix(5,5); 
    }
}

/********* NEWTON , SUBDIVIDE AND INITIAL ITERATIONS *********************/
int RCFTSTLDBeamColumn3D::update()
{
#ifdef COMPOSITE_DEBUG
  ofstream geom;
  geom.open("geom.dat",ios::app);

  ofstream FS;
  FS.open("FS.dat",ios::app);

  ofstream unbal;
  unbal.open("unbal.dat",ios::app);

  ofstream mpls;
  mpls.open("mpls.dat",ios::app);

  ofstream lstiff;
  lstiff.open("lstiff.dat",ios::app);
#endif

  int i,j,k;

  this->calcDeformedLength();

  double L = getDeformedLength();

  double oneOverL  = 1.0/L;

  if( initialFlag == 2 )
  this->revertToLastCommit();

  const Vector &dub = getBasicIncrDeltaDisp();
  
  ub = ub + dub;

  const Vector du = getLocalIncrDeltaDisp();

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double wt[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wt);

  /* FIRST CALCULATE THE INCREMENT IN AXIAL FORCE AND ADD THIS TO THE     */
  /* GEOMETRIC STIFFNESS MATRIX BEFORE RECOVERING MOMENTS                 */
  /* ONLY IF GEOMETRICALLY NONLINEAR                                      */
  
#ifdef COMPOSITE_DEBUG
  FS<<"nd1[i]"<<endl;
#endif
  for ( i = 0; i < numSections; i++ ){
     nldhat[i] = this->getNld_hat(i, ub);
     dhat[i] = getd_hat(i, ub);
     nd1[i] = this->getNd1(i, ub);
#ifdef COMPOSITE_DEBUG
	 FS>>nd1[i];
#endif
     for( j = 0; j < 3; j++ ){
         for( k = 0; k < 5; k++ ){
               nd1T[i](k,j) = nd1[i](j,k);
         }
     }
     int res = sections[i]->setTrialSectionDeformation(dhat[i]);
     ksa[i] = sections[i]->getSectionTangent();
     invertMatrix(4,ksa[i],fsa[i]);
     for( j = 0; j < 3; j++ ){
       for( k = 0; k< 3; k++ ){
         fs[i](j,k) = fsa[i](j,k);
         ks[i](j,k) = ksa[i](j,k);
       }
     }
     DQ[i] = ks[i] * dhat[i];
  }

  Vector temp_f(6);
  Matrix nldhatT(6,3);

  double tempfint2[6];

  for( i = 0; i < 6; i++ ){
       tempfint2[i] = fint2(i);
  }

  fint2.Zero();

  for( i = 0; i < numSections; i++ ){
         for( j = 0; j < 6; j++ ){
                for( k = 0; k < 3; k++ ){
                        nldhatT(j,k) = nldhat[i](k,j);
                }
         }
         temp_f = nldhatT * DQ[i];
         for( j = 0; j < 6; j++ ){
                 fint2(j) = fint2(j) + L * wt[i] * temp_f(j);
         }
  }

  Vector dfint(6);

  for( i = 0; i < 6; i++ ){
       dfint(i) = fint2(i) - tempfint2[i];
  }

  df_i(0) = - dfint(0);
  df_i(1) = ( dfint(1) + dfint(3) - fk(0) * (du(7) - du(1)) ) / L;
  df_i(2) = ( dfint(2) + dfint(4) - fk(0) * (du(8) - du(2)) ) / L;
  df_i(3) =  du(3) * ksa[0](3,3) / L - du(9) * ksa[0](3,3) / L;
  df_i(4) = - dfint(2);
  df_i(5) = dfint(1);
  df_i(6) = - df_i(0);
  df_i(7) = - df_i(1);
  df_i(8) = - df_i(2);
  df_i(9) = - du(3) * ksa[1](3,3) / L + du(9) * ksa[1](3,3) / L;
  df_i(10) = - dfint(4);
  df_i(11) = dfint(3);

  fk_incr(0) = - df_i(0);
  fk_incr(1) = df_i(4);
  fk_incr(2) = - df_i(5);
  fk(0)     -= df_i(0);
  fk(1)     += df_i(4);
  fk(2)     -= df_i(5);
  fk(3)      = 0.0;
  fk(4)      = 0.0;
  fk(5)      = 0.0;
  fk_incr(6) = df_i(6);
  fk_incr(7) = - df_i(10);
  fk_incr(8) = df_i(11);
  fk(6)     += df_i(6);
  fk(7)     -= df_i(10);
  fk(8)     += df_i(11);
  fk(9)      = 0.0;
  fk(10)     = 0.0;
  fk(11)     = 0.0;

  //UPDATE TOTAL ELEMENT FORCES WITH RESPECT TO LOCAL COORDINATES

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

  /********************************************************/
  /*   GET CROSS-SECTION STIFFNESS AND FORCES FOR THE     */
  /*			  I-END AND J-END	   	  */
  /********************************************************/
  
  Matrix &kcrsi = ks[0];
  Matrix &kcrsj = ks[1];

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
  
  p 	= ( fk(0) + fk(6) ) / 2.0 ;
  pa 	= fk(0);
  pb 	= fk(6);
  my 	= fk(1);
  mya = fk(1);
  myb = fk(7);
  mz 	= fk(2);
  mza = fk(2);
  mzb = fk(8);
  lyy = fk(3);
  lzz = fk(4);
  lyz = fk(5);
  dmy = myb - mya;
  dmz = mzb - mza;
  dlyy= fk(9) - fk(3);
  dlzz= fk(10) - fk(4);
  dlyz= fk(11) - fk(5);

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
  kt(1,1)  += ( 4.0 * lzz / L ) + ( dlzz / L );
  kt(1,2)  += ( 4.0 * lyz / L ) + ( dlyz / L );
  kt(1,3)  += ( 2.0 * lzz / L ) + ( dlzz / L );
  kt(1,4)  += ( 2.0 * lyz / L ) + ( dlyz / L );
  kt(1,5)  += - dmy / 12.0;
  kt(1,6)  += - ( 4.0 * mz / L ) - ( 4.0 * dmz / ( 3.0 * L ) );
  kt(2,2)  += ( 4.0 * lyy / L ) + ( dlyy / L );
  kt(2,3)  += ( 2.0 * lyz / L ) + ( dlyz / L );
  kt(2,4)  += ( 2.0 * lyy / L ) + ( dlyy / L );
  kt(2,5)  += dmz / 12.0;
  kt(2,6)  -= ( 4.0 * my / L ) + ( 4.0 * dmy / ( 3.0 * L ) );
  kt(3,3)  += ( 4.0 * lzz / L ) + ( 3.0 * dlzz / L );
  kt(3,4)  += ( 4.0 * lyz / L ) + ( 3.0 * dlyz / L );
  kt(3,5)  += dmy / 12.0;
  kt(3,6)  += - ( 4.0 * mz / L ) - ( 8.0 * dmz / ( 3.0 * L ) );
  kt(4,4)  += ( 4.0 * lyy / L ) + ( 3.0 * dlyy / L );
  kt(4,5)  += - dmz / 12.0;
  kt(4,6)  -= ( 4.0 * my / L ) + ( 8.0 * dmy / ( 3.0 * L ) );
  kt(5,5)  += ( ( lzz + lyy ) / L ) + ( ( dlzz + dlyy ) / ( 2.0 * L ) );
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
RCFTSTLDBeamColumn3D::getMass(void)
{
  theMatrix.Zero();

  return theMatrix;
}



void
RCFTSTLDBeamColumn3D::zeroLoad(void)
{
	
}

int
RCFTSTLDBeamColumn3D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  return 0;
}


int
RCFTSTLDBeamColumn3D::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}

const Vector &
RCFTSTLDBeamColumn3D::getResistingForceIncInertia()
{
  // Compute the current resisting force
  theVector = this->getResistingForce();

  // add the damping forces if rayleigh damping
  if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
     theVector += this->getRayleighDampingForces();

  return theVector;
}



bool
RCFTSTLDBeamColumn3D::isSubdomain(void)
{
    return false;
}


void
RCFTSTLDBeamColumn3D::Print(OPS_Stream &s, int flag)
{
   if (flag == 1)
   {
      s << "\nElement: " << this->getTag() << " Type: RCFTSTLDBeamColumn3D ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << numSections;
      s << "\tMass density: " << rho;
      for (int i = 0; i < numSections; i++)
         s << "\nSection "<<i<<" :" << *sections[i];
    }
   else
   {
      s << "\nElement: " << this->getTag() << " Type: RCFTSTLDBeamColumn3D ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << numSections;
      s << "\tMass density: " << rho << endln;
   }
}


OPS_Stream &operator<<(OPS_Stream &s, RCFTSTLDBeamColumn3D &E)
{
    E.Print(s);
    return s;
}


int
RCFTSTLDBeamColumn3D::displaySelf(Renderer &theViewer, int displayMode, float fact)
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
RCFTSTLDBeamColumn3D::setResponse(char **argv, int argc, Information &eleInformation)
{
    // global force -
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
	|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
	return new ElementResponse(this, 1, theVector);
    else
      return 0;
}

int
RCFTSTLDBeamColumn3D::getResponse(int responseID, Information &eleInfo)
{
    switch (responseID) {
      case 1:  // global forces
         return eleInfo.setVector(this->getResistingForce());
      default:
         return -1;
  }
}

int
RCFTSTLDBeamColumn3D::setParameter (const char **argv, int argc, Parameter &param)
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
		cerr << "RCFTSTLDBeamColumn3D::setParameter() - could not set parameter. " << endl;
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
RCFTSTLDBeamColumn3D::updateParameter (int parameterID, Information &info)
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
		cerr << "RCFTSTLDBeamColumn3D::updateParameter() - could not update parameter. " << endl;
		return ok;
	}
	else {
		return ok;
	}
     }
     else {
	cerr << "RCFTSTLDBeamColumn3D::updateParameter() - could not update parameter. " << endl;
	return -1;
     }
}

void
RCFTSTLDBeamColumn3D::setSectionPointers(int numSec, RCFTSTLFiberSection3D **secPtrs)
{
  if (numSec > maxNumSections) {
    opserr << "Error: RCFTSTLDBeamColumn3D::setSectionPointers -- max number of sections exceeded";
  }

  numSections = numSec;

  if (secPtrs == 0) {
    opserr << "Error: RCFTSTLDBeamColumn3D::setSectionPointers -- invalid section pointer";
  }

  sections = new RCFTSTLFiberSection3D *[numSections];
  if (sections == 0) {
    opserr << "Error: RCFTSTLDBeamColumn3D::setSectionPointers -- could not allocate section pointers";
  }

  for (int i = 0; i < numSections; i++) {

    if (secPtrs[i] == 0) {
      opserr << "Error: RCFTSTLDBeamColumn3D::setSectionPointers -- null section pointer " << i << endln;
    }

    sections[i] = (RCFTSTLFiberSection3D*) secPtrs[i]->getCopy();

    if (sections[i] == 0) {
      opserr << "Error: RCFTSTLDBeamColumn3D::setSectionPointers -- could not create copy of section " << i << endln;
    }

  }

  ks  = new Matrix [numSections];
  fs  = new Matrix [numSections];
  ksa = new Matrix [numSections];
  fsa = new Matrix [numSections];
  dhat = new Vector [numSections];
  DQ = new Vector [numSections];
  nd1 = new Matrix [numSections];
  nldhat = new Matrix [numSections];
  nd1T = new Matrix [numSections];
  nd1Tf = new Matrix [numSections];
  nd1Tfnd1 = new Matrix [numSections];

  if (ks == 0) {
    opserr << "RCFTSTLDBeamColumn3D::setSectionPointers -- failed to allocate ks array";
  }

  if (fs == 0) {
    opserr << "RCFTSTLDBeamColumn3D::setSectionPointers -- failed to allocate fs array";
  }

  if (ksa == 0) {
    opserr << "RCFTSTLDBeamColumn3D::setSectionPointers -- failed to allocate ksa array";
  }
  
  if (fsa == 0) {
    opserr << "RCFTSTLDBeamColumn3D::setSectionPointers -- failed to allocate fsa array";
  }
  
  if (nldhat == 0) {
    opserr << "RCFTSTLDBeamColumn3D::setSectionPointers -- failed to allocate nldhat array";
  }
  
  if (dhat == 0) {
    opserr << "RCFTSTLDBeamColumn3D::setSectionPointers -- failed to allocate dhat array";
  }

  if (DQ == 0) {
    opserr << "RCFTSTLDBeamColumn3D::setSectionPointers -- failed to allocate dhat array";
  }

  if (nd1 == 0) {
    opserr << "RCFTSTLDBeamColumn3D::setSectionPointers -- failed to allocate nd1 array";
  }
  
  if (nd1T == 0) {
    opserr <<  "RCFTSTLDBeamColumn3D::setSectionPointers -- failed to allocate nd1T array";
  }

  if (nd1Tf == 0) {
    opserr <<  "RCFTSTLDBeamColumn3D::setSectionPointers -- failed to allocate nd1Tf array";
  }

  if (nd1Tfnd1 == 0) {
    opserr <<  "RCFTSTLDBeamColumn3D::setSectionPointers -- failed to allocate nd1Tfnd1 array";
  }

}

Vector
RCFTSTLDBeamColumn3D::getLocalIncrDeltaDisp(void)
{
#ifdef COMPOSITE_DEBUG
    ofstream newton;
    newton.open("newton.dat",ios::app);
#endif

    const Vector &disp1 = theNodes[0]->getIncrDeltaDisp();
    const Vector &disp2 = theNodes[1]->getIncrDeltaDisp();

    double ug[12];
    for (int i = 0; i < 6; i++) {
        ug[i]   = disp1(i);
        ug[i+6] = disp2(i);
    }

#ifdef COMPOSITE_DEBUG
	newton<<"\n global displ. \n"<<endl;
    
    for(int i = 0; i < 12; i++) {
        newton<<ug[i]<<endl;
    }
#endif

    //double L = crdTransf->getInitialLength();
  
    double L = getDeformedLength();
    
    double oneOverL = 1.0/L;

    Vector ul(12);

#ifdef COMPOSITE_DEBUG
	newton<<"\n rotation matrix \n"<<endl;
    
    for(int i = 0; i < 3; i++) {
	 newton<<R[i][0]<<"   "<<R[i][1]<<"   "<<R[i][2]<<endl;
    }
#endif

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

#ifdef COMPOSITE_DEBUG
	newton<<"\n local displacement \n"<<endl;
	
	newton>>ul;
#endif

    return ul;
}

const Vector& 
RCFTSTLDBeamColumn3D::getBasicIncrDisp(void)
{
    return ub;									    
}	


Vector
RCFTSTLDBeamColumn3D::getBasicIncrDeltaDisp(void)
{
#ifdef COMPOSITE_DEBUG
   ofstream output;
   output.open("localaxes.dat",ios::app);
   
   ofstream newton;
   newton.open("newton101.dat",ios::app);
   
   ofstream basic;
   basic.open("basic.dat",ios::app);
   
   ofstream mpls;
   mpls.open("mpls.dat",ios::app);
#endif

   const Vector &disp1 = theNodes[0]->getIncrDeltaDisp();
   const Vector &disp2 = theNodes[1]->getIncrDeltaDisp();

#ifdef COMPOSITE_DEBUG
   output<<"\nnode 1 incredeltadisp\n";
   output>>disp1;

   output<<"\nnode 2 incredeltadisp\n";
   output>>disp2;
#endif

   const Vector &disp3 = theNodes[0]->getIncrDisp();
   const Vector &disp4 = theNodes[1]->getIncrDisp();

#ifdef COMPOSITE_DEBUG
   output<<"\nnode 1 incrdisp\n";
   output>>disp3;

   output<<"\nnode 2 incrdisp\n";
   output>>disp4;
#endif

   const Vector &disp5 = theNodes[0]->getTrialDisp();
   const Vector &disp6 = theNodes[1]->getTrialDisp();

#ifdef COMPOSITE_DEBUG
   output<<"\nnode 1 trialdisp\n";
   output>>disp5;

   output<<"\nnode 2 trialdisp\n";
   output>>disp6;
#endif

   double ug[12];
   for (int i = 0; i < 6; i++) {
       ug[i]   = disp1(i);
       ug[i+6] = disp2(i);
   }

#ifdef COMPOSITE_DEBUG
   mpls<<"\n global disp"<<endl;
   mpls<<ug[0]<<"  "<<ug[1]<<"   "<<ug[2]<<"   "<<ug[3]<<"   "<<ug[4]<<"   "<<ug[5]<<"   "<<ug[6]
       <<"  "<<ug[7]<<"   "<<ug[8]<<"   "<<ug[9]<<"   "<<ug[10]<<"   "<<ug[11]<<endl;
#endif

   double L = getDeformedLength();

   double oneOverL = 1.0/L;

   Vector dub(7);

   double ul[12];

#ifdef COMPOSITE_DEBUG
   basic<<R[0][0]<<"   "<<R[0][1]<<"   "<<R[0][2]<<endl;
   basic<<R[1][0]<<"   "<<R[1][1]<<"   "<<R[1][2]<<endl;
   basic<<R[2][0]<<"   "<<R[2][1]<<"   "<<R[2][2]<<endl;
#endif

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

#ifdef COMPOSITE_DEBUG
   mpls<<"\n local disp"<<endl;
   mpls<<ul[0]<<"  "<<ul[1]<<"   "<<ul[2]<<"   "<<ul[3]<<"   "<<ul[4]<<"   "<<ul[5]<<"   "<<ul[6]
       <<"  "<<ul[7]<<"   "<<ul[8]<<"   "<<ul[9]<<"   "<<ul[10]<<"   "<<ul[11]<<endl;
#endif

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

#ifdef COMPOSITE_DEBUG
   mpls<<"\n $$$$$$$$ sr $$$$$$$$"<<endl;

   mpls>>sr;

   mpls>>ss;
#endif

   for ( int i = 0; i < 6; i++ ){
       dub(6) -= sr(0,i) * dub(i) * ss(0,0);
   }

   Vector dub2(6);

   dub2(0) = dub(0);
   dub2(1) = dub(1);
   dub2(2) = dub(2);
   dub2(3) = dub(3);
   dub2(4) = dub(4);
   dub2(5) = dub(6);

   return dub2;									    
}

Vector
RCFTSTLDBeamColumn3D::getd_hat(int sec, const Vector &v)
{
#ifdef COMPOSITE_DEBUG
   ofstream newton;
   newton.open("dhat.dat",ios::app);
#endif

   //double L = crdTransf->getInitialLength();

   double L = getDeformedLength();
   double oneOverL  = 1.0/L;

   double xi[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double wt[maxNumSections];
   beamIntegr->getSectionWeights(numSections, L, wt);

   double temp_x, temp_A, temp_B, temp_C, temp_D, temp_E, temp_F;

   Vector D_hat(3);

   temp_x = L * xi[sec];
   temp_A = 1 - 4 * ( temp_x / L ) + 3 * pow ( ( temp_x / L ) , 2 );
   temp_B = - 2 * ( temp_x / L ) + 3 * pow ( ( temp_x / L ) , 2 );
   temp_C = - 1 / L + 4 * temp_x / ( L * L );
   temp_D = - 8 * temp_x / ( L * L ) + 4 / L;
   temp_E = - 4 / L + 6 * temp_x / ( L * L );
   temp_F =  - 2 / L + 6 * temp_x / ( L * L );

#ifdef COMPOSITE_DEBUG
   newton<<"\n INSIDE GETD_HAT "<<sec<<endl;

   newton>>v;

   newton<<"\n intermediate variables"<<endl;

   newton<<" temp_A "<<temp_A<<"  "<<"temp_B "<<temp_B<<"   "<<"temp_C "<<temp_C
	   <<"   "<<"temp_D "<<temp_D<<"   "<<"temp_E "<<temp_E<<endl;
#endif

   //D_hat(0) =  temp_C * v(0) +
   //            0.5 * ( temp_C * temp_C * v(0) + temp_C * temp_D * v(6) ) * v(0) +
   //            0.5 * ( temp_A * temp_A * v(1) + temp_A * temp_B * v(3) ) * v(1) +
   //            0.5 * ( temp_A * temp_A * v(2) + temp_A * temp_B * v(4) ) * v(2) +
   //            0.5 * ( temp_A * temp_B * v(1) + temp_B * temp_B * v(3) ) * v(3) +
   //            0.5 * ( temp_A * temp_B * v(2) + temp_B * temp_B * v(4) ) * v(4) +
   //            0.5 * ( temp_D * temp_D * v(6) + temp_C * temp_D * v(0) ) * v(6) +
   //            temp_D * v(6);

   D_hat(0) =  temp_C*v(0) + temp_D*v(5) +
               v(1)*v(1)/15 - v(3)*v(1)/30 + v(3)*v(3)/15 +
               v(2)*v(2)/15 - v(4)*v(2)/30 + v(4)*v(4)/15;

   D_hat(1) =  temp_E * v(1) + temp_F * v(3);

   D_hat(2) =  temp_E * v(2) + temp_F * v(4);

#ifdef COMPOSITE_DEBUG
   newton<<"\n D_hat values"<<endl;

   newton>>D_hat;
#endif

   return D_hat;
}

Matrix
RCFTSTLDBeamColumn3D::getNld_hat(int sec, const Vector &v)
{
#ifdef COMPOSITE_DEBUG
   ofstream newton;
   newton.open("mpls.dat",ios::app);
#endif

   int i,j;

   double L = getDeformedLength();

   double oneOverL  = 1.0/L;

   double xi[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double wt[maxNumSections];
   beamIntegr->getSectionWeights(numSections, L, wt);

   double temp_x, temp_A, temp_B, temp_C, temp_D, temp_E, temp_F;

   temp_x = L * xi[sec];
   temp_A = 1 - 4 * ( temp_x / L ) + 3 * pow ( ( temp_x / L ) , 2 );
   temp_B = - 2 * ( temp_x / L ) + 3 * pow ( ( temp_x / L ) , 2 );
   temp_C = - 1 / L + 4 * temp_x / ( L * L );
   temp_D = - 8 * temp_x / ( L * L ) + 4 / L;
   temp_E = - 4 / L + 6 * temp_x / ( L * L );
   temp_F =  - 2 / L + 6 * temp_x / ( L * L );

   Matrix Nld_hat(3,6);

   for(i = 0; i < 3; i++){
       for(j = 0; j < 6; j++){
               Nld_hat(i,j) = 0.0;
       }
   }

   Nld_hat(0,0) = temp_C + temp_C * temp_C * v(0) + temp_C * temp_D * v(5);
   Nld_hat(0,1) = ( temp_A * temp_A * v(1) + temp_A * temp_B * v(3) );
   Nld_hat(0,2) = ( temp_A * temp_A * v(2) + temp_A * temp_B * v(4) );
   Nld_hat(0,3) = ( temp_A * temp_B * v(1) + temp_B * temp_B * v(3) );
   Nld_hat(0,4) = ( temp_A * temp_B * v(2) + temp_B * temp_B * v(4) );
   Nld_hat(0,5) = temp_D + temp_D * temp_D * v(5) + temp_C * temp_D * v(0);
   Nld_hat(1,1) = temp_E;
   Nld_hat(1,3) = temp_F;
   Nld_hat(2,2) = temp_E;
   Nld_hat(2,4) = temp_F;

   return Nld_hat;
}

Matrix
RCFTSTLDBeamColumn3D::getNd1(int sec, const Vector &v)
{
#ifdef COMPOSITE_DEBUG
	ofstream newton;
    newton.open("newton.dat",ios::app);
#endif

    double L = getDeformedLength();
    double oneOverL  = 1.0/L;

    double xi[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, xi);

    double wt[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wt);

    double temp_x, temp_A, temp_B;

    temp_x = L * xi[sec];

    temp_A = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[1]
            + L * ( - pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[3];

    temp_B = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[2]
            + L * ( - pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[4];

    Matrix Nd1(3,5);

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 5; j++){
            Nd1(i,j) = 0.0;
        }
    }

    Nd1(0,0)   = 1.0;
    Nd1(1,0)   = temp_A;
    Nd1(1,1)   = temp_x / L - 1.0;
    Nd1(1,3)   = temp_x / L;
    Nd1(2,0)   = temp_B;
    Nd1(2,2)   = temp_x / L - 1.0;
    Nd1(2,4)   = temp_x / L;

    return Nd1;
}
	
void RCFTSTLDBeamColumn3D::calcDeformedLength(void)
{
#ifdef COMPOSITE_DEBUG
   ofstream uiuc;
   uiuc.open("uiuc.dat",ios::app);
#endif

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
 
#ifdef COMPOSITE_DEBUG
   uiuc<<"\n inside deformed length "<<endl;
   
   uiuc>>dispi;
   uiuc>>dispj;
   uiuc>>crdi;
   uiuc>>crdj;
#endif

   deflength = sqrt((ix-jx)*(ix-jx)+(iy-jy)*(iy-jy)+(iz-jz)*(iz-jz));
}


double RCFTSTLDBeamColumn3D::getDeformedLength(void)
{
   return deflength; 
}

int
RCFTSTLDBeamColumn3D::sendSelf(int commitTag, Channel &theChannel)
{  
  return 0;
}
    
int
RCFTSTLDBeamColumn3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}






    
									    

