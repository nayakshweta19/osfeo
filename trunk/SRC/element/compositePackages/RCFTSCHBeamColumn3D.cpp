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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/element/RCFTSCHBeamColumn3D.cpp,v $

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <iomanip>

#include <Information.h>
#include <Parameter.h>
#include <RCFTSCHBeamColumn3D.h>
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
#define  NND   9         // number of nodal dof's
#define  NEGD  18        // number of element global dof's
#define  NEBD  12        // number of element dof's in the basic system

//using std::endl;
//using std::ios;
//using std::ifstream;
//using std::ofstream;
using namespace std;

Matrix RCFTSCHBeamColumn3D::theMatrix(18,18);
Vector RCFTSCHBeamColumn3D::theVector(18);
double RCFTSCHBeamColumn3D::workArea[400];

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
// CONSTRUCTOR FOR PARALLEL PROCESSING

RCFTSCHBeamColumn3D::RCFTSCHBeamColumn3D():
Element(0,ELE_TAG_RCFTSCHBeamColumn3D), connectedExternalNodes(2),
beamIntegr(0), numSections(0), sections(0), crdTransf(0),
rho(0.0), initialFlag(0), kv(NEBD,NEBD), Sg(NEGD), Sglobal(NEGD), kvcommit(NEBD,NEBD),
ksa(0), dhat(0), DSQa(0), Ki(0), XAxis(3), YAxis(3), ZAxis(3), 
sr(2,12), ss(2,2), fk(24), ub(12), df_i(18), slp_strn(0), str_f4(0), f4(0), str_f4inv(0), d4(0), Tagg(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;

  Sg.Zero();
  kv.Zero();
  sr.Zero();
  ss.Zero();
  fk.Zero();
  ub.Zero();
  df_i.Zero();
  deflength = 0.0;

  p_si = 0.0;
  p_ci = 0.0;
  my_si = 0.0;
  mz_si = 0.0;
  my_ci = 0.0;
  mz_ci = 0.0;

  p_sj = 0.0;
  p_cj = 0.0;
  my_sj = 0.0;
  mz_sj = 0.0;
  my_cj = 0.0;
  mz_cj = 0.0;
}


// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points.
// allocates the necessary space needed by each object
RCFTSCHBeamColumn3D::RCFTSCHBeamColumn3D (int tag, int nodeI, int nodeJ,
                                      int numSec, RCFTAggregator **sec,
				      BeamIntegration &bi,
                                      RCFTCrdTransf3D &coordTransf, double massDensPerUnitLength):
Element(tag,ELE_TAG_RCFTSCHBeamColumn3D), connectedExternalNodes(2),
beamIntegr(0), numSections(0), sections(0), crdTransf(0),
rho(massDensPerUnitLength),initialFlag(0),kv(NEBD,NEBD), Sg(NEGD), Sglobal(NEGD), kvcommit(NEBD,NEBD),
ksa(0), dhat(0), DSQa(0), Ki(0), XAxis(3), YAxis(3), ZAxis(3),
sr(2,12), ss(2,2), fk(24), ub(12), df_i(18), slp_strn(0.0), str_f4(0), f4(0), str_f4inv(0), d4(0), Tagg(tag)
{

   theNodes[0] = 0;
   theNodes[1] = 0;

   connectedExternalNodes(0) = nodeI;
   connectedExternalNodes(1) = nodeJ;

   // get copy of the sections

   beamIntegr = bi.getCopy();
   if(beamIntegr == 0) {
     opserr<<"Error: RCFTSCHBeamColumn3D::RCFTSCHBeamColumn3D: could not create copy of beam integration object" << endln;
     exit(-1);
   }

   // get copy of the transformation object
   
   crdTransf = (RCFTCrdTransf3D*) coordTransf.getCopy();

   //deflength = crdTransf->getInitialLength();
   
   if (crdTransf == 0) {
      opserr << "Error: RCFTSCHBeamColumn3D::RCFTSCHBeamColumn3D: could not create copy of coordinate transformation object" << endln;
      exit(-1);
   }

   this->setSectionPointers(numSec,sec);

   for(int i = 0; i < numSec; i++){
       dhat[i].Zero();
       ksa[i].Zero();
       DSQa[i].Zero();
       str_f4[i].Zero();
       f4[i].Zero();
       str_f4inv[i].Zero();
       d4[i].Zero();
   }

   //Fill in transformation matrix
   R[0][0]   = 0.0; R[0][1] = 0.0; R[0][2] = 0.0;
   R[1][0]   = 0.0; R[1][1] = 0.0; R[1][2] = 0.0;
   R[2][0]   = 0.0; R[2][1] = 0.0; R[2][2] = 0.0;

   ss.Zero();
   sr.Zero();
   
   kv.Zero();
   fk.Zero();
   ub.Zero();
   df_i.Zero();

   p_si = 0.0;
   p_ci = 0.0;
   my_si = 0.0;
   mz_si = 0.0;
   my_ci = 0.0;
   mz_ci = 0.0;

   p_sj = 0.0;
   p_cj = 0.0;
   my_sj = 0.0;
   mz_sj = 0.0;
   my_cj = 0.0;
   mz_cj = 0.0;
}

// ~RCFT_BeamColumn3D():
// 	destructor
//  delete must be invoked on any objects created by the object
RCFTSCHBeamColumn3D::~RCFTSCHBeamColumn3D()
{

   if (sections) {
      for (int i=0; i < numSections; i++)
         if (sections[i])
            delete sections[i];
      delete [] sections;
   }

   if(dhat != 0)
     delete [] dhat;

   if(ksa != 0)
     delete [] ksa;

   if (crdTransf != 0)
     delete crdTransf;

   if(beamIntegr != 0)
     delete beamIntegr;

   if (Ki != 0)
     delete Ki;

   if(DSQa != 0)
     delete [] DSQa;

   if(str_f4 != 0)
     delete [] str_f4;

   if(f4 != 0)
     delete [] f4;

   if(str_f4inv != 0)
     delete [] str_f4inv;
  
   if(d4 != 0)
     delete [] d4;

}


int
RCFTSCHBeamColumn3D::getNumExternalNodes(void) const
{
   return 2;
}


const ID &
RCFTSCHBeamColumn3D::getExternalNodes(void)
{
   return connectedExternalNodes;
}


Node **
RCFTSCHBeamColumn3D::getNodePtrs(void)
{
   return theNodes;
}

int
RCFTSCHBeamColumn3D::getNumDOF(void)
{
   return NEGD;
}


void
RCFTSCHBeamColumn3D::setDomain(Domain *theDomain)
{
   // check Domain is not null - invoked when object removed from a domain
 
   if (theDomain == 0)
   {
      theNodes[0] = 0;
      theNodes[1] = 0;

      opserr << "RCFTSCHBeamColumn3D::setDomain:  theDomain = 0 ";
      exit(0);
   }

   // get pointers to the nodes

   int Nd1 = connectedExternalNodes(0);
   int Nd2 = connectedExternalNodes(1);

   theNodes[0] = theDomain->getNode(Nd1);
   theNodes[1] = theDomain->getNode(Nd2);

   if (theNodes[0] == 0)
   {
      opserr << "RCFTSCHBeamColumn3D::setDomain: Nd1: ";
      opserr << Nd1 << "does not exist in model\n";
      exit(0);
   }

   if (theNodes[1] == 0)
   {
      opserr << "RCFTSCHBeamColumn3D::setDomain: Nd2: ";
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
      opserr << "RCFTSCHBeamColumn3D::setDomain(): Nd2 or Nd1 incorrect dof ";
      exit(0);
   }

   // initialize the transformation
   if (crdTransf->initialize(theNodes[0], theNodes[1]))
   {
      opserr << "RCFTSCHBeamColumn3D::setDomain(): Error initializing coordinate transformation";
      exit(0);
   }

   deflength = crdTransf->getInitialLength();

   Li = crdTransf->getInitialLength();

   crdTransf->getLocalAxes(XAxis, YAxis, ZAxis);

   R[0][0] = XAxis(0); R[0][1] = XAxis(1); R[0][2] = XAxis(2);
   R[1][0] = YAxis(0); R[1][1] = YAxis(1); R[1][2] = YAxis(2);
   R[2][0] = ZAxis(0); R[2][1] = ZAxis(1); R[2][2] = ZAxis(2);

   // get element length
   double L = crdTransf->getInitialLength();
   if (L == 0.0)
   {
      opserr << "RCFTSCHBeamColumn3D::setDomain(): Zero element length:" << this->getTag();
      exit(0);
   }

   if (initialFlag == 0)
     this->initializeSectionHistoryVariables();

}

int
RCFTSCHBeamColumn3D::commitState()
{
   ofstream dunhat;
   dunhat.open("dunhat.dat",ios::app); 

   ofstream output;
   output.open("stlcon.dat",ios::app);

   ofstream cont1;
   cont1.open("cont1.dat",ios::app);

   ofstream cont2;
   cont2.open("cont2.dat",ios::app);
   
   int err = 0;
   int i = 0;
   int j = 0;

   // call element commitState to do any base class stuff
   if ((err = this->Element::commitState()) != 0) {
     opserr << "RCFTSCHBeamColumn3D::commitState () - failed in base class";
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

   ub.Zero();
   for( i = 0; i < numSections; i++){
        dhat[i].Zero();
   }
   slp_strn = 0.0;

   if( Tagg == 1 ){
    cont1<<DSQa[0](0)<<"  "<<DSQa[1](0)<<"  ";
   }
   if( Tagg == 2 ){
    cont1<<DSQa[0](0)<<"  "<<DSQa[1](0)<<"  ";
   }
   if( Tagg == 3 ){
    cont1<<DSQa[0](0)<<"  "<<DSQa[1](0)<<"  ";
   }
   if( Tagg == 4 ){
    cont1<<DSQa[0](0)<<"  "<<DSQa[1](0)<<"  ";
   }
   if( Tagg == 5 ){
    cont1<<DSQa[0](0)<<"  "<<DSQa[1](0)<<"  ";
   }
   if( Tagg == 6 ){
    cont1<<DSQa[0](0)<<"  "<<DSQa[1](0)<<"  ";
   }
   if( Tagg == 7 ){
    cont1<<DSQa[0](0)<<"  "<<DSQa[1](0)<<"  ";
   }
   if( Tagg == 8 ){
    cont1<<DSQa[0](0)<<"  "<<DSQa[1](0)<<"  ";
   }
   if( Tagg == 9 ){
    cont1<<DSQa[0](0)<<"  "<<DSQa[1](0)<<"  ";
   }
   if( Tagg == 10 ){
    cont1<<DSQa[0](0)<<"  "<<DSQa[1](0)<<endl;
   }

   if( Tagg == 1 ){
    cont2<<DSQa[0](6)<<"  "<<DSQa[1](6)<<"  ";
   }
   if( Tagg == 2 ){
    cont2<<DSQa[0](6)<<"  "<<DSQa[1](6)<<"  ";
   }
   if( Tagg == 3 ){
    cont2<<DSQa[0](6)<<"  "<<DSQa[1](6)<<"  ";
   }
   if( Tagg == 4 ){
    cont2<<DSQa[0](6)<<"  "<<DSQa[1](6)<<"  ";
   }
   if( Tagg == 5 ){
    cont2<<DSQa[0](6)<<"  "<<DSQa[1](6)<<"  ";
   }
   if( Tagg == 6 ){
    cont2<<DSQa[0](6)<<"  "<<DSQa[1](6)<<"  ";
   }
   if( Tagg == 7 ){
    cont2<<DSQa[0](6)<<"  "<<DSQa[1](6)<<"  ";
   }
   if( Tagg == 8 ){
    cont2<<DSQa[0](6)<<"  "<<DSQa[1](6)<<"  ";
   }
   if( Tagg == 9 ){
    cont2<<DSQa[0](6)<<"  "<<DSQa[1](6)<<"  ";
   }
   if( Tagg == 10 ){
    cont2<<DSQa[0](6)<<"  "<<DSQa[1](6)<<endl;
   }

   return err;
}


int RCFTSCHBeamColumn3D::revertToLastCommit()
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
   kv   = kvcommit;

   initialFlag = 0;
   // this->update();

   return err;
}


int RCFTSCHBeamColumn3D::revertToStart()
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
   kvcommit.Zero();

   kv.Zero();

   initialFlag = 0;
   // this->update();
   return err;
}

const Matrix &
RCFTSCHBeamColumn3D::getInitialStiff(void)
{
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
  const Matrix &kcrsj = sections[numSections-1]->getInitialTangent();

  Matrix kvInit(12, 12);

  kvInit.Zero();

  Matrix kt(14,14);
  Matrix temp_sr(2,12);
  Matrix temp_rr(12,12);
  Matrix rs(12,2);
  Matrix rr(12,12);

  /* RIGIDITY TERMS										 */
  double	gj 	= 0.0;	/* TORSIONAL RIGIDITY						 */
  double	kslip 	= 0.0;	/* SLIP LAYER STIFFNESS						 */
  double	dkslip 	= 0.0;	/* SLIP LAYER STIFFNESS						 */

  double	ea_c 	= 0.0;	/* AXIAL RIGIDITY IN CONCRETE					 */
  double	eqy_c	= 0.0;	/* FIRST MOMENT RIGIDITY IN CONCRETE			 	 */
  double	eqz_c	= 0.0;	/* FIRST MOMENT RIGIDITY IN CONCRETE			 	 */
  double	eiyy_c	= 0.0;	/* MOM. OF INERTIA IN CONCRETE					 */
  double	eizz_c	= 0.0;	/* MOM. OF INERTIA IN CONCRETE			 		 */
  double	eiyz_c	= 0.0;	/* MOM. OF INERTIA IN CONCRETE					 */
  double	dea_c 	= 0.0;	/* DIFF. BETWEEN i & j AXIAL RIGIDITY IN CONCRETE		 */
  double	deqy_c	= 0.0;	/* DIFF IN FIRST MOMENT RIGIDITY IN CONCRETE	 		 */
  double	deqz_c	= 0.0;	/* DIFF IN FIRST MOMENT RIGIDITY IN CONCRETE	 		 */
  double	deiyy_c	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN CONCRETE			 	 */
  double	deizz_c	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN CONCRETE			 	 */
  double	deiyz_c	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN CONCRETE			 	 */

  double	ea_s 	= 0.0;	/* AXIAL RIGIDITY IN STEEL					 */
  double	eqy_s	= 0.0;	/* FIRST MOMENT RIGIDITY IN STEEL				 */
  double	eqz_s	= 0.0;	/* FIRST MOMENT RIGIDITY IN STEEL				 */
  double	eiyy_s	= 0.0;	/* MOM. OF INERTIA IN STEEL					 */
  double	eizz_s	= 0.0;	/* MOM. OF INERTIA IN STEEL					 */
  double	eiyz_s	= 0.0;	/* MOM. OF INERTIA IN STEEL					 */
  double	dea_s 	= 0.0;	/* DIFF. BETWEEN i & j AXIAL RIGIDITY IN STEEL	 		 */
  double	deqy_s	= 0.0;	/* DIFF IN FIRST MOMENT RIGIDITY IN STEEL		 	 */
  double	deqz_s	= 0.0;	/* DIFF IN FIRST MOMENT RIGIDITY IN STEEL		 	 */
  double	deiyy_s	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN STEEL			 	 */
  double	deizz_s	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN STEEL			 	 */
  double	deiyz_s	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN STEEL			 	 */

  /************************************************************************/
  /* CALCULATE THE VALUES OF THE LOCAL VARIABLES			  */
  /************************************************************************/
  gj 	= kcrsi(6,6);
  kslip = kcrsi(12,12) * 2 * ( kcrsi(8,8) + kcrsi(7,7) );
  dkslip= kcrsj(12,12) * 2 * ( kcrsj(8,8) + kcrsj(7,7) ) - kcrsi(12,12) * 2 * ( kcrsi(8,8) + kcrsi(7,7) );

  ea_c 	 = kcrsi(0,0);
  eqy_c	 = -kcrsi(0,2);
  eqz_c	 = -kcrsi(1,0);
  eiyy_c = kcrsi(2,2);
  eizz_c = kcrsi(1,1);
  eiyz_c = kcrsi(1,2);
  dea_c  = kcrsj(0,0) - kcrsi(0,0);
  deqy_c = kcrsj(0,2) - kcrsi(0,2);
  deqz_c = kcrsj(1,0) - kcrsi(1,0);
  deiyy_c= kcrsj(2,2) - kcrsi(2,2);
  deizz_c= kcrsj(1,1) - kcrsi(1,1);
  deiyz_c= kcrsj(1,2) - kcrsi(1,2);

  ea_s 	 = kcrsi(3,3);
  eqy_s	 = -kcrsi(3,5);
  eqz_s	 = -kcrsi(3,4);
  eiyy_s = kcrsi(5,5);
  eizz_s = kcrsi(4,4);
  eiyz_s = kcrsi(4,5);
  dea_s  = kcrsj(3,3) - kcrsi(3,3);
  deqy_s = kcrsj(3,5) - kcrsi(3,5);
  deqz_s = kcrsj(3,4) - kcrsi(3,4);
  deiyy_s= kcrsj(5,5) - kcrsi(5,5);
  deizz_s= kcrsj(4,4) - kcrsi(4,4);
  deiyz_s= kcrsj(4,5) - kcrsi(4,5);

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
  /* CONCRETE TERMS */
  kt(1,1)       = 7.0 * ea_c / ( 3.0 * L );
  kt(1,1)       = kt(1,1) + 11.0 * dea_c / ( 6.0 * L );
  kt(1,2) 	= ( eqz_c / L ) + ( 2.0 * deqz_c / ( 3.0 * L ) );
  kt(1,3) 	= ( eqy_c / L ) + ( 2.0 * deqy_c / ( 3.0 * L ) );
  kt(1,4) 	= ( 3.0 * eqz_c / L ) + ( 7.0 * deqz_c / ( 3.0 * L ) );
  kt(1,5) 	= ( 3.0 * eqy_c / L ) + ( 7.0 * deqy_c / ( 3.0 * L ) );
  kt(1,12)	= - ( 8.0 * ea_c / ( 3.0 * L ) ) - ( 2.0 * dea_c / L );
  kt(2,2) 	= ( 4.0 * eizz_c / L ) + ( deizz_c / L );
  kt(2,3) 	= ( 4.0 * eiyz_c / L ) + ( deiyz_c / L );
  kt(2,4) 	= ( 2.0 * eizz_c / L ) + ( deizz_c / L );
  kt(2,5) 	= ( 2.0 * eiyz_c / L ) + ( deiyz_c / L );
  kt(2,12)      = - ( 4.0 * eqz_c / L ) - ( 4.0 * deqz_c / ( 3.0 * L ) );
  kt(3,3) 	= ( 4.0 * eiyy_c / L ) + ( deiyy_c / L );
  kt(3,4) 	= ( 2.0 * eiyz_c / L ) + ( deiyz_c / L );
  kt(3,5) 	= ( 2.0 * eiyy_c / L ) + ( deiyy_c / L );
  kt(3,12)      = - ( 4.0 * eqy_c / L ) - ( 4.0 * deqy_c / ( 3.0 * L ) );
  kt(4,4) 	= ( 4.0 * eizz_c / L ) + ( 3.0 * deizz_c / L );
  kt(4,5) 	= ( 4.0 * eiyz_c / L ) + ( 3.0 * deiyz_c / L );
  kt(4,12)      = - ( 4.0 * eqz_c / L ) - ( 8.0 * deqz_c / ( 3.0 * L ) );
  kt(5,5) 	= ( 4.0 * eiyy_c / L ) + ( 3.0 * deiyy_c / L );
  kt(5,12)      = - ( 4.0 * eqy_c / L ) - ( 8.0 * deqy_c / ( 3.0 * L ) );
  kt(12,12)     = ( 16.0 * ea_c / ( 3.0 * L ) ) + ( 8.0 * dea_c / ( 3.0 * L ) );

  /* STEEL TERMS */
  kt(6,6)	= ( 7.0 * ea_s / ( 3.0 * L ) ) + ( 11.0 * dea_s / ( 6.0 * L ) );
  kt(6,7) 	= ( eqz_s / L ) + ( 2.0 * deqz_s / ( 3.0 * L ) );
  kt(6,8) 	= ( eqy_s / L ) + ( 2.0 * deqy_s / ( 3.0 * L ) );
  kt(6,9)       = ( 3.0 * eqz_s / L ) + ( 7.0 * deqz_s / ( 3.0 * L ) );
  kt(6,10)      = ( 3.0 * eqy_s / L ) + ( 7.0 * deqy_s / ( 3.0 * L ) );
  kt(6,13)	= - ( 8.0 * ea_s / ( 3.0 * L ) ) - ( 2.0 * dea_s / L );
  kt(7,7) 	= ( 4.0 * eizz_s / L ) + ( deizz_s / L );
  kt(7,8) 	= ( 4.0 * eiyz_s / L ) + ( deiyz_s / L );
  kt(7,9)       = ( 2.0 * eizz_s / L ) + ( deizz_s / L );
  kt(7,10)      = ( 2.0 * eiyz_s / L ) + ( deiyz_s / L );
  kt(7,13)      = - ( 4.0 * eqz_s / L ) - ( 4.0 * deqz_s / ( 3.0 * L ) );
  kt(8,8) 	= ( 4.0 * eiyy_s / L ) + ( deiyy_s / L );
  kt(8,9)       = ( 2.0 * eiyz_s / L ) + ( deiyz_s / L );
  kt(8,10)      = ( 2.0 * eiyy_s / L ) + ( deiyy_s / L );
  kt(8,13)      = - ( 4.0 * eqy_s / L ) - ( 4.0 * deqy_s / ( 3.0 * L ) );
  kt(9,9)       = ( 4.0 * eizz_s / L ) + ( 3.0 * deizz_s / L );
  kt(9,10)      = ( 4.0 * eiyz_s / L ) + ( 3.0 * deiyz_s / L );
  kt(9,13)      = - ( 4.0 * eqz_s / L ) - ( 8.0 * deqz_s / ( 3.0 * L ) );
  kt(10,10)     = ( 4.0 * eiyy_s / L ) + ( 3.0 * deiyy_s / L );
  kt(10,13)     = - ( 4.0 * eqy_s / L ) - ( 8.0 * deqy_s / ( 3.0 * L ) );
  kt(11,11)     = gj / L;
  kt(13,13)     = ( 16.0 * ea_s / ( 3.0 * L ) ) + ( 8.0 * dea_s / ( 3.0 * L ) );

  /* SPRING TERMS */
  kt(0,0) 	+= ( L * kslip ) + ( L * dkslip / 2.0 );
  kt(0,1) 	+= ( ( L * kslip ) + ( L * dkslip ) ) / 6.0;
  kt(0,6) 	+= - ( ( L * kslip ) + ( L * dkslip ) ) / 6.0;
  kt(0,12)      += ( ( 2.0 * L * kslip ) + ( L * dkslip ) ) / 3.0;
  kt(0,13)      += - ( ( 2.0 * L * kslip ) + ( L * dkslip ) ) / 3.0;
  kt(1,1) 	+= ( 2.0 * L * kslip / 15.0 ) + ( 7.0 * L * dkslip / 60.0 );
  kt(1,6) 	+= - ( 2.0 * L * kslip / 15.0 ) - ( 7.0 * L * dkslip / 60.0 );
  kt(1,12)      += ( L * kslip / 15.0 ) + ( L * dkslip / 15.0 );
  kt(1,13)      += - ( L * kslip / 15.0 ) - ( L * dkslip / 15.0 );
  kt(6,6) 	+= ( 2.0 * L * kslip / 15.0 ) + ( 7.0 * L * dkslip / 60.0 );
  kt(6,12)      += - ( L * kslip / 15.0 ) - ( L * dkslip / 15.0 );
  kt(6,13)      += ( L * kslip / 15.0 ) + ( L * dkslip / 15.0 );
  kt(12,12)     += ( 8.0 * L * kslip / 15.0 ) + ( 4.0 * L * dkslip / 15.0 );
  kt(12,13)     += - ( 8.0 * L * kslip / 15.0 ) - ( 4.0 * L * dkslip / 15.0 );
  kt(13,13)     += ( 8.0 * L * kslip / 15.0 ) + ( 4.0 * L * dkslip / 15.0 );

  /****************************************************/
  /* GENERATE TERMS IN LOWER TRIANGULAR PORTION OF kt */
  /****************************************************/
  for ( i = 0; i < 13; i++ ){
	for ( j = ( i + 1 ); j < 14; j++ ){
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
  /*************************************************************/
  /* GENERATE TERMS IN TWO_BY_TWO MATRIX FOR CONDENSATION      */
  /*************************************************************/
  double ssdet;
  
  ssdet = kt(12,12) * kt(13,13) - pow( kt(12,13) , 2 );

  ss(0,0) = kt(13,13) / ssdet;
  ss(1,1) = kt(12,12) / ssdet;
  ss(0,1) = - kt(12,13) / ssdet;
  ss(1,0) = ss(0,1);

  /*************************************************************/
  /* GENERATE TERMS IN OTHER MATRICIES FOR CONDENSATION	       */
  /*************************************************************/
  for ( i = 0; i < 12; i++ ){
	for ( j = 0; j < 2; j++ ){
		rs(i,j) = kt(i,12 + j);
		sr(j,i) = rs(i,j);
	}
	for ( j = 0; j < 12; j++ ){
		rr(i,j) =  kt(i,j);
	}
  }

  /*********************************************/
  /* GENERATE CONDENSED NATURAL ELEMENT MATRIX */
  /*********************************************/
  /* MATRIX MULTIPLICATION */
  for ( i = 0; i < 2; i++ ){
	for ( k = 0; k < 12; k++ ){
		for ( j = 0; j < 2; j++ ){
			temp_sr(i,k) += ss(i,j) * sr(j,k);
		}
	}
  }

  /* MATRIX MULTIPLICATION */
  for ( i = 0; i < 12; i++ ){
	for ( k = 0; k < 12; k++ ){
		for ( j = 0; j < 2; j++ ){
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
RCFTSCHBeamColumn3D::getTangentStiff(void)
{
  int i;  
  crdTransf->update();  // Will remove once we clean up the corotational 2d transformation -- MHS
  const Matrix &KV = crdTransf->getGlobalStiffMatrix(kv,fk);
	
  return KV;
}

const Vector &
RCFTSCHBeamColumn3D::getResistingForce(void)
{
  return Sglobal;
}

void RCFTSCHBeamColumn3D::calcResistingForce(void)
{
  crdTransf->update();
  Vector p0(18);
  p0.Zero();
  Sg = crdTransf->getGlobalResistingForce(df_i, p0);
  Sglobal  = Sglobal + Sg;
}

void
RCFTSCHBeamColumn3D::initializeSectionHistoryVariables (void)
{
    for (int i = 0; i < numSections; i++)
    {
        ksa[i]       = Matrix(12,12);
	dhat[i]      = Vector(6);
        DSQa[i]      = Vector(13);
        str_f4[i]    = Matrix(4,4);
        str_f4inv[i] = Matrix(4,4);
        f4[i]        = Vector(4);
        d4[i]        = Vector(4); 
    }
}

/********* NEWTON , SUBDIVIDE AND INITIAL ITERATIONS *********************/
int RCFTSCHBeamColumn3D::update()
{
  //ofstream geom;
  //geom.open("geom.dat",ios::app);
  //ofstream unbal;
  //unbal.open("unbal.dat",ios::app);
  //ofstream mpls;
  //mpls.open("mpls.dat",ios::app);
  //ofstream lstiff;
  //lstiff.open("lstiff.dat",ios::app);
  //ofstream newton19;
  //newton19.open("newton19.dat",ios::app);
  //ofstream mpls1;
  //mpls1.open("mpls1.dat",ios::app);
  //ofstream mpls2;
  //mpls2.open("mpls2.dat",ios::app);
  //ofstream FS;
  //FS.open("FS.dat",ios::app);
  int i,j,k;
  this->calcDeformedLength();
  double L = crdTransf->getInitialLength();
  double oneOverL  = 1.0/L;
  if( initialFlag == 2 )
  this->revertToLastCommit();
  const Vector du = getLocalIncrDeltaDisp();
  //mpls<<"incremental local displacement"<<endl;
  //mpls>>du;
  const Matrix &KL = crdTransf->getLocalStiffMatrix(kv,fk);
  df_i = KL * du; 
  Vector df_nat(12);

  /************************************************************************/
  /* THIS SECTION COMPUTES THE ELONGATION OF THE STEEL (NATURAL DOF 7)    */
  /* THEN THE ELONGATION OF THE CONCRETE (NATURAL DOF 2), THEN THE DIFF   */
  /* IN THE STEEL AND CONCRETE AT END i (NATURAL DOF 1), THEN THE NATURAL */
  /* ROTATION OF EACH END OF THE ELEMENT ENDS.  THESE CALCULATIONS FOLLOW */
  /* MORALES, p. 50, Yang and Kuo, p. 179, Cook, p. 353, etc.             */
  /************************************************************************/

  /* CALCULATE THE INCREMENTAL ELONGATIONS IN THE X, Y, AND Z GLOBAL SPACE*/

  double d9_0, d10_1, d11_2, d15_6, d16_7, d17_8;

  d9_0 = du(9) - du(0);   
  d10_1 = du(10) - du(1);
  d11_2 = du(11) - du(2);
  d15_6 = du(15) - du(6);
  d16_7 = du(16) - du(7);
  d17_8 = du(17) - du(8);

  /************************************************************************/
  /* CALCULATE THE ELONGATION OF THE STEEL AND THE CONCRETE               */
  /* STEEL                                                                */
  /************************************************************************/
  Vector dub(14);

  /* STEEL */
  dub(6) = ( ( ( 2.0 * L + d9_0 ) * d9_0 ) +
          pow( d10_1, 2 ) + pow( d11_2 , 2 ) ) / ( L + L );

  /* CONCRETE */
  dub(1) = ( ( ( 2.0 * L + d15_6 ) * d15_6 ) +
          pow( d16_7, 2 ) + pow( d17_8 , 2 ) ) / ( L + L );

  /************************************************************************/
  /* CALCULATE THE RIGID BODY ROTATION OF THE ELEMENT, ASSUMING THAT THE  */
  /* MOVEMENT OF THE STEEL DESCRIBES THE CONCRETE ALSO.                   */
  /************************************************************************/
  double theta_rigid_z, theta_rigid_y;

  theta_rigid_z = ( du(10) - du(1) ) / L;

  theta_rigid_y = ( du(11) - du(2) ) / L;

  /************************************************************************/
  /* CALCULATE THE NATURAL ROTATIONS OF THE ELEMENT                       */
  /* STEEL ROTATIONS                                                      */
  /************************************************************************/
  dub(2) = du(5) - theta_rigid_z;
  dub(3) = - du(4) - theta_rigid_y;
  dub(4) = du(14) - theta_rigid_z;
  dub(5) = - du(13) - theta_rigid_y;

  /************************************************************************/
  /* CONCRETE ROTATIONS                                                   */
  /************************************************************************/
  dub(7) = dub(2);
  dub(8) = dub(3);
  dub(9) = dub(4);
  dub(10) = dub(5);

  dub(11) = du(12) - du(3);

  /* i-END DIFFERENCE BETWEEN STEEL AND CONCRETE                          */
  dub(0) = du(6) - du(0);

  slp_strn = dub(0);

  /* FIRST CALCULATE THE INCREMENT IN AXIAL FORCE AND ADD THIS TO THE     */
  /* GEOMETRIC STIFFNESS MATRIX BEFORE RECOVERING MOMENTS                 */

  /* SPRING FORCE */
  df_nat(0)    = 0.0;
  for ( int ctr2 = 0; ctr2 < 12; ctr2++ ){
        df_nat(0) += dub(ctr2) * kv(0,ctr2);
  }

  /* CONCRETE FORCE */
  df_nat(1) = 0.0;
  for ( int ctr2 = 0; ctr2 < 12; ctr2++ ){
        df_nat(1) += dub(ctr2) * kv(1,ctr2);
  }

  /* STEEL FORCE */
  df_nat(6)    = 0.0;
  for ( int ctr2 = 0; ctr2 < 12; ctr2++ ){
        df_nat(6) += dub(ctr2) * kv(6,ctr2);
  }

  /* COMPUTE THE AVERAGE STEEL AND CONCRETE AXIAL FORCES TO USE IN Kg     */
  double p_c = df_nat(1) - ( df_nat(0) / 2.0 );
  double pa_c = df_nat(1) - df_nat(0);
  double pb_c = df_nat(1);
  double p_s = df_nat(6) + ( df_nat(0) / 2.0 );
  double pa_s = df_nat(6) + df_nat(0);
  double pb_s = df_nat(6);

  /* UPDATE THE NATURAL GEOMETRIC STIFFNESS MATRIX                        */

  kv(2,2)   += ( 3.0 * pa_c + pb_c ) * ( L / 30.0 );
  kv(2,4)   += - p_c * L / 30.0;
  kv(3,3)   += ( 3.0 * pa_c + pb_c ) * ( L / 30.0 );
  kv(3,5)   += - p_c * L / 30.0;
  kv(4,2)   += - p_c * L / 30.0;
  kv(4,4)   += ( 3.0 * pb_c + pa_c ) * ( L / 30.0 );
  kv(5,3)   += - p_c * L / 30.0;
  kv(5,5)   += ( 3.0 * pb_c + pa_c ) * ( L / 30.0 );
  kv(7,7)   += ( 3.0 * pa_s + pb_s ) * ( L / 30.0 );
  kv(7,9)   += - p_s * L / 30.0;
  kv(8,8)   += ( 3.0 * pa_s + pb_s ) * ( L / 30.0 );
  kv(8,10)  += - p_s * L / 30.0;
  kv(9,7)   += - p_s * L / 30.0;
  kv(9,9)   += ( 3.0 * pb_s + pa_s ) * ( L / 30.0 );
  kv(10,8)  += - p_s * L / 30.0;
  kv(10,10) += ( 3.0 * pb_s + pa_s ) * ( L / 30.0 );

  /* CALCULATE THE NATURAL ELEMENT FORCES FROM THE NATURAL DEFORMATIONS   */
  for ( int ctr1 = 0; ctr1 < 12; ctr1++ ){
        df_nat(ctr1)  = 0.0;
        for ( int ctr2 = 0; ctr2 < 12; ctr2++ ){
                df_nat(ctr1) += dub(ctr2) * kv(ctr1,ctr2);
        }
  }

  //mpls<<"incrmental natural forces"<<endl;
  //mpls>>df_nat;

  /************************************************************************/
  /* COMPUTE THE NATURAL DEFORMATION OF THE MIDPOINT DEGREES OF FREEDOM.  */
  /* THIS IS THE ELONGATION OF THE MEMBER BETWEEN THE MIDPOINT AND THE    */
  /* i-end OF THE MEMBER.  THIS IS ACCOMPLISHED BY UNDOING THE STATIC     */
  /* CONDENSATION THAT WAS DONE TO THE NATURAL STIFFNESS MATRIX IN THE    */
  /* FUNCTION b_cft_tangent_k.c                                           */
  /*                                                                      */
  /* THE NATURAL DOFS ARE:                                                */
  /*      dub(12):  CONCRETE ELONGATION BETWEEN MIDPOINT NODE AND         */
  /*                      i-END OF THE MEMBER                             */
  /*      dub(13):  STEEL ELONGATION BETWEEN MIDPOINT NODE AND            */
  /*                      i-END OF THE MEMBER                             */
  /************************************************************************/
  
  dub(12) = 0.0;
  dub(13) = 0.0;

  double temp_ab[2];

  temp_ab[0] = 0.0;
  temp_ab[1] = 0.0;

  int ctr1;

  for ( ctr1 = 0; ctr1 < 12; ctr1++ ){
         temp_ab[0] +=  sr(0,ctr1) * dub(ctr1);
         temp_ab[1] +=  sr(1,ctr1) * dub(ctr1);
  }
  for ( ctr1 = 0; ctr1 < 2; ctr1++ ){
         dub(12) -= ss(0,ctr1) * temp_ab[ctr1];
         dub(13) -= ss(1,ctr1) * temp_ab[ctr1];
  }

  //mpls<<"incremental natural displacement"<<endl;
  //mpls>>dub;

  df_i(1) = ( df_nat(7) + df_nat(9) ) / L - ( fk(1) + fk(13) ) / 2.0 * theta_rigid_z;
  df_i(2) = ( df_nat(8) + df_nat(10) ) / L - ( fk(1) + fk(13) ) / 2.0 * theta_rigid_y;
  df_i(7) = ( df_nat(2) + df_nat(4) ) / L - ( fk(0) + fk(12) ) / 2.0 * theta_rigid_z;
  df_i(8) = ( df_nat(3) + df_nat(5) ) / L - ( fk(0) + fk(12) ) / 2.0 * theta_rigid_y;
  df_i(0) = - df_nat(6) - df_nat(0);
  df_i(4) = - df_nat(3) - df_nat(8);
  df_i(5) = df_nat(2) + df_nat(7);
  df_i(6) = - df_nat(1) + df_nat(0);
  df_i(9) = df_nat(6);
  df_i(10) = - df_i(1);
  df_i(11) = - df_i(2);
  df_i(13) = - df_nat(5) - df_nat(10);
  df_i(14) = df_nat(4) + df_nat(9);
  df_i(15) = df_nat(1);
  df_i(16) = - df_i(7);
  df_i(17) = - df_i(8);

  //FS<<"new Element"<<endl;
  for ( i = 0; i < numSections; i++){
    //FS<<"cross_section stiffness"<<endl;
    ksa[i] = sections[i]->getSectionTangent();
    //FS>>ksa[i];
    str_f4[i](0,0) = ksa[i](0,0);
    str_f4[i](0,1) = 0.0;
    str_f4[i](0,2) = ksa[i](0,1);
    str_f4[i](0,3) = ksa[i](0,2);

    str_f4[i](1,0) = 0.0;
    str_f4[i](1,1) = ksa[i](3,3);
    str_f4[i](1,2) = ksa[i](3,4);
    str_f4[i](1,3) = ksa[i](3,5);

    str_f4[i](2,0) = ksa[i](0,1);
    str_f4[i](2,1) = ksa[i](3,4);
    str_f4[i](2,2) = ksa[i](4,4) + ksa[i](1,1);
    str_f4[i](2,3) = ksa[i](4,5) + ksa[i](1,2);

    str_f4[i](3,0) = ksa[i](0,2);
    str_f4[i](3,1) = ksa[i](3,5);
    str_f4[i](3,2) = ksa[i](4,5) + ksa[i](1,2);
    str_f4[i](3,3) = ksa[i](5,5) + ksa[i](2,2);
    //mpls>>str_f4[i];
    invertMatrix(3,str_f4[i],str_f4inv[i]);
    if( i == 0 ){
    f4[i](0) = - df_i(6);
    f4[i](1) = - df_i(0);
    f4[i](2) = - df_i(5);
    f4[i](3) =   df_i(4);
    }
    else if( i == 1 ){
    f4[i](0) =  df_i(15);
    f4[i](1) =  df_i(9);
    f4[i](2) =  df_i(14);
    f4[i](3) = - df_i(13);
    }
    //mpls<<"force"<<endl;
    //mpls>>f4[i];
    d4[i] = str_f4inv[i] * f4[i];
    //mpls<<"strain"<<endl;
    //mpls>>d4[i];
    dhat[i](0) = dhat[i](0) + d4[i](0);
    dhat[i](1) = dhat[i](1) + d4[i](2);
    dhat[i](2) = dhat[i](2) + d4[i](3);
    dhat[i](3) = dhat[i](3) + d4[i](1); 
    dhat[i](4) = dhat[i](4) + d4[i](2);
    dhat[i](5) = dhat[i](5) + d4[i](3); 
    Vector aggdhat(13);
    for( int j = 0; j < 6; j++ ){
       aggdhat(j) = dhat[i](j);
    } 
    aggdhat(6) = 0.0;
    aggdhat(7) = 0.0;
    aggdhat(8) = 0.0;
    aggdhat(9) = 0.0;
    aggdhat(10) = 0.0;
    aggdhat(11) = 0.0;
    aggdhat(12) = slp_strn;
    int res = sections[i]->setTrialSectionDeformation(aggdhat);
    DSQa[i] = sections[i]->getStressResultant();    
    //mpls<<"DSQa[i]"<<endl;
    //mpls>>DSQa[i];
    //if( i == 0 ){
    //  mpls1>>DSQa[i];
    //}
    //if( i == 1 ){
    //  mpls2>>DSQa[i];
    //}
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

  //mpls<<"transformation vectors"<<endl;
  //mpls>>XAxis;
  //mpls>>YAxis;
  //mpls>>ZAxis;

  fk(0) = DSQa[0](0);
  fk(1) = DSQa[0](6);
  fk(2) = DSQa[0](2);
  fk(3) = DSQa[0](1);
  fk(4) = DSQa[0](8);
  fk(5) = DSQa[0](7);
  fk(6) = DSQa[0](3);
  fk(7) = DSQa[0](4);
  fk(8) = DSQa[0](5);
  fk(9) = DSQa[0](9);
  fk(10)= DSQa[0](10);
  fk(11)= DSQa[0](11);
  fk(12)= DSQa[1](0);
  fk(13)= DSQa[1](6);
  fk(14)= DSQa[1](2);
  fk(15)= DSQa[1](1);
  fk(16)= DSQa[1](8);
  fk(17)= DSQa[1](7);
  fk(18)= DSQa[1](3);
  fk(19)= DSQa[1](4);
  fk(20)= DSQa[1](5);
  fk(21)= DSQa[1](9);
  fk(22)= DSQa[1](10);
  fk(23)= DSQa[1](11);

  //mpls<<"incremental nodal forces"<<endl;
  //mpls>>df_i;

  /********************************************************/
  /*   GET CROSS-SECTION STIFFNESS AND FORCES FOR THE     */
  /*			  I-END AND J-END	   	  */
  /********************************************************/
  
  Matrix &kcrsi = ksa[0];
  Matrix &kcrsj = ksa[numSections-1];

  Matrix kt(14,14);
  Matrix temp_sr(2,12);
  Matrix temp_rr(12,12);

  /************************************************************************/
  /* TO FACILITATE SIMPLER SYNTAX IN THE EQUATIONS, DEFINE SOME LOCAL	  */
  /* VARIABLES TO REPRESENT RIGIDITY TERMS				  */
  /************************************************************************/
  /* FORCE TERMS										 */
  double	pc_diff	= 0.0;	/* AVERAGE AXIAL FORCE IN CONCRETE		 	 	 */
  double	my_c 	= 0.0;	/* Y_AXIS MOMENT IN CONCRETE					 */
  double	mya_c 	= 0.0;	/* Y_AXIS MOMENT IN CONCRETE					 */
  double	myb_c 	= 0.0;	/* Y_AXIS MOMENT IN CONCRETE					 */
  double	mz_c 	= 0.0;	/* Z_AXIS MOMENT IN CONCRETE					 */
  double	mza_c 	= 0.0;	/* Z_AXIS MOMENT IN CONCRETE					 */
  double	mzb_c 	= 0.0;	/* Z_AXIS MOMENT IN CONCRETE					 */
  double	lzz_c 	= 0.0;	/* ZZ_AXIS SECOND MOMENT OF THE AREA IN CONCRETE 		 */
  double	lyy_c 	= 0.0;	/* YY_AXIS SECOND MOMENT OF THE AREA IN CONCRETE 		 */
  double	lyz_c 	= 0.0;	/* YZ_AXIS SECOND MOMENT OF THE AREA IN CONCRETE 		 */
  double	dmy_c 	= 0.0;	/* Y_AXIS MOMENT IN CONCRETE					 */
  double	dmz_c 	= 0.0;	/* Z_AXIS MOMENT IN CONCRETE					 */
  double	dlzz_c 	= 0.0;	/* ZZ_AXIS SECOND MOMENT OF THE AREA IN CONCRETE 		 */
  double	dlyy_c 	= 0.0;	/* YY_AXIS SECOND MOMENT OF THE AREA IN CONCRETE 		 */
  double	dlyz_c 	= 0.0;	/* YZ_AXIS SECOND MOMENT OF THE AREA IN CONCRETE 		 */

  double	ps_diff	= 0.0;	/* AVERAGE AXIAL FORCE IN STEEL					 */
  double	my_s 	= 0.0;	/* Y_AXIS MOMENT IN STEEL					 */
  double	mya_s 	= 0.0;	/* Y_AXIS MOMENT IN STEEL					 */
  double	myb_s 	= 0.0;	/* Y_AXIS MOMENT IN STEEL					 */
  double	mz_s 	= 0.0;	/* Z_AXIS MOMENT IN STEEL					 */
  double	mza_s 	= 0.0;	/* Z_AXIS MOMENT IN STEEL					 */
  double	mzb_s 	= 0.0;	/* Z_AXIS MOMENT IN STEEL					 */
  double	lzz_s 	= 0.0;	/* ZZ_AXIS SECOND MOMENT OF THE AREA IN STEEL	 		 */
  double	lyy_s 	= 0.0;	/* YY_AXIS SECOND MOMENT OF THE AREA IN STEEL	 	 	 */
  double	lyz_s 	= 0.0;	/* YZ_AXIS SECOND MOMENT OF THE AREA IN STEEL	 		 */
  double	dmy_s 	= 0.0;	/* Y_AXIS MOMENT IN STEEL					 */
  double	dmz_s 	= 0.0;	/* Z_AXIS MOMENT IN STEEL					 */
  double	dlzz_s 	= 0.0;	/* ZZ_AXIS SECOND MOMENT OF THE AREA IN STEEL	 		 */
  double	dlyy_s 	= 0.0;	/* YY_AXIS SECOND MOMENT OF THE AREA IN STEEL	 		 */
  double	dlyz_s 	= 0.0;	/* YZ_AXIS SECOND MOMENT OF THE AREA IN STEEL	 		 */

  /* RIGIDITY TERMS										 */
  double	gj 	= 0.0;	/* TORSIONAL RIGIDITY						 */
  double	kslip 	= 0.0;	/* SLIP LAYER STIFFNESS						 */
  double	dkslip 	= 0.0;	/* SLIP LAYER STIFFNESS						 */

  double	ea_c 	= 0.0;	/* AXIAL RIGIDITY IN CONCRETE					 */
  double	eqy_c	= 0.0;	/* FIRST MOMENT RIGIDITY IN CONCRETE			 	 */
  double	eqz_c	= 0.0;	/* FIRST MOMENT RIGIDITY IN CONCRETE			 	 */
  double	eiyy_c	= 0.0;	/* MOM. OF INERTIA IN CONCRETE					 */
  double	eizz_c	= 0.0;	/* MOM. OF INERTIA IN CONCRETE			 		 */
  double	eiyz_c	= 0.0;	/* MOM. OF INERTIA IN CONCRETE					 */
  double	dea_c 	= 0.0;	/* DIFF. BETWEEN i & j AXIAL RIGIDITY IN CONCRETE		 */
  double	deqy_c	= 0.0;	/* DIFF IN FIRST MOMENT RIGIDITY IN CONCRETE	 		 */
  double	deqz_c	= 0.0;	/* DIFF IN FIRST MOMENT RIGIDITY IN CONCRETE	 		 */
  double	deiyy_c	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN CONCRETE			 	 */
  double	deizz_c	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN CONCRETE			 	 */
  double	deiyz_c	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN CONCRETE			 	 */

  double	ea_s 	= 0.0;	/* AXIAL RIGIDITY IN STEEL					 */
  double	eqy_s	= 0.0;	/* FIRST MOMENT RIGIDITY IN STEEL				 */
  double	eqz_s	= 0.0;	/* FIRST MOMENT RIGIDITY IN STEEL				 */
  double	eiyy_s	= 0.0;	/* MOM. OF INERTIA IN STEEL					 */
  double	eizz_s	= 0.0;	/* MOM. OF INERTIA IN STEEL					 */
  double	eiyz_s	= 0.0;	/* MOM. OF INERTIA IN STEEL					 */
  double	dea_s 	= 0.0;	/* DIFF. BETWEEN i & j AXIAL RIGIDITY IN STEEL	 		 */
  double	deqy_s	= 0.0;	/* DIFF IN FIRST MOMENT RIGIDITY IN STEEL		 	 */
  double	deqz_s	= 0.0;	/* DIFF IN FIRST MOMENT RIGIDITY IN STEEL		 	 */
  double	deiyy_s	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN STEEL			 	 */
  double	deizz_s	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN STEEL			 	 */
  double	deiyz_s	= 0.0;	/* DIFF. IN MOM. OF INERTIA IN STEEL			 	 */

  /************************************************************************/
  /* CALCULATE THE VALUES OF THE LOCAL VARIABLES			  */
  /************************************************************************/
  /* FORCE TERMS							  */
  /************************************************************************/
  /* NOTE THAT THESE TERMS ARE INTERNAL FORCES THAT FOLLOW THE INTERNAL	  */
  /* FORCE POSITIVE SIGN CONVENTION.  THEY DO NOT REPRESENT THE POSITIVE  */
  /* ELEMENT END FORCE SIGN CONVENTION.					  */
  /************************************************************************/
  Matrix rs(12,2);
  Matrix rr(12,12);
  
  p_c 	= ( fk(0) + fk(12) ) / 2.0 ;
  pa_c 	= fk(0);
  pb_c 	= fk(12);
  my_c 	= fk(2);
  mya_c = fk(2);
  myb_c = fk(14);
  mz_c 	= fk(3);
  mza_c = fk(3);
  mzb_c = fk(15);
  lyy_c = fk(6);
  lzz_c = fk(7);
  lyz_c = fk(8);
  dmy_c = myb_c - mya_c;
  dmz_c = mzb_c - mza_c;
  dlyy_c= fk(18) - fk(6);
  dlzz_c= fk(19) - fk(7);
  dlyz_c= fk(20) - fk(8);

  p_s 	= ( fk(1) + fk(13) ) / 2.0 ;
  pa_s 	= fk(1);
  pb_s 	= fk(13);
  my_s 	= fk(4);
  mya_s = fk(4);
  myb_s = fk(16);
  mz_s 	= fk(5);
  mza_s = fk(5);
  mzb_s = fk(17);
  lyy_s = fk(9);
  lzz_s = fk(10);
  lyz_s = fk(11);
  dmy_s = myb_s - mya_s;
  dmz_s = mzb_s - mza_s;
  dlyy_s= fk(21) - fk(9);
  dlzz_s= fk(22) - fk(10);
  dlyz_s= fk(23) - fk(11);

  gj 	= kcrsi(6,6);
  kslip = kcrsi(12,12) * 2 * ( kcrsi(7,7) + kcrsi(8,8) );
  dkslip = kcrsj(12,12) * 2 * ( kcrsj(7,7) + kcrsj(8,8) ) - kcrsi(12,12) * 2 * ( kcrsi(7,7) + kcrsi(8,8) );

  ea_c 	 = kcrsi(0,0);
  eqy_c	 = kcrsi(0,2);
  eqz_c	 = kcrsi(1,0);
  eiyy_c = kcrsi(2,2);
  eizz_c = kcrsi(1,1);
  eiyz_c = kcrsi(1,2);
  dea_c  = kcrsj(0,0) - kcrsi(0,0);
  deqy_c = (kcrsj(0,2) - kcrsi(0,2));
  deqz_c = kcrsj(1,0) - kcrsi(1,0);
  deiyy_c= kcrsj(2,2) - kcrsi(2,2);
  deizz_c= kcrsj(1,1) - kcrsi(1,1);
  deiyz_c= (kcrsj(1,2) - kcrsi(1,2));

  ea_s 	 = kcrsi(3,3);
  eqy_s	 = kcrsi(3,5);
  eqz_s	 = kcrsi(3,4);
  eiyy_s = kcrsi(5,5);
  eizz_s = kcrsi(4,4);
  eiyz_s = kcrsi(4,5);
  dea_s  = kcrsj(3,3) - kcrsi(3,3);
  deqy_s = (kcrsj(3,5) - kcrsi(3,5));
  deqz_s = kcrsj(3,4) - kcrsi(3,4);
  deiyy_s= kcrsj(5,5) - kcrsi(5,5);
  deizz_s= kcrsj(4,4) - kcrsi(4,4);
  deiyz_s= (kcrsj(4,5) - kcrsi(4,5));

  //lstiff<<"concrete"<<endl;
  //
  //lstiff<<"\nea_c "<<ea_c<<" eqy_c "<<eqy_c<<" eqz_c "<<eqz_c<<" eiyy_c "<<eiyy_c<<" eizz_c "
  //      <<eizz_c<<" eiyz_c "<<eiyz_c<<" dea_c "<<dea_c<<" deqy_c "<<deqy_c<<" deqz_c "<<deqz_c<<" deiyy_c "<<deiyy_c<<" deizz_c "<<
  //      deizz_c<<" deiyz_c "<<deiyz_c<<"  "<<kcrsj(1,0)<<"  "<<kcrsi(1,0)<<endl;
  //
  //lstiff<<"steel"<<endl;
  //
  //lstiff<<"\nea_s "<<ea_s<<" eqy_s "<<eqy_s<<" eqz_s "<<eqz_s<<" eiyy_s "<<eiyy_s<<" eizz_s "
  //      <<eizz_s<<" eiyz_s "<<eiyz_s<<" dea_s "<<dea_s<<" deqy_s "<<deqy_s<<" deqz_s "<<deqz_s<<" deiyy_s "<<deiyy_s<<" deizz_s "<<
  //      deizz_s<<" deiyz_s "<<deiyz_s<<endl;



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
  /* CONCRETE TERMS */

  kt(1,1)       = 7.0 * ea_c / ( 3.0 * L );
  kt(1,1)       = kt(1,1) + 11.0 * dea_c / ( 6.0 * L );
  kt(1,2) 	= ( eqz_c / L ) + ( 2.0 * deqz_c / ( 3.0 * L ) );
  kt(1,3) 	= ( eqy_c / L ) + ( 2.0 * deqy_c / ( 3.0 * L ) );
  kt(1,4) 	= ( 3.0 * eqz_c / L ) + ( 7.0 * deqz_c / ( 3.0 * L ) );
  kt(1,5) 	= ( 3.0 * eqy_c / L ) + ( 7.0 * deqy_c / ( 3.0 * L ) );
  kt(1,12)	= - ( 8.0 * ea_c / ( 3.0 * L ) ) - ( 2.0 * dea_c / L );
  kt(2,2) 	= ( 4.0 * eizz_c / L ) + ( deizz_c / L );
  kt(2,3) 	= ( 4.0 * eiyz_c / L ) + ( deiyz_c / L );
  kt(2,4) 	= ( 2.0 * eizz_c / L ) + ( deizz_c / L );
  kt(2,5) 	= ( 2.0 * eiyz_c / L ) + ( deiyz_c / L );
  kt(2,12)      = - ( 4.0 * eqz_c / L ) - ( 4.0 * deqz_c / ( 3.0 * L ) );
  kt(3,3) 	= ( 4.0 * eiyy_c / L ) + ( deiyy_c / L );
  kt(3,4) 	= ( 2.0 * eiyz_c / L ) + ( deiyz_c / L );
  kt(3,5) 	= ( 2.0 * eiyy_c / L ) + ( deiyy_c / L );
  kt(3,12)      = - ( 4.0 * eqy_c / L ) - ( 4.0 * deqy_c / ( 3.0 * L ) );
  kt(4,4) 	= ( 4.0 * eizz_c / L ) + ( 3.0 * deizz_c / L );
  kt(4,5) 	= ( 4.0 * eiyz_c / L ) + ( 3.0 * deiyz_c / L );
  kt(4,12)      = - ( 4.0 * eqz_c / L ) - ( 8.0 * deqz_c / ( 3.0 * L ) );
  kt(5,5) 	= ( 4.0 * eiyy_c / L ) + ( 3.0 * deiyy_c / L );
  kt(5,12)      = - ( 4.0 * eqy_c / L ) - ( 8.0 * deqy_c / ( 3.0 * L ) );
  kt(12,12)     = ( 16.0 * ea_c / ( 3.0 * L ) ) + ( 8.0 * dea_c / ( 3.0 * L ) );

  /* STEEL TERMS */
  kt(6,6)	= ( 7.0 * ea_s / ( 3.0 * L ) ) + ( 11.0 * dea_s / ( 6.0 * L ) );
  kt(6,7) 	= ( eqz_s / L ) + ( 2.0 * deqz_s / ( 3.0 * L ) );
  kt(6,8) 	= ( eqy_s / L ) + ( 2.0 * deqy_s / ( 3.0 * L ) );
  kt(6,9)       = ( 3.0 * eqz_s / L ) + ( 7.0 * deqz_s / ( 3.0 * L ) );
  kt(6,10)      = ( 3.0 * eqy_s / L ) + ( 7.0 * deqy_s / ( 3.0 * L ) );
  kt(6,13)	= - ( 8.0 * ea_s / ( 3.0 * L ) ) - ( 2.0 * dea_s / L );
  kt(7,7) 	= ( 4.0 * eizz_s / L ) + ( deizz_s / L );
  kt(7,8) 	= ( 4.0 * eiyz_s / L ) + ( deiyz_s / L );
  kt(7,9)       = ( 2.0 * eizz_s / L ) + ( deizz_s / L );
  kt(7,10)      = ( 2.0 * eiyz_s / L ) + ( deiyz_s / L );
  kt(7,13)      = - ( 4.0 * eqz_s / L ) - ( 4.0 * deqz_s / ( 3.0 * L ) );
  kt(8,8) 	= ( 4.0 * eiyy_s / L ) + ( deiyy_s / L );
  kt(8,9)       = ( 2.0 * eiyz_s / L ) + ( deiyz_s / L );
  kt(8,10)      = ( 2.0 * eiyy_s / L ) + ( deiyy_s / L );
  kt(8,13)      = - ( 4.0 * eqy_s / L ) - ( 4.0 * deqy_s / ( 3.0 * L ) );
  kt(9,9)       = ( 4.0 * eizz_s / L ) + ( 3.0 * deizz_s / L );
  kt(9,10)      = ( 4.0 * eiyz_s / L ) + ( 3.0 * deiyz_s / L );
  kt(9,13)      = - ( 4.0 * eqz_s / L ) - ( 8.0 * deqz_s / ( 3.0 * L ) );
  kt(10,10)     = ( 4.0 * eiyy_s / L ) + ( 3.0 * deiyy_s / L );
  kt(10,13)     = - ( 4.0 * eqy_s / L ) - ( 8.0 * deqy_s / ( 3.0 * L ) );
  kt(11,11)     = gj / L;
  kt(13,13)     = ( 16.0 * ea_s / ( 3.0 * L ) ) + ( 8.0 * dea_s / ( 3.0 * L ) );

  /* SPRING TERMS */
  kt(0,0) 	+= ( L * kslip ) + ( L * dkslip / 2.0 );
  kt(0,1) 	+= ( ( L * kslip ) + ( L * dkslip ) ) / 6.0;
  kt(0,6) 	+= - ( ( L * kslip ) + ( L * dkslip ) ) / 6.0;
  kt(0,12)      += ( ( 2.0 * L * kslip ) + ( L * dkslip ) ) / 3.0;
  kt(0,13)      += - ( ( 2.0 * L * kslip ) + ( L * dkslip ) ) / 3.0;
  kt(1,1) 	+= ( 2.0 * L * kslip / 15.0 ) + ( 7.0 * L * dkslip / 60.0 );
  kt(1,6) 	+= - ( 2.0 * L * kslip / 15.0 ) - ( 7.0 * L * dkslip / 60.0 );
  kt(1,12)      += ( L * kslip / 15.0 ) + ( L * dkslip / 15.0 );
  kt(1,13)      += - ( L * kslip / 15.0 ) - ( L * dkslip / 15.0 );
  kt(6,6) 	+= ( 2.0 * L * kslip / 15.0 ) + ( 7.0 * L * dkslip / 60.0 );
  kt(6,12)      += - ( L * kslip / 15.0 ) - ( L * dkslip / 15.0 );
  kt(6,13)      += ( L * kslip / 15.0 ) + ( L * dkslip / 15.0 );
  kt(12,12)     += ( 8.0 * L * kslip / 15.0 ) + ( 4.0 * L * dkslip / 15.0 );
  kt(12,13)     += - ( 8.0 * L * kslip / 15.0 ) - ( 4.0 * L * dkslip / 15.0 );
  kt(13,13)     += ( 8.0 * L * kslip / 15.0 ) + ( 4.0 * L * dkslip / 15.0 );

  /**************************************************************/
  /* GENERATE GEOMETRIC TERMS IN UPPER TRIANGULAR PORTION OF kt */
  /**************************************************************/
  kt(2,2)	+= ( 3.0 * pa_c + pb_c ) * ( L / 30.0 );
  kt(2,4)	+= - ( pa_c + pb_c ) * ( L / 60.0 );
  kt(3,3)	+= ( 3.0 * pa_c + pb_c ) * ( L / 30.0 );
  kt(3,5)	+= - ( pa_c + pb_c ) * ( L / 60.0 );
  kt(4,4)	+= ( 3.0 * pb_c + pa_c ) * ( L / 30.0 );
  kt(5,5)	+= ( 3.0 * pb_c + pa_c ) * ( L / 30.0 );
  kt(7,7)	+= ( 3.0 * pa_s + pb_s ) * ( L / 30.0 );
  kt(7,9)	+= - ( pa_s + pb_s ) * ( L / 60.0 );
  kt(8,8)	+= ( 3.0 * pa_s + pb_s ) * ( L / 30.0 );
  kt(8,10)	+= - ( pa_s + pb_s ) * ( L / 30.0 );
  kt(9,9)       += ( 3.0 * pb_s + pa_s ) * ( L / 30.0 );
  kt(10,10)     += ( 3.0 * pb_s + pa_s ) * ( L / 30.0 );

  //mpls<<"pa_c "<<pa_c<<" pb_c "<<pb_c<<" pa_s "<<pa_s<<" pb_s "<<pb_s<<endl;
  
  //mpls<<( 3.0 * pa_c + pb_c ) * ( L / 30.0 )<<"  "<<- ( pa_c + pb_c ) * ( L / 60.0 )<<"   "<<
  //    ( 3.0 * pa_c + pb_c ) * ( L / 30.0 )<<"   "<<- ( pa_c + pb_c ) * ( L / 60.0 )<<"   "<<
  //    ( 3.0 * pb_c + pa_c ) * ( L / 30.0 )<<"   "<<( 3.0 * pb_c + pa_c ) * ( L / 30.0 )<<"   "<<
  //    ( 3.0 * pa_s + pb_s ) * ( L / 30.0 )<<"   "<<- ( pa_s + pb_s ) * ( L / 60.0 )<<"   "<<
  //    ( 3.0 * pa_s + pb_s ) * ( L / 30.0 )<<"   "<<- ( pa_s + pb_s ) * ( L / 30.0 )<<"   "<<
  //    ( 3.0 * pb_s + pa_s ) * ( L / 30.0 )<<"   "<< ( 3.0 * pb_s + pa_s ) * ( L / 30.0 )<<endl; 

  /* HIGHER ORDER GEOMETRIC STIFFNESS TERMS */
  /* CONCRETE TERMS							*/
  kt(1,1)	+= ( 7.0 * p_c / ( 3.0 * L ) );
  kt(1,2) 	+= ( mz_c / L ) + ( 2.0 * dmz_c / ( 3.0 * L ) );
  kt(1,3) 	+= ( my_c / L ) + ( 2.0 * dmy_c / ( 3.0 * L ) );
  kt(1,4) 	+= ( 3.0 * mz_c / L ) + ( 7.0 * dmz_c / ( 3.0 * L ) );
  kt(1,5) 	+= ( 3.0 * my_c / L ) + ( 7.0 * dmy_c / ( 3.0 * L ) );
  kt(1,12)	+= - ( 8.0 * p_c / ( 3.0 * L ) );
  kt(2,2) 	+= ( 4.0 * lzz_c / L ) + ( dlzz_c / L );
  kt(2,3) 	+= ( 4.0 * lyz_c / L ) + ( dlyz_c / L );
  kt(2,4) 	+= ( 2.0 * lzz_c / L ) + ( dlzz_c / L );
  kt(2,5) 	+= ( 2.0 * lyz_c / L ) + ( dlyz_c / L );
  kt(2,12)      += - ( 4.0 * mz_c / L ) - ( 4.0 * dmz_c / ( 3.0 * L ) );
  kt(3,3) 	+= ( 4.0 * lyy_c / L ) + ( dlyy_c / L );
  kt(3,4) 	+= ( 2.0 * lyz_c / L ) + ( dlyz_c / L );
  kt(3,5) 	+= ( 2.0 * lyy_c / L ) + ( dlyy_c / L );
  kt(3,12)      -= ( 4.0 * my_c / L ) + ( 4.0 * dmy_c / ( 3.0 * L ) );
  kt(4,4) 	+= ( 4.0 * lzz_c / L ) + ( 3.0 * dlzz_c / L );
  kt(4,5) 	+= ( 4.0 * lyz_c / L ) + ( 3.0 * dlyz_c / L );
  kt(4,12)      += - ( 4.0 * mz_c / L ) - ( 8.0 * dmz_c / ( 3.0 * L ) );
  kt(5,5) 	+= ( 4.0 * lyy_c / L ) + ( 3.0 * dlyy_c / L );
  kt(5,12)      -= ( 4.0 * my_c / L ) + ( 8.0 * dmy_c / ( 3.0 * L ) );
  kt(12,12)     += ( 16.0 * p_c / ( 3.0 * L ) );

  /* STEEL TERMS							*/
  kt(6,6)	+= ( 7.0 * p_s / ( 3.0 * L ) );
  kt(6,7) 	+= ( mz_s / L ) + ( 2.0 * dmz_s / ( 3.0 * L ) );
  kt(6,8) 	+= ( my_s / L ) + ( 2.0 * dmy_s / ( 3.0 * L ) );
  kt(6,9)       += ( 3.0 * mz_s / L ) + ( 7.0 * dmz_s / ( 3.0 * L ) );
  kt(6,10)      += ( 3.0 * my_s / L ) + ( 7.0 * dmy_s / ( 3.0 * L ) );
  kt(6,13)	+= - ( 8.0 * p_s / ( 3.0 * L ) );
  kt(7,7) 	+= ( 4.0 * lzz_s / L ) + ( dlzz_s / L );
  kt(7,8) 	+= ( 4.0 * lyz_s / L ) + ( dlyz_s / L );
  kt(7,9)       += ( 2.0 * lzz_s / L ) + ( dlzz_s / L );
  kt(7,10)      += ( 2.0 * lyz_s / L ) + ( dlyz_s / L );
  kt(7,13)      += - ( 4.0 * mz_s / L ) - ( 4.0 * dmz_s / ( 3.0 * L ) );
  kt(8,8) 	+= ( 4.0 * lyy_s / L ) + ( dlyy_s / L );
  kt(8,9)       += ( 2.0 * lyz_s / L ) + ( dlyz_s / L );
  kt(8,10)      += ( 2.0 * lyy_s / L ) + ( dlyy_s / L );
  kt(8,13)      -= ( 4.0 * my_s / L ) + ( 4.0 * dmy_s / ( 3.0 * L ) );
  kt(9,9)       += ( 4.0 * lzz_s / L ) + ( 3.0 * dlzz_s / L );
  kt(9,10)      += ( 4.0 * lyz_s / L ) + ( 3.0 * dlyz_s / L );
  kt(9,13)      += - ( 4.0 * mz_s / L ) - ( 8.0 * dmz_s / ( 3.0 * L ) );
  kt(10,10)     += ( 4.0 * lyy_s / L ) + ( 3.0 * dlyy_s / L );
  kt(10,13)     -= ( 4.0 * my_s / L ) + ( 8.0 * dmy_s / ( 3.0 * L ) );
  kt(13,13)     += ( 16.0 * p_s / ( 3.0 * L ) );

  kt(2,11)      += ( mya_c - myb_c ) / 12.0;
  kt(3,11)      += ( mzb_c - mza_c ) / 12.0;
  kt(4,11)      += ( myb_c - mya_c ) / 12.0;
  kt(5,11)      += ( mza_c - mzb_c ) / 12.0;
  kt(7,11)      += ( mya_s - myb_s ) / 12.0;
  kt(8,11)      += ( mzb_s - mza_s ) / 12.0;
  kt(9,11)      += ( myb_s - mya_s ) / 12.0;
  kt(10,11)     += ( mza_s - mzb_s ) / 12.0;
  kt(11,11)     += ( ( lyy_c + lyy_s + lzz_c + lzz_s ) / L ) + ( ( dlyy_c + dlyy_s + dlzz_c + dlzz_s ) / ( 2.0 * L ) ) ;

  /****************************************************/
  /* GENERATE TERMS IN LOWER TRIANGULAR PORTION OF kt */
  /****************************************************/
  for ( i = 0; i < 13; i++ ){
	for ( j = ( i + 1 ); j < 14; j++ ){
		kt(j,i) = kt(i,j);
	}
  }

  //lstiff<<"natural stiffness matrix"<<endl;
  //lstiff>>kt;

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
  /*************************************************************/
  /* GENERATE TERMS IN TWO_BY_TWO MATRIX FOR CONDENSATION      */
  /*************************************************************/
  double ssdet;
  ssdet = kt(12,12) * kt(13,13) - pow( kt(12,13) , 2 );

  ss(0,0) = kt(13,13) / ssdet;
  ss(1,1) = kt(12,12) / ssdet;
  ss(0,1) = - kt(12,13) / ssdet;
  ss(1,0) = ss(0,1);

  /*************************************************************/
  /* GENERATE TERMS IN OTHER MATRICIES FOR CONDENSATION	       */
  /*************************************************************/
  for ( i = 0; i < 12; i++ ){
	for ( j = 0; j < 2; j++ ){
		rs(i,j) = kt(i,12 + j);
		sr(j,i) = rs(i,j);
	}
	for ( j = 0; j < 12; j++ ){
		rr(i,j) =  kt(i,j);
	}
  }

  /*********************************************/
  /* GENERATE CONDENSED NATURAL ELEMENT MATRIX */
  /*********************************************/
  /* MATRIX MULTIPLICATION */
  for ( i = 0; i < 2; i++ ){
	for ( k = 0; k < 12; k++ ){
		for ( j = 0; j < 2; j++ ){
			temp_sr(i,k) += ss(i,j) * sr(j,k);
		}
	}
  }

  /* MATRIX MULTIPLICATION */
  for ( i = 0; i < 12; i++ ){
	for ( k = 0; k < 12; k++ ){
		for ( j = 0; j < 2; j++ ){
			temp_rr(i,k) += rs(i,j) * temp_sr(j,k);
		}
		/* FINAL FORM OF THE CONDENSED NATURAL ELEMENT STIFFNESS MATRIX */
		kv(i,k) = kt(i,k) - temp_rr(i,k);
	}
  }

  Li = L;

  return 0;
}


const Matrix &
RCFTSCHBeamColumn3D::getMass(void)
{
  theMatrix.Zero();

  return theMatrix;
}



void
RCFTSCHBeamColumn3D::zeroLoad(void)
{
	
}

int
RCFTSCHBeamColumn3D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  return 0;
}


int
RCFTSCHBeamColumn3D::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}

const Vector &
RCFTSCHBeamColumn3D::getResistingForceIncInertia()
{
  // Compute the current resisting force
  theVector = this->getResistingForce();

  // add the damping forces if rayleigh damping
  if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
     theVector += this->getRayleighDampingForces();

  return theVector;
}



bool
RCFTSCHBeamColumn3D::isSubdomain(void)
{
    return false;
}


void
RCFTSCHBeamColumn3D::Print(OPS_Stream &s, int flag)
{
   if (flag == 1)
   {
      s << "\nElement: " << this->getTag() << " Type: RCFTSCHBeamColumn3D ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << numSections;
      s << "\tMass density: " << rho;
      for (int i = 0; i < numSections; i++)
         s << "\nSection "<<i<<" :" << *sections[i];
    }
   else
   {
      s << "\nElement: " << this->getTag() << " Type: RCFTSCHBeamColumn3D ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << numSections;
      s << "\tMass density: " << rho << endln;
   }
}


OPS_Stream &operator<<(OPS_Stream &s, RCFTSCHBeamColumn3D &E)
{
    E.Print(s);
    return s;
}


int
RCFTSCHBeamColumn3D::displaySelf(Renderer &theViewer, int displayMode, float fact)
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
RCFTSCHBeamColumn3D::setResponse(char **argv, int argc, Information &eleInformation)
{
    // global force -
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
	|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
	return new ElementResponse(this, 1, theVector);
    else
      return 0;
}

int
RCFTSCHBeamColumn3D::getResponse(int responseID, Information &eleInfo)
{
    switch (responseID) {
      case 1:  // global forces
         return eleInfo.setVector(this->getResistingForce());
      default:
         return -1;
  }
}

int
RCFTSCHBeamColumn3D::setParameter (const char **argv, int argc, Parameter &param)
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
		cerr << "RCFTSCHBeamColumn3D::setParameter() - could not set parameter. " << endl;
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
RCFTSCHBeamColumn3D::updateParameter (int parameterID, Information &info)
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
		cerr << "RCFTSCHBeamColumn3D::updateParameter() - could not update parameter. " << endl;
		return ok;
	}
	else {
		return ok;
	}
     }
     else {
	cerr << "RCFTSCHBeamColumn3D::updateParameter() - could not update parameter. " << endl;
	return -1;
     }
}

void
RCFTSCHBeamColumn3D::setSectionPointers(int numSec, RCFTAggregator **secPtrs)
{
  if (numSec > maxNumSections) {
    opserr << "Error: RCFTSCHBeamColumn3D::setSectionPointers -- max number of sections exceeded";
  }

  numSections = numSec;

  if (secPtrs == 0) {
    opserr << "Error: RCFTSCHBeamColumn3D::setSectionPointers -- invalid section pointer";
  }

  sections = new RCFTAggregator *[numSections];
  if (sections == 0) {
    opserr << "Error: RCFTSCHBeamColumn3D::setSectionPointers -- could not allocate section pointers";
  }

  for (int i = 0; i < numSections; i++) {

    if (secPtrs[i] == 0) {
      opserr << "Error: RCFTSCHBeamColumn3D::setSectionPointers -- null section pointer " << i << endln;
    }

    sections[i] = (RCFTAggregator*) secPtrs[i]->getCopy();

    if (sections[i] == 0) {
      opserr << "Error: RCFTSCHBeamColumn3D::setSectionPointers -- could not create copy of section " << i << endln;
    }

  }

  ksa  = new Matrix [numSections];
  if (ksa == 0) {
    opserr << "RCFTSCHBeamColumn3D::setSectionPointers -- failed to allocate ks array";
  }
	  
  dhat  = new Vector [numSections];
  if (dhat == 0) {
    opserr << "RCFTSCHBeamColumn3D::setSectionPointers -- failed to allocate dhat array";
  }
      
  DSQa  = new Vector [numSections];
  if (DSQa == 0) {
    opserr << "RCFTSCHBeamColumn3D::setSectionPointers -- failed to allocate DSQa array";
  }

  str_f4  = new Matrix [numSections];
  if (str_f4 == 0) {
    opserr << "RCFTSCHBeamColumn3D::setSectionPointers -- failed to allocate str_f4 array";
  }

  str_f4inv = new Matrix [numSections];
  if (str_f4inv == 0) {
    opserr << "RCFTSCHBeamColumn3D::setSectionPointers -- failed to allocate str_f4inv array";
  }

  f4 = new Vector [numSections];
  if (f4 == 0) {
    opserr << "RCFTSCHBeamColumn3D::setSectionPointers -- failed to allocate f4 array";
  }

  d4 = new Vector [numSections];
  if (d4 == 0) {
    opserr << "RCFTSCHBeamColumn3D::setSectionPointers -- failed to allocate d4 array";
  }
}

Vector
RCFTSCHBeamColumn3D::getLocalIncrDeltaDisp(void)
{
    const Vector &disp1 = theNodes[0]->getIncrDeltaDisp();
    const Vector &disp2 = theNodes[1]->getIncrDeltaDisp();

    double ug[18];
    for (int i = 0; i < 9; i++) {
        ug[i]   = disp1(i);
        ug[i+9] = disp2(i);
    }

    double L = crdTransf->getInitialLength();
    
    double oneOverL = 1.0/L;

    Vector ul(18);

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

    ul(12) = R[0][0]*ug[12] + R[0][1]*ug[13] + R[0][2]*ug[14];
    ul(13) = R[1][0]*ug[12] + R[1][1]*ug[13] + R[1][2]*ug[14];
    ul(14) = R[2][0]*ug[12] + R[2][1]*ug[13] + R[2][2]*ug[14];

    ul(15) = R[0][0]*ug[15] + R[0][1]*ug[16] + R[0][2]*ug[17];
    ul(16) = R[1][0]*ug[15] + R[1][1]*ug[16] + R[1][2]*ug[17];
    ul(17) = R[2][0]*ug[15] + R[2][1]*ug[16] + R[2][2]*ug[17];

    return ul;
}

void RCFTSCHBeamColumn3D::calcDeformedLength(void)
{
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

    deflength = sqrt((ix-jx)*(ix-jx)+(iy-jy)*(iy-jy)+(iz-jz)*(iz-jz));
}


double RCFTSCHBeamColumn3D::getDeformedLength(void)
{
   return deflength; 
}

int
RCFTSCHBeamColumn3D::sendSelf(int commitTag, Channel &theChannel)
{  
  return 0;
}
    
int
RCFTSCHBeamColumn3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}
