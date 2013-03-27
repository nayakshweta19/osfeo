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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/element/RCFTSTLGMBeamColumn3D.cpp,v $



#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <iomanip>

#include <Information.h>
#include <Parameter.h>
#include <RCFTSTLGMBeamColumn3D.h>
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
//#include <iomanip.h>

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

Matrix RCFTSTLGMBeamColumn3D::theMatrix(12,12);
Vector RCFTSTLGMBeamColumn3D::theVector(12);
double RCFTSTLGMBeamColumn3D::workArea[400];

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
// CONSTRUCTOR FOR PARALLEL PROCESSING

RCFTSTLGMBeamColumn3D::RCFTSTLGMBeamColumn3D():
Element(0,ELE_TAG_RCFTSTLGMBeamColumn3D), connectedExternalNodes(2),
beamIntegr(0), numSections(0), sections(0), crdTransf(0),
rho(0.0), initialFlag(0), kv(NEBD,NEBD), Se(NEBD), Sg(NEGD), df_i(NEGD), f_i(NEGD), Cf_i(NEGD), Sglobal(NEGD), CSglobal(NEGD), kvcommit(NEBD,NEBD), Secommit(NEBD),
ks(0), Ki(0), XAxis(3), YAxis(3), ZAxis(3), sr(1,6), ss(1,1), T(5), V(5), V2(5), G(5,5), GT(5,5), G2(5,5), G2T(5,5), H(5,5), H12(5,5), H22(5,5), Md(5,5), Hinv(5,5),Kg(5,5), fint2(5), fnat2(5), Tfnat2(5), fk(5), Cfk(5), fk_incr(5), ub(5), Tagg(0), itr(0), cnvg(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;

  Sg.Zero();
  T.Zero();
  G.Zero();
  GT.Zero();
  G2.Zero();
  G2T.Zero();
  H.Zero();
  H12.Zero();
  H22.Zero();
  Md.Zero();
  Kg.Zero();
  Hinv.Zero();
  fint2.Zero();
  fnat2.Zero();
  Tfnat2.Zero();
  kv.Zero();
  sr.Zero();
  ss.Zero();
  fk.Zero();
  Cfk.Zero();
  fk_incr.Zero();
  ub.Zero();
  df_i.Zero();
  f_i.Zero();
  Cf_i.Zero();
  V.Zero();
  V2.Zero();
  deflength = 0.0;

}


// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points.
// allocates the necessary space needed by each object
RCFTSTLGMBeamColumn3D::RCFTSTLGMBeamColumn3D (int tag, int nodeI, int nodeJ,
                                      int numSec, RCFTSTLFiberSection3D **sec,
                                      BeamIntegration &bi,
                                      CrdTransf &coordTransf, double massDensPerUnitLength):
Element(tag,ELE_TAG_RCFTSTLGMBeamColumn3D), connectedExternalNodes(2),
beamIntegr(0), numSections(0), sections(0), crdTransf(0),
rho(massDensPerUnitLength),
initialFlag(0),
kv(NEBD,NEBD), Se(NEBD), Sg(NEGD), df_i(NEGD), f_i(NEGD), Cf_i(NEGD), Sglobal(NEGD), CSglobal(NEGD), kvcommit(NEBD,NEBD), Secommit(NEBD),
ks(0), Ki(0), XAxis(3), YAxis(3), ZAxis(3), sr(1,6), ss(1,1), T(5), V(5), V2(5), G(5,5), GT(5,5), G2(5,5), G2T(5,5), H(5,5), H12(5,5), H22(5,5), Md(5,5), Hinv(5,5),Kg(5,5), fint2(5), fnat2(5), Tfnat2(5), fk(5), Cfk(5), fk_incr(5), ub(5), Tagg(tag), itr(0), cnvg(0)
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

   for(int i = 0; i < numSec; i++){
       dhat[i].Zero();
       duhat[i].Zero();
       sduhat[i].Zero();
   }

   //Fill in transformation matrix
   R[0][0]   = 0.0; R[0][1] = 0.0; R[0][2] = 0.0;
   R[1][0]   = 0.0; R[1][1] = 0.0; R[1][2] = 0.0;
   R[2][0]   = 0.0; R[2][1] = 0.0; R[2][2] = 0.0;

   CR[0][0]   = 0.0; CR[0][1] = 0.0; CR[0][2] = 0.0;
   CR[1][0]   = 0.0; CR[1][1] = 0.0; CR[1][2] = 0.0;
   CR[2][0]   = 0.0; CR[2][1] = 0.0; CR[2][2] = 0.0;

   ss.Zero();
   sr.Zero();

   T.Zero();
   G.Zero();
   GT.Zero();
   G2.Zero();
   G2T.Zero();
   H.Zero();
   H12.Zero();
   H22.Zero();
   Md.Zero();
   Hinv.Zero();
   fint2.Zero();
   fnat2.Zero();
   Tfnat2.Zero();

   kv.Zero();
   Kg.Zero();
   fk.Zero();
   Cfk.Zero();
   df_i.Zero();
   f_i.Zero();
   Cf_i.Zero();
   fk_incr.Zero();
   ub.Zero();
   V.Zero();
   V2.Zero();
   deflength = 0.0;
}

RCFTSTLGMBeamColumn3D::~RCFTSTLGMBeamColumn3D()
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

   if(nldhatT != 0)
     delete [] nldhatT;

   if(dhat != 0)
     delete [] dhat;

   if(duhat != 0)
     delete [] duhat;

   if(sduhat != 0)
     delete [] sduhat;

   if(nd1 != 0)
     delete [] nd1;

   if(nd2 != 0)
     delete [] nd2;

   if(nd1T != 0)
     delete [] nd1T;

   if(nd2T != 0)
     delete [] nd2T;

   if(nd1Tf != 0)
     delete [] nd1Tf;

   if(nd1Tfnd1 != 0)
     delete [] nd1Tfnd1;
  
   if(nd1Tfnd2 != 0)
     delete [] nd1Tfnd2;

   if(DQ != 0)
     delete [] DQ;

   if(DSQ != 0)
     delete [] DSQ;

   if(DSQa != 0)
     delete [] DSQa;

   if(CDSQa != 0)
     delete [] CDSQa;

   if (gd_delta != 0)
     delete [] gd_delta;

   if (crdTransf)
     delete crdTransf;

   if(beamIntegr != 0)
     delete beamIntegr;

   if (Ki != 0)
     delete Ki;

}


int
RCFTSTLGMBeamColumn3D::getNumExternalNodes(void) const
{
   return 2;
}


const ID &
RCFTSTLGMBeamColumn3D::getExternalNodes(void)
{
   return connectedExternalNodes;
}


Node **
RCFTSTLGMBeamColumn3D::getNodePtrs(void)
{
   return theNodes;
}

int
RCFTSTLGMBeamColumn3D::getNumDOF(void)
{
   return NEGD;
}


void
RCFTSTLGMBeamColumn3D::setDomain(Domain *theDomain)
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

   CR[0][0] = XAxis(0); CR[0][1] = XAxis(1); CR[0][2] = XAxis(2);
   CR[1][0] = YAxis(0); CR[1][1] = YAxis(1); CR[1][2] = YAxis(2);
   CR[2][0] = ZAxis(0); CR[2][1] = ZAxis(1); CR[2][2] = ZAxis(2);

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
RCFTSTLGMBeamColumn3D::commitState()
{
#ifdef COMPOSITE_DEBUG
   ofstream dunhat;
   dunhat.open("dunhat.dat",ios::app); 
   
   ofstream output;
   output.open("stlcon.dat",ios::app);
   
   ofstream crv11;
   crv11.open("crrv11.dat",ios::app);
   
   ofstream crv12;
   crv12.open("crrv12.dat",ios::app);
   
   ofstream crv13;
   crv13.open("crrv13.dat",ios::app);
   
   ofstream crv14;
   crv14.open("crrv14.dat",ios::app);
   
   ofstream crv21;
   crv21.open("crrv21.dat",ios::app);
   
   ofstream crv22;
   crv22.open("crrv22.dat",ios::app);
   
   ofstream crv23;
   crv23.open("crrv23.dat",ios::app);
   
   ofstream crv24;
   crv24.open("crrv24.dat",ios::app);
   
   ofstream crv31;
   crv31.open("crrv31.dat",ios::app);
   
   ofstream crv32;
   crv32.open("crrv32.dat",ios::app);
   
   ofstream crv33;
   crv33.open("crrv33.dat",ios::app);
   
   ofstream crv34;
   crv34.open("crrv34.dat",ios::app);
   
   ofstream crv41;
   crv41.open("crrv41.dat",ios::app);
   
   ofstream crv42;
   crv42.open("crrv42.dat",ios::app);
   
   ofstream crv43;
   crv43.open("crrv43.dat",ios::app);
   
   ofstream crv44;
   crv44.open("crrv44.dat",ios::app);
   
   ofstream mdc;
   mdc.open("mdc.dat",ios::app);
   
   ofstream mdrc;
   mdrc.open("mdrc.dat",ios::app);
   
   mdrc<<"commit_State"<<endl;
#endif

   int err = 0;
   int i = 0;
   int j = 0;

   cnvg = 0;

   // call element commitState to do any base class stuff
   if ((err = this->Element::commitState()) != 0) {
     opserr << "RCFTBeamColumn3D::commitState () - failed in base class";
     return err;
   }

   do {
      err = sections[i++]->commitState();
#ifdef COMPOSITE_DEBUG
	  if( i == 1 ){
         mdrc>>sections[0]->getSectionTangent();
      }
      if( i == 2 ){
         mdrc>>sections[1]->getSectionTangent();
      }
      if( i == 3 ){
         mdrc>>sections[2]->getSectionTangent();
      }
      if( i == 4 ){
         mdrc>>sections[3]->getSectionTangent();
      }
#endif
   } while (err == 0 && i < numSections);

   if (err)
      return err;

   // commit the transformation between coord. systems
   if ((err = crdTransf->commitState()) != 0)
      return err;

   // commit the element variables state
   kvcommit = kv;
   Secommit = Se;
   CSglobal = Sglobal;

   Cf_i = f_i;
   Cfk = fk;
   ub.Zero();
   fnat2.Zero();
   fint2.Zero();
   for( i = 0; i < numSections; i++){
        sduhat[i] = sduhat[i] + duhat[i];
        dhat[i].Zero();
        duhat[i].Zero();
        DSQ[i].Zero();
        DQ[i].Zero();
        CDSQa[i] = DSQa[i];
   }

#ifdef COMPOSITE_DEBUG
   if( Tagg == 1 ){
   crv11<<sduhat[0](1)<<"  "<<DSQa[0](1)<<endl;
   crv12<<sduhat[1](1)<<"  "<<DSQa[1](1)<<endl;
   crv13<<sduhat[2](1)<<"  "<<DSQa[2](1)<<endl;
   crv14<<sduhat[3](1)<<"  "<<DSQa[3](1)<<endl;
   }
   if( Tagg == 2 ){
   crv21<<sduhat[0](1)<<"  "<<DSQa[0](1)<<endl;
   crv22<<sduhat[1](1)<<"  "<<DSQa[1](1)<<endl;
   crv23<<sduhat[2](1)<<"  "<<DSQa[2](1)<<endl;
   crv24<<sduhat[3](1)<<"  "<<DSQa[3](1)<<endl;
   }
   if( Tagg == 3 ){
   crv31<<sduhat[0](1)<<"  "<<DSQa[0](1)<<endl;
   crv32<<sduhat[1](1)<<"  "<<DSQa[1](1)<<endl;
   crv33<<sduhat[2](1)<<"  "<<DSQa[2](1)<<endl;
   crv34<<sduhat[3](1)<<"  "<<DSQa[3](1)<<endl;
   }
   if( Tagg == 4 ){
   crv41<<sduhat[0](1)<<"  "<<DSQa[0](1)<<endl;
   crv42<<sduhat[1](1)<<"  "<<DSQa[1](1)<<endl;
   crv43<<sduhat[2](1)<<"  "<<DSQa[2](1)<<endl;
   crv44<<sduhat[3](1)<<"  "<<DSQa[3](1)<<endl;
   }
#endif

   //this->calcDeformedLength();

   //crdTransf->update();

   //crdTransf->getLocalAxes(XAxis, YAxis, ZAxis);

   CR[0][0]   = R[0][0];
   CR[0][1]   = R[0][1];
   CR[0][2]   = R[0][2];

   CR[1][0]   = R[1][0];
   CR[1][1]   = R[1][1];
   CR[1][2]   = R[1][2];

   CR[2][0]   = R[2][0];
   CR[2][1]   = R[2][1];
   CR[2][2]   = R[2][2];

   itr = 0;
   return err;
}


int RCFTSTLGMBeamColumn3D::revertToLastCommit()
{
#ifdef COMPOSITE_DEBUG
   ofstream mdc;
   mdc.open("mdc.dat",ios::app);
#endif

   int err;
   int i = 0;

   cnvg = 1;

#ifdef COMPOSITE_DEBUG
   mdc<<"revert_To_LastCommit"<<endl;
#endif

   do {
      err = sections[i]->revertToLastCommit();

      ksa[i] = sections[i]->getSectionTangent();

#ifdef COMPOSITE_DEBUG
	  mdc>>ksa[i];
#endif

      invertMatrix(4,ksa[i],fsa[i]);

      for( int j = 0; j < 3; j++ ){
        for( int k = 0; k< 3; k++ ){
          fs[i](j,k) = fsa[i](j,k);
          ks[i](j,k) = ksa[i](j,k);
        }
      }

      DSQa[i] = CDSQa[i];

      DQ[i].Zero();

      DSQ[i].Zero();

      dhat[i].Zero();

      duhat[i].Zero();


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
   Sglobal = CSglobal;
   fk   = Cfk; 

   fint2.Zero();
   fnat2.Zero();
   ub.Zero();
   f_i = Cf_i;

   initialFlag = 0;

   R[0][0]   = CR[0][0];
   R[0][1]   = CR[0][1];
   R[0][2]   = CR[0][2];

   R[1][0]   = CR[1][0];
   R[1][1]   = CR[1][1];
   R[1][2]   = CR[1][2];

   R[2][0]   = CR[2][0];
   R[2][1]   = CR[2][1];
   R[2][2]   = CR[2][2];

   // this->update();

   return err;
}


int RCFTSTLGMBeamColumn3D::revertToStart()
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
RCFTSTLGMBeamColumn3D::getInitialStiff(void)
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
RCFTSTLGMBeamColumn3D::getTangentStiff(void){
#ifdef COMPOSITE_DEBUG
  ofstream stf;
  stf.open("stf.dat",ios::app);
#endif
  //  int i,j;  
  crdTransf->update();  // Will remove once we clean up the corotational 2d transformation -- MHS
  const Matrix &KV = crdTransf->getGlobalStiffMatrix(kv,fk);
#ifdef COMPOSITE_DEBUG
  for( i = 0; i < 12; i++ ){
     for( j = 0; j < 12; j++ ){
     //stf<<setprecision(15)<<KV(i,j)<<"  ";
     stf<<KV(i,j)<<"  ";
     }
     stf<<"\n";
  }
#endif
  return KV;
}

const Vector &
RCFTSTLGMBeamColumn3D::getResistingForce(void){
#ifdef COMPOSITE_DEBUG
  ofstream stf;
  stf.open("stf.dat",ios::app);
  stf>>Sglobal;
#endif
  return Sglobal;
}

void RCFTSTLGMBeamColumn3D::calcResistingForce(void){
  crdTransf->update();
  Vector p0(12);
  p0.Zero();
  Sg = crdTransf->getGlobalResistingForce(f_i, p0);
  Sglobal  =  Sg;
}

void
RCFTSTLGMBeamColumn3D::initializeSectionHistoryVariables (void){
  for (int i = 0; i < numSections; i++){
        ksa[i]       = Matrix(4,4);
        fsa[i]       = Matrix(4,4);
        ks[i]        = Matrix(3,3);
        fs[i]        = Matrix(3,3);
        dhat[i]      = Vector(3);
        duhat[i]     = Vector(3);
        sduhat[i]    = Vector(3); 
        nldhat[i]    = Matrix(3,5);
        nldhatT[i]   = Matrix(5,3);
        nd1[i]       = Matrix(3,5);
        nd2[i]       = Matrix(3,5);
        nd2T[i]      = Matrix(5,3);
        nd1T[i]      = Matrix(5,3);
        nd1Tf[i]     = Matrix(5,3);
        nd1Tfnd1[i]  = Matrix(5,5);
        nd1Tfnd2[i]  = Matrix(5,5); 
        DQ[i]        = Vector(3);
        DSQ[i]       = Vector(3);
        DSQa[i]      = Vector(3);
        CDSQa[i]     = Vector(3);
        gd_delta[i]  = Vector(3);
  }
}

/********* NEWTON , SUBDIVIDE AND INITIAL ITERATIONS *********************/
int RCFTSTLGMBeamColumn3D::update()
{
#ifdef COMPOSITE_DEBUG
  ofstream geom;
  geom.open("geom.dat",ios::app);
  
  ofstream unbal;
  unbal.open("unbal.dat",ios::app);
  
  ofstream mpls;
  mpls.open("mpls.dat",ios::app);
  
  ofstream lstiff;
  lstiff.open("lstiff.dat",ios::app);
  
  ofstream RQ;
  RQ.open("RQ.dat",ios::app);
  
  ofstream FS;
  FS.open("FS.dat",ios::app);
  
  ofstream eig;
  eig.open("eig.dat",ios::app);
#endif

  if( cnvg == 1 ){
    cnvg = 0;
    return 0;
  }

  int i,j,k;

  itr = itr + 1;

  this->calcDeformedLength();

  double L = getDeformedLength();

#ifdef COMPOSITE_DEBUG
  mpls<<"Length"<<endl;
  mpls<<L;
#endif

  double oneOverL  = 1.0/L;

  if( initialFlag == 2 )
  this->revertToLastCommit();

  const Vector &dub = getBasicIncrDeltaDisp();
  
  ub = ub + dub;

#ifdef COMPOSITE_DEBUG
  mpls<<"ub"<<endl;
  mpls>>ub;

  mpls<<"dub"<<endl;
  mpls>>dub;
#endif

  const Vector du = getLocalIncrDeltaDisp();

#ifdef COMPOSITE_DEBUG
  mpls<<"du"<<endl;
  mpls>>du;
#endif

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double wt[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wt);

  double temp_x, temp_A, temp_B;
//  double d_d[3];

#ifdef COMPOSITE_DEBUG
  mpls<<"dhat[i]"<<endl;
#endif
  for ( i = 0; i < numSections; i++ ){
     nldhat[i] = this->getNld_hat(i, ub);
     dhat[i] = this->getd_hat(i, ub);
#ifdef COMPOSITE_DEBUG
	 mpls>>dhat[i];
#endif
     nd1[i] = this->getNd1(i, ub);
     nd2[i] = this->getNd2(i);
     for( j = 0; j < 3; j++ ){
         for( k = 0; k < 5; k++ ){
               nd1T[i](k,j) = nd1[i](j,k);
         }
     }
     for( j = 0; j < 3; j++ ){
         for( k = 0; k < 5; k++ ){
               nd2T[i](k,j) = nd2[i](j,k);
               nldhatT[i](k,j) = nldhat[i](j,k);
         }
     }
  }

  V.Zero();
  for( i = 0; i < numSections; i++ ){
     V = V + L * wt[i] * nd1T[i] * (dhat[i] - duhat[i] - ( fs[i] * ( DQ[i] - DSQ[i] ) ) );
  }

  H.Zero();
  for( i = 0; i < numSections; i++ ){
     H = H + L * wt[i] *  nd1T[i] * fs[i] * nd1[i];
  }

  invertMatrix(5, H, Hinv);

  G.Zero();

  for( i = 0; i < numSections; i++ ){
       G = G + L * wt[i] * nd1T[i] * nldhat[i];
  }

  G2.Zero();

  for( i = 0; i < numSections; i++ ){
       G2 = G2 + L * wt[i] * nd2T[i] * nldhat[i];
  }

  for( i = 0; i < 5; i++ ){
     for( j = 0; j < 5; j++ ){
       G2T(i,j) = G2(j,i);
     }
  }

  Vector dfnat(5);

  dfnat  =  Hinv * V;

  fnat2   = fnat2 + dfnat;

#ifdef COMPOSITE_DEBUG
  mpls<<"fnat2"<<endl;
  mpls>>fnat2;
#endif

  Tfnat2 = Tfnat2 + fnat2;

#ifdef COMPOSITE_DEBUG
  mpls<<"duhat[i]"<<endl;
#endif
  for ( i = 0; i < numSections; i++){

     DQ[i] = nd1[i] * fnat2;

     Vector d_delta(3);

     Vector aggduhat(4);

     d_delta = fs[i] * ( DQ[i] - DSQ[i] );

     duhat[i] =  duhat[i] + d_delta;
  
     for( int m = 0; m < 3; m++ ){
         aggduhat(m) = duhat[i](m);
     }

     aggduhat(3) = 0.0;
   
#ifdef COMPOSITE_DEBUG
	 int res = sections[i]->setTrialSectionDeformation(aggduhat);

     DSQa[i] = sections[i]->getStressResultant();

     mpls>>DSQa[i];
#endif

     ksa[i] = sections[i]->getSectionTangent();

     invertMatrix(4,ksa[i],fsa[i]);

     for( j = 0; j < 3; j++ ){
       for( k = 0; k< 3; k++ ){
         fs[i](j,k) = fsa[i](j,k);
         ks[i](j,k) = ksa[i](j,k);
       }
     }

     DSQ[i] = ks[i] * duhat[i];

     DSQa[i](0) = CDSQa[i](0) + DSQ[i](0);
     DSQa[i](1) = CDSQa[i](1) + DSQ[i](1);
     DSQa[i](2) = CDSQa[i](2) + DSQ[i](2);

     gd_delta[i] = fs[i] * ( DQ[i] - DSQ[i] );  

  } // for ( i = 0; i < numSections; i++)

  V.Zero();

  for( i = 0; i < numSections; i++ ){
     V = V + L * wt[i] * nd1T[i] * (dhat[i] - duhat[i] - gd_delta[i]);
  }

  V2.Zero();

  for( i = 0; i < numSections; i++ ){
     V2 = V2 + L * wt[i] * nd2T[i] * (dhat[i] - duhat[i]);
  }

  H.Zero();

  for( i = 0; i < numSections; i++ ){
     H = H + L * wt[i] * nd1T[i] * fs[i] * nd1[i];
  }

  invertMatrix(5, H, Hinv);

  H12.Zero();

  for( i = 0; i < numSections; i++ ){
     H12 = H12 + L * wt[i] * nd1T[i] * fs[i] * nd2[i];
  }

  H22.Zero();

  for( i = 0; i < numSections; i++ ){
     H22 = H22 + L * wt[i] * nd2T[i] * fs[i] * nd2[i];
  }

  Matrix temp_Md(5,5);
  Md.Zero();

  /***************************************************/
  /* CALCULATE THE G and T MATRIX                    */
  /***************************************************/
  for( i = 0; i < numSections; i++ ){
     temp_x = L * xi[i];
     temp_A =  ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * L;
     temp_B =  ( - pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * L;
     temp_Md.Zero();
     temp_Md(0,1) = temp_A * ( dhat[i](1) - duhat[i](1) );
     temp_Md(0,2) = temp_A * ( dhat[i](2) - duhat[i](2) );
     temp_Md(0,3) = temp_B * ( dhat[i](1) - duhat[i](1) );
     temp_Md(0,4) = temp_B * ( dhat[i](2) - duhat[i](2) );
     Md = Md + L * wt[i] * temp_Md;
  }

  /***************************************************/
  /* CALCULATE THE G and T MATRIX                    */
  /***************************************************/
 
  Matrix GMH(5,5);

  Matrix GMHT(5,5);

  Matrix GT(5,5);

  GMH = G + Md - H12;

  for( i = 0; i < 5; i++ ){
     for( j = 0; j < 5; j++ ){
        GMHT(i,j) = GMH(j,i);
        GT(i,j) = G(j,i);
     }
  }

  Kg.Zero();

  for( i = 0; i < numSections; i++ ){
      Kg = Kg + L * wt[i] * this->getKg(i);
  }

  //mpls<<"Kg"<<endl;
  //mpls>>Kg;

  /***************************************************/
  /* CALCULATE STIFFNESS MATRIX WITHOUT TORSION TERM */
  /***************************************************/

  Matrix K_temp(5,5);

  K_temp = ( Kg + G2 + G2T - H22 ) + GMHT * Hinv * GMH;

  //K_temp =  Kg + GT * Hinv * G;

  //mpls<<"GMHT * Hinv * GMH"<<endl;
  //mpls>>GMHT * Hinv * GMH;

  //Matrix K_k(5,5);

  //K_k.Zero(); 

  //for( i = 0; i < numSections; i++ ){
  //   K_k = K_k + L * wt[i] * nldhatT[i] * ks[i] * nldhat[i];
  //}

  Matrix kt(6,6);

  kt.Zero();

  for( i = 0; i < 5; i++ ){
     for( j = 0; j< 5; j++ ){
         kt(i,j) = K_temp(i,j);
     }
  }

  kt(5,5) =  ksa[0](3,3) / L;

  kv = kt;

  double tempfint2[5];

  for( i = 0; i < 5; i++ ){
       tempfint2[i] = fint2(i);
  }

  fint2 = GT * fnat2 + V2 + GMHT * Hinv * V;

  //fint2 = GT * fnat2 +  GT * Hinv * V;

  //fint2.Zero();

  //for( i = 0; i < numSections; i++ ){
  //   fint2 = fint2 + L * wt[i] * nldhatT[i] * DSQ[i];
  // }

  //mpls<<"fint2"<<endl;
  //mpls>>fint2;

  Vector dfint(5);

  for( i = 0; i < 5; i++ ){
       dfint(i) = fint2(i) - tempfint2[i];
  }

  fk(0)  = fk(0) + dfint(0);
  fk(1)  = fk(1) + dfint(1);  
  fk(2)  = fk(2) + dfint(2);
  fk(3)  = fk(3) + dfint(3);
  fk(4)  = fk(4) + dfint(4);

  df_i(0) = - dfint(0);
  df_i(1) = (dfint(1) + dfint(3))/L; 
  df_i(2) = (dfint(2) + dfint(4))/L; 
  df_i(3) =  du(3) * ksa[0](3,3) / L - du(9) * ksa[0](3,3) / L;
  df_i(4) = - dfint(2);
  df_i(5) = dfint(1);
  df_i(6) = - df_i(0);
  df_i(7) = - df_i(1);
  df_i(8) = - df_i(2);
  df_i(9) = - du(3) * ksa[1](3,3) / L + du(9) * ksa[1](3,3) / L;
  df_i(10) = - dfint(4);
  df_i(11) = dfint(3);

  //mpls<<"df_i"<<endl;
  //mpls>>df_i;

  f_i(0) = f_i(0) + df_i(0);
  f_i(1) = f_i(1) + df_i(1);
  f_i(2) = f_i(2) + df_i(2);
  f_i(3) = f_i(3) + df_i(3);
  f_i(4) = f_i(4) + df_i(4);
  f_i(5) = f_i(5) + df_i(5);
  f_i(6) = f_i(6) + df_i(6);
  f_i(7) = f_i(7) + df_i(7);
  f_i(8) = f_i(8) + df_i(8);
  f_i(9) = f_i(9) + df_i(9);
  f_i(10) = f_i(10) + df_i(10);
  f_i(11) = f_i(11) + df_i(11);

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

  return 0;
}


const Matrix &
RCFTSTLGMBeamColumn3D::getMass(void)
{
  theMatrix.Zero();

  return theMatrix;
}



void
RCFTSTLGMBeamColumn3D::zeroLoad(void)
{
	
}

int
RCFTSTLGMBeamColumn3D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  return 0;
}


int
RCFTSTLGMBeamColumn3D::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}

const Vector &
RCFTSTLGMBeamColumn3D::getResistingForceIncInertia()
{
  // Compute the current resisting force
  theVector = this->getResistingForce();

  // add the damping forces if rayleigh damping
  if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
     theVector += this->getRayleighDampingForces();

  return theVector;
}



bool
RCFTSTLGMBeamColumn3D::isSubdomain(void)
{
    return false;
}


void
RCFTSTLGMBeamColumn3D::Print(OPS_Stream &s, int flag)
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


OPS_Stream &operator<<(OPS_Stream &s, RCFTSTLGMBeamColumn3D &E)
{
    E.Print(s);
    return s;
}


int
RCFTSTLGMBeamColumn3D::displaySelf(Renderer &theViewer, int displayMode, float fact)
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
RCFTSTLGMBeamColumn3D::setResponse(char **argv, int argc, Information &eleInformation)
{
    // global force -
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
	|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
	return new ElementResponse(this, 1, theVector);
    else
      return 0;
}

int
RCFTSTLGMBeamColumn3D::getResponse(int responseID, Information &eleInfo)
{
    switch (responseID) {
      case 1:  // global forces
         return eleInfo.setVector(this->getResistingForce());
      default:
         return -1;
  }
}

int
RCFTSTLGMBeamColumn3D::setParameter (const char **argv, int argc, Parameter &param)
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
		cerr << "RCFTSTLGMBeamColumn3D::setParameter() - could not set parameter. " << endl;
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
RCFTSTLGMBeamColumn3D::updateParameter (int parameterID, Information &info)
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
		cerr << "RCFTSTLGMBeamColumn3D::updateParameter() - could not update parameter. " << endl;
		return ok;
	}
	else {
		return ok;
	}
     }
     else {
	cerr << "RCFTSTLGMBeamColumn3D::updateParameter() - could not update parameter. " << endl;
	return -1;
     }
}

void
RCFTSTLGMBeamColumn3D::setSectionPointers(int numSec, RCFTSTLFiberSection3D **secPtrs)
{
  if (numSec > maxNumSections) {
    opserr << "Error: RCFTSTLGMBeamColumn3D::setSectionPointers -- max number of sections exceeded";
  }

  numSections = numSec;

  if (secPtrs == 0) {
    opserr << "Error: RCFTSTLGMBeamColumn3D::setSectionPointers -- invalid section pointer";
  }

  sections = new RCFTSTLFiberSection3D *[numSections];
  if (sections == 0) {
    opserr << "Error: RCFTSTLGMBeamColumn3D::setSectionPointers -- could not allocate section pointers";
  }

  for (int i = 0; i < numSections; i++) {

    if (secPtrs[i] == 0) {
      opserr << "Error: RCFTSTLGMBeamColumn3D::setSectionPointers -- null section pointer " << i << endln;
    }

    sections[i] = (RCFTSTLFiberSection3D*) secPtrs[i]->getCopy();

    if (sections[i] == 0) {
      opserr << "Error: RCFTSTLGMBeamColumn3D::setSectionPointers -- could not create copy of section " << i << endln;
    }

  }

  ks  = new Matrix [numSections];
  fs  = new Matrix [numSections];
  ksa = new Matrix [numSections];
  fsa = new Matrix [numSections];
  dhat = new Vector [numSections];
  duhat = new Vector [numSections];
  sduhat = new Vector [numSections];
  nd1 = new Matrix [numSections];
  nd2 = new Matrix [numSections];
  nldhat = new Matrix [numSections];
  nldhatT = new Matrix [numSections];
  nd1T = new Matrix [numSections];
  nd2T = new Matrix [numSections];
  nd1Tf = new Matrix [numSections];
  nd1Tfnd1 = new Matrix [numSections];
  nd1Tfnd2 = new Matrix [numSections];
  DQ = new Vector [numSections];
  DSQa = new Vector [numSections];
  DSQ = new Vector [numSections];
  CDSQa = new Vector [numSections];
  gd_delta = new Vector [numSections];

  if (ks == 0) {
    opserr << "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate ks array";
  }

  if (fs == 0) {
    opserr << "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate fs array";
  }

  if (ksa == 0) {
    opserr << "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate ksa array";
  }
  
  if (fsa == 0) {
    opserr << "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate fsa array";
  }
  
  if (nldhat == 0) {
    opserr << "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate nldhat array";
  }

  if (nldhatT == 0) {
    opserr << "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate nldhatT array";
  }
  
  if (dhat == 0) {
    opserr << "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate dhat array";
  }

  if (nd1 == 0) {
    opserr << "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate nd1 array";
  }
  
  if (nd1T == 0) {
    opserr <<  "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate nd1T array";
  }

  if (nd1Tf == 0) {
    opserr <<  "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate nd1Tf array";
  }

  if (nd1Tfnd1 == 0) {
    opserr <<  "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate nd1Tfnd1 array";
  }

  if (DSQ == 0) {
    opserr <<  "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate DSQ array";
  }

  if (DQ == 0) {
    opserr <<  "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate DQ array";
  }

  if (DSQa == 0) {
    opserr <<  "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate DSQa array";
  }

  if (CDSQa == 0) {
    opserr <<  "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate CDSQa array";
  }

  if (nd2T == 0) {
    opserr <<  "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate nd2T array";
  }

  if (nd2 == 0) {
    opserr <<  "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate nd2 array";
  }

  if (nd1Tfnd2 == 0) {
    opserr <<  "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate nd1Tfnd2 array";
  }

  if (gd_delta == 0) {
    opserr <<  "RCFTSTLGMBeamColumn3D::setSectionPointers -- failed to allocate gd_delta array";
  }

}

Vector
RCFTSTLGMBeamColumn3D::getLocalIncrDeltaDisp(void)
{
#ifdef COMPOSITE_DEBUG
	ofstream mpls;
    mpls.open("mpls.dat",ios::app);

    mpls<<"nodal_disp"<<endl;
#endif

    const Vector &disp1 = theNodes[0]->getIncrDeltaDisp();
    const Vector &disp2 = theNodes[1]->getIncrDeltaDisp();

#ifdef COMPOSITE_DEBUG
	mpls>>disp1;
    mpls>>disp2;
#endif

    double ug[12];
    for (int i = 0; i < 6; i++) {
        ug[i]   = disp1(i);
        ug[i+6] = disp2(i);
    }

    double L = getDeformedLength();
    
    double oneOverL = 1.0/L;

    Vector ul(12);

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

    return ul;
}

const Vector& 
RCFTSTLGMBeamColumn3D::getBasicIncrDisp(void)
{
    return ub;									    
}	


Vector
RCFTSTLGMBeamColumn3D::getBasicIncrDeltaDisp(void)
{
   const Vector &disp1 = theNodes[0]->getIncrDeltaDisp();
   const Vector &disp2 = theNodes[1]->getIncrDeltaDisp();

   double ug[12];
   for (int i = 0; i < 6; i++) {
       ug[i]   = disp1(i);
       ug[i+6] = disp2(i);
   }

   double L = getDeformedLength();

   double oneOverL = 1.0/L;

   Vector dub(5);

   double ul[12];

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
    
   //dub(0) = d6_0;
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

   return dub;									    
}

Vector
RCFTSTLGMBeamColumn3D::getd_hat(int sec, const Vector &v)
{
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
   temp_C = 1 / L;
   temp_D = - 8 * temp_x / ( L * L ) + 4 / L;
   temp_E = - 4 / L + 6 * temp_x / ( L * L );
   temp_F =  - 2 / L + 6 * temp_x / ( L * L );

   D_hat(0) =  temp_C * v(0) +
               0.5 * ( temp_C * temp_C * v(0) ) * v(0) +
               0.5 * ( temp_A * temp_A * v(1) + temp_A * temp_B * v(3) ) * v(1) +
               0.5 * ( temp_A * temp_A * v(2) + temp_A * temp_B * v(4) ) * v(2) +
               0.5 * ( temp_A * temp_B * v(1) + temp_B * temp_B * v(3) ) * v(3) +
               0.5 * ( temp_A * temp_B * v(2) + temp_B * temp_B * v(4) ) * v(4); 

   //D_hat(0) =  temp_C * v(0) +
   //            0.5 * ( temp_A * temp_A * v(1) + temp_A * temp_B * v(3) ) * v(1) +
   //            0.5 * ( temp_A * temp_A * v(2) + temp_A * temp_B * v(4) ) * v(2) +
   //            0.5 * ( temp_A * temp_B * v(1) + temp_B * temp_B * v(3) ) * v(3) +
   //            0.5 * ( temp_A * temp_B * v(2) + temp_B * temp_B * v(4) ) * v(4);

   D_hat(1) =  temp_E * v(1) + temp_F * v(3);

   D_hat(2) =  temp_E * v(2) + temp_F * v(4);

   return D_hat;
}

Matrix
RCFTSTLGMBeamColumn3D::getKg(int sec)
{
   double L = getDeformedLength();
   double oneOverL  = 1.0/L;

   double xi[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double wt[maxNumSections];
   beamIntegr->getSectionWeights(numSections, L, wt);

   double temp_x, temp_A, temp_B;

   temp_x = L * xi[sec];

   Matrix kg(5,5);

   for(int i = 0; i < 5; i++){
       for(int j = 0; j < 5; j++){
           kg(i,j) = 0.0;
       }
   }

   temp_A = 1 - 4 * temp_x / L + 3 * ( temp_x * temp_x ) / ( L * L );
   temp_B = - 2 * temp_x / L + 3 * ( temp_x * temp_x ) / ( L * L );

   double Ps = DSQa[sec](0);

   kg(1,1) = Ps * temp_A * temp_A;
   kg(1,3) = Ps * temp_A * temp_B;
   kg(2,2) = Ps * temp_A * temp_A;
   kg(2,4) = Ps * temp_A * temp_B;
   kg(3,1) = Ps * temp_A * temp_B;
   kg(3,3) = Ps * temp_B * temp_B;
   kg(4,2) = Ps * temp_A * temp_B;
   kg(4,4) = Ps * temp_B * temp_B;

   return kg;
}


Matrix
RCFTSTLGMBeamColumn3D::getNld_hat(int sec, const Vector &v){
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
   temp_C =  1 / L;
   temp_D = - 8 * temp_x / ( L * L ) + 4 / L;
   temp_E = - 4 / L + 6 * temp_x / ( L * L );
   temp_F =  - 2 / L + 6 * temp_x / ( L * L );
 
   Matrix Nld_hat(3,5);

   for(i = 0; i < 3; i++){
       for(j = 0; j < 5; j++){
               Nld_hat(i,j) = 0.0;
       }
   }

   Nld_hat(0,0) = temp_C + temp_C * temp_C * v(0);
   //Nld_hat(0,0) = temp_C;
   Nld_hat(0,1) =  temp_A * temp_A * v(1) + temp_A * temp_B * v(3);
   Nld_hat(0,2) =  temp_A * temp_A * v(2) + temp_A * temp_B * v(4);
   Nld_hat(0,3) =  temp_A * temp_B * v(1) + temp_B * temp_B * v(3);
   Nld_hat(0,4) =  temp_A * temp_B * v(2) + temp_B * temp_B * v(4);
   Nld_hat(1,1) = temp_E;
   Nld_hat(1,3) = temp_F;
   Nld_hat(2,2) = temp_E;
   Nld_hat(2,4) = temp_F;

   return Nld_hat;
}

Matrix
RCFTSTLGMBeamColumn3D::getNd2(int sec){
   double L = getDeformedLength();
   double oneOverL  = 1.0/L;

   double xi[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double wt[maxNumSections];
   beamIntegr->getSectionWeights(numSections, L, wt);

   double temp_x, temp_A, temp_B;
//, temp_C, temp_D
   temp_x = L * xi[sec];

   Matrix Nd2(3,6);

   for(int i = 0; i < 5; i++){
       for(int j = 0; j < 5; j++){
           Nd2(i,j) = 0.0;
       }
   }

   temp_A = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) );
   temp_B = L * ( -pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) );

   double P = DSQa[sec](0);

   Nd2(1,2) = P * temp_A;
   Nd2(1,4) = P * temp_B;
   Nd2(2,1) = P * temp_A;
   Nd2(2,3) = P * temp_B;

   return Nd2;
}


Matrix
RCFTSTLGMBeamColumn3D::getNd1(int sec, const Vector &v){
   double L = getDeformedLength();
   double oneOverL  = 1.0/L;

   double xi[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double wt[maxNumSections];
   beamIntegr->getSectionWeights(numSections, L, wt);

   double temp_x, temp_A, temp_B;

   temp_x = L * xi[sec];

   temp_A = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[2]
           + L * ( - pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[4];

   temp_B = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[1]
           + L * ( - pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[3];

   Matrix Nd1(3,5);

   for(int i = 0; i < 3; i++){
       for(int j = 0; j < 5; j++){
           Nd1(i,j) = 0.0;
       }
   }

   Nd1(0,0)   = 1.0;
   Nd1(1,0)   = temp_A;
   Nd1(1,1)   = - temp_x / L + 1.0;
   Nd1(1,3)   = temp_x / L;
   Nd1(2,0)   = temp_B;
   Nd1(2,2)   = - temp_x / L + 1.0;
   Nd1(2,4)   = temp_x / L;

   return Nd1;
}
	
void RCFTSTLGMBeamColumn3D::calcDeformedLength(void){
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


double RCFTSTLGMBeamColumn3D::getDeformedLength(void){
   return deflength; 
}

int
RCFTSTLGMBeamColumn3D::sendSelf(int commitTag, Channel &theChannel){  
   return 0;
}
    
int
RCFTSTLGMBeamColumn3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker){
   return 0;
}

double 
RCFTSTLGMBeamColumn3D::getminEigenValue(int n, double *b)
{
   //this funtion dtermines the minimum eigen value of the global stiffness matrix
   //global stiffness matrix should be entered in vector form and starting index is zero
   int i, j, p, q, u, w, t, s;
   double eig, fm, cn, sn, omega, x, y, d;
   double eps = 1.e-3;

   double *a = new double[ n * n ];
   for ( i = 0; i < n * n; i++ ) *(a+i) = * (b+i);

   while (1){
      fm = 0.0;
      for ( i = 0; i < n; i++ )
          for ( j = 0; j < n; j++ ){
             d = fabs(a[i*n+j]);
             if ( ( i != j ) && ( d > fm ) ) { fm = d; p = i; q= j; }
          }
          if ( fm < eps ){
             eig = a[0];
             for ( i = 1; i < n; i++ ) if ( a[i*n+i] < eig ) eig = a[i*n+i];
             return (eig);
          }
          u = p * n + q; w = p * n + p;
          t = q * n + p; s = q * n + q;

          x = -a[u]; y = (a[s]-a[w])/2.0;
          omega = x/sqrt(x*x+y*y);

          if ( y < 0.0 ) omega = -omega;
             sn = 1.0 + sqrt(1.0-omega*omega);
             sn = omega / sqrt(2.0*sn);
             cn = sqrt(1.0-sn*sn);
             fm = a[w];
             a[w] = fm*cn*cn + a[s]*sn*sn + a[u]*omega;
             a[s] = fm*sn*sn + a[s]*cn*cn - a[u]*omega;
             a[u] = 0.0; a[t] = 0.0;
             for ( j = 0; j < n; j++ )
                 if ( ( j != p ) && ( j != q ) ){
                     u = p * n + j; w = q * n + j;
                     fm = a[u];
                     a[u] = fm * cn + a[w] * sn;
                     a[w] = -fm * sn + a[w] * cn;
                 }
             for ( i = 0; i < n; i++ )
                 if ( ( i != p ) && ( i != q ) ){
                     u = i * n + p; w = i * n + q;
                     fm = a[u];
                     a[u] = fm * cn + a[w] * sn;
                     a[w] = -fm * sn + a[w] * cn;
                 }
   }
   eig = a[0];
   for ( i = 1; i < n; i++ ) if ( a[i*n+i] < eig ) eig = a[i*n+i];
    delete a;
   return (eig);
}






    
									    

