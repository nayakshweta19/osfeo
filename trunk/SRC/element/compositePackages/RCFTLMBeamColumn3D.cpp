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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/element/RCFTLMBeamColumn3D.cpp,v $



#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <iomanip>

#include <Information.h>
#include <Parameter.h>
#include <RCFTLMBeamColumn3D.h>
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

Matrix RCFTLMBeamColumn3D::theMatrix(18,18);
Vector RCFTLMBeamColumn3D::theVector(18);
double RCFTLMBeamColumn3D::workArea[400];

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
// CONSTRUCTOR FOR PARALLEL PROCESSING

RCFTLMBeamColumn3D::RCFTLMBeamColumn3D():
Element(0,ELE_TAG_RCFTLMBeamColumn3D), connectedExternalNodes(2),
beamIntegr(0), numSections(0), sections(0), crdTransf(0),
rho(0.0), 
initialFlag(0),
kv(NEBD,NEBD), Se(NEBD), Sg(NEGD), Sglobal(NEGD), Sgb(NEGD), kvcommit(NEBD,NEBD), Secommit(NEBD),
fs(0), ks(0), fsa(0), ksa(0), nldhat(0), nldhatT(0), nd1(0), nd2(0), dhat(0), nd1T(0), nd2T(0), nd1Tf(0), 
nd1Tfnd1(0), nd1Tfnd2(0), duhat(0), sduhat(0), gd_delta(0), DQ(0), DSQ(0), DSQa(0), CDSQa(0), 
Ki(0), XAxis(3), YAxis(3), ZAxis(3), V(12), V2(12), G(12,12), Gsc(12,12), H(12,12),
H2(12,12), Hinv(12,12), fnat2(12), Tfnat2(12), fint2(12), Md(12,12), sr(2,12), ss(2,2), fk(24), fk_incr(24), ub(12), Tagg(0),
ugti(9), ugtj(9), T(12), S(12), df_i(18), itr(0), str_f4(0), f4(0), sf4(0), str_f4inv(0), d4(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;

  Sg.Zero();
  Sgb.Zero();
  Sglobal.Zero();
  kv.Zero();
  sr.Zero();
  ss.Zero();
  V.Zero();
  T.Zero();
  S.Zero();
  V2.Zero();
  G.Zero();
  Gsc.Zero();
  H.Zero();
  H2.Zero();
  Hinv.Zero();
  Md.Zero();
  fk.Zero();
  fk_incr.Zero();
  fint2.Zero();
  fnat2.Zero();
  Tfnat2.Zero();
  ub.Zero();
  ugti.Zero();
  ugtj.Zero();
  df_i.Zero();
  deflength = 0.0;

}


// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points.
// allocates the necessary space needed by each object
RCFTLMBeamColumn3D::RCFTLMBeamColumn3D (int tag, int nodeI, int nodeJ,
                                      int numSec, RCFTAggregator **sec,
                                      BeamIntegration &bi,
                                      CrdTransf &coordTransf, double massDensPerUnitLength):
Element(tag,ELE_TAG_RCFTLMBeamColumn3D), connectedExternalNodes(2),
beamIntegr(0), numSections(0), sections(0), crdTransf(0),
rho(massDensPerUnitLength),
initialFlag(0),
kv(NEBD,NEBD), Se(NEBD), Sg(NEGD), Sglobal(NEGD), Sgb(NEGD), kvcommit(NEBD,NEBD), Secommit(NEBD),
fs(0), ks(0), fsa(0), ksa(0), nldhat(0), nldhatT(0), nd1(0), nd2(0), dhat(0), nd1T(0), nd2T(0), nd1Tf(0),
nd1Tfnd1(0), nd1Tfnd2(0), duhat(0), sduhat(0), gd_delta(0), DQ(0), DSQ(0), CDSQa(0),  
Ki(0), XAxis(3), YAxis(3), ZAxis(3), V(12), V2(12), G(12,12), Gsc(12,12), H(12,12),
H2(12,12), Hinv(12,12), fnat2(12), Tfnat2(12), fint2(12), Md(12,12), sr(2,12), ss(2,2), fk(24), fk_incr(24), ub(12), Tagg(tag),
ugti(9), ugtj(9), T(12), S(12), df_i(18), itr(0), str_f4(0), sf4(0), f4(0), str_f4inv(0), d4(0)
{

   theNodes[0] = 0;
   theNodes[1] = 0;

   connectedExternalNodes(0) = nodeI;
   connectedExternalNodes(1) = nodeJ;

   // get copy of the sections

   beamIntegr = bi.getCopy();
   if(beamIntegr == 0) {
     opserr<<"Error: RCFTLMBeamColumn3D::RCFTLMBeamColumn3D: could not create copy of beam integration object" << endln;
     exit(-1);
   }

   // get copy of the transformation object
   
   crdTransf = coordTransf.getCopy3d();

   //deflength = crdTransf->getInitialLength();
   
   if (crdTransf == 0) {
      opserr << "Error: RCFTLMBeamColumn3D::RCFTLMBeamColumn3D: could not create copy of coordinate transformation object" << endln;
      exit(-1);
   }

   this->setSectionPointers(numSec,sec);

   for(int i = 0; i < numSec; i++){
       dhat[i].Zero();
   }

   //Fill in transformation matrix
   R[0][0]   = 0.0; R[0][1] = 0.0; R[0][2] = 0.0;
   R[1][0]   = 0.0; R[1][1] = 0.0; R[1][2] = 0.0;
   R[2][0]   = 0.0; R[2][1] = 0.0; R[2][2] = 0.0;

   ss.Zero();
   sr.Zero();
   
   kv.Zero();
   V.Zero();
   T.Zero();
   S.Zero();
   Sgb.Zero();
   Sglobal.Zero();
   V2.Zero();
   G.Zero();
   Gsc.Zero();
   H.Zero();
   H2.Zero();
   Hinv.Zero();
   Md.Zero();
   fk.Zero();
   fk_incr.Zero();
   ub.Zero();
   ugti.Zero();
   ugtj.Zero();
   df_i.Zero();
}

// ~RCFT_BeamColumn3D():
// 	destructor
//  delete must be invoked on any objects created by the object
RCFTLMBeamColumn3D::~RCFTLMBeamColumn3D()
{

   if (sections) {
      for (int i=0; i < numSections; i++)
         if (sections[i])
            delete sections[i];
      delete [] sections;
   }

   if(nd1 != 0)
     delete [] nd1;

   if(nldhat != 0)
     delete [] nldhat;

   if(nldhatT != 0)
     delete [] nldhatT;

   if(nd2 != 0)
     delete [] nd2;

   if(nd2T != 0)
     delete [] nd2T;

   if(nd1T != 0)
     delete [] nd1T;

   if(nd1Tf != 0)
     delete [] nd1Tf;

   if(nd1Tfnd1 != 0)
     delete [] nd1Tfnd1;

   if(nd1Tfnd2 != 0)
     delete [] nd1Tfnd2;

   if(dhat != 0)
     delete [] dhat;

   if(fsa != 0)
     delete [] fsa;

   if(ksa != 0)
     delete [] ksa;

   if(fs != 0)
     delete [] fs;

   if(ks != 0)
     delete [] ks;

   if(duhat != 0)
     delete [] duhat;

   if(sduhat != 0)
     delete [] sduhat;

   if(gd_delta != 0)
     delete [] gd_delta;
  
   if (crdTransf)
     delete crdTransf;

   if(beamIntegr != 0)
     delete beamIntegr;

   if (Ki != 0)
     delete Ki;

   if (DQ != 0)
     delete [] DQ;

   if(DSQ != 0)
     delete [] DSQ;

   if(DSQa != 0)
     delete [] DSQa;

   if(CDSQa != 0)
     delete [] CDSQa;

   if(str_f4 != 0)
     delete [] str_f4;

   if(f4 != 0)
     delete [] f4;

   if(sf4 != 0)
     delete [] sf4;

   if(str_f4inv != 0)
     delete [] str_f4inv;

   if(d4 != 0)
     delete [] d4;
}


int
RCFTLMBeamColumn3D::getNumExternalNodes(void) const
{
   return 2;
}


const ID &
RCFTLMBeamColumn3D::getExternalNodes(void)
{
   return connectedExternalNodes;
}


Node **
RCFTLMBeamColumn3D::getNodePtrs(void)
{
   return theNodes;
}

int
RCFTLMBeamColumn3D::getNumDOF(void)
{
   return NEGD;
}


void
RCFTLMBeamColumn3D::setDomain(Domain *theDomain)
{
   // check Domain is not null - invoked when object removed from a domain
 
   if (theDomain == 0)
   {
      theNodes[0] = 0;
      theNodes[1] = 0;

      opserr << "RCFTLMBeamColumn3D::setDomain:  theDomain = 0 ";
      exit(0);
   }

   // get pointers to the nodes

   int Nd1 = connectedExternalNodes(0);
   int Nd2 = connectedExternalNodes(1);

   theNodes[0] = theDomain->getNode(Nd1);
   theNodes[1] = theDomain->getNode(Nd2);

   if (theNodes[0] == 0)
   {
      opserr << "RCFTLMBeamColumn3D::setDomain: Nd1: ";
      opserr << Nd1 << "does not exist in model\n";
      exit(0);
   }

   if (theNodes[1] == 0)
   {
      opserr << "RCFTLMBeamColumn3D::setDomain: Nd2: ";
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
      opserr << "RCFTLMBeamColumn3D::setDomain(): Nd2 or Nd1 incorrect dof ";
      exit(0);
   }

   // initialize the transformation
   if (crdTransf->initialize(theNodes[0], theNodes[1]))
   {
      opserr << "RCFTLMBeamColumn3D::setDomain(): Error initializing coordinate transformation";
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
      opserr << "RCFTLMBeamColumn3D::setDomain(): Zero element length:" << this->getTag();
      exit(0);
   }

   if (initialFlag == 0)
     this->initializeSectionHistoryVariables();

}



int
RCFTLMBeamColumn3D::commitState()
{
   //ofstream dunhat;
   //dunhat.open("dunhat.dat",ios::app); 
   //
   //ofstream output;
   //output.open("stlcon.dat",ios::app);
   //
   //ofstream CDQ;
   //CDQ.open("CDQ.dat",ios::app);
   //
   //ofstream crv11;
   //crv11.open("crrv11.dat",ios::app);
   //
   //ofstream crv12;
   //crv12.open("crrv12.dat",ios::app);
   //
   //ofstream crv13;
   //crv13.open("crrv13.dat",ios::app);
   //
   //ofstream crv14;
   //crv14.open("crrv14.dat",ios::app);
   //
   //ofstream crv21;
   //crv21.open("crrv21.dat",ios::app);
   //
   //ofstream crv22;
   //crv22.open("crrv22.dat",ios::app);
   //
   //ofstream crv23;
   //crv23.open("crrv23.dat",ios::app);
   //
   //ofstream crv24;
   //crv24.open("crrv24.dat",ios::app);
   //
   //ofstream crv31;
   //crv31.open("crrv31.dat",ios::app);
   //
   //ofstream crv32;
   //crv32.open("crrv32.dat",ios::app);
   //
   //ofstream crv33;
   //crv33.open("crrv33.dat",ios::app);
   //
   //ofstream crv34;
   //crv34.open("crrv34.dat",ios::app);

   int err = 0;
   int i = 0;
   int j = 0;

   // call element commitState to do any base class stuff
   if ((err = this->Element::commitState()) != 0) {
     opserr << "RCFTLMBeamColumn3D::commitState () - failed in base class";
     return err;
   }

   do {
      //output<<"section #"<<i<<endl;	   
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

   fnat2.Zero();
   fint2.Zero();
   ub.Zero();
   
   for( i = 0; i < numSections; i++){
        //CDQ<<"section "<<i<<" th"<<endl;
        //CDQ>>DQ[i];
        //CDQ>>DSQ[i];
        //CDQ<<(DQ[i]-DSQ[i]).Norm()<<endl;
        sduhat[i] = sduhat[i] + duhat[i];
		duhat[i].Zero();
        dhat[i].Zero();
		DSQ[i].Zero();
        DQ[i].Zero();
        f4[i].Zero();
        sf4[i].Zero();
        d4[i].Zero(); 
		CDSQa[i] = DSQa[i];
   }
   //if( Tagg == 1 ){
   //crv11<<sduhat[0](1)<<"  "<<DSQa[0](1)<<"   "<<sduhat[0](4)<<"  "<<DSQa[0](7)<<"  "<<ks[0](1,1)<<endl;
   //crv12<<sduhat[1](1)<<"  "<<DSQa[1](1)<<"   "<<sduhat[1](4)<<"  "<<DSQa[1](7)<<"  "<<ks[1](1,1)<<endl;
   //crv13<<sduhat[2](1)<<"  "<<DSQa[2](1)<<"   "<<sduhat[2](4)<<"  "<<DSQa[2](7)<<"  "<<ks[2](1,1)<<endl; 
   //crv14<<sduhat[3](1)<<"  "<<DSQa[3](1)<<"   "<<sduhat[3](4)<<"  "<<DSQa[3](7)<<"  "<<ks[3](1,1)<<endl;
   //}
   //if( Tagg == 2 ){
   //crv21<<sduhat[0](1)<<"  "<<DSQa[0](1)<<"   "<<sduhat[0](4)<<"  "<<DSQa[0](7)<<"  "<<ks[0](1,1)<<endl;
   //crv22<<sduhat[1](1)<<"  "<<DSQa[1](1)<<"   "<<sduhat[1](4)<<"  "<<DSQa[1](7)<<"  "<<ks[1](1,1)<<endl;
   //crv23<<sduhat[2](1)<<"  "<<DSQa[2](1)<<"   "<<sduhat[2](4)<<"  "<<DSQa[2](7)<<"  "<<ks[2](1,1)<<endl;
   //crv24<<sduhat[3](1)<<"  "<<DSQa[3](1)<<"   "<<sduhat[3](4)<<"  "<<DSQa[3](7)<<"  "<<ks[3](1,1)<<endl;
   //}
   //if( Tagg == 3 ){
   //crv31<<sduhat[0](1)<<"  "<<DSQa[0](1)<<"   "<<sduhat[0](4)<<"  "<<DSQa[0](7)<<"  "<<ks[0](1,1)<<endl;
   //crv32<<sduhat[1](1)<<"  "<<DSQa[1](1)<<"   "<<sduhat[1](4)<<"  "<<DSQa[1](7)<<"  "<<ks[1](1,1)<<endl;
   //crv33<<sduhat[2](1)<<"  "<<DSQa[2](1)<<"   "<<sduhat[2](4)<<"  "<<DSQa[2](7)<<"  "<<ks[2](1,1)<<endl;
   //crv34<<sduhat[3](1)<<"  "<<DSQa[3](1)<<"   "<<sduhat[3](4)<<"  "<<DSQa[3](7)<<"  "<<ks[3](1,1)<<endl;
   //}
   itr = 0;
   return err;
}


int RCFTLMBeamColumn3D::revertToLastCommit()
{
   int err;
   int i = 0;

   do {
      err = sections[i]->revertToLastCommit();

      fs[i]  = sections[i]->getSectionFlexibility();

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


int RCFTLMBeamColumn3D::revertToStart()
{
   // revert the sections state to start
   int err;
   int i = 0;

   do {
       fs[i].Zero();
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
RCFTLMBeamColumn3D::getInitialStiff(void)
{
  //ofstream newton; 
  //newton.open("newton.dat",ios::app);

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

  /************************************************************************/
  /* TO FACILITATE SIMPLER SYNTAX IN THE EQUATIONS, DEFINE SOME LOCAL	  */
  /* VARIABLES TO REPRESENT RIGIDITY TERMS				  */
  /************************************************************************/
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
  gj 	= kcrsi(6,6);
  kslip = kcrsi(12,12) * 2 * ( kcrsi(8,8) + kcrsi(7,7) );
  dkslip= kcrsj(12,12) * 2 * ( kcrsj(8,8) + kcrsj(7,7) ) - kcrsi(12,12) * 2 * ( kcrsi(8,8) + kcrsi(7,7) );

  ea_c 	 = kcrsi(0,0);
  eqy_c	 = kcrsi(0,2);
  eqz_c	 = kcrsi(1,0);
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
  eqy_s	 = kcrsi(3,5);
  eqz_s	 = kcrsi(3,4);
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

  //newton<<"\n INITIAL GLOBAL STIFFNESS MATRIX"<<Tagg<<endl;
  
  //newton>>(*Ki);
  
  return *Ki;

}

const Matrix &
RCFTLMBeamColumn3D::getTangentStiff(void)
{
  //ofstream mpls;
  //mpls.open("mpls.dat",ios::app);
  int i;  
  crdTransf->update();  // Will remove once we clean up the corotational 2d transformation -- MHS
  const Matrix &KV = crdTransf->getGlobalStiffMatrix(kv,fk);
	
  //mpls>>kv;
  return KV;
}

const Vector &
RCFTLMBeamColumn3D::getResistingForce(void)
{
  return Sglobal;
}

void RCFTLMBeamColumn3D::calcResistingForce(void)
{
  //ofstream unbal;
  //unbal.open("iforce.dat",ios::app);
  crdTransf->update();
  Vector p0(18);
  p0.Zero();
  Sgb = Sgb + df_i;
  Sg = crdTransf->getGlobalResistingForce(df_i, p0);
  Sglobal  = Sglobal + Sg;
  //unbal>>Sg;
  //unbal>>Sglobal;
}

void
RCFTLMBeamColumn3D::initializeSectionHistoryVariables (void){
    for (int i = 0; i < numSections; i++){
	fs[i]       = Matrix(6,6);
	ks[i]       = Matrix(6,6);
        fsa[i]      = Matrix(12,12);
        ksa[i]      = Matrix(12,12);
	dhat[i]     = Vector(6);
	nldhat[i]   = Matrix(6,12);
        nldhatT[i]   = Matrix(12,6);
        dhat[i]     = Vector(6);
        duhat[i]    = Vector(6);
        sduhat[i]   = Vector(6);
        sduhat[i].Zero();
        nd1[i]      = Matrix(6,12);
        nd2[i]      = Matrix(6,12);
        nd1T[i]     = Matrix(12,6);
        nd2T[i]     = Matrix(12,6);
        nd1Tf[i]    = Matrix(12,6);
        nd1Tfnd1[i] = Matrix(12,12);
        nd1Tfnd2[i] = Matrix(12,12);
        DQ[i]       = Vector(6);
        DSQ[i]      = Vector(6);
        DSQa[i]     = Vector(13);
        CDSQa[i]    = Vector(13);
        gd_delta[i] = Vector(6);
        str_f4[i]   = Matrix(4,4);
        f4[i]       = Vector(4);
        sf4[i]      = Vector(4);
        str_f4inv[i]= Matrix(4,4);
        d4[i]       = Vector(4);
    }
}

/********* NEWTON , SUBDIVIDE AND INITIAL ITERATIONS *********************/
int RCFTLMBeamColumn3D::update(){
  //ofstream geom;
  //geom.open("geom.dat",ios::app);
  //
  //ofstream unbal;
  //unbal.open("unbal.dat",ios::app);
  //
  //ofstream mpls;
  //mpls.open("mpls.dat",ios::app);
  //
  //ofstream lstiff;
  //lstiff.open("lstiff.dat",ios::app);
  //
  //ofstream newton19;
  //newton19.open("newton19.dat",ios::app);
  //
  //ofstream FS;
  //FS.open("FS.dat",ios::app);
  //
  //ofstream FN;
  //FN.open("FN.dat",ios::app);
  // 
  //ofstream VV;
  //VV.open("VV.dat",ios::app);
  //
  //ofstream dq;
  //dq.open("DQ.dat",ios::app);
  //
  //ofstream DH;
  //DH.open("DH.dat",ios::app);
  //
  //ofstream DU;
  //DU.open("DU.dat",ios::app);
  //
  //ofstream check3;
  //check3.open("check3.dat",ios::app);

  int i,j,k;

  itr = itr + 1;

  //mpls<<"ELEMENT NUMBER "<<Tagg<<" ITERATION "<<itr<<endl;

  //dq<<"ELEMENT NUMBER "<<Tagg<<endl;

  double L = crdTransf->getInitialLength();

  double oneOverL  = 1.0/L;

  if( initialFlag == 2 )
  this->revertToLastCommit();

  const Vector &dub = getBasicIncrDeltaDisp();
  
  ub = ub + dub;

  //DH>>ub;

  const Vector du = getLocalIncrDeltaDisp();

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double wt[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wt);

  int flag = 1;

  double tempfint2[12];

  if( flag == 0 ){

  /********************************************************************************/
  /*                                                                              */
  /* CALCULATE INCREMENTAL STRAIN VECTOR FROM INCREMENTAL NATURAL DEFORMATIONS    */
  /*                                                                              */
  /********************************************************************************/

  for ( i = 0; i < numSections; i++ ){
      nldhat[i] = this->getNld_hat(i, ub);
      dhat[i] = this->getd_hat(i, ub);
      nd1[i] = this->getNd1(i, ub);
      for( j = 0; j < 6; j++ ){
          for( k = 0; k < 12; k++ ){
                nd1T[i](k,j) = nd1[i](j,k);
          }
      }
  }

  V.Zero();
  
  double temp_x, temp_A, temp_B, temp_C, temp_D;
  double d_d[6];

  for( i = 0; i < numSections; i++ ){
       V = V + L * wt[i] * nd1T[i] * ( dhat[i] - duhat[i] -  ( fs[i] * ( DQ[i] - DSQ[i] ) ) );
  }

  H.Zero(); 

  for( i = 0; i < numSections; i++ ){
       nd1Tf[i] = nd1T[i] * fs[i];
       nd1Tfnd1[i] = nd1Tf[i] * nd1[i];
       H = H + L * wt[i] * nd1Tfnd1[i];
  }

  invertMatrix(12, H, Hinv);

  Matrix temp_G(12,12);

  G.Zero();

  for( i = 0; i < numSections; i++ ){
       temp_G = nd1T[i] * nldhat[i];
       G = G + L * wt[i] * temp_G;
  }

  Vector dfnat(12);

  dfnat  =  Hinv * V;

  fnat2   = fnat2 + dfnat;

  Tfnat2 = Tfnat2 + dfnat;

  for ( i = 0; i < numSections; i++){
     double slp_strn;
     double slp_force;

     DQ[i] = nd1[i] * fnat2;

     f4[i](0) = DQ[i](0);
     f4[i](1) = DQ[i](3);
     f4[i](2) = DQ[i](1) + DQ[i](4); 
     f4[i](3) = DQ[i](2) + DQ[i](5);

     str_f4[i](0,0) = ks[i](0,0);
     str_f4[i](1,1) = ks[i](3,3);
     str_f4[i](0,1) = str_f4[i](1,0) = 0.0;
     str_f4[i](0,2) = str_f4[i](2,0) = ks[i](0,1);
     str_f4[i](0,3) = str_f4[i](3,0) = ks[i](0,2);
     str_f4[i](1,2) = str_f4[i](2,1) = ks[i](3,4);
     str_f4[i](1,3) = str_f4[i](3,1) = ks[i](3,5);
     str_f4[i](2,2) = ks[i](1,1) + ks[i](4,4);
     str_f4[i](2,3) = str_f4[i](3,2) = ks[i](2,1) + ks[i](5,4);
     str_f4[i](3,3) = ks[i](2,2) + ks[i](5,5); 
          
     invertMatrix(4,str_f4[i],str_f4inv[i]);

     d4[i] =  str_f4inv[i] * ( f4[i] - sf4[i] );

     temp_x = L * xi[i];

     temp_A = - temp_x/L + 2 * temp_x * temp_x / ( L * L );

     temp_B = - 4 * temp_x * temp_x / ( L * L ) + 4 * temp_x / L;

     slp_strn = temp_A * ub(5) + temp_B * ub(11) - temp_A * ub(0) - temp_B * ub(10);

     Vector d_delta(6);

     Vector aggduhat(13);

     d_delta(0) = d4[i](0);
     d_delta(1) = d4[i](2); 
     d_delta(2) = d4[i](3); 
     d_delta(3) = d4[i](1); 
     d_delta(4) = d4[i](2);
     d_delta(5) = d4[i](3);

     duhat[i] =  duhat[i] + d_delta;

     for( int m = 0; m < 6; m++ ){
         aggduhat(m) = duhat[i](m);
     }

     aggduhat(6) = 0.0;
     aggduhat(7) = 0.0;
     aggduhat(8) = 0.0;
     aggduhat(9) = 0.0;
     aggduhat(10) = 0.0;
     aggduhat(11) = 0.0;
     aggduhat(12) = slp_strn;

     int res = sections[i]->setTrialSectionDeformation(aggduhat);

     DSQa[i] = sections[i]->getStressResultant();

     ksa[i] = sections[i]->getSectionTangent();

     invertMatrix(12,ksa[i],fsa[i]);

     for( j = 0; j < 6; j++ ){
       for( k = 0; k< 6; k++ ){
         fs[i](j,k) = fsa[i](j,k);
         ks[i](j,k) = ksa[i](j,k);
       }
     }

     str_f4[i](0,0) = ks[i](0,0);
     str_f4[i](1,1) = ks[i](3,3);
     str_f4[i](0,1) = str_f4[i](1,0) = 0.0;
     str_f4[i](0,2) = str_f4[i](2,0) = ks[i](0,1);
     str_f4[i](0,3) = str_f4[i](3,0) = ks[i](0,2);
     str_f4[i](1,2) = str_f4[i](2,1) = ks[i](3,4);
     str_f4[i](1,3) = str_f4[i](3,1) = ks[i](3,5);
     str_f4[i](2,2) = ks[i](1,1) + ks[i](4,4);
     str_f4[i](2,3) = str_f4[i](3,2) = ks[i](2,1) + ks[i](5,4);
     str_f4[i](3,3) = ks[i](2,2) + ks[i](5,5);

     invertMatrix(4,str_f4[i],str_f4inv[i]);

     if( i == 0 ){
         slp_force = 2 * ( ksa[i](7,7) + ksa[i](8,8) ) * ( xi[i+1] * L ) * 0.5 * DSQa[i](12);
     }
     else if( i == numSections - 1 ){
         slp_force = 2 * ( ksa[i](7,7) + ksa[i](8,8) ) * ( xi[numSections-1] * L 
	 	  - xi[numSections-2] * L ) * 0.5 * DSQa[i](12);
     }
     else{
         slp_force = 2 * ( ksa[i](7,7) + ksa[i](8,8) ) * ( ( xi[i] * L  - xi[i-1] * L ) * 0.5
                  + ( xi[i+1] * L  - xi[i] * L ) * 0.5 ) * DSQa[i](12);
     }

     slp_force =  2 * ( ksa[i](7,7) + ksa[i](8,8) ) * DSQa[i](12);
 
     DSQ[i] = ks[i]*duhat[i];

     //dq<<"DQ[i]"<<endl;
     //dq>>DQ[i];
     //dq<<"DSQ[i]"<<endl;
     //dq>>DSQ[i];

     sf4[i](0) = DSQ[i](0);
     sf4[i](1) = DSQ[i](3);
     sf4[i](2) = DSQ[i](1) + DSQ[i](4);
     sf4[i](3) = DSQ[i](2) + DSQ[i](5);
  
     Vector gd4(4);
     gd4 = str_f4inv[i] * ( f4[i] - sf4[i] );
     gd_delta[i](0) = gd4(0);
     gd_delta[i](1) = gd4(2);
     gd_delta[i](2) = gd4(3);
     gd_delta[i](3) = gd4(1);
     gd_delta[i](4) = gd4(2);
     gd_delta[i](5) = gd4(3);

  } // for ( i = 0; i < numSections; i++)

  V.Zero();

  for( i = 0; i < numSections; i++ ){
     V = V + L * wt[i] * nd1T[i] * (dhat[i] - duhat[i] - gd_delta[i]);
  }

  H.Zero();

  for( i = 0; i < numSections; i++ ){
     nd1Tf[i] = nd1T[i] * fs[i];
     nd1Tfnd1[i] = nd1Tf[i] * nd1[i];
     H = H + L * wt[i] * nd1Tfnd1[i];
  }

  invertMatrix(12, H, Hinv);


  Matrix temp_Gsc(12,12);

  Gsc.Zero();

  for( i = 0; i < numSections; i++ ){
       for( j = 0; j < 12; j++ ){
            for( k = 0; k < 12; k++ ){
                 temp_Gsc(j,k) = 0.0;
            }
       }
       //ksa[i] = sections[i]->getSectionTangent();
       temp_Gsc(0,0) = ksa[i](12,12);
       temp_Gsc(0,5) = - ksa[i](12,12);
       temp_Gsc(5,0) = - ksa[i](12,12);
       temp_Gsc(5,5) = ksa[i](12,12);
       for( j = 0; j < 12; j++ ){
           for( k = 0; k < 12; k++ ){
                Gsc(j,k) = Gsc(j,k) + L * wt[i] * temp_Gsc(j,k);
           }
       }
  }

  for( i = 0; i < 12; i++ ){
       tempfint2[i] = fint2(i);
  }

  Matrix GT(12,12);

  Matrix GscT(12,12);

  for( i = 0; i < 12; i++ ){
     for( j = 0; j < 12; j++ ){
        GT(i,j) = G(j,i);
        GscT(i,j) = Gsc(j,i);
     }
  } 
 
  /**************************************/
  /* CALCULATE THE INTERNAL LOAD VECTOR */
  /**************************************/

  fint2 = GT * fnat2 + GscT * ub + ( GT * Hinv ) * V;

  //fint2 = GT * fnat2 + V2 + ( GMHT * Hinv ) * V;

  }
  else if( flag == 1 ){
      for ( i = 0; i < numSections; i++ ){
      	 nldhat[i] = this->getNld_hat(i, ub);  
         for( j = 0; j < 12; j++ ){
            for( k = 0; k < 6; k++ ){
               nldhatT[i](j,k) = nldhat[i](k,j);     
            }
         }
     	 dhat[i] = this->getd_hat(i, ub);
      }
      for ( i = 0; i < numSections; i++ ){
         double temp_x = L * xi[i];
         double temp_A = - temp_x/L + 2 * temp_x * temp_x / ( L * L );
         double temp_B = - 4 * temp_x * temp_x / ( L * L ) + 4 * temp_x / L;
         double slp_strn = temp_A * ub(5) + temp_B * ub(11) - temp_A * ub(0) - temp_B * ub(10);
         Vector aggdhat(13);
         for( int m = 0; m < 6; m++ ){
            aggdhat(m) = dhat[i](m);
         }
         aggdhat(6) = 0.0;
         aggdhat(7) = 0.0;
         aggdhat(8) = 0.0;
         aggdhat(9) = 0.0;
         aggdhat(10) = 0.0;
         aggdhat(11) = 0.0;
         aggdhat(12) = slp_strn;
         int res = sections[i]->setTrialSectionDeformation(aggdhat);
         ksa[i] = sections[i]->getSectionTangent();
         invertMatrix(12,ksa[i],fsa[i]);
         for( j = 0; j < 6; j++ ){
            for( k = 0; k< 6; k++ ){
                fs[i](j,k) = fsa[i](j,k);
                ks[i](j,k) = ksa[i](j,k);
            }
         }
         DSQ[i] = ks[i] * dhat[i]; 
         //dq>>DSQ[i];
     }

     Vector Z(12);

     Z.Zero();

     for( i = 0; i < numSections; i++ ){
         Z = Z + L * wt[i] * nldhatT[i] * DSQ[i];
     }

     Matrix temp_Gsc(12,12);

     Gsc.Zero();

     for( i = 0; i < numSections; i++ ){
        for( j = 0; j < 12; j++ ){
            for( k = 0; k < 12; k++ ){
                 temp_Gsc(j,k) = 0.0;
            }
        }
        temp_Gsc(0,0) = ksa[i](12,12);
        temp_Gsc(0,5) = - ksa[i](12,12);
        temp_Gsc(5,0) = - ksa[i](12,12);
        temp_Gsc(5,5) = ksa[i](12,12);
        for( j = 0; j < 12; j++ ){
           for( k = 0; k < 12; k++ ){
                Gsc(j,k) = Gsc(j,k) + L * wt[i] * temp_Gsc(j,k);
           }
        }
      }
      for( i = 0; i < 12; i++ ){
         tempfint2[i] = fint2(i);
      }
      Matrix GscT(12,12);
      for( i = 0; i < 12; i++ ){
         for( j = 0; j < 12; j++ ){
            GscT(i,j) = Gsc(j,i);
         }
      }
      fint2 = Z;// + GscT * ub;
     // dq>>fint2;   
  }
  
  /***********************************************/
  /* CALCULATE INCREMENTAL INTERNAL LOAD VECTOR  */
  /***********************************************/
  Vector dfint(12);

  for( i = 0; i < 12; i++ ){
       dfint(i) = fint2(i) - tempfint2[i];
  }

  double theta_rigid_z, theta_rigid_y;

  theta_rigid_z = ( du(10) - du(1) ) / L;
  theta_rigid_y = ( du(11) - du(2) ) / L;

  /************************************************************************/
  /* NOW UPDATE THE INCREMENTAL ELEMENT FORCES INCLUDING SHEAR            */
  /************************************************************************/
  df_i(0)  = - dfint(5);
  df_i(1)  = (dfint(6) + dfint(8)) / L;
  df_i(2)  = (dfint(7) + dfint(9)) / L;
  df_i(3)  = du(3) * ksa[0](6,6) / L - du(12) * ksa[0](6,6) / L;
  df_i(4)  = - dfint(2) - dfint(7);
  df_i(5)  = dfint(1) + dfint(6);
  df_i(6)  = - dfint(0);
  df_i(7)  = (dfint(1) + dfint(3)) / L;
  df_i(8)  = (dfint(2) + dfint(4)) / L;
  df_i(9)  = dfint(5);
  df_i(10) = - df_i(1);
  df_i(11) = - df_i(2);
  df_i(12) = du(3) * ksa[0](6,6) / L - du(12) * ksa[0](6,6) / L;    
  df_i(13) = - dfint(4) - dfint(9);
  df_i(14) = dfint(3) + dfint(8);
  df_i(15) = dfint(0);
  df_i(16) = - df_i(7);
  df_i(17) = - df_i(8);

  this->calcResistingForce();

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

  fk_incr(0) = dfint(0);
  fk_incr(1) = dfint(5);
  fk_incr(2) = - dfint(2);
  fk_incr(3) = - dfint(1);
  fk_incr(4) = - dfint(7);
  fk_incr(5) = - dfint(6);
  fk_incr(6) = 0.0;
  fk_incr(8) = 0.0;
  fk_incr(7) = 0.0;
  fk_incr(9) = 0.0;
  fk_incr(10)= 0.0;
  fk_incr(11)= 0.0;
  fk(0) += fk_incr(0);
  fk(1) += fk_incr(1);
  fk(2) += fk_incr(2);
  fk(3) += fk_incr(3);
  fk(4) += fk_incr(4);
  fk(5) += fk_incr(5);
  fk(6) = DSQa[0](3);
  fk(7) = DSQa[0](5);
  fk(8) = DSQa[0](4);
  fk(9) = DSQa[0](9);
  fk(10)= DSQa[0](10);
  fk(11)= DSQa[0](11);

  fk_incr(12) = dfint(0);
  fk_incr(13) = dfint(5);
  fk_incr(14) = dfint(4);
  fk_incr(15) = dfint(3);
  fk_incr(16) = dfint(9);
  fk_incr(17) = dfint(8);
  fk_incr(18) = 0.0;
  fk_incr(20) = 0.0;
  fk_incr(19) = 0.0;
  fk_incr(21) = 0.0;
  fk_incr(23) = 0.0;
  fk_incr(22) = 0.0;
  fk(12) += fk_incr(12);
  fk(13) += fk_incr(13);
  fk(14) += fk_incr(14);
  fk(15) += fk_incr(15);
  fk(16) += fk_incr(16);
  fk(17) += fk_incr(17);
  fk(18) = DSQa[numSections-1](3);
  fk(19) = DSQa[numSections-1](4);
  fk(20) = DSQa[numSections-1](5);
  fk(21) = DSQa[numSections-1](9);
  fk(22) = DSQa[numSections-1](10);
  fk(23) = DSQa[numSections-1](11);

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

  for ( i = 0; i < 13; i++ ){
        for ( j = ( i + 1 ); j < 14; j++ ){
                kt(j,i) = kt(i,j);
        }
  }

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

  //DH>>kt;

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
RCFTLMBeamColumn3D::getMass(void)
{
  theMatrix.Zero();

  return theMatrix;
}



void
RCFTLMBeamColumn3D::zeroLoad(void)
{
	
}

int
RCFTLMBeamColumn3D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  return 0;
}


int
RCFTLMBeamColumn3D::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}

const Vector &
RCFTLMBeamColumn3D::getResistingForceIncInertia()
{
  // Compute the current resisting force
  theVector = this->getResistingForce();

  // add the damping forces if rayleigh damping
  if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
     theVector += this->getRayleighDampingForces();

  return theVector;
}



bool
RCFTLMBeamColumn3D::isSubdomain(void)
{
    return false;
}


void
RCFTLMBeamColumn3D::Print(OPS_Stream &s, int flag)
{
   if (flag == 1)
   {
      s << "\nElement: " << this->getTag() << " Type: RCFTLMBeamColumn3D ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << numSections;
      s << "\tMass density: " << rho;
      for (int i = 0; i < numSections; i++)
         s << "\nSection "<<i<<" :" << *sections[i];
    }
   else
   {
      s << "\nElement: " << this->getTag() << " Type: RCFTLMBeamColumn3D ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << numSections;
      s << "\tMass density: " << rho << endln;
   }
}


OPS_Stream &operator<<(OPS_Stream &s, RCFTLMBeamColumn3D &E)
{
    E.Print(s);
    return s;
}


int
RCFTLMBeamColumn3D::displaySelf(Renderer &theViewer, int displayMode, float fact)
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
RCFTLMBeamColumn3D::setResponse(char **argv, int argc, Information &eleInformation)
{
    // global force -
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
	|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
	return new ElementResponse(this, 1, theVector);
    else
      return 0;
}

int
RCFTLMBeamColumn3D::getResponse(int responseID, Information &eleInfo)
{
    switch (responseID) {
      case 1:  // global forces
         return eleInfo.setVector(this->getResistingForce());
      default:
         return -1;
  }
}

int
RCFTLMBeamColumn3D::setParameter (const char **argv, int argc, Parameter &param)
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
		cerr << "RCFTLMBeamColumn3D::setParameter() - could not set parameter. " << endl;
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
RCFTLMBeamColumn3D::updateParameter (int parameterID, Information &info)
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
		cerr << "RCFTLMBeamColumn3D::updateParameter() - could not update parameter. " << endl;
		return ok;
	}
	else {
		return ok;
	}
     }
     else {
	cerr << "RCFTLMBeamColumn3D::updateParameter() - could not update parameter. " << endl;
	return -1;
     }
}

void
RCFTLMBeamColumn3D::setSectionPointers(int numSec, RCFTAggregator **secPtrs)
{
  if (numSec > maxNumSections) {
    opserr << "Error: RCFTLMBeamColumn3D::setSectionPointers -- max number of sections exceeded";
  }

  numSections = numSec;

  if (secPtrs == 0) {
    opserr << "Error: RCFTLMBeamColumn3D::setSectionPointers -- invalid section pointer";
  }

  sections = new RCFTAggregator *[numSections];
  if (sections == 0) {
    opserr << "Error: RCFTLMBeamColumn3D::setSectionPointers -- could not allocate section pointers";
  }

  for (int i = 0; i < numSections; i++) {

    if (secPtrs[i] == 0) {
      opserr << "Error: RCFTLMBeamColumn3D::setSectionPointers -- null section pointer " << i << endln;
    }

    sections[i] = (RCFTAggregator*) secPtrs[i]->getCopy();

    if (sections[i] == 0) {
      opserr << "Error: RCFTLMBeamColumn3D::setSectionPointers -- could not create copy of section " << i << endln;
    }

  }

  // allocate section flexibility matrices and section deformation vectors
  fsa  = new Matrix [numSections];
  if (fsa == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate fs array";
  }

  ksa  = new Matrix [numSections];
  if (ksa == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate ks array";
  }
	  
  fs  = new Matrix [numSections];
  if (fs == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate fs array";
  }

  ks  = new Matrix [numSections];
  if (ks == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate ks array";
  }

  nldhat  = new Matrix [numSections];
  if (nldhat == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate nldhat array";
  }

  nldhatT  = new Matrix [numSections];
  if (nldhatT == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate nldhatT array";
  }

  dhat  = new Vector [numSections];
  if (dhat == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate dhat array";
  }
      
  duhat  = new Vector [numSections];
  if (duhat == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate duhat array";
  }

  sduhat  = new Vector [numSections];
  if (sduhat == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate sduhat array";
  }

  nd1  = new Matrix [numSections];
  if (nd1 == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate nd1 array";
  }

  nd2  = new Matrix [numSections];
  if (nd2 == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate nd2 array";
  }

  nd1T  = new Matrix [numSections];
  if (nd1T == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate nd1T array";
  }

  nd2T  = new Matrix [numSections];
  if (nd2T == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate nd1T array";
  }
 
  nd1Tf  = new Matrix [numSections];
  if (nd1Tf == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate nd1Tf array";
  }

  nd1Tfnd1 = new Matrix [numSections];
  if (nd1Tfnd1 == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate nd1Tfnd1 array";
  }

  nd1Tfnd2  = new Matrix [numSections];
  if (nd1Tfnd2 == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate nd1Tfnd2 array";
  }

  DQ  = new Vector [numSections];
  if (DQ == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate DQ array";
  }

  DSQ  = new Vector [numSections];
  if (DSQ == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate DSQ array";
  }

  DSQa  = new Vector [numSections];
  if (DSQa == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate DSQa array";
  }

  CDSQa  = new Vector [numSections];
  if (CDSQa == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate CDSQa array";
  }

  gd_delta  = new Vector [numSections];
  if (gd_delta == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate gd_delta array";
  }

  str_f4  = new Matrix [numSections];
  if (str_f4 == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate str_f4 array";
  }

  f4  = new Vector [numSections];
  if (f4 == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate f4 array";
  }

  sf4  = new Vector [numSections];
  if (sf4 == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate sf4 array";
  }

  str_f4inv  = new Matrix [numSections];
  if (str_f4inv == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate str_f4inv array";
  }

  d4  = new Vector [numSections];
  if (d4 == 0) {
    opserr << "RCFTLMBeamColumn3D::setSectionPointers -- failed to allocate d4 array";
  }
  
}

Vector
RCFTLMBeamColumn3D::getLocalIncrDeltaDisp(void){
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

const Vector& 
RCFTLMBeamColumn3D::getBasicIncrDisp(void){
   return ub;									    
}	


Vector
RCFTLMBeamColumn3D::getBasicIncrDeltaDisp(void){
   const Vector &disp1 = theNodes[0]->getIncrDeltaDisp();
   const Vector &disp2 = theNodes[1]->getIncrDeltaDisp();

   double ug[18];
   for (int i = 0; i < 9; i++) {
       ug[i]   = disp1(i);
       ug[i+9] = disp2(i);
   }

   double L = crdTransf->getInitialLength();

   double oneOverL = 1.0/L;

   Vector dub(12);

   double ul[18];

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

   ul[12] = R[0][0]*ug[12] + R[0][1]*ug[13] + R[0][2]*ug[14];
   ul[13] = R[1][0]*ug[12] + R[1][1]*ug[13] + R[1][2]*ug[14];
   ul[14] = R[2][0]*ug[12] + R[2][1]*ug[13] + R[2][2]*ug[14];

   ul[15] = R[0][0]*ug[15] + R[0][1]*ug[16] + R[0][2]*ug[17];
   ul[16] = R[1][0]*ug[15] + R[1][1]*ug[16] + R[1][2]*ug[17];
   ul[17] = R[2][0]*ug[15] + R[2][1]*ug[16] + R[2][2]*ug[17];

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

   double d9_0, d10_1, d11_2, d15_6, d16_7, d17_8;

   d9_0   = ul[9] - ul[0];
   d10_1  = ul[10] - ul[1];
   d11_2  = ul[11] - ul[2];
   d15_6  = ul[15] - ul[6];
   d16_7  = ul[16] - ul[7];
   d17_8  = ul[17] - ul[8];

   /************************************************************************/
   /* CALCULATE THE ELONGATION OF THE STEEL AND THE CONCRETE               */
   /* STEEL                                                                */
   /************************************************************************/
   /* STEEL */
   dub(5) = d9_0; 

   /* CONCRETE */
   dub(0) = d15_6; 

   /************************************************************************/
   /* CALCULATE THE RIGID BODY ROTATION OF THE ELEMENT, ASSUMING THAT THE  */
   /* MOVEMENT OF THE STEEL DESCRIBES THE CONCRETE ALSO.                   */
   /************************************************************************/
   double theta_rigid_z, theta_rigid_y;

   theta_rigid_z = ( ul[10] - ul[1] ) / L;

   theta_rigid_y = ( ul[11] - ul[2] ) / L;

   /************************************************************************/
   /* CALCULATE THE NATURAL ROTATIONS OF THE ELEMENT                       */
   /* STEEL ROTATIONS                                                      */
   /************************************************************************/
   dub(6) = ul[5] - theta_rigid_z;
   dub(7) = - ul[4] - theta_rigid_y;
   dub(8) = ul[14] - theta_rigid_z;
   dub(9) = - ul[13] - theta_rigid_y;

   /************************************************************************/
   /* CONCRETE ROTATIONS                                                   */
   /************************************************************************/
   dub(1) = dub(6);
   dub(2) = dub(7);
   dub(3) = dub(8);
   dub(4) = dub(9);

   /************************************************************************/
   /* COMPUTE THE NATURAL DEFORMATION OF THE MIDPOINT DEGREES OF FREEDOM.  */
   /* THIS IS THE ELONGATION OF THE MEMBER BETWEEN THE MIDPOINT AND THE    */
   /* i-end OF THE MEMBER.  THIS IS ACCOMPLISHED BY UNDOING THE STATIC     */
   /* CONDENSATION THAT WAS DONE TO THE NATURAL STIFFNESS MATRIX IN THE    */
   /* FUNCTION b_cft_tangent_k.c                                           */
   /*                                                                      */
   /* THE NATURAL DOFS ARE:                                                */
   /*      u_nat[12]:  CONCRETE ELONGATION BETWEEN MIDPOINT NODE AND       */
   /*                      i-END OF THE MEMBER                             */
   /*      u_nat[13]:  STEEL ELONGATION BETWEEN MIDPOINT NODE AND          */
   /*                      i-END OF THE MEMBER                             */
   /************************************************************************/
   double u_nat[14];

   u_nat[0] = ul[6] - ul[0];
   u_nat[1] = dub(0);
   u_nat[2] = dub(1);
   u_nat[3] = dub(2);
   u_nat[4] = dub(3);
   u_nat[5] = dub(4);
   u_nat[6] = dub(5);
   u_nat[7] = dub(6);
   u_nat[8] = dub(7);
   u_nat[9] = dub(8);
   u_nat[10] = dub(9);
   u_nat[11] = ul[12] - ul[3];

   u_nat[12] = 0.0;
   u_nat[13] = 0.0;

   double temp_sr[2];

   temp_sr[0] = 0.0;
   temp_sr[1] = 0.0;

   for ( int ctr1 = 0; ctr1 < 12; ctr1++ ){
          temp_sr[0] +=  sr(0,ctr1) * u_nat[ctr1];
          temp_sr[1] +=  sr(1,ctr1) * u_nat[ctr1];
   }
   for ( int ctr1 = 0; ctr1 < 2; ctr1++ ){
          u_nat[12] -= ss(0,ctr1) * temp_sr[ctr1];
          u_nat[13] -= ss(1,ctr1) * temp_sr[ctr1];
   }

   dub(10) = u_nat[12];
   dub(11) = u_nat[13];

   return dub;									    
}	

Vector 
RCFTLMBeamColumn3D::getd_hat(int sec, const Vector &v)
{
   double L = crdTransf->getInitialLength();

   double oneOverL  = 1.0/L;

   double xi[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double wt[maxNumSections];
   beamIntegr->getSectionWeights(numSections, L, wt);

   double temp_x, temp_C, temp_D, temp_E, temp_F;

   Vector D_hat(6);

   temp_x = L * xi[sec];
   temp_C = - 1 / L + 4 * temp_x / ( L * L );
   temp_D = - 8 * temp_x / ( L * L ) + 4 / L;
   temp_E = - 4 / L + 6 * temp_x / ( L * L );
   temp_F = - 2 / L + 6 * temp_x / ( L * L );

   D_hat(0) = temp_C * v(0) + temp_D * v(10);

   D_hat(1) = temp_E * v(1) + temp_F * v(3);

   D_hat(2) = temp_E * v(2) + temp_F * v(4);

   D_hat(3) = temp_C * v(5) + temp_D * v(11);

   D_hat(4) = temp_E * v(6) + temp_F * v(8);

   D_hat(5) = temp_E * v(7) + temp_F * v(9);

   return D_hat;     
}

Matrix
RCFTLMBeamColumn3D::getNld_hat(int sec, const Vector &v) 
{
   int i,j;

   double L = crdTransf->getInitialLength();

   double oneOverL  = 1.0/L;

   double xi[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double wt[maxNumSections];
   beamIntegr->getSectionWeights(numSections, L, wt);

   double temp_x, temp_C, temp_D, temp_E, temp_F;

   temp_x = L * xi[sec];
   temp_C = - 1 / L + 4 * temp_x / ( L * L );
   temp_D = - 8 * temp_x / ( L * L ) + 4 / L;
   temp_E = - 4 / L + 6 * temp_x / ( L * L );
   temp_F =  - 2 / L + 6 * temp_x / ( L * L );     

   Matrix Nld_hat(6,12);
    
   for(i = 0; i < 6; i++){
	for(j = 0; j < 12; j++){
		Nld_hat(i,j) = 0.0;
	}
   }

   Nld_hat(0,0)  = temp_C;
   Nld_hat(0,10) = temp_D;
   Nld_hat(1,1)  = temp_E;
   Nld_hat(1,3)  = temp_F;
   Nld_hat(2,2)  = temp_E;
   Nld_hat(2,4)  = temp_F;
   Nld_hat(3,5)  = temp_C;
   Nld_hat(3,11) = temp_D;
   Nld_hat(4,6)  = temp_E;
   Nld_hat(4,8)  = temp_F;
   Nld_hat(5,7)  = temp_E;
   Nld_hat(5,9)  = temp_F;

   return Nld_hat;
}

Matrix
RCFTLMBeamColumn3D::getNd2(int sec)
{
   double L = crdTransf->getInitialLength(); 
   double oneOverL  = 1.0/L;

   double xi[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double wt[maxNumSections];
   beamIntegr->getSectionWeights(numSections, L, wt);

   double temp_x, temp_A, temp_B, temp_C, temp_D;

   temp_x = L * xi[sec];

   Matrix Nd2(6,12);

   for(int i = 0; i < 6; i++){
       for(int j = 0; j < 12; j++){
           Nd2(i,j) = 0.0;
       }
   }

   temp_A = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) );
   temp_B = L * ( -pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) );

   Nd2(1,1) = Tfnat2(0) * temp_A;
   Nd2(1,3) = Tfnat2(0) * temp_B;
   Nd2(2,2) = Tfnat2(0) * temp_A;
   Nd2(2,4) = Tfnat2(0) * temp_B;
   Nd2(4,6) = Tfnat2(6) * temp_A;
   Nd2(4,8) = Tfnat2(6) * temp_B;
   Nd2(5,7) = Tfnat2(6) * temp_A;
   Nd2(5,9) = Tfnat2(6) * temp_B;

   return Nd2;
}

Matrix
RCFTLMBeamColumn3D::getNd1(int sec, const Vector &v)
{
   double L = crdTransf->getInitialLength();

   double oneOverL  = 1.0/L;

   double xi[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double wt[maxNumSections];
   beamIntegr->getSectionWeights(numSections, L, wt);

   double temp_x, temp_A, temp_B, temp_C, temp_D;

   temp_x = L * xi[sec];

   Matrix Nd1(6,12);
    
   for(int i = 0; i < 6; i++){
      for(int j = 0; j < 12; j++){
          Nd1(i,j) = 0.0;
      }
   }
		
   Nd1(0,0)   = temp_x / L - 1.0;
   Nd1(0,1)   = temp_x / L;
   Nd1(1,2)   = temp_x / L - 1.0;
   Nd1(1,4)   = temp_x / L;
   Nd1(2,3)   = temp_x / L - 1.0;
   Nd1(2,5)   = temp_x / L;
   Nd1(3,6)   = temp_x / L - 1.0;
   Nd1(3,7)   = temp_x / L;
   Nd1(4,8)   = temp_x / L - 1.0;
   Nd1(4,10)  = temp_x / L;
   Nd1(5,9)   = temp_x / L - 1.0;
   Nd1(5,11)  = temp_x / L;

   return Nd1;
}

void RCFTLMBeamColumn3D::calcDeformedLength(void)
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


double RCFTLMBeamColumn3D::getDeformedLength(void)
{
   return deflength; 
}

int
RCFTLMBeamColumn3D::sendSelf(int commitTag, Channel &theChannel)
{  
   return 0;
}
    
int
RCFTLMBeamColumn3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
   return 0;
}






    
									    

