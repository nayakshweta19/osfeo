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

// $Revision: 1.3 $
// $Date: 2009-01-07 23:09:33 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/element/RCFTMMBeamColumn3D.cpp,v $



#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <iomanip>

#include <Information.h>
#include <Parameter.h>
#include <RCFTMMBeamColumn3D.h>
#include <MatrixUtil.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <math.h>
#include <OPS_Globals.h>
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

Matrix RCFTMMBeamColumn3D::theMatrix(18,18);
Vector RCFTMMBeamColumn3D::theVector(18);
double RCFTMMBeamColumn3D::workArea[400];

int RCFTMMBeamColumn3D::recoveryflag = 0;

//GaussQuadRule1d01 RCFTMMBeamColumn3D::quadRule;

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
// CONSTRUCTOR FOR PARALLEL PROCESSING

RCFTMMBeamColumn3D::RCFTMMBeamColumn3D():
Element(0,ELE_TAG_RCFTMMBeamColumn3D), connectedExternalNodes(2),
beamIntegr(0), numSections(0), sections(0), crdTransf(0),
rho(0.0), 
initialFlag(0),cnvg(0),
kv(12,12), kvInit(12,12), Se(NEBD), Sg(NEGD), Sglobal(NEGD), CSglobal(NEGD), Sgb(NEGD), kvcommit(12,12), Secommit(NEBD),
fs(0), ks(0), fsa(0), ksa(0), nldhat(0), nldhatT(0), nldhatsc(0), nldhatscT(0), nd1(0), nd2(0), dhat(0), nd1T(0), nd2T(0), nd1Tf(0), 
nd1Tfnd1(0), nd1Tfnd2(0), duhat(0), sduhat(0), duhatcommit(0), gd_delta(0), DQ(0), DSQ(0), DSQa(0), CDSQa(0), 
Ki(0), XAxis(3), YAxis(3), ZAxis(3), V(12), V2(12), G(12,13), G2(13,13), G2T(13,13), Gsc(13,13), H(12,12),
H12(12,13), H22(13,13), Ksc(13,13), Kk(13,13), Hinv(12,12), fnat2(12), Cfnat2(12), Tfnat2(12), fint2(13), Cfint2(13), Md(12,13), sr(2,12), Csr(2,12), Css(2,2), ss(2,2), fk(24), Cfk(24), fk_incr(24), ub(13), ubcommit(13), ubf(12), s(12), Tagg(0), ugti(9), ugtj(9), T(12), S(12), df_i(18), f_i(18), Cf_i(18), itr(0), str_f4(0), f4(0), sf4(0), str_f4inv(0), d4(0), intgr(0), Tfint2(13), Flag(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;

  Sg.Zero();
  Sgb.Zero();
  Sglobal.Zero();
  kv.Zero();
  kvInit.Zero();
  sr.Zero();
  ss.Zero();
  Csr.Zero();
  Css.Zero();
  V.Zero();
  T.Zero();
  S.Zero();
  V2.Zero();
  G.Zero();
  G2.Zero();
  G2T.Zero(); 
  Gsc.Zero();
  H.Zero();
  H22.Zero();
  Hinv.Zero();
  Md.Zero();
  Kk.Zero();
  Ksc.Zero();
  Kg.Zero();
  fk.Zero();
  fk_incr.Zero();
  fint2.Zero();
  Cfint2.Zero();
  fnat2.Zero();
  Cfnat2.Zero();
  Tfnat2.Zero();
  ub.Zero();
  ubcommit.Zero();
  ubf.Zero();
  s.Zero(); 
  ugti.Zero();
  ugtj.Zero();
  df_i.Zero();
  f_i.Zero();
  Cf_i.Zero();
  Tfint2.Zero();
  deflength = 0.0;

}


// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points.
// allocates the necessary space needed by each object
RCFTMMBeamColumn3D::RCFTMMBeamColumn3D (int tag, int nodeI, int nodeJ,
                                      int numSec, RCFTAggregator **sec,
                                      BeamIntegration &bi,
                                      CrdTransf &coordTransf, double massDensPerUnitLength):
Element(tag,ELE_TAG_RCFTMMBeamColumn3D), connectedExternalNodes(2),
beamIntegr(0), numSections(0), sections(0), crdTransf(0),
rho(massDensPerUnitLength),
initialFlag(0), cnvg(0),
kv(12,12), kvInit(12,12), Se(NEBD), Sg(NEGD), Sglobal(NEGD), Sgb(NEGD), kvcommit(12,12), Secommit(NEBD),
fs(0), ks(0), fsa(0), ksa(0), nldhat(0), nldhatT(0), nldhatsc(0), nldhatscT(0), nd1(0), nd2(0), dhat(0), nd1T(0), nd2T(0), nd1Tf(0),
nd1Tfnd1(0), nd1Tfnd2(0), duhat(0), sduhat(0), duhatcommit(0), gd_delta(0), DQ(0), DSQ(0), CDSQa(0),  
Ki(0), XAxis(3), YAxis(3), ZAxis(3), V(12), V2(13), G(12,13), G2(13,13), G2T(13,13), Gsc(13,13), H(12,12),
H12(12,13), H22(13,13), Kk(13,13), Ksc(13,13), Kg(13,13), Hinv(12,12), fnat2(12), Cfnat2(12), Tfnat2(12), fint2(13), Cfint2(13), Md(12,13), sr(2,12), Csr(2,12), ss(2,2), Css(2,2), fk(24), fk_incr(24), ub(13), ubcommit(13), ubf(12), s(12), Tagg(tag), ugti(9), ugtj(9), T(12), S(12), df_i(18), f_i(18), Cf_i(18), itr(0), str_f4(0), sf4(0), f4(0), str_f4inv(0), d4(0), intgr(0), Tfint2(13), Flag(0)
{

   theNodes[0] = 0;
   theNodes[1] = 0;

   connectedExternalNodes(0) = nodeI;
   connectedExternalNodes(1) = nodeJ;

   // get copy of the sections

   beamIntegr = bi.getCopy();
   if(beamIntegr == 0) {
     opserr<<"Error: RCFTMMBeamColumn3D::RCFTMMBeamColumn3D: could not create copy of beam integration object" << endln;
     exit(-1);
   }

   // get copy of the transformation object
   
   crdTransf = coordTransf.getCopy3d();

   //deflength = crdTransf->getInitialLength();
   
   if (crdTransf == 0) {
      opserr << "Error: RCFTMMBeamColumn3D::RCFTMMBeamColumn3D: could not create copy of coordinate transformation object" << endln;
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

   CR[0][0]   = 0.0; CR[0][1] = 0.0; CR[0][2] = 0.0;
   CR[1][0]   = 0.0; CR[1][1] = 0.0; CR[1][2] = 0.0;
   CR[2][0]   = 0.0; CR[2][1] = 0.0; CR[2][2] = 0.0;

   ss.Zero();
   sr.Zero();

   Css.Zero();
   Csr.Zero();
   
   kv.Zero();
   kvInit.Zero();
   V.Zero();
   T.Zero();
   S.Zero();
   Sgb.Zero();
   Sglobal.Zero();
   V2.Zero();
   G.Zero();
   G2.Zero();
   G2T.Zero();
   Gsc.Zero();
   H.Zero();
   H12.Zero();
   H22.Zero();
   Kk.Zero();
   Ksc.Zero();
   Kg.Zero();
   Hinv.Zero();
   Md.Zero();
   fk.Zero();
   fk_incr.Zero();
   ub.Zero();
   ubcommit.Zero();
   ubf.Zero();
   s.Zero();
   ugti.Zero();
   ugtj.Zero();
   df_i.Zero();
   f_i.Zero();
   Cf_i.Zero();
   Tfint2.Zero(); 
   Cfint2.Zero();
}

// ~RCFT_BeamColumn3D():
// 	destructor
//  delete must be invoked on any objects created by the object
RCFTMMBeamColumn3D::~RCFTMMBeamColumn3D()
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

   if(nldhatsc != 0)
     delete [] nldhatsc;

   if(nldhatscT != 0)
     delete [] nldhatscT;

   if(nd2 != 0)
     delete [] nd2;

   if(nd2T != 0)
     delete [] nd2T;

   if(nd1T != 0)
     delete [] nd1T;

   if(nd1Tf != 0)
     delete [] nd1Tf;

   if(nd2Tf != 0)
     delete [] nd2Tf;

   if(nd1Tfnd1 != 0)
     delete [] nd1Tfnd1;

   if(nd1Tfnd2 != 0)
     delete [] nd1Tfnd2;
 
   if(nd2Tfnd2 != 0)
     delete [] nd2Tfnd2;

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

   if(duhatcommit != 0)
     delete [] duhatcommit;

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
RCFTMMBeamColumn3D::getNumExternalNodes(void) const
{
   return 2;
}


const ID &
RCFTMMBeamColumn3D::getExternalNodes(void)
{
   return connectedExternalNodes;
}


Node **
RCFTMMBeamColumn3D::getNodePtrs(void)
{
   return theNodes;
}

int
RCFTMMBeamColumn3D::getNumDOF(void)
{
   return NEGD;
}


void
RCFTMMBeamColumn3D::setDomain(Domain *theDomain)
{
   // check Domain is not null - invoked when object removed from a domain
 
   if (theDomain == 0)
   {
      theNodes[0] = 0;
      theNodes[1] = 0;

      opserr << "RCFTMMBeamColumn3D::setDomain:  theDomain = 0 ";
      exit(0);
   }

   // get pointers to the nodes

   int Nd1 = connectedExternalNodes(0);
   int Nd2 = connectedExternalNodes(1);

   theNodes[0] = theDomain->getNode(Nd1);
   theNodes[1] = theDomain->getNode(Nd2);

   if (theNodes[0] == 0)
   {
      opserr << "RCFTMMBeamColumn3D::setDomain: Nd1: ";
      opserr << Nd1 << "does not exist in model\n";
      exit(0);
   }

   if (theNodes[1] == 0)
   {
      opserr << "RCFTMMBeamColumn3D::setDomain: Nd2: ";
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
      opserr << "RCFTMMBeamColumn3D::setDomain(): Nd2 or Nd1 incorrect dof ";
      exit(0);
   }

   // initialize the transformation
   if (crdTransf->initialize(theNodes[0], theNodes[1]))
   {
      opserr << "RCFTMMBeamColumn3D::setDomain(): Error initializing coordinate transformation";
      exit(0);
   }

   deflength = crdTransf->getInitialLength();

   Li = crdTransf->getInitialLength();

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
      opserr << "RCFTMMBeamColumn3D::setDomain(): Zero element length:" << this->getTag();
      exit(0);
   }

   if (initialFlag == 0)
     this->initializeSectionHistoryVariables();

}



int
RCFTMMBeamColumn3D::commitState()
{
   int err = 0;
   int i = 0;
   int j = 0;

   cnvg = 0;

#ifdef COMPOSITE_DEBUG
   ofstream crv11;
   crv11.open("crv11.dat",ios::app);
   
   ofstream crv12;
   crv12.open("crv12.dat",ios::app);
   
   ofstream crv13;
   crv13.open("crv13.dat",ios::app);
   
   ofstream crv21;
   crv21.open("crv21.dat",ios::app);
   
   ofstream crv22;
   crv22.open("crv22.dat",ios::app);
    
   ofstream crv23;
   crv23.open("crv23.dat",ios::app);
   
   ofstream crv31;
   crv31.open("crv31.dat",ios::app);
   
   ofstream crv32;
   crv32.open("crv32.dat",ios::app);
   
   ofstream crv33;
   crv33.open("crv33.dat",ios::app);
#endif

   // call element commitState to do any base class stuff
   if ((err = this->Element::commitState()) != 0) {
     opserr << "RCFTMMBeamColumn3D::commitState () - failed in base class";
     return err;
   }

   do {
      err = sections[i++]->commitState();
   } while (err == 0 && i < numSections);

   if (err)
      return err;

   // commit the transformation between coord. systems
   if ((err = crdTransf->commitState()) != 0)
      return err;

   // commit the element variables state
   kvcommit = kv;
   Csr = sr; 
   Css = ss;
   Cf_i = f_i;
   Secommit = Se;
   ubcommit = ub;
   Cfnat2 = fnat2;
   Cfint2 = fint2;
   fnat2.Zero();
   fint2.Zero();
   ub.Zero();
   ubf.Zero();
   s.Zero(); 
   CSglobal = Sglobal;
   Cfk = fk;
   CR[0][0]   = R[0][0];
   CR[0][1]   = R[0][1];
   CR[0][2]   = R[0][2];

   CR[1][0]   = R[1][0];
   CR[1][1]   = R[1][1];
   CR[1][2]   = R[1][2];

   CR[2][0]   = R[2][0];
   CR[2][1]   = R[2][1];
   CR[2][2]   = R[2][2];

   //flag = 0;
   
   for( i = 0; i < numSections; i++){
        sduhat[i] = sduhat[i] + duhat[i];
        duhatcommit[i] = duhat[i];
	    duhat[i].Zero();
        dhat[i].Zero();
	    DSQ[i].Zero();
        DQ[i].Zero();
        f4[i].Zero();
        sf4[i].Zero();
        d4[i].Zero(); 
	    CDSQa[i] = DSQa[i];
#ifdef COMPOSITE_DEBUG
		if( ( Tagg == 1 ) && ( i == 0 ) ){
           crv11<<sduhat[i](4)<<"  "<<DSQa[i](1)+DSQa[i](7)<<"  "<<sduhat[i](5)<<"  "<<DSQa[i](2)+DSQa[i](8)<<endl;
        } 
        if( ( Tagg == 1 ) && ( i == 1 ) ){
           crv12<<sduhat[i](4)<<"  "<<DSQa[i](1)+DSQa[i](7)<<"  "<<sduhat[i](5)<<"  "<<DSQa[i](2)+DSQa[i](8)<<endl;
        }
        if( ( Tagg == 1 ) && ( i == 2 ) ){
           crv13<<sduhat[i](4)<<"  "<<DSQa[i](1)+DSQa[i](7)<<"  "<<sduhat[i](5)<<"  "<<DSQa[i](2)+DSQa[i](8)<<endl;
        }
        if( ( Tagg == 2 ) && ( i == 0 ) ){
           crv21<<sduhat[i](4)<<"  "<<DSQa[i](1)+DSQa[i](7)<<"  "<<sduhat[i](5)<<"  "<<DSQa[i](2)+DSQa[i](8)<<endl;
        }
        if( ( Tagg == 2 ) && ( i == 1 ) ){
           crv22<<sduhat[i](4)<<"  "<<DSQa[i](1)+DSQa[i](7)<<"  "<<sduhat[i](5)<<"  "<<DSQa[i](2)+DSQa[i](8)<<endl;
        }
        if( ( Tagg == 2 ) && ( i == 2 ) ){
           crv23<<sduhat[i](4)<<"  "<<DSQa[i](1)+DSQa[i](7)<<"  "<<sduhat[i](5)<<"  "<<DSQa[i](2)+DSQa[i](8)<<endl;
        }
        if( ( Tagg == 3 ) && ( i == 0 ) ){
           crv31<<sduhat[i](4)<<"  "<<DSQa[i](1)+DSQa[i](7)<<"  "<<sduhat[i](5)<<"  "<<DSQa[i](2)+DSQa[i](8)<<endl;
        }
        if( ( Tagg == 3 ) && ( i == 1 ) ){
           crv32<<sduhat[i](4)<<"  "<<DSQa[i](1)+DSQa[i](7)<<"  "<<sduhat[i](5)<<"  "<<DSQa[i](2)+DSQa[i](8)<<endl;
        }
        if( ( Tagg == 3 ) && ( i == 2 ) ){
           crv33<<sduhat[i](4)<<"  "<<DSQa[i](1)+DSQa[i](7)<<"  "<<sduhat[i](5)<<"  "<<DSQa[i](2)+DSQa[i](8)<<endl;
        }
#endif
   }

   itr = 0;
   return err;
}


int RCFTMMBeamColumn3D::revertToLastCommit()
{
   int err;
   int i = 0;
   int j,k;

   cnvg = 1;

   do {
      err = sections[i]->revertToLastCommit();

      ksa[i] = sections[i]->getSectionTangent();

      invertMatrix(12,ksa[i],fsa[i]);

      for( j = 0; j < 6; j++ ){
        for( k = 0; k< 6; k++ ){
          fs[i](j,k) = fsa[i](j,k);
          ks[i](j,k) = ksa[i](j,k);
        }
      }

      DSQa[i] = CDSQa[i];

      DQ[i].Zero();

      DSQ[i].Zero();

      dhat[i].Zero();

      duhat[i].Zero();

      sf4[i].Zero();

      i++;
   } while (err == 0 && i < numSections);


   if (err)
      return err;

   // revert the transformation to last commit
   if ((err = crdTransf->revertToLastCommit()) != 0)
      return err;

   // revert the element state to last commit
   Se    = Secommit;
   kv    = kvcommit;
   ss    = Css;
   sr    = Csr;
   f_i   = Cf_i;
   fint2.Zero();
   fnat2.Zero();
   ub.Zero();
   fk = Cfk;
   Sglobal = CSglobal;
   R[0][0]   = CR[0][0];
   R[0][1]   = CR[0][1];
   R[0][2]   = CR[0][2];

   R[1][0]   = CR[1][0];
   R[1][1]   = CR[1][1];
   R[1][2]   = CR[1][2];

   R[2][0]   = CR[2][0];
   R[2][1]   = CR[2][1];
   R[2][2]   = CR[2][2];

   initialFlag = 0;
   // this->update();

   return err;
}


int RCFTMMBeamColumn3D::revertToStart()
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
RCFTMMBeamColumn3D::getInitialStiff(void)
{
  if (Ki != 0)
    return *Ki;

  int i, j, k;

  double L = crdTransf->getInitialLength();
  double oneOverL  = 1.0/L;

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wt);

  for ( i = 0; i < numSections; i++ ){
      nldhat[i] = this->getNld_hat(i, ub);
      nldhatsc[i] = this->getNld_hatsc(i, ub);
      nd1[i] = this->getNd1(i, ub);
      for( j = 0; j < 6; j++ ){
          for( k = 0; k < 12; k++ ){
                nd1T[i](k,j) = nd1[i](j,k);
          }
      }
      for( j = 0; j < 6; j++ ){
          for( k = 0; k < 13; k++ ){
                nldhatT[i](k,j) = nldhat[i](j,k);
          }
      }
      for( j = 0; j < 1; j++ ){
          for( k = 0; k < 13; k++ ){
                nldhatscT[i](k,j) = nldhatsc[i](j,k);
          }
      }
      ksa[i] = sections[i]->getInitialTangent();
      invertMatrix(12,ksa[i],fsa[i]);
      for( j = 0; j < 6; j++ ){
        for( k = 0; k< 6; k++ ){
          fs[i](j,k) = fsa[i](j,k);
          ks[i](j,k) = ksa[i](j,k);
        }
      }
  }
 
  G.Zero();

  for( i = 0; i < numSections; i++ ){
       G = G + L * wt[i] * nd1T[i] * nldhat[i];
  }

  Matrix GT(13,12);

  for( i = 0; i < 13; i++ ){
     for( j = 0; j < 12; j++ ){
        GT(i,j) = G(j,i);
     }
  }

  H.Zero();

  for( i = 0; i < numSections; i++ ){
       nd1Tf[i] = nd1T[i] * fs[i];
       nd1Tfnd1[i] = nd1Tf[i] * nd1[i];
       H = H + L * wt[i] * nd1Tfnd1[i];
  }

  invertMatrix(12, H, Hinv);
 
  Ksc.Zero();

  for( i = 0; i < numSections; i++ ){
     Ksc = Ksc + L * wt[i] * nldhatscT[i] * ksa[i](12,12) * 2 * ( ksa[i](7,7) + ksa[i](8,8) ) * nldhatsc[i];
  }

  Matrix K_temp(13,13);

  K_temp = Ksc + GT * Hinv * G;

  /************************************************/
  /* CALCULATE STIFFNESS MATRIX WITH TORSION TERM */
  /************************************************/

  Matrix kt(14,14);

  kt.Zero();

  for( i = 0; i < 11; i++ ){
     for( j = 0; j< 11; j++ ){
         kt(i,j) = K_temp(i,j);
     }
  }

  kt(11,11) =  ksa[0](6,6) / L;

  for( i = 11; i < 13; i++ ){
     for( j = 0; j < 11; j++ ){
        kt(i+1,j) = K_temp(i,j);
     }
  }

  for( i = 0; i < 11; i++ ){
     for( j = 11; j < 13; j++ ){
        kt(i,j+1) = K_temp(i,j);
     }
  }

  for( i = 11; i < 13; i++ ){
     for( j = 11; j < 13; j++ ){
        kt(i+1,j+1) = K_temp(i,j);
     }
  }

  Matrix temp_sr(2,12);
  Matrix temp_rr(12,12);

  Matrix rs(12,2);
  Matrix rr(12,12);

  /************************************************************************/
  /*                                                                      */
  /* STATIC CONDENSATION PROCEDURE - SEE DON WHITES THESIS                */
  /*                                                                      */
  /*    | |rr| {rs} |   { D  }          {F}                               */
  /*    |           |             =                                       */
  /*    | {sr} |ss| |   { Dm }          {0}                               */
  /*                                                                      */
  /* so condensed set of equations can be written as:                     */
  /*                                                                      */
  /* |k_condensed| { D } = { F }                                          */
  /*                                                                      */
  /* where:                                                               */
  /*                                                                      */
  /* |k_condensed| = |rr| - {rs} inv|ss| {sr}                             */
  /*                                                                      */
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
  /* GENERATE TERMS IN OTHER MATRICIES FOR CONDENSATION        */
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
RCFTMMBeamColumn3D::getTangentStiff(void)
{
//  int i;  
  crdTransf->update();  // Will remove once we clean up the corotational 2d transformation -- MHS
  const Matrix &KV = crdTransf->getGlobalStiffMatrix(kv,fk);
  return KV;
}

const Vector &
RCFTMMBeamColumn3D::getResistingForce(void)
{
  return Sglobal;
}

void RCFTMMBeamColumn3D::calcResistingForce(void)
{
  crdTransf->update();
  Vector p0(18);
  p0.Zero();
  Sgb = Sgb + df_i;
  Sg = crdTransf->getGlobalResistingForce(f_i, p0);
  Sglobal = Sg;
}

void
RCFTMMBeamColumn3D::initializeSectionHistoryVariables (void)
{
    for (int i = 0; i < numSections; i++)
    {
	    fs[i]       = Matrix(6,6);
	    ks[i]       = Matrix(6,6);
        fsa[i]      = Matrix(12,12);
        ksa[i]      = Matrix(12,12);
	    dhat[i]     = Vector(6);
	    nldhat[i]   = Matrix(6,13);
        nldhatT[i]  = Matrix(13,6);
        nldhatsc[i] = Matrix(1,13);
        nldhatscT[i]= Matrix(13,1);
        dhat[i]     = Vector(6);
        duhat[i]    = Vector(6);
        duhatcommit[i] = Vector(6);
        sduhat[i]   = Vector(6);
        sduhat[i].Zero();
        nd1[i]      = Matrix(6,12);
        nd2[i]      = Matrix(6,13);
        nd1T[i]     = Matrix(12,6);
        nd2T[i]     = Matrix(13,6);
        nd1Tf[i]    = Matrix(12,6);
        nd2Tf[i]    = Matrix(13,6);
        nd1Tfnd1[i] = Matrix(12,12);
        nd1Tfnd2[i] = Matrix(12,13);
        nd2Tfnd2[i] = Matrix(13,13);
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
int RCFTMMBeamColumn3D::update()
{
  if( cnvg == 1 ){
    cnvg = 0;
    return 0;
  }

  int i,j,k;
 
  int linear = 0;
 
  double tempfint2[13];
 
  itr = itr + 1;

  intgr = 0;

  this->calcDeformedLength();

  double L = getDeformedLength();

  double oneOverL  = 1.0/L;

  if( initialFlag == 2 )
  this->revertToLastCommit();

  const Vector &dubwt = getBasicIncrDeltaDisp();
  
  Vector dub(13);

  Vector dubf(12);

  for( int i = 0; i < 11; i++ ){
      dub(i) = dubwt(i);
  }
  dub(11) = dubwt(12);
  dub(12) = dubwt(13);

  dubf(0) = dub(1);
  dubf(1) = dub(1);
  dubf(2) = dub(2);
  dubf(3) = dub(3);
  dubf(4) = dub(4);
  dubf(5) = dub(5);
  dubf(6) = dub(6);
  dubf(7) = dub(6);
  dubf(8) = dub(7);
  dubf(9) = dub(8);
  dubf(10) = dub(9);
  dubf(11) = dub(10);
   
  ub = ub + dub;

  ubf = ubf + dubf;

  const Vector du = getLocalIncrDeltaDisp();

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wt);
  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
  /********************************************************************************/
  /*                                                                              */
  /* CALCULATE INCREMENTAL STRAIN VECTOR FROM INCREMENTAL NATURAL DEFORMATIONS    */
  /*                                                                              */
  /********************************************************************************/

  //if( recoveryflag == 0 ){

  //Flag = 1;

  //recoveryflag = 1;

  if( recoveryflag == 0 ){

  for ( i = 0; i < numSections; i++ ){
      nldhat[i] = this->getNld_hat(i, ub);
      nldhatsc[i] = this->getNld_hatsc(i, ub);
      dhat[i] = this->getd_hat(i, ub);
      nd1[i] = this->getNd1(i, ub);
      nd2[i] = this->getNd2(i);
      for( j = 0; j < 6; j++ ){
          for( k = 0; k < 12; k++ ){
                nd1T[i](k,j) = nd1[i](j,k);
          }
      }
      for( j = 0; j < 6; j++ ){
          for( k = 0; k < 13; k++ ){
                nldhatT[i](k,j) = nldhat[i](j,k);
                nd2T[i](k,j) = nd2[i](j,k);
          }
      }
      for( j = 0; j < 1; j++ ){
          for( k = 0; k < 13; k++ ){
                nldhatscT[i](k,j) = nldhatsc[i](j,k);
          }
      }
      //ksa[i] = sections[i]->getSectionTangent();

      //invertMatrix(12,ksa[i],fsa[i]);

      for( j = 0; j < 6; j++ ){
        for( k = 0; k< 6; k++ ){
          fs[i](j,k) = fsa[i](j,k);
          ks[i](j,k) = ksa[i](j,k);
        }
      }
  }

  V.Zero();
  
  double temp_x, temp_A, temp_B;
//  double d_d[6];, temp_C, temp_D

  for( i = 0; i < numSections; i++ ){
       V = V + L * wt[i] * nd1T[i] * (dhat[i] - duhat[i] -  ( fs[i] * ( DQ[i] - DSQ[i] ) ) );
  }

  G.Zero();

  for( i = 0; i < numSections; i++ ){
       G = G + L * wt[i] * nd1T[i] * nldhat[i];
  }

  G2.Zero();

  for( i = 0; i < numSections; i++ ){
       G2 = G2 + L * wt[i] * nd2T[i] * nldhat[i];
  }

  for( i = 0; i < 13; i++ ){
     for( j = 0; j < 13; j++ ){
       G2T(i,j) = G2(j,i);
     }
  }

  H.Zero();

  for( i = 0; i < numSections; i++ ){
       nd1Tf[i] = nd1T[i] * fs[i];
       nd1Tfnd1[i] = nd1Tf[i] * nd1[i];
       H = H + L * wt[i] * nd1Tfnd1[i];
  }

  invertMatrix(12, H, Hinv);

  fnat2   = fnat2 + Hinv * V;

  for ( i = 0; i < numSections; i++){
     double slp_strn;

     temp_x = L * xi[i];

     temp_A = - temp_x/L + 2 * temp_x * temp_x / ( L * L );

     temp_B = - 4 * temp_x * temp_x / ( L * L ) + 4 * temp_x / L;

     slp_strn = temp_A * ub(6) + temp_B * ub(12) - temp_A * ub(1) - temp_B * ub(11) - ub(0);

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

     DSQ[i](0) = DSQa[i](0) - CDSQa[i](0);
     DSQ[i](1) = DSQa[i](1) - CDSQa[i](1);
     DSQ[i](2) = DSQa[i](2) - CDSQa[i](2);
     DSQ[i](3) = DSQa[i](6) - CDSQa[i](6);
     DSQ[i](4) = DSQa[i](7) - CDSQa[i](7);
     DSQ[i](5) = DSQa[i](8) - CDSQa[i](8);

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

  V2.Zero();

  for( i = 0; i < numSections; i++ ){
     V2 = V2 + L * wt[i] * nd2T[i] * (dhat[i] - duhat[i]);
  }

  H.Zero();

  for( i = 0; i < numSections; i++ ){
     nd1Tf[i] = nd1T[i] * fs[i];
     nd1Tfnd1[i] = nd1Tf[i] * nd1[i];
     H = H + L * wt[i] * nd1Tfnd1[i];
  }

  invertMatrix(12, H, Hinv);

  H12.Zero();

  for( i = 0; i < numSections; i++ ){
     nd1Tf[i] = nd1T[i] * fs[i];
     nd1Tfnd2[i] = nd1Tf[i] * nd2[i];
     H12 = H12 + L * wt[i] * nd1Tfnd2[i];
  }

  H22.Zero();

  for( i = 0; i < numSections; i++ ){
     nd2Tf[i] = nd2T[i] * fs[i];
     nd2Tfnd2[i] = nd2Tf[i] * nd2[i];
     H22 = H22 + L * wt[i] * nd2Tfnd2[i];
  }

  Ksc.Zero();

  for( i = 0; i < numSections; i++ ){
     Ksc = Ksc + L * wt[i] * nldhatscT[i] * ksa[i](12,12) * 2 * ( ksa[i](7,7) + ksa[i](8,8) ) * nldhatsc[i];
  }

  //FORMATION OF THE GEOMETRIC STIFFNESS MATRIX

  Kg.Zero();

  for( i = 0; i < numSections; i++ ){
     Kg = Kg + L * wt[i] * this->getKg(i);
  }

  Matrix temp_Md(12,13);

  Md.Zero();

  /***************************************************/
  /* CALCULATE THE G and T MATRIX                    */
  /***************************************************/
  for( i = 0; i < numSections; i++ ){
     	temp_x = L * xi[i];
      	temp_A =  ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * L;
      	temp_B =  ( - pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * L;
     	temp_Md.Zero();
      	temp_Md(0,2) = temp_A * ( dhat[i](1) - duhat[i](1) );
      	temp_Md(0,3) = temp_A * ( dhat[i](2) - duhat[i](2) );
      	temp_Md(0,4) = temp_B * ( dhat[i](1) - duhat[i](1) );
      	temp_Md(0,5) = temp_B * ( dhat[i](2) - duhat[i](2) );
      	temp_Md(6,7) = temp_A * ( dhat[i](4) - duhat[i](4) );
      	temp_Md(6,8) = temp_A * ( dhat[i](5) - duhat[i](5) );
      	temp_Md(6,9) = temp_B * ( dhat[i](4) - duhat[i](4) );
      	temp_Md(6,10) = temp_B * ( dhat[i](5) - duhat[i](5) );
        Md = Md + L * wt[i] * temp_Md;
  }

  Matrix temp_Gsc(13,13);

  Gsc.Zero();

  for( i = 0; i < numSections; i++ ){
       Gsc = Gsc + L * wt[i] * nldhatscT[i] * ksa[i](12,12) * 2 * ( ksa[i](7,7) + ksa[i](8,8) ) * nldhatsc[i];
  }

  for( i = 0; i < 13; i++ ){
       tempfint2[i] = fint2(i);
  }

  Matrix GMH(12,13);

  Matrix GMHT(13,12);

  Matrix GT(13,12);

  Matrix GscT(13,13);

  GMH = G + Md - H12;

  for( i = 0; i < 13; i++ ){
     for( j = 0; j < 12; j++ ){
        GMHT(i,j) = GMH(j,i); 
        GT(i,j) = G(j,i);
     }
  } 

  for( i = 0; i < 13; i++ ){
     for( j = 0; j < 13; j++ ){
        GscT(i,j) = Gsc(j,i);
     }
  }

  /***************************************************/
  /* CALCULATE STIFFNESS MATRIX WITHOUT TORSION TERM */
  /***************************************************/

  Matrix K_temp(13,13);

  K_temp = ( Kg + Ksc + G2 + G2T - H22 ) + GMHT * Hinv * GMH;

  //K_temp = ( Kg + Ksc ) + GMHT * Hinv * GMH;

  Matrix kt(14,14);

  kt.Zero();

  for( i = 0; i < 11; i++ ){
     for( j = 0; j< 11; j++ ){
         kt(i,j) = K_temp(i,j);
     }
  }

  kt(11,11) =  ksa[0](6,6) / L;

  for( i = 11; i < 13; i++ ){
     for( j = 0; j < 11; j++ ){
        kt(i+1,j) = K_temp(i,j);
     }
  }

  for( i = 0; i < 11; i++ ){
     for( j = 11; j < 13; j++ ){
        kt(i,j+1) = K_temp(i,j);
     }
  }

  for( i = 11; i < 13; i++ ){
     for( j = 11; j < 13; j++ ){
        kt(i+1,j+1) = K_temp(i,j);
     }
  }

  Matrix temp_sr(2,12);
  Matrix temp_rr(12,12);

  Matrix rs(12,2);
  Matrix rr(12,12);
  /************************************************************************/
  /*                                                                      */
  /* STATIC CONDENSATION PROCEDURE - SEE DON WHITES THESIS                */
  /*                                                                      */
  /*    | |rr| {rs} |   { D  }          {F}                               */
  /*    |           |             =                                       */
  /*    | {sr} |ss| |   { Dm }          {0}                               */
  /*                                                                      */
  /* so condensed set of equations can be written as:                     */
  /*                                                                      */
  /* |k_condensed| { D } = { F }                                          */
  /*                                                                      */
  /* where:                                                               */
  /*                                                                      */
  /* |k_condensed| = |rr| - {rs} inv|ss| {sr}                             */
  /*                                                                      */
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
  /* GENERATE TERMS IN OTHER MATRICIES FOR CONDENSATION        */
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
 
  /**************************************/
  /* CALCULATE THE INTERNAL LOAD VECTOR */
  /**************************************/

  fint2 = GT * fnat2 + GscT * ub + V2 + ( GMHT * Hinv ) * V;

  //if( fint2.Norm()/Cfint2.Norm() < 0.1 | fint2.Norm()/Cfint2.Norm() > 10.0 ){
  //    linear = 1;
  //    cout<<"linear trigger"<<endl;
  //}

  }//if(flag == 0)
  if( ( recoveryflag == 1 ) | ( linear == 1 ) ){

  for ( i = 0; i < numSections; i++ ){
      nldhat[i] = this->getNld_hat(i, ub);
      nldhatsc[i] = this->getNld_hatsc(i, ub);
      dhat[i] = this->getd_hat(i, ub);
      for( j = 0; j < 6; j++ ){
          for( k = 0; k < 13; k++ ){
                nldhatT[i](k,j) = nldhat[i](j,k);
          }
      }
      for( j = 0; j < 1; j++ ){
          for( k = 0; k < 13; k++ ){
                nldhatscT[i](k,j) = nldhatsc[i](j,k);
          }
      }
  }
  for ( i = 0; i < numSections; i++)
  {
      double temp_x,temp_A,temp_B;
      double slp_strn;
//      double slp_force;
      temp_x = L * xi[i];
      temp_A = - temp_x/L + 2 * temp_x * temp_x / ( L * L );
      temp_B = - 4 * temp_x * temp_x / ( L * L ) + 4 * temp_x / L;
      slp_strn = temp_A * ub(6) + temp_B * ub(12) - temp_A * ub(1) - temp_B * ub(11) - ub(0);
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
      DSQa[i] = sections[i]->getStressResultant();
      DSQ[i](0) = DSQa[i](0) - CDSQa[i](0);
      DSQ[i](1) = DSQa[i](1) - CDSQa[i](1);
      DSQ[i](2) = DSQa[i](2) - CDSQa[i](2);
      DSQ[i](3) = DSQa[i](6) - CDSQa[i](6);
      DSQ[i](4) = DSQa[i](7) - CDSQa[i](7);
      DSQ[i](5) = DSQa[i](8) - CDSQa[i](8);
      ksa[i] = sections[i]->getSectionTangent();
      for( j = 0; j < 6; j++ ){
    	  for( k = 0; k< 6; k++ ){
        	ks[i](j,k) = ksa[i](j,k);
       	  }
      }
      double tempgk[9];

      for(j = 0; j < 3; j++){
          for(k = 0; k < 3; k++){
              tempgk[j*3+k]=ks[i](j,k);
          }
      }

  } // for ( i = 0; i < numSections; i++)

  Vector Z(13);

  Z.Zero();

  for( i = 0; i < numSections; i++ ){
      Z = Z + L * wt[i] * nldhatT[i] * DSQ[i];
  }

  Gsc.Zero();

  for( i = 0; i < numSections; i++ ){
      Gsc = Gsc + L * wt[i] * nldhatscT[i] * ksa[i](12,12) * 2 * ( ksa[i](7,7) + ksa[i](8,8) ) * nldhatsc[i];
  }

  for( i = 0; i < 13; i++ ){
       tempfint2[i] = fint2(i);
  }

  Matrix GscT(13,13);

  for( i = 0; i < 13; i++ ){
     for( j = 0; j < 13; j++ ){
        GscT(i,j) = Gsc(j,i);
     }
  }

  Kk.Zero();

  for( i = 0; i < numSections; i++ ){
      Kk = Kk + L * wt[i] * nldhatT[i] * ks[i] * nldhat[i];
  }

  Ksc.Zero();

  for( i = 0; i < numSections; i++ ){
      Ksc = Ksc + L * wt[i] * nldhatscT[i] * ksa[i](12,12) * 2 * ( ksa[i](7,7) + ksa[i](8,8) ) * nldhatsc[i];
  }

  Kg.Zero();
  
  for( i = 0; i < numSections; i++ ){
      Kg = Kg + L * wt[i] * this->getKg(i);
  }
 
  /***************************************************/
  /* CALCULATE STIFFNESS MATRIX WITHOUT TORSION TERM */
  /***************************************************/

  Matrix K_temp(13,13);;

  K_temp = Kg + Ksc + Kk; 

  /************************************************/
  /* CALCULATE STIFFNESS MATRIX WITH TORSION TERM */
  /************************************************/

  Matrix kt(14,14);

  kt.Zero();

  for( i = 0; i < 11; i++ ){
     for( j = 0; j< 11; j++ ){
         kt(i,j) = K_temp(i,j);
     }
  }

  kt(11,11) =  ksa[0](6,6) / L;

  for( i = 11; i < 13; i++ ){
     for( j = 0; j < 11; j++ ){
        kt(i+1,j) = K_temp(i,j);
     }
  }

  for( i = 0; i < 11; i++ ){
     for( j = 11; j < 13; j++ ){
        kt(i,j+1) = K_temp(i,j);
     }
  }

  for( i = 11; i < 13; i++ ){
     for( j = 11; j < 13; j++ ){
        kt(i+1,j+1) = K_temp(i,j);
     }
  }

  Matrix temp_sr(2,12);
  Matrix temp_rr(12,12);

  Matrix rs(12,2);
  Matrix rr(12,12);
  /************************************************************************/
  /*                                                                      */
  /* STATIC CONDENSATION PROCEDURE - SEE DON WHITES THESIS                */
  /*                                                                      */
  /*    | |rr| {rs} |   { D  }          {F}                               */
  /*    |           |             =                                       */
  /*    | {sr} |ss| |   { Dm }          {0}                               */
  /*                                                                      */
  /* so condensed set of equations can be written as:                     */
  /*                                                                      */
  /* |k_condensed| { D } = { F }                                          */
  /*                                                                      */
  /* where:                                                               */
  /*                                                                      */
  /* |k_condensed| = |rr| - {rs} inv|ss| {sr}                             */
  /*                                                                      */
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
  /* GENERATE TERMS IN OTHER MATRICIES FOR CONDENSATION        */
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
  
  /**************************************/
  /* CALCULATE THE INTERNAL LOAD VECTOR */
  /**************************************/
  fint2 = Z + GscT * ub;

  }

  /***********************************************/
  /* CALCULATE INCREMENTAL INTERNAL LOAD VECTOR  */
  /***********************************************/
  Vector dfint(13);

  for( i = 0; i < 13; i++ ){
       dfint(i) = fint2(i) - tempfint2[i];
  }

  /************************************************************************/
  /* NOW UPDATE THE INCREMENTAL ELEMENT FORCES INCLUDING SHEAR            */
  /************************************************************************/
  //df_i(0)  = - dfint(6);
  df_i(0)  = - dfint(6) - dfint(0);
  df_i(1)  = (dfint(7) + dfint(9)) / L;
  df_i(2)  = (dfint(8) + dfint(10)) / L;
  df_i(3)  = du(3) * ksa[0](6,6) / L - du(12) * ksa[0](6,6) / L;
  df_i(4)  = - dfint(3) - dfint(8);
  df_i(5)  = dfint(2) + dfint(7);
  //df_i(6)  = - dfint(1);
  df_i(6)  = - dfint(1) + dfint(0);
  df_i(7)  = (dfint(2) + dfint(4)) / L;
  df_i(8)  = (dfint(3) + dfint(5)) / L;
  df_i(9)  = dfint(6);
  df_i(10) = - df_i(1);
  df_i(11) = - df_i(2);
  df_i(12) = du(3) * ksa[0](6,6) / L - du(12) * ksa[0](6,6) / L;    
  df_i(13) = - dfint(5) - dfint(10);
  df_i(14) = dfint(4) + dfint(9);
  df_i(15) = dfint(1);
  df_i(16) = - df_i(7);
  df_i(17) = - df_i(8);

  f_i(0)  =  f_i(0) + df_i(0);
  f_i(1)  =  f_i(1) + df_i(1);
  f_i(2)  =  f_i(2) + df_i(2);
  f_i(3)  =  f_i(3) + df_i(3);
  f_i(4)  =  f_i(4) + df_i(4);
  f_i(5)  =  f_i(5) + df_i(5);
  f_i(6)  =  f_i(6) + df_i(6);
  f_i(7)  =  f_i(7) + df_i(7);
  f_i(8)  =  f_i(8) + df_i(8);
  f_i(9)  =  f_i(9) + df_i(9);
  f_i(10) =  f_i(10) + df_i(10);
  f_i(11) =  f_i(11) + df_i(11);
  f_i(12) =  f_i(12) + df_i(12);
  f_i(13) =  f_i(13) + df_i(13);
  f_i(14) =  f_i(14) + df_i(14);
  f_i(15) =  f_i(15) + df_i(15);
  f_i(16) =  f_i(16) + df_i(16);
  f_i(17) =  f_i(17) + df_i(17);
  
  fk_incr(0) = dfint(1);
  fk_incr(1) = dfint(6);
  fk_incr(2) = - dfint(3);
  fk_incr(3) = - dfint(2);
  fk_incr(4) = - dfint(8);
  fk_incr(5) = - dfint(7);
  fk_incr(6) = 0.0;
  fk_incr(8) = 0.0;
  fk_incr(7) = 0.0;
  fk_incr(9) = 0.0;
  fk_incr(11)= 0.0;
  fk_incr(10)= 0.0;
  fk(1) += fk_incr(1);
  fk(0) += fk_incr(0);
  fk(4) += fk_incr(4);
  fk(5) += fk_incr(5);
  fk(2) += fk_incr(2);
  fk(3) += fk_incr(3);
  fk(6) = DSQa[0](3);
  fk(8) = DSQa[0](4);
  fk(7) = DSQa[0](5);
  fk(9) = DSQa[0](9);
  fk(11)= DSQa[0](11);
  fk(10)= DSQa[0](10);

  fk_incr(12) = dfint(1);
  fk_incr(13) = dfint(6);
  fk_incr(14) = dfint(5);
  fk_incr(15) = dfint(4);
  fk_incr(16) = dfint(10);
  fk_incr(17) = dfint(9);
  fk_incr(18) = 0.0;
  fk_incr(20) = 0.0;
  fk_incr(19) = 0.0;
  fk_incr(21) = 0.0;
  fk_incr(23) = 0.0;
  fk_incr(22) = 0.0;
  fk(13) += fk_incr(13);
  fk(12) += fk_incr(12);
  fk(16) += fk_incr(16);
  fk(17) += fk_incr(17);
  fk(14) += fk_incr(14);
  fk(15) += fk_incr(15);
  fk(18) = DSQa[numSections-1](3);
  fk(20) = DSQa[numSections-1](5);
  fk(19) = DSQa[numSections-1](4);
  fk(21) = DSQa[numSections-1](9);
  fk(23) = DSQa[numSections-1](11);
  fk(22) = DSQa[numSections-1](10);

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
RCFTMMBeamColumn3D::getMass(void)
{
  theMatrix.Zero();

  return theMatrix;
}



void
RCFTMMBeamColumn3D::zeroLoad(void)
{
	
}

int
RCFTMMBeamColumn3D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  return 0;
}


int
RCFTMMBeamColumn3D::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}

const Vector &
RCFTMMBeamColumn3D::getResistingForceIncInertia()
{
  // Compute the current resisting force
  theVector = this->getResistingForce();

  // add the damping forces if rayleigh damping
  if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
     theVector += this->getRayleighDampingForces();

  return theVector;
}



bool
RCFTMMBeamColumn3D::isSubdomain(void)
{
    return false;
}


void
RCFTMMBeamColumn3D::Print(OPS_Stream &s, int flag)
{
   if (flag == 1)
   {
      s << "\nElement: " << this->getTag() << " Type: RCFTMMBeamColumn3D ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << numSections;
      s << "\tMass density: " << rho;
      for (int i = 0; i < numSections; i++)
         s << "\nSection "<<i<<" :" << *sections[i];
    }
   else
   {
      s << "\nElement: " << this->getTag() << " Type: RCFTMMBeamColumn3D ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << numSections;
      s << "\tMass density: " << rho << endln;
   }
}


OPS_Stream &operator<<(OPS_Stream &s, RCFTMMBeamColumn3D &E)
{
    E.Print(s);
    return s;
}


int
RCFTMMBeamColumn3D::displaySelf(Renderer &theViewer, int displayMode, float fact)
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
RCFTMMBeamColumn3D::setResponse(char **argv, int argc, Information &eleInformation)
{
    // global force -
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
	|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
	return new ElementResponse(this, 1, theVector);
    else
      return 0;
}

int
RCFTMMBeamColumn3D::getResponse(int responseID, Information &eleInfo)
{
    switch (responseID) {
      case 1:  // global forces
         return eleInfo.setVector(this->getResistingForce());
      default:
         return -1;
  }
}

int
RCFTMMBeamColumn3D::setParameter (const char **argv, int argc, Parameter &param)
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
		 return param.addObject(1, this);

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
		cerr << "RCFTMMBeamColumn3D::setParameter() - could not set parameter. " << endl;
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
RCFTMMBeamColumn3D::updateParameter (int parameterID, Information &info)
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
		cerr << "RCFTMMBeamColumn3D::updateParameter() - could not update parameter. " << endl;
		return ok;
	}
	else {
		return ok;
	}
     }
     else {
	cerr << "RCFTMMBeamColumn3D::updateParameter() - could not update parameter. " << endl;
	return -1;
     }
}

void
RCFTMMBeamColumn3D::setSectionPointers(int numSec, RCFTAggregator **secPtrs)
{
  if (numSec > maxNumSections) {
    opserr << "Error: RCFTMMBeamColumn3D::setSectionPointers -- max number of sections exceeded";
  }

  numSections = numSec;

  if (secPtrs == 0) {
    opserr << "Error: RCFTMMBeamColumn3D::setSectionPointers -- invalid section pointer";
  }

  sections = new RCFTAggregator *[numSections];
  if (sections == 0) {
    opserr << "Error: RCFTMMBeamColumn3D::setSectionPointers -- could not allocate section pointers";
  }

  for (int i = 0; i < numSections; i++) {

    if (secPtrs[i] == 0) {
      opserr << "Error: RCFTMMBeamColumn3D::setSectionPointers -- null section pointer " << i << endln;
    }

    sections[i] = (RCFTAggregator*) secPtrs[i]->getCopy();

    if (sections[i] == 0) {
      opserr << "Error: RCFTMMBeamColumn3D::setSectionPointers -- could not create copy of section " << i << endln;
    }

  }

  // allocate section flexibility matrices and section deformation vectors
  fsa  = new Matrix [numSections];
  if (fsa == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate fs array";
  }

  ksa  = new Matrix [numSections];
  if (ksa == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate ks array";
  }
	  
  fs  = new Matrix [numSections];
  if (fs == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate fs array";
  }

  ks  = new Matrix [numSections];
  if (ks == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate ks array";
  }

  nldhat  = new Matrix [numSections];
  if (nldhat == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate nldhat array";
  }

  nldhatT  = new Matrix [numSections];
  if (nldhatT == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate nldhatT array";
  }

  nldhatscT  = new Matrix [numSections];
  if (nldhatscT == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate nldhatscT array";
  }
  
  nldhatsc  = new Matrix [numSections];
  if (nldhatsc == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate nldhatsc array";
  }

  dhat  = new Vector [numSections];
  if (dhat == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate dhat array";
  }
      
  duhat  = new Vector [numSections];
  if (duhat == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate duhat array";
  }

  duhatcommit  = new Vector [numSections];
  if (duhatcommit == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate duhatcommit array";
  }

  sduhat  = new Vector [numSections];
  if (sduhat == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate sduhat array";
  }

  nd1  = new Matrix [numSections];
  if (nd1 == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate nd1 array";
  }

  nd2  = new Matrix [numSections];
  if (nd2 == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate nd2 array";
  }

  nd1T  = new Matrix [numSections];
  if (nd1T == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate nd1T array";
  }

  nd2T  = new Matrix [numSections];
  if (nd2T == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate nd1T array";
  }
 
  nd1Tf  = new Matrix [numSections];
  if (nd1Tf == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate nd1Tf array";
  }

  nd2Tf  = new Matrix [numSections];
  if (nd2Tf == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate nd2Tf array";
  }

  nd1Tfnd1 = new Matrix [numSections];
  if (nd1Tfnd1 == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate nd1Tfnd1 array";
  }

  nd1Tfnd2  = new Matrix [numSections];
  if (nd1Tfnd2 == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate nd1Tfnd2 array";
  }

  nd2Tfnd2  = new Matrix [numSections];
  if (nd2Tfnd2 == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate nd2Tfnd2 array";
  }

  DQ  = new Vector [numSections];
  if (DQ == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate DQ array";
  }

  DSQ  = new Vector [numSections];
  if (DSQ == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate DSQ array";
  }

  DSQa  = new Vector [numSections];
  if (DSQa == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate DSQa array";
  }

  CDSQa  = new Vector [numSections];
  if (CDSQa == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate CDSQa array";
  }

  gd_delta  = new Vector [numSections];
  if (gd_delta == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate gd_delta array";
  }

  str_f4  = new Matrix [numSections];
  if (str_f4 == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate str_f4 array";
  }

  f4  = new Vector [numSections];
  if (f4 == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate f4 array";
  }

  sf4  = new Vector [numSections];
  if (sf4 == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate sf4 array";
  }

  str_f4inv  = new Matrix [numSections];
  if (str_f4inv == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate str_f4inv array";
  }

  d4  = new Vector [numSections];
  if (d4 == 0) {
    opserr << "RCFTMMBeamColumn3D::setSectionPointers -- failed to allocate d4 array";
  }
  
}

Vector
RCFTMMBeamColumn3D::getLocalIncrDeltaDisp(void)
{

    const Vector &disp1 = theNodes[0]->getIncrDeltaDisp();
    const Vector &disp2 = theNodes[1]->getIncrDeltaDisp();

    double ug[18];
    for (int i = 0; i < 9; i++) {
        ug[i]   = disp1(i);
        ug[i+9] = disp2(i);
    }

    double L = getDeformedLength();
    
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
RCFTMMBeamColumn3D::getBasicIncrDisp(void)
{
    return ub;									    
}	


Vector
RCFTMMBeamColumn3D::getBasicIncrDeltaDisp(void)
{
   const Vector &disp1 = theNodes[0]->getIncrDeltaDisp();
   const Vector &disp2 = theNodes[1]->getIncrDeltaDisp();

   double ug[18];
   for (int i = 0; i < 9; i++) {
       ug[i]   = disp1(i);
       ug[i+9] = disp2(i);
   }

   double L = this->getDeformedLength(); 

   double oneOverL = 1.0/L;

   Vector dub(14);

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

   d9_0 = ul[9] - ul[0];
   d10_1 = ul[10] - ul[1];
   d11_2 = ul[11] - ul[2];
   d15_6 = ul[15] - ul[6];
   d16_7 = ul[16] - ul[7];
   d17_8 = ul[17] - ul[8];

   /************************************************************************/
   /* CALCULATE THE ELONGATION OF THE STEEL AND THE CONCRETE               */
   /* STEEL                                                                */
   /************************************************************************/
   /* STEEL */
   dub(0) = ul[6] - ul[0];

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

   theta_rigid_z = ( ul[10] - ul[1] ) / L;

   theta_rigid_y = ( ul[11] - ul[2] ) / L;

   /************************************************************************/
   /* CALCULATE THE NATURAL ROTATIONS OF THE ELEMENT                       */
   /* STEEL ROTATIONS                                                      */
   /************************************************************************/
   dub(7) = ul[5] - theta_rigid_z;
   dub(8) = - ul[4] - theta_rigid_y;
   dub(9) = ul[14] - theta_rigid_z;
   dub(10) = - ul[13] - theta_rigid_y;

   /************************************************************************/
   /* CONCRETE ROTATIONS                                                   */
   /************************************************************************/
   dub(2) = dub(7);
   dub(3) = dub(8);
   dub(4) = dub(9);
   dub(5) = dub(10);

   dub(11) = ul[12] - ul[3];

   /************************************************************************/
   /* COMPUTE THE NATURAL DEFORMATION OF THE MIDPOINT DEGREES OF FREEDOM.  */
   /* THIS IS THE ELONGATION OF THE MEMBER BETWEEN THE MIDPOINT AND THE    */
   /* i-end OF THE MEMBER.  THIS IS ACCOMPLISHED BY UNDOING THE STATIC     */
   /* CONDENSATION THAT WAS DONE TO THE NATURAL STIFFNESS MATRIX IN THE    */
   /* FUNCTION                                           		   */
   /*                                                                      */
   /* THE NATURAL DOFS ARE:                                                */
   /*      u_nat[12]:  CONCRETE ELONGATION BETWEEN MIDPOINT NODE AND       */
   /*                      i-END OF THE MEMBER                             */
   /*      u_nat[13]:  STEEL ELONGATION BETWEEN MIDPOINT NODE AND          */
   /*                      i-END OF THE MEMBER                             */
   /************************************************************************/
   dub(12) = 0.0;
   dub(13) = 0.0;

   double temp_sr[2];

   temp_sr[0] = 0.0;
   temp_sr[1] = 0.0;

   int ctr1;

   for ( ctr1 = 0; ctr1 < 12; ctr1++ ){
          temp_sr[0] +=  sr(0,ctr1) * dub(ctr1);
          temp_sr[1] +=  sr(1,ctr1) * dub(ctr1);
   }
   for ( ctr1 = 0; ctr1 < 2; ctr1++ ){
          dub(12) -= ss(0,ctr1) * temp_sr[ctr1];
          dub(13) -= ss(1,ctr1) * temp_sr[ctr1];
   }

   return dub;			
}	

Vector 
RCFTMMBeamColumn3D::getd_hat(int sec, const Vector &v)
{
    double L = getDeformedLength();
    double oneOverL  = 1.0/L;

    double xi[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, xi);

    double wt[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wt);

    double temp_x, temp_A, temp_B, temp_C, temp_D, temp_E, temp_F;

    Vector D_hat(6);

    temp_x = L * xi[sec];

    temp_A = 1 - 4 * ( temp_x / L ) + 3 * pow ( ( temp_x / L ) , 2 );
    temp_B = - 2 * ( temp_x / L ) + 3 * pow ( ( temp_x / L ) , 2 );
    temp_C = - 1 / L + 4 * temp_x / ( L * L );
    temp_D = - 8 * temp_x / ( L * L ) + 4 / L;
    temp_E = - 4 / L + 6 * temp_x / ( L * L );
    temp_F = - 2 / L + 6 * temp_x / ( L * L );

     D_hat(0) = temp_C * v(1) +
        0.5 * ( temp_C * temp_C * v(1) + temp_C * temp_D * v(11) ) * v(1)+
        0.5 * ( temp_A * temp_A * v(2) + temp_A * temp_B * v(4) ) * v(2) +
        0.5 * ( temp_A * temp_A * v(3) + temp_A * temp_B * v(5) ) * v(3) +
        0.5 * ( temp_A * temp_B * v(2) + temp_B * temp_B * v(4) ) * v(4) +
        0.5 * ( temp_A * temp_B * v(3) + temp_B * temp_B * v(5) ) * v(5) +
        0.5 * ( temp_D * temp_D * v(11) + temp_C * temp_D * v(1) ) * v(11) +
        temp_D * v(11);

    D_hat(1) = temp_E * v(2) + temp_F * v(4);

    D_hat(2) = temp_E * v(3) + temp_F * v(5);

    D_hat(3) = temp_C * v(6) +
        0.5 * ( temp_C * temp_C * v(6) + temp_C * temp_D * v(12) ) * v(1)+
        0.5 * ( temp_A * temp_A * v(7) + temp_A * temp_B * v(9) ) * v(7) +
        0.5 * ( temp_A * temp_A * v(8) + temp_A * temp_B * v(10) ) * v(8) +
        0.5 * ( temp_A * temp_B * v(7) + temp_B * temp_B * v(9) ) * v(9) +
        0.5 * ( temp_A * temp_B * v(8) + temp_B * temp_B * v(10) ) * v(10)+
        0.5 * ( temp_D * temp_D * v(12) + temp_C * temp_D * v(6) ) * v(12) +
        temp_D * v(12);

    D_hat(4) = temp_E * v(7) + temp_F * v(9);

    D_hat(5) = temp_E * v(8) + temp_F * v(10);

    return D_hat;     
}

Matrix
RCFTMMBeamColumn3D::getNld_hat(int sec, const Vector &v) 
{
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

    Matrix Nld_hat(6,13);
    
    for(i = 0; i < 6; i++){
	for(j = 0; j < 13; j++){
		Nld_hat(i,j) = 0.0;
	}
    }

    Nld_hat(0,1)  = temp_C + ( temp_C * temp_C * v(1) + temp_C * temp_D * v(11) );
    Nld_hat(0,2)  = ( temp_A * temp_A * v(2) + temp_A * temp_B * v(4) );
    Nld_hat(0,3)  = ( temp_A * temp_A * v(3) + temp_A * temp_B * v(5) );
    Nld_hat(0,4)  = ( temp_A * temp_B * v(2) + temp_B * temp_B * v(4) );
    Nld_hat(0,5)  = ( temp_A * temp_B * v(3) + temp_B * temp_B * v(5) );
    Nld_hat(0,11) = temp_D + ( temp_D * temp_D * v(11) + temp_C * temp_D * v(1) );
    Nld_hat(1,2)  = temp_E;
    Nld_hat(1,4)  = temp_F;
    Nld_hat(2,3)  = temp_E;
    Nld_hat(2,5)  = temp_F;
    Nld_hat(3,6)  = temp_C + ( temp_C * temp_C * v(6) + temp_C * temp_D * v(12) );
    Nld_hat(3,7)  = ( temp_A * temp_A * v(7) + temp_A * temp_B * v(9) );
    Nld_hat(3,8)  = ( temp_A * temp_A * v(8) + temp_A * temp_B * v(10) );
    Nld_hat(3,9)  = ( temp_A * temp_B * v(7) + temp_B * temp_B * v(9) );
    Nld_hat(3,10)  = ( temp_A * temp_B * v(8) + temp_B * temp_B * v(10) );
    Nld_hat(3,12) = temp_D + ( temp_D * temp_D * v(12) + temp_C * temp_D * v(6) );
    Nld_hat(4,7)  = temp_E;
    Nld_hat(4,9)  = temp_F;
    Nld_hat(5,8)  = temp_E;
    Nld_hat(5,10)  = temp_F;

    return Nld_hat;
}

Matrix
RCFTMMBeamColumn3D::getNld_hatsc(int sec, const Vector &v)
{
    int i,j;

    double L = getDeformedLength();

    double oneOverL  = 1.0/L;

    double xi[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, xi);

    double wt[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wt);

    double temp_x, temp_A, temp_B;

    temp_x = L * xi[sec];
    temp_A = ( temp_x / L ) - 2 * pow ( ( temp_x / L ) , 2 );
    temp_B =  - 4 * ( temp_x / L ) + 4 * pow ( ( temp_x / L ) , 2 );

    Matrix Nld_hatsc(1,13);

    for(i = 0; i < 1; i++){
        for(j = 0; j < 13; j++){
                Nld_hatsc(i,j) = 0.0;
        }
    }

    Nld_hatsc(0,0)   = -1;
    Nld_hatsc(0,1)   = temp_A;
    Nld_hatsc(0,6)   = -temp_A;
    Nld_hatsc(0,11)  = temp_B;
    Nld_hatsc(0,12)  = -temp_B;

    return Nld_hatsc;
}

Matrix
RCFTMMBeamColumn3D::getNd2(int sec)
{
    double L = getDeformedLength();
    double oneOverL  = 1.0/L;

    double xi[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, xi);

    double wt[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wt);

    double temp_x, temp_A, temp_B;
//, temp_C, temp_D
    temp_x = L * xi[sec];

    Matrix Nd2(6,13);

    for(int i = 0; i < 6; i++){
        for(int j = 0; j < 13; j++){
            Nd2(i,j) = 0.0;
        }
    }

    temp_A = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) );
    temp_B = L * ( -pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) );

    double Pc = DSQa[sec](0);
    double Ps = DSQa[sec](6);

    Nd2(1,2) =  Pc * temp_A;
    Nd2(1,4) =  Pc * temp_B;
    Nd2(2,3) =  Pc * temp_A;
    Nd2(2,5) =  Pc * temp_B;
    Nd2(4,7) =  Ps * temp_A;
    Nd2(4,9) =  Ps * temp_B;
    Nd2(5,8) =  Ps * temp_A;
    Nd2(5,10) = Ps * temp_B;

    return Nd2;
}

Matrix
RCFTMMBeamColumn3D::getKg(int sec)
{
    double L = getDeformedLength();
    double oneOverL  = 1.0/L;

    double xi[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, xi);

    double wt[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wt);

    double temp_x, temp_A, temp_B;

    temp_x = L * xi[sec];

    Matrix kg(13,13);

    for(int i = 0; i < 13; i++){
        for(int j = 0; j < 13; j++){
            kg(i,j) = 0.0;
        }
    }

    temp_A = 1 - 4 * temp_x / L + 3 * ( temp_x * temp_x ) / ( L * L );
    temp_B = - 2 * temp_x / L + 3 * ( temp_x * temp_x ) / ( L * L );

    double Pc = DSQa[sec](0);
    double Ps = DSQa[sec](6);

    kg(2,2)   = Pc * temp_A * temp_A;
    kg(2,4)   = Pc * temp_A * temp_B;
    kg(3,3)   = Pc * temp_A * temp_A;
    kg(3,5)   = Pc * temp_A * temp_B;
    kg(4,2)   = Pc * temp_A * temp_B;
    kg(4,4)   = Pc * temp_B * temp_B;
    kg(5,3)   = Pc * temp_A * temp_B;
    kg(5,5)   = Pc * temp_B * temp_B;
    kg(7,7)   = Ps * temp_A * temp_A;
    kg(7,9)   = Ps * temp_A * temp_B;
    kg(8,8)   = Ps * temp_A * temp_A;
    kg(8,10)  = Ps * temp_A * temp_B;
    kg(9,7)   = Ps * temp_A * temp_B;
    kg(9,9)   = Ps * temp_B * temp_B;
    kg(10,8)  = Ps * temp_A * temp_B;
    kg(10,10) = Ps * temp_B * temp_B;  
    
    return kg;
}

Matrix
RCFTMMBeamColumn3D::getNd1(int sec, const Vector &v)
{
    double L = getDeformedLength();
    double oneOverL  = 1.0/L;

    double xi[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, xi);

    double wt[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wt);

    double temp_x, temp_A, temp_B, temp_C, temp_D;

    temp_x = L * xi[sec];

    temp_A = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[2]
            + L * ( - pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[4];

    temp_B = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[3]
            + L * ( - pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[5];

    temp_C = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[7]
            + L * ( - pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[9];

    temp_D = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[8]
            + L * ( - pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) ) * v[10];

    Matrix Nd1(6,12);
    
    for(int i = 0; i < 6; i++){
	for(int j = 0; j < 12; j++){
	    Nd1(i,j) = 0.0;
	}
    }
		
    Nd1(0,0)   = -temp_x / L + 1.0;
    Nd1(0,1)   = temp_x / L;
    Nd1(1,0)   = temp_A;
    Nd1(1,2)   = -temp_x / L + 1.0;
    Nd1(1,4)   = temp_x / L;
    Nd1(2,0)   = temp_B;
    Nd1(2,3)   = -temp_x / L + 1.0;
    Nd1(2,5)   = temp_x / L;
    Nd1(3,6)   = -temp_x / L + 1.0;
    Nd1(3,7)   = temp_x / L;
    Nd1(4,6)   = temp_C;
    Nd1(4,8)   = -temp_x / L + 1.0;
    Nd1(4,10)  = temp_x / L;
    Nd1(5,6)   = temp_D;
    Nd1(5,9)   = -temp_x / L + 1.0;
    Nd1(5,11)  = temp_x / L;

    return Nd1;
}

void RCFTMMBeamColumn3D::calcDeformedLength(void)
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

double RCFTMMBeamColumn3D::getminEigenValue(int n, double *b)
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

double RCFTMMBeamColumn3D::getDeformedLength(void)
{
   return deflength; 
}

int
RCFTMMBeamColumn3D::sendSelf(int commitTag, Channel &theChannel)
{  
  return 0;
}
    
int
RCFTMMBeamColumn3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}






    
									    

