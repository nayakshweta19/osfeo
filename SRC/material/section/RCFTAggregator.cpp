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
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.2 $
// $Date: 2008-07-03 18:08:38 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/RCFTAggregator.cpp,v $
                                                                        
                                                                        
// File: ~/section/RCFTAggregator.C
//
// Written: Cenk Tort
// Created: February 2005
// Revision: A
//
// Purpose: This file contains the implementation for the SectionAggregator class. 

#include <stdlib.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <classTags.h>
#include <RCFTAggregator.h>
#include <MaterialResponse.h>
#include <ID.h>
#include <iostream>

#include <classTags.h>

#define maxOrder 15 //TORT

//using std::endl;
//using std::ios;
//using std::ifstream;
//using std::ofstream;
using namespace std;

// Assumes section order is less than or equal to maxOrder.
// Can increase if needed!!!
double RCFTAggregator::workArea[2*maxOrder*(maxOrder+1)];
int    RCFTAggregator::codeArea[maxOrder];

// constructors:
RCFTAggregator::RCFTAggregator (int tag, SectionForceDeformation &theSec,
				      int numAdds, UniaxialMaterial **theAdds,
				      const ID &addCodes): 
  SectionForceDeformation(tag, SEC_TAG_Aggregator), 
  theSection(0), theAdditions(0), matCodes(0), numMats(numAdds),
  e(0), s(0), ks(0), ksCommit(0), fs(0), fsCommit(0), theCode(0),
  otherDbTag(0)
{
    RCFTFiberSection3D &rcftsec = dynamic_cast<RCFTFiberSection3D&>(theSec);

    theSection = rcftsec.getCopy();

    if (!theSection) {
      opserr << "SectionAggregator::SectionAggregator -- failed to get copy of section\n";
      exit(-1);
    }

    if (!theAdds) {
      opserr << "SectionAggregator::SectionAggregator -- null uniaxial material array passed\n";
      exit(-1);
    }
    theAdditions = new UniaxialMaterial *[numMats];

    if (!theAdditions) {
      opserr << "SectionAggregator::SectionAggregator -- failed to allocate pointers\n";
      exit(-1);
    }    
    int i;
    
    for (i = 0; i < numMats; i++) {
      if (!theAdds[i]) {
	opserr << "SectionAggregator::SectionAggregator -- null uniaxial material pointer passed\n";
	exit(-1);
      }	
      theAdditions[i] = theAdds[i]->getCopy(this);
      
      if (!theAdditions[i]) {
	opserr << "SectionAggregator::SectionAggregator -- failed to copy uniaxial material\n";
	exit(-1);
      }
    }

    int order = rcftsec.getOrder()+numAdds;

    if (order > maxOrder) {
      opserr << "SectionAggregator::SectionAggregator -- order too big, need to modify the #define in SectionAggregator.cpp to " <<
	order << endln;
      exit(-1);
    }

    theCode = new ID(codeArea, order);
    e = new Vector(workArea, order);
    s = new Vector(&workArea[maxOrder], order);
    ks = new Matrix(&workArea[2*maxOrder], order, order);
    ksCommit = new Matrix(&workArea[2*maxOrder], order, order);
    fs = new Matrix(&workArea[maxOrder*(maxOrder+2)], order, order);
    fsCommit = new Matrix(&workArea[order*(order+2)], order, order);

    matCodes = new ID(addCodes);

    if (theCode == 0 || e == 0 || s == 0 || ks == 0 || fs == 0 || ksCommit == 0 || fsCommit == 0 || matCodes == 0) {
      opserr << "SectionAggregator::SectionAggregator -- out of memory\n";
      exit(-1);
    }	

}

RCFTAggregator::RCFTAggregator (int tag, int numAdds,
				      UniaxialMaterial **theAdds,
				      const ID &addCodes): 
  SectionForceDeformation(tag, SEC_TAG_Aggregator), 
  theSection(0), theAdditions(0), matCodes(0), numMats(numAdds),
  e(0), s(0), ks(0), ksCommit(0), fs(0), fsCommit(0), theCode(0),
  otherDbTag(0)
{
  if (!theAdds) {
    opserr << "SectionAggregator::SectionAggregator -- null uniaxial material array passed\n";
    exit(-1);
  }

  theAdditions = new UniaxialMaterial *[numMats];

  if (!theAdditions) {
    opserr << "SectionAggregator::SectionAggregator -- failed to allocate pointers\n";
    exit(-1);
  }    
    int i;
    
    for (i = 0; i < numMats; i++) {
      if (!theAdds[i]) {
	opserr << "SectionAggregator::SectionAggregator -- null uniaxial material pointer passed\n";
	exit(-1);
      }		
	
      theAdditions[i] = theAdds[i]->getCopy(this);
      
      if (!theAdditions[i]) {
	opserr << "SectionAggregator::SectionAggregator -- failed to copy uniaxial material\n";
	exit(-1);
      }
    }

    int order = numAdds;

    if (order > maxOrder) {
      opserr << "SectionAggregator::SectionAggregator -- order too big, need to modify the #define in SectionAggregator.cpp to %d\n";
      exit(-1);
    }

    theCode = new ID(codeArea, order);
    e = new Vector(workArea, order);
    s = new Vector(&workArea[maxOrder], order);
    ks = new Matrix(&workArea[2*maxOrder], order, order);
    ksCommit = new Matrix(&workArea[2*maxOrder],order, order);
    fs = new Matrix(&workArea[maxOrder*(maxOrder+2)], order, order);
    fsCommit = new Matrix(&workArea[order*(order+2)], order, order);

    matCodes = new ID(addCodes);

    if (theCode == 0 || e == 0 || s == 0 || ks == 0 || fs == 0 || matCodes == 0) {
      opserr << "SectionAggregator::SectionAggregator -- out of memory\n";
      exit(-1);
    }
}

RCFTAggregator::RCFTAggregator (int tag, SectionForceDeformation &theSec,
				      UniaxialMaterial &theAddition, int c) :
  SectionForceDeformation(tag, SEC_TAG_Aggregator),
  theSection(0), theAdditions(0), matCodes(0), numMats(1),
  e(0), s(0), ks(0), ksCommit(0), fs(0), fsCommit(0), theCode(0),
  otherDbTag(0)
{
  theSection = theSec.getCopy();
  
  if (!theSection) {
    opserr << "SectionAggregator::SectionAggregator -- failed to get copy of section\n";
    exit(-1);
  }

  theAdditions = new UniaxialMaterial *[1];
  
  theAdditions[0] = theAddition.getCopy(this);
  
  if (!theAdditions[0]) {
    opserr << "SectionAggregator::SectionAggregator -- failed to copy uniaxial material\n";
    exit(-1);
  }
    
  matCodes = new ID(1);
  (*matCodes)(0) = c;
  
  int order = theSec.getOrder()+1;
  
  if (order > maxOrder) {
    opserr << "SectionAggregator::SectionAggregator -- order too big, need to modify the #define in SectionAggregator.cpp to %d\n";
    exit(-1);
  }
  
  theCode = new ID(codeArea, order);
  e = new Vector(workArea, order);
  s = new Vector(&workArea[maxOrder], order);
  ks = new Matrix(&workArea[2*maxOrder], order, order);
  ksCommit = new Matrix(&workArea[2*maxOrder], order, order);
  fs = new Matrix(&workArea[maxOrder*(maxOrder+2)], order, order);
  fsCommit = new Matrix(&workArea[order*(order+2)], order, order);

  int i;

  if (theCode == 0 || e == 0 || s == 0 || ks == 0 || fs == 0 || matCodes == 0) {
    opserr << "SectionAggregator::SectionAggregator -- out of memory\n";
    exit(-1);
  }
}

// constructor for blank object that recvSelf needs to be invoked upon
RCFTAggregator::RCFTAggregator():
  SectionForceDeformation(0, SEC_TAG_Aggregator),
  theSection(0), theAdditions(0), matCodes(0), numMats(0), 
  e(0), s(0), ks(0), ksCommit(0), fs(0), fsCommit(0), theCode(0),
  otherDbTag(0)
{

}

// destructor:
RCFTAggregator::~RCFTAggregator()
{
   int i;

   if (theSection)
       delete theSection;

   for (i = 0; i < numMats; i++)
       if (theAdditions[i])
	   delete theAdditions[i];

   if (theAdditions)
       delete [] theAdditions;
   
   if (e != 0)
     delete e;

   if (s != 0)
     delete s;

   if (ks != 0)
     delete ks;

   if (ksCommit != 0)
     delete ksCommit;

   if (fs != 0)
     delete fs;

   if (fsCommit != 0)
     delete fsCommit;

   if (theCode != 0)
     delete theCode;

   if (matCodes != 0)
     delete matCodes;
}

int RCFTAggregator::setTrialSectionDeformation (const Vector &def)
{
  int ret = 0;
  int i = 0;

  int theSectionOrder = 0;

  if (theSection) {
    theSectionOrder = theSection->getOrder();
    Vector v(workArea, theSectionOrder);
    
    for (i = 0; i < theSectionOrder; i++)
      v(i) = def(i);
    
    ret = theSection->setTrialSectionDeformation(v);
  }

  int order = theSectionOrder + numMats;
  
  for ( ; i < order; i++)
    ret += theAdditions[i-theSectionOrder]->setTrialStrain(def(i));
  
  return ret;
}

const Vector &
RCFTAggregator::getSectionDeformation(void)
{
  int i = 0;

  int theSectionOrder = 0;
    
  if (theSection) {
    const Vector &eSec = theSection->getSectionDeformation();
    theSectionOrder = theSection->getOrder();
    
    for (i = 0; i < theSectionOrder; i++)
      (*e)(i) = eSec(i);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*e)(i) = theAdditions[i-theSectionOrder]->getStrain();

  return *e;
}

const Matrix &
RCFTAggregator::getSectionTangent(void)
{
  int i = 0;
  int theSectionOrder = 0;

  // Zero before assembly
  ks->Zero();

  if (theSection) {
    const Matrix &kSec = theSection->getSectionTangent();
    theSectionOrder = theSection->getOrder();
    //output<<"\n before slip is added \n"<<endl;
    //for( i=0;i<12;i++){
    //    output<<kSec(i,0)<<"  "<<kSec(i,1)<<"  "<<kSec(i,2)<<"  "<<kSec(i,3)<<"  "<<kSec(i,4)<<"  "<<
    //        kSec(i,5)<<"  "<<kSec(i,6)<<"  "<<kSec(i,7)<<"  "<<kSec(i,8)<<"  "<<kSec(i,9)<<"  "<<
    //        kSec(i,10)<<"  "<<kSec(i,11)<<endl;
    //}
    for (i = 0; i < theSectionOrder; i++)
      for (int j = 0; j < theSectionOrder; j++)
	(*ks)(i,j) = kSec(i,j);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*ks)(i,i) = theAdditions[i-theSectionOrder]->getTangent();

  //output<<"\n after slip is added \n"<<endl;
  //for( i=0;i<13;i++){
  //  output<<(*ks)(i,0)<<"  "<<(*ks)(i,1)<<"  "<<(*ks)(i,2)<<"  "<<(*ks)(i,3)<<"  "<<(*ks)(i,4)<<"  "<<
  //          (*ks)(i,5)<<"  "<<(*ks)(i,6)<<"  "<<(*ks)(i,7)<<"  "<<(*ks)(i,8)<<"  "<<(*ks)(i,9)<<"  "<<
  //          (*ks)(i,10)<<"  "<<(*ks)(i,11)<<"  "<<(*ks)(i,12)<<endl;
  //}    
	  
  return *ks;
}

const Matrix &
RCFTAggregator::getSectionCommitTangent(void)
{
  int i = 0;
  int theSectionOrder = 0;

  RCFTSlipMaterial **theslipAdditions = (RCFTSlipMaterial **) theAdditions;

  // Zero before assembly
  ksCommit->Zero();

  //cout<<"\n BEFOREBEFORE CASTING \n"<<endl;

  RCFTFiberSection3D *rcftsec = (RCFTFiberSection3D*) theSection;

  //cout<<"\n BEFORE SLIP IS ADDED AGGREGATOR \n"<<endl;

  if (theSection) {	

      //output<<"\n INSIDELOOP \n";	  
      const Matrix &kSec = rcftsec->getSectionCommitTangent();
      //output>>kSec;
      theSectionOrder = theSection->getOrder();

      cout<<theSectionOrder<<endl;
       
      //output>>kSec;
      
      for (i = 0; i < theSectionOrder; i++)
           for (int j = 0; j < theSectionOrder; j++)
               (*ksCommit)(i,j) = kSec(i,j);
  }

  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
     (*ksCommit)(i,i) = theslipAdditions[i - theSectionOrder]->getCommitTangent();

  //output<<"\n AFTER SLIP IS ADDED \n"<<endl;

  //for( i=0;i<13;i++){
  //    output<<(*ksCommit)(i,0)<<"  "<<(*ksCommit)(i,1)<<"  "<<(*ksCommit)(i,2)<<"  "<<(*ksCommit)(i,3)<<"  "<<(*ksCommit)(i,4)<<"  "<<
  //            (*ksCommit)(i,5)<<"  "<<(*ksCommit)(i,6)<<"  "<<(*ksCommit)(i,7)<<"  "<<(*ksCommit)(i,8)<<"  "<<(*ksCommit)(i,9)<<"  "<<
  //            (*ksCommit)(i,10)<<"  "<<(*ksCommit)(i,11)<<"  "<<(*ksCommit)(i,12)<<endl;
  //}

  return *ksCommit;
}

const Matrix &
RCFTAggregator::getInitialTangent(void)
{
  int i = 0;

  int theSectionOrder = 0;

  // Zero before assembly
  ks->Zero();

  if (theSection) {
    const Matrix &kSec = theSection->getInitialTangent();
    theSectionOrder = theSection->getOrder();

    for (i = 0; i < theSectionOrder; i++)
      for (int j = 0; j < theSectionOrder; j++)
	(*ks)(i,j) = kSec(i,j);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*ks)(i,i) = theAdditions[i-theSectionOrder]->getInitialTangent();

  //for( i=0;i<13;i++){
  //  output<<(*ks)(i,0)<<"  "<<(*ks)(i,1)<<"  "<<(*ks)(i,2)<<"  "<<(*ks)(i,3)<<"  "<<(*ks)(i,4)<<"  "<<
  //          (*ks)(i,5)<<"  "<<(*ks)(i,6)<<"  "<<(*ks)(i,7)<<"  "<<(*ks)(i,8)<<"  "<<(*ks)(i,9)<<"  "<<
  //          (*ks)(i,10)<<"  "<<(*ks)(i,11)<<"  "<<(*ks)(i,12)<<endl;
  //}
  
  return *ks;
}

const Matrix &
RCFTAggregator::getSectionFlexibility(void)
{
  int i = 0;
    
  int theSectionOrder = 0;

  // Zero before assembly
  fs->Zero();
  
  if (theSection) {
    const Matrix &fSec = theSection->getSectionFlexibility();
    theSectionOrder = theSection->getOrder();

    for (i = 0; i < theSectionOrder; i++)
      for (int j = 0; j < theSectionOrder; j++)
	(*fs)(i,j) = fSec(i,j);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++) {
    double k = theAdditions[i-theSectionOrder]->getTangent();
    if (k == 0.0) {
      opserr << "SectionAggregator::getSectionFlexibility -- singular section stiffness\n";
			       
      (*fs)(i,i) = 1.e14;
    }
    else
      (*fs)(i,i) = 1/k;
  }	
  
  return *fs;
}

const Matrix &
RCFTAggregator::getInitialFlexibility(void)
{
  int i = 0;
    
  int theSectionOrder = 0;

  // Zero before assembly
  fs->Zero();

  if (theSection) {
    const Matrix &fSec = theSection->getInitialFlexibility();
    theSectionOrder = theSection->getOrder();

    for (i = 0; i < theSectionOrder; i++)
      for (int j = 0; j < theSectionOrder; j++)
	(*fs)(i,j) = fSec(i,j);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++) {
    double k = theAdditions[i-theSectionOrder]->getInitialTangent();
    (*fs)(i,i) = 1.0/k;
  }	
  
  return *fs;
}

const Vector &
RCFTAggregator::getStressResultant(void)
{
  int i = 0;

  int theSectionOrder = 0;
    
  if (theSection) {
    const Vector &sSec = theSection->getStressResultant();
    theSectionOrder = theSection->getOrder();
    
    for (i = 0; i < theSectionOrder; i++)
      (*s)(i) = sSec(i);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*s)(i) = theAdditions[i-theSectionOrder]->getStress();
  
  return *s;
}

SectionForceDeformation *
RCFTAggregator::getCopy(void)
{
  RCFTAggregator *theCopy = 0;

  if (theSection)
    theCopy = new RCFTAggregator(this->getTag(), *theSection,
				    numMats, theAdditions, *matCodes);
  else
    theCopy = new RCFTAggregator(this->getTag(), numMats,
				    theAdditions, *matCodes);
  
  if (theCopy == 0) {
    opserr << "SectionAggregator::getCopy -- failed to allocate copy\n";
    exit(-1);
  }
			  
  
  return theCopy;
}

const ID&
RCFTAggregator::getType ()
{
  int i = 0;

  int theSectionOrder = 0;
    
  if (theSection) {
    const ID &secType = theSection->getType();
    theSectionOrder = theSection->getOrder();
    
    for (i = 0; i < theSectionOrder; i++)
      (*theCode)(i) = secType(i);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*theCode)(i) = (*matCodes)(i-theSectionOrder);

  return *theCode;
}

int
RCFTAggregator::getOrder () const
{
  int order = numMats;

  if (theSection != 0)
    order += theSection->getOrder();

  return order;
}

int
RCFTAggregator::commitState(void)
{
  ofstream slip;
  slip.open("slip.dat",ios::app);
  int err = 0;

  if (theSection)
    err += theSection->commitState();

  for (int i = 0; i < numMats; i++){
    err += theAdditions[i]->commitState();
    slip<<theAdditions[i]->getStrain()<<"  "<<theAdditions[i]->getStress()<<"  "<<theAdditions[i]->getTangent()<<"  "<<0<<endl;
  }

  return err;
}

int
RCFTAggregator::revertToLastCommit(void)
{
  int err = 0;
  
  int i = 0;
  
  // Revert the section
  if (theSection)
    err += theSection->revertToLastCommit();
  
  // Do the same for the uniaxial materials
  for (i = 0; i < numMats; i++)
    err += theAdditions[i]->revertToLastCommit();
  
  return err;
}	

int
RCFTAggregator::revertToStart(void)
{
  int err = 0;
  
  // Revert the section
  if (theSection)
    err += theSection->revertToStart();
  
  // Do the same for the uniaxial materials
  for (int i = 0; i < numMats; i++)
    err += theAdditions[i]->revertToStart();
  
  return err;
}

int
RCFTAggregator::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

  // Need otherDbTag since classTags ID and data ID may be the same size
  if (otherDbTag == 0) 
    otherDbTag = theChannel.getDbTag();
  
  // Create ID for tag and section order data
  static ID data(5);
  
  int order = this->getOrder();
  
  data(0) = this->getTag();
  data(1) = otherDbTag;
  data(2) = order;
  data(3) = (theSection != 0) ? theSection->getOrder() : 0;
  data(4) = numMats;

  // Send the tag and section order data
  res += theChannel.sendID(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "SectionAggregator::sendSelf -- could not send data ID\n";
			    
    return res;
  }
  
  // Determine how many classTags there are and allocate ID vector
  // for the tags and section code
  int numTags = (theSection == 0) ? numMats : numMats + 1;
  ID classTags(2*numTags + numMats);
  
  // Loop over the UniaxialMaterials filling in class and db tags
  int i, dbTag;
  for (i = 0; i < numMats; i++) {
    classTags(i) = theAdditions[i]->getClassTag();
    
    dbTag = theAdditions[i]->getDbTag();
    
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
	theAdditions[i]->setDbTag(dbTag);
    }
    
    classTags(i+numTags) = dbTag;
  }
  
  // Put the Section class and db tags into the ID vector
  if (theSection != 0) {
    classTags(numTags-1) = theSection->getClassTag();
    
    dbTag = theSection->getDbTag();
    
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
	theSection->setDbTag(dbTag);
    }
    
    classTags(2*numTags-1) = dbTag;
  }
  
  // Put the UniaxialMaterial codes into the ID vector
  int j = 2*numTags;
  for (i = 0; i < numMats; i++, j++)
    classTags(j) = (*matCodes)(i);
  
  // Send the material class and db tags and section code
  res += theChannel.sendID(otherDbTag, cTag, classTags);
  if (res < 0) {
    opserr << "SectionAggregator::sendSelf -- could not send classTags ID\n";
    return res;
  }

  // Ask the UniaxialMaterials to send themselves
  for (i = 0; i < numMats; i++) {
    res += theAdditions[i]->sendSelf(cTag, theChannel);
    if (res < 0) {
      opserr << "SectionAggregator::sendSelf -- could not send UniaxialMaterial, i = " << i << endln;
      return res;
    }
  }
  
  // Ask the Section to send itself
  if (theSection != 0) {
    res += theSection->sendSelf(cTag, theChannel);
    if (res < 0) {
      opserr << "SectionAggregator::sendSelf -- could not send SectionForceDeformation\n";
      return res;
    }
  }
  
  return res;
}


int
RCFTAggregator::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  // Create an ID and receive tag and section order
  static ID data(5);
  res += theChannel.recvID(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "SectionAggregator::recvSelf -- could not receive data ID\n";
    return res;
  }
  
  this->setTag(data(0));
  otherDbTag = data(1);
  int order = data(2);
  int theSectionOrder = data(3);
  numMats = data(4);

  if (order > 0) {
    if (e == 0 || e->Size() != order) {
      if (e != 0) {
	delete e;
	delete s;
	delete ks;
	delete fs;
	delete theCode;
      }
      e = new Vector(workArea, order);
      s = new Vector(&workArea[maxOrder], order);
      ks = new Matrix(&workArea[2*maxOrder], order, order);
      fs = new Matrix(&workArea[maxOrder*(maxOrder+2)], order, order);
      theCode = new ID(codeArea, order);
    }
  }

  if (numMats > 0) {
    if (matCodes == 0 || matCodes->Size() != numMats) {
      if (matCodes != 0)
	delete matCodes;

      matCodes = new ID(numMats);
    }
  }

  // Determine how many classTags there are and allocate ID vector
  int numTags = (theSectionOrder == 0) ? numMats : numMats + 1;
  ID classTags(numTags*2 + numMats);
  
  // Receive the material class and db tags
  res += theChannel.recvID(otherDbTag, cTag, classTags);
  if (res < 0) {
    opserr << "SectionAggregator::recvSelf -- could not receive classTags ID\n";
    return res;
  }

  // Check if null pointer, allocate if so
  if (theAdditions == 0) {
    theAdditions = new UniaxialMaterial *[numMats];
    if (theAdditions == 0) {
      opserr << "SectionAggregator::recvSelf -- could not allocate UniaxialMaterial array\n";
      return -1;
    }
    // Set pointers to null ... will get allocated by theBroker
    for (int j = 0; j < numMats; j++)
      theAdditions[j] = 0;
  }
  
  // Loop over the UniaxialMaterials
  int i, classTag;
  for (i = 0; i < numMats; i++) {
    classTag = classTags(i);
    
    // Check if the UniaxialMaterial is null; if so, get a new one
    if (theAdditions[i] == 0)
      theAdditions[i] = theBroker.getNewUniaxialMaterial(classTag);
    
    // Check that the UniaxialMaterial is of the right type; if not, delete
    // the current one and get a new one of the right type
    else if (theAdditions[i]->getClassTag() != classTag) {
      delete theAdditions[i];
      theAdditions[i] = theBroker.getNewUniaxialMaterial(classTag);
    }
    
    // Check if either allocation failed
    if (theAdditions[i] == 0) {
      opserr << "SectionAggregator::recvSelf -- could not get UniaxialMaterial, i = " << i << endln;
      return -1;
    }
    
    // Now, receive the UniaxialMaterial
    theAdditions[i]->setDbTag(classTags(i+numTags));
    res += theAdditions[i]->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "SectionAggregator::recvSelf -- could not receive UniaxialMaterial, i = " << i << endln;
      return res;
    }
  }

  // If there is no Section to receive, return
  if (theSectionOrder == 0)
    return res;
  
  classTag = classTags(numTags-1);
  
  // Check if the Section is null; if so, get a new one
  if (theSection == 0)
    theSection = theBroker.getNewSection(classTag);
  
  // Check that the Section is of the right type; if not, delete
  // the current one and get a new one of the right type
  else if (theSection->getClassTag() != classTag) {
    delete theSection;
    theSection = theBroker.getNewSection(classTag);
  }
  
  // Check if either allocation failed
  if (theSection == 0) {
    opserr << "SectionAggregator::recvSelf -- could not get a SectionForceDeformation\n";
    return -1;
  }

  // Now, receive the Section
  theSection->setDbTag(classTags(2*numTags-1));
  res += theSection->recvSelf(cTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "SectionAggregator::recvSelf -- could not receive SectionForceDeformation\n";
    return res;
  }
  
  // Fill in the section code
  int j = 2*numTags;
  for (i = 0; i < numMats; i++, j++)
    (*matCodes)(i) = classTags(j);

  return res;
}

void
RCFTAggregator::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {
      theSection->Print(s, flag);
  } else {
    s << "\nSection Aggregator, tag: " << this->getTag() << endln;
    if (theSection) {
      s << "\tSection, tag: " << theSection->getTag() << endln;
      theSection->Print(s, flag);
    }
    s << "\tUniaxial Additions" << endln;
    for (int i = 0; i < numMats; i++)
      s << "\t\tUniaxial Material, tag: " << theAdditions[i]->getTag() << endln;
    s << "\tUniaxial codes " << *matCodes << endln;
  }
}

Response*
RCFTAggregator::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  // See if the response is one of the defaults
  Response *res = SectionForceDeformation::setResponse(argv, argc, output);
  if (res != 0)
    return res;
  
  // If not, forward the request to the section (need to do this to get fiber response)
  // CURRENTLY NOT SENDING ANYTHING OFF TO THE UniaxialMaterials ... Probably
  // don't need anything more from them than stress, strain, and stiffness, 
  // which are covered in base class method ... can change if need arises
  else if (theSection != 0)
    return theSection->setResponse(argv, argc, output);
  
  else
    return 0;
}

int
RCFTAggregator::getResponse(int responseID, Information &info)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  return SectionForceDeformation::getResponse(responseID, info);
}

int
RCFTAggregator::setVariable(const char *argv)
{
  // Axial strain
  if (strcmp(argv,"axialStrain") == 0)
    return 1;
  // Curvature about the section z-axis
  else if (strcmp(argv,"curvatureZ") == 0)
    return 2;
  // Curvature about the section y-axis
  else if (strcmp(argv,"curvatureY") == 0)
    return 3;
  else
    return -1;
}

int
RCFTAggregator::getVariable(int variableID, double &info)
{
  int i;

  info = 0.0;

  int order = numMats;
  if (theSection != 0)
    order += theSection->getOrder();

  const Vector &e = this->getSectionDeformation();
  const ID &code  = this->getType();

  switch (variableID) {
  case 1:	// Axial strain
    // Series model ... add all sources of deformation
    for (i = 0; i < order; i++)
      if (code(i) == SECTION_RESPONSE_P)
	info += e(i);
    return 0;
  case 2:	// Curvature about the section z-axis
    for (i = 0; i < order; i++)
      if (code(i) == SECTION_RESPONSE_MZ)
	info += e(i);
    return 0;
  case 3:	// Curvature about the section y-axis
    for (i = 0; i < order; i++)
      if (code(i) == SECTION_RESPONSE_MY)
	info += e(i);
    return 0;
  default:
    return -1;
  }
}
