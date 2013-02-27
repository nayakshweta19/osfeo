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
// $Date: 2008-07-03 17:58:49 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/analysis/integrator/ArcLengthy.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/ArcLengthw.C
// 
// Written: FMK and Cenk Tort
// Created: 02/2006
// Revision: A
//


#include <ArcLengthy.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>

//using std::ofstream;
//using std::ios;
//using std::endl;
using namespace std;

int ArcLengthy::iflag = 0;
int ArcLengthy::sign= 1;
double ArcLengthy::normdeltaUhat11 = 0;
Vector *ArcLengthy::deltaUhat1 = 0;

ArcLengthy::ArcLengthy(double lambda1)
:StaticIntegrator(INTEGRATOR_TAGS_ArcLengthy),
 deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0), 
 phat(0), deltaLambdaStep(0.0), currentLambda(0.0), Lambda1(lambda1) //, deltaUhat1(0) 
{

}

ArcLengthy::~ArcLengthy()
{
    // delete any vector object created
    if (deltaUhat != 0)
	delete deltaUhat;
    if (deltaU != 0)
	delete deltaU;
    if (deltaUstep != 0)
	delete deltaUstep;
    if (deltaUbar != 0)
	delete deltaUbar;
    if (phat != 0)
	delete phat;
}

int
ArcLengthy::newStep(void)
{
    ofstream factor;
    factor.open("factor.dat",ios::app);
    // get pointers to AnalysisModel and LinearSOE
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING ArcLengthw::newStep() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    // get the current load factor
    currentLambda = theModel->getCurrentDomainTime();
 
    //factor<<"currentLambda"<<endl;
    //factor<<currentLambda<<endl;

    //factor<<"Lambda1"<<endl;
    //factor<<Lambda1<<endl;

    // determine dUhat
    this->formTangent();
    theLinSOE->setB(*phat);
    theLinSOE->solve();
    (*deltaUhat) = theLinSOE->getX();
    Vector &dUhat = *deltaUhat;

    factor<<"dUhat"<<endl;
    //factor>>dUhat;

    //factor<<"deltaUhat1"<<endl;
    //factor>>(*deltaUhat1);

    int size = dUhat.Size();
    int i = 0;
    double dLambda = 0.0;
    double GSP; 
    //factor<<"*phat"<<endl;
    //factor>>*phat;
    //factor<<"dUhat"<<endl;
    //factor>>dUhat;
    //factor<<"iflag"<<endl;
    //factor<<iflag<<endl;
    if( iflag == 0 ){
       //factor<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;	    
       normdeltaUhat11 = (*deltaUhat) ^ (*deltaUhat);
       GSP = 1.0;    
       dLambda = Lambda1; 
       iflag = 1;  
    }
    else if( iflag == 1 ){
       GSP = normdeltaUhat11 / ( (*deltaUhat1) ^ (*deltaUhat) );	 
       if( GSP < 0.0 ){
	  sign = - sign;     
          //dLambda = -Lambda1 * sqrt(fabs(GSP));
       }
       //else if( GSP > 0.0 ){
       //   dLambda = Lambda1 * sqrt(fabs(GSP));
       //} 
       dLambda = sign * Lambda1 * sqrt(fabs(GSP));  
    }
    //factor<<"deltaUhat1"<<endl;
    //factor>>*deltaUhat1;
    //factor<<"dektaUhat"<<endl;
    //factor>>(*deltaUhat);
    //factor<<"normdeltaUhat11"<<endl;
    //factor<<normdeltaUhat11<<endl;
    //factor<<"iflag"<<endl;
    //factor<<iflag<<endl;
    //factor<<"GSP"<<endl;
    //factor<<GSP<<endl;    

    deltaLambdaStep = dLambda;
    currentLambda += dLambda;

    (*deltaU) = dUhat;
    (*deltaU) *= dLambda;
    (*deltaUstep) = (*deltaU);

    // update model with delta lambda and delta U
    theModel->incrDisp(*deltaU);    
    theModel->applyLoadDomain(currentLambda);    
    theModel->updateDomain();
    (*deltaUhat1) = (*deltaUhat);
    //factor<<"deltaUhat1"<<endl;
    //factor>>(*deltaUhat1);
    return 0;
}

int
ArcLengthy::update(const Vector &dU)
{
    ofstream factor;
    factor.open("factor.dat",ios::app);

    //factor<<"insideupdate"<<endl;
    //factor>>dU;
    //factor<<"insideupdate1"<<endl;

    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING ArcLengthy::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    (*deltaUbar) = dU; // have to do this as the SOE is gonna change

    // determine dUhat    
    theLinSOE->setB(*phat);
    theLinSOE->solve();
    (*deltaUhat) = theLinSOE->getX();    

    double dLambda = -((*deltaUhat1)^(*deltaUbar))/((*deltaUhat1)^(*deltaUhat));

    //factor<<"dLambda"<<endl;
    //factor<<dLambda<<endl;
    
    (*deltaU) = (*deltaUbar);    
    deltaU->addVector(1.0, *deltaUhat,dLambda);
    
    // update dU and dlambda
    (*deltaUstep) += *deltaU;
    deltaLambdaStep += dLambda;
    currentLambda += dLambda;

    //factor<<"currentLambda"<<endl;
    //factor<<currentLambda<<endl;

    // update the model
    theModel->incrDisp(*deltaU);    
    theModel->applyLoadDomain(currentLambda);    
    theModel->updateDomain();

    //factor<<"updatedomain"<<endl;
    
    // set the X soln in linearSOE to be deltaU for convergence Test
    theLinSOE->setX(*deltaU);

    return 0;
}



int 
ArcLengthy::domainChanged(void)
{
    // we first create the Vectors needed
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING ArcLengthy::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }    
    int size = theModel->getNumEqn(); // ask model in case N+1 space

    if (deltaUhat == 0 || deltaUhat->Size() != size) { // create new Vector
	if (deltaUhat != 0)
	    delete deltaUhat;   // delete the old
	deltaUhat = new Vector(size);
	if (deltaUhat == 0 || deltaUhat->Size() != size) { // check got it
	    opserr << "FATAL ArcLengthy::domainChanged() - ran out of memory for";
	    opserr << " deltaUhat Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUhat1 == 0 || deltaUhat1->Size() != size) { // create new Vector
	if (deltaUhat1 != 0)
	    delete deltaUhat1;   // delete the old
	deltaUhat1 = new Vector(size);
	if (deltaUhat1 == 0 || deltaUhat1->Size() != size) { // check got it
	    opserr << "FATAL ArcLengthy::domainChanged() - ran out of memory for";
	    opserr << " deltaUhat1 Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUbar == 0 || deltaUbar->Size() != size) { // create new Vector
	if (deltaUbar != 0)
	    delete deltaUbar;   // delete the old
	deltaUbar = new Vector(size);
	if (deltaUbar == 0 || deltaUbar->Size() != size) { // check got it
	    opserr << "FATAL ArcLengthy::domainChanged() - ran out of memory for";
	    opserr << " deltaUbar Vector of size " << size << endln;
	    exit(-1);
	}
    }

    
    if (deltaU == 0 || deltaU->Size() != size) { // create new Vector
	if (deltaU != 0)
	    delete deltaU;   // delete the old
	deltaU = new Vector(size);
	if (deltaU == 0 || deltaU->Size() != size) { // check got it
	    opserr << "FATAL ArcLengthy::domainChanged() - ran out of memory for";
	    opserr << " deltaU Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	if (deltaUstep != 0)
	    delete deltaUstep;  
	deltaUstep = new Vector(size);
	if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	    opserr << "FATAL ArcLengthy::domainChanged() - ran out of memory for";
	    opserr << " deltaUstep Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (phat == 0 || phat->Size() != size) { 
	if (phat != 0)
	    delete phat;  
	phat = new Vector(size);
	if (phat == 0 || phat->Size() != size) { 
	    opserr << "FATAL ArcLengthy::domainChanged() - ran out of memory for";
	    opserr << " phat Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUhat1 == 0)
       deltaUhat1  = new Vector (size);

    // now we have to determine phat
    // do this by incrementing lambda by 1, applying load
    // and getting phat from unbalance.
    currentLambda = theModel->getCurrentDomainTime();
    currentLambda += 1.0;
    theModel->applyLoadDomain(currentLambda);    
    this->formUnbalance(); // NOTE: this assumes unbalance at last was 0
    (*phat) = theLinSOE->getB();
    currentLambda -= 1.0;
    theModel->setCurrentDomainTime(currentLambda);    
    
    return 0;
}

int
ArcLengthy::sendSelf(int cTag,
		    Channel &theChannel)
{
  Vector data(3);
  data(0) = Lambda1;
  data(1) = deltaLambdaStep;
  data(2) = currentLambda;

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "ArcLengthw::sendSelf() - failed to send the data\n";
      return -1;
  }
  return 0;
}


int
ArcLengthy::recvSelf(int cTag,
		    Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  Vector data(3);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "ArcLengthw::sendSelf() - failed to send the data\n";
      return -1;
  }      

  // set the data
  Lambda1 = data(0);
  deltaLambdaStep = data(1);
  currentLambda = data(2);
  return 0;
}

void
ArcLengthy::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
	double cLambda = theModel->getCurrentDomainTime();
    } else 
	s << "\t ArcLengthw - no associated AnalysisModel\n";
}








