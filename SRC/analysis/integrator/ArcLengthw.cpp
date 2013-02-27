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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/analysis/integrator/ArcLengthw.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/ArcLengthw.C
// 
// Written: FMK and Cenk Tort
// Created: 02/2006
// Revision: A
//


#include <ArcLengthw.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <EquiSolnAlgo.h>
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

ArcLengthw::ArcLengthw(double ifactor, int jd)
:StaticIntegrator(INTEGRATOR_TAGS_ArcLengthw),
 iFactor(ifactor), Jd(jd), dWi(0.0), deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0), 
 phat(0), deltaLambdaStep(0.0), currentLambda(0.0), iflag(0)
{
}

ArcLengthw::~ArcLengthw()
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
ArcLengthw::newStep(void)
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
 
    factor<<"currentLambda"<<endln;
    factor<<currentLambda<<endln;

    // determine dUhat
    this->formTangent();
    theLinSOE->setB(*phat);
    theLinSOE->solve();
    (*deltaUhat) = theLinSOE->getX();
    Vector &dUhat = *deltaUhat;

    factor<<"dUhat"<<endln;
    //factor>>dUhat;

    int size = dUhat.Size();
    int i = 0;
    double sum = 0.0;
    int Ji_1 = 0;
    double dLambda = 0.0;
    factor<<"dWibefore"<<endln;
    factor<<dWi<<endln;
    factor<<"*phat"<<endln;
    //factor>>*phat;
    factor<<"dUhat"<<endln;
    //factor>>dUhat;
    factor<<"iFactor"<<endln;
    factor<<iFactor<<endln;
    factor<<"Jd"<<endln;
    factor<<Jd<<endln;
    factor<<"iflag"<<endln;
    factor<<iflag<<endln;
    double dJd = Jd;
    double dJi_1 = 1.0;
    if( iflag == 0 ){
       dWi = ( (*phat) ^ dUhat ) * iFactor * iFactor;
       dLambda = iFactor; 
       iflag = 1;  
    }
    else if( iflag == 1 ){
	   Ji_1 = 10; //theAlgo->getNumIteration();
       dJi_1 = Ji_1;
       dWi = dWi * pow(( dJd / dJi_1 ),0.01);
       dLambda = dWi / ( (*phat)^(dUhat) );
    }
    if( Ji_1 >0){
    factor<<"Jd/Ji-1"<<endln;
    factor<<dJd/dJi_1<<endln;
    }
    factor<<"iflag"<<endln;
    factor<<iflag<<endln;

    factor<<"Ji_1"<<endln;
    factor<<Ji_1<<endln;

    factor<<"dWi"<<endln;
    factor<<dWi<<endln;
    
    deltaLambdaStep = dLambda;
    currentLambda += dLambda;

    (*deltaU) = dUhat;
    (*deltaU) *= dLambda;
    (*deltaUstep) = (*deltaU);

    // update model with delta lambda and delta U
    theModel->incrDisp(*deltaU);    
    theModel->applyLoadDomain(currentLambda);    
    theModel->updateDomain();

    return 0;
}

int
ArcLengthw::update(const Vector &dU)
{
    ofstream factor;
    factor.open("FS.dat",ios::app);

    factor<<"insideupdate"<<endln;
    //factor>>dU;
    factor<<"insideupdate1"<<endln;

    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING ArcLengthw::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    (*deltaUbar) = dU; // have to do this as the SOE is gonna change

    // determine dUhat    
    theLinSOE->setB(*phat);
    theLinSOE->solve();
    (*deltaUhat) = theLinSOE->getX();    

    double dLambda = -((*phat)^(*deltaUbar))/((*phat)^(*deltaUhat));

    (*deltaU) = (*deltaUbar);    
    deltaU->addVector(1.0, *deltaUhat,dLambda);
    
    // update dU and dlambda
    (*deltaUstep) += *deltaU;
    deltaLambdaStep += dLambda;
    currentLambda += dLambda;

    // update the model
    theModel->incrDisp(*deltaU);    
    theModel->applyLoadDomain(currentLambda);    
    theModel->updateDomain();
    
    // set the X soln in linearSOE to be deltaU for convergence Test
    theLinSOE->setX(*deltaU);

    return 0;
}



int 
ArcLengthw::domainChanged(void)
{
    // we first create the Vectors needed
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING ArcLengthw::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }    
    int size = theModel->getNumEqn(); // ask model in case N+1 space

    if (deltaUhat == 0 || deltaUhat->Size() != size) { // create new Vector
	if (deltaUhat != 0)
	    delete deltaUhat;   // delete the old
	deltaUhat = new Vector(size);
	if (deltaUhat == 0 || deltaUhat->Size() != size) { // check got it
	    opserr << "FATAL ArcLengthw::domainChanged() - ran out of memory for";
	    opserr << " deltaUhat Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUbar == 0 || deltaUbar->Size() != size) { // create new Vector
	if (deltaUbar != 0)
	    delete deltaUbar;   // delete the old
	deltaUbar = new Vector(size);
	if (deltaUbar == 0 || deltaUbar->Size() != size) { // check got it
	    opserr << "FATAL ArcLengthw::domainChanged() - ran out of memory for";
	    opserr << " deltaUbar Vector of size " << size << endln;
	    exit(-1);
	}
    }

    
    if (deltaU == 0 || deltaU->Size() != size) { // create new Vector
	if (deltaU != 0)
	    delete deltaU;   // delete the old
	deltaU = new Vector(size);
	if (deltaU == 0 || deltaU->Size() != size) { // check got it
	    opserr << "FATAL ArcLengthw::domainChanged() - ran out of memory for";
	    opserr << " deltaU Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	if (deltaUstep != 0)
	    delete deltaUstep;  
	deltaUstep = new Vector(size);
	if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	    opserr << "FATAL ArcLengthw::domainChanged() - ran out of memory for";
	    opserr << " deltaUstep Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (phat == 0 || phat->Size() != size) { 
	if (phat != 0)
	    delete phat;  
	phat = new Vector(size);
	if (phat == 0 || phat->Size() != size) { 
	    opserr << "FATAL ArcLengthw::domainChanged() - ran out of memory for";
	    opserr << " phat Vector of size " << size << endln;
	    exit(-1);
	}
    }    

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
ArcLengthw::sendSelf(int cTag,
		    Channel &theChannel)
{
  Vector data(4);
  data(0) = iFactor;
  data(1) = Jd;
  data(2) = deltaLambdaStep;
  data(3) = currentLambda;

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "ArcLengthw::sendSelf() - failed to send the data\n";
      return -1;
  }
  return 0;
}


int
ArcLengthw::recvSelf(int cTag,
		    Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  Vector data(4);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "ArcLengthw::sendSelf() - failed to send the data\n";
      return -1;
  }      

  // set the data
  iFactor = data(0);
  Jd = data(1);
  deltaLambdaStep = data(2);
  currentLambda = data(3);
  return 0;
}

void
ArcLengthw::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
	double cLambda = theModel->getCurrentDomainTime();
	s << "\t ArcLengthw - currentLambda: " << cLambda;
	s << "  Initial Load Factor: " << iFactor <<  "  Jd: ";
	s << Jd << endln;
    } else 
	s << "\t ArcLengthw - no associated AnalysisModel\n";
}








