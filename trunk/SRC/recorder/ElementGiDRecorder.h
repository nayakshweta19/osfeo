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
                                                                        
// $Revision: 1.14 $
// $Date: 2009/04/14 21:14:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/ElementGiDRecorder.h,v $
                                                                        
                                                                        
#ifndef ElementGiDRecorder_h
#define ElementGiDRecorder_h


// Written:
// Created:
// Revision:
//
// Description: This file contains the class definition for ElementGiDRecorder.
// A ElementGiDRecorder is used to obtain a response from an element during 
// the analysis.
//
// What: "@(#) ElementGiDRecorder.h, revA"

#include <Recorder.h>
#include <Information.h>
#include <ID.h>

class Domain;
class Vector;
class Matrix;
class Element;
class Response;
class FE_Datastore;

class ElementGiDRecorder: public Recorder
{
  public:
    ElementGiDRecorder();
    ElementGiDRecorder(ID &eleIDs, 
		    const char **argv, 
		    int argc,
		    bool echoTime, 
		    Domain &theDomain, 
		    OPS_Stream &theOutputHandler,
		    double deltaT = 0.0);

    ~ElementGiDRecorder();

    int record(int commitTag, double timeStamp);
    int restart(void);    

    int setDomain(Domain &theDomain);
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
  protected:

    
  private:	
    int initialize(void);

    int numEle;
    ID *eleID;

    Response **theResponses;

    Domain *theDomain;
    OPS_Stream *theOutputHandler;

    bool echoTimeFlag;             // flag indicating if pseudo time also printed

    double deltaT;
    double nextTimeStampToRecord;
	int stepN;
	bool hasLinear;
	bool hasTri3;
	bool hasQuad4;
	bool hasQuad8;
	bool hasQuad9;
	bool hasBrick;

    //Vector *data;
    bool initializationDone;
    char **responseArgs;
    int numArgs;

    int addColumnInfo;
};


#endif
