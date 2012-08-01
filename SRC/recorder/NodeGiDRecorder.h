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
                                                                        
// $Revision: 1.17 $
// $Date: 2010/04/23 22:47:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/NodeGiDRecorder.h,v $
                                                                        
#ifndef NodeGiDRecorder_h
#define NodeGiDRecorder_h

// Written: fmk 
// Created: 09/00
// Revision: A
//
// Description: This file contains the class definition for 
// NodeGiDRecorder. A NodeGiDRecorder is used to store the specified nodal dof responses
// for the specified nodes in a file.
//
// What: "@(#) NodeGiDRecorder.h, revA"


#include <Recorder.h>
#include <ID.h>
#include <Vector.h>

class Domain;
class FE_Datastore;
class Node;
class String;

class NodeGiDRecorder: public Recorder
{
  public:
    NodeGiDRecorder();
    NodeGiDRecorder(const ID &theDof, 
		 const ID *theNodes, 
		 const char *dataToStore,
		 Domain &theDomain,
		 OPS_Stream &theOutputHandler,
		 double deltaT = 0.0,
		 bool echoTimeFlag = true); 
    
    ~NodeGiDRecorder();

    int record(int commitTag, double timeStamp);

    int domainChanged(void);    
    int setDomain(Domain &theDomain);
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

  protected:

  private:	
    int initialize(void);

    ID *theDofs;
    ID *theNodalTags;
    Node **theNodes;
    Vector response;

    Domain *theDomain;
    OPS_Stream *theOutputHandler;

    bool echoTimeFlag;   // flag indicating whether time to be included in o/p
    int dataFlag;        // flag indicating what it is to be stored in recorder

    double deltaT;
    double nextTimeStampToRecord;

    bool initializationDone;
    int numValidNodes;

    int addColumnInfo;
};

#endif
