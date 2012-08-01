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
                                                                        
// $Revision: 1.41 $
// $Date: 2010-06-01 23:42:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/NodeGiDRecorder.cpp,v $
                                                                        
// Written: fmk 
//
// Description: This file contains the class definition for NodeGiDRecorder.
// A NodeGiDRecorder is used to record the specified dof responses 
// at a collection of nodes over an analysis. (between commitTag of 0 and
// last commitTag).
//
// What: "@(#) NodeGiDRecorder.C, revA"

#include <NodeGiDRecorder.h>
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <FE_Datastore.h>
#include <FEM_ObjectBroker.h>

#include <string.h>
#include <stdlib.h>

#define RECORDER_TAGS_NodeGiDRecorder 21

NodeGiDRecorder::NodeGiDRecorder()
:Recorder(RECORDER_TAGS_NodeGiDRecorder),
 theDofs(0), theNodalTags(0), theNodes(0), response(0), 
 theDomain(0), theOutputHandler(0),
 echoTimeFlag(true), dataFlag(0), 
 deltaT(0), nextTimeStampToRecord(0.0), 
 initializationDone(false), numValidNodes(0), addColumnInfo(0)
{

}

NodeGiDRecorder::NodeGiDRecorder(const ID &dofs, 
			   const ID *nodes, 
			   const char *dataToStore,
			   Domain &theDom,
			   OPS_Stream &theOutputHandler,
			   double dT,
			   bool timeFlag)
:Recorder(RECORDER_TAGS_NodeGiDRecorder),
 theDofs(0), theNodalTags(0), theNodes(0), response(0), 
 theDomain(&theDom), theOutputHandler(&theOutputHandler),
 echoTimeFlag(timeFlag), dataFlag(0), 
 deltaT(dT), nextTimeStampToRecord(0.0), 
 initializationDone(false), numValidNodes(0), addColumnInfo(0)
{

  //
  // store copy of dof's to be recorder, verifying dof are valid, i.e. >= 0
  //

  int numDOF = dofs.Size();

  if (numDOF != 0) {
    
    theDofs = new ID(numDOF);
    
    int count = 0;
    int i;
    for (i=0; i<numDOF; i++) {
      int dof = dofs(i);
      if (dof >= 0) {
	(*theDofs)[count] = dof;
	count++;
      } else {
	opserr << "NodeGiDRecorder::NodeGiDRecorder - invalid dof  " << dof;
	opserr << " will be ignored\n";
      }
    }
  }

  // 
  // create memory to hold nodal ID's
  //

  if (nodes != 0) {
    int numNode = nodes->Size();
    if (numNode != 0) {
      theNodalTags = new ID(*nodes);
      if (theNodalTags == 0 || theNodalTags->Size() != nodes->Size()) {
	opserr << "NodeGiDRecorder::NodeGiDRecorder - out of memory\n";
      }
    }
  } 

  //
  // set the data flag used as a switch to get the response in a record
  //

  if (dataToStore == 0 || (strcmp(dataToStore, "disp") == 0)) {
    dataFlag = 0;
  } else if ((strcmp(dataToStore, "vel") == 0)) {
    dataFlag = 1;
  } else if ((strcmp(dataToStore, "accel") == 0)) {
    dataFlag = 2;
  } else if ((strcmp(dataToStore, "incrDisp") == 0)) {
    dataFlag = 3;
  } else if ((strcmp(dataToStore, "incrDeltaDisp") == 0)) {
    dataFlag = 4;
  } else if ((strcmp(dataToStore, "unbalance") == 0)) {
    dataFlag = 5;
  } else if ((strcmp(dataToStore, "unbalanceInclInertia") == 0) ||
	     (strcmp(dataToStore, "unbalanceIncInertia") == 0) ||
	     (strcmp(dataToStore, "unbalanceIncludingInertia") == 0))  {
    dataFlag = 6;
  } else if ((strcmp(dataToStore, "reaction") == 0)) {
    dataFlag = 7;
  } else if (((strcmp(dataToStore, "reactionIncInertia") == 0))
	     || ((strcmp(dataToStore, "reactionInclInertia") == 0))
	     || ((strcmp(dataToStore, "reactionIncludingInertia") == 0))) {
    dataFlag = 8;
  } else if (((strcmp(dataToStore, "rayleighForces") == 0))
	     || ((strcmp(dataToStore, "rayleighDampingForces") == 0))) {
    dataFlag = 9;
  } else if ((strncmp(dataToStore, "eigen",5) == 0)) {
    int mode = atoi(&(dataToStore[5]));
    if (mode > 0)
      dataFlag = 10 + mode;
    else
      dataFlag = 10;

  } else {
    dataFlag = 10;
    opserr << "NodeGiDRecorder::NodeGiDRecorder - dataToStore " << dataToStore;
    opserr << "not recognized (disp, vel, accel, incrDisp, incrDeltaDisp)\n";
  }
}


NodeGiDRecorder::~NodeGiDRecorder()
{
  if (theOutputHandler != 0) {
    theOutputHandler->endTag(); // Data
    delete theOutputHandler;
  }

  if (theDofs != 0)
    delete theDofs;

  if (theNodalTags != 0)
    delete theNodalTags;

  if (theNodes != 0)
    delete [] theNodes;
}

int 
NodeGiDRecorder::record(int commitTag, double timeStamp)
{
  if (theDomain == 0 || theDofs == 0) {
    return 0;
  }

  if (theOutputHandler == 0) {
    opserr << "NodeGiDRecorder::record() - no DataOutputHandler has been set\n";
    return -1;
  }

  if (initializationDone != true) 
    if (this->initialize() != 0) {
      opserr << "NodeGiDRecorder::record() - failed in initialize()\n";
      return -1;
    }

  int numDOF = theDofs->Size();
  
  if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {

    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;

    //
    // if need nodal reactions get the domain to calculate them
    // before we iterate over the nodes
    //

    if (dataFlag == 7)
      theDomain->calculateNodalReactions(0);
    else if (dataFlag == 8)
      theDomain->calculateNodalReactions(1);
    if (dataFlag == 9)
      theDomain->calculateNodalReactions(2);
    //
    // add time information if requested
    //

    int timeOffset = 0;
    //if (echoTimeFlag == true) {
    //  timeOffset = 1;
    //  response(0) = timeStamp;
    //}
    //
    // now we go get the responses from the nodes & place them in disp vector
    //

	char dataType[10];
	char outputData[100];
    if (dataFlag == 0) {
      strcpy(dataType,"Disp");
    } else if (dataFlag == 1) {
      strcpy(dataType,"Vel");
    } else if (dataFlag == 2) {
      strcpy(dataType,"Accel");
    } else if (dataFlag == 3) {
      strcpy(dataType,"dD");
    } else if (dataFlag == 4) {
      strcpy(dataType,"ddD");
    } else if (dataFlag == 5) {
      strcpy(dataType,"U");
    } else if (dataFlag == 6) {
      strcpy(dataType,"U");
    } else if (dataFlag == 7) {
      strcpy(dataType,"R");
    } else if (dataFlag == 8) {
      strcpy(dataType,"R");
    } else if (dataFlag > 10) {
      sprintf(dataType,"E%d", dataFlag-10);
    } else
      strcpy(dataType,"Unknown");

	if (numDOF == 2) {
	  // Result " Nodal Disp ," Analysis"      1.00000 Vector OnNodes
	  sprintf(outputData,"Result \" Nodal %s\" ,\" Analysis\" %12.5f Vector OnNodes\n",dataType,timeStamp);
	  theOutputHandler->write(outputData,strlen(outputData));
	  // ComponentNames "X-Disp"  "Y-Disp"
	  sprintf(outputData, "ComponentNames \"X-%s\"  \"Y-%s\"\n",dataType,dataType);
	  theOutputHandler->write(outputData,strlen(outputData));
	  // Values
	  theOutputHandler->write("Values\n",8);

	} else if (numDOF == 3) {
      // Result " Nodal Disp ," Analysis"      1.00000 Vector OnNodes
      sprintf(outputData,"Result \" Nodal %s\" ,\" Analysis\" %12.5f Vector OnNodes\n",dataType,timeStamp);
      theOutputHandler->write(outputData,strlen(outputData));
      //ComponentNames "X-Disp"  "Y-Disp"
      sprintf(outputData, "ComponentNames \"X-%s\"  \"Y-%s\"  \"Z-%s\"\n",dataType,dataType,dataType);
      theOutputHandler->write(outputData,strlen(outputData));
	  // Values
	  theOutputHandler->write("Values\n",8);
	}

    int cnt;

    if (dataFlag != 10) {

      for (int i=0; i<numValidNodes; i++) {

	cnt = 0;//i*numDOF + timeOffset; 
	response(0) = theNodes[i]->getTag();
	cnt ++;
	Node *theNode = theNodes[i];
	if (dataFlag == 0) {
	  const Vector &theResponse = theNode->getTrialDisp();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	  response(cnt) = theResponse(dof);
	    }
	    else {
	  response(cnt) = 0.0;
	    }
	    cnt++;
	  }
	}else if (dataFlag == 1) {
	  const Vector &theResponse = theNode->getTrialVel();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof);    
	    } else 
	      response(cnt) = 0.0;    
	    
	    cnt++;
	  }
	} else if (dataFlag == 2) {
	  const Vector &theResponse = theNode->getTrialAccel();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof);    
	    } else 
	      response(cnt) = 0.0;    
	    
	    cnt++;
	  }
	} else if (dataFlag == 3) {
	  const Vector &theResponse = theNode->getIncrDisp();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof);    
	    } else 
	      response(cnt) = 0.0;    
	    
	    cnt++;
	  }
	} else if (dataFlag == 4) {
	  const Vector &theResponse = theNode->getIncrDeltaDisp();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof);    
	    } else 
	      response(cnt) = 0.0;    
	    
	    cnt++;
	  }
	} else if (dataFlag == 5) {
	  const Vector &theResponse = theNode->getUnbalancedLoad();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof);
	    } else 
	      response(cnt) = 0.0;
	    
	    cnt++;
	  }
	  
	} else if (dataFlag == 6) {
	  const Vector &theResponse = theNode->getUnbalancedLoadIncInertia();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof);
	    } else 
	      response(cnt) = 0.0;
	    
	    cnt++;
	  }
	  
	  
	} else if (dataFlag == 7 || dataFlag == 8 || dataFlag == 9) {
	  const Vector &theResponse = theNode->getReaction();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof);
	    } else 
	      response(cnt) = 0.0;
	    cnt++;
	  }
	  
	} else if (10 <= dataFlag  && dataFlag < 1000) {
	  int mode = dataFlag - 10;
	  int column = mode - 1;
	  
	  const Matrix &theEigenvectors = theNode->getEigenvectors();
	  if (theEigenvectors.noCols() > column) {
	    int noRows = theEigenvectors.noRows();
	    for (int j=0; j<numDOF; j++) {
	      int dof = (*theDofs)(j);
	      if (noRows > dof) {
		response(cnt) = theEigenvectors(dof,column);
	      } else 
		response(cnt) = 0.0;
	      cnt++;		
	    }
	  }
	}
	
	else {
	  // unknown response
	  for (int j=0; j<numDOF; j++) {
	    response(cnt) = 0.0;
	  }
	}	
	
	//this->initStep();
	// insert the data into the database
    theOutputHandler->write(response);
      }
	  // End Values
	  theOutputHandler->write("End Values\n",12);
    } else {
	  opserr << "NodeGiDRecorder::record() --- wrong response type! " << endln;
    }
  }

  return 0;
}

int 
NodeGiDRecorder::setDomain(Domain &theDom)
{
  theDomain = &theDom;
  return 0;
}

int 
NodeGiDRecorder::sendSelf(int commitTag, Channel &theChannel)
{
  addColumnInfo = 1;

  if (theChannel.isDatastore() == 1) {
    opserr << "NodeGiDRecorder::sendSelf() - does not send data to a datastore\n";
    return -1;
  }

  static ID idData(7); 
  idData.Zero();
  if (theDofs != 0)
    idData(0) = theDofs->Size();
  if (theNodalTags != 0)
    idData(1) = theNodalTags->Size();
  if (theOutputHandler != 0) {
    idData(2) = theOutputHandler->getClassTag();
  }
  
  if (echoTimeFlag == true)
    idData(3) = 1;
  else
    idData(3) = 0;

  idData(4) = dataFlag;
  idData(5) = 0;

  idData(6) = this->getTag();

  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "NodeGiDRecorder::sendSelf() - failed to send idData\n";
    return -1;
  }

  if (theDofs != 0) 
    if (theChannel.sendID(0, commitTag, *theDofs) < 0) {
      opserr << "NodeGiDRecorder::sendSelf() - failed to send dof id's\n";
      return -1;
    }

  if (theNodalTags != 0)
    if (theChannel.sendID(0, commitTag, *theNodalTags) < 0) {
      opserr << "NodeGiDRecorder::sendSelf() - failed to send nodal tags\n";
      return -1;
    }

  static Vector data(2);
  data(0) = deltaT;
  data(1) = nextTimeStampToRecord;
  if (theChannel.sendVector(0, commitTag, data) < 0) {
    opserr << "NodeGiDRecorder::sendSelf() - failed to send data\n";
    return -1;
  }

  if (theOutputHandler->sendSelf(commitTag, theChannel) < 0) {
    opserr << "NodeGiDRecorder::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }

  return 0;
}

int 
NodeGiDRecorder::recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  addColumnInfo = 1;

  if (theChannel.isDatastore() == 1) {
    opserr << "NodeGiDRecorder::sendSelf() - does not send data to a datastore\n";
    return -1;
  }

  static ID idData(7); 
  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "NodeGiDRecorder::recvSelf() - failed to send idData\n";
    return -1;
  }


  int numDOFs = idData(0);
  int numNodes = idData(1);

  this->setTag(idData(6));

  if (idData(3) == 1)
    echoTimeFlag = true;
  else
    echoTimeFlag = false;    

  dataFlag = idData(4);
  //sensitivity = idData(5);

  //
  // get the DOF ID data
  //

  if (theDofs == 0 || theDofs->Size() != numDOFs) {
    if (theDofs != 0)
      delete theDofs;

    if (numDOFs != 0) {
      theDofs = new ID(numDOFs);
      if (theDofs == 0 || theDofs->Size() != numDOFs) {
	opserr << "NodeGiDRecorder::recvSelf() - out of memory\n";
	return -1;
      }	
    }
  }
  if (theDofs != 0)
    if (theChannel.recvID(0, commitTag, *theDofs) < 0) {
      opserr << "NodeGiDRecorder::recvSelf() - failed to recv dof data\n";
      return -1;
    } 

  //
  // get the NODAL tag data
  //

  if (theNodalTags == 0 || theNodalTags->Size() != numNodes) {
    if (theNodalTags != 0)
      delete theNodalTags;

    if (numNodes != 0) {
      theNodalTags = new ID(numNodes);
      if (theNodalTags == 0 || theNodalTags->Size() != numNodes) {
	opserr << "NodeGiDRecorder::recvSelf() - out of memory\n";
	return -1;
      }	
    }
  }
  if (theNodalTags != 0)
    if (theChannel.recvID(0, commitTag, *theNodalTags) < 0) {
      opserr << "NodeGiDRecorder::recvSelf() - failed to recv dof data\n";
      return -1;
    } 


  static Vector data(2);
  if (theChannel.recvVector(0, commitTag, data) < 0) {
    opserr << "NodeGiDRecorder::sendSelf() - failed to receive data\n";
    return -1;
  }
  deltaT = data(0);
  nextTimeStampToRecord = data(1);


  if (theOutputHandler != 0)
    delete theOutputHandler;

  theOutputHandler = theBroker.getPtrNewStream(idData(2));
  if (theOutputHandler == 0) {
    opserr << "NodeGiDRecorder::sendSelf() - failed to get a data output handler\n";
    return -1;
  }

  if (theOutputHandler->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "NodeGiDRecorder::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }

  return 0;
}

int
NodeGiDRecorder::domainChanged(void)
{
  return 0;
}

int
NodeGiDRecorder::initialize(void)
{
  if (theDofs == 0 || theDomain == 0) {
    opserr << "NodeGiDRecorder::initialize() - either nodes, dofs or domain has not been set\n";
    return -1;
  }

  //
  // create & set nodal array pointer
  //

  if (theNodes != 0) 
    delete [] theNodes;

  numValidNodes = 0;

  if (theNodalTags != 0) {

    int numNode = theNodalTags->Size();
    theNodes = new Node *[numNode];
    if (theNodes == 0) {
      opserr << "NodeGiDRecorder::domainChanged - out of memory\n";
      return -1;
    }

    for (int i=0; i<numNode; i++) {
      int nodeTag = (*theNodalTags)(i);
      Node *theNode = theDomain->getNode(nodeTag);
      if (theNode != 0) {
	theNodes[numValidNodes] = theNode;
	numValidNodes++;
      }
    }
  } else {

    int numNodes = theDomain->getNumNodes();
    theNodes = new Node *[numNodes];
    opserr << "NodeGiDRecorder::initialize - numNodes: " << numNodes << endln;
    if (theNodes == 0) {
      opserr << "NodeGiDRecorder::domainChanged - out of memory\n";
      return -1;
    }
    NodeIter &theDomainNodes = theDomain->getNodes();
    Node *theNode;
    numValidNodes = 0;
    while (((theNode = theDomainNodes()) != 0) && (numValidNodes < numNodes)) {
      theNodes[numValidNodes] = theNode;
      numValidNodes++;
    }
  }

  //
  // resize the response vector
  //

  int timeOffset = 0;
  if (echoTimeFlag == true)
    timeOffset = 1;

  int numValidResponse = theDofs->Size()+1; //+ timeOffset;numValidNodes*
  response.resize(numValidResponse);
  response.Zero();

  //ID orderResponse(numValidResponse);

  //
  // need to create the data description, i.e. what each column of data is
  //
  
  //char outputData[32];
  //char dataType[10];
  //
  //if (dataFlag == 0) {
  //  strcpy(dataType,"Disp");
  //} else if (dataFlag == 1) {
  //  strcpy(dataType,"Vel");
  //} else if (dataFlag == 2) {
  //  strcpy(dataType,"Accel");
  //} else if (dataFlag == 3) {
  //  strcpy(dataType,"dD");
  //} else if (dataFlag == 4) {
  //  strcpy(dataType,"ddD");
  //} else if (dataFlag == 5) {
  //  strcpy(dataType,"U");
  //} else if (dataFlag == 6) {
  //  strcpy(dataType,"U");
  //} else if (dataFlag == 7) {
  //  strcpy(dataType,"R");
  //} else if (dataFlag == 8) {
  //  strcpy(dataType,"R");
  //} else if (dataFlag > 10) {
  //  sprintf(dataType,"E%d", dataFlag-10);
  //} else
  //  strcpy(dataType,"Unknown");
  //
  //int numDOF = theDofs->Size();
  
  // write out info to handler if parallel execution
  //  

  //ID xmlOrder(numValidNodes);
  //
  //if (echoTimeFlag == true)  
  //  xmlOrder.resize(numValidNodes+1);
  //
  //if (theNodalTags != 0 && addColumnInfo == 1) {
  //
  //  int numNode = theNodalTags->Size();
  //  int count = 0;
  //  int nodeCount = 0;
  //
  //  if (echoTimeFlag == true)  {
  //    orderResponse(count++) = 0;
  //    xmlOrder(nodeCount++) = 0;
  //  }
  //  
  //  for (int i=0; i<numNode; i++) {
  //    int nodeTag = (*theNodalTags)(i);
  //    Node *theNode = theDomain->getNode(nodeTag);
  //    if (theNode != 0) {
  //xmlOrder(nodeCount++) = i+1;
  //for (int j=0; j<numDOF; j++)
  //  orderResponse(count++) = i+1;
  //    }
  //  }
  //
  //  theOutputHandler->setOrder(xmlOrder);
  //}

  //char nodeCrdData[20];
  //sprintf(nodeCrdData,"coord");

  //if (echoTimeFlag == true) {
  //  if (theNodalTags != 0 && addColumnInfo == 1) {
  //    theOutputHandler->tag("TimeOutput");
  //    theOutputHandler->tag("ResponseType", "time");
  //    theOutputHandler->endTag();
  //  }
  //}
  //
  //for (int i=0; i<numValidNodes; i++) {
  //  int nodeTag = theNodes[i]->getTag();
  //  const Vector &nodeCrd = theNodes[i]->getCrds();
  //  int numCoord = nodeCrd.Size();
  //
  //
  //  theOutputHandler->tag("NodeOutput");
  //  theOutputHandler->attr("nodeTag", nodeTag);
  //
  //  for (int j=0; j<3; j++) {
  //    sprintf(nodeCrdData,"coord%d",j+1);
  //    if (j < numCoord)
	//theOutputHandler->attr(nodeCrdData, nodeCrd(j));      
  //    else
	//theOutputHandler->attr(nodeCrdData, 0.0);      
  //  }
  //
  //  for (int k=0; k<theDofs->Size(); k++) {
  //    sprintf(outputData, "%s%d", dataType, k+1);
  //    theOutputHandler->tag("ResponseType",outputData);
  //  }
  //
  //  theOutputHandler->endTag();
  //}

  //if (theNodalTags != 0 && addColumnInfo == 1) {
  //  theOutputHandler->setOrder(orderResponse);
  //}

  //theOutputHandler->tag("Data");
  initializationDone = true;

  return 0;
}