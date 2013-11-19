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
                                                                        
// $Revision: 1.34 $
// $Date: 2009/04/30 23:25:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/ElementGiDRecorder.cpp,v $
                                                                        
// Written:
// Created:
//
// Description: This file contains the class implementations of ElementGiDRecorder.
//
// What: "@(#) ElementGiDRecorder.C, revA"

#include <ElementGiDRecorder.h>
#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <string.h>
#include <Response.h>
#include <FE_Datastore.h>
#include <OPS_Globals.h>
#include <Message.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <packages.h>
#include <elementAPI.h>

ElementGiDRecorder::ElementGiDRecorder()
:Recorder(RECORDER_TAGS_ElementGiDRecorder),
 numEle(0), eleID(0), theResponses(0), 
 theDomain(0), theOutputHandler(0),
 echoTimeFlag(true), deltaT(0), nextTimeStampToRecord(0.0),
 initializationDone(false), responseArgs(0), numArgs(0), addColumnInfo(0),
 hasLinear(false), hasTri3(false), hasQuad4(false), hasQuad8(false), hasQuad9(false), hasBrick(false)
{

}

ElementGiDRecorder::ElementGiDRecorder(ID &eleIDs,
				 const char **argv, 
				 int argc,
				 bool echoTime, 
				 Domain &theDom, 
				 OPS_Stream &theOutputHandler,
				 double dT)
:Recorder(RECORDER_TAGS_ElementGiDRecorder),
 numEle(0), eleID(0), theResponses(0), 
 theDomain(&theDom), theOutputHandler(&theOutputHandler),
 echoTimeFlag(echoTime), deltaT(dT), nextTimeStampToRecord(0.0),
 initializationDone(false), responseArgs(0), numArgs(0), addColumnInfo(0),
 hasLinear(false), hasTri3(false), hasQuad4(false), hasQuad8(false), hasQuad9(false), hasBrick(false)
{

  if (eleIDs != 0) {
    numEle = eleIDs.Size();
    eleID = new ID(eleIDs);
    if (eleID == 0 || eleID->Size() != numEle)
      opserr << "ElementGiDRecorder::ElementGiDRecorder() - out of memory\n";
  } 

  //
  // create a copy of the response request
  //

  responseArgs = new char *[argc];
  if (responseArgs == 0) {
    opserr << "ElementGiDRecorder::ElementGiDRecorder() - out of memory\n";
    numEle = 0;
  }
  
  for (int i=0; i<argc; i++) {
    responseArgs[i] = new char[strlen(argv[i])+1];
    if (responseArgs[i] == 0) {
      delete [] responseArgs;
      opserr << "ElementGiDRecorder::ElementGiDRecorder() - out of memory\n";
      numEle = 0;
    }
    strcpy(responseArgs[i], argv[i]);
  }
  
  numArgs = argc;
}

ElementGiDRecorder::~ElementGiDRecorder()
{
  theOutputHandler->endTag(); // Data

  if (theOutputHandler != 0)
    delete theOutputHandler;

  //
  // invoke the destructor on the response objects
  //

  if (eleID != 0)
    delete eleID;

  if (theResponses != 0) {
    for (int i = 0; i < numEle; i++)
      delete theResponses[i];
    delete [] theResponses;
  }

  //if (data != 0)
  //  delete data;
  
  // 
  // invoke destructor on response args
  //
  
  for (int i=0; i<numArgs; i++)
    delete [] responseArgs[i];
  delete [] responseArgs;

}

int 
ElementGiDRecorder::record(int commitTag, double timeStamp)
{
  // 
  // check that initialization has been done
  //

  if (initializationDone == false) {
    if (this->initialize() != 0) {
      opserr << "ElementGiDRecorder::record() - failed to initialize\n";
      return -1;
    }
  }
  
  int result = 0;
  if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {
	stepN ++;

    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;

    int loc = 0;
    //if (echoTimeFlag == true) 
    //  (*data)(loc++) = timeStamp;
    //int ndf = OPS_GetNDF(); 
	//int ndm = OPS_GetNDM();
	char outputData[100];
	// **** Linear Elements - 2 Nodes
    if (hasLinear == 1) {
      if (stricmp(responseArgs[numArgs-1], "force") == 0 && numArgs >1 ) { // force on section
	    // Result " Name of Results " "Analysis"      1.00000 Vector OnNodes
        sprintf(outputData, "Result \"Element_%s\" \"Loading_Analysis\"\t%i", responseArgs[numArgs-1], stepN);
        theOutputHandler->write(outputData,80);
	    theOutputHandler->write("\tMatrix OnGaussPoints \"Linear_GuassPoint_Set\"\n",80);
        // ComponentNames 
        theOutputHandler->write(" ComponentNames \"Mz\"  \"P\"  \"Vy\"  \"My\"  \"Vz\"  \"T\"\n",80);
	  }
	  else if (stricmp(responseArgs[numArgs-1], "force") == 0 && numArgs == 1) { // force on element node
		// Result " Name of Results " "Analysis"      1.00000 Vector OnNodes
		sprintf(outputData, "Result \"Element_%s\" \"Loading_Analysis\"\t%i", responseArgs[numArgs-1], stepN);
		theOutputHandler->write(outputData,80);
		theOutputHandler->write("\tVector OnNodes \n",80);
		// ComponentNames 
		if (OPS_GetNDM() == 2) {
		  theOutputHandler->write(" ComponentNames \"Mz\"  \"P\"  \"Vy\"\n",80);
		}
		else if (OPS_GetNDM() == 3) {
		  theOutputHandler->write(" ComponentNames \"Mz\"  \"P\"  \"Vy\"  \"My\"  \"Vz\"  \"T\"\n",80);
		}

	  }
	  else if (stricmp(responseArgs[numArgs-1], "deformation") == 0 ) { // deformation on section
		// Result " Name of Results " "Analysis"      1.00000 Vector OnNodes
        sprintf(outputData, "Result \"Element_%s\" \"Loading_Analysis\"\t%i", responseArgs[numArgs-1], stepN);
        theOutputHandler->write(outputData,80);
	    theOutputHandler->write("\tMatrix OnGaussPoints \"Linear_GuassPoint_Set\"\n",80);
        // ComponentNames 
        theOutputHandler->write(" ComponentNames \"kappaZ\"  \"eps\"  \"gammaY\"  \"kappaY\"  \"gammaZ\"  \"theta\"\n",80);
	  }

	  // Values
	  theOutputHandler->write(" Values\n",8);
	    
      // for each element if responses exist, put them in response vector
      for (int i=0; i< numEle; i++) {
	    sprintf(outputData, "%i\t", (*eleID)(i));
	    theOutputHandler->write(outputData,20);
        if (theResponses[i] != 0) {
	  // ask the element for the response
	  int res, dateLength;
	  loc = 0;
	  theOutputHandler->setPrecision(12);
	  if (( res = theResponses[i]->getResponse()) < 0) {
	    result += res;
	    dateLength = 0;
	  }else {
	    Information &eleInfo = theResponses[i]->getInformation();
	    const Vector &eleData = eleInfo.getData();
	    // for each integration points
	    int nIP = theOutputHandler->getEleGPs();
	    for (int k=0; k<nIP; k++) {
	  dateLength = eleData.Size()/nIP;
	  for (int j=0; j<dateLength; j++) {
	    double temp = eleData(k*dateLength+j);
	    theOutputHandler->write(temp);
	  }
	  // send the response vector to the output handler for o/p
	  //theOutputHandler->write(*data);
	  theOutputHandler->write("\n",2);
	    }
	  }
        }
      }
	  // End Values
	  theOutputHandler->write("End Values\n\n",14);
    }
    // **** Quadrilateral Elements - 4 Nodes
    if (hasQuad4 == 1) {
	  // Print HEADER
	  if (stricmp(responseArgs[0], "material") == 0  ) { // integration of material point
        opserr << "4 Points quad element integrPoint record has not been implemented,yet." << endln;
	  }
	  else if (stricmp(responseArgs[numArgs-1], "stresses") == 0 ) { // stress on one gauss point
		// Result " Name of Results " "Analysis"      1.00000 Vector OnNodes
        sprintf(outputData, "Result \"Element_%s\" \"Loading_Analysis\"\t%i", responseArgs[numArgs-1], stepN);
        theOutputHandler->write(outputData,80);
	    theOutputHandler->write("\tVector OnGaussPoints \"Quadrilateral_4_GuassPoint_Set\"\n",80);
        // ComponentNames 
        theOutputHandler->write(" ComponentNames \"Sxx\"  \"Syy\"  \"Sxy\"\n",80);
	  }
	  else if (stricmp(responseArgs[numArgs-1], "force") == 0) {
		opserr << "4 Points quad element force record has not been implemented,yet." << endln;
	  }
	  // Values
	  theOutputHandler->write(" Values\n",8);
	    
      // for each element if responses exist, put them in response vector
      for (int i=0; i< numEle; i++) {
	    sprintf(outputData, "%i\t", (*eleID)(i));
	    theOutputHandler->write(outputData,20);
        if (theResponses[i] != 0) {
	  // ask the element for the response
	  int res, dateLength;
	  loc = 0;
	  theOutputHandler->setPrecision(12);
	  if (( res = theResponses[i]->getResponse()) < 0) {
	    result += res;
	    dateLength = 0;
	  }else {
	    Information &eleInfo = theResponses[i]->getInformation();
	    const Vector &eleData = eleInfo.getData();
	    // for each integration points
		dateLength = theOutputHandler->getEleGPs();
	    int nIP = eleData.Size()/dateLength;
	    for (int k=0; k<nIP; k++) {
	  for (int j=0; j<dateLength; j++) {
	    double temp = eleData(k*dateLength+j);
	    theOutputHandler->write(temp);
	  }
	  // send the response vector to the output handler for o/p
	  //theOutputHandler->write(*data);
	  theOutputHandler->write("\n",2);
	    }
	  }
        }
      }
	  // End Values
	  theOutputHandler->write("End Values\n\n",14);
    }
    // **** Triangular Elements - 3 Nodes
    if (hasTri3 == 1) {
      // Print HEADER
	  if (stricmp(responseArgs[0], "material") == 0  ) { // integration of material point
	    opserr << "3 Points element integrPoint record has not been implemented,yet." << endln;
	  }
	  else if (stricmp(responseArgs[numArgs-1], "stresses") == 0 ) { // stress on one gauss point
		// Result " Name of Results " "Analysis"      1.00000 Vector OnNodes
        sprintf(outputData, "Result \"Element_%s\" \"Loading_Analysis\"\t%i", responseArgs[numArgs-1], stepN);
        theOutputHandler->write(outputData,80);
	    theOutputHandler->write("\tVector OnGaussPoints \"Triangle_GuassPoint_Set\"\n",80);
        // ComponentNames 
        theOutputHandler->write(" ComponentNames \"Sxx\"  \"Syy\"  \"Sxy\"\n",80);
	  }
	  else if (stricmp(responseArgs[numArgs-1], "force") == 0) {
		opserr << "3 Points element force record has not been implemented,yet." << endln;
	  }
	  // Values
	  theOutputHandler->write(" Values\n",8);
	    
      // for each element if responses exist, put them in response vector
      for (int i=0; i< numEle; i++) {
	    sprintf(outputData, "%i\t", (*eleID)(i));
	    theOutputHandler->write(outputData,20);
        if (theResponses[i] != 0) {
	  // ask the element for the response
	  int res, dateLength;
	  loc = 0;
	  theOutputHandler->setPrecision(12);
	  if (( res = theResponses[i]->getResponse()) < 0) {
	    result += res;
	    dateLength = 0;
	  }else {
	    Information &eleInfo = theResponses[i]->getInformation();
	    const Vector &eleData = eleInfo.getData();
	    // for each integration points
		dateLength = theOutputHandler->getEleGPs();
		int nIP = eleData.Size()/dateLength;
	    for (int k=0; k<nIP; k++) {
	  for (int j=0; j<dateLength; j++) {
	    double temp = eleData(k*dateLength+j);
	    theOutputHandler->write(temp);
	  }
	  // send the response vector to the output handler for o/p
	  //theOutputHandler->write(*data);
	  theOutputHandler->write("\n",2);
	    }
	  }
        }
      }
	  // End Values
	  theOutputHandler->write("End Values\n\n",14);
    }
    // **** Quadrilateral Elements - 9 Nodes
    if (hasQuad9 == 1) {
	  // Print HEADER
	  if (stricmp(responseArgs[0], "material") == 0  ) { // integration of material point
		opserr << "9 Points quad element integrPoint record has not been implemented,yet." << endln;
	  }
	  else if (stricmp(responseArgs[numArgs-1], "stresses") == 0 ) { // stress on one gauss point
		// Result " Name of Results " "Analysis"      1.00000 Vector OnNodes
        sprintf(outputData, "Result \"Element_%s\" \"Loading_Analysis\"\t%i", responseArgs[numArgs-1], stepN);
        theOutputHandler->write(outputData,80);
	    theOutputHandler->write("\tMatrix OnGaussPoints \"Quadrilateral_9_GuassPoint_Set\"\n",90);
        // ComponentNames 
        theOutputHandler->write(" ComponentNames \"p11\"  \"p22\"  \"p12\"  \"m11\"  \"m22\"  \"m12\"  \"q1\"  \"q2\"\n",80);
	  }
	  else if (stricmp(responseArgs[numArgs-1], "force") == 0) {
		opserr << "9 Points quad element force record has not been implemented,yet." << endln;
	  }
	  // Values
	  theOutputHandler->write(" Values\n",8);
	    
      // for each element if responses exist, put them in response vector
      for (int i=0; i< numEle; i++) {
	    sprintf(outputData, "%i\t", (*eleID)(i));
	    theOutputHandler->write(outputData,20);
        if (theResponses[i] != 0) {
	  // ask the element for the response
	  int res, dateLength;
	  loc = 0;
	  theOutputHandler->setPrecision(12);
	  if (( res = theResponses[i]->getResponse()) < 0) {
	    result += res;
	    dateLength = 0;
	  }else {
	    Information &eleInfo = theResponses[i]->getInformation();
	    const Vector &eleData = eleInfo.getData();
	    // for each integration points
		dateLength = theOutputHandler->getEleGPs();
		int nIP = eleData.Size()/dateLength;
	    for (int k=0; k<nIP; k++) {
	  for (int j=0; j<dateLength; j++) {
	    double temp = eleData(k*dateLength+j);
	    theOutputHandler->write(temp);
	  }
	  // send the response vector to the output handler for o/p
	  //theOutputHandler->write(*data);
	  theOutputHandler->write("\n",2);
	    }
	  }
        }
      }
	  // End Values
	  theOutputHandler->write("End Values\n\n",14);
    }
    // **** Hexahedra Elements - 8 Nodes
    if (hasBrick == 1) {
	  // Print HEADER
	  if (stricmp(responseArgs[0], "material") == 0  ) { // integration of material point
        opserr << "8 Points brick element integrPoint record has not been implemented,yet." << endln;
      }
	  else if (stricmp(responseArgs[numArgs-1], "stresses") == 0 ) { // stress on one gauss point
	  	// Result " Name of Results " "Analysis"      1.00000 Vector OnNodes
	  	sprintf(outputData, "Result \"Element_%s\" \"Loading_Analysis\"\t%i", responseArgs[numArgs-1], stepN);
	  	theOutputHandler->write(outputData,80);
	  	theOutputHandler->write("\tMatrix OnGaussPoints \"Hexahedra_GuassPoint_Set\"\n",90);
	  	// ComponentNames 
	  	theOutputHandler->write(" ComponentNames \"S11\"  \"S22\"  \"S33\"  \"S12\"  \"S23\"  \"S31\"\n",80);
		
	  }
	  else if (stricmp(responseArgs[numArgs-1], "force") == 0) {
	  	opserr << "8 Points brick element force record has not been implemented,yet." << endln;
	  }

	  // Values
	  theOutputHandler->write(" Values\n",8);
	  
      // for each element if responses exist, put them in response vector
      for (int i=0; i< numEle; i++) {
	    sprintf(outputData, "%i\t", (*eleID)(i));
	    theOutputHandler->write(outputData,20);
        if (theResponses[i] != 0) {
	  // ask the element for the response
	  int res, dateLength;
	  loc = 0;
	  theOutputHandler->setPrecision(12);
	  if (( res = theResponses[i]->getResponse()) < 0) {
	    result += res;
	    dateLength = 0;
	  }else {
	    Information &eleInfo = theResponses[i]->getInformation();
	    const Vector &eleData = eleInfo.getData();
	    // for each integration points
		dateLength = theOutputHandler->getEleGPs();
		int nIP = eleData.Size()/dateLength;
	    for (int k=0; k<nIP; k++) {
	  for (int j=0; j<dateLength; j++) {
	    double temp = eleData(k*dateLength+j);
	    theOutputHandler->write(temp);
	  }
	  // send the response vector to the output handler for o/p
	  //theOutputHandler->write(*data);
	  theOutputHandler->write("\n",2);
	    }
	  }
        }
      }
	  // End Values
	  theOutputHandler->write("End Values\n\n",14);
    }
  }

  // successfully completion - return 0
  return result;
}

int
ElementGiDRecorder::restart(void)
{
  //if (data != 0)
  //  data->Zero();
  return 0;
}

int 
ElementGiDRecorder::setDomain(Domain &theDom)
{
  theDomain = &theDom;
  return 0;
}

int
ElementGiDRecorder::sendSelf(int commitTag, Channel &theChannel)
{
  addColumnInfo = 1;

  if (theChannel.isDatastore() == 1) {
    opserr << "ElementGiDRecorder::sendSelf() - does not send data to a datastore\n";
    return -1;
  }

  //
  // into an ID, place & send (*eleID) size, numArgs and length of all responseArgs
  //

  static ID idData(6);
  if (eleID != 0)
    idData(0) = eleID->Size();
  else
    idData(0) = 0;

  idData(1) = numArgs;

  int msgLength = 0;
  for (int i=0; i<numArgs; i++) 
    msgLength += strlen(responseArgs[i])+1;

  idData(2) = msgLength;

  if (theOutputHandler != 0) {
    idData(3) = theOutputHandler->getClassTag();
  } else 
    idData(3) = 0;

  if (echoTimeFlag == true)
    idData(4) = 1;
  else
    idData(4) = 0;


  idData(5) = this->getTag();

  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "ElementGiDRecorder::sendSelf() - failed to send idData\n";
    return -1;
  }

  static Vector dData(1);
  dData(1) = deltaT;
  if (theChannel.sendVector(0, commitTag, dData) < 0) {
    opserr << "ElementGiDRecorder::sendSelf() - failed to send dData\n";
    return -1;
  }
  
  //
  // send the eleID
  //

  if (eleID != 0)
    if (theChannel.sendID(0, commitTag, *eleID) < 0) {
      opserr << "ElementGiDRecorder::sendSelf() - failed to send idData\n";
      return -1;
    }

  //
  // create a single char array holding all strings
  //    will use string terminating character to differentiate strings on other side
  //

  if (msgLength ==  0) {
    opserr << "ElementGiDRecorder::sendSelf() - no data to send!!\n";
    return -1;
  }

  char *allResponseArgs = new char[msgLength];
  if (allResponseArgs == 0) {
    opserr << "ElementGiDRecorder::sendSelf() - out of memory\n";
    return -1;
  }

  char *currentLoc = allResponseArgs;
  for (int j=0; j<numArgs; j++) {
    strcpy(currentLoc, responseArgs[j]);
    currentLoc += strlen(responseArgs[j]);
    currentLoc++;
  }

  //
  // send this single char array
  //


  Message theMessage(allResponseArgs, msgLength);
  if (theChannel.sendMsg(0, commitTag, theMessage) < 0) {
    opserr << "ElementGiDRecorder::sendSelf() - failed to send message\n";
    return -1;
  }

  //
  // invoke sendSelf() on the output handler
  //

  if (theOutputHandler == 0 || theOutputHandler->sendSelf(commitTag, theChannel) < 0) {
    opserr << "ElementGiDRecorder::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }

  //
  // clean up & return success
  //

  delete [] allResponseArgs;
  return 0;
}

int 
ElementGiDRecorder::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
  addColumnInfo = 1;

  if (theChannel.isDatastore() == 1) {
    opserr << "ElementGiDRecorder::recvSelf() - does not recv data to a datastore\n";
    return -1;
  }

  if (responseArgs != 0) {
    for (int i=0; i<numArgs; i++)
      delete [] responseArgs[i];
  
    delete [] responseArgs;
  }

  //
  // into an ID of size 2 recv eleID size and length of all responseArgs
  //

  static ID idData(6);
  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "ElementGiDRecorder::recvSelf() - failed to recv idData\n";
    return -1;
  }

  int eleSize = idData(0);
  numArgs = idData(1);
  int msgLength = idData(2);

  this->setTag(idData(5));

  if (idData(4) == 1)
    echoTimeFlag = true;
  else
    echoTimeFlag = false;    

  numEle = eleSize;

  static Vector dData(1);
  if (theChannel.recvVector(0, commitTag, dData) < 0) {
    opserr << "ElementGiDRecorder::sendSelf() - failed to send dData\n";
    return -1;
  }
  deltaT = dData(1);


  //
  // resize & recv the eleID
  //

  if (eleSize != 0) {
    eleID = new ID(eleSize);
    if (eleID == 0) {
      opserr << "ElementGiDRecorder::recvSelf() - failed to recv idData\n";
      return -1;
    }
    if (theChannel.recvID(0, commitTag, *eleID) < 0) {
      opserr << "ElementGiDRecorder::recvSelf() - failed to recv idData\n";
      return -1;
    }
  }

  //
  // recv the single char array of element response args
  //

  if (msgLength == 0) {
    opserr << "ElementGiDRecorder::recvSelf() - 0 sized string for responses\n";
    return -1;
  }

  char *allResponseArgs = new char[msgLength];
  if (allResponseArgs == 0) {
    opserr << "ElementGiDRecorder::recvSelf() - out of memory\n";
    return -1;
  }

  Message theMessage(allResponseArgs, msgLength);
  if (theChannel.recvMsg(0, commitTag, theMessage) < 0) {
    opserr << "ElementGiDRecorder::recvSelf() - failed to recv message\n";
    return -1;
  }

  //
  // now break this single array into many
  // 

  responseArgs = new char *[numArgs];
  if (responseArgs == 0) {
    opserr << "ElementGiDRecorder::recvSelf() - out of memory\n";
    return -1;
  }

  char *currentLoc = allResponseArgs;
  for (int j=0; j<numArgs; j++) {

    int argLength = strlen(currentLoc)+1;

    responseArgs[j] = new char[argLength];
    if (responseArgs[j] == 0) {
      opserr << "ElementGiDRecorder::recvSelf() - out of memory\n";
      return -1;
    }

    strcpy(responseArgs[j], currentLoc);
    currentLoc += argLength;
  }

  //
  // create a new handler object and invoke recvSelf() on it
  //

  if (theOutputHandler != 0)
    delete theOutputHandler;

  theOutputHandler = theBroker.getPtrNewStream(idData(3));
  if (theOutputHandler == 0) {
    opserr << "NodeRecorder::sendSelf() - failed to get a data output handler\n";
    return -1;
  }

  if (theOutputHandler->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "NodeRecorder::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }

  //
  // clean up & return success
  //

  delete [] allResponseArgs;
  return 0;
}

int 
ElementGiDRecorder::initialize(void)
{
  if (theDomain == 0)
    return 0;

  if (theResponses != 0) {
    for (int i = 0; i < numEle; i++)
      delete theResponses[i];
    delete [] theResponses;
  }

  int numDbColumns = 0;

  // Set the response objects:
  //   1. create an array of pointers for them
  //   2. iterate over the elements invoking setResponse() to get the new objects & determine size of data
  //

  int i =0;
  //ID xmlOrder(0,64);
  //ID responseOrder(0,64);

  if (eleID != 0) {

    //
    // if we have an eleID we know Response size so allocate Response holder & loop over & ask each element
    //

    //int eleCount = 0;
    //int responseCount = 0;

    if (echoTimeFlag == true && addColumnInfo == 1) {
    //  xmlOrder[0] = 0;
    //  responseOrder[0] = 0;
    //  eleCount = 1;
    //  responseCount =1;
    }

    // loop over ele & set Responses
    for (i=0; i<numEle; i++) {
      Element *theEle = theDomain->getElement((*eleID)(i));
      if (theEle != 0) {
	//xmlOrder[eleCount] = i+1;
	//eleCount++;

	// Cycle over Elements to understand what type of elements are there
	int tag = theEle->getTag();
    // Check type of Element with Number of Nodes
    int nNode = theEle->getNumExternalNodes();
    if (nNode == 2) {
      hasLinear = 1;
    } else if (nNode == 4) {
      hasQuad4 = 1;
    } else if (nNode == 3) {
      hasTri3 = 1;
    } else if (nNode == 9) {
      hasQuad9 = 1;
    } else if (nNode == 8) {
      const char *name = theEle->getClassType();
      if (strcmp(name,"Brick") == 0) {
        hasBrick = 1;
      } else {
        hasQuad8 = 1;
      }
    }
      }
    }

    //theOutputHandler->setOrder(xmlOrder);

    //
    // do time
    //

    if (echoTimeFlag == true) {
    //  theOutputHandler->tag("TimeOutput");
    //  theOutputHandler->tag("ResponseType", "time");
    //  theOutputHandler->endTag(); // TimeOutput
    //  numDbColumns += 1;
    }

    //
    // if we have an eleID we know Response size so allocate Response holder
	// & loop over & ask each element
    // allocate memory for Responses & set to 0
    theResponses = new Response *[numEle];
    if (theResponses == 0) {
      opserr << "ElementGiDRecorder::initialize() - out of memory\n";
      return -1;
    }

    for (int k=0; k<numEle; k++)
      theResponses[k] = 0;

    // loop over ele & set Responses
    for (i=0; i<numEle; i++) {
      Element *theEle = theDomain->getElement((*eleID)(i));
      if (theEle == 0) {
	theResponses[i] = 0;
      } else {
	theResponses[i] = theEle->setResponse((const char **)responseArgs, numArgs, *theOutputHandler);
	if (theResponses[i] != 0) {
	  // from the response type determine no of cols for each
	  Information &eleInfo = theResponses[i]->getInformation();
	  const Vector &eleData = eleInfo.getData();
	  int dataSize = eleData.Size();
	  numDbColumns = (numDbColumns > dataSize ? numDbColumns:dataSize); //neallee need revised
	  //if (addColumnInfo == 1) {
	  //for (int j=0; j<dataSize; j++)
	  //  responseOrder[responseCount++] = i+1;
	  //}
	}
      }
    }

    //theOutputHandler->setOrder(responseOrder);

  } else {

    //if (echoTimeFlag == true) {
    //  theOutputHandler->tag("TimeOutput");
    //  theOutputHandler->tag("ResponseType", "time");
    //  theOutputHandler->endTag(); // TimeOutput
    //  numDbColumns += 1;
    //}

    //
    // if no eleID we don't know response size so make initial guess & loop over & ask ele
    // if guess to small, we enlarge
    //

    // initial size & allocation
    int numResponse = 0;
    numEle = 12;
    theResponses = new Response *[numEle];

    if (theResponses == 0) {
      opserr << "ElementGiDRecorder::initialize() - out of memory\n";
      return -1;
    }

    for (int k=0; k<numEle; k++)
      theResponses[k] = 0;

    // loop over ele & set Responses
    ElementIter &theElements = theDomain->getElements();
    Element *theEle;

    while ((theEle = theElements()) != 0) {
      Response *theResponse = theEle->setResponse((const char **)responseArgs, numArgs, *theOutputHandler);
      if (theResponse != 0) {
	if (numResponse == numEle) {
	  // Why is this created locally and not used? -- MHS
	  Response **theNextResponses = new Response *[numEle*2];
	  if (theNextResponses != 0) {
	    for (i=0; i<numEle; i++)
	      theNextResponses[i] = theResponses[i];
	    for (int j=numEle; j<2*numEle; j++)
	      theNextResponses[j] = 0;
	  }
	  numEle = 2*numEle;
	  delete [] theNextResponses;
	}
	theResponses[numResponse] = theResponse;

	// from the response type determine no of cols for each
	Information &eleInfo = theResponses[numResponse]->getInformation();
	const Vector &eleData = eleInfo.getData();
	
	int temp = eleData.Size();
	numDbColumns = (numDbColumns > temp ? numDbColumns:temp); //neallee need revised

	numResponse++;

      }

	  // Cycle over Elements to understand what type of elements are there
	  int tag = theEle->getTag();
      // Check type of Element with Number of Nodes
      int nNode = theEle->getNumExternalNodes();
      if (nNode == 2) {
        hasLinear = 1;
	  //  LinearEle = theElement;
      } else if (nNode == 4) {
        hasQuad4 = 1;
      } else if (nNode == 3) {
        hasTri3 = 1;
      } else if (nNode == 9) {
        hasQuad9 = 1;
      } else if (nNode == 8) {
        const char *name = theEle->getClassType();
        if (strcmp(name,"Brick") == 0) {
          hasBrick = 1;
        } else {
          hasQuad8 = 1;
        }
      }
    }
    numEle = numResponse;
  }

  // create the vector to hold the data
  //data = new Vector(numDbColumns);

  //if (data == 0) {
  //  opserr << "ElementGiDRecorder::initialize() - out of memory\n";
  //  return -1;
  //}
  
  theOutputHandler->write("GiD Post Results File 1.0\n\n",40);
  int nIP = numDbColumns/theOutputHandler->getEleGPs();
  // **** Linear Elements - 2 Nodes
  if (hasLinear == 1) {
  	// Print HEADER
  	theOutputHandler->write("GaussPoints \"Linear_GuassPoint_Set\" ElemType Linear\n", 56); //defined for all meshtype
	theOutputHandler->write("Number Of Gauss Points: ", 25);
	theOutputHandler->write(nIP); //nIP
	theOutputHandler->write("\nNodes included\n",18); // for Lobbato Integeration assumption
	theOutputHandler->write("Natural Coordinates: internal\n", 34);
	theOutputHandler->write("End GaussPoints\n\n", 20);
  }
  // **** Quadrilateral Elements - 4 Nodes
  if (hasQuad4 == 1) {
  	// Print HEADER
  	theOutputHandler->write("GaussPoints \"Quadrilateral_4_GuassPoint_Set\" ElemType Quadrilateral\n", 70); //defined for all meshtype
	theOutputHandler->write("Number Of Gauss Points: ", 25);
	theOutputHandler->write(nIP); //nIP for quad = 4, for SSPquad = 1
	theOutputHandler->write("\nNatural Coordinates: internal\n", 34);
	theOutputHandler->write("End GaussPoints\n\n", 20);
  }
  // **** Triangular Elements - 3 Nodes
  if (hasTri3 == 1) {
  	// Print HEADER
  	theOutputHandler->write("GaussPoints \"Triangle_GuassPoint_Set\" ElemType Triangle\n", 60); //defined for all meshtype
	theOutputHandler->write("Number Of Gauss Points: ", 25);
	theOutputHandler->write(1); //nIP tri31 nIP = 1
	theOutputHandler->write("\nNatural Coordinates: internal\n", 34);
	theOutputHandler->write("End GaussPoints\n\n", 20);
  }
  // **** Quadrilateral Elements - 9 Nodes
  if (hasQuad9 == 1) {
  	// Print HEADER
  	theOutputHandler->write("GaussPoints \"Quadrilateral_9_GuassPoint_Set\" ElemType Quadrilateral\n", 72); //defined for all meshtype
	theOutputHandler->write("Number Of Gauss Points: ", 25);
	theOutputHandler->write(9); //nIP
	theOutputHandler->write("\nNatural Coordinates: internal\n", 34);
	theOutputHandler->write("End GaussPoints\n\n", 20);
  }
  // **** Hexahedra Elements - 8 Nodes
  if (hasBrick == 1) {
  	// Print HEADER
	theOutputHandler->write("GaussPoints \"Hexahedra_GuassPoint_Set\" ElemType Hexahedra\n", 62); //defined for all meshtype
	theOutputHandler->write("Number Of Gauss Points: ", 25);
	theOutputHandler->write(nIP); //nIP, for brick = 8
	theOutputHandler->write("\nNatural Coordinates: internal\n", 34);
	theOutputHandler->write("End GaussPoints\n\n", 20);
  }
  
  initializationDone = true;

  stepN =0;

  return 0;
}
