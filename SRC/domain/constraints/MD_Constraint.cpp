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
// $Date: 2008-07-03 18:00:48 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/domain/constraints/MD_Constraint.cpp,v $
                                                                        
                                                                        
// File: ~/domain/constraints//MD_Constraint.C
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the implementation of class MD_Constraint.
//
// The class MD_Constraint interface:
//

#include <MD_Constraint.h>

#include <stdlib.h>
#include <math.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <iostream>
#include <fstream>

static int numMDs = 0;
static int nextTag = 0;

// constructor for FEM_ObjectBroker			// Cenk
MD_Constraint::MD_Constraint()		
:MP_Constraint(CNSTRNT_TAG_MD_Constraint),thisDomain(0), 
 nodeRetained(0),nodeConstrained(0),commit_constraint(0), constraint(0),constrDOF(0),retainDOF(0), ConstrainedNode(0),
 dbTag1(0), dbTag2(0), dbTag3(0), global_vector(0), commit_global_vector(0), cnvg(0)
{
    
}

// general constructor for ModelBuilder
MD_Constraint::MD_Constraint(Domain *theDomain, int tag, int nodeRetain, int nodeConstr, 
			     ID &constrainedDOF, ID &retainedDOF, Vector &gv)
:MP_Constraint(CNSTRNT_TAG_MD_Constraint), 
 nodeRetained(nodeRetain), nodeConstrained(nodeConstr), 
 commit_constraint(0), constraint(0), constrDOF(0), retainDOF(0), ConstrainedNode(0), dbTag1(0), dbTag2(0), dbTag3(0), cnvg(0)
{
    
    constrDOF = new ID(constrainedDOF);
    retainDOF = new ID(retainedDOF);    
    global_vector = new Vector(gv);
    commit_global_vector = new Vector(gv);
    ConstrainedNode = theDomain->getNode(nodeConstrained);
    if(ConstrainedNode == NULL){
	opserr<<"MD_Constrained::MD_Constrained: nodeConstrained: ";
	opserr<<nodeConstrained<<"does not exist in model\n";
	exit(0);
    }
    
    if (constrDOF == 0 || constrainedDOF.Size() != constrDOF->Size() ||
	retainDOF == 0 || retainedDOF.Size() != retainDOF->Size()) { 
	opserr << "MD_Constraint::MD_Constraint - ran out of memory 1\n";
	exit(-1);
    }    
    
    constraint = new Matrix(9,9);
    if (constraint == 0 || constraint->noCols() != constraint->noCols()) { 
	opserr << "MD_Constraint::MD_Constraint - ran out of memory 2\n";
	exit(-1);
    }        
    (*constraint).Zero();

    commit_constraint = new Matrix(9,9);
    if (commit_constraint == 0 || commit_constraint->noCols() != commit_constraint->noCols()) {
        opserr << "MD_Constraint::MD_Constraint - ran out of memory 2\n";
        exit(-1);
    }
    (*commit_constraint).Zero();

}


MD_Constraint::~MD_Constraint()
{
    // invoke the destructor on the matrix and the two ID objects
    if (constraint != 0)
	  delete constraint;
    if (commit_constraint != 0)
      delete commit_constraint;
    if (constrDOF != 0)
	  delete constrDOF;
    if (retainDOF != 0)
	  delete retainDOF;  
    if(global_vector != 0) 
	  delete global_vector;  
    if(commit_global_vector != 0)
      delete commit_global_vector;  
}


int
MD_Constraint::getNodeRetained(void) const
{
    // return id of retained node
    return nodeRetained;
}

int
MD_Constraint::getNodeConstrained(void) const
{
    // return id of constrained node    
    return nodeConstrained;
}

int
MD_Constraint::commitState(void)
{
    (*commit_global_vector) = (*global_vector);
    (*commit_constraint) = (*constraint);
    cnvg = 0;
    return 0;
}


int
MD_Constraint::revertToLastCommit(void)
{
    (*global_vector) = (*commit_global_vector); 
    (*constraint) = (*commit_constraint);
    cnvg = 1;
    return 0;
}


const ID &
MD_Constraint::getConstrainedDOFs(void) const
{
    if (constrDOF == 0) {
	opserr << "MD_Constraint::getConstrainedDOF - no ID was set, ";
	opserr << "was recvSelf() ever called? or subclass incorrect?\n";	
	exit(-1);
    }

    // return the ID corresponding to constrained DOF of Ccr
    return *constrDOF;    
}

const ID &
MD_Constraint::getRetainedDOFs(void) const
{
    if (retainDOF == 0) {
	opserr << "MD_Constraint::getRetainedDOFs - no ID was set\n ";
	opserr << "was recvSelf() ever called? or subclass incorrect?\n";		
	exit(-1);
    }

    // return the ID corresponding to retained DOF of Ccr
    return *retainDOF;    
}

int 
MD_Constraint::applyConstraint(double timeStamp)
{
    /* INITIALIZE THE ROTATION MATRIX */
    int i,j;
    double r[3][3];
    for ( i = 0; i < 3; i++ ){
        for ( j = 0; j < 3; j++ ){
            r[i][j] = 0.0;
        }
    }
    /***********************************************************/
    /* RETRIEVE THE GLOBAL INCREMENTAL ROTATIONS AT THIS JOINT */
    /***********************************************************/
    const Vector &dq = ConstrainedNode->getIncrDeltaDisp();

    double alpha, beta, gamma;

    alpha = dq(3);
    beta  = dq(4);
    gamma = dq(5);

    /**************************************************************/
    /* DETERMINE THE SINE AND COSINE OF THE INCREMENTAL ROTATIONS */
    /**************************************************************/
    double sin_alpha, sin_beta, sin_gamma, cos_alpha, cos_beta, cos_gamma;
    sin_alpha = sin(alpha );
    sin_beta  = sin(beta );
    sin_gamma = sin(gamma );

    cos_alpha = cos(alpha);
    cos_beta  = cos(beta);
    cos_gamma = cos(gamma);

    /**************************************************************/
    /* COMPUTE THE TERMS OF THE ROTATION MATRIX FOR THIS EQUATION */
    /**************************************************************/
    r[0][0] = cos_gamma * cos_beta;
    r[0][1] = sin_gamma * cos_beta;
    r[0][2] = sin_beta;
    r[1][0] = - cos_gamma * sin_alpha * sin_beta - sin_gamma * cos_alpha;
    r[1][1] = - sin_gamma * sin_alpha * sin_beta + cos_gamma * cos_alpha;
    r[1][2] = sin_alpha * cos_beta;
    r[2][0] = - cos_gamma * cos_alpha * sin_beta + sin_gamma * sin_alpha;
    r[2][1] = - sin_gamma * cos_alpha * sin_beta - cos_gamma * sin_alpha;
    r[2][2] = cos_alpha * cos_beta;
    /*******************************************************/
    /* COMPUTE A NEW RESTRAINT VECTOR IN THE NEW DIRECTION */
    /*******************************************************/
    double x_length, y_length, z_length;
    x_length = 0.0;
    y_length = 0.0;
    z_length = 0.0;
    
    for ( i = 0; i < 3; i++ ){
	   x_length += (*global_vector)(i) * r[0][i];
	   y_length += (*global_vector)(i) * r[1][i];
	   z_length += (*global_vector)(i) * r[2][i];
    }
    /******************************************************************/
    /* MAKE THE VECTOR A UNIT VECTOR AND UPDATE THE LAGRANGE EQUATION */
    /******************************************************************/
    double length;
    length = sqrt( pow( x_length, 2 ) + pow( y_length, 2 ) + pow( z_length, 2 ) );

    (*global_vector)(0) = x_length / length;
    (*global_vector)(1) = y_length / length;
    (*global_vector)(2) = z_length / length;

    (*constraint).Zero();

    /* COMPUTE n*nT*/
    for( i = 0; i < 3; i++ ){
         for( j = 0; j < 3; j++ ){
               (*constraint)(i,j) =  (*global_vector)(i) * (*global_vector)(j);
               (*constraint)(i+6,j+6) = (*global_vector)(i) * (*global_vector)(j);
               (*constraint)(i+6,j) = - (*global_vector)(i) * (*global_vector)(j);
               (*constraint)(i,j+6) = - (*global_vector)(i) * (*global_vector)(j);
	 }
    }

    return 0;
}

bool
MD_Constraint::isTimeVarying(void) const
{
    return true;
}


int
MD_Constraint::update(void) 
{
    /* INITIALIZE THE ROTATION MATRIX */
    if( cnvg == 1 ){
       cnvg = 0;
       return 0;
    }

    int i,j;
    double r[3][3];
    for ( i = 0; i < 3; i++ ){
        for ( j = 0; j < 3; j++ ){
            r[i][j] = 0.0;
        }
    }
    /***********************************************************/
    /* RETRIEVE THE GLOBAL INCREMENTAL ROTATIONS AT THIS JOINT */
    /***********************************************************/
    const Vector &dq = ConstrainedNode->getIncrDeltaDisp();

    double alpha, beta, gamma;

    alpha = dq(3);
    beta  = dq(4);
    gamma = dq(5);

    alpha = 0.0;
    beta  = 0.0;
    gamma = 0.0;

    /**************************************************************/
    /* DETERMINE THE SINE AND COSINE OF THE INCREMENTAL ROTATIONS */
    /**************************************************************/
    double sin_alpha, sin_beta, sin_gamma, cos_alpha, cos_beta, cos_gamma;
    sin_alpha = sin(alpha );
    sin_beta  = sin(beta );
    sin_gamma = sin(gamma );

    cos_alpha = cos(alpha);
    cos_beta  = cos(beta);
    cos_gamma = cos(gamma);

    /**************************************************************/
    /* COMPUTE THE TERMS OF THE ROTATION MATRIX FOR THIS EQUATION */
    /**************************************************************/
    r[0][0] = cos_gamma * cos_beta;
    r[0][1] = sin_gamma * cos_beta;
    r[0][2] = sin_beta;
    r[1][0] = - cos_gamma * sin_alpha * sin_beta - sin_gamma * cos_alpha;
    r[1][1] = - sin_gamma * sin_alpha * sin_beta + cos_gamma * cos_alpha;
    r[1][2] = sin_alpha * cos_beta;
    r[2][0] = - cos_gamma * cos_alpha * sin_beta + sin_gamma * sin_alpha;
    r[2][1] = - sin_gamma * cos_alpha * sin_beta - cos_gamma * sin_alpha;
    r[2][2] = cos_alpha * cos_beta;
    /*******************************************************/
    /* COMPUTE A NEW RESTRAINT VECTOR IN THE NEW DIRECTION */
    /*******************************************************/
    double x_length, y_length, z_length;
    x_length = 0.0;
    y_length = 0.0;
    z_length = 0.0;
    
    for ( i = 0; i < 3; i++ ){
	   x_length += (*global_vector)(i) * r[0][i];
	   y_length += (*global_vector)(i) * r[1][i];
	   z_length += (*global_vector)(i) * r[2][i];
    }
    /******************************************************************/
    /* MAKE THE VECTOR A UNIT VECTOR AND UPDATE THE LAGRANGE EQUATION */
    /******************************************************************/
    double length;
    length = sqrt( pow( x_length, 2 ) + pow( y_length, 2 ) + pow( z_length, 2 ) );

    (*global_vector)(0) = x_length / length;
    (*global_vector)(1) = y_length / length;
    (*global_vector)(2) = z_length / length;

    (*constraint).Zero();

    /* COMPUTE n*nT*/
    for( i = 0; i < 3; i++ ){
         for( j = 0; j < 3; j++ ){
               (*constraint)(i,j) =  (*global_vector)(i) * (*global_vector)(j);
               (*constraint)(i+6,j+6) = (*global_vector)(i) * (*global_vector)(j);
               (*constraint)(i+6,j) = - (*global_vector)(i) * (*global_vector)(j);
               (*constraint)(i,j+6) = - (*global_vector)(i) * (*global_vector)(j);
	 }
    }

    return 0;
}


const Matrix &
MD_Constraint::getConstraint(void)
{
	if (constraint == 0) {
	opserr << "MD_Constraint::getConstraint - no Matrix was set\n";
	exit(-1);
    }

    // return the constraint matrix Ccr
    //fstream mdc;
    //mdc.open("mdc.dat",ios::app);
    //mdc>>(*constraint);
    return *constraint; 
}

int 
MD_Constraint::sendSelf(int cTag, Channel &theChannel)
{
        static ID data(10);
    int dataTag = this->getDbTag();

    data(0) = this->getTag(); 
    data(1) = nodeRetained;
    data(2) = nodeConstrained;
    if (constraint == 0) data(3) = 0; else data(3) = constraint->noRows();
    if (constraint == 0) data(4) = 0; else data(4) = constraint->noCols();    
    if (constrDOF == 0) data(5) = 0; else data(5) = constrDOF->Size();    
    if (retainDOF == 0) data(6) = 0; else data(6) = retainDOF->Size();        
    
    // need two database tags for ID objects
    if (constrDOF != 0 && dbTag1 == 0) 
      dbTag1 = theChannel.getDbTag();
    if (retainDOF != 0 && dbTag2 == 0) 
      dbTag2 = theChannel.getDbTag();

    data(7) = dbTag1;
    data(8) = dbTag2;
	data(9) = nextTag;

    int result = theChannel.sendID(dataTag, cTag, data);
    if (result < 0) {
	opserr << "WARNING MD_Constraint::sendSelf - error sending ID data\n";
	return result;  
    }    
    
    if (constraint != 0 && constraint->noRows() != 0) {
	int result = theChannel.sendMatrix(dataTag, cTag, *constraint);
	if (result < 0) {
	    opserr << "WARNING MD_Constraint::sendSelf ";
	    opserr << "- error sending Matrix data\n"; 
	    return result;  
	}
    }

    if (constrDOF != 0 && constrDOF->Size() != 0) {
	int result = theChannel.sendID(dbTag1, cTag, *constrDOF);
	if (result < 0) {
	    opserr << "WARNING MD_Constraint::sendSelf ";
	    opserr << "- error sending constrained data\n"; 
	    return result;  
	}
    }

    if (retainDOF != 0 && retainDOF->Size() != 0) {
	int result = theChannel.sendID(dbTag2, cTag, *retainDOF);
	if (result < 0) {
	    opserr << "WARNING MD_Constraint::sendSelf ";
	    opserr << "- error sending retained data\n"; 
	    return result;  
	}
    }
    
    return 0;
}


int 
MD_Constraint::recvSelf(int cTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    static ID data(10);
    int result = theChannel.recvID(dataTag, cTag, data);
    if (result < 0) {
	opserr << "WARNING MD_Constraint::recvSelf - error receiving ID data\n";
	return result;  
    }    

    this->setTag(data(0));
    nodeRetained = data(1);
    nodeConstrained = data(2);
    int numRows = data(3); 
    int numCols = data(4);
    dbTag1 = data(7);
    dbTag2 = data(8);
    nextTag = data(9);

    if (numRows != 0 && numCols != 0) {
	constraint = new Matrix(numRows,numCols);
	
	int result = theChannel.recvMatrix(dataTag, cTag, *constraint);
	if (result < 0) {
	    opserr << "WARNING MD_Constraint::recvSelf ";
	    opserr << "- error receiving Matrix data\n"; 
	    return result;  
	}
    }    
    int size = data(5);
    if (size != 0) {
	constrDOF = new ID(size);
	int result = theChannel.recvID(dbTag1, cTag, *constrDOF);
	if (result < 0) {
	    opserr << "WARNING MD_Constraint::recvSelf ";
	    opserr << "- error receiving constrained data\n"; 
	    return result;  
	}	
    }
    
    size = data(6);
    if (size != 0) {
	retainDOF = new ID(size);
	int result = theChannel.recvID(dbTag2, cTag, *retainDOF);
	if (result < 0) {
	    opserr << "WARNING MP_Retainaint::recvSelf ";
	    opserr << "- error receiving retained data\n"; 
	    return result;  
	}	
    }    
    
    return 0;
}



void
MD_Constraint::Print(OPS_Stream &s, int flag)
{     
    s << "MD_Constraint: " << this->getTag() << "\n";
    s << "\tNode Constrained: " << nodeConstrained;
    s << " node Retained: " << nodeRetained ;
    if (constrDOF != 0)
	s << " constrained dof: " << *constrDOF;    
    if (retainDOF != 0)
	s << " retained dof: " << *retainDOF;        
    if (constraint != 0)
	s << " constraint matrix: " << *constraint << "\n";
}


