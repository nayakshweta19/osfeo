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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2007-09-21 15:28:12 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/domain/constraints/MD_Constraint.h,v $
                                                                        
                                                                        
#ifndef MD_Constraint_h
#define MD_Constraint_h

#include <DomainComponent.h>
#include <bool.h>
#include <MP_Constraint.h>
#include <Node.h>
#include <Domain.h>

class Matrix;
class ID;

class MD_Constraint : public MP_Constraint
{
  public:
    // constructors        
    MD_Constraint();		//Cenk 

    MD_Constraint(Domain *theDomain,
		  int tag,
		  int nodeRetain, 
		  int nodeConstr, 
		  ID &constrainedDOF,
    		  ID &retainedDOF,
		  Vector &gv);

    // destructor    
    ~MD_Constraint();

    // method to get information about the constraint
    int getNodeRetained(void) const;
    int getNodeConstrained(void) const;    
    const ID &getConstrainedDOFs(void) const;        
    const ID &getRetainedDOFs(void) const;            
    int applyConstraint(double pseudoTime);
    bool isTimeVarying(void) const;
    const Matrix &getConstraint(void);    
    int commitState(void);
    int revertToLastCommit(void);
    int update(void);
    // methods for output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag =0);

  protected:
    
  private:
    int cnvg;
    int nodeRetained;        // to identify the retained node
    int nodeConstrained;     // to identify  the constrained node
    Matrix *constraint;      // pointer to the constraint matrix
    Matrix *commit_constraint;
    ID *constrDOF;           // ID of constrained DOF at constrained node
    ID *retainDOF;           // ID of related DOF at retained node
    Node *ConstrainedNode;
    Vector *global_vector;
    Vector *commit_global_vector;
    Domain *thisDomain;
    
    int dbTag1, dbTag2, dbTag3;      // need a dbTag for the two ID's
};

#endif

