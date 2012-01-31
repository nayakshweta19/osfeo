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
// $Date: 2007-09-21 15:28:07 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/analysis/fe_ele/penalty/PenaltyMD_FE.h,v $
                                                                        
                                                                        
#ifndef PenaltyMD_FE_h
#define PenaltyMD_FE_h

#include <FE_Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>

class Element;
class Integrator;
class IncrementalIntegrator;
class AnalysisModel;
class Domain;
class MP_Constraint;
class Node;

class PenaltyMD_FE: public FE_Element
{
  public:
    PenaltyMD_FE(int tag, Domain &theDomain, MP_Constraint &theMP, double alpha);
    virtual ~PenaltyMD_FE();    

    // public methods
    virtual int  setID(void);
    virtual const Matrix &getTangent(Integrator *theIntegrator);
    virtual const Vector &getResidual(Integrator *theIntegrator);
    virtual const Vector &getTangForce(const Vector &x, double fact = 1.0);
 
    virtual const Vector &getK_Force(const Vector &x, double fact = 1.0);
    virtual const Vector &getC_Force(const Vector &x, double fact = 1.0);
    virtual const Vector &getM_Force(const Vector &x, double fact = 1.0);
   
  protected:
    
  private:
    void determineTangent(void);
    
    MP_Constraint *theMP;
    Node *theConstrainedNode;
    Node *theRetainedNode;    

    Matrix *tang;
    Vector *resid;
    double alpha;
	
    
};

#endif


