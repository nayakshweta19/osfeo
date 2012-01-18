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

// $Revision: 1.3 $
// $Date: 2009/05/19 22:10:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/NewmarkExplicit.h,v $


#ifndef NewmarkExplicit_h
#define NewmarkExplicit_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net) 
// Created: 02/05
// Revision: A
//
// Description: This file contains the class definition for NewmarkExplicit.
// NewmarkExplicit is an algorithmic class for performing a transient analysis
// using the explicit Newmark integration scheme (beta = 0).
//
// What: "@(#) NewmarkExplicit.h, revA"

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class NewmarkExplicit : public TransientIntegrator
{
public:
    // constructors
    NewmarkExplicit();
    NewmarkExplicit(double gamma,
        bool updDomFlag = false);
    NewmarkExplicit(double gamma,
        double alphaM, double betaK, double betaKi, double betaKc,
        bool updDomFlag = false); 
    
    // destructor
    ~NewmarkExplicit();
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);        
    
    int domainChanged(void);    
    int newStep(double deltaT);    
    int revertToLastStep(void);            
    int update(const Vector &aiPlusOne);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);        
    
protected:
    
private:
    double gamma;
    bool updDomFlag;    // a flag indicating if updateDomain() is called

    // rayleigh damping factors
    double alphaM;
    double betaK;
    double betaKi;
    double betaKc;

    int updateCount;                // method should only have one update per step
    double c2, c3;                  // some constants we need to keep
    Vector *Ut, *Utdot, *Utdotdot;  // response quantities at time t
    Vector *U, *Udot, *Udotdot;     // response quantities at time t+deltaT
};

#endif
