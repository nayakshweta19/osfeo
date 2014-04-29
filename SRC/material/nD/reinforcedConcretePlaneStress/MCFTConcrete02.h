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
// $Date: 2007/06/08 00:38:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nd/reinforcedConcretePlaneStress/MCFTConcrete02.h,v $
                                                                      
// Written: Li Ning (Neallee@tju.edu.cn)
// Created: 07/11
//
// Description: This file contains the class definition for 
// MCFTConcrete02. MCFTConcrete02 is based on an f2c of the FEDEAS material

#ifndef MCFTConcrete02_h
#define MCFTConcrete02_h

#include <UniaxialMaterial.h>

class MCFTConcrete02 : public UniaxialMaterial
{
  public:
    MCFTConcrete02(int tag, double _fc, double _epsc0, double _fcu,
	     double _epscu, double _rat, double _ft, double _Ets);

    MCFTConcrete02(void);

    virtual ~MCFTConcrete02();

    const char *getClassType(void) const {return "MCFTConcrete02";};    
    double getInitialTangent(void);
    UniaxialMaterial *getCopy(void);

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);      
    double getStress(void);
    double getTangent(void);
	double getSecant(void);
    double getPD(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    

	Response *setResponse (const char **argv, int argc, 
		OPS_Stream &theOutputStream);
	int getResponse (int responseID, Information &matInformation);   
    
	void Print(OPS_Stream &s, int flag =0);

    int getVariable(const char *variable, Information &);
    
 protected:
    
 private:
    void Tens_Envlp (double epsc, double &sigc, double &Ect);
    void Compr_Envlp (double epsc, double &sigc, double &Ect);

	// Calculate the plastic strain offsets
	void determineCompEpscp(double eps);
	void determineTensEpscp(double eps);

    // matpar : Concrete FIXED PROPERTIES
    double fc;    // concrete compression strength           : mp(1)
    double epsc0; // strain at compression strength          : mp(2)
    double fcu;   // stress at ultimate (crushing) strain    : mp(3)
    double epscu; // ultimate (crushing) strain              : mp(4)       
    double rat;   // ratio between unloading slope at epscu and original slope : mp(5)
    double ft;    // concrete tensile strength               : mp(6)
    double Ets;   // tension stiffening slope                : mp(7)

    // hstvP : Concrete HISTORY VARIABLES last committed step
    double ecminP;  //  hstP(1)
    double deptP;   //  hstP(2)
    double epsP;  //  = strain at previous converged step
    double sigP;  //  = stress at previous converged step
    double eP;    //   stiffness modulus at last converged step;

    // hstv : Concrete HISTORY VARIABLES  current step
    double ecmin;  
    double dept;   
    double sig;   
    double e;     
    double eps;   

	double epscp; // plastic strain offset
};


#endif

