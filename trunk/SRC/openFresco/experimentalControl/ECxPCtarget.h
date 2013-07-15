/* ****************************************************************** **
**    OpenFRESCO - Open Framework                                     **
**                 for Experimental Setup and Control                 **
**                                                                    **
**                                                                    **
** Copyright (c) 2006, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited. See    **
** file 'COPYRIGHT_UCB' in main directory for information on usage    **
** and redistribution, and for a DISCLAIMER OF ALL WARRANTIES.        **
**                                                                    **
** Developed by:                                                      **
**   Andreas Schellenberg (andreas.schellenberg@gmx.net)              **
**   Yoshikazu Takahashi (yos@catfish.dpri.kyoto-u.ac.jp)             **
**   Gregory L. Fenves (fenves@berkeley.edu)                          **
**   Stephen A. Mahin (mahin@berkeley.edu)                            **
**                                                                    **
** ****************************************************************** */

// $Revision: 314 $
// $Date: 2011-05-23 05:17:07 +0800 (星期一, 23 五月 2011) $
// $URL: svn://opensees.berkeley.edu/usr/local/svn/OpenFresco/trunk/SRC/experimentalControl/ECxPCtarget.h $

#ifndef ECxPCtarget_h
#define ECxPCtarget_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 09/06
// Revision: A
//
// Description: This file contains the class definition for ECxPCtarget.
// ECxPCtarget is a controller class for communicating with a xPC Target
// digital signal processor.

#include "ExperimentalControl.h"

class ECxPCtarget : public ExperimentalControl
{
public:
    // constructors
    ECxPCtarget(int tag, int pcType, char *ipAddress,
        char *ipPort, char *appName, char *appPath = 0,
        int timeOut = 10);
    ECxPCtarget(const ECxPCtarget &ec);
    
    // destructor
    virtual ~ECxPCtarget();
    
    // method to get class type
    const char *getClassType() const {return "ECxPCtarget";};
    
    // public methods to set and to get response
    virtual int setup();
    virtual int setSize(ID sizeT, ID sizeO);
    
    virtual int setTrialResponse(const Vector* disp, 
        const Vector* vel,
        const Vector* accel,
        const Vector* force,
        const Vector* time);
    virtual int getDaqResponse(Vector* disp,
        Vector* vel,
        Vector* accel,
        Vector* force,
        Vector* time);
    
    virtual int commitState();
    
    virtual ExperimentalControl *getCopy();
    
    // public methods for experimental control recorder
    virtual Response *setResponse(const char **argv, int argc,
        OPS_Stream &output);
    virtual int getResponse(int responseID, Information &info);
    
    // public methods for output
    void Print(OPS_Stream &s, int flag = 0);    
    
protected:
    // protected methods to set and to get response
    virtual int control();
    virtual int acquire();
    
private:
    void sleep(const clock_t wait);
    
    int pcType, port, timeOut;
    char *ipAddress, *ipPort, *appName, *appPath;
    char errMsg[256];

    double newTarget, switchPC, atTarget;
    double *ctrlDisp, *ctrlVel, *ctrlAccel;
    double *daqDisp, *daqForce;
    
    int newTargetId, switchPCId, atTargetId;
    int ctrlDispId, ctrlVelId, ctrlAccelId;
    int *daqDispId, *daqForceId;
};

#endif