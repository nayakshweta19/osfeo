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
**   Andreas Schellenberg (andreas.schellenberg@gmail.com)            **
**   Yoshikazu Takahashi (yos@catfish.dpri.kyoto-u.ac.jp)             **
**   Gregory L. Fenves (fenves@berkeley.edu)                          **
**   Stephen A. Mahin (mahin@berkeley.edu)                            **
**                                                                    **
** ****************************************************************** */

// $Revision: 358 $
// $Date: 2014-01-31 10:16:15 +0800 (星期五, 31 一月 2014) $
// $URL: svn://opensees.berkeley.edu/usr/local/svn/OpenFresco/trunk/SRC/experimentalControl/ECSCRAMNetGT.h $

#ifndef ECSCRAMNetGT_h
#define ECSCRAMNetGT_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 08/13
// Revision: A
//
// Description: This file contains the class definition for ECSCRAMNetGT.
// ECSCRAMNetGT is a controller class for communicating with a shared
// common RAM network (SCRAMNet GT).

#include "ExperimentalControl.h"

extern "C" {
#include <scgtapi.h>
}

class ECSCRAMNetGT : public ExperimentalControl
{
public:
    // constructors
    ECSCRAMNetGT(int tag, int memOffset, int numActCh);
    ECSCRAMNetGT(const ECSCRAMNetGT &ec);
    
    // destructor
    virtual ~ECSCRAMNetGT();
    
    // method to get class type
    const char *getClassType() const {return "ECSCRAMNetGT";};
    
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
    
    scgtHandle gtHandle;
    scgtInterrupt interrupt;

    scgtDeviceInfo *deviceInfo;
    
    const int memOffset;
    const int numActCh;

    const int *memPtrBASE;
	float *memPtrOPF;

    unsigned int *newTarget, *switchPC, *atTarget;
    float *ctrlDisp, *ctrlVel, *ctrlAccel, *ctrlForce, *ctrlTime;
    float *daqDisp, *daqVel, *daqAccel, *daqForce, *daqTime;

    unsigned int flag;
};

#endif