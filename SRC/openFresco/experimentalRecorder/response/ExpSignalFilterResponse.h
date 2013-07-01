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

// $Revision: 329 $
// $Date: 2012-01-21 13:23:43 +0800 (星期六, 21 一月 2012) $
// $URL: svn://opensees.berkeley.edu/usr/local/svn/OpenFresco/trunk/SRC/experimentalRecorder/response/ExpSignalFilterResponse.h $

#ifndef ExpSignalFilterResponse_h
#define ExpSignalFilterResponse_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 08/08
// Revision: A
//
// Description: This file contains the class definition for
// ExpSignalFilterResponse.

#include <Response.h>

class ExperimentalSignalFilter;
class ID;
class Vector;
class Matrix;
class Tensor;

class ExpSignalFilterResponse : public Response
{
public:
    ExpSignalFilterResponse(ExperimentalSignalFilter *esf, int id);
    ExpSignalFilterResponse(ExperimentalSignalFilter *esf, int id, int val);
    ExpSignalFilterResponse(ExperimentalSignalFilter *esf, int id, double val);
    ExpSignalFilterResponse(ExperimentalSignalFilter *esf, int id, const ID &val);
    ExpSignalFilterResponse(ExperimentalSignalFilter *esf, int id, const Vector &val);
    ExpSignalFilterResponse(ExperimentalSignalFilter *esf, int id, const Matrix &val);
    ~ExpSignalFilterResponse();
    
    int getResponse();
    
private:
    ExperimentalSignalFilter *theFilter;
    int responseID;
};

#endif
