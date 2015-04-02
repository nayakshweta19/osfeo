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

// $Revision: 387 $
// $Date: 2015-03-21 13:07:49 +0800 (星期六, 21 三月 2015) $
// $URL: svn://opensees.berkeley.edu/usr/local/svn/OpenFresco/trunk/SRC/experimentalControl/ECSimulation.cpp $

// Written: Yoshi (yos@catfish.dpri.kyoto-u.ac.jp)
// Created: 09/06
// Revision: A
//
// Description: This file contains the implementation of
// ECSimulation.

#include "ECSimulation.h"


ECSimulation::ECSimulation(int tag)
    : ExperimentalControl(tag)
{
    // does nothing
}


ECSimulation::ECSimulation(const ECSimulation& ec)
    : ExperimentalControl(ec)
{
    // does nothing
}


ECSimulation::~ECSimulation()
{
    // does nothing
}
