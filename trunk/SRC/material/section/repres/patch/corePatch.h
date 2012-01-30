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
**   Cenk Tort (tort0008@umn.edu)				      **
**   Jerome F. Hajjar (hajjar@struc.ce.umn.edu)                       **
**								      **
**   University of Minnesota - Twin Cities			      **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2007-09-21 15:28:40 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/repres/patch/corePatch.h,v $
                                                                        
                                                                        
// File: corePatch.h
// Written by Cenk Tort
// September 2004

#ifndef corePatch_h 
#define corePatch_h 

#include <Patch.h>

class Cell;
class Matrix;

class corePatch: public Patch
{
  public:

    corePatch();
    corePatch(int materialID, double D, double B, double T, int CFL_D, int CFL_B, int CWL_D, int CWL_B, int CC_D, int CC_B);
        
    ~corePatch();
    
    // edition functions

    void setMaterialID(int materialID);
    void setDiscretization(double D, double B, double T, int CFL_D, int CFL_B, int CWL_D, int CWL_B, int CC_D, int CC_B);
   
    // reinforcing bar inquiring functions
    
    int     getMaterialID         (void) const; 
    int     getNumCells           (void) const;
    Cell  **getCells              (void) const;
    Patch  *getCopy               (void) const;

    void   getDiscretization(double &D, double &B, double &T, int &CFL_D, int &CFL_B, int &CWL_D, int &CWL_B, int &CC_D,  int &CC_B) const;

    void Print(OPS_Stream &s, int flag =0) const;   
    friend OPS_Stream &operator<<(OPS_Stream &s, corePatch &core_Patch);    
    
  protected:
    
  private:
    int    matID;
    double d;
    double b;
    double t;
    int    cfl_d;
    int    cfl_b;
    int    cwl_d;
    int    cwl_b;
    int    cc_d;
    int    cc_b;
};


#endif
 
