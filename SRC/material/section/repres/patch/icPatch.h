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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/repres/patch/icPatch.h,v $
                                                                        
                                                                        
// File: ic_Patch.h
// Written by Cenk Tort
// September 2004

#ifndef icPatch_h 
#define icPatch_h 

#include <Patch.h>

class Cell;
class Matrix;

class icPatch: public Patch
{
  public:

    icPatch();
    icPatch(int materialID, double D, double B, double TF, int SFL_D, int SFL_B);
        
    ~icPatch();
    
    // edition functions

    void setMaterialID(int materialID);
    void setDiscretization(double D, double B, double TF, int SFL_D, int SFL_B);
   
    // reinforcing bar inquiring functions
    
    int     getMaterialID         (void) const; 
    int     getNumCells           (void) const;
    Cell  **getCells              (void) const;
    Patch  *getCopy               (void) const;

    void getDiscretization(double &D, double &B, double &TF, int &SFL_D, int &SFL_B) const;

    void Print(OPS_Stream &s, int flag =0) const;   
    friend OPS_Stream &operator<<(OPS_Stream &s, icPatch &ic_Patch);    
    
  protected:
    
  private:
    double d;
    double b;
    double tf;
    int    matID;
    int    sfl_d;
    int    sfl_b;
};


#endif
 
