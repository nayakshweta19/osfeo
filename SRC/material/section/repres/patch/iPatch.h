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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/repres/patch/iPatch.h,v $
                                                                        
                                                                        
// File: i_Patch.h
// Written by Cenk Tort
// September 2004

#ifndef iPatch_h 
#define iPatch_h 

#include <Patch.h>

class Cell;
class Matrix;

class iPatch: public Patch
{
  public:

    iPatch();
    iPatch(int materialID, double D, double B, double TF, double TW, double FILLET, int SFL_D, int SFL_B, 
	   int SWL_D, int SWL_B);
        
    ~iPatch();
    
    // edition functions

    void setMaterialID(int materialID);
    void setDiscretization(double D, double B, double TF, double TW, double FILLET, int SFL_D, int SFL_B, 
           int SWL_D, int SWL_B);
   
    // reinforcing bar inquiring functions
    
    int     getMaterialID         (void) const; 
    int     getNumCells           (void) const;
    Cell  **getCells              (void) const;
    Patch  *getCopy               (void) const;

    void getDiscretization(double &D, double &B, double &TF, double &TW, double &FILLET, int &SFL_D, int &SFL_B, 
	   int &SWL_D, int &SWL_B) const;

    void Print(OPS_Stream &s, int flag =0) const;   
    friend OPS_Stream &operator<<(OPS_Stream &s, iPatch &i_Patch);    
    
  protected:
    
  private:
    double d;
    double b;
    double tf;
    double tw;
    double fillet;
    int    matID;
    int    sfl_d;
    int    sfl_b;
    int    swl_d;
    int    swl_b;
};


#endif
 
