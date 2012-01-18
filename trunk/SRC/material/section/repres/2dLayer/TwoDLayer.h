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
// $Date: 2003/02/14 23:01:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/2dLayer/TwoDLayer.h,v $
                                                                        
                                                                        
// File: TwoDLayer.h
// Written by Remo M. de Souza
// December 1998

#ifndef TwoDLayer_h 
#define TwoDLayer_h 

#include <OPS_Globals.h>

class ReinfBar;

class TwoDLayer
{
  public:

    TwoDLayer();
    virtual ~TwoDLayer();
    
    // edition functions

    virtual void setNumReinfBars     (int numReinfBars)        = 0;
    virtual void setMaterialID       (int materialID)          = 0;
    virtual void setReinfBarDiameter (double reinfBarDiemater) = 0;
    virtual void setReinfBarArea     (double reinfBarArea)     = 0;

    // reinforcing layer inquiring functions
    
    virtual int         getNumReinfBars     (void) const = 0;
    virtual int         getMaterialID       (void) const = 0; 
    virtual double      getReinfBarDiameter (void) const = 0;
    virtual double      getReinfBarArea     (void) const = 0;
    virtual TwoDLayer *getCopy             (void) const = 0;
    virtual ReinfBar   *getReinfBars        (void) const = 0;     
   
    virtual void Print(OPS_Stream &s, int flag =0) const = 0;   
    friend OPS_Stream &operator<<(OPS_Stream &s, const TwoDLayer &TwoDLayer);    
    
  protected:
    
  private:
};


#endif

