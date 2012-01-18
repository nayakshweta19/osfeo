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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/2dLayer/Straight2dLayer.h,v $
                                                                        
                                                                        
// File: Straight2dLayer.h
// Written by Remo M. de Souza
// December 1998


#ifndef Straight2dLayer_h 
#define Straight2dLayer_h 

#include <ReinfLayer.h>

class ReinfBar;

class Straight2dLayer : public ReinfLayer
{
  public:

    Straight2dLayer();
    Straight2dLayer(int materialID, int numReinfBars, double  reinfBarArea,
                       const Vector &initialPosition, 
                       const Vector &finalPosition);

    ~Straight2dLayer();
    
    // edition functions

    void setNumReinfBars     (int numReinfBars);
    void setMaterialID       (int materialID);  
    void setReinfBarDiameter (double reinfBarDiameter);
    void setReinfBarArea     (double reinfBarArea);

    void setInitialPosition (const Vector &initialPosition);
    void setFinalPosition   (const Vector &finalPosition);

    // inquiring functions

    int           getNumReinfBars     (void) const;
    int           getMaterialID       (void) const;
    double        getReinfBarDiameter (void) const;
    double        getReinfBarArea     (void) const;
    ReinfBar     *getReinfBars        (void) const;

  
    ReinfLayer   *getCopy             (void) const;
    const Vector &getInitialPosition  (void) const;
    const Vector &getFinalPosition    (void) const;

    void Print(OPS_Stream &s, int flag =0) const;   
    friend OPS_Stream &operator<<(OPS_Stream &s, const Straight2dLayer &Straight2dLayer);    
    
  protected:
    
  private:
    int    nReinfBars;
    int    matID;
    double barDiam;
    double area;
    Vector initPosit;
    Vector finalPosit;
};


#endif

