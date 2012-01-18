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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/2dLayer/Circ2dLayer.cpp,v $
                                                                        
                                                                        
// File: Circ2dLayer.C 
// Written by Remo M. de Souza 
// December 1998

#include <math.h>
#include <Matrix.h>
#include <Vector.h>

#include <ReinfBar.h>
#include <Circ2dLayer.h>


Circ2dLayer::Circ2dLayer(void):
                        nReinfBars(0), matID(0), barDiam(0.0),
                        area(0.0), centerPosit(2), arcRad(0.0),
                        initAng(0.0), finalAng(0.0)          
{
}


Circ2dLayer::Circ2dLayer(int materialID, int numReinfBars, 
                               double reinfBarArea,
                               const Vector &centerPosition,
                               double arcRadius, double initialAngle,
                               double finalAngle):
                               nReinfBars(numReinfBars),
                               matID(materialID), area(reinfBarArea),
                               barDiam(0.0),centerPosit(centerPosition),
                               arcRad(arcRadius),initAng(initialAngle), 
                               finalAng(finalAngle)
{
}

Circ2dLayer::Circ2dLayer(int materialID, int numReinfBars, double  reinfBarArea,
							   const Vector &centerPosition, double radius):
nReinfBars(numReinfBars), matID(materialID), area(reinfBarArea),
barDiam(0.0), centerPosit(centerPosition), arcRad(radius),
initAng(0.0), finalAng(0.0)
{
	// Figure out final angle so that complete circle does not put
	// two bars at the same location
	if (nReinfBars > 0)
		finalAng = 360.0 - 360.0/nReinfBars;
}

Circ2dLayer::~Circ2dLayer()
{

}


void Circ2dLayer::setNumReinfBars(int numReinfBars)
{
   nReinfBars = numReinfBars;
}

void Circ2dLayer::setMaterialID (int materialID)
{
   matID = materialID;
}

void Circ2dLayer::setReinfBarDiameter (double reinfBarDiameter)
{
   barDiam = reinfBarDiameter;
   double pi = acos(-1.0);
   area = pi * barDiam*barDiam/4.0;
}

void Circ2dLayer::setReinfBarArea(double reinfBarArea)
{
   area = reinfBarArea;
}


int Circ2dLayer::getNumReinfBars (void) const
{
   return nReinfBars;
}

int Circ2dLayer::getMaterialID (void) const
{
   return matID;
}

double Circ2dLayer::getReinfBarDiameter (void) const
{
   return barDiam;
}

double Circ2dLayer::getReinfBarArea (void) const
{
   return area;
}

ReinfBar * 
Circ2dLayer::getReinfBars (void) const
{
   double theta, dtheta;
   Vector barPosit(2);
   int i;
   ReinfBar *reinfBars;
   double pi = acos(-1.0);
   double initAngRad, finalAngRad;

   if (nReinfBars > 1)
   {
      initAngRad  = pi * initAng  / 180.0;
      finalAngRad = pi * finalAng / 180.0;
 
      dtheta = (finalAngRad - initAngRad) /(nReinfBars - 1);

      reinfBars = new ReinfBar [nReinfBars];

      for (i = 0; i < nReinfBars; i++)
      {
         theta = initAngRad + dtheta * i;
         barPosit(0) = centerPosit(0) + arcRad*cos(theta);
         barPosit(1) = centerPosit(1) + arcRad*sin(theta);

         reinfBars[i].setPosition(barPosit);
         reinfBars[i].setArea(this->area);
      }
   }
   else
      return 0;

   return reinfBars;         
}


const Vector & 
Circ2dLayer::getCenterPosition(void) const
{
   return centerPosit;
}

double Circ2dLayer::getArcRadius(void) const 
{
   return arcRad;
}

double Circ2dLayer::getInitAngle(void) const 
{
   return initAng;
}

double Circ2dLayer::getFinalAngle(void) const 
{
   return finalAng;
}


ReinfLayer * 
Circ2dLayer::getCopy (void) const
{
   Circ2dLayer *theCopy = new Circ2dLayer (matID, nReinfBars, area,
                                                 centerPosit, arcRad,
                                                 initAng, finalAng);
   return theCopy;
}



void Circ2dLayer::Print(OPS_Stream &s, int flag) const
{
   s << "\nReinforcing Layer type:  Circ";
   s << "\nMaterial ID: " << matID;
   s << "\nReinf. bar diameter: " << barDiam;
   s << "\nReinf. bar area: " << area;
   s << "\nCenter Position: " << centerPosit;
   s << "\nArc Radius: " << arcRad;
   s << "\nInitial angle: " << initAng;
   s << "\nFinal angle: " << finalAng;
}


OPS_Stream &operator<<(OPS_Stream &s, const Circ2dLayer &Circ2dLayer)
{  
   Circ2dLayer.Print(s);
   return s;
}
 
