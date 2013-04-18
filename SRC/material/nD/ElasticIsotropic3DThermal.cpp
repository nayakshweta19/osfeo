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
                                                                        
#include <ElasticIsotropic3DThermal.h>           
#include <Channel.h>

Vector ElasticIsotropic3DThermal::sigma(6);
Matrix ElasticIsotropic3DThermal::D(6,6);

ElasticIsotropic3DThermal::ElasticIsotropic3DThermal
(int tag, double E, double nu, double rho) :
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropic3DThermal, E, nu, rho),
 epsilon(6), Cepsilon(6),E0T(E),ThermalElongation(0)
{
  E=E0T;
  epsilon.Zero();
  Cepsilon.Zero();
}

ElasticIsotropic3DThermal::ElasticIsotropic3DThermal():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropic3DThermal, 0.0, 0.0),
 epsilon(6), Cepsilon(6),E0T(E),ThermalElongation(0)
{
  E=E0T;
  epsilon.Zero();
  Cepsilon.Zero();
}

ElasticIsotropic3DThermal::~ElasticIsotropic3DThermal ()
{

}

int
ElasticIsotropic3DThermal::setTrialStrain (const Vector &strain)
{
  epsilon = strain;
  return 0;
}

int
ElasticIsotropic3DThermal::setTrialStrain (const Vector &strain, const Vector &rate)
{
  epsilon = strain;
  return 0;
}

int
ElasticIsotropic3DThermal::setTrialStrainIncr (const Vector &strain)
{
  epsilon += strain;
  return 0;
}

int
ElasticIsotropic3DThermal::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  epsilon += strain;
  return 0;
}

const Matrix&
ElasticIsotropic3DThermal::getTangent (void)
{
  double mu2 = E/(1.0+v);
  double lam = v*mu2/(1.0-2.0*v);
  double mu = 0.50*mu2;
  mu2 += lam;

  D(0,0) = D(1,1) = D(2,2) = mu2;
  D(0,1) = D(1,0) = D(0,2) = D(2,0) = D(1,2) = D(2,1) = lam;
  D(3,3) = mu;
  D(4,4) = mu;
  D(5,5) = mu;

  return D;
}

const Matrix&
ElasticIsotropic3DThermal::getInitialTangent (void)
{
  //  return this->getTangent();
  double mu2 = E/(1.0+v);
  double lam = v*mu2/(1.0-2.0*v);
  double mu = 0.50*mu2;
  mu2 += lam;

  D(0,0) = D(1,1) = D(2,2) = mu2;
  D(0,1) = D(1,0) = D(0,2) = D(2,0) = D(1,2) = D(2,1) = lam;
  D(3,3) = mu;
  D(4,4) = mu;
  D(5,5) = mu;

  return D;
}

const Vector&
ElasticIsotropic3DThermal::getStress (void)
{
  double mu2 = E/(1.0+v);
  double lam = v*mu2/(1.0-2.0*v);
  double mu = 0.50*mu2;
  mu2 += lam;

  double eps0 = epsilon(0);
  double eps1 = epsilon(1);
  double eps2 = epsilon(2);

  D(0,0) = D(1,1) = D(2,2) = mu2;
  D(0,1) = D(1,0) = D(0,2) = D(2,0) = D(1,2) = D(2,1) = lam;
  D(3,3) = mu;
  D(4,4) = mu;
  D(5,5) = mu;

  //sigma = D*epsilon;
  sigma(0) = mu2*eps0 + lam*(eps1+eps2);
  sigma(1) = mu2*eps1 + lam*(eps0+eps2);
  sigma(2) = mu2*eps2 + lam*(eps0+eps1);
  sigma(3) = mu*epsilon(3);
  sigma(4) = mu*epsilon(4);
  sigma(5) = mu*epsilon(5);
	
  return sigma;
}
double 
ElasticIsotropic3DThermal::setThermalTangentAndElongation(double &TempT, double&ET, double&Elong)
{

	// EN 1992 pt 1-2-1. Class N hot rolled  reinforcing steel at elevated temperatures
  if (TempT <= 100) {
		E = E0T*(1 - (TempT - 20)*0.375/80);
  }
  else if (TempT <= 200) {
      E = E0T*(0.625 - (TempT - 100)*0.193/100);
  }
  else if (TempT <= 300) {
      E = E0T*(0.432 - (TempT - 200)*0.128/100);
  }
  else if (TempT <= 400) {
		E = E0T*(0.304 - (TempT - 300)*0.1165/100);
  }
  else if (TempT <= 500) {
	  E = E0T*(0.1875 - (TempT - 400)*0.0875/100);
  }
  else if (TempT <= 600) {
	  E = E0T*(0.1 - (TempT - 500)*0.055/100);
  }
  else if (TempT <= 700) {
	  E = E0T*(0.045 - (TempT - 600)*0.015/100);
  }
  else if (TempT <= 800) {
	  E = E0T*(0.03 - (TempT - 700)*0.015/100);
  }
  else if (TempT <= 900) {
	  E = E0T*(0.015 - (TempT - 800)*0.007/100);
  }
  else if (TempT <= 1000) {
	  E = E0T*(0.008 - (TempT - 900)*0.004/100);
  }
  else if (TempT <= 1100) {
	  E = E0T*(0.004 - (TempT - 1000)*0.003/100);
  }

  else  {
      opserr << "the temperature is invalid\n"; 
  } 

  // Calculate thermal elongation 
  // concrete
     if (TempT <= 20) {
      ThermalElongation = 0.0;
  }
  else if (TempT <= 700) {
      ThermalElongation = -1.8e-4 + 9e-6 *TempT + 2.3e-11 *TempT*TempT*TempT;
  }
  else if (TempT <= 1200) {
      ThermalElongation = 14e-3;
  }
  else {
	  opserr << "the temperature is invalid\n";
  }

  //ET = E;  
  ET =E;
  Elong = ThermalElongation;

  return 0;
}

const Vector&
ElasticIsotropic3DThermal::getStrain (void)
{
  return epsilon;
}

int
ElasticIsotropic3DThermal::commitState (void)
{
  Cepsilon=epsilon;
  return 0;
}

int
ElasticIsotropic3DThermal::revertToLastCommit (void)
{
  epsilon=Cepsilon;
  return 0;
}

int
ElasticIsotropic3DThermal::revertToStart (void)
{
  epsilon.Zero();
  Cepsilon.Zero();
  return 0;
}

NDMaterial*
ElasticIsotropic3DThermal::getCopy (void)
{
  ElasticIsotropic3DThermal *theCopy =
    new ElasticIsotropic3DThermal (this->getTag(), E, v, rho);
  
  theCopy->epsilon = epsilon;
  theCopy->Cepsilon = Cepsilon;
  
  return theCopy;
}

const char*
ElasticIsotropic3DThermal::getType (void) const
{
  return "ThreeDimensionalimensional";
}

int
ElasticIsotropic3DThermal::getOrder (void) const
{
  return 6;
}

int 
ElasticIsotropic3DThermal::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(10);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = v;
  data(3) = rho;
  data(4) = Cepsilon(0);
  data(5) = Cepsilon(1);
  data(6) = Cepsilon(2);
  data(7) = Cepsilon(3);
  data(8) = Cepsilon(4);
  data(9) = Cepsilon(5);
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ElasticIsotropic3DThermal::sendSelf -- could not send Vector\n";
    return res;
  }
  
  return res;
}

int 
ElasticIsotropic3DThermal::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(10);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ElasticIsotropic3DThermal::sendSelf -- could not send Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  E = data(1);
  v = data(2);
  rho = data(3);
  epsilon(0)=data(4);
  epsilon(1)=data(5);
  epsilon(2)=data(6);
  epsilon(3)=data(7);
  epsilon(4)=data(8);
  epsilon(5)=data(9);

  Cepsilon = epsilon;

  return res;
}


