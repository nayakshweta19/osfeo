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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2007-09-21 15:28:41 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/uniaxial/RCFT_concMaterial.h,v $

// Written: MHS
// Created: Aug 2001

// Modified by : Cenk Tort - tort0008@umn.edu
// University of Minnesota - Civil Engineering Department
// Date : Wed Jul 23 17:32:23 EDT 2003
// Description: This file contains the class definition for
// cyclic uniaxial concrete stress-strain relationship for
// RCFT members.


#ifndef RCFT_concMaterial_h
#define RCFT_concMaterial_h

#include <UniaxialMaterial.h>

class RCFT_concMaterial : public UniaxialMaterial
{
 public:
  RCFT_concMaterial(int tag, double fc, double depth, double thickness, double fy, double es);
  RCFT_concMaterial();    
  ~RCFT_concMaterial();
  
  int setTrialStrain(double strain, double strainRate = 0.0); 
  double getStrain(void);
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void);
  
  int commitState(void);
  int revertToLastCommit(void);    
  int revertToStart(void);        
  
  UniaxialMaterial *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);
  
 private:
  double Ttangent;
  double Tstrain;
  double Tstress;
  double Trule;
  double Frule, CFrule;
  int Cstart_flag;
  int Tstart_flag;
  int crack, Ccrack;
  int direction, Cdirection;
  int cp;
  int itr;
  double Fc;
  double Fc_n;
  double Fc_p;
  double nn;
  double np;
  double Ec;
  double ec_n;
  double ec_p;
  double Ec_sft;
  double Fc_res_n;
  double ec_res_n;
  double D;
  double t;
  double Fy;
  double Es;
  double eo_n, Ceo_n;
  double eo_p, Ceo_p;
  double funld_n, Cfunld_n;
  double eunld_n, Ceunld_n;
  double Esec_n, CEsec_n;
  double Epl_n, CEpl_n;
  double df_n, Cdf_n;
  double de_n, Cde_n;
  double epl_n, Cepl_n;
  double fnew_n, Cfnew_n;
  double Enew_n, CEnew_n;
  double ere_n, Cere_n;
  double fre_n, Cfre_n;
  double Ere_n, CEre_n;
  double fnew_str_n, Cfnew_str_n;
  double Enew_str_n, CEnew_str_n;
  double ere_str_n, Cere_str_n;
  double fre_str_n, Cfre_str_n;
  double Ere_str_n, CEre_str_n;
  double deo_n, Cdeo_n;
  double funld_p, Cfunld_p;
  double eunld_p, Ceunld_p;
  double Esec_p, CEsec_p;
  double Epl_p, CEpl_p;
  double df_p, Cdf_p;
  double de_p, Cde_p;
  double epl_p, Cepl_p;
  double fnew_p, Cfnew_p;
  double Enew_p, CEnew_p;
  double ere_p, Cere_p;
  double fre_p, Cfre_p;
  double Ere_p, CEre_p;
  double xu_n, Cxu_n;
  double xu_p, Cxu_p;
  double deo_p, Cdeo_p;
  double eunld_p8, funld_p8, Eunld_p8;
  double Ceunld_p8, Cfunld_p8, CEunld_p8;
  double eunld_n7, funld_n7, Eunld_n7;
  double Ceunld_n7, Cfunld_n7, CEunld_n7;
  double ea, Cea;
  double fa, Cfa;
  double Ea, CEa;
  double er, Cer;
  double fr, Cfr;
  double Er, CEr;
  double eb, Ceb;
  double fb, Cfb;
  double Eb, CEb;
  double fnew_str_p, Cfnew_str_p;
  double Enew_str_p, CEnew_str_p;
  double ere_str_p, Cere_str_p;
  double fre_str_p, Cfre_str_p;
  double Ere_str_p, CEre_str_p;
  double Cstrain;
  double Cstress;
  double Crule;
  double Ctangent;
  double r_n;
  double r_p;
  double xcr_n;
  double xcr_p;
  double ecr;
  int history_n, Chistory_n;
  int history_p, Chistory_p;
  double Cea1, ea1, Ceb1, eb1;
  double Cfa1, fa1, Cfb1, fb1; 
  double CEa1, Ea1, CEb1, Eb1;
  double Cer1, er1, Cfr1, fr1;
  double CEr1, Er1;
  int Urule; 
  int CUrule;
  double stress_envlp_p(double strn);
  double tangent_envlp_p(double strn);
  double stress_envlp_n(double strn);
  double tangent_envlp_n(double strn);
  double stress_tran(double strn, double ei, double fi, double Ei, double ef, double ff, double Ef);
  double tangent_tran(double strn, double ei, double fi, double Ei, double ef, double ff, double Ef);
  double stress_tran3(double strn, double ei, double fi, double Ei, double ef, double ff, double Ef);
  double tangent_tran3(double strn, double ei, double fi, double Ei, double ef, double ff, double Ef);
  double stress_tran_p(double strn, double ei, double fi, double Ei, double ef, double ff, double Ef);
  double tangent_tran_p(double strn, double ei, double fi, double Ei, double ef, double ff, double Ef);
  void commitStatevar(void);
  void backtocommitStatevar(void);
};
#endif



