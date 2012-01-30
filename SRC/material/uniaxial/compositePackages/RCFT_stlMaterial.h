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
// $Date: 2007-09-21 15:28:43 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/uniaxial/RCFT_stlMaterial.h,v $

// Written: MHS
// Created: Aug 2001

// Modified by : Cenk Tort - tort0008@umn.edu
// University of Minnesota - Civil Engineering Department
// Date : Wed Jul 23 17:32:23 EDT 2003
// Description: This file contains the class definition for
// cyclic uniaxial concrete stress-strain relationship for
// RCFT members.


#ifndef RCFT_stlMaterial_h
#define RCFT_stlMaterial_h

#include <UniaxialMaterial.h>

class RCFT_stlMaterial : public UniaxialMaterial
{
 public:
  RCFT_stlMaterial(int tag, double Fy, double Fu, double es, double depth, double thickness, double Eplo);
  RCFT_stlMaterial();    
  ~RCFT_stlMaterial();
  
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
  void commitStatevar(void);
  void backtocommitStatevar(void);
  double fy, fu, Ee, D, t;
  double Ksft, eplbf, Rbso, Epoi, alfa, a, bb, c, w, ksi, e, fE;
  double Rlso, CRlso, Rls, CRls, Tls_p, Cls_p, Tls_n, Cls_n, Tbs_p, Cbs_p;
  double Tbs_n, Cbs_n, Tmem_p, Cmem_p, Tmem_n, Cmem_n, Tvbs_p, Cvbs_p, Tvbs_n, Cvbs_n;  
  int elastic, Celastic, plastic, Cplastic, memory, Cmemory, lb, Clb, elb, Celb, repeat, ld, Cld;
  double ep, Cep, ep_ref, Cep_ref, elb_ref, Celb_ref, epmin, Cepmin, epmax, Cepmax;
  double ebar_p, Cebar_p, Ep, CEp, Epo, CEpo, W, CW, delta, Cdelta, delta_p, Cdelta_p;
  double delta_y, Cdelta_y, delta_in, Cdelta_in, deltap_in, Cdeltap_in, h, Ch, cbs, Ccbs;
  double elbf, Celbf, slbf, Cslbf, Trule, Crule, Ttangent, Ctangent, Tstress, Cstress;
  double Tstrain, Cstrain, strs, strn, tgnt;
  double Cstrs, Ctgnt;
  double epo;
  double strain_inc;
  double eunld, funld;
  double Ceunld, Cfunld;
  double eresp_n, fres_n;
  double Cmax_strs, Tmax_strs;
  int buckled;
  double lbstrn;
  double lbstrs;
  double lbW;
};

#endif



