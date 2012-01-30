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
                                                                        
// $Revision: 1.2 $
// $Date: 2008-09-16 20:34:23 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/uniaxial/CCFT_stlMaterial.h,v $

// Written: MHS
// Created: Aug 2001

// Modified by : Cenk Tort - tort0008@umn.edu
// University of Minnesota - Civil Engineering Department
// Date : Wed Jul 23 17:32:23 EDT 2003
// Description: This file contains the class definition for
// cyclic uniaxial concrete stress-strain relationship for
// RCFT members.
// Further Modified by: Mark Denavit
// University of Illinois at Urbana-Champaign
// for circular CFT (CCFT) members

#ifndef CCFT_stlMaterial_h
#define CCFT_stlMaterial_h

#include <UniaxialMaterial.h>

class CCFT_stlMaterial : public UniaxialMaterial
{
 public:
  CCFT_stlMaterial(int tag, double Fy, double Fu, double es, double depth, double thickness, double Eplo);
  CCFT_stlMaterial();    
  ~CCFT_stlMaterial();
  
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
  void commitStatevar(void); //function to commit all state variables en masse
  void backtocommitStatevar(void); //function to revert all state variables to last commit en masse
  int elastic, Celastic; //flag for type of loading: 1=elastic 0=plastic
  int plastic, Cplastic; //flag for type of loading: 1=plastic 0=elastic
  int memory, Cmemory; //flag if memory surface has been breached: 1=NOT breached 0=IS breached
  int lb, Clb; //flag if local buckling is occuring: 1=IS occuring 0=NOT occuring
  int elb, Celb; //does not appear to be used
  int repeat; //does not appear to be used
  int ld, Cld; //flag for direction of loading: 1=positive 0=negative 
  int buckled; //flag if buckled: 1=?? 0=??
  double fy, fu, Ee; //user input variables: steel yield strength, steel ultimate strength, steel elastic modulus
  double D, t; //user input variables: tube outer diamter, tube thickness
  double epo; //user input variable: inital plastic strain
  double Ksft; //softening slope of steel 
  double eplbf; //plastic strain at initiation of local buckling
  double Rbso, Epoi, alfa, a, bb, c, w, ksi, e, fE; //parameters from the Japanese material model
  double Rlso, CRlso; //initial size of the loading surface
  double Rls, CRls; //current size of the loading surface
  double Tls_p, Cls_p, Tls_n, Cls_n; //current positive and negative values of the loading surface
  double Tbs_p, Cbs_p, Tbs_n, Cbs_n; //current positive and negative values of the bounding surface
  double Tmem_p, Cmem_p, Tmem_n, Cmem_n; //current positive and negative values of the memory surface
  double Tvbs_p, Cvbs_p, Tvbs_n, Cvbs_n; //current positive and negative values of the virtual bounding surface
  double ep, Cep; //plastic strain
  double ep_ref, Cep_ref; //reference plastic strain for local buckling
  double elb_ref, Celb_ref; //does not appear to be used
  double epmin, Cepmin, epmax, Cepmax; //minimum and maximum plastic strains experienced by the fiber
  double ebar_p, Cebar_p; //range of plastic strain experienced by the fiber (epmax-epmin)
  double Ep, CEp; //current plastic stiffness
  double Epo, CEpo; //slope of the current bounding line
  double W, CW; //accumulated plastic work
  double delta, Cdelta; //distance between bounding line and loading point
  double delta_p, Cdelta_p; //does not appear to be used
  double delta_y, Cdelta_y; //distance between bounding line and virtual bounding line
  double delta_in, Cdelta_in; //value of delta at the inital yield state in the current loading path
  double deltap_in, Cdeltap_in; //value of distance between virtual bounding line and loading point at the inital yield state in the current loading path
  double h, Ch; //shape parameter for determination of Ep
  double cbs, Ccbs; //does not appear to be used
  double elbf, Celbf; //does not appear to be used
  double slbf, Cslbf; //does not appear to be used
  double Trule, Crule; //flag for type of loading: 1=elastic 2=plastic
  double Ttangent, Ctangent; //current tangent modulus of the fiber (both elastic and plastic)
  double Tstress, Cstress; //current stress in the fiber
  double Tstrain, Cstrain; //current total strain in the fiber
  double strs, Cstrs; //appears to be identical to Tstress, Cstress 
  double strn; //does not appear to be used
  double tgnt, Ctgnt; //appears to be identical to Ttangent, Ctangent 
  double strain_inc; //increment in strain for the current step
  double eunld, Ceunld, funld, Cfunld; //unloading stress and strain
  double eresp_n; //strain at initiation of constant residual strength
  double fres_n; //constant residual strength after local buckling
  double Cmax_strs, Tmax_strs; //maximum stress attained by the material (for size of memory surface)
  double lbstrn, lbstrs; //strain and stress at initiation of local buckling
  double lbW; //parameter for determination of the bounding surface in the local buckled range
};

#endif



