C
C  This file is part of MUMPS 4.8.4, built on Mon Dec 15 15:31:38 UTC 2008
C
C
C  This version of MUMPS is provided to you free of charge. It is public
C  domain, based on public domain software developed during the Esprit IV
C  European project PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL.
C  Since this first public domain version in 1999, the developments are
C  supported by the following institutions: CERFACS, ENSEEIHT-IRIT, and
C  INRIA.
C
C  Main contributors are Patrick Amestoy, Iain Duff, Abdou Guermouche,
C  Jacko Koster, Jean-Yves L'Excellent, and Stephane Pralet.
C
C  Up-to-date copies of the MUMPS package can be obtained
C  from the Web pages:
C  http://mumps.enseeiht.fr/  or  http://graal.ens-lyon.fr/MUMPS
C
C
C   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
C   EXPRESSED OR IMPLIED. ANY USE IS AT YOUR OWN RISK.
C
C
C  User documentation of any code that uses this software can
C  include this complete notice. You can acknowledge (using
C  references [1], [2], and [3]) the contribution of this package
C  in any scientific publication dependent upon the use of the
C  package. You shall use reasonable endeavours to notify
C  the authors of the package of this publication.
C
C   [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent,
C   Multifrontal parallel distributed symmetric and unsymmetric solvers,
C   in Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000).
C
C   [2] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
C   A fully asynchronous multifrontal solver using distributed dynamic
C   scheduling, SIAM Journal of Matrix Analysis and Applications,
C   Vol 23, No 1, pp 15-41 (2001).
C
C   [3] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
C   S. Pralet, Hybrid scheduling for the parallel solution of linear
C   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).
C
      INTEGER XXI, XXR, XXS, XXN, XXP
      PARAMETER(XXI=0,XXR=1,XXS=2,XXN=3,XXP=4)
      INTEGER XXVADDRL, XXWRITTENL
      PARAMETER(XXVADDRL=5,XXWRITTENL=7)
      INTEGER XXNDIAG2W  
      PARAMETER(XXNDIAG2W=8)
      INTEGER XXVADDRU, XXWRITTENU
      PARAMETER(XXVADDRU=9,XXWRITTENU=11)
      INTEGER XSIZE_IC, XSIZE_OOC_SYM, XSIZE_OOC_UNSYM
      INTEGER XSIZE_OOC_NOPANEL 
      PARAMETER (XSIZE_IC=5,XSIZE_OOC_SYM=9,XSIZE_OOC_UNSYM=12,
     *           XSIZE_OOC_NOPANEL=7)
      INTEGER IXSZ
      PARAMETER(IXSZ= 222)    
      INTEGER S_CB1COMP
      PARAMETER (S_CB1COMP=314)
      INTEGER S_ACTIVE, S_ALL, S_NOLCBCONTIG,
     *        S_NOLCBNOCONTIG, S_NOLCLEANED,
     *        S_NOLCBNOCONTIG38, S_NOLCBCONTIG38,
     *        S_NOLCLEANED38, C_FINI
      PARAMETER(S_ACTIVE=400, S_ALL=401, S_NOLCBCONTIG=402,
     *          S_NOLCBNOCONTIG=403, S_NOLCLEANED=404,
     *          S_NOLCBNOCONTIG38=405, S_NOLCBCONTIG38=406,
     *          S_NOLCLEANED38=407,C_FINI=1)
      INTEGER S_FREE, S_NOTFREE
      PARAMETER(S_FREE=54321,S_NOTFREE=-123456)
      INTEGER TOP_OF_STACK
      PARAMETER(TOP_OF_STACK=-999999)
      INTEGER XTRA_SLAVES_SYM, XTRA_SLAVES_UNSYM
      PARAMETER(XTRA_SLAVES_SYM=3, XTRA_SLAVES_UNSYM=1)
         INTEGER S_ROOT2SON_CALLED, S_REC_CONTSTATIC, 
     &  S_ROOTBAND_INIT
         PARAMETER(S_ROOT2SON_CALLED=-341,S_REC_CONTSTATIC=1,
     &             S_ROOTBAND_INIT=0)
