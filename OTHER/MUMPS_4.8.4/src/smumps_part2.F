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
      SUBROUTINE SMUMPS_267( BUFR, LBUFR, LBUFR_BYTES,
     *                           N, ID, MYID, STEP,
     *                           PTLUST_S, 
     *                           K10, K28, SLAVEF, IW, LIW,
     *                           A, LA, COMM, XSIZE)
      IMPLICIT NONE
      INTEGER LBUFR, LBUFR_BYTES, XSIZE
      INTEGER K10, K28, LA, LIW, COMM
      INTEGER ID, MYID, SLAVEF, N 
      INTEGER BUFR( LBUFR ) 
      INTEGER STEP(N), PTLUST_S( K28 ), IW( LIW )
      REAL A( LA )
      INCLUDE 'mpif.h'
      INTEGER POSITION, IOLDPS, IERR
      INTEGER NELIM, LCONT, NROW, NPIV, IAPOS, NSLAVES, ISLAVE
      INTEGER INODE, NCOL_L, NROW_L, NCB, BLSIZE, LASTBL
      INTEGER SIZELUMASTER, SIZELUSLAVE, NFRONT
      INCLUDE 'mumps_headers.h'
      ISLAVE = ID - MYID
      IF ( ISLAVE .LT. 1 ) ISLAVE = ISLAVE + SLAVEF 
      POSITION = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 INODE, 1, MPI_INTEGER,
     *                 COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 NCOL_L, 1, MPI_INTEGER,
     *                 COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 NROW_L, 1, MPI_INTEGER,
     *                 COMM, IERR )
      IOLDPS = PTLUST_S(STEP( INODE ))
      LCONT  = IW( IOLDPS +XSIZE)
      NELIM  = IW( IOLDPS + 1 +XSIZE)
      NROW   = IW( IOLDPS + 2 +XSIZE)
      NPIV   = IW( IOLDPS + 3 +XSIZE)
      IAPOS  = IW( IOLDPS + 4 +XSIZE)
      NSLAVES= IW( IOLDPS + 5 +XSIZE)
      NFRONT = LCONT + NPIV
      NCB    = NFRONT - NROW
      IF ( NSLAVES .NE. min(NCB/K10,SLAVEF-1)) THEN
        WRITE(*,*) 'Process_factor: error1. Wrong number of slaves'
        WRITE(*,*) 'INODE,NFRONT,NROW,NSLAVES=',
     * INODE,NFRONT,NCB,NSLAVES
        stop
      END IF
      BLSIZE = NCB / NSLAVES
      LASTBL = mod(NCB,NSLAVES) + BLSIZE
      SIZELUMASTER = ( LCONT + NROW ) * NPIV
      SIZELUSLAVE  = NPIV * BLSIZE
         IF ( NROW_L .NE. LASTBL .AND. NROW_L .NE. BLSIZE ) THEN
           WRITE(*,*) 'Process_factor: error2. wrong block size'
           CALL MUMPS_ABORT()
         END IF
         IF ( NSLAVES .eq. 0 ) THEN
           WRITE(*,*) 'Process_factor: error3'
           CALL MUMPS_ABORT()
         END IF
         IF ( NPIV .NE. NCOL_L ) THEN
           WRITE(*,*) 'Process_factor: error4'
           CALL MUMPS_ABORT()
         END IF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 A(IAPOS + SIZELUMASTER + ( ISLAVE - 1 )
     *                 * SIZELUSLAVE),
     *                 NCOL_L * NROW_L,
     *                 MPI_REAL,
     *                 COMM, IERR )
      RETURN
      END SUBROUTINE SMUMPS_267
      SUBROUTINE SMUMPS_152(SSARBR, MYID, N, IPOSBLOCK,
     *       RPOSBLOCK,
     *       IW, LIW,
     *       LRLU, LRLUS, IPTRLU,
     *       IWPOSCB, LA, KEEP, KEEP8, IN_PLACE_STATS
     $     )
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INTEGER IPOSBLOCK, RPOSBLOCK,
     *         LIW, LRLU, LRLUS,
     *        IPTRLU, IWPOSCB, N, LA
      LOGICAL IN_PLACE_STATS
      INTEGER IW( LIW ), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER MYID
      LOGICAL SSARBR
      INTEGER SIZFI_BLOCK, SIZFR_BLOCK, SIZFI, SIZFR
      INTEGER IPOSSHIFT
      INTEGER SIZFR_BLOCK_EFF, SIZEHOLE, MEM_INC
      INCLUDE 'mumps_headers.h'
      IPOSSHIFT = IPOSBLOCK + KEEP(IXSZ)
      SIZFI_BLOCK=IW(IPOSBLOCK+XXI)
      SIZFR_BLOCK=IW(IPOSBLOCK+XXR)
      IF (KEEP(216).eq.3) THEN
        SIZFR_BLOCK_EFF=SIZFR_BLOCK
      ELSE
        CALL SMUMPS_628( IW(IPOSBLOCK),
     *                     LIW-IPOSBLOCK+1,
     *                     SIZEHOLE, KEEP(IXSZ))
        SIZFR_BLOCK_EFF=SIZFR_BLOCK-SIZEHOLE
      ENDIF
      IF ( IPOSBLOCK .eq. IWPOSCB + 1 ) THEN
         IPTRLU  = IPTRLU  + SIZFR_BLOCK
         IWPOSCB = IWPOSCB + SIZFI_BLOCK
         LRLU    = LRLU  + SIZFR_BLOCK
         IF (.NOT. IN_PLACE_STATS) THEN
           LRLUS   = LRLUS + SIZFR_BLOCK_EFF
         ENDIF
      MEM_INC = -SIZFR_BLOCK_EFF
      IF (IN_PLACE_STATS) THEN
        MEM_INC= 0
      ENDIF
      CALL SMUMPS_471(SSARBR,.FALSE.,
     *         LA-LRLUS,0,MEM_INC,KEEP,KEEP8,LRLU)
 90      IF ( IWPOSCB .eq. LIW ) GO TO 100
         IPOSSHIFT = IWPOSCB + KEEP(IXSZ)
         SIZFI = IW( IWPOSCB+1+XXI )
         SIZFR = IW( IWPOSCB+1+XXR )
         IF ( IW( IWPOSCB+1+XXS ) .EQ. S_FREE ) THEN
              IPTRLU  = IPTRLU + SIZFR
              LRLU    = LRLU + SIZFR
              IWPOSCB = IWPOSCB + SIZFI
              GO TO 90
         ENDIF
 100     CONTINUE
         IW( IWPOSCB+1+XXP)=TOP_OF_STACK
      ELSE
         IW( IPOSBLOCK +XXS)=S_FREE
         IF (.NOT. IN_PLACE_STATS) LRLUS = LRLUS + SIZFR_BLOCK_EFF
      CALL SMUMPS_471(SSARBR,.FALSE.,
     *            LA-LRLUS,0,-SIZFR_BLOCK_EFF,KEEP,KEEP8,LRLU)
      END IF
      RETURN
      END SUBROUTINE SMUMPS_152
      SUBROUTINE SMUMPS_144( COMM_LOAD, ASS_IRECV, 
     *           N, INODE, FPERE, IW, LIW, A, LA,
     *           UU, NOFFW,
     *           NPVW,
     *             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, COMP,
     *             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP, PIMASTER,
     *             PAMASTER,
     *             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     *             OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL,
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE, AVOID_DELAYED,
     *              DKEEP,PIVNUL_LIST,LPN_LIST)
      USE SMUMPS_OOC
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, INODE, FPERE, LIW, LA, NOFFW, NPVW
      INTEGER IW( LIW )
      REAL A( LA )
      REAL UU, SEUIL
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER COMM, MYID, LBUFR, LBUFR_BYTES
      INTEGER LPTRAR, NELT
      INTEGER ICNTL(40), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER NBFIN, SLAVEF, 
     *        IFLAG, IERROR, LEAF, LPOOL
      INTEGER POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, 
     *        COMP 
      INTEGER FRTPTR( N + 1 ), FRTELT( NELT )
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL),
     *        ITLOC(N), FILS(N),
     *        PTRARW(LPTRAR), PTRAIW(LPTRAR),
     *        ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER  INTARR(max(1,KEEP(14)))
      INTEGER PTRIST(KEEP(28)), PTRAST(KEEP(28)),
     *        PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     *        STEP(N),
     *        
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28)),
     *        NSTK_S(KEEP(28)), NBPROCFILS(KEEP(28)),
     *        PROCNODE_STEPS(KEEP(28))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      DOUBLE PRECISION OPASSW, OPELIW
      REAL DBLARR(max(1,KEEP(13)))
      LOGICAL AVOID_DELAYED
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      REAL DKEEP(30)
      INTEGER INOPV, IFINB, NFRONT, NPIV, IBEGKJI, NBOLKJ,
     &        NBTLKJ, IBEG_BLOCK
      INTEGER NASS, NEL1, IEND, IOLDPS, POSELT, dummy, allocok
      LOGICAL LASTBL 
      REAL UUTEMP
      INTEGER , ALLOCATABLE, DIMENSION ( : ) :: IPIV
      INTEGER LAFAC, LIWFAC, STRAT, TYPEFile, LNextPiv2beWritten, 
     &        UNextPiv2beWritten, IFLAG_OOC,
     &        PP_FIRST2SWAP_L, PP_FIRST2SWAP_U,
     &        PP_LastPIVRPTRFilled_L, PP_LastPIVRPTRFilled_U 
      TYPE(IO_BLOCK) :: MonBloc 
      INCLUDE 'mumps_headers.h'
      EXTERNAL SMUMPS_224, SMUMPS_233, 
     *         SMUMPS_225, SMUMPS_232,
     *         SMUMPS_294,
     *         SMUMPS_44
      LOGICAL STATICMODE
      REAL SEUIL_LOC
      INOPV = 0
      SEUIL_LOC = SEUIL
      IF(KEEP(97) .EQ. 0) THEN
         STATICMODE = .FALSE.
      ELSE
         STATICMODE = .TRUE.
      ENDIF
      IF (AVOID_DELAYED) THEN
         STATICMODE = .TRUE.
         UUTEMP=UU
         SEUIL_LOC = max(SEUIL,epsilon(SEUIL))
      ELSE
         UUTEMP=UU
      ENDIF
      IBEG_BLOCK=1
      dummy  = 0
      IOLDPS = PTLUST_S(STEP( INODE ))
      POSELT = PTRAST(STEP( INODE ))
      NFRONT = IW(IOLDPS+KEEP(IXSZ))
      NASS   = iabs(IW(IOLDPS+2+KEEP(IXSZ)))
      IF (NASS .GT. KEEP(3)) THEN
        NBOLKJ = min( KEEP(6), NASS )
      ELSE
        NBOLKJ = min( KEEP(5),NASS )
      ENDIF
      NBTLKJ = NBOLKJ
      ALLOCATE( IPIV( NASS ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(*,*) MYID,' : FACTO_NIV2 :failed to allocate ',NASS,
     * ' integers'
        IFLAG  = -13
        IERROR =NASS
        GO TO 490
      END IF
      IF (KEEP(201).EQ.1) THEN 
          LAFAC     = IW(IOLDPS+XXR)
          LIWFAC    = IW(IOLDPS+XXI)
          TYPEFile  = TYPEF_U
          LNextPiv2beWritten = 1 
          UNextPiv2beWritten = 1 
          PP_FIRST2SWAP_L = LNextPiv2beWritten 
          PP_FIRST2SWAP_U = UNextPiv2beWritten 
          MonBloc%LastPanelWritten_L = 0 
          MonBloc%LastPanelWritten_U = 0 
          MonBloc%INODE    = INODE
          MonBloc%MASTER   = .TRUE.
          MonBloc%Typenode = 2
          MonBloc%NROW     = NASS
          MonBloc%NCOL     = NFRONT
          MonBloc%NFS      = NASS
          MonBloc%Last     = .FALSE.   
          MonBloc%LastPiv  = -68877    
          NULLIFY(MonBloc%INDICES)
      ENDIF
 50   CONTINUE
      IBEGKJI = IBEG_BLOCK
      CALL SMUMPS_224(NFRONT,NASS,IBEGKJI, NASS, IPIV,
     *                N,INODE,IW,LIW,A,LA,INOPV,NOFFW,
     *                IFLAG,IOLDPS,POSELT,UU,SEUIL_LOC,KEEP,KEEP8,
     *                 DKEEP(1),PIVNUL_LIST(1),LPN_LIST, 
     *     PP_FIRST2SWAP_L,  MonBloc%LastPanelWritten_L,
     *     PP_LastPIVRPTRFilled_L,
     *     PP_FIRST2SWAP_U,  MonBloc%LastPanelWritten_U,
     *     PP_LastPIVRPTRFilled_U)
      IF (IFLAG.LT.0) GOTO 490   
      IF (INOPV.EQ.1) THEN
         IF(STATICMODE) THEN
            INOPV = -1
            GOTO 50
         ENDIF
      ENDIF
      IF (INOPV.GE.1) THEN
          LASTBL = (INOPV.EQ.1)
          IEND = IW(IOLDPS+1+KEEP(IXSZ))
          CALL SMUMPS_294( COMM_LOAD, ASS_IRECV, 
     *             N, INODE, FPERE, IW, LIW, 
     *             IOLDPS, POSELT, A, LA, NFRONT,
     *             IBEGKJI, IEND, IPIV, NASS,LASTBL, dummy, 
     *
     *             COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, COMP,
     *             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     *             PIMASTER, PAMASTER,
     *             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     *             OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) GOTO 500
      ENDIF
      IF (INOPV.EQ.1) GO TO 500
      IF (INOPV.EQ.2) THEN
         CALL SMUMPS_233(IBEG_BLOCK,NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *            IOLDPS,POSELT,NBOLKJ,NBTLKJ,KEEP(4),KEEP(IXSZ))
         GOTO 50
      ENDIF
      NPVW = NPVW + 1
      IF (NASS.LE.1) THEN
        IFINB = -1
      ELSE
         CALL SMUMPS_225(IBEG_BLOCK,
     *             NFRONT, NASS, N,INODE,IW,LIW,A,LA,
     *             IOLDPS,POSELT,IFINB,
     *             NBTLKJ,KEEP(4),KEEP(IXSZ))
      ENDIF
      IW(IOLDPS+1+KEEP(IXSZ)) = IW(IOLDPS+1+KEEP(IXSZ)) + 1
       IF (IFINB.EQ.0) GOTO 50
       IF ((IFINB.EQ.1).OR.(IFINB.EQ.-1)) THEN
          LASTBL = (IFINB.EQ.-1) 
          IEND = IW(IOLDPS+1+KEEP(IXSZ))
          CALL SMUMPS_294(COMM_LOAD, ASS_IRECV, 
     *             N, INODE, FPERE, IW, LIW, 
     *             IOLDPS, POSELT, A, LA, NFRONT, 
     *             IBEGKJI, IEND, IPIV, NASS, LASTBL, dummy, 
     *
     *             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, COMP,
     *             PTRIST, PTRAST, PTLUST_S, PTRFAC,
     *             STEP, PIMASTER, PAMASTER,
     *             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     *             OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) GOTO 500
       ENDIF
       IF (IFINB.EQ.(-1)) GOTO 500
       NPIV   = IW(IOLDPS+1+KEEP(IXSZ))
       NEL1   = NASS - NPIV
      CALL SMUMPS_232(A,LA,
     *           NFRONT,NPIV,NASS,POSELT,NBTLKJ)
      IF (KEEP(201).EQ.1) THEN 
          STRAT            = STRAT_TRY_WRITE
          MonBloc%LastPiv  = NPIV
          TYPEFile         = TYPEF_BOTH_LU  
          CALL SMUMPS_688
     &          ( STRAT, TYPEFile, 
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC )
          IF (IFLAG_OOC < 0 ) IFLAG=IFLAG_OOC
          IF (IFLAG<0) RETURN
      ENDIF
      GO TO 50
 490  CONTINUE
      CALL SMUMPS_44( MYID, SLAVEF, COMM )
 500  CONTINUE
      DEALLOCATE( IPIV )
      IF (KEEP(201).EQ.1) THEN 
          STRAT        = STRAT_WRITE_MAX   
          MonBloc%Last = .TRUE.
          MonBloc%LastPiv  = IW(IOLDPS+1+KEEP(IXSZ))
          TYPEFile     = TYPEF_BOTH_LU  
          CALL SMUMPS_688
     &          ( STRAT, TYPEFile, 
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC )
          IF (IFLAG_OOC < 0 ) IFLAG=IFLAG_OOC
          IF (IFLAG<0) RETURN
          CALL SMUMPS_644(IWPOS, 
     &      IOLDPS, IW, LIW, MonBloc , NFRONT, KEEP)
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_144
      SUBROUTINE SMUMPS_176( COMM_LOAD, ASS_IRECV, 
     *    root, FRERE, IROOT, 
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     *    PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      USE SMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER IROOT
      INTEGER ICNTL( 40 ), KEEP( 500 )
      INTEGER*8 KEEP8(150)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER POSFAC,IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS
      INTEGER N, LIW, LA
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     * PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S( KEEP(28) ), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF, NBFIN
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER ITLOC( N ), FILS( N )
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND(KEEP(28)), FRERE(KEEP(28))
      INTEGER INTARR( max(1,KEEP(14)) )
      REAL DBLARR( max(1,KEEP(13)) )
      INTEGER NELIM, NB_CONTRI_GLOBAL, NUMORG, 
     *        NFRONT, IROW, JCOL, PDEST, HF, IOLDPS,
     *        IN, DEB_ROW, ILOC_ROW, IFSON, ILOC_COL,
     *        IPOS_SON, NELIM_SON, NSLAVES_SON, HS,
     *        IROW_SON, ICOL_SON, ISLAVE, IERR, 
     *        NELIM_SENT, IPOS_STATREC
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      INCLUDE 'mumps_headers.h'
      INCLUDE 'mumps_tags.h'
      NB_CONTRI_GLOBAL = KEEP(41)
      NUMORG    = root%ROOT_SIZE
      NELIM     = KEEP(42)
      NFRONT    = NUMORG + KEEP(42)
      DO IROW = 0, root%NPROW - 1
        DO JCOL = 0, root%NPCOL - 1
            PDEST = IROW * root%NPCOL + JCOL
          IF ( PDEST .NE. MYID ) THEN
           CALL SMUMPS_73(NFRONT, 
     *     NB_CONTRI_GLOBAL, PDEST, COMM, IERR)
              if (IERR.lt.0) then
                write(6,*) ' error detected by ',
     *          'SMUMPS_73'
                stop
               endif
           ENDIF
        END DO
      END DO
      CALL  SMUMPS_270( NFRONT,
     *    NB_CONTRI_GLOBAL, root,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     *    PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM, COMM_LOAD,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND )
       IF (IFLAG < 0 ) RETURN
      HF = 6 + KEEP(IXSZ)
      IOLDPS = PTLUST_S(STEP(IROOT))
      IN = IROOT
      DEB_ROW = IOLDPS + HF
      ILOC_ROW    = DEB_ROW
      DO WHILE (IN.GT.0) 
       IW(ILOC_ROW)           = IN
       IW(ILOC_ROW+NFRONT)    = IN
       ILOC_ROW = ILOC_ROW + 1
       IN = FILS(IN)
      END DO
      IFSON = -IN
      ILOC_ROW    = IOLDPS + HF + NUMORG
      ILOC_COL    = ILOC_ROW + NFRONT
      IF ( NELIM.GT.0 ) THEN
        IN = IFSON
        DO WHILE (IN.GT.0)
          IPOS_SON  = PIMASTER(STEP(IN))
          IF (IPOS_SON .EQ. 0) GOTO 100
          NELIM_SON   = IW(IPOS_SON+1+KEEP(IXSZ))
              if (NELIM_SON.eq.0) then
                write(6,*) ' error 1 in process_last_rtnelind'
                CALL MUMPS_ABORT()
              endif
          NSLAVES_SON = IW(IPOS_SON+5+KEEP(IXSZ))
          HS          = 6 + NSLAVES_SON + KEEP(IXSZ)
          IROW_SON    = IPOS_SON + HS
          ICOL_SON    = IROW_SON + NELIM_SON
          IW(ILOC_ROW:ILOC_ROW+NELIM_SON-1) = 
     *            IW(IROW_SON:IROW_SON+NELIM_SON-1) 
          IW(ILOC_COL:ILOC_COL+NELIM_SON-1) = 
     *            IW(ICOL_SON:ICOL_SON+NELIM_SON-1) 
          NELIM_SENT = ILOC_ROW - IOLDPS - HF + 1
          DO ISLAVE = 0,NSLAVES_SON
            IF (ISLAVE.EQ.0) THEN
             PDEST= MUMPS_275(STEP(IN),PROCNODE_STEPS,SLAVEF)
            ELSE
             PDEST = IW(IPOS_SON + 5 + ISLAVE+KEEP(IXSZ))
            ENDIF
            IF (PDEST.NE.MYID) THEN
             CALL SMUMPS_74(IN, NELIM_SENT,
     *        PDEST, COMM, IERR )
               if (IERR.lt.0) then
                write(6,*) ' error detected by ',
     *          'SMUMPS_73'
                CALL MUMPS_ABORT()
               endif
            ELSE
             CALL SMUMPS_271( COMM_LOAD, ASS_IRECV,
     *       IN, NELIM_SENT, root,
     *
     *       BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *       IWPOS, IWPOSCB, IPTRLU,
     *       LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST, 
     *       PTLUST_S, PTRFAC,
     *       PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *       IFLAG, IERROR, COMM,
     *       NBPROCFILS,
     *       IPOOL, LPOOL, LEAF,
     *       NBFIN, MYID, SLAVEF,
     *
     *       OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *       INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *       LPTRAR, NELT, FRTPTR, FRTELT, 
     *       ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
             IF ( ISLAVE .NE. 0 ) THEN
               IF (KEEP(50) .EQ. 0) THEN
                IPOS_STATREC = PTRIST(STEP(IN))+6+KEEP(IXSZ)
               ELSE
                IPOS_STATREC = PTRIST(STEP(IN))+8+KEEP(IXSZ)
               ENDIF
               IF (IW(IPOS_STATREC).EQ. S_REC_CONTSTATIC) THEN
                  IW(IPOS_STATREC) = S_ROOT2SON_CALLED
               ELSE
                CALL SMUMPS_626( N, IN, PTRIST, PTRAST,
     *          IW, LIW, A, LA, LRLU, LRLUS, IWPOSCB,
     *          IPTRLU, STEP, MYID, KEEP
     *        )
               ENDIF
             ENDIF
             IPOS_SON  = PIMASTER(STEP(IN))
            ENDIF
          END DO
          CALL  SMUMPS_152( .FALSE.,MYID,N, IPOS_SON, 
     *       PTRAST(STEP(IN)),
     *       IW, LIW,
     *       LRLU, LRLUS, IPTRLU,
     *       IWPOSCB, LA, KEEP,KEEP8, .FALSE.
     $         )
          ILOC_ROW = ILOC_ROW + NELIM_SON
          ILOC_COL = ILOC_COL + NELIM_SON
 100      CONTINUE
          IN = FRERE(STEP(IN))
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_176
      SUBROUTINE SMUMPS_268(MYID,BUFR, LBUFR, 
     *     LBUFR_BYTES,
     *     PROCNODE_STEPS, SLAVEF,
     *     IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS,
     *     N, IW, LIW, A, LA,
     *     PTRIST, PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S,
     *     COMP,
     *     IFLAG, IERROR, COMM, COMM_LOAD, NBPROCFILS,
     *     IPOOL, LPOOL, LEAF, KEEP,KEEP8, ND, FILS, FRERE, ITLOC, 
     *
     *     ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR
      INTEGER MYID
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER SLAVEF
      INTEGER IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS
      INTEGER N, LIW, LA
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER PTRIST(KEEP(28)), PTRAST(KEEP(28)), STEP(N),
     *        
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER PROCNODE_STEPS( KEEP(28) ), ITLOC( N )
      INTEGER COMP
      INTEGER NSTK_S( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM, COMM_LOAD
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER ND(KEEP(28)), FILS( N ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER POSITION, IFATH, ISON, NROW, NCOL, NELIM, 
     *        NSLAVES
      INTEGER NOREAL, NOINT, INIV2, NCOL_EFF
      DOUBLE PRECISION FLOP1
      INTEGER NBROWS_ALREADY_SENT, NBROWS_PACKET
      INTEGER NOREAL_PACKET
      LOGICAL PERETYPE2
      INCLUDE 'mumps_headers.h'
      INTEGER  MUMPS_330
      EXTERNAL MUMPS_330
      POSITION = 0
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *          IFATH, 1, MPI_INTEGER
     *        , COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        ISON , 1, MPI_INTEGER, 
     *        COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        NSLAVES, 1,
     *        MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *          NROW , 1, MPI_INTEGER
     *        , COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *          NCOL , 1, MPI_INTEGER
     *        , COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &          NBROWS_ALREADY_SENT, 1,
     &          MPI_INTEGER, COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &          NBROWS_PACKET, 1,
     &          MPI_INTEGER, COMM, IERR)
      IF ( NSLAVES .NE. 0 .and. KEEP(50).ne.0 ) THEN
        NCOL_EFF = NROW
      ELSE
        NCOL_EFF = NCOL
      ENDIF
      NOREAL_PACKET = NBROWS_PACKET * NCOL_EFF
      IF (NBROWS_ALREADY_SENT .EQ. 0) THEN
        NOINT = 6 + NROW + NCOL + NSLAVES + KEEP(IXSZ)
        NOREAL= NROW * NCOL_EFF
        CALL SMUMPS_22(.FALSE.,0,.FALSE.,.FALSE.,
     *   MYID,N,KEEP,KEEP8,IW,LIW,A,LA,
     *   LRLU, IPTRLU,IWPOS,IWPOSCB,
     *   PTRIST,PTRAST,STEP, PIMASTER, PAMASTER, ITLOC,
     *   NOINT, NOREAL, ISON, S_NOTFREE, .TRUE.,
     *   COMP, LRLUS, IFLAG, IERROR
     $     )
        IF ( IFLAG .LT. 0 ) THEN
          RETURN
        ENDIF
        PIMASTER(STEP( ISON )) = IWPOSCB + 1
        PAMASTER(STEP( ISON )) = IPTRLU  + 1
        IW( IWPOSCB + 1 + KEEP(IXSZ) ) = NCOL
        NELIM = NROW
        IW( IWPOSCB + 2 + KEEP(IXSZ) ) = NELIM
        IW( IWPOSCB + 3 + KEEP(IXSZ) ) = NROW
        IF ( NSLAVES .NE. 0 .and. KEEP(50).ne.0 ) THEN
          IW( IWPOSCB + 4 + KEEP(IXSZ) ) = NROW - NCOL
          IF ( NROW - NCOL .GE. 0 ) THEN
            WRITE(*,*) 'Error in PROCESS_MAITRE2:',NROW,NCOL
            CALL MUMPS_ABORT()
          END IF
        ELSE
          IW( IWPOSCB + 4 + KEEP(IXSZ) ) = 0
        END IF
        IW( IWPOSCB + 5 + KEEP(IXSZ) ) = 1
        IW( IWPOSCB + 6 + KEEP(IXSZ) ) = NSLAVES
        IF (NSLAVES.GT.0) THEN
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 IW( IWPOSCB + 7 + KEEP(IXSZ) ),
     *                 NSLAVES, MPI_INTEGER, COMM, IERR )
        ENDIF
        CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        IW(IWPOSCB + 7 + KEEP(IXSZ) + NSLAVES),
     *        NROW, MPI_INTEGER, COMM, IERR)
        CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        IW(IWPOSCB + 7 + KEEP(IXSZ) + NROW + NSLAVES),
     *        NCOL, MPI_INTEGER, COMM, IERR)
        IF ( ( KEEP(48).NE. 0 ).AND.(NSLAVES .GT. 0 )) THEN
          INIV2 = ISTEP_TO_INIV2 ( STEP(ISON) )
          CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        TAB_POS_IN_PERE(1,INIV2),
     *        NSLAVES+1, MPI_INTEGER, COMM, IERR)
          TAB_POS_IN_PERE(SLAVEF+2,INIV2) = NSLAVES
        ENDIF
      ENDIF
      IF (NOREAL_PACKET.GT.0) THEN
        CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        A(PAMASTER(STEP(ISON)) +
     *          NBROWS_ALREADY_SENT * NCOL_EFF),
     *        NOREAL_PACKET, MPI_REAL, COMM, IERR)
      ENDIF
      IF ( NBROWS_ALREADY_SENT + NBROWS_PACKET .EQ. NROW ) THEN
        PERETYPE2 = ( MUMPS_330(STEP(IFATH),PROCNODE_STEPS,
     *              SLAVEF) .EQ. 2 )
        NSTK_S( STEP(IFATH ))       = NSTK_S( STEP(IFATH) ) - 1
        IF ( NSTK_S( STEP(IFATH)) .EQ. 0 ) THEN
          CALL SMUMPS_507(N, IPOOL, LPOOL, PROCNODE_STEPS,
     *         SLAVEF, KEEP(28), KEEP(76), KEEP(80), KEEP(47),
     *         STEP, IFATH )
          IF (KEEP(47) .GE. 3) THEN
             CALL SMUMPS_500(
     $            IPOOL, LPOOL, 
     *            PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     *            MYID, STEP, N, ND, FILS )
          ENDIF
          CALL MUMPS_137( IFATH, N, PROCNODE_STEPS,
     *                            SLAVEF, ND,
     *                            FILS,FRERE, STEP, PIMASTER,
     *                            KEEP(28), KEEP(50),
     *                            FLOP1,IW, LIW, KEEP(IXSZ) )
          IF (IFATH.NE.KEEP(20))
     *    CALL SMUMPS_190(1, .FALSE., FLOP1, KEEP,KEEP8)
        END IF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_268
      SUBROUTINE SMUMPS_242(DATA, LDATA, MPITYPE, ROOT, COMMW, TAG,
     *SLAVEF)
      USE SMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR
      INTEGER LDATA, ROOT, COMMW, TAG, MPITYPE, SLAVEF
      INTEGER DEST
      INTEGER DATA(LDATA)
      DO 10 DEST = 0, SLAVEF - 1
        IF (DEST .NE. ROOT) THEN
          IF ( LDATA .EQ. 1 .and. MPITYPE .EQ. MPI_INTEGER ) THEN
            CALL SMUMPS_62( DATA(1), DEST, TAG, 
     *                                COMMW, IERR )
          ELSE
            WRITE(*,*) 'Error : bad argument to SMUMPS_242'
            CALL MUMPS_ABORT()
          END IF
        ENDIF
   10 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_242
      SUBROUTINE SMUMPS_44( MYID, SLAVEF, COMM )
      INTEGER MYID, SLAVEF, COMM
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER DUMMY (1)
      CALL SMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID,
     *                 COMM, TERREUR, SLAVEF )
      RETURN
      END SUBROUTINE SMUMPS_44
      SUBROUTINE SMUMPS_464( K34, K35 )
      INTEGER K34, K35
#if defined(t3e)
      K34 = 8
#else
      K34 = 4
#endif
      K35 = 4
#if defined(t3e)
      K35 = K35 * 2
#endif
      RETURN
      END SUBROUTINE SMUMPS_464
      SUBROUTINE SMUMPS_20( NPROCS, CNTL, ICNTL,
     *                    KEEP,KEEP8,
     *                    INFO, INFOG, RINFO, RINFOG, SYM, PAR,
     *                    DKEEP)
      REAL    DKEEP(30)
      REAL    CNTL(15), RINFO(20), RINFOG(20)
      INTEGER ICNTL(40), KEEP(500), SYM, PAR, NPROCS
      INTEGER INFO(40), INFOG(40)
      INTEGER*8 KEEP8(150)
C     Let $A_{preproc}$ be the preprocessed matrix to be factored (see
      KEEP(1:500) = 0
      KEEP8(1:150)= 0_8
      INFO(1:40)  = 0
      INFOG(1:40) = 0
      ICNTL(1:40) = 0
      RINFO(1:20) = 0.0E0
      RINFOG(1:20)= 0.0E0
      CNTL(1:15)  = 0.0E0
      DKEEP(1:30) = 0.0E0
      KEEP( 50 ) = SYM
      IF ( KEEP(50).NE.1 .and. KEEP(50).NE.2 ) KEEP( 50 ) = 0
      IF ( KEEP(50) .NE. 1 ) THEN
        CNTL(1)   = 0.01E0
      ELSE
        CNTL(1)   = 0.0E0
      END IF
      CNTL(2) = sqrt(epsilon(0.0E0))
      CNTL(3) = 0.0E0
      CNTL(4) = -1.0E0
      CNTL(5) = 0.0E0
      CNTL(6) = -1.0E0
      KEEP(46) = PAR
      IF ( KEEP(46) .NE. 0 .AND.
     *     KEEP(46) .NE. 1 ) THEN
           KEEP(46) = 1
      END IF
      ICNTL(1)  = 6
      ICNTL(2)  = 0
      ICNTL(3)  = 6
      ICNTL(4)  = 2
      ICNTL(5)  = 0
      IF (SYM.NE.1) THEN
       ICNTL(6)  = 7
      ELSE
       ICNTL(6)  = 0
      ENDIF
      ICNTL(7) = 7
      ICNTL(8)  = 77
      ICNTL(9)  = 1
      ICNTL(10)  = 0
      ICNTL(11)  = 0
      IF(SYM .EQ. 2) THEN
         ICNTL(12)  = 0
      ELSE
         ICNTL(12)  = 1
      ENDIF
      ICNTL(13) = 0
      IF (SYM.eq.1.AND.NPROCS.EQ.1) THEN 
        ICNTL(14) = 5  
      ELSE IF (NPROCS .GT. 4) THEN
        ICNTL(14) = 30
      ELSE
        ICNTL(14) = 20
      END IF
      ICNTL(15) = 0
      ICNTL(16) = 0
      ICNTL(17) = 0
      ICNTL(18) = 0
      ICNTL(19) = 0
      ICNTL(20) = 0
      ICNTL(21) = 0
      ICNTL(22) = 0
      ICNTL(23) = 0
      ICNTL(24) = 0
      ICNTL(27)= -8
      ICNTL(39) = 1
      ICNTL(40)  = 0 
      KEEP(12) = 0
      KEEP(11) = 2147483646
      KEEP(24) = 18
      KEEP(68) = 0
      KEEP(36) = 1
      KEEP(1) = 8
      KEEP(7)  = 150
      KEEP(8)  = 120
      KEEP(57) = 500
      KEEP(58) = 250
#if defined(t3e) || defined(sgi)
      KEEP(9)  = 400
      IF ( SYM .eq. 0 ) THEN
        KEEP(4)  = 32
        KEEP(3)  = 96
        KEEP(5)  = 16
        KEEP(6)  = 32
        KEEP(85) =  160
        KEEP(62) =  50
        IF (NPROCS.GE.128) KEEP(62)=200
        IF (NPROCS.GE.128) KEEP(9)=500
        IF (NPROCS.GE.256) KEEP(9)=600
      ELSE
        KEEP(4)  = 24
        KEEP(3) =  96
        KEEP(5) =  16
        KEEP(6) =  48
        KEEP(85) =  80
        KEEP(62) = 100
        IF (NPROCS.GE.128) KEEP(62)=150
        IF (NPROCS.GE.64) KEEP(9)=500
        IF (NPROCS.GE.128) KEEP(9)=600
        IF (NPROCS.GE.256) KEEP(9)=700
      END IF
#elif defined(SP_)
      IF ( SYM .eq. 0 ) THEN
        KEEP(4)  = 32
        KEEP(3)  = 96
        KEEP(5)  = 16
        KEEP(6)  = 32
        KEEP(9)  = 800
        KEEP(85) = 400
      ELSE
        KEEP(4)  = 24
        KEEP(3) =  96
        KEEP(5) =  16
        KEEP(6) =  48
        KEEP(9)  = 400
        KEEP(85) = 100
      END IF
      KEEP(62) = 150
#else
      IF ( SYM .eq. 0 ) THEN
        KEEP(4)  = 32
        KEEP(3)  = 96
        KEEP(5)  = 16
        KEEP(6)  = 32
        KEEP(9)  = 700
        KEEP(85) =  300
        KEEP(62) =  50
        IF (NPROCS.GE.128) KEEP(62)=200
        IF (NPROCS.GE.128) KEEP(9)=800
        IF (NPROCS.GE.256) KEEP(9)=900
      ELSE
        KEEP(4)  = 24 
        KEEP(3)  = 96  
        KEEP(5)  = 16
        KEEP(6)  = 48
        KEEP(9)  = 400
        KEEP(85) = 100
        KEEP(62) = 100
        IF (NPROCS.GE.128) KEEP(62)=150
        IF (NPROCS.GE.64) KEEP(9)=800
        IF (NPROCS.GE.128) KEEP(9)=900
      END IF
#endif
      KEEP(63) = 60
      KEEP(48) = 5
      KEEP(17) = 0
      CALL SMUMPS_464( KEEP(34), KEEP(35))
      KEEP(37) = max(800, 2*NPROCS)
      IF ( NPROCS > 256 ) THEN
        KEEP(39) = 10000
      ELSEIF ( NPROCS > 128 ) THEN
        KEEP(39) = 20000
      ELSEIF ( NPROCS > 64 ) THEN
        KEEP(39) = 40000
      ELSEIF ( NPROCS > 16 ) THEN
        KEEP(39) = 80000
      ELSE
        KEEP(39) = 160000
      END IF
      KEEP(40) = -1 - 456789
      KEEP(45) = 0
      KEEP(47) = 2
#if defined(SP_)
      KEEP( 51 )  = 70
#elif defined(valuesV416sgi)
      KEEP( 51 )  = 48
#elif defined(t3e) || defined(sgi)
      KEEP( 51 )  = 24
#else
      KEEP( 51 )  = 48
#endif
      KEEP(64) = 10
      KEEP(69) = 4
      KEEP(75) = 1
      KEEP(76) = 2
      KEEP(77) = 30
      IF (NPROCS.GT.4) THEN
          KEEP(78)=max(
     &       int(log(real(NPROCS))/log(real(2))) - 2 
     &       , 0         )
      ENDIF
      KEEP(210) = 2 
      KEEP(79) = -10
      KEEP(80) = 1
      KEEP(81) = 0
      KEEP(82) = 5
      KEEP(83) = min(8,NPROCS/4)
      KEEP(83) = max(min(4,NPROCS),max(KEEP(83),1))
      KEEP(86)=1
      KEEP(87)=0
      KEEP(88)=0
      KEEP(90)=1
      KEEP(91)=min(8, NPROCS)
      KEEP(91) = max(min(4,NPROCS),min(KEEP(83),KEEP(91)))
      IF(NPROCS.LT.48)THEN
         KEEP(102)=150
      ELSEIF(NPROCS.LT.128)THEN
         KEEP(102)=150
      ELSEIF(NPROCS.LT.256)THEN
         KEEP(102)=200
      ELSEIF(NPROCS.LT.512)THEN
         KEEP(102)=300
      ELSEIF(NPROCS.GE.512)THEN
         KEEP(102)=400
      ENDIF
#if defined(OLD_OOC_NOPANEL)
      KEEP(99)=0  
#else
      KEEP(99)=4  
#endif
      KEEP(100)=0
      KEEP(204)=0
      KEEP(205)=0
      KEEP(209)=-1 
      KEEP(104) = 16
      KEEP(107)=0
      KEEP(211)=2
      KEEP(213)=201
      KEEP(217)=0
      KEEP(215)=0
      KEEP(216)=1
      KEEP(218)=50
      KEEP(219)=1
      IF (KEEP(50).EQ.2) THEN
        KEEP(227)= max(2,32)
      ELSE
        KEEP(227)= max(1,32)
      ENDIF
      KEEP(231) = 1
      KEEP(232) = 3
      KEEP(233) = 0
      KEEP(239) = 1
      KEEP(240) = 10
      KEEP(241) = 0
      DKEEP(4) = -1.0E0
      DKEEP(5) = -1.0E0
      IF(NPROCS.LE.8)THEN
         KEEP(238)=12
      ELSE
         KEEP(238)=7
      ENDIF
      KEEP(234)= 1
      DKEEP(3)=-5.0D0
      RETURN
      END SUBROUTINE SMUMPS_20
      SUBROUTINE SMUMPS_195(N, NZ, IRN, ICN, LIW, IKEEP, PTRAR,
     *     IORD, NFSIZ, FILS, FRERE, LISTVAR_SCHUR, SIZE_SCHUR,
     *     ICNTL, INFO, KEEP,KEEP8, NSLAVES, SYM, PIV, id)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER N,NZ,LIW,IORD,SIZE_SCHUR, NSLAVES, SYM
      INTEGER PTRAR(N,4), NFSIZ(N), FILS(N), FRERE(N)
      INTEGER IKEEP(N,3)
      INTEGER  LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER INFO(40), ICNTL(40), KEEP(500)
      INTEGER*8 KEEP8(150)
      TYPE (SMUMPS_STRUC) :: id      
      INTEGER IRN(NZ), ICN(NZ)  
      INTEGER, DIMENSION(:), ALLOCATABLE :: IW
      INTEGER IERR
      INTEGER K,I,L1,L2,IWFR,NCMPA,LLIW, IN, IFSON
      INTEGER NEMIN, LP, MP, LDIAG, ITEMP, symmetry
      INTEGER MedDens, NBQD, AvgDens
      LOGICAL PROK, COMPRESS_SCHUR
      INTEGER NBBUCK
      INTEGER, DIMENSION(:), ALLOCATABLE :: HEAD
      INTEGER NUMFLAG
      INTEGER OPT_METIS_SIZE
      INTEGER, DIMENSION(:), ALLOCATABLE :: OPTIONS_METIS
      REAL, DIMENSION(:), ALLOCATABLE :: COLSCA_TEMP
      INTEGER THRESH, IVersion
      LOGICAL AGG6
      INTEGER MINSYM
      PARAMETER (MINSYM=50)
      INTEGER K79REF
      PARAMETER(K79REF=12000000)
      INTEGER PIV(N)
      INTEGER MTRANS, COMPRESS,NCMP,IERROR,J,JPERM,NCST
      INTEGER TOTEL
      LOGICAL IDENT,SPLITROOT
      EXTERNAL MUMPS_197, SMUMPS_198, 
     *     SMUMPS_199, SMUMPS_351,
     *     SMUMPS_557, SMUMPS_201
#if defined(OLDDFS)
      EXTERNAL SMUMPS_200
#endif
      EXTERNAL SMUMPS_623
      EXTERNAL SMUMPS_547, SMUMPS_550,
     *     SMUMPS_556
      ALLOCATE( IW ( LIW ), stat = IERR )
      IF ( IERR .GT. 0 ) THEN
         INFO( 1 ) = -7
         INFO( 2 ) = LIW
         RETURN
      ENDIF
      LLIW = LIW - 2*N - 1
      L1 = LLIW + 1
      L2 = L1 + N
      LP    = ICNTL(1)
      MP    = ICNTL(3)
      PROK  = (MP.GT.0)
      LDIAG = ICNTL(4)
      COMPRESS_SCHUR = .FALSE.
      IF (KEEP(1).LT.0) KEEP(1) = 0
      NEMIN = KEEP(1)
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
         WRITE (MP,99999) N, NZ, LIW, INFO(1)
         K = min0(10,NZ)
         IF (LDIAG.EQ.4) K = NZ
         IF (K.GT.0) WRITE (MP,99998) (IRN(I),ICN(I),I=1,K)
         K = min0(10,N)
         IF (LDIAG.EQ.4) K = N
         IF (IORD.EQ.1 .AND. K.GT.0) THEN
            WRITE (MP,99997) (IKEEP(I,1),I=1,K)
         ENDIF
      ENDIF
      NCMP    = N   
      IF (KEEP(60).NE.0) THEN
         IF ((SIZE_SCHUR.LE.0 ).OR.
     &        (SIZE_SCHUR.GE.N) ) GOTO 90
      ENDIF
#if defined(metis)
      IF  ( ( KEEP(60).NE.0).AND.(SIZE_SCHUR.GT.0)
     &     .AND.
     &     ((IORD.EQ.7).OR.(IORD.EQ.5))
     &     )THEN
         COMPRESS_SCHUR=.TRUE.
         NCMP          = N-SIZE_SCHUR
         CALL SMUMPS_623(N,NCMP,NZ,IRN, ICN, IW(1), LLIW, 
     *        IW(L2), PTRAR(1,2),
     *        PTRAR, IW(L1), IWFR, KEEP(113), KEEP(114),
     *        INFO(1), INFO(2), ICNTL, symmetry, 
     *        SYM, MedDens, NBQD, AvgDens, 
     *        LISTVAR_SCHUR, SIZE_SCHUR, 
     *        FRERE,FILS)
         IORD = 5
         KEEP(95) = 1
         NBQD     = 0           
      ELSE
#endif
         CALL SMUMPS_351(N,NZ,IRN, ICN, IW(1), LLIW, 
     *        IW(L2), PTRAR(1,2),
     *        PTRAR, IW(L1), IWFR, KEEP(113), KEEP(114),
     *        INFO(1), INFO(2), ICNTL, symmetry, 
     *        SYM, MedDens, NBQD, AvgDens)
#if defined(metis)
      ENDIF
#endif
      INFO(8) = symmetry
      IF(NBQD .GT. 0) THEN
         IF( KEEP(50) .EQ. 2 .AND. ICNTL(12) .LE. 1 ) THEN
            IF(KEEP(95) .NE. 1) THEN
               IF ( PROK ) 
     *              WRITE( MP,*) 
     *              'Compressed/constrained ordering set OFF'
               KEEP(95) = 1   
            ENDIF
         ENDIF
      ENDIF
      IF ( (KEEP(60).NE.0) .AND. (IORD.GT.1) .AND.
     &     .NOT. COMPRESS_SCHUR ) THEN
         IORD = 0               
      ENDIF 
      CALL SMUMPS_701( N, SYM, NSLAVES, IORD, 
     &     symmetry, MedDens, NBQD, AvgDens,
     &     PROK, MP )
      IF(SYM .EQ. 2) THEN
         IF(KEEP(95) .EQ. 3 .AND. IORD .NE. 2) THEN
            IF (PROK) WRITE(MP,*)
     *           'WARNING: ANAL_F constrained ordering not available ', 
     *           'with selected ordering'
            KEEP(95) = 2
         ENDIF
         IF(KEEP(95) .EQ. 2 .AND. IORD .EQ. 0) THEN 
            IF (PROK) WRITE(MP,*)
     *           'WARNING: ANAL_F AMD not available with ', 
     *           'compressed ordering -> move to QAMD'
            IORD = 6
         ENDIF
      ELSE
         KEEP(95) = 1
      ENDIF
      MTRANS = KEEP(23)
      COMPRESS = KEEP(95) - 1
      IF(COMPRESS .GT. 0 .AND. KEEP(52) .EQ. -2) THEN
         IF(id%CNTL(4) .GE. 0.0E0) THEN
            IF (KEEP(1).LE.8) THEN
               KEEP(1) = 16   
            ELSE
               KEEP(1) = 2*KEEP(1)
            ENDIF
            IF (PROK) 
     *           WRITE(MP,*) 'Setting static pivoting ON, COMPRESS =',
     *           COMPRESS
         ENDIF
      ENDIF
      IF(MTRANS .GT. 0 .AND. KEEP(50) .EQ. 2) THEN
         KEEP(23) = 0
      ENDIF
      IF(COMPRESS .EQ. 2) THEN
         IF (IORD.NE.2) THEN
            WRITE(*,*) "IORD not compatible with COMPRESS:",
     *           IORD, COMPRESS
            CALL MUMPS_ABORT()
         ENDIF
         CALL  SMUMPS_556(
     *        N,PIV,FRERE,FILS,NFSIZ,IKEEP,
     *        NCST,KEEP,KEEP8,id)
      ENDIF
      IF ( IORD .NE. 1 ) THEN
         IF(COMPRESS .GE. 1) THEN
            CALL SMUMPS_547(
     *           N,NZ, IRN, ICN, PIV,
     *           NCMP, IW(1), LLIW, IW(L2),PTRAR(1,2),PTRAR, 
     *           IW(L1), FILS, IWFR,
     *           IERROR, KEEP,KEEP8, ICNTL)
            symmetry = 100
         ENDIF
         IF ( (symmetry.LT.MINSYM).AND.(SYM.EQ.0) ) THEN
            IF(KEEP(23) .EQ. 7 ) THEN
               KEEP(23) = -5
               DEALLOCATE (IW)
               RETURN
            ELSE IF(KEEP(23) .EQ. -9876543) THEN
               IDENT = .TRUE.
               KEEP(23) = 5
               IF (PROK) WRITE(MP,'(A)')
     *              ' ... Apply column permutation (already computed)'
               DO J=1,N
                  JPERM = PIV(J)
                  FILS(JPERM) = J 
                  IF (JPERM.NE.J) IDENT = .FALSE.
               ENDDO
               IF (.NOT.IDENT) THEN
                  DO K=1,NZ
                     J = ICN(K)
                     IF ((J.LE.0).OR.(J.GT.N)) CYCLE
                     ICN(K) = FILS(J)
                  ENDDO
                  ALLOCATE(COLSCA_TEMP(N), stat=IERR)
                  IF ( IERR > 0 ) THEN
                     INFO( 1 ) = -7
                     INFO( 2 ) = LIW
                     RETURN
                  ENDIF
                  DO J = 1, N
                     COLSCA_TEMP(J)=id%COLSCA(J)
                  ENDDO
                  DO J=1, N
                     id%COLSCA(FILS(J))=COLSCA_TEMP(J)
                  ENDDO
                  DEALLOCATE(COLSCA_TEMP)
                  IF (MP.GT.0 .AND. ICNTL(4).GE.2)
     *                 WRITE(MP,'(/A)')
     *                 ' WARNING input matrix data modified'
                  CALL SMUMPS_351
     *                 (N,NZ,IRN, ICN, IW(1), LLIW, IW(L2), PTRAR(1,2), 
     *                 PTRAR, IW(L1), IWFR, KEEP(113), KEEP(114), 
     *                 INFO(1), INFO(2), ICNTL, symmetry, SYM,
     *                 MedDens, NBQD, AvgDens)
                  INFO(8) = symmetry
                  NCMP = N
               ELSE
                  KEEP(23) = 0
               ENDIF
            ENDIF
         ELSE IF (KEEP(23) .EQ. 7 .OR. KEEP(23) .EQ. -9876543 ) THEN
            IF (PROK) WRITE(MP,'(A)')
     *           ' ... No column permutation'
            KEEP(23) = 0
         ENDIF
      ENDIF                     
      IF (IORD.NE.1 .AND. IORD.NE.5) THEN
         IF (PROK) THEN
            IF (IORD.EQ.2) THEN
               WRITE(MP,'(A)') ' Ordering based on AMF '
#if defined(scotch)
            ELSE IF (IORD.EQ.3) THEN
               WRITE(MP,'(A)') ' Ordering based on SCOTCH '
#endif
#if defined(pord)
            ELSE IF (IORD.EQ.4) THEN
               WRITE(MP,'(A)') ' Ordering based on PORD '
#endif
            ELSE IF (IORD.EQ.6) THEN
               WRITE(MP,'(A)') ' Ordering based on QAMD '
            ELSE
               WRITE(MP,'(A)') ' Ordering based on AMD '
            ENDIF
         ENDIF
         IF ( KEEP(60) .NE. 0 ) THEN
            CALL MUMPS_162(N, LLIW, IW(L2), IWFR, PTRAR(1,2), IW(1), 
     *           IW(L1), IKEEP, 
     *           IKEEP(1,2), NCMPA, FILS, IKEEP(1,3), PTRAR, PTRAR(1,3),
     *           LISTVAR_SCHUR, SIZE_SCHUR)
            IF (KEEP(60)==1) THEN
               KEEP(20) = LISTVAR_SCHUR(1)
            ELSE
               KEEP(38) = LISTVAR_SCHUR(1)
            ENDIF
         ELSE
            IF ( .FALSE. ) THEN
#if defined(pord)
            ELSEIF (IORD .EQ. 4) THEN
               IF(COMPRESS .EQ. 1) THEN
                  DO I=L1,L1-1+KEEP(93)/2
                     IW(I) = 2
                  ENDDO
                  DO I=L1+KEEP(93)/2,L1+NCMP-1
                     IW(I) = 1
                  ENDDO
                  CALL MUMPS_PORDF_WND(NCMP, IWFR-1, IW(L2), IW, 
     *                 IW(L1), NCMPA, N)
                  CALL SMUMPS_549(NCMP,IW(L2),IKEEP(1,1),
     *                 FRERE,PTRAR(1,1))
                  DO I=1,NCMP
                     IKEEP(IKEEP(I,1),2)=I
                  ENDDO
               ELSE
                  CALL MUMPS_PORDF(NCMP, IWFR-1, IW(L2), IW(1), 
     *                 IW(L1), NCMPA)
               ENDIF
               IF ( NCMPA .NE. 0 ) THEN
                  write(6,*) ' Out PORD, NCMPA=', NCMPA
                  INFO( 1 ) = -9999
                  INFO( 2 ) = 4
                  RETURN
               ENDIF
#endif    
#if defined(scotch)
            ELSEIF (IORD .EQ. 3) THEN
               CALL MUMPS_SCOTCH(NCMP, LLIW, IW(L2), IWFR,
     $              PTRAR(1,2), IW(1), IW(L1), IKEEP, 
     *              IKEEP(1,2), NCMPA)
               IF ( NCMPA .NE. 0 ) THEN
                  write(6,*) ' Out SCTOCH, NCMPA=', NCMPA
                  INFO( 1 ) = -9999
                  INFO( 2 ) = 3
                  RETURN
               ENDIF
#endif
            ELSEIF (IORD .EQ. 2) THEN
               NBBUCK = 2*N
               ALLOCATE( HEAD ( 0: NBBUCK + 1), stat = IERR )
               IF ( IERR .GT. 0 ) THEN
                  INFO( 1 ) = -7
                  INFO( 2 ) = NBBUCK+2
                  RETURN
               ENDIF
               IF(COMPRESS .GE. 1) THEN
                  DO I=L1,L1-1+KEEP(93)/2
                     IW(I) = 2
                  ENDDO
                  DO I=L1+KEEP(93)/2,L1+NCMP-1
                     IW(I) = 1
                  ENDDO
               ELSE
                  IW(L1) = -1
               ENDIF
               IF(COMPRESS .LE. 1) THEN
                  CALL MUMPS_337(NCMP, NBBUCK, LLIW, IW(L2),
     *                 IWFR, PTRAR(1,2),
     *                 IW(1), IW(L1), IKEEP, IKEEP(1,2), NCMPA, FILS,
     *                 IKEEP(1,3), PTRAR, PTRAR(1,3), HEAD)
               ELSE
                  IF(PROK) WRITE(MP,'(A)') 
     *                 ' Constrained Ordering based on AMF'
                  CALL MUMPS_560(NCMP, NBBUCK, LLIW, IW(L2),
     *                 IWFR, PTRAR(1,2), 
     *                 IW(1), IW(L1), IKEEP, IKEEP(1,2), NCMPA, FILS, 
     *                 IKEEP(1,3), PTRAR, PTRAR(1,3), HEAD,
     *                 NFSIZ, FRERE)
               ENDIF
               DEALLOCATE(HEAD)
            ELSEIF (IORD .EQ. 6) THEN
               ALLOCATE( HEAD ( N ), stat = IERR )
               IF ( IERR .GT. 0 ) THEN
                  INFO( 1 ) = -7
                  INFO( 2 ) = N
                  RETURN
               ENDIF
               THRESH = 1
               IVersion = 2
               IF(COMPRESS .EQ. 1) THEN
                  DO I=L1,L1-1+KEEP(93)/2
                     IW(I) = 2
                  ENDDO
                  DO I=L1+KEEP(93)/2,L1+NCMP-1
                     IW(I) = 1
                  ENDDO
                  TOTEL = KEEP(93)+KEEP(94)
               ELSE
                  IW(L1) = -1
                  TOTEL = N
               ENDIF
               CALL MUMPS_421(TOTEL,IVersion, THRESH, HEAD,
     *              NCMP, LLIW, IW(L2), IWFR, PTRAR(1,2), IW(1),
     *              IW(L1), IKEEP, IKEEP(1,2), NCMPA, FILS,
     $              IKEEP(1,3), PTRAR, PTRAR(1,3))
               DEALLOCATE(HEAD)
            ELSE
               CALL MUMPS_197(NCMP, LLIW, IW(L2), IWFR, PTRAR(1,2),
     $              IW(1), IW(L1), IKEEP, IKEEP(1,2), NCMPA, FILS,
     $              IKEEP(1,3), PTRAR, PTRAR(1,3))
            ENDIF
         ENDIF
         IF(COMPRESS .GE. 1) THEN
            CALL SMUMPS_550(N,NCMP,KEEP(94),KEEP(93),
     *           PIV,IKEEP(1,1),IKEEP(1,2))
            COMPRESS = -1
         ENDIF
      ENDIF  
#if defined(metis)
      IF (IORD.EQ.5) THEN
         IF (PROK) THEN
            WRITE(MP,'(A)') ' Ordering based on METIS '
         ENDIF
         NUMFLAG = 1
         OPT_METIS_SIZE = 8
         ALLOCATE( OPTIONS_METIS (OPT_METIS_SIZE ), stat = IERR )
         IF ( IERR .GT. 0 ) THEN
            INFO( 1 ) = -7
            INFO( 2 ) = OPT_METIS_SIZE
            RETURN
         ENDIF
         OPTIONS_METIS(1) = 0
         IF (COMPRESS .EQ. 1) THEN
            DO I=1,KEEP(93)/2
               FILS(I) = 2
            ENDDO
            DO I=KEEP(93)/2+1,NCMP
               FILS(I) = 1
            ENDDO
            CALL METIS_NODEWND(NCMP, IW(L2), IW(1),FILS, 
     *           NUMFLAG, OPTIONS_METIS,
     *           IKEEP(1,2), IKEEP(1,1) )
         ELSE
            CALL METIS_NODEND(NCMP, IW(L2), IW(1), NUMFLAG, 
     *           OPTIONS_METIS,
     *           IKEEP(1,2), IKEEP(1,1) )
         ENDIF
         DEALLOCATE (OPTIONS_METIS)
         IF ( COMPRESS_SCHUR ) THEN
            CALL SMUMPS_622(
     *           N, NCMP, IKEEP(1,1),IKEEP(1,2),
     *           LISTVAR_SCHUR, SIZE_SCHUR, FILS)
            COMPRESS = -1       
         ENDIF
         IF (COMPRESS .EQ. 1) THEN
            CALL SMUMPS_550(N,NCMP,KEEP(94),
     *           KEEP(93),PIV,IKEEP(1,1),IKEEP(1,2))
            COMPRESS = -1       
         ENDIF
      ENDIF                     
#endif
      IF (PROK) THEN
         IF (IORD.EQ.1) THEN
            WRITE(MP,'(A)') ' Ordering given is used'
         ENDIF
      ENDIF
      IF ((IORD.EQ.1) 
     $     ) THEN
         DO K=1,N
            PTRAR(K,1) = 0
         ENDDO
         DO K=1,N
            IF ((IKEEP(K,1).LE.0).OR.(IKEEP(K,1).GT.N)) 
     *           GO TO 40
            IF (PTRAR(IKEEP(K,1),1).EQ.1) THEN
               GOTO 40
            ELSE
               PTRAR(IKEEP(K,1),1) = 1
            ENDIF
         ENDDO
      ENDIF
      IF (IORD.EQ.1 .OR. IORD.EQ.5 .OR. COMPRESS.EQ.-1) THEN
         IF (KEEP(106)==1) THEN
            IF ( COMPRESS .EQ. -1 ) THEN
               CALL SMUMPS_351(N,NZ,IRN, ICN, IW(1), LLIW,
     *              IW(L2), PTRAR(1,2),
     *              PTRAR, IW(L1), IWFR, KEEP(113), KEEP(114),
     *              INFO(1), INFO(2), ICNTL, symmetry, SYM,
     *              MedDens, NBQD, AvgDens)
               INFO(8) = symmetry
            ENDIF
            COMPRESS = 0
            ALLOCATE( HEAD ( 2*N ), stat = IERR )
            IF ( IERR .GT. 0 ) THEN
               INFO( 1 ) = -7
               INFO( 2 ) = 2*N
               RETURN
            ENDIF
            THRESH = -1
            IF (KEEP(60) == 0) THEN
               ITEMP = 0 
            ELSE
               ITEMP = SIZE_SCHUR
               IF (KEEP(60)==1) THEN
                  KEEP(20) = LISTVAR_SCHUR(1)
               ELSE
                  KEEP(38) = LISTVAR_SCHUR(1)
               ENDIF
            ENDIF
            AGG6 = ( NSLAVES == 1 )
            CALL MUMPS_422(THRESH, HEAD,
     *           N, LLIW, IW(L2), IWFR, PTRAR(1,2), IW,
     *           IW(L1), HEAD(N+1),
     *           IKEEP(1,2), NCMPA, FILS, IKEEP(1,3), PTRAR, PTRAR(1,3), 
     *           IKEEP(1,1), LISTVAR_SCHUR, ITEMP, AGG6)
            DEALLOCATE(HEAD)
         ELSE
            CALL SMUMPS_198(N, NZ, IRN, ICN, IKEEP, IW(1),
     *           LLIW, IW(L2),
     *           PTRAR(1,2), IW(L1), IWFR,
     *           INFO(1),INFO(2), KEEP(11), MP)
            IF (KEEP(60) .EQ. 0) THEN
               ITEMP = 0 
               CALL SMUMPS_199(N, IW(L2), IW, LLIW, IWFR, IKEEP,
     *              IKEEP(1,2), IW(L1),
     *              PTRAR, NCMPA, ITEMP)
            ELSE
               CALL SMUMPS_199(N, IW(L2), IW, LLIW, IWFR, IKEEP,
     *              IKEEP(1,2), IW(L1),
     *              PTRAR, NCMPA, SIZE_SCHUR)
               IF (KEEP(60) .EQ. 1) THEN
                  KEEP(20) = LISTVAR_SCHUR(1)
               ELSE
                  KEEP(38) = LISTVAR_SCHUR(1)
               ENDIF
            ENDIF
         ENDIF                  
      ENDIF                     
#if defined(OLDDFS)
      CALL SMUMPS_200
     *     (N, IW(L2), IW(L1), IKEEP(1,1), IKEEP(1,2), IKEEP(1,3),
     *     NFSIZ, INFO(6), FILS, FRERE, PTRAR(1,3), NEMIN, KEEP(60))
#else
      CALL SMUMPS_557
     *     (N, IW(L2), IW(L1), IKEEP(1,1), IKEEP(1,2), IKEEP(1,3),
     *     NFSIZ, PTRAR, INFO(6), FILS, FRERE, 
     *     PTRAR(1,3), NEMIN, PTRAR(1,4), KEEP(60),
     *     KEEP(20),KEEP(38),PTRAR(1,2),KEEP(104),IW(1),KEEP(50))
#endif
      IF (KEEP(60).NE.0)  THEN
         IF (KEEP(60)==1) THEN
            IN = KEEP(20)
         ELSE
            IN = KEEP(38)
         ENDIF
         DO WHILE (IN.GT.0) 
            IN = FILS (IN)
         END DO
         IFSON = -IN
         IF (KEEP(60)==1) THEN
            IN = KEEP(20)
         ELSE
            IN = KEEP(38)
         ENDIF
         DO I=2,SIZE_SCHUR
            FILS(IN) = LISTVAR_SCHUR (I)
            IN       = FILS(IN)
            FRERE (IN) = N+1
         ENDDO
         FILS(IN) = -IFSON
      ENDIF
      CALL SMUMPS_201(IKEEP(1,2),
     *     PTRAR(1,3), INFO(6),
     *     INFO(5), KEEP(2), KEEP(50),
     *     KEEP(101),KEEP(108),KEEP(5),
     *     KEEP(6), KEEP(226))
      IF ( KEEP(53) .NE. 0 ) THEN
         CALL MUMPS_209( N, FRERE, FILS, NFSIZ, KEEP(20) )
      END IF
      IF (  (KEEP(48) == 4 .AND. KEEP(10).GT.0)
     *     .OR.
     *     (KEEP (48)==5 .AND. KEEP(10) .GT. 0 )
     *     .OR.
     *     (KEEP(24).NE.0.AND.KEEP(10).GT.0) ) THEN 
         CALL SMUMPS_510(KEEP(10), KEEP(2),
     *        KEEP(48), KEEP(50), NSLAVES)
      END IF
      IF (KEEP(210).LT.0.OR.KEEP(210).GT.2) KEEP(210)=0
      IF (KEEP(210).EQ.0.AND.KEEP(201).GT.0) KEEP(210)=1 
      IF (KEEP(210).EQ.0.AND.KEEP(201).EQ.0) KEEP(210)=2 
      IF (KEEP(210).EQ.2) KEEP(79)=huge(KEEP(79))
      IF (KEEP(210).EQ.1.AND.KEEP(79).LE.0) THEN
         IF ( huge(KEEP(79)) / K79REF + 1 .GE. NSLAVES ) THEN
            KEEP(79)=huge(KEEP(79))
         ELSE
            KEEP(79)=K79REF * NSLAVES
         ENDIF
      ENDIF
      IF (KEEP(210).EQ.1) THEN
         SPLITROOT = .FALSE. 
         IF ( KEEP(62).GE.1) THEN
            CALL SMUMPS_97(N, FRERE, FILS, NFSIZ,INFO(6),
     *           NSLAVES, KEEP,KEEP8, SPLITROOT,
     *           MP, LDIAG,INFO(1),INFO(2))
            IF (INFO(1).LT.0) RETURN
         ENDIF
      ENDIF
      SPLITROOT = ((ICNTL(13).GT.0 .AND. NSLAVES.GT.ICNTL(13)) .OR.
     &     ICNTL(13).EQ.-1 )
     &     .AND. (KEEP(60).EQ.0)
      IF (SPLITROOT) THEN
         CALL SMUMPS_97(N, FRERE, FILS, NFSIZ,INFO(6),
     *        NSLAVES, KEEP,KEEP8, SPLITROOT,
     *        MP, LDIAG,INFO(1),INFO(2))
         IF (INFO(1).LT.0) RETURN
      ENDIF
C     $        IKEEP(1,3), IKEEP(1,2), IRN, ICN, PTRAR, KEEP)
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
         K = min0(10,N)
         IF (LDIAG.EQ.4) K = N
         IF (K.GT.0) WRITE (MP,99997) (IKEEP(I,1),I=1,K)
         IF (K.GT.0) WRITE (MP,99991) (IKEEP(I,2),I=1,K)
         IF (K.GT.0) WRITE (MP,99990) (IKEEP(I,3),I=1,K)
         IF (K.GT.0) WRITE (MP,99986) (PTRAR(I,1),I=1,K)
         IF (K.GT.0) WRITE (MP,99985) (PTRAR(I,2),I=1,K)
         IF (K.GT.0) WRITE (MP,99984) (PTRAR(I,3),I=1,K)
         IF (K.GT.0) WRITE (MP,99987) (NFSIZ(I),I=1,K)
         IF (K.GT.0) WRITE (MP,99989) (FILS(I),I=1,K)
         IF (K.GT.0) WRITE (MP,99988) (FRERE(I),I=1,K)
      ENDIF
      GO TO 90
 40   INFO(1) = -4
      INFO(2) = K
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) WRITE (LP,99996) INFO(1)
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) WRITE (LP,99982) INFO(2)
      GOTO 90
 90   CONTINUE 
      DEALLOCATE(IW)
      RETURN
99999 FORMAT (/'Entering analysis phase with ...'/
     *     '                N         NZ         LIW       INFO(1)'/,
     *     9X, I8, I11, I12, I14)
99998 FORMAT ('Matrix entries:    IRN()   ICN()'/
     *     (I12, I7, I12, I7, I12, I7))
99997 FORMAT ('IKEEP(.,1)=', 10I6/(12X, 10I6))
99996 FORMAT (/'** Error return ** from Analysis *  INFO(1)=', I3)
99991 FORMAT ('IKEEP(.,2)=', 10I6/(12X, 10I6))
99990 FORMAT ('IKEEP(.,3)=', 10I6/(12X, 10I6))
99989 FORMAT ('FILS (.)  =', 10I6/(12X, 10I6))
99988 FORMAT ('FRERE(.)  =', 10I6/(12X, 10I6))
99987 FORMAT ('NFSIZ(.)  =', 10I6/(12X, 10I6))
99986 FORMAT ('PTRAR(.,1)=', 10I6/(12X, 10I6))
99985 FORMAT ('PTRAR(.,2)=', 10I6/(12X, 10I6))
99984 FORMAT ('PTRAR(.,3)=', 10I6/(12X, 10I6))
99982 FORMAT ('Error in permutation array KEEP   INFO(2)=', I3)
      END SUBROUTINE SMUMPS_195
#if defined(NOT_USED)
      SUBROUTINE SMUMPS_196(N,NZ, IRN, ICN, IW, LW, IPE, LEN,
     *  IQ, FLAG, IWFR,
     * NRORM, NIORM, IFLAG,IERROR, ICNTL, 
     * symmetry, SYM)
      INTEGER N,NZ,LW,IFLAG,IERROR,NRORM,NIORM,IWFR
      INTEGER symmetry, SYM
      INTEGER ICNTL(40)
      INTEGER  IRN(NZ), ICN(NZ) 
      INTEGER LEN(N), IQ(N)
      INTEGER IPE(N)
      INTEGER FLAG(N), IW(LW)
      INTEGER MP, MPG
      INTEGER I,K,J,N1,LAST,NDUP,K1,K2,L
      INTEGER NBERR
      INTEGER NZOFFA, NDIAGA
      REAL RSYM
      INTRINSIC nint
      MP = ICNTL(2)
      MPG = ICNTL(3)
      NIORM  = 3*N
      NDIAGA = 0
      IERROR = 0
      DO 10 I=1,N
        IPE(I) = 0
   10 CONTINUE
      DO 50 K=1,NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     *                          .OR.(J.LT.1)) THEN
           IERROR = IERROR + 1
        ELSE
          IF (I.NE.J) THEN
           IPE(I) = IPE(I) + 1
           IPE(J) = IPE(J) + 1
           NIORM  = NIORM + 1
          ELSE
           NDIAGA = NDIAGA + 1
          ENDIF
        ENDIF
   50 CONTINUE
      NZOFFA  = NIORM - 3*N
      IF (IERROR.GE.1) THEN
         NBERR  = 0
         IF (mod(IFLAG,2).EQ.0) IFLAG  = IFLAG+1
         IF ((MP.GT.0).AND.(ICNTL(4).GE.2))  THEN 
          WRITE (MP,99999) 
          DO 70 K=1,NZ
           I = IRN(K)
           J = ICN(K)
           IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     *                            .OR.(J.LT.1)) THEN
            NBERR = NBERR + 1
            IF (NBERR.LE.10)  THEN
               IF (mod(K,10).GT.3 .OR. mod(K,10).EQ.0 .OR.
     *             (10.LE.K .AND. K.LE.20)) THEN
                 WRITE (MP,'(I8,A,I8,A,I8,A)')
     *             K,'th entry (in row',I,' and column',J,') ignored'
               ELSE
                 IF (mod(K,10).EQ.1) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'st entry (in row',I,' and column',J,') ignored'
                 IF (mod(K,10).EQ.2) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'nd entry (in row',I,' and column',J,') ignored'
                 IF (mod(K,10).EQ.3) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'rd entry (in row',I,' and column',J,') ignored'
               ENDIF
            ELSE
               GO TO 100
            ENDIF
           ENDIF
   70     CONTINUE
         ENDIF
      ENDIF
  100 NRORM = NIORM - 2*N
      IQ(1) = 1
      N1 = N - 1
      IF (N1.GT.0) THEN
        DO 110 I=1,N1
            IQ(I+1) = IPE(I) + IQ(I) 
  110   CONTINUE
      ENDIF
      LAST = max(IPE(N)+IQ(N)-1,IQ(N))
      DO 115 I = 1,N
         FLAG(I) = 0
         IPE(I)  = IQ(I)
  115 CONTINUE
      DO 130 K=1,LAST
        IW(K) = 0
  130 CONTINUE
      IWFR = LAST + 1
      DO 200 K=1,NZ
         I = IRN(K)
         J = ICN(K)
         IF (I.NE.J) THEN
          IF (I.LT.J) THEN
            IF ((I.GE.1).AND.(J.LE.N)) THEN
             IW(IQ(I)) = -J
             IQ(I)     = IQ(I) + 1 
            ENDIF
          ELSE
            IF ((J.GE.1).AND.(I.LE.N)) THEN
             IW(IQ(J)) = -I
             IQ(J)     = IQ(J) + 1
            ENDIF
          ENDIF
         ENDIF
  200 CONTINUE
      NDUP = 0
      DO 260 I=1,N
        K1 = IPE(I) 
        K2 = IQ(I) -1
        IF (K1.GT.K2) THEN
         LEN(I) = 0
         IPE(I) = 0
         IQ(I)  = 0
        ELSE
         DO 240 K=K1,K2
           J     = -IW(K)
           IF (J.LE.0) GO TO 250
           L     = IQ(J) 
           IQ(J) = L + 1
           IF (FLAG(J).EQ.I) THEN
            NDUP = NDUP + 1
            IW(L) = 0
            IW(K) = 0
           ELSE
            IW(L)   = I
            IW(K)   = J
            FLAG(J) = I
           ENDIF
  240    CONTINUE
  250    IQ(I) = IQ(I) - IPE(I)
         IF (NDUP.EQ.0) LEN(I) = IQ(I)
        ENDIF
  260 CONTINUE
      IF (NDUP.NE.0) THEN
       IWFR = 1
       DO 280 I=1,N
         K1 = IPE(I) 
         IF (K1.EQ.0) GO TO 280
         K2 = K1 + IQ(I) - 1
         L = IWFR
         IPE(I) = IWFR
         DO 270 K=K1,K2
           IF (IW(K).NE.0) THEN
            IW(IWFR) = IW(K)
            IWFR     = IWFR + 1
           ENDIF
  270    CONTINUE
         LEN(I) = IWFR - L 
  280  CONTINUE
      ENDIF
      IF (SYM.EQ.0) THEN
      RSYM =  real(NDIAGA+2*NZOFFA - (IWFR-1))/
     *            real(NZOFFA+NDIAGA) 
      symmetry = nint (100.0*RSYM)
      ELSE
       symmetry = 100
      ENDIF
      RETURN
99999 FORMAT (/'*** Warning message from analysis routine ***')
      END SUBROUTINE SMUMPS_196
#endif
      SUBROUTINE SMUMPS_199(N,IPE,IW, LW, IWFR, IPS, IPV, NV, FLAG,
     *                  NCMPA, SIZE_SCHUR)
      INTEGER N,LW,IWFR,NCMPA,SIZE_SCHUR
      INTEGER FLAG(N)
      INTEGER IPS(N), IPV(N)
      INTEGER IW(LW), NV(N), IPE(N)
      INTEGER I,J,ML,MS,ME,IP,MINJS,IE,KDUMMY,JP
      INTEGER LN,JP1,JS,LWFR,JP2,JE
      DO 10 I=1,N
        FLAG(I) = 0
        NV(I) = 0
        J = IPS(I)
        IPV(J) = I
   10 CONTINUE
      NCMPA = 0
      DO 100 ML=1,N-SIZE_SCHUR 
        MS = IPV(ML)
        ME = MS
        FLAG(MS) = ME
        IP = IWFR
        MINJS = N
        IE = ME
        DO 70 KDUMMY=1,N
          JP = IPE(IE)
          LN = 0
          IF (JP.LE.0) GO TO 60
          LN = IW(JP)
          DO 50 JP1=1,LN
            JP = JP + 1
            JS = IW(JP)
            IF (FLAG(JS).EQ.ME) GO TO 50
            FLAG(JS) = ME
            IF (IWFR.LT.LW) GO TO 40
            IPE(IE) = JP
            IW(JP) = LN - JP1
            CALL SMUMPS_194(N, IPE, IW, IP-1, LWFR,NCMPA)
            JP2 = IWFR - 1
            IWFR = LWFR
            IF (IP.GT.JP2) GO TO 30
            DO 20 JP=IP,JP2
              IW(IWFR) = IW(JP)
              IWFR = IWFR + 1
   20       CONTINUE
   30       IP = LWFR
            JP = IPE(IE)
   40       IW(IWFR) = JS
            MINJS = min0(MINJS,IPS(JS)+0)
            IWFR = IWFR + 1
   50     CONTINUE
   60     IPE(IE) = -ME
          JE = NV(IE)
          NV(IE) = LN + 1
          IE = JE
          IF (IE.EQ.0) GO TO 80
   70   CONTINUE
   80   IF (IWFR.GT.IP) GO TO 90
        IPE(ME) = 0
        NV(ME) = 1
        GO TO 100
   90   MINJS = IPV(MINJS)
        NV(ME) = NV(MINJS)
        NV(MINJS) = ME
        IW(IWFR) = IW(IP)
        IW(IP) = IWFR - IP
        IPE(ME) = IP
        IWFR = IWFR + 1
  100 CONTINUE
      IF (SIZE_SCHUR == 0) RETURN
      DO ML = N-SIZE_SCHUR+1,N
        ME = IPV(ML)
        IE = ME
        DO KDUMMY=1,N
          JP = IPE(IE)
          LN = 0
          IF (JP.LE.0) GO TO 160
          LN = IW(JP)
  160     IPE(IE) = -IPV(N-SIZE_SCHUR+1)
          JE = NV(IE)
          NV(IE) = LN + 1
          IE = JE
          IF (IE.EQ.0) GO TO 190
        ENDDO
  190   NV(ME) = 0
        IPE(ME) = -IPV(N-SIZE_SCHUR+1)
      ENDDO
      ME = IPV(N-SIZE_SCHUR+1)
      IPE(ME) = 0
      NV(ME) = SIZE_SCHUR
      RETURN
      END SUBROUTINE SMUMPS_199
      SUBROUTINE SMUMPS_198(N, NZ, IRN, ICN, PERM,
     * IW, LW, IPE, IQ, FLAG,
     * IWFR, IFLAG, IERROR, IOVFLO, MP)
      INTEGER N,NZ,LW,IWFR,IFLAG,IERROR
      INTEGER PERM(N)
      INTEGER IQ(N)
      INTEGER IRN(NZ), ICN(NZ) 
      INTEGER IPE(N), IW(LW), FLAG(N)
      INTEGER MP
      INTEGER IOVFLO
      INTEGER I,J,K,LBIG,L,ID,IN,LEN,JDUMMY,K1,K2
      IERROR = 0
      DO 10 I=1,N
        IQ(I) = 0
   10 CONTINUE
      DO 80 K=1,NZ
        I = IRN(K)
        J = ICN(K)
        IW(K) = -I
        IF (I.EQ.J) GOTO 40
        IF (I.GT.J) GOTO 30
        IF (I.GE.1 .AND. J.LE.N) GO TO 60
        GO TO 50
   30   IF (J.GE.1 .AND. I.LE.N) GO TO 60
        GO TO 50
   40   IW(K) = 0
        IF (I.GE.1 .AND. I.LE.N) GO TO 80
   50   IERROR = IERROR + 1
        IW(K) = 0
        IF (IERROR.LE.1 .AND. MP.GT.0) WRITE (MP,99999) 
        IF (IERROR.LE.10 .AND. MP.GT.0) WRITE (MP,99998) K, I, J
        GO TO 80
   60   IF (PERM(J).GT.PERM(I)) GO TO 70
        IQ(J) = IQ(J) + 1
        GO TO 80
   70   IQ(I) = IQ(I) + 1
   80 CONTINUE
      IF (IERROR.GE.1) THEN
        IF (mod(IFLAG,2) .EQ. 0) IFLAG = IFLAG+1
      ENDIF
      IWFR = 1
      LBIG = 0
      DO 100 I=1,N
        L = IQ(I)
        LBIG = MAX0(L,LBIG)
        IWFR = IWFR + L
        IPE(I) = IWFR - 1
  100 CONTINUE
      DO 140 K=1,NZ
        I = -IW(K)
        IF (I.LE.0) GO TO 140
        L = K
        IW(K) = 0
        DO 130 ID=1,NZ
          J = ICN(L)
          IF (PERM(I).LT.PERM(J)) GO TO 110
          L = IPE(J)
          IPE(J) = L - 1
          IN = IW(L)
          IW(L) = I
          GO TO 120
  110     L = IPE(I)
          IPE(I) = L - 1
          IN = IW(L)
          IW(L) = J
  120     I = -IN
          IF (I.LE.0) GO TO 140
  130   CONTINUE
  140 CONTINUE
      K = IWFR - 1
      L = K + N
      IWFR = L + 1
      DO 170 I=1,N
        FLAG(I) = 0
        J = N + 1 - I
        LEN = IQ(J)
        IF (LEN.LE.0) GO TO 160
        DO 150 JDUMMY=1,LEN
          IW(L) = IW(K)
          K = K - 1
          L = L - 1
  150   CONTINUE
  160   IPE(J) = L
        L = L - 1
  170 CONTINUE
      IF (LBIG.GE.IOVFLO) GO TO 190
      DO 180 I=1,N
        K = IPE(I)
        IW(K) = IQ(I)
        IF (IQ(I).EQ.0) IPE(I) = 0
  180 CONTINUE
      GO TO 230
  190 IWFR = 1
      DO 220 I=1,N
        K1 = IPE(I) + 1
        K2 = IPE(I) + IQ(I)
        IF (K1.LE.K2) GO TO 200
        IPE(I) = 0
        GO TO 220
  200   IPE(I) = IWFR
        IWFR = IWFR + 1
        DO 210 K=K1,K2
          J = IW(K)
          IF (FLAG(J).EQ.I) GO TO 210
          IW(IWFR) = J
          IWFR = IWFR + 1
          FLAG(J) = I
  210   CONTINUE
        K = IPE(I)
        IW(K) = IWFR - K - 1
  220 CONTINUE
  230 RETURN
99999 FORMAT (' *** WARNING MESSAGE FROM SMUMPS_198 ***' )
99998 FORMAT (I6, ' NON-ZERO (IN ROW, I6, 11H AND COLUMN ', I6,
     * ') IGNORED')
      END SUBROUTINE SMUMPS_198
      SUBROUTINE SMUMPS_194(N, IPE, IW, LW, IWFR,NCMPA)
      INTEGER N,LW,IWFR,NCMPA
      INTEGER IPE(N)
      INTEGER   IW(LW)
      INTEGER I,K1,LWFR,IR,K,K2
      NCMPA = NCMPA + 1
      DO 10 I=1,N
        K1 = IPE(I)
        IF (K1.LE.0) GO TO 10
        IPE(I) = IW(K1)
        IW(K1) = -I
   10 CONTINUE
      IWFR = 1
      LWFR = IWFR
      DO 60 IR=1,N
        IF (LWFR.GT.LW) GO TO 70
        DO 20 K=LWFR,LW
          IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
        GO TO 70
   30   I = -IW(K)
        IW(IWFR) = IPE(I)
        IPE(I) = IWFR
        K1 = K + 1
        K2 = K + IW(IWFR)
        IWFR = IWFR + 1
        IF (K1.GT.K2) GO TO 50
        DO 40 K=K1,K2
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN
      END SUBROUTINE SMUMPS_194
#if defined(OLDDFS)
      SUBROUTINE SMUMPS_200(N, IPE, NV, IPS, NE, NA, NFSIZ, 
     *                  NSTEPS,
     *                  FILS, FRERE,NDD,NEMIN, KEEP60)
      INTEGER N,NSTEPS
      INTEGER NDD(N)
      INTEGER FILS(N), FRERE(N)
      INTEGER IPS(N), NE(N), NA(N), NFSIZ(N)
      INTEGER IPE(N), NV(N)
      INTEGER NEMIN, KEEP60
      INTEGER I,IF,IS,NR,NR1,INS,INL,INB,INF,INFS,INSW
      INTEGER K,L,ISON,IN,INP,IFSON,INC,INO
      INTEGER INOS,IB,IL
      DO 10 I=1,N
        IPS(I) = 0
        NE(I) = 0
   10 CONTINUE
      DO 20 I=1,N
        IF (NV(I).GT.0) GO TO 20
        IF = -IPE(I)
        IS = -IPS(IF)
        IF (IS.GT.0) IPE(I) = IS
        IPS(IF) = -I
   20 CONTINUE
      NR = N + 1
      DO 50 I=1,N
        IF (NV(I).LE.0) GO TO 50
        IF = -IPE(I)
        IF (IF.NE.0) THEN
         IS = -IPS(IF)
         IF (IS.GT.0) IPE(I) = IS
         IPS(IF) = -I
        ELSE
         NR = NR - 1
         NE(NR) = I
        ENDIF
   50 CONTINUE
      DO 999 I=1,N
       FILS(I) = IPS(I)
 999  CONTINUE
      NR1 = NR
      INS = 0
 1000 IF (NR1.GT.N) GO TO 1151
      INS = NE(NR1)
      NR1 = NR1 + 1
 1070 INL = FILS(INS)
      IF (INL.LT.0) THEN
       INS = -INL
       GO TO 1070
      ENDIF
 1080 IF (IPE(INS).LT.0) THEN
       INS       = -IPE(INS)
       FILS(INS) = 0
       GO TO 1080
      ENDIF
      IF (IPE(INS).EQ.0) THEN
       INS = 0
       GO TO 1000
      ENDIF
      INB = IPE(INS)
      IF (NV(INB).EQ.0) THEN
       INS = INB
       GO TO 1070
      ENDIF
      IF (NV(INB).GE.NV(INS)) THEN
       INS = INB
       GO TO 1070
      ENDIF
      INF = INB
 1090 INF = IPE(INF)
      IF (INF.GT.0) GO TO 1090
      INF  = -INF
      INFS = -FILS(INF)
      IF (INFS.EQ.INS) THEN
       FILS(INF) = -INB
       IPS(INF)  = -INB
       IPE(INS)  = IPE(INB)
       IPE(INB)  = INS
       INS       = INB
       GO TO 1070
      ENDIF
      INSW = INFS
 1100 INFS = IPE(INSW)
      IF (INFS.NE.INS) THEN
       INSW = INFS
       GO TO 1100
      ENDIF
      IPE(INS) = IPE(INB)
      IPE(INB) = INS
      IPE(INSW)= INB
      INS      =INB
      GO TO 1070
 1151 CONTINUE
      DO 51 I=1,N
       FRERE(I) = IPE(I)
       FILS(I) = IPS(I)
 51   CONTINUE
      IS = 1
      I  = 0
      IL = 0
      DO 160 K=1,N
        IF (I.GT.0) GO TO 60
        I = NE(NR)
        NE(NR) = 0
        NR = NR + 1
        IL = N
        NA(N) = 0
   60   DO 70 L=1,N
          IF (IPS(I).GE.0) GO TO 80
          ISON = -IPS(I)
          IPS(I) = 0
          I = ISON
          IL = IL - 1
          NA(IL) = 0
   70   CONTINUE
   80   IPS(I) = K
        NE(IS) = NE(IS) + 1
        IF (NV(I).GT.0) GO TO 89
      IN = I
 81   IN =  FRERE(IN)
      IF (IN.GT.0) GO TO 81
      IF = -IN
      IN = IF
 82   INL = IN
      IN = FILS(IN)
      IF (IN.GT.0) GO TO 82
      IFSON = -IN
      FILS(INL) = I
      IN = I
 83   INP = IN
      IN = FILS(IN)
      IF (IN.GT.0) GO TO 83
      IF (IFSON .EQ. I) GO TO 86
      FILS(INP) = -IFSON
      IN = IFSON
 84   INC =IN
      IN = FRERE(IN)
      IF (IN.NE.I) GO TO 84
      FRERE(INC) = FRERE(I)
      GO TO 120
 86   IF (FRERE(I).LT.0) FILS(INP) = 0
      IF (FRERE(I).GT.0) FILS(INP) = -FRERE(I)
      GO TO 120
   89   IF (IL.LT.N) NA(IL+1) = NA(IL+1) + 1
        NA(IS) = NA(IL)
        NDD(IS) = NV(I)
        NFSIZ(I) = NV(I)
        IF (NA(IS).LT.1) GO TO 110
        IF (   (KEEP60.NE.0).AND.
     &         (NE(IS).EQ.NDD(IS)) ) GOTO 110
        IF (NDD(IS-1)-NE(IS-1).EQ.NDD(IS)) GO TO 100
        IF ((NE(IS-1).GE.NEMIN).AND.
     &         (NE(IS).GE.NEMIN) ) GO TO 110
        IF (2*NE(IS-1)*(NDD(IS)-NDD(IS-1)+NE(IS-1)).GE.
     &    ((NDD(IS)+NE(IS-1))*
     &    (NDD(IS)+NE(IS-1))*NEMIN/100)) GO TO 110
  100   NA(IS-1) = NA(IS-1) + NA(IS) - 1
        NDD(IS-1) = NDD(IS) + NE(IS-1)
        NE(IS-1) = NE(IS) + NE(IS-1)
        NE(IS) = 0
      IN=I
 101  INL = IN
      IN = FILS(IN)
      IF (IN.GT.0) GO TO 101
      IFSON = -IN
      IN = IFSON
 102  INO = IN
      IN =  FRERE(IN)
      IF (IN.GT.0) GO TO 102
      FILS(INL) = INO
      NFSIZ(I) = NDD(IS-1)
      IN = INO
 103  INP = IN
      IN = FILS(IN)
      IF (IN.GT.0) GO TO 103
      INOS = -IN
      IF (IFSON.EQ.INO) GO TO 107
      IN = IFSON
      FILS(INP) = -IFSON
 105  INS = IN
      IN =  FRERE(IN)
      IF (IN.NE.INO) GO TO 105
      IF (INOS.EQ.0) FRERE(INS) = -I
      IF (INOS.NE.0) FRERE(INS) =  INOS
      IF (INOS.EQ.0) GO TO 109
 107  IN = INOS
      IF (IN.EQ.0) GO TO 109
 108  INT = IN
      IN =  FRERE(IN)
      IF (IN.GT.0) GO TO 108
      FRERE(INT) = -I
 109  CONTINUE
        GO TO 120
  110   IS = IS + 1
  120   IB = IPE(I)
        IF (IB.LT.0) GOTO 150
        IF (IB.EQ.0) GOTO 140
        NA(IL) = 0
  140   I = IB
        GO TO 160
  150   I = -IB
        IL = IL + 1
  160 CONTINUE
      NSTEPS = IS - 1
      DO 170 I=1,N
        K = FILS(I)
        IF (K.GT.0) THEN
          FRERE(K)  = N + 1
          NFSIZ(K)  = 0
        ENDIF
 170  CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_200
#else
      SUBROUTINE SMUMPS_557(N, IPE, NV, IPS, NE, NA, NFSIZ,
     *               NODE, NSTEPS,
     *               FILS, FRERE, ND, NEMIN, SUBORD, KEEP60, 
     *               KEEP20, KEEP38, NAMALG,NAMALGMAX,
     *               CUMUL,KEEP50)
      IMPLICIT NONE
      INTEGER N, NSTEPS, KEEP60, KEEP20, KEEP38, KEEP50
      INTEGER ND(N), NFSIZ(N)
      INTEGER IPE(N), FILS(N), FRERE(N), SUBORD(N)
      INTEGER NV(N), IPS(N), NE(N), NA(N), NODE(N)
      INTEGER NEMIN,AMALG_COUNT
      INTEGER NAMALG(N),NAMALGMAX, CUMUL(N)
      DOUBLE PRECISION ACCU, FLOPS_FATHER, FLOPS_SON,
     *                  FLOPS_AVANT, FLOPS_APRES
      INTEGER I,IF,IS,NR,NR1,INS,INL,INB,INF,INFS,INSW
      INTEGER K,L,ISON,IN,IFSON,INO
      INTEGER INOS,IB,IL,INT
      INTEGER IPERM
#if ! defined(NOAMALGTOFATHER)
      INTEGER DADI
      LOGICAL AMALG_TO_father_OK
#endif
      AMALG_COUNT = 0
      DO 10 I=1,N
        CUMUL(I)= 0
        IPS(I)  = 0
        NE(I)   = 0
        NODE(I) = 1
        SUBORD(I) = 0
        NAMALG(I) = 0
   10 CONTINUE
      FRERE(1:N) = IPE(1:N)
      NR = N + 1
      DO 50 I=1,N
        IF = -FRERE(I)
        IF (NV(I).EQ.0) THEN
          IF (SUBORD(IF).NE.0) SUBORD(I) = SUBORD(IF)
          SUBORD(IF) = I
          NODE(IF) = NODE(IF)+1
        ELSE
          IF (IF.NE.0) THEN
            IS = -IPS(IF)
            IF (IS.GT.0) FRERE(I) = IS
            IPS(IF) = -I
          ELSE
            NR = NR - 1
            NE(NR) = I
          ENDIF
        ENDIF
   50 CONTINUE
      DO 999 I=1,N
       FILS(I) = IPS(I)
 999  CONTINUE
      NR1 = NR
      INS = 0
 1000 IF (NR1.GT.N) GO TO 1151
      INS = NE(NR1)
      NR1 = NR1 + 1
 1070 INL = FILS(INS)
      IF (INL.LT.0) THEN
       INS = -INL
       GO TO 1070
      ENDIF
 1080 IF (FRERE(INS).LT.0) THEN
       INS       = -FRERE(INS)
       FILS(INS) = 0
       GO TO 1080
      ENDIF
      IF (FRERE(INS).EQ.0) THEN
       INS = 0
       GO TO 1000
      ENDIF
      INB = FRERE(INS)
#if ! defined(NOAMALGTOFATHER)
      IF (NV(INB).GE.NV(INS)) THEN
#else
      IF (NV(INB).GE.NV(INS)) THEN
#endif
       INS = INB
       GO TO 1070
      ENDIF
      INF = INB
 1090 INF = FRERE(INF)
      IF (INF.GT.0) GO TO 1090
      INF  = -INF
      INFS = -FILS(INF)
      IF (INFS.EQ.INS) THEN
        FILS(INF) = -INB
        IPS(INF)  = -INB
        FRERE(INS)  = FRERE(INB)
        FRERE(INB)  = INS
      ELSE
        INSW = INFS
 1100   INFS = FRERE(INSW)
        IF (INFS.NE.INS) THEN
          INSW = INFS
          GO TO 1100
        ENDIF
        FRERE(INS) = FRERE(INB)
        FRERE(INB) = INS
        FRERE(INSW)= INB
      ENDIF
        INS      = INB
        GO TO 1070
 1151 DO 51 I=1,N
       FILS(I) = IPS(I)
 51   CONTINUE
      IS = 1
      I = 0
      IPERM = 1
      DO 160 K=1,N
        IF (I.LE.0) THEN
         IF (NR.GT.N) EXIT
         I = NE(NR)
         NE(NR) = 0
         NR = NR + 1
         IL = N
         NA(N) = 0
        ENDIF
        DO 70 L=1,N
          IF (IPS(I).GE.0) EXIT
          ISON = -IPS(I)
          IPS(I) = 0
          I = ISON
          IL = IL - 1
          NA(IL) = 0
   70   CONTINUE
#if ! defined(NOAMALGTOFATHER)
        DADI = -IPE(I)  
        IF ( (DADI.NE.0) .AND.
     *      (
     *       (KEEP60.EQ.0).OR.
     *       ( (KEEP20.NE.DADI).AND.(KEEP38.NE.DADI) )
     *      ) 
     *     ) THEN
           ACCU = 
     &     ( dble(20000)*
     &       dble(NODE(I))*dble(NV(DADI)-NV(I)+NODE(I)) 
     &      )
     &      /
     &          ( dble(NV(DADI)+NODE(I))*
     &              dble(NV(DADI)+NODE(I))  )
           ACCU = ACCU + CUMUL(I) 
           AMALG_TO_father_OK = ( (NODE(I).LE.NEMIN).OR.
     &         (NODE(DADI).LE.NEMIN) )
           AMALG_TO_father_OK = ( AMALG_TO_father_OK .AND.
     &       ( 
     &        ( dble(2*(NODE(I)))*
     &         dble((NV(DADI)-NV(I)+NODE(I)))
     &        ) .LT.
     &        (  dble(NV(DADI)+NODE(I))*
     &         dble(NV(DADI)+NODE(I))*dble(NEMIN)/dble(100)
     &        ) 
     &       ) )
           AMALG_TO_father_OK = ( AMALG_TO_father_OK .AND.
     &     ( ACCU .LE. dble(NEMIN)*dble(100) )
     &                           )
           IF (AMALG_TO_father_OK) THEN
              CALL MUMPS_511(NV(I),NODE(I),NODE(I),
     *                                  KEEP50,1,FLOPS_SON)
              CALL MUMPS_511(NV(DADI),NODE(DADI),
     *                             NODE(DADI),
     *                             KEEP50,1,FLOPS_FATHER)
              FLOPS_AVANT = FLOPS_FATHER+FLOPS_SON
     *                      + max(dble(200.0) * dble(NV(I)-NODE(I))
     *                            * dble(NV(I)-NODE(I)),
     *                            dble(10000.0))
              CALL MUMPS_511(NV(DADI)+NODE(I),
     *                             NODE(DADI)+NODE(I),
     *                             NODE(DADI)+NODE(I),
     *                             KEEP50,1,FLOPS_APRES)
              IF(FLOPS_APRES .GT. FLOPS_AVANT) THEN
                 AMALG_TO_father_OK = .FALSE.
              ENDIF
           ENDIF
           IF (AMALG_TO_father_OK) THEN
              NAMALG(DADI) = NAMALG(DADI) +  NAMALG(I) + 1
              CUMUL(DADI)  =  CUMUL(DADI)+nint(ACCU)
           ELSE
              AMALG_TO_father_OK = (AMALG_TO_father_OK .OR.
     &                       ( NV(I)-NODE(I).EQ.NV(DADI)) )
              IF (AMALG_TO_father_OK) THEN
                NAMALG(DADI) = NAMALG(DADI) +  NAMALG(I)
                CUMUL(DADI)  = CUMUL(DADI)+nint(ACCU)
              ENDIF
           ENDIF
           IF (AMALG_TO_father_OK) THEN
             AMALG_COUNT = AMALG_COUNT+1
             IN = DADI
 75          IF (SUBORD(IN).EQ.0) GOTO 76
               IN = SUBORD(IN)
               GOTO 75
 76          CONTINUE
             SUBORD(IN) = I
             NV(I)      = 0
             IFSON = -FILS(DADI)
             IF (IFSON.EQ.I) THEN
              IF (FILS(I).LT.0) THEN
                FILS(DADI) =  FILS(I)
                GOTO 78
              ELSE
                IF (FRERE(I).GT.0) THEN
                  FILS(DADI) = -FRERE(I)  
                ELSE
                  FILS(DADI) = 0
                ENDIF
                GOTO 90
              ENDIF
             ENDIF
             IN = IFSON
  77         INS = IN
             IN = FRERE(IN)
             IF (IN.NE.I) GOTO 77
             IF (FILS(I) .LT.0) THEN
               FRERE(INS) = -FILS(I)
             ELSE
               FRERE(INS) = FRERE(I)  
               GOTO 90
             ENDIF
  78         CONTINUE
             IN = -FILS(I)
  79         INO = IN
             IN = FRERE(IN)
             IF (IN.GT.0) GOTO 79
             FRERE(INO) = FRERE(I)
  90         CONTINUE
             NODE(DADI) = NODE(DADI)+ NODE(I) 
             NV(DADI)   = NV(DADI) +  NODE(I) 
             NA(IL+1)   = NA(IL+1) + NA(IL)
             GOTO 120
           ENDIF
        ENDIF
#endif
        NE(IS) = NE(IS) + NODE(I) 
        IF (IL.LT.N) NA(IL+1) = NA(IL+1) + 1
        NA(IS) = NA(IL)
        ND(IS) = NV(I)
        NODE(I) = IS
        IPS(I) = IPERM
        IPERM = IPERM + 1
        IN = I
  777   IF (SUBORD(IN).EQ.0) GO TO 778
          IN = SUBORD(IN)
          NODE(IN) = IS
          IPS(IN) = IPERM
          IPERM = IPERM + 1
          GO TO 777
  778   IF (NA(IS).LE.0) GO TO 110
#if defined(NOAMALGTOFATHER)
        IF (   (KEEP60.NE.0).AND.
     &         (NE(IS).EQ.ND(IS)) ) GOTO 110
        IF (ND(IS-1)-NE(IS-1).EQ.ND(IS)) THEN
           GO TO 100
        ENDIF
        IF(NAMALG(IS-1) .GE. NAMALGMAX) THEN
           GOTO 110
        ENDIF
        IF ((NE(IS-1).GE.NEMIN).AND.
     &         (NE(IS).GE.NEMIN) ) GO TO 110
        IF (2*NE(IS-1)*(ND(IS)-ND(IS-1)+NE(IS-1)).GE.
     &    ((ND(IS)+NE(IS-1))*
     &    (ND(IS)+NE(IS-1))*NEMIN/100)) GO TO 110
        NAMALG(IS-1) = NAMALG(IS-1)+1
  100   NA(IS-1) = NA(IS-1) + NA(IS) - 1
        ND(IS-1) = ND(IS) + NE(IS-1)
        NE(IS-1) = NE(IS) + NE(IS-1)
        NE(IS) = 0
        NODE(I) = IS-1
        IFSON = -FILS(I)
        IN = IFSON
 102    INO = IN
        IN =  FRERE(IN)
        IF (IN.GT.0) GO TO 102
        NV(INO) = 0
        IN = I
  888   IF (SUBORD(IN).EQ.0) GO TO 889
        IN = SUBORD(IN)
        GO TO 888
  889   SUBORD(IN) = INO
      INOS = -FILS(INO)
      IF (IFSON.EQ.INO) THEN 
         FILS(I) = -INOS
         GO TO 107
      ENDIF
      IN = IFSON
 105  INS = IN
      IN =  FRERE(IN)
      IF (IN.NE.INO) GO TO 105
        IF (INOS.EQ.0) THEN
          FRERE(INS) = -I
          GO TO 120
        ELSE
          FRERE(INS) =  INOS
        ENDIF
 107    IN = INOS
        IF (IN.EQ.0) GO TO 120
 108    INT = IN
        IN =  FRERE(IN)
        IF (IN.GT.0) GO TO 108
        FRERE(INT) = -I
        GO TO 120
#endif
  110   IS = IS + 1
  120   IB = FRERE(I)
        IF (IB.GE.0) THEN
          IF (IB.GT.0) NA(IL) = 0
          I = IB
        ELSE
          I = -IB
          IL = IL + 1
        ENDIF
  160 CONTINUE
      NSTEPS = IS - 1
      DO I=1, N
        IF (NV(I).EQ.0) THEN
          FRERE(I) = N+1
        ELSE
          NFSIZ(I) = ND(NODE(I))
          IF (SUBORD(I) .NE.0) THEN
           INOS = -FILS(I)  
           INO = I
           DO WHILE (SUBORD(INO).NE.0) 
             IS = SUBORD(INO)
             FILS(INO) = IS
             INO = IS
           END DO
           FILS(INO) = -INOS
          ENDIF
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_557
#endif
      SUBROUTINE SMUMPS_201(NE, ND, NSTEPS,
     * MAXFR, MAXELIM, K50, MAXFAC, MAXNPIV,
     * K5,K6,PANEL_SIZE)
      IMPLICIT NONE
      INTEGER NSTEPS,MAXNPIV
      INTEGER MAXFR, MAXELIM, K50, MAXFAC
      INTEGER K5,K6,PANEL_SIZE
      INTEGER NE(NSTEPS), ND(NSTEPS)
      INTEGER ITREE, NFR, NELIM
      INTEGER LKJIB
      LKJIB   = max(k5,k6)
      MAXFR   = 0
      MAXFAC  = 0
      MAXELIM = 0
      MAXNPIV = 0
      PANEL_SIZE = 0
      DO ITREE=1,NSTEPS
        NELIM = NE(ITREE)
        NFR = ND(ITREE)
        IF (NFR.GT.MAXFR)         MAXFR   = NFR
        IF (NFR-NELIM.GT.MAXELIM) MAXELIM = NFR - NELIM
        IF (NELIM .GT. MAXNPIV) THEN
           IF(NFR .NE. NELIM) MAXNPIV = NELIM
        ENDIF
        IF (K50.EQ.0) THEN
          MAXFAC = max(MAXFAC, (2*NFR - NELIM)*NELIM )
          PANEL_SIZE = max(PANEL_SIZE, NFR*LKJIB)
        ELSE
         MAXFAC = max(MAXFAC, NFR * NELIM)
         PANEL_SIZE = max(PANEL_SIZE, NELIM*LKJIB)
         PANEL_SIZE = max(PANEL_SIZE, (NFR-NELIM)*LKJIB)
        ENDIF
      END DO
      RETURN
      END SUBROUTINE SMUMPS_201
      SUBROUTINE SMUMPS_348( N, FILS, FRERE,
     * NSTK, NA )
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N
      INTEGER, INTENT(IN)    :: FILS(N), FRERE(N)
      INTEGER, INTENT(INOUT) ::  NSTK(N), NA(N) 
      INTEGER NBROOT, NBLEAF, ILEAF, I, IN, ISON
      NA   = 0
      NSTK = 0
      NBROOT  = 0
      ILEAF   = 1
      DO 11 I=1,N
         IF (FRERE(I).EQ. N+1) CYCLE
         IF (FRERE(I).EQ.0) NBROOT = NBROOT + 1
         IN = I
 12      IN = FILS(IN)
         IF (IN.GT.0) GO TO 12
         IF (IN.EQ.0) THEN 
            NA(ILEAF) = I
            ILEAF     = ILEAF + 1
            CYCLE
         ENDIF
         ISON = -IN
 13      NSTK(I) = NSTK(I) + 1
         ISON = FRERE(ISON)
         IF (ISON.GT.0) GO TO 13
 11   CONTINUE
      NBLEAF = ILEAF-1
      IF (N.GT.1) THEN
         IF (NBLEAF.GT.N-2) THEN
            IF (NBLEAF.EQ.N-1) THEN
               NA(N-1) = -NA(N-1)-1
               NA(N)   = NBROOT
            ELSE
               NA(N) = -NA(N)-1
            ENDIF
         ELSE
            NA(N-1) = NBLEAF
            NA(N)   = NBROOT
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_348
      SUBROUTINE SMUMPS_202( N, NZ, PERM, IRN, ICN, PTRAR, KEEP )
      IMPLICIT NONE
      INTEGER, INTENT(IN)      :: N, NZ
      INTEGER, INTENT(IN)      :: KEEP(500)
      INTEGER, TARGET          :: PTRAR(N,2)
      INTEGER, INTENT(IN)      :: IRN(NZ), ICN(NZ), PERM(N)
      INTEGER I
      INTEGER IOLD, K, JOLD, INEW, JNEW, ISHIFT
      DO 50 IOLD=1,N
         PTRAR(IOLD,1) = 0
         PTRAR(IOLD,2) = 0
 50   CONTINUE
      DO 70 K=1,NZ
         IOLD = IRN(K)
         JOLD = ICN(K)
         IF ( (IOLD.GT.N).OR.(JOLD.GT.N).OR.(IOLD.LT.1)
     *        .OR.(JOLD.LT.1) ) GOTO 70
         IF (IOLD.NE.JOLD) THEN
            INEW = PERM(IOLD)
            JNEW = PERM(JOLD)
            IF ( KEEP( 50 ) .EQ. 0 ) THEN
               IF (INEW.LT.JNEW) THEN
                  PTRAR(IOLD,2) = PTRAR(IOLD,2) + 1
               ELSE
                  PTRAR(JOLD,1) = PTRAR(JOLD,1) + 1
               ENDIF
            ELSE
               IF ( INEW .LT. JNEW ) THEN
                  PTRAR( IOLD, 1 ) = PTRAR( IOLD, 1 ) + 1
               ELSE 
                  PTRAR( JOLD, 1 ) = PTRAR( JOLD, 1 ) + 1
               END IF
            ENDIF
         ENDIF
 70   CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_202
      SUBROUTINE SMUMPS_203( N, NZ, MTRANS, PERM,
     &     id, ICNTL, INFO)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (SMUMPS_STRUC) :: id
      INTEGER N, NZ, LIWG
      INTEGER PERM(N)
      INTEGER MTRANS
      INTEGER ICNTL(40), INFO(40)
      INTEGER  allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IW
      REAL, ALLOCATABLE, DIMENSION(:) :: S2
      TARGET :: S2
      INTEGER LS2,LSC
      INTEGER ICNTL64(10), INFO64(10)
      INTEGER ICNTL_SYM_MWM(10),INFO_SYM_MWM(10)
      REAL CNTL64(10)
      INTEGER LDW, LDWMIN
      INTEGER MPRINT,LP, MP, IPIW, LIW, LIWMIN
      INTEGER JPERM
      INTEGER NUMNZ, I, J, JPOS, K, NZREAL
      INTEGER PLENR, IP, IRNW,RSPOS,CSPOS
      LOGICAL PROK, IDENT, DUPPLI
      INTEGER NZTOT, K50, KER_SIZE, NZER_DIAG, MTRANSLOC,RZ_DIAG
      LOGICAL SCALINGLOC
      INTEGER,POINTER,DIMENSION(:) :: ZERODIAG
      INTEGER,POINTER,DIMENSION(:) :: STR_KER
      INTEGER,POINTER,DIMENSION(:) :: MARKED
      INTEGER,POINTER,DIMENSION(:) :: FLAG
      INTEGER,POINTER,DIMENSION(:) :: PIV_OUT
      REAL THEMIN, THEMAX, COLNORM,MAXDBL
      REAL ZERO,TWO,ONE
      PARAMETER(ZERO = 0.0E0,TWO = 2.0E0,ONE = 1.0E0)
      MPRINT = ICNTL(3)
      LP     = ICNTL(1)
      MP     = ICNTL(2)
      PROK = (MPRINT.GT.0)
      IF (PROK) WRITE(MPRINT,101)
 101  FORMAT(/'****** Preprocessing of original matrix '/)
      K50 = id%KEEP(50)
      SCALINGLOC = .FALSE.
      IF(id%KEEP(52) .EQ. -2) THEN
         IF(.not.associated(id%A)) THEN
            INFO(1) = -22
            INFO(2) = 4
            GOTO 500
         ELSE
            SCALINGLOC = .TRUE.
         ENDIF
      ELSE IF(id%KEEP(52) .EQ. 77) THEN
         SCALINGLOC = .TRUE.
         IF(K50 .NE. 2) THEN
            IF( MTRANS .NE. 5 .AND. MTRANS .NE. 6 
     *           .AND. MTRANS .NE. 7) THEN
               SCALINGLOC = .FALSE.
               IF (PROK) 
     *              WRITE(MPRINT,*) 'Analysis: auto scaling set OFF'
            ENDIF
         ENDIF
         IF(.not.associated(id%A)) THEN
            SCALINGLOC = .FALSE.
            IF (PROK) 
     *           WRITE(MPRINT,*) 'Analysis: auto scaling set OFF'
         ENDIF
      ENDIF
      IF(SCALINGLOC) THEN
         IF (PROK) WRITE(MPRINT,*) 
     *        'Scaling will be computed during analysis'
      ENDIF
      MTRANSLOC = MTRANS
      IF (MTRANS.LT.0 .OR. MTRANS.GT.7) GO TO 500
      IF (K50 .EQ. 0) THEN
         IF(.NOT. SCALINGLOC .AND. MTRANS .EQ. 7) THEN 
            GO TO 500
         ENDIF
         IF(SCALINGLOC) THEN
            MTRANSLOC = 5
         ENDIF
      ELSE
         IF (MTRANS .EQ. 7) MTRANSLOC = 5
      ENDIF
      IF(SCALINGLOC .AND. MTRANSLOC .NE. 5 .AND.
     *     MTRANSLOC .NE. 6 ) THEN
         IF (PROK) WRITE(MPRINT,*)
     *        'WARNING scaling required: set MTRANS option to 5'
         MTRANSLOC = 5
      ENDIF
      IF (N.EQ.1) THEN
        MTRANS=0
        GO TO 500
      ENDIF
      IF(K50 .EQ. 2) THEN
         NZTOT = 2*NZ+N
      ELSE
         NZTOT = NZ
      ENDIF
      ZERODIAG => id%IS1(N+1:2*N)
      STR_KER => id%IS1(2*N+1:3*N)
      CALL SMUMPS_448(ICNTL64,CNTL64)
      ICNTL64(1) = ICNTL(1)
      ICNTL64(2) = ICNTL(2)
      ICNTL64(3) = ICNTL(2)
      ICNTL64(4) = -1
      IF (ICNTL(4).EQ.3) ICNTL64(4) = 0
      IF (ICNTL(4).EQ.4) ICNTL64(4) = 1
      ICNTL64(5) = -1
      IF (PROK) THEN
         WRITE(MPRINT,'(A,I3)')
     *     'Compute maximum matching (Maximum Transversal):',
     *        MTRANSLOC
         IF (MTRANSLOC.EQ.1)
     *   WRITE(MPRINT,'(A,I3)')' ... JOB =',MTRANSLOC
         IF (MTRANSLOC.EQ.2)
     *   WRITE(MPRINT,'(A,I3,A)')
     *     ' ... JOB =',MTRANSLOC,': BOTTLENECK THESIS'
         IF (MTRANSLOC.EQ.3)
     *   WRITE(MPRINT,'(A,I3,A)')
     *     ' ... JOB =',MTRANSLOC,': BOTTLENECK SIMAX'
         IF (MTRANSLOC.EQ.4)
     *   WRITE(MPRINT,'(A,I3,A)')
     *     ' ... JOB =',MTRANSLOC,': MAXIMIZE SUM DIAGIONAL'
         IF (MTRANSLOC.EQ.5 .OR. MTRANSLOC.EQ.6)
     *   WRITE(MPRINT,'(A,I3,A)')
     *     ' ... JOB =',MTRANSLOC,
     *     ': MAXIMIZE PRODUCT DIAGONAL AND SCALE'
      ENDIF
      id%INFOG(23) = MTRANSLOC
      CNTL64(2) = huge(CNTL64(2))
      IRNW = 1
      IP = IRNW + NZTOT
      PLENR = IP + N + 1
      IPIW = PLENR
      IF (MTRANSLOC.EQ.1) LIWMIN = 5*N
      IF (MTRANSLOC.EQ.2) LIWMIN = 4*N
      IF (MTRANSLOC.EQ.3) LIWMIN = 10*N + NZTOT
      IF (MTRANSLOC.EQ.4) LIWMIN = 5*N
      IF (MTRANSLOC.EQ.5) LIWMIN = 5*N
      IF (MTRANSLOC.EQ.6) LIWMIN = 5*N + NZTOT
      LIW = LIWMIN
      LIWG  = LIW + (NZTOT + N + 1)
      ALLOCATE(IW(LIWG), stat=allocok)
      IF (allocok .GT. 0 ) GOTO 410
      IF (MTRANSLOC.EQ.1) THEN
       LDWMIN = N+3
      ENDIF
      IF (MTRANSLOC.EQ.2) LDWMIN = max(N+NZTOT,N+3)
      IF (MTRANSLOC.EQ.3) LDWMIN = max(NZTOT+1,N+3)
      IF (MTRANSLOC.EQ.4) LDWMIN = 2*N + max(NZTOT,N+3)
      IF (MTRANSLOC.EQ.5) LDWMIN = 3*N + NZTOT
      IF (MTRANSLOC.EQ.6) LDWMIN = 4*N + NZTOT
      LDW   = LDWMIN
      ALLOCATE(S2(LDW), stat=allocok)
      IF(MTRANSLOC .NE. 1) LDW = LDW-NZTOT
      RSPOS = NZTOT
      CSPOS = RSPOS+N
      IF (allocok .GT. 0 ) GOTO 430
      NZREAL = 0
      DO 5 J=1,N
        IW(PLENR+J-1) = 0
  5   CONTINUE
      IF(K50 .EQ. 0) THEN
         DO 10 K=1,NZ
            I = id%IRN(K)
            J = id%JCN(K)
            IF ( (J.LE.N).AND.(J.GE.1).AND.
     *           (I.LE.N).AND.(I.GE.1) ) THEN
               IW(PLENR+J-1) = IW(PLENR+J-1) + 1
               NZREAL = NZREAL + 1
            ENDIF
 10      CONTINUE
      ELSE
         ZERODIAG = 0
         NZER_DIAG = N
         RZ_DIAG = 0
         DO K=1,NZ
            I = id%IRN(K)
            J = id%JCN(K)
            IF ( (J.LE.N).AND.(J.GE.1).AND.
     *           (I.LE.N).AND.(I.GE.1) ) THEN
               IW(PLENR+J-1) = IW(PLENR+J-1) + 1
               NZREAL = NZREAL + 1
               IF(I .NE. J) THEN
                  IW(PLENR+I-1) = IW(PLENR+I-1) + 1
                  NZREAL = NZREAL + 1
               ELSE
                  IF(ZERODIAG(I) .EQ. 0) THEN
                     ZERODIAG(I) = K
                     IF(associated(id%A)) THEN
                        IF(abs(id%A(K)) .EQ. real(0.0E0)) THEN
                           RZ_DIAG = RZ_DIAG + 1
                        ENDIF
                     ENDIF
                     NZER_DIAG = NZER_DIAG - 1                     
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         IF(MTRANSLOC .GE. 4) THEN
            DO I =1, N
               IF(ZERODIAG(I) .EQ. 0) THEN
                  IW(PLENR+I-1) = IW(PLENR+I-1) + 1
                  NZREAL = NZREAL + 1
               ENDIF
            ENDDO
         ENDIF
      ENDIF
      IW(IP)   = 1
      DO 20 J=1,N
        IW(IP+J)   = IW(IP+J-1)+IW(PLENR+J-1)
  20  CONTINUE
      DO 25 J=1, N
        IW(PLENR+J-1 ) = IW(IP+J-1 )
  25  CONTINUE
      IF(K50 .EQ. 0) THEN
         IF (MTRANSLOC.EQ.1) THEN
            DO 30 K=1,NZ
               I = id%IRN(K)
               J = id%JCN(K)
               IF ( (J.LE.N).AND.(J.GE.1) .AND.
     *              (I.LE.N).AND.(I.GE.1)) THEN
                  JPOS            = IW(PLENR+J-1)
                  IW(IRNW+JPOS-1) = I
                  IW(PLENR+J-1)   = IW(PLENR+J-1) + 1
              ENDIF
 30         CONTINUE
         ELSE
            IF ( .not.associated(id%A)) THEN
               INFO(1) = -22
               INFO(2) = 4
               GOTO 500
            ENDIF
            DO 35 K=1,NZ
               I = id%IRN(K)
               J = id%JCN(K)
               IF ( (J.LE.N).AND.(J.GE.1) .AND.
     *              (I.LE.N).AND.(I.GE.1)) THEN
                  JPOS            = IW(PLENR+J-1)
                  IW(IRNW+JPOS-1) = I
                  S2(JPOS)         = abs(id%A(K))
                  IW(PLENR+J-1)   = IW(PLENR+J-1) + 1
               ENDIF
 35         CONTINUE
         ENDIF
      ELSE
         IF (MTRANSLOC.EQ.1) THEN
            DO K=1,NZ
               I = id%IRN(K)
               J = id%JCN(K)
               IF ( (J.LE.N).AND.(J.GE.1) .AND.
     *              (I.LE.N).AND.(I.GE.1)) THEN
                  JPOS            = IW(PLENR+J-1)
                  IW(IRNW+JPOS-1) = I
                  IW(PLENR+J-1)   = IW(PLENR+J-1) + 1
                  IF(I.NE.J) THEN
                     JPOS            = IW(PLENR+I-1)
                     IW(IRNW+JPOS-1) = J
                     IW(PLENR+I-1)   = IW(PLENR+I-1) + 1
                  ENDIF
               ENDIF
            ENDDO
         ELSE
            IF ( .not.associated(id%A)) THEN
               INFO(1) = -22
               INFO(2) = 4
               GOTO 500
            ENDIF
            K = 1
            THEMIN = ZERO
            DO
               IF(THEMIN .NE. ZERO) EXIT
               THEMIN = abs(id%A(K))
               K = K+1
            ENDDO
            THEMAX = THEMIN
            DO K=1,NZ
               I = id%IRN(K)
               J = id%JCN(K)
               IF ( (J.LE.N).AND.(J.GE.1) .AND.
     *              (I.LE.N).AND.(I.GE.1)) THEN
                  JPOS            = IW(PLENR+J-1)
                  IW(IRNW+JPOS-1) = I
                  S2(JPOS)         = abs(id%A(K))
                  IW(PLENR+J-1)   = IW(PLENR+J-1) + 1
                  IF(abs(id%A(K)) .GT. THEMAX) THEN
                     THEMAX = abs(id%A(K))
                  ELSE IF(abs(id%A(K)) .LT. THEMIN 
     *                    .AND. abs(id%A(K)).GT. ZERO) THEN
                     THEMIN = abs(id%A(K))
                  ENDIF
                  IF(I.NE.J) THEN
                     JPOS            = IW(PLENR+I-1)
                     IW(IRNW+JPOS-1) = J
                     S2(JPOS)         = abs(id%A(K))
                     IW(PLENR+I-1)   = IW(PLENR+I-1) + 1
                  ENDIF
               ENDIF
            ENDDO
            DO I =1, N
               IF(ZERODIAG(I) .EQ. 0) THEN
                  JPOS            = IW(PLENR+I-1)
                  IW(IRNW+JPOS-1) = I
                  S2(JPOS)         = ZERO
                  IW(PLENR+I-1)   = IW(PLENR+I-1) + 1
               ENDIF
            ENDDO
            CNTL64(2) = (log(THEMAX/THEMIN))*(real(N))
     *           - log(THEMIN) + ONE
         ENDIF
      ENDIF
      DUPPLI = .FALSE.
      I = NZREAL
      FLAG => id%IS1(3*N+1:4*N)
      IF(MTRANSLOC.NE.1) THEN
         CALL SMUMPS_563(N,NZREAL,IW(IP),IW(IRNW),S2,
     *        PERM,FLAG)
      ELSE
         CALL SMUMPS_562(N,NZREAL,IW(IP),IW(IRNW),
     *        PERM,FLAG)
      ENDIF
      IF(NZREAL .NE. I) DUPPLI = .TRUE.
      LS2 = NZTOT
      IF ( MTRANSLOC .EQ. 1 ) THEN
         LS2 = 1
         LDW = 1
      ENDIF
      CALL SMUMPS_559(MTRANSLOC ,N, N, NZREAL, 
     *     IW(IP), IW(IRNW), S2(1), LS2,
     *     NUMNZ, PERM, LIW, IW(IPIW), LDW, S2(LS2+1),
     *     ICNTL64, CNTL64, INFO64)
      IF (INFO64(1).LT.0) THEN
         IF (LP.GT.0 .AND. ICNTL(4).GE.1)
     *        WRITE(LP,'(A,I5)')
     *   ' INTERNAL ERROR in MAXTRANS INFO(1)=',INFO64(1)
         INFO(1) = -9964
         INFO(2) = INFO64(1)
         GO TO 500
      ENDIF
      IF (INFO64(1).GT.0) THEN
         IF (MP.GT.0 .AND. ICNTL(4).GE.2)
     *        WRITE(MP,'(A,I5)')
     *        ' WARNING in MAXTRANS INFO(1)=',INFO64(1)
      ENDIF
      KER_SIZE = 0
      IF(K50 .EQ. 2) THEN
         DO I=1,N
            IF(ZERODIAG(I) .EQ. 0) THEN
               IF(PERM(I) .EQ. I) THEN
                  KER_SIZE = KER_SIZE + 1
                  PERM(I) = -I
                  STR_KER(KER_SIZE) = I
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      IF (NUMNZ.LT.N) GO TO 400
      IF(KER_SIZE .GT. 0) GOTO 400
      IF(K50 .EQ. 0) THEN
         IDENT = .TRUE.
         IF (MTRANS .EQ. 0 ) GOTO 102
         DO 80 J=1,N
            JPERM = PERM(J)
            IW(PLENR+JPERM-1) = J
            IF (JPERM.NE.J) IDENT = .FALSE.
 80      CONTINUE
         IF(IDENT) THEN 
            MTRANS = 0
         ELSE
            IF(MTRANS .EQ. 7) THEN
               MTRANS = -9876543
               GOTO 102
            ENDIF
            IF (PROK) WRITE(MPRINT,'(A)')
     *           ' ... Apply column permutation'
            DO 100 K=1,NZ
               J = id%JCN(K)
               IF ((J.LE.0).OR.(J.GT.N)) GO TO 100
               id%JCN(K) = IW(PLENR+J-1)
 100        CONTINUE
            IF (MP.GT.0 .AND. ICNTL(4).GE.2)
     *           WRITE(MP,'(/A)')
     *           ' WARNING input matrix data modified'
         ENDIF
 102     CONTINUE
         IF (SCALINGLOC) THEN
            IF ( associated(id%COLSCA))
     *           DEALLOCATE( id%COLSCA )
            IF ( associated(id%ROWSCA))
     *           DEALLOCATE( id%ROWSCA )
            ALLOCATE( id%COLSCA(N), stat=allocok)
            IF (allocok .GT.0) THEN
               id%INFO(1)=-5
               id%INFO(2)=N
               IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
                  WRITE (LP,'(/A)') '** Error in SMUMPS_203'
                  WRITE (LP,'(A)') 
     *                 '** Failure during allocation of COLSCA'
                  GOTO 500
               ENDIF
            ENDIF
            ALLOCATE( id%ROWSCA(N), stat=allocok)
            IF (allocok .GT.0) THEN
               id%INFO(1)=-5
               id%INFO(2)=N
               IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
                  WRITE (LP,'(/A)') '** Error in SMUMPS_203'
                  WRITE (LP,'(A)') 
     *                 '** Failure during allocation of ROWSCA'
                  GOTO 500
               ENDIF
            ENDIF
            id%KEEP(52) = -2
            id%KEEP(74) = 1
            MAXDBL = log(huge(MAXDBL))
            DO J=1,N
               IF(S2(RSPOS+J) .GT. MAXDBL) THEN
                  S2(RSPOS+J) = ZERO
               ENDIF
               IF(S2(CSPOS+J) .GT. MAXDBL) THEN
                  S2(CSPOS+J)= ZERO
               ENDIF
            ENDDO
            DO 105 J=1,N
               id%ROWSCA(J) = exp(S2(RSPOS+J))
               IF(id%ROWSCA(J) .EQ. ZERO) THEN
                  id%ROWSCA(J) = ONE
               ENDIF
               IF ( MTRANS .EQ.  -9876543 .OR. MTRANS.EQ. 0 ) THEN
                 id%COLSCA(J)= exp(S2(CSPOS+J))
                 IF(id%COLSCA(J) .EQ. ZERO) THEN
                   id%COLSCA(J) = ONE
                 ENDIF
               ELSE
                 id%COLSCA(IW(PLENR+J-1))= exp(S2(CSPOS+J))
                 IF(id%COLSCA(IW(PLENR+J-1)) .EQ. ZERO) THEN
                   id%COLSCA(IW(PLENR+J-1)) = ONE
                 ENDIF
               ENDIF
 105        CONTINUE
         ENDIF
      ELSE
         IDENT = .FALSE.         
         IF(SCALINGLOC) THEN
            IF ( associated(id%COLSCA)) DEALLOCATE( id%COLSCA )
            IF ( associated(id%ROWSCA)) DEALLOCATE( id%ROWSCA )
            ALLOCATE( id%COLSCA(N), stat=allocok)
            IF (allocok .GT.0) THEN
               id%INFO(1)=-5
               id%INFO(2)=N
               IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
                  WRITE (LP,'(/A)') '** Error in SMUMPS_203'
                  WRITE (LP,'(A)') 
     *                 '** Failure during allocation of COLSCA'
                  GOTO 500
               ENDIF
            ENDIF
            ALLOCATE( id%ROWSCA(N), stat=allocok)
            IF (allocok .GT.0) THEN
               id%INFO(1)=-5
               id%INFO(2)=N
               IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
                  WRITE (LP,'(/A)') '** Error in SMUMPS_203'
                  WRITE (LP,'(A)') 
     *                 '** Failure during allocation of ROWSCA'
                  GOTO 500
               ENDIF
            ENDIF
            id%KEEP(52) = -2
            id%KEEP(74) = 1
            MAXDBL = log(huge(MAXDBL))
            DO J=1,N
               IF(S2(RSPOS+J)+S2(CSPOS+J) .GT. MAXDBL) THEN
                  S2(RSPOS+J) = ZERO
                  S2(CSPOS+J)= ZERO
               ENDIF
            ENDDO
            DO J=1,N
               IF(PERM(J) .GT. 0) THEN
                  id%ROWSCA(J) = 
     *                 exp((S2(RSPOS+J)+S2(CSPOS+J))/TWO)
                  IF(id%ROWSCA(J) .EQ. ZERO) THEN
                     id%ROWSCA(J) = ONE
                  ENDIF
                  id%COLSCA(J)= id%ROWSCA(J)
               ENDIF
            ENDDO
            DO JPOS=1,KER_SIZE
               I = STR_KER(JPOS)
               COLNORM = ZERO
               DO J = IW(IP+I-1),IW(IP+J) - 1
                  COLNORM = max(COLNORM,S2(J))
               ENDDO
               COLNORM = exp(COLNORM) 
               id%ROWSCA(I) = ONE / COLNORM
               id%COLSCA(I) = id%ROWSCA(I)
            ENDDO
         ENDIF
         IF(MTRANS .EQ. 7 .OR. id%KEEP(95) .EQ. 0) THEN
            IF( (NZER_DIAG+RZ_DIAG) .LT. (N/10) 
     *           .AND. id%KEEP(95) .EQ. 0) THEN
               MTRANS = 0
               id%KEEP(95) = 1
               GOTO 390
            ELSE
               IF(id%KEEP(95) .EQ. 0) THEN
                 IF(SCALINGLOC) THEN
                  id%KEEP(95) = 3
                 ELSE
                  id%KEEP(95) = 2   
                 ENDIF
               ENDIF
               IF(MTRANS .EQ. 7) MTRANS = 5
            ENDIF
         ENDIF
         IF(MTRANS .EQ. 0) GOTO 390
         ICNTL_SYM_MWM = 0
         INFO_SYM_MWM = 0
         IF(MTRANS .EQ. 5 .OR. MTRANS .EQ. 6 .OR.
     *        MTRANS .EQ. 7) THEN
            ICNTL_SYM_MWM(1) = 0
            ICNTL_SYM_MWM(2) = 1
         ELSE IF(MTRANS .EQ. 4) THEN
            ICNTL_SYM_MWM(1) = 2
            ICNTL_SYM_MWM(2) = 1
         ELSE
            ICNTL_SYM_MWM(1) = 0
            ICNTL_SYM_MWM(2) = 1
         ENDIF
         MARKED => id%IS1(2*N+1:3*N)
         FLAG => id%IS1(3*N+1:4*N)
         PIV_OUT => id%IS1(4*N+1:5*N)
         IF(MTRANSLOC .LT. 4) THEN
            LSC = 1
         ELSE
            LSC = 2*N
         ENDIF
         CALL SMUMPS_551(
     *        N, NZREAL, IW(IP), IW(IRNW), S2(1),LSC, PERM, ZERODIAG,
     *        ICNTL_SYM_MWM, S2(LSC+1),MARKED,FLAG,
     *        PIV_OUT, INFO_SYM_MWM)
         IF(INFO_SYM_MWM(1) .NE. 0) THEN
            WRITE(*,*) '** Error in SMUMPS_203'
            RETURN
         ENDIF
         IF(INFO_SYM_MWM(3) .EQ. N) THEN
            IDENT = .TRUE.
         ELSEIF( (N-INFO_SYM_MWM(4)-INFO_SYM_MWM(3)) .GT. N/10
     *           ) THEN
            IDENT = .TRUE.
            id%KEEP(95) = 1
         ELSE
            DO I=1,N
               PERM(I) = PIV_OUT(I)
            ENDDO
         ENDIF
         id%KEEP(93) = INFO_SYM_MWM(4)
         id%KEEP(94) = INFO_SYM_MWM(3)
         IF (IDENT) MTRANS=0
      ENDIF
 390  IF(MTRANS .EQ. 0) THEN
         id%KEEP(95) = 1 
         IF (PROK) THEN
           WRITE (MPRINT,'(A)')
     *  ' ... Column permutation not used'
         ENDIF
      ENDIF
      GO TO 500
 400  IF ((LP.GE.0).AND.(ICNTL(4).GE.1))
     *   WRITE (LP,'(/A)') '** Error: Matrix is structurally singular'
      INFO(1) = -6
      INFO(2) = NUMNZ
      GOTO 500
 410  IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
       WRITE (LP,'(/A)') '** Error in SMUMPS_203'
       WRITE (LP,'(A,I9)')
     & '** Failure during allocation of INTEGER array of size ',
     & LIWG
      ENDIF
      INFO(1) = -5
      INFO(2) = LIWG
      GOTO 500
 430  IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
       WRITE (LP,'(/A)') '** Error in SMUMPS_203'
       WRITE (LP,'(A)') '** Failure during allocation of S2'
      ENDIF
      INFO(1) = -5
      INFO(2) = LDW
 500  CONTINUE
      IF (allocated(IW)) DEALLOCATE(IW)
      IF (allocated(S2)) DEALLOCATE(S2)
      RETURN
      END SUBROUTINE SMUMPS_203
      SUBROUTINE SMUMPS_21( N, MYID,PROCNODE, ND,
     *                         SSARBR, NBSA, FRERE, FILS,
     *                         SLAVEF,
     *                         KEEP,KEEP8, ICNTL, INFO )
      IMPLICIT NONE
      INTEGER KEEP(500), ICNTL(40), INFO(40)
      INTEGER*8 KEEP8(150)
      INTEGER N, NBSA, SLAVEF, MYID
      INTEGER PROCNODE( N ), ND( N )
      INTEGER SSARBR( NBSA ), FRERE( N ), FILS( N )
      LOGICAL  MUMPS_167, MUMPS_283
      INTEGER  MUMPS_330, MUMPS_275, MUMPS_497
      EXTERNAL MUMPS_330, MUMPS_275,
     *         MUMPS_167, MUMPS_283, MUMPS_497
      INTEGER I, INODE, INODEROOT, IN, NPIV, NFRONT, NFR, NB_NIV2
      INTEGER IROOTTREE, SIZEROOT, PROCVAL, NSLAVES, KMAX
      INTEGER NNODE, ICT, TPN, MP, allocok
      LOGICAL ONE_LEVEL2
      INTEGER, ALLOCATABLE :: NUMNODE( : ), COST( : ),
     *                        MEMUSED( : ), NODEMAP(:)
      MP = ICNTL(3)
      SIZEROOT  = -1
      IROOTTREE = -1
      NB_NIV2  = 0
      IF ( SLAVEF .eq. 1 ) THEN
       DO INODE = 1, N
          PROCNODE(INODE) = PROCNODE( INODE ) - SLAVEF
       END DO
      ELSE
       DO I = 1, NBSA
         INODEROOT = SSARBR( I )
         PROCNODE( INODEROOT ) = PROCNODE( INODEROOT )
     *                         - SLAVEF
         INODE = INODEROOT
 30      CONTINUE
         IN = INODE
         DO WHILE ( IN .NE. 0 )
           INODE = IN
           DO WHILE ( IN .GT. 0 )
             IN = FILS( IN )
           END DO
           IF ( IN .LT. 0 ) IN = - IN
         END DO
  10     CONTINUE
         IF ( INODE .EQ. INODEROOT ) GOTO 100
         PROCNODE( INODE ) = PROCNODE( INODE ) - 2 * SLAVEF
         IN = FRERE( INODE )
         INODE = abs( IN )
         IF ( IN .LT. 0 ) THEN
             GO TO 10
         ELSE
             GO TO 30
         END IF
 100     CONTINUE
        END DO
      END IF
      ONE_LEVEL2 = .FALSE.
      DO INODE = 1, N
        IF ( FRERE( INODE ) .LT. N + 1 .AND.
     *     .NOT. MUMPS_167( INODE, PROCNODE, SLAVEF ) ) THEN
          IN = INODE
          NPIV = 0
          DO WHILE ( IN.GT.0 )
            IN = FILS( IN )
            NPIV = NPIV + 1
          END DO
          NFRONT = ND( INODE )
          IF ( (MUMPS_283( INODE,PROCNODE,SLAVEF).AND.
     &     (IN.NE.0.OR.NFRONT.NE.NPIV )) .OR. SLAVEF.EQ.1 ) GOTO 991
          IF ( NFRONT-NPIV > KEEP(9) .AND.
     *         ( NPIV > KEEP(4) .OR. KEEP(45) .EQ. 0 ) 
     *         .AND. ICNTL(40) .EQ. 0 ) THEN
            PROCNODE( INODE ) = PROCNODE( INODE ) + SLAVEF
            ONE_LEVEL2 = .TRUE.
            NB_NIV2 = NB_NIV2 + 1
          END IF
          IF ( NFRONT .EQ. NPIV ) THEN
            IF ( NFRONT .GT. SIZEROOT ) THEN
              IROOTTREE = INODE
              SIZEROOT  = NFRONT
            END IF
          END IF
        END IF
 991    CONTINUE
        IF ( ( SLAVEF .eq. 1 .OR.
     *         MUMPS_283( INODE, PROCNODE, SLAVEF ) )
     *         .AND. KEEP(53) .NE. 0 .AND.
     *         FRERE( INODE ) .LT. N + 1 ) THEN
          IN = INODE
          NPIV = 0
          DO WHILE ( IN.GT.0 )
            IN = FILS( IN )
            NPIV = NPIV + 1
          END DO
          NFRONT = ND( INODE )
          IF ( NFRONT .EQ. NPIV ) THEN
            IF ( NFRONT .GT. SIZEROOT ) THEN
              IROOTTREE = INODE
              SIZEROOT  = NFRONT
            END IF
          END IF
        END IF
      END DO
      KEEP ( 56 ) = NB_NIV2
      IF ( SIZEROOT .GT. KEEP( 37 )
     *       .and. (ICNTL(13) .LE. 0 )
     *       .and. KEEP(60) .EQ. 0
#ifndef null_space_old
     *       .and. KEEP(53) .eq. 0
#endif
     *  ) THEN
         IF ( MP .GT. 0 )
     *   WRITE( MP, * ) 'A root of estimated size ',SIZEROOT,
     *   ' has been selected for Scalapack'
        KEEP(38) = IROOTTREE
        PROCNODE( IROOTTREE ) = 1 + 2 * SLAVEF
      ELSE IF (KEEP(60) .EQ. 2 .OR. KEEP(60).EQ.3) THEN
        PROCNODE( KEEP(38) ) = 1 + 2 * SLAVEF
        IF ( MP .GT. 0 )
     *  WRITE( MP, * ) ' Largest root node of size ',SIZEROOT
      ELSE
        KEEP(38) = 0
        IF ( MP .GT. 0 )
     *  WRITE( MP, * ) ' Largest root node of size ',SIZEROOT
      END IF
      DO INODE = 1, N
        IF ( FRERE( INODE ) .LT. N + 1 ) THEN
          IN = FILS( INODE )
          DO WHILE ( IN > 0 )
            PROCNODE( IN ) = PROCNODE( INODE )
            IN = FILS( IN )
          END DO
        END IF
      END DO
      IF ( KEEP( 45 ) .NE. 0 ) RETURN
      NNODE = 0
      DO INODE = 1, N
        IF ( FRERE( INODE ) .LT. N + 1 ) THEN
          TPN = MUMPS_330( INODE, PROCNODE, SLAVEF )
          IF ( TPN .NE. 3 .AND.
     *         .not. MUMPS_283( INODE, PROCNODE, SLAVEF ) .AND.
     *         .not. MUMPS_167( INODE, PROCNODE, SLAVEF ) ) THEN
               NNODE = NNODE + 1
          END IF
        END IF
      END DO
      IF ( NNODE .GT. 0 ) THEN
        ALLOCATE( NUMNODE( NNODE  ), stat = allocok )
        if (allocok .gt. 0) THEN
          INFO(1) = -7
          INFO(2) = NNODE
          RETURN
        endif
        ALLOCATE( NODEMAP( NNODE  ), stat = allocok )
        if (allocok .gt. 0) THEN
          INFO(1) = -7
          INFO(2) = NNODE
          RETURN
        endif
        ALLOCATE( COST   ( NNODE  ), stat = allocok )
        if (allocok .gt. 0) THEN
          INFO(1) = -7
          INFO(2) = NNODE
          RETURN
        endif
        ALLOCATE( MEMUSED( SLAVEF ), stat = allocok )
        if (allocok .gt. 0) THEN
          INFO(1) = -7
          INFO(2) = SLAVEF
          RETURN
        endif
        MEMUSED( 1: SLAVEF ) = 0
        ICT = 0
        DO INODE = 1, N
        IF ( FRERE( INODE ) .LT. N + 1 ) THEN
          TPN = MUMPS_330( INODE, PROCNODE, SLAVEF )
          IF ( TPN .NE. 3 .AND.
     *       .not. MUMPS_283( INODE, PROCNODE, SLAVEF ) .AND.
     *       .not. MUMPS_167( INODE, PROCNODE, SLAVEF ) ) THEN
             ICT            = ICT + 1
             NUMNODE( ICT ) = INODE
             NFR  = ND( INODE )
             NPIV = 0
             IN   = INODE
             DO WHILE ( IN .GT. 0 )
               NPIV = NPIV + 1
               IN = FILS( IN )
             END DO
             IF ( TPN .eq. 1 ) THEN
               IF ( KEEP( 50 ) .eq. 0 ) THEN
                 COST( ICT ) = NPIV * ( 2 * NFR - NPIV )
               ELSE
                 COST( ICT ) = NPIV * NFR
               END IF
             ELSE
               KMAX  = MUMPS_497( KEEP(10), NFR - NPIV)
               NSLAVES     = max( min((NFR-NPIV)/KMAX, SLAVEF-1),
     *                            1 )
               IF ( KEEP(50) .eq. 0 ) THEN
                 COST( ICT ) = NPIV * NFR - NPIV *
     *           ( ( NFR - NPIV ) / NSLAVES )
               ELSE
                 COST( ICT ) = NPIV * NPIV - NPIV *
     *           ( ( NFR - NPIV ) / NSLAVES )
               END IF
             END IF
          ELSE IF ( MUMPS_283( INODE, PROCNODE, SLAVEF ) .OR.
     *              MUMPS_167( INODE, PROCNODE, SLAVEF )  ) THEN
               NFR  = ND( INODE )
               NPIV = 0
               IN   = INODE
               DO WHILE ( IN .GT. 0 )
                 NPIV = NPIV + 1
                 IN = FILS( IN )
               END DO
               IF ( KEEP( 50 ) .eq. 0 ) THEN
                 MEMUSED(1+MUMPS_275( INODE, PROCNODE, SLAVEF))
     *           = MEMUSED(1+MUMPS_275( INODE, PROCNODE, SLAVEF))
     *           + NPIV * ( 2 * NFR - NPIV )
               ELSE
                 MEMUSED(1+MUMPS_275( INODE, PROCNODE, SLAVEF))
     *           = MEMUSED(1+MUMPS_275( INODE, PROCNODE, SLAVEF))
     *           + NPIV * NFR
               END IF
          END IF
        END IF
        END DO
        CALL SMUMPS_212( N, NNODE, NUMNODE, COST, MEMUSED,
     *                               PROCNODE, SLAVEF, NODEMAP )
        DO I = 1, NNODE
          INODE   = NUMNODE( I )
          PROCVAL = PROCNODE( INODE ) + ( NODEMAP( I ) - 1 )
     *            - MUMPS_275( INODE, PROCNODE, SLAVEF )
          IN = INODE
          DO WHILE ( IN > 0 )
            PROCNODE( IN ) = PROCVAL
            IN = FILS( IN )
          END DO
        END DO
        DEALLOCATE( NUMNODE )
        DEALLOCATE( NODEMAP )
        DEALLOCATE( COST    )
        DEALLOCATE( MEMUSED )
      END IF
      RETURN
      END SUBROUTINE SMUMPS_21
      SUBROUTINE SMUMPS_212( N, NNODE, NUMNODE, COST,
     *       MEMUSED, PROCNODE, SLAVEF, NODEMAP )
      IMPLICIT NONE
      INTEGER N, NNODE, SLAVEF
      INTEGER PROCNODE( N )
      INTEGER NUMNODE( NNODE ), COST( NNODE ), NODEMAP( NNODE )
      INTEGER MEMUSED( SLAVEF )
      INTEGER IMAXLOAD, IMINLOAD, IPROC, INODE
      CALL SMUMPS_149
     *       ( MEMUSED, SLAVEF, IMINLOAD, IMAXLOAD )
      CALL SMUMPS_260( NNODE, NUMNODE, COST )
      DO INODE = 1, NNODE
        IF ( COST( INODE ) > 0 ) THEN
          IPROC = IMINLOAD
        ELSE
          IPROC = IMAXLOAD
        END IF
        MEMUSED( IPROC ) = MEMUSED( IPROC ) + COST( INODE )
        NODEMAP( INODE ) = IPROC
        CALL SMUMPS_149
     *       ( MEMUSED, SLAVEF, IMINLOAD, IMAXLOAD )
      END DO
      RETURN
      END SUBROUTINE SMUMPS_212
      SUBROUTINE SMUMPS_149
     *           ( MEMUSED, SLAVEF, IMINLOAD, IMAXLOAD )
      INTEGER SLAVEF, IMINLOAD, IMAXLOAD
      INTEGER MEMUSED( SLAVEF )
      INTEGER MINLOAD, MAXLOAD, IPROC
      MINLOAD = MEMUSED( 1 )
      MAXLOAD = MEMUSED( 1 )
      IMINLOAD = 1
      IMAXLOAD = 1
      DO IPROC = 2, SLAVEF
        IF ( MEMUSED( IPROC ) .GT. MAXLOAD ) THEN
          MAXLOAD  = MEMUSED( IPROC )
          IMAXLOAD = IPROC
        END IF
        IF ( MEMUSED( IPROC ) .LT. MINLOAD ) THEN
          MINLOAD  = MEMUSED( IPROC )
          IMINLOAD = IPROC
        END IF
      END DO  
      RETURN
      END SUBROUTINE SMUMPS_149
      SUBROUTINE SMUMPS_260( NNODE, NUMNODE, COST )
      IMPLICIT NONE
      INTEGER NNODE
      INTEGER NUMNODE( NNODE ), COST( NNODE )
      LOGICAL FIN
      INTEGER I, LOC_COST, LOC_NUM
      FIN = .FALSE.
      DO WHILE ( .NOT. FIN )
        FIN = .TRUE.
        DO I = 1, NNODE - 1
          IF ( abs( COST( I ) ) .LT. abs( COST( I + 1 ) ) ) THEN
            LOC_COST      = COST( I )
            COST( I )     = COST( I + 1 )
            COST( I + 1 ) = LOC_COST
            LOC_NUM          = NUMNODE( I )
            NUMNODE( I )     = NUMNODE( I + 1 )
            NUMNODE( I + 1 ) = LOC_NUM
            FIN = .FALSE.
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE SMUMPS_260
      SUBROUTINE SMUMPS_100
     *( MYID, COMM, KEEP,KEEP8, INFO, INFOG, RINFO, RINFOG, ICNTL )
      IMPLICIT NONE
      INTEGER COMM, MYID, KEEP(500), INFO(40), ICNTL(40), INFOG(40)
      INTEGER*8 KEEP8(150)
      REAL RINFO(20), RINFOG(20)
      INCLUDE 'mpif.h'
      INTEGER MASTER, IERR_MPI, MPG
      PARAMETER( MASTER = 0 )
      MPG = ICNTL(3)
      IF ( MYID.eq.MASTER.and.MPG.GT.0) THEN
       WRITE(MPG, 99992) INFO(1), INFO(2),
     *  KEEP8(109), KEEP8(111), INFOG(4),
     *  INFOG(5), KEEP(28), INFOG(7), KEEP(23), ICNTL(7), 
     *  KEEP(12), KEEP(56), KEEP(61), RINFOG(1)
       IF (KEEP(95).GT.1)             
     *                    WRITE(MPG, 99993) KEEP(95)  
       IF (KEEP(54).GT.0) WRITE(MPG, 99994) KEEP(54)  
       IF (KEEP(60).GT.0) WRITE(MPG, 99995) KEEP(60)  
      ENDIF
      RETURN
99992 FORMAT(/'Leaving analysis phase with  ...'/
     1       'INFOG(1)                                       =',I16/
     2       'INFOG(2)                                       =',I16/
     3       ' -- (20) Number of entries in factors (estim.) =',I16/
     3       ' --  (3) Storage of factors  (REAL, estimated) =',I16/
     3       ' --  (4) Storage of factors  (INT , estimated) =',I16/
     5       ' --  (5) Maximum frontal size      (estimated) =',I16/
     6       ' --  (6) Number of nodes in the tree           =',I16/
     6       ' --  (7) Ordering option effectively used      =',I16/
     7       'ICNTL(6) Maximum transversal option            =',I16/
     8       'ICNTL(7) Pivot order option                    =',I16/
     8       'Percentage of memory relaxation (effective)    =',I16/
     8       'Number of level 2 nodes                        =',I16/
     8       'Number of split nodes                          =',I16/
     9   'RINFO(1) Operations during elimination (estim) =  ',1PD10.3)
99993 FORMAT('Ordering compressed/constrained (ICNTL(12))    =',I16)
99994 FORMAT('Distributed matrix entry format (ICNTL(18))    =',I16)
99995 FORMAT('Effective Schur option (ICNTL(19))             =',I16)
      END SUBROUTINE SMUMPS_100
      SUBROUTINE SMUMPS_97
     *           ( N, FRERE, FILS, NFSIZ, NSTEPS, NSLAVES, 
     *             KEEP, KEEP8, SPLITROOT, MP, LDIAG, INFO1, INFO2 )
      IMPLICIT NONE
      INTEGER N, NSTEPS, NSLAVES, KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER FRERE( N ), FILS( N ), NFSIZ( N )
      LOGICAL SPLITROOT
      INTEGER MP, LDIAG  
      INTEGER INFO1, INFO2
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPOOL 
      INTEGER INODE, DEPTH, I, IBEG, IEND, IIPOOL, NROOT
      INTEGER MAX_DEPTH, ISON, TOT_CUT, MAX_CUT, STRAT, K79
      INTEGER NFRONT, K82, allocok
      K79  = KEEP(79)
      K82  = abs(KEEP(82))  
      STRAT=KEEP(62)
      IF (KEEP(210).EQ.1) THEN
        MAX_DEPTH = 2*NSLAVES*K82
        STRAT     = STRAT/4
      ELSE
        IF (( NSLAVES .eq. 1 ).AND. (.NOT. SPLITROOT) ) RETURN
        MAX_DEPTH = int( log( real( NSLAVES - 1 ) ) 
     &                 / log(2.0E0) )
      ENDIF
      ALLOCATE(IPOOL(NSTEPS+1), stat=allocok)
      IF (allocok.GT.0) THEN
        INFO1= -7
        INFO2= NSTEPS+1
        RETURN
      ENDIF
      NROOT = 0
      DO INODE = 1, N
        IF ( FRERE(INODE) .eq. 0 ) THEN
          NROOT = NROOT + 1
          IPOOL( NROOT ) = INODE
        END IF
      END DO
      IBEG = 1
      IEND = NROOT
      IIPOOL   = NROOT + 1
      IF (SPLITROOT) MAX_DEPTH=1
      DO DEPTH = 1, MAX_DEPTH
        DO I = IBEG, IEND
          INODE = IPOOL( I )
          ISON = INODE
          DO WHILE ( ISON .GT. 0 )
            ISON = FILS( ISON )
          END DO
          ISON = - ISON
          DO WHILE ( ISON .GT. 0 )
            IPOOL( IIPOOL ) = ISON
            IIPOOL = IIPOOL + 1
            ISON = FRERE( ISON )
          END DO
        END DO
        IPOOL( IBEG ) = -IPOOL( IBEG )
        IBEG = IEND + 1
        IEND = IIPOOL - 1
      END DO
      IPOOL( IBEG ) = -IPOOL( IBEG )
      TOT_CUT = 0
      IF (SPLITROOT) THEN
        MAX_CUT = NROOT*max(K82,2)
        INODE = abs(IPOOL(1))
        NFRONT = NFSIZ( INODE )
        K79     = max(real(NFRONT*NFRONT)/
     &            real((K82+1)*(K82+1)), real(1))
      ELSE
         MAX_CUT = 2 * NSLAVES
         IF (KEEP(210).EQ.1) THEN
            MAX_CUT = 4 * (MAX_CUT + 4)
         ENDIF
      ENDIF
      DEPTH   = -1
      DO I = 1, IIPOOL - 1
        INODE = IPOOL( I )
        IF ( INODE .LT. 0 ) THEN
          INODE = -INODE
          DEPTH = DEPTH + 1
        END IF
        CALL SMUMPS_313
     *           ( INODE, N, FRERE, FILS, NFSIZ, NSTEPS, NSLAVES,
     *             KEEP,KEEP8, TOT_CUT, STRAT, DEPTH, 
     *             K79, SPLITROOT, MP, LDIAG )
        IF ( TOT_CUT > MAX_CUT )  EXIT
      END DO
      KEEP(61) = TOT_CUT
      DEALLOCATE(IPOOL)
      RETURN
      END SUBROUTINE SMUMPS_97
      RECURSIVE SUBROUTINE SMUMPS_313
     * ( INODE, N, FRERE, FILS, NFSIZ, NSTEPS, NSLAVES, KEEP,KEEP8,
     *   TOT_CUT, STRAT, DEPTH, K79, SPLITROOT, MP, LDIAG )
      IMPLICIT NONE
      INTEGER INODE, N, NSTEPS, NSLAVES, KEEP(500), STRAT, 
     &        DEPTH, TOT_CUT, K79, MP, LDIAG
      INTEGER*8 KEEP8(150)
      INTEGER FRERE( N ), FILS( N ), NFSIZ( N )
      LOGICAL SPLITROOT
      INTEGER I, IN, NPIV, NFRONT, NSLAVES_ESTIM
      REAL WK_SLAVE, WK_MASTER
      INTEGER INODE_SON, INODE_FATH, IN_SON, IN_FATH, IN_GRANDFATH
      INTEGER NPIV_SON, NPIV_FATH, MAX_CUT
      INTEGER NCB, KAVERAGE, NSLAVESMIN, NSLAVESMAX
      INTEGER  MUMPS_50,
     *         MUMPS_52
      EXTERNAL  MUMPS_50,
     *         MUMPS_52
      IF  ( (KEEP(210).EQ.1.AND.KEEP(60).EQ.0) .OR.
     &       (SPLITROOT) ) THEN
        IF ( FRERE ( INODE ) .eq. 0 ) THEN 
          NFRONT = NFSIZ( INODE )
          NPIV = NFRONT
          NCB = 0
          IF (NFRONT*NFRONT.GT.K79) THEN 
           GOTO 333
          ELSE IF ( (SPLITROOT).AND.
     &      (MP.GT.0 .AND. LDIAG.GE.2)) THEN
             WRITE(MP,'(A,I11)') 
     &        ' Order of root node after splitting : ',
     &        NFRONT
          ENDIF
        ENDIF
      ENDIF
      IF ( FRERE ( INODE ) .eq. 0 ) RETURN
      NFRONT = NFSIZ( INODE )
      IN = INODE
      NPIV = 0
      DO WHILE( IN > 0 )
        IN = FILS( IN )
        NPIV = NPIV + 1
      END DO
      NCB = NFRONT - NPIV
      IF ( (NFRONT - (NPIV/2)) .LE. KEEP(9)) RETURN
      IF ((KEEP(50) == 0.and.NFRONT * NPIV > K79 ) .OR.
     &(KEEP(50) .NE.0.and.NPIV * NPIV > K79 )) GOTO 333
      IF (KEEP(210).EQ.1) THEN
        NSLAVESMIN    = 1   
        NSLAVESMAX    = 64  
        NSLAVES_ESTIM = 32+NSLAVES
      ELSE
        NSLAVESMIN = MUMPS_50 
     *         ( NSLAVES, KEEP(48), KEEP(10), KEEP(50),
     *         NFRONT, NCB)
        NSLAVESMAX = MUMPS_52 
     *        ( NSLAVES, KEEP(48), KEEP(10), KEEP(50),
     *          NFRONT, NCB)
        NSLAVES_ESTIM = max (1, 
     *   nint( real(NSLAVESMAX-NSLAVESMIN)/real(3) )
     *                    )
        NSLAVES_ESTIM = min (NSLAVES_ESTIM, NSLAVES-1)
      ENDIF
      IF ( KEEP(50) .eq. 0 ) THEN
       WK_MASTER = 0.6667E0 * 
     *                real(NPIV)*real(NPIV)*real(NPIV) +
     *                real(NPIV)*real(NPIV)*real(NCB)
       WK_SLAVE  = real( NPIV ) * real( NCB ) *
     *         ( 2.0E0 * real(NFRONT) - real(NPIV) )
     *         / real(NSLAVES_ESTIM)
      ELSE
       WK_MASTER = real(NPIV)*real(NPIV)*real(NPIV) / real(3)
       WK_SLAVE  = 
     *           (real(NPIV)*real(NCB)*real(NFRONT)) 
     *           / real(NSLAVES_ESTIM)
      ENDIF
      IF (KEEP(210).EQ.1) THEN
        IF ( ( 100 + STRAT )
     &        * WK_SLAVE / 100 .GE. WK_MASTER ) RETURN
      ELSE
        IF ( ( 100 + STRAT * max( DEPTH-1, 1 ) )
     &        * WK_SLAVE / 100 .GE. WK_MASTER ) RETURN
      ENDIF
 333  CONTINUE
      IF (NPIV .LE. 1 ) RETURN
       NSTEPS  = NSTEPS + 1
       TOT_CUT = TOT_CUT + 1
       NPIV_SON  = max(NPIV/2,1)
       NPIV_FATH = NPIV - NPIV_SON
       INODE_SON = INODE
       IN_SON = INODE
       DO I = 1, NPIV_SON - 1
         IN_SON = FILS( IN_SON )
       END DO
       INODE_FATH = FILS( IN_SON )
       IF ( INODE_FATH .LT. 0 ) THEN
       write(*,*) 'Error: INODE_FATH < 0 ', INODE_FATH
       END IF
       IN_FATH = INODE_FATH
       DO WHILE ( FILS( IN_FATH ) > 0 )
         IN_FATH = FILS( IN_FATH )
       END DO
       FRERE( INODE_FATH ) = FRERE( INODE_SON )
       FRERE( INODE_SON  ) = - INODE_FATH
       FILS ( IN_SON     ) = FILS( IN_FATH )
       FILS ( IN_FATH    ) = - INODE_SON
       IN = FRERE( INODE_FATH )
       DO WHILE ( IN > 0 )
           IN = FRERE( IN )
       END DO
       IF ( IN .eq. 0 )  GO TO 10
       IN = -IN
       DO WHILE ( FILS( IN ) > 0 )
           IN = FILS( IN )
       END DO
       IN_GRANDFATH = IN
       IF ( FILS( IN_GRANDFATH ) .eq. - INODE_SON ) THEN
           FILS( IN_GRANDFATH ) = -INODE_FATH
       ELSE
           IN = IN_GRANDFATH
           IN = - FILS ( IN )
           DO WHILE ( FRERE( IN ) > 0 )
             IF ( FRERE( IN ) .eq. INODE_SON ) THEN
               FRERE( IN ) = INODE_FATH
               GOTO 10
             END IF
             IN = FRERE( IN )
           END DO
           WRITE(*,*) 'ERROR 2 in SPLIT NODE',
     *          IN_GRANDFATH, IN, FRERE(IN)
       END IF
 10    CONTINUE
       NFSIZ(INODE_SON) = NFRONT
       NFSIZ(INODE_FATH) = NFRONT - NPIV_SON
       KEEP(2) = max( KEEP(2), NFRONT - NPIV_SON )
      CALL SMUMPS_313
     *  ( INODE_FATH, N, FRERE, FILS, NFSIZ, NSTEPS,
     *   NSLAVES, KEEP,KEEP8, TOT_CUT, STRAT, DEPTH, 
     *   K79, SPLITROOT, MP, LDIAG )
      IF (.NOT. SPLITROOT) THEN
        CALL SMUMPS_313
     *   ( INODE_SON, N, FRERE, FILS, NFSIZ, NSTEPS,
     *   NSLAVES, KEEP,KEEP8, TOT_CUT, STRAT, DEPTH, 
     *   K79, SPLITROOT, MP, LDIAG )
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_313
      SUBROUTINE SMUMPS_351
     * (N,NZ, IRN, ICN, IW, LW, IPE, LEN,
     * IQ, FLAG, IWFR,
     * NRORM, NIORM, IFLAG,IERROR, ICNTL, 
     * symmetry, SYM, MedDens, NBQD, AvgDens)
      INTEGER N,NZ,LW,IFLAG,IERROR,NRORM,NIORM,IWFR
      INTEGER symmetry, SYM
      INTEGER MedDens, NBQD, AvgDens
      INTEGER ICNTL(40)
      INTEGER  IRN(NZ), ICN(NZ) 
      INTEGER LEN(N)
      INTEGER IPE(N+1)
      INTEGER FLAG(N), IW(LW)
      INTEGER IQ(N)
      INTEGER MP, MPG
      INTEGER I,K,J,N1,LAST,NDUP,K1,K2,L
      INTEGER NBERR, THRESH
      INTEGER NZOFFA, NDIAGA
      REAL RSYM
      INTRINSIC nint
      MP = ICNTL(2)
      MPG= ICNTL(3)
      NIORM  = 3*N
      NDIAGA = 0
      IERROR = 0
      DO 10 I=1,N
        IPE(I) = 0
   10 CONTINUE
      DO 50 K=1,NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     *                          .OR.(J.LT.1)) THEN
           IERROR = IERROR + 1
        ELSE
          IF (I.NE.J) THEN
           IPE(I) = IPE(I) + 1
           IPE(J) = IPE(J) + 1
           NIORM  = NIORM + 1
          ELSE
           NDIAGA = NDIAGA + 1
          ENDIF
        ENDIF
   50 CONTINUE
      NZOFFA  = NIORM - 3*N
      IF (IERROR.GE.1) THEN
         NBERR  = 0
         IF (mod(IFLAG,2) .EQ. 0) IFLAG = IFLAG+1
         IF ((MP.GT.0).AND.(ICNTL(4).GE.2))  THEN 
          WRITE (MP,99999) 
          DO 70 K=1,NZ
           I = IRN(K)
           J = ICN(K)
           IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     *                            .OR.(J.LT.1)) THEN
            NBERR = NBERR + 1
            IF (NBERR.LE.10)  THEN
               IF (mod(K,10).GT.3 .OR. mod(K,10).EQ.0 .OR.
     *             (10.LE.K .AND. K.LE.20)) THEN
                 WRITE (MP,'(I8,A,I8,A,I8,A)')
     *             K,'th entry (in row',I,' and column',J,') ignored'
               ELSE
                 IF (mod(K,10).EQ.1) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'st entry (in row',I,' and column',J,') ignored'
                 IF (mod(K,10).EQ.2) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'nd entry (in row',I,' and column',J,') ignored'
                 IF (mod(K,10).EQ.3) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'rd entry (in row',I,' and column',J,') ignored'
               ENDIF
            ELSE
               GO TO 100
            ENDIF
           ENDIF
   70     CONTINUE
         ENDIF
      ENDIF
  100 NRORM = NIORM - 2*N
      IQ(1) = 1
      N1 = N - 1
      IF (N1.GT.0) THEN
        DO 110 I=1,N1
            IQ(I+1) = IPE(I) + IQ(I) 
  110   CONTINUE
      ENDIF
      LAST = max(IPE(N)+IQ(N)-1,IQ(N))
      FLAG(1:N) = 0
      IPE(1:N)  = IQ(1:N)
      IW(1:LAST) = 0
      IWFR = LAST + 1
      DO 200 K=1,NZ
         I = IRN(K)
         J = ICN(K)
         IF (I.NE.J) THEN
          IF (I.LT.J) THEN
            IF ((I.GE.1).AND.(J.LE.N)) THEN
             IW(IQ(I)) = -J
             IQ(I)     = IQ(I) + 1 
            ENDIF
          ELSE
            IF ((J.GE.1).AND.(I.LE.N)) THEN
             IW(IQ(J)) = -I
             IQ(J)     = IQ(J) + 1
            ENDIF
          ENDIF
         ENDIF
  200 CONTINUE
      NDUP = 0
      DO 260 I=1,N
        K1 = IPE(I) 
        K2 = IQ(I) -1
        IF (K1.GT.K2) THEN
         LEN(I) = 0
         IQ(I)  = 0
        ELSE
         DO 240 K=K1,K2
           J     = -IW(K)
           IF (J.LE.0) GO TO 250
           L     = IQ(J) 
           IQ(J) = L + 1
           IF (FLAG(J).EQ.I) THEN
            NDUP = NDUP + 1
            IW(L) = 0
            IW(K) = 0
           ELSE
            IW(L)   = I
            IW(K)   = J
            FLAG(J) = I
           ENDIF
  240    CONTINUE
  250    IQ(I) = IQ(I) - IPE(I)
         IF (NDUP.EQ.0) LEN(I) = IQ(I)
        ENDIF
  260 CONTINUE
      IF (NDUP.NE.0) THEN
       IWFR = 1
       DO 280 I=1,N
         IF (IQ(I).EQ.0) THEN
             LEN(I) = 0
            IPE(I) = IWFR
            GOTO 280
         ENDIF
         K1 = IPE(I) 
         K2 = K1 + IQ(I) - 1
         L = IWFR
         IPE(I) = IWFR
         DO 270 K=K1,K2
           IF (IW(K).NE.0) THEN
            IW(IWFR) = IW(K)
            IWFR     = IWFR + 1
           ENDIF
  270    CONTINUE
         LEN(I) = IWFR - L 
  280  CONTINUE
      ENDIF
      IPE(N+1) = IPE(N) + LEN(N)
      IWFR = IPE(N+1)
      IF (SYM.EQ.0) THEN
      RSYM =  real(NDIAGA+2*NZOFFA - (IWFR-1))/
     *            real(NZOFFA+NDIAGA) 
      symmetry = nint (100.0*RSYM)
         IF (MPG .GT. 0)
     &  write(MPG,'(A,I5)') 
     &  ' ... Structural symmetry (in percent)=', symmetry
        IF (MP.GT.0 .AND. MPG.NE.MP)
     &  write(MP,'(A,I5)') 
     &  ' ... Structural symmetry (in percent)=', symmetry
      ELSE
       symmetry = 100
      ENDIF
      AvgDens = nint(real(IWFR-1)/real(N))
      THRESH  = AvgDens*50 - AvgDens/10 + 1
      NBQD    = 0
      IF (N.GT.2) THEN
        IQ(1:N) = 0
        DO I= 1, N
          K = max(LEN(I),1)
          IQ(K) = IQ(K) + 1
          IF (K.GT.THRESH) NBQD = NBQD+1
        ENDDO
        K = 0
        MedDens = 0
        DO WHILE (K .LT. (N/2))
         MedDens = MedDens + 1
         K       = K+IQ(MedDens)
        ENDDO
      ELSE
        MedDens = AvgDens
      ENDIF
         IF (MPG .GT. 0)
     &  write(MPG,'(A,3I5)') 
     &  ' Density: NBdense, Average, Median   =',
     &  NBQD, AvgDens, MedDens
        IF (MP.GT.0 .AND. MPG.NE.MP)
     &  write(MP,'(A,3I5)') 
     &  ' Density: NBdense, Average, Median   =',
     &  NBQD, AvgDens, MedDens
      RETURN
99999 FORMAT (/'*** Warning message from analysis routine ***')
      END SUBROUTINE SMUMPS_351
      SUBROUTINE SMUMPS_701(N, SYM, NPROCS, IORD,
     *                       symmetry,MedDens, NBQD, AvgDens,
     *                       PROK, MP)
      IMPLICIT NONE
      INTEGER, intent(in)    :: N, NPROCS, SYM
      INTEGER, intent(in)    :: symmetry,MedDens, NBQD, AvgDens, MP
      LOGICAL, intent(in)    :: PROK
      INTEGER, intent(inout)   :: IORD
      INTEGER MAXQD
      PARAMETER (MAXQD=2)
      INTEGER SMALLSYM, SMALLUNS
      PARAMETER (SMALLUNS=5000, SMALLSYM=10000)
#if ! defined(metis)
      IF ( IORD .EQ. 5 ) THEN
        IF (PROK) WRITE(MP,*)
     *  'WARNING: METIS not available. Ordering set to default.'
        IORD = 7
      END IF
#endif
#if ! defined(pord)
      IF ( IORD .EQ. 4 ) THEN
        IF (PROK) WRITE(MP,*)
     *  'WARNING: PORD not available. Ordering set to default.'
        IORD = 7
      END IF
#endif
#if ! defined(scotch)
      IF ( IORD .EQ. 3 ) THEN
        IF (PROK) WRITE(MP,*)
     *  'WARNING: SCOTCH not available. Ordering set to default.'
        IORD = 7
      END IF
#endif
      IF (IORD.EQ.7) THEN
        IF (SYM.NE.0) THEN
          IF ( N.LE.SMALLSYM ) THEN 
             IF (NBQD.GE.MAXQD) THEN
               IORD = 6         
             ELSE
               IORD = 2         
             ENDIF
          ELSE
             IF (NBQD.GE.MedDens*NPROCS) THEN
               IORD = 6      
               RETURN
             ENDIF
#if  defined(metis)
             IORD = 5
#else
#if defined(pord)
               IORD = 4
#else
               IORD = 6
#endif
#endif
          ENDIF
        ELSE
          IF ( N.LE.SMALLUNS ) THEN
            IF (NBQD.GE.MAXQD) THEN
              IORD = 6  
            ELSE
              IORD = 2   
            ENDIF
          ELSE
            IF (NBQD.GE.MedDens*NPROCS) THEN
              IORD = 6      
              RETURN
            ENDIF
#if  defined(metis)
            IORD = 5
#else
#if defined(pord)
              IORD = 4
#else
              IORD = 6
#endif
#endif
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_701
      SUBROUTINE SMUMPS_510
     *     (KEEP10, KEEP2, KEEP48 ,KEEP50, NSLAVES)
      IMPLICIT NONE
      INTEGER NSLAVES, KEEP10, KEEP2, KEEP48, KEEP50
      KEEP10 = max(KEEP10*KEEP2,1)
#if defined(t3e) 
      KEEP10 = min(1500000, KEEP10)
#elif defined(SP_)
      KEEP10 = min(3000000, KEEP10)
#else
      KEEP10 = min(2000000, KEEP10)
#endif
#if defined(t3e) 
      IF (NSLAVES .GT. 64) THEN
         KEEP10 = 
     &        min((8*(KEEP2*KEEP2)/NSLAVES)+1, KEEP10)
      ELSE
         KEEP10 = 
     &        min((4*(KEEP2*KEEP2)/NSLAVES)+1, KEEP10)
      ENDIF 
#else
      IF (NSLAVES.GT.64) THEN
         KEEP10 = 
     &        min((6*(KEEP2*KEEP2)/NSLAVES)+1, KEEP10)
      ELSE
         KEEP10 = 
     &        min((4*(KEEP2*KEEP2)/NSLAVES)+1, KEEP10)
      ENDIF
#endif
         IF (KEEP50 .EQ. 0 ) THEN
            KEEP10 = max(KEEP10,((7*KEEP2*KEEP2)/
     *           4/max(NSLAVES-1,1)) + KEEP2)
         ELSE
            KEEP10 = max(KEEP10,((7*KEEP2*KEEP2)/
     *           4/(max(NSLAVES-1,1))) +KEEP2)
         ENDIF
      IF (KEEP50 .EQ. 0 ) THEN
#if defined(t3e)
         KEEP10 = max(KEEP10,200000)
#else 
         KEEP10 = max(KEEP10,300000)
#endif
      ELSE
#if defined(t3e)
         KEEP10 = max(KEEP10,40000)
#else 
         KEEP10 = max(KEEP10,80000)
#endif
      ENDIF
      KEEP10 = -KEEP10 
      RETURN
      END SUBROUTINE SMUMPS_510
      SUBROUTINE SMUMPS_559(JOB,M,N,NE,
     *     IP,IRN,A,LA,NUM,PERM,LIW,IW,LDW,DW,
     *     ICNTL,CNTL,INFO)
      IMPLICIT NONE
      INTEGER NICNTL, NCNTL, NINFO
      PARAMETER (NICNTL=10, NCNTL=10, NINFO=10)
      INTEGER JOB,M,N,NE,NUM,LIW,LDW
      INTEGER IP(N+1),IRN(NE),PERM(M),IW(LIW)
      INTEGER ICNTL(NICNTL),INFO(NINFO)
      INTEGER LA
      REAL A(LA)
      REAL DW(LDW),CNTL(NCNTL)
      INTEGER I,J,K,WARN1,WARN2,WARN4
      REAL FACT,ZERO,ONE,RINF,RINF2,RINF3
      PARAMETER (ZERO=0.0E+00,ONE=1.0E+0)
      EXTERNAL SMUMPS_457,SMUMPS_444,SMUMPS_451,
     &         SMUMPS_452,SMUMPS_454
      INTRINSIC abs,log
      RINF = CNTL(2)
      RINF2 = huge(RINF2)/real(2*N)
      RINF3 = 0.0E0
      WARN1 = 0
      WARN2 = 0
      WARN4 = 0
      IF (JOB.LT.1 .OR. JOB.GT.6) THEN
         INFO(1) = -1
         INFO(2) = JOB
         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
         GO TO 99
      ENDIF
      IF (M.LT.1 .OR. M.LT.N) THEN
         INFO(1) = -2
         INFO(2) = M
         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'M',M
         GO TO 99
      ENDIF
      IF (N.LT.1) THEN
         INFO(1) = -2
         INFO(2) = N
         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
         GO TO 99
      ENDIF
      IF (NE.LT.1) THEN
         INFO(1) = -3
         INFO(2) = NE
         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NE',NE
         GO TO 99
      ENDIF
      IF (JOB.EQ.1) K = 4*N +   M
      IF (JOB.EQ.2) K = 2*N + 2*M
      IF (JOB.EQ.3) K = 8*N + 2*M + NE
      IF (JOB.EQ.4) K = 3*N + 2*M
      IF (JOB.EQ.5) K = 3*N + 2*M
      IF (JOB.EQ.6) K = 3*N + 2*M + NE
      IF (LIW.LT.K) THEN
         INFO(1) = -4
         INFO(2) = K
         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
         GO TO 99
      ENDIF
      IF (JOB.GT.1) THEN
         IF (JOB.EQ.2) K =       M
         IF (JOB.EQ.3) K = 1
         IF (JOB.EQ.4) K =     2*M
         IF (JOB.EQ.5) K = N + 2*M
         IF (JOB.EQ.6) K = N + 3*M
         IF (LDW.LT.K) THEN
            INFO(1) = -5
            INFO(2) = K
            IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
            GO TO 99
         ENDIF
      ENDIF
      IF (ICNTL(5).EQ.0) THEN
         DO 3 I = 1,M
            IW(I) = 0
 3       CONTINUE
         DO 6 J = 1,N
            DO 4 K = IP(J),IP(J+1)-1
               I = IRN(K)
               IF (I.LT.1 .OR. I.GT.M) THEN
                  INFO(1) = -6
                  INFO(2) = J
                  IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
                  GO TO 99
               ENDIF
               IF (IW(I).EQ.J) THEN
                  INFO(1) = -7
                  INFO(2) = J
                  IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I 
                  GO TO 99
               ELSE
                  IW(I) = J
               ENDIF
 4          CONTINUE
 6       CONTINUE
      ENDIF
      IF (ICNTL(3).GE.0) THEN
         IF (ICNTL(4).EQ.0 .OR. ICNTL(4).EQ.1) THEN
            WRITE(ICNTL(3),9020) JOB,M,N,NE
            IF (ICNTL(4).EQ.0) THEN
               WRITE(ICNTL(3),9021) (IP(J),J=1,min(10,N+1))
               WRITE(ICNTL(3),9022) (IRN(J),J=1,min(10,NE))
               IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(J),J=1,min(10,NE))
            ELSEIF (ICNTL(4).EQ.1) THEN
               WRITE(ICNTL(3),9021) (IP(J),J=1,N+1)
               WRITE(ICNTL(3),9022) (IRN(J),J=1,NE)
               IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(J),J=1,NE)
            ENDIF
            WRITE(ICNTL(3),9024) (ICNTL(J),J=1,NICNTL)
            WRITE(ICNTL(3),9025) (CNTL(J),J=1,NCNTL)
         ENDIF
      ENDIF
      DO 8 I=1,NINFO
         INFO(I) = 0
    8 CONTINUE
      IF (JOB.EQ.1) THEN
         DO 10 J = 1,N
            IW(J) = IP(J+1) - IP(J)
 10      CONTINUE
         CALL SMUMPS_457(M,N,IRN,NE,IP,IW(1),PERM,NUM,
     &        IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1))
         GO TO 90
      ENDIF
      IF (JOB.EQ.2) THEN
         DW(1) = max(ZERO,CNTL(1))
         CALL SMUMPS_444(M,N,NE,IP,IRN,A,PERM,NUM,
     &        IW(1),IW(N+1),IW(2*N+1),IW(2*N+M+1),DW,RINF2)
         GO TO 90
      ENDIF
      IF (JOB.EQ.3) THEN
         DO 20 K = 1,NE
            IW(K) = IRN(K)
 20      CONTINUE
         CALL SMUMPS_451(N,NE,IP,IW,A)
         FACT = max(ZERO,CNTL(1))
         CALL SMUMPS_452(M,N,NE,IP,IW(1),A,PERM,NUM,IW(NE+1),
     &        IW(NE+N+1),IW(NE+2*N+1),IW(NE+3*N+1),IW(NE+4*N+1),
     &        IW(NE+5*N+1),IW(NE+5*N+M+1),FACT,RINF2)
         GO TO 90
      ENDIF
      IF (JOB.EQ.4) THEN
         DO 50 J = 1,N
            FACT = ZERO
            DO 30 K = IP(J),IP(J+1)-1
               IF (abs(A(K)).GT.FACT) FACT = abs(A(K))
 30         CONTINUE
            IF(FACT .GT. RINF3) RINF3 = FACT
            DO 40 K = IP(J),IP(J+1)-1
               A(K) = FACT - abs(A(K))
 40         CONTINUE
 50      CONTINUE
         DW(1) = max(ZERO,CNTL(1))
         DW(2) = RINF3
         IW(1) = JOB
         CALL SMUMPS_454(M,N,NE,IP,IRN,A,PERM,NUM,
     &        IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1),
     &        DW(1),DW(M+1),RINF2)
         GO TO 90
      ENDIF
      IF (JOB.EQ.5 .or. JOB.EQ.6) THEN
         RINF3=ONE
         IF (JOB.EQ.5) THEN
            DO 75 J = 1,N
               FACT = ZERO
               DO 60 K = IP(J),IP(J+1)-1
                  IF (A(K).GT.FACT) FACT = A(K)
 60            CONTINUE
               DW(2*M+J) = FACT
               IF (FACT.NE.ZERO) THEN
                  FACT = log(FACT)
                  IF(FACT .GT. RINF3) RINF3=FACT
                  DO 70 K = IP(J),IP(J+1)-1
                     IF (A(K).NE.ZERO) THEN
                        A(K) = FACT - log(A(K))
                        IF(A(K) .GT. RINF3) RINF3=A(K)
                     ELSE
                        A(K) = FACT + RINF
                     ENDIF
 70               CONTINUE
               ELSE
                  DO 71 K = IP(J),IP(J+1)-1
                     A(K) = ONE
 71               CONTINUE
               ENDIF
 75         CONTINUE
         ENDIF
         IF (JOB.EQ.6) THEN
            DO 175 K = 1,NE
               IW(3*N+2*M+K) = IRN(K)
 175        CONTINUE
            DO 61 I = 1,M
               DW(2*M+N+I) = ZERO
 61         CONTINUE
            DO 63 J = 1,N
               DO 62 K = IP(J),IP(J+1)-1
                  I = IRN(K)
                  IF (A(K).GT.DW(2*M+N+I)) THEN
                     DW(2*M+N+I) = A(K)
                  ENDIF
 62            CONTINUE
 63         CONTINUE
            DO 64 I = 1,M
               IF (DW(2*M+N+I).NE.ZERO) THEN
                  DW(2*M+N+I) = 1/DW(2*M+N+I)
               ENDIF
 64         CONTINUE
            DO 66 J = 1,N
               DO 65 K = IP(J),IP(J+1)-1
                  I = IRN(K)
                  A(K) = DW(2*M+N+I) * A(K)
 65            CONTINUE
 66         CONTINUE
            CALL SMUMPS_451(N,NE,IP,IW(3*N+2*M+1),A)
            DO 176 J = 1,N
               IF (IP(J).NE.IP(J+1)) THEN
                  FACT = A(IP(J))
               ELSE
                  FACT = ZERO
               ENDIF
               DW(2*M+J) = FACT
               IF (FACT.NE.ZERO) THEN
                  FACT = log(FACT)
                  DO 170 K = IP(J),IP(J+1)-1
                     IF (A(K).NE.ZERO) THEN
                        A(K) = FACT - log(A(K))
                        IF(A(K) .GT. RINF3) RINF3=A(K)
                     ELSE
                        A(K) = FACT + RINF
                     ENDIF
 170              CONTINUE
               ELSE
                  DO 171 K = IP(J),IP(J+1)-1
                     A(K) = ONE
 171              CONTINUE
               ENDIF
 176        CONTINUE
         ENDIF
         DW(1) = max(ZERO,CNTL(1))
         RINF3 = RINF3+ONE
         DW(2) = RINF3
         IW(1) = JOB
         IF (JOB.EQ.5) THEN
            CALL SMUMPS_454(M,N,NE,IP,IRN,A,PERM,NUM,
     &           IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1),
     &           DW(1),DW(M+1),RINF2)
         ENDIF
         IF (JOB.EQ.6) THEN
            CALL SMUMPS_454(M,N,NE,IP,IW(3*N+2*M+1),A,PERM,NUM,
     &           IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1),
     &           DW(1),DW(M+1),RINF2)
         ENDIF
         IF (JOB.EQ.6) THEN
            DO 79 I = 1,M
               IF (DW(2*M+N+I).NE.0) THEN
                  DW(I) = DW(I) + log(DW(2*M+N+I))
               ENDIF
 79         CONTINUE
         ENDIF
         IF (NUM.EQ.N) THEN
            DO 80 J = 1,N
               IF (DW(2*M+J).NE.ZERO) THEN
                  DW(M+J) = DW(M+J) - log(DW(2*M+J))
               ELSE
                  DW(M+J) = ZERO
               ENDIF
 80         CONTINUE
         ENDIF
         FACT = 0.5*log(RINF2)
         DO 86 I = 1,M
            IF (DW(I).LT.FACT) GO TO 86
            WARN2 = 2
            GO TO 90
 86      CONTINUE 
         DO 87 J = 1,N
            IF (DW(M+J).LT.FACT) GO TO 87
            WARN2 = 2
            GO TO 90
 87      CONTINUE 
      ENDIF
 90   IF (NUM.LT.N) WARN1 = 1
      IF (JOB.EQ.4 .OR. JOB.EQ.5 .OR. JOB.EQ.6) THEN 
         IF (CNTL(1).LT.ZERO) WARN4 = 4
      ENDIF
      IF (INFO(1).EQ.0) THEN
         INFO(1) = WARN1 + WARN2 + WARN4
         IF (INFO(1).GT.0 .AND. ICNTL(2).GT.0) THEN
            WRITE(ICNTL(2),9010) INFO(1)
            IF (WARN1.EQ.1) WRITE(ICNTL(2),9011)
            IF (WARN2.EQ.2) WRITE(ICNTL(2),9012)
            IF (WARN4.EQ.4) WRITE(ICNTL(2),9014)
         ENDIF
      ENDIF
      IF (ICNTL(3).GE.0) THEN
         IF (ICNTL(4).EQ.0 .OR. ICNTL(4).EQ.1) THEN
            WRITE(ICNTL(3),9030) (INFO(J),J=1,2)
            WRITE(ICNTL(3),9031) NUM
            IF (ICNTL(4).EQ.0) THEN
               WRITE(ICNTL(3),9032) (PERM(J),J=1,min(10,M))
               IF (JOB.EQ.5 .OR. JOB.EQ.6) THEN
                  WRITE(ICNTL(3),9033) (DW(J),J=1,min(10,M))
                  WRITE(ICNTL(3),9034) (DW(M+J),J=1,min(10,N))
               ENDIF
            ELSEIF (ICNTL(4).EQ.1) THEN
               WRITE(ICNTL(3),9032) (PERM(J),J=1,M)
               IF (JOB.EQ.5 .OR. JOB.EQ.6) THEN
                  WRITE(ICNTL(3),9033) (DW(J),J=1,M)
                  WRITE(ICNTL(3),9034) (DW(M+J),J=1,N)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
 99   RETURN
 9001 FORMAT (' ****** Error in SMUMPS_443. INFO(1) = ',I2,
     &     ' because ',(A),' = ',I10)
 9004 FORMAT (' ****** Error in SMUMPS_443. INFO(1) = ',I2/
     &     '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in SMUMPS_443. INFO(1) = ',I2/
     &     '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in SMUMPS_443. INFO(1) = ',I2/
     &     '        Column ',I8,
     &     ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in SMUMPS_443. INFO(1) = ',I2/
     &     '        Column ',I8,
     &     ' contains two or more entries with row index ',I8)
 9008 FORMAT (' ****** Error in SMUMPS_443. INFO(1) = ',I2/
     &     '        CNTL(2) is negative or zero.')
 9010 FORMAT (' ****** Warning from SMUMPS_443. INFO(1) = ',I2)
 9011 FORMAT ('        - The matrix is structurally singular.')
 9012 FORMAT ('        - Some scaling factors may be too large.')
 9014 FORMAT ('        - CNTL(1) is negative and was treated as zero.')
 9020 FORMAT (' ****** Input parameters for SMUMPS_443:'/
     &     ' JOB =',I10/' M   =',I10/' N   =',I10/' NE  =',I10)
 9021 FORMAT (' IP(1:N+1)   = ',8I8/(15X,8I8))
 9022 FORMAT (' IRN(1:NE)   = ',8I8/(15X,8I8))
 9023 FORMAT (' A(1:NE)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9024 FORMAT (' ICNTL(1:10) = ',8I8/(15X,2I8))
 9025 FORMAT (' CNTL(1:10)  = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for SMUMPS_443:'/
     &     ' INFO(1:2)   = ',2I8)
 9031 FORMAT (' NUM         = ',I8)
 9032 FORMAT (' PERM(1:M)   = ',8I8/(15X,8I8))
 9033 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
 9034 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END SUBROUTINE SMUMPS_559
      SUBROUTINE SMUMPS_563(N,NZ,IP,IRN,A,FLAG,POSI)
      IMPLICIT NONE
      INTEGER N,NZ
      INTEGER IP(N+1),IRN(NZ)
      REAL A(NZ)
      INTEGER WR_POS,BEG_COL,ROW,COL,K,SV_POS
      INTEGER FLAG(N), POSI(N)
      INTEGER ierr
      FLAG = 0
      WR_POS = 1
      DO COL=1,N
         BEG_COL = WR_POS
         DO K=IP(COL),IP(COL+1)-1
            ROW = IRN(K)
            IF(FLAG(ROW) .NE. COL) THEN
               IRN(WR_POS) = ROW
               A(WR_POS) = A(K)
               FLAG(ROW) = COL
               POSI(ROW) = WR_POS
               WR_POS = WR_POS+1
            ELSE
               SV_POS = POSI(ROW)
               A(SV_POS) = A(SV_POS) + A(K)
            ENDIF
         ENDDO
         IP(COL) = BEG_COL
      ENDDO
      IP(N+1) = WR_POS
      NZ = WR_POS-1
      RETURN
      END SUBROUTINE SMUMPS_563
      SUBROUTINE SMUMPS_562(N,NZ,IP,IRN,FLAG,POSI)
      IMPLICIT NONE
      INTEGER N,NZ
      INTEGER IP(N+1),IRN(NZ)
      INTEGER WR_POS,BEG_COL,ROW,COL,K,SV_POS
      INTEGER FLAG(N), POSI(N)
      INTEGER ierr
      FLAG = 0
      WR_POS = 1
      DO COL=1,N
         BEG_COL = WR_POS
         DO K=IP(COL),IP(COL+1)-1
            ROW = IRN(K)
            IF(FLAG(ROW) .NE. COL) THEN
               IRN(WR_POS) = ROW
               FLAG(ROW) = COL
               POSI(ROW) = WR_POS
               WR_POS = WR_POS+1
            ENDIF
         ENDDO
         IP(COL) = BEG_COL
      ENDDO
      IP(N+1) = WR_POS
      NZ = WR_POS-1
      RETURN
      END SUBROUTINE SMUMPS_562
      SUBROUTINE SMUMPS_181( N, POOL, NSTK,
     &          PERM, FILS, FRERE, WORK )
      IMPLICIT NONE
      INTEGER, INTENT(IN)  ::  N
      INTEGER, INTENT(IN)  ::  FILS( N ), FRERE( N )
      INTEGER, INTENT(IN), TARGET :: POOL(N), NSTK(N)
      INTEGER, INTENT(OUT), TARGET :: WORK(N,2)
      INTEGER, INTENT(OUT) ::  PERM( N )
      INTEGER  :: IPERM, INODE, IN, I
      INTEGER  :: INBLEAF, INBROOT
      INTEGER, POINTER  :: IPOOL(:), INSTK(:)
      IPOOL => WORK(1:N,1)
      INSTK => WORK(1:N,2)
      DO I=1, N
         IPOOL(I) = POOL(I)
         INSTK(I) = NSTK(I)
      END DO
      IF (N.GT.1) THEN
         IF(POOL(N-1) .LT. 0) THEN
            INBLEAF = N-1
            INBROOT = IPOOL(N)
            IPOOL(N-1) = -IPOOL(N-1)-1
         ELSE IF ((POOL(N-1) .GT. 0) .AND.
     $           (POOL(N) .LT. 0)) THEN
            INBLEAF = N
            INBROOT = N
            IPOOL(N) = -IPOOL(N)-1
         ELSE
            INBLEAF = POOL(N-1)
            INBROOT = POOL(N)
         END IF
      END IF
      IPERM = 1
      DO WHILE ( INBLEAF .NE. 0 )
        INODE = IPOOL( INBLEAF )
        INBLEAF = INBLEAF - 1
        IN = INODE
        DO WHILE ( IN .GT. 0 )
          PERM ( IN ) = IPERM
          IPERM = IPERM + 1
          IN = FILS( IN )
        END DO
        IN = FRERE( INODE )
        DO WHILE( IN .GT. 0 )
          IN = FRERE( IN )
        END DO
        IF ( IN .eq. 0 ) THEN
          INBROOT = INBROOT - 1
        ELSE
          IN = - IN
          INSTK( IN ) = INSTK( IN ) - 1
          IF ( INSTK( IN ) .eq. 0 ) THEN
            INBLEAF = INBLEAF + 1
            IPOOL( INBLEAF ) = IN
          END IF
        END IF
      END DO
      RETURN
      END SUBROUTINE SMUMPS_181
      SUBROUTINE SMUMPS_448(ICNTL,CNTL)
      IMPLICIT NONE
      INTEGER NICNTL, NCNTL
      PARAMETER (NICNTL=10, NCNTL=10)
      INTEGER ICNTL(NICNTL)
      REAL CNTL(NCNTL)
      INTEGER I
      ICNTL(1) =  6
      ICNTL(2) =  6
      ICNTL(3) = -1
      ICNTL(4) = -1
      ICNTL(5) =  0
      DO 10 I = 6,NICNTL
        ICNTL(I) = 0
   10 CONTINUE
      CNTL(1) = 0.0E0
      CNTL(2) = 0.0E0
      DO 20 I = 3,NCNTL
        CNTL(I) = 0.0E0
   20 CONTINUE
      RETURN
      END
      SUBROUTINE SMUMPS_444
     &           (M,N,NE,IP,IRN,A,IPERM,NUM,JPERM,PR,Q,L,D,RINF) 
      IMPLICIT NONE
      INTEGER M,N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(M),JPERM(N),PR(N),Q(M),L(M)
      REAL A(NE)
      REAL D(M), RINF
      INTEGER I,II,J,JJ,JORD,Q0,QLEN,IDUM,JDUM,ISP,JSP,
     &        K,KK,KK1,KK2,I0,UP,LOW
      REAL    CSP,DI,DNEW,DQ0,AI,A0,BV,TBV,RLX
      REAL    ZERO,MINONE
      PARAMETER (ZERO=0.0D+0,MINONE=-1.0D+0)
      INTRINSIC abs,min
      EXTERNAL SMUMPS_445, SMUMPS_446, SMUMPS_447, SMUMPS_455
      RLX = D(1)
      NUM = 0
      BV = RINF
      DO 10 K = 1,N
        JPERM(K) = 0
        PR(K) = IP(K)
   10 CONTINUE
      DO 12 K = 1,M
        IPERM(K) = 0
        D(K) = ZERO
   12 CONTINUE
      DO 30 J = 1,N
        A0 = MINONE
        DO 20 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = abs(A(K))
          IF (AI.GT.D(I)) D(I) = AI
          IF (JPERM(J).NE.0) GO TO 20
          IF (AI.GE.BV) THEN
            A0 = BV
            IF (IPERM(I).NE.0) GO TO 20
            JPERM(J) = I 
            IPERM(I) = J
            NUM = NUM + 1
          ELSE
            IF (AI.LE.A0) GO TO 20
            A0 = AI
            I0 = I
          ENDIF
   20   CONTINUE
        IF (A0.NE.MINONE .AND. A0.LT.BV) THEN
          BV = A0
          IF (IPERM(I0).NE.0) GO TO 30
          IPERM(I0) = J
          JPERM(J) = I0
          NUM = NUM + 1
        ENDIF
   30 CONTINUE
      IF (M.EQ.N) THEN
        DO 35 I = 1,M
          BV = min(BV,D(I))
   35   CONTINUE
      ENDIF
      IF (NUM.EQ.N) GO TO 1000
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        DO 50 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = abs(A(K))
          IF (AI.LT.BV) GO TO 50
          IF (IPERM(I).EQ.0) GO TO 90
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 50
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).NE.0) GO TO 70
            IF (abs(A(KK)).GE.BV) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   50   CONTINUE
        GO TO 95
   80   JPERM(JJ) = II
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = I
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 99 I = 1,M
        D(I) = MINONE
        L(I) = 0
   99 CONTINUE
      TBV = BV * (1-RLX)
      DO 100 JORD = 1,N
        IF (JPERM(JORD).NE.0) GO TO 100
        QLEN = 0
        LOW = M + 1
        UP = M + 1
        CSP = MINONE
        J = JORD
        PR(J) = -1
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = abs(A(K))
          IF (CSP.GE.DNEW) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = I
            JSP = J
            IF (CSP.GE.TBV) GO TO 160
          ELSE
            D(I) = DNEW
            IF (DNEW.GE.TBV) THEN
              LOW = LOW - 1
              Q(LOW) = I
            ELSE
              QLEN = QLEN + 1
              L(I) = QLEN
              CALL SMUMPS_445(I,M,Q,D,L,1)
            ENDIF
            JJ = IPERM(I)
            PR(JJ) = J
          ENDIF
  115   CONTINUE
        DO 150 JDUM = 1,NUM
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (CSP.GE.D(I)) GO TO 160
            BV = D(I)
            TBV = BV * (1-RLX)
            DO 152 IDUM = 1,M
              CALL SMUMPS_446(QLEN,M,Q,D,L,1)
              L(I) = 0
              LOW = LOW - 1
              Q(LOW) = I
              IF (QLEN.EQ.0) GO TO 153
              I = Q(1)
              IF (D(I).LT.TBV) GO TO 153
  152       CONTINUE
          ENDIF
  153     UP = UP - 1
          Q0 = Q(UP)
          DQ0 = D(Q0)
          L(Q0) = UP
          J = IPERM(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
            IF (L(I).GE.UP) GO TO 155
            DNEW = min(DQ0,abs(A(K)))
            IF (CSP.GE.DNEW) GO TO 155
            IF (IPERM(I).EQ.0) THEN
              CSP = DNEW
              ISP = I
              JSP = J
              IF (CSP.GE.TBV) GO TO 160
            ELSE
              DI = D(I)
              IF (DI.GE.TBV .OR. DI.GE.DNEW) GO TO 155
              D(I) = DNEW
              IF (DNEW.GE.TBV) THEN
                IF (DI.NE.MINONE) THEN
                  CALL SMUMPS_447(L(I),QLEN,M,Q,D,L,1)
                ENDIF
                L(I) = 0
                LOW = LOW - 1
                Q(LOW) = I
              ELSE
                IF (DI.EQ.MINONE) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL SMUMPS_445(I,M,Q,D,L,1)
              ENDIF
              JJ = IPERM(I)
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE
  160   IF (CSP.EQ.MINONE) GO TO 190
        BV = min(BV,CSP)
        TBV = BV * (1-RLX)
        NUM = NUM + 1
        I = ISP
        J = JSP
        DO 170 JDUM = 1,NUM+1
          I0 = JPERM(J)
          JPERM(J) = I
          IPERM(I) = J
          J = PR(J)
          IF (J.EQ.-1) GO TO 190
          I = I0
  170   CONTINUE
  190   DO 191 KK = UP,M
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  191   CONTINUE 
        DO 192 KK = LOW,UP-1
          I = Q(KK)
          D(I) = MINONE
  192   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  193   CONTINUE
  100 CONTINUE
 1000 IF (M.EQ.N .and. NUM.EQ.N) GO TO 2000
      CALL SMUMPS_455(M,N,IPERM,L,JPERM)
 2000 RETURN
      END
      SUBROUTINE SMUMPS_445(I,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER I,N,IWAY
      INTEGER Q(N),L(N)
      REAL D(N)
      INTEGER IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      REAL DI
      POS = L(I)
      IF (POS.LE.1) GO TO 20
      DI = D(I)
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20 
          Q(POS) = QK
          L(QK) = POS 
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
      ELSE
        DO 15 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   15   CONTINUE
      ENDIF
   20 Q(POS) = I
      L(I) = POS
      RETURN
      END
      SUBROUTINE SMUMPS_446(QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER QLEN,N,IWAY
      INTEGER Q(N),L(N)
      REAL D(N)
      INTEGER I,IDUM,K,POS,POSK
      PARAMETER (K=2)
      REAL DK,DR,DI
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = 1
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 20
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   10   CONTINUE
      ELSE
        DO 15 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 20
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   15   CONTINUE
      ENDIF
   20 Q(POS) = I
      L(I) = POS
      RETURN
      END
      SUBROUTINE SMUMPS_447(POS0,QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER POS0,QLEN,N,IWAY
      INTEGER Q(N),L(N)
      REAL D(N)
      INTEGER I,IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      REAL DK,DR,DI
      IF (QLEN.EQ.POS0) THEN
        QLEN = QLEN - 1
        RETURN
      ENDIF
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = POS0
      IF (IWAY.EQ.1) THEN
        IF (POS.LE.1) GO TO 20
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20 
          Q(POS) = QK
          L(QK) = POS 
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
   20   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 30 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   30   CONTINUE
      ELSE
        IF (POS.LE.1) GO TO 34
        DO 32 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 34 
          Q(POS) = QK
          L(QK) = POS 
          POS = POSK
          IF (POS.LE.1) GO TO 34
   32   CONTINUE
   34   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 36 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   36   CONTINUE
      ENDIF
   40 Q(POS) = I
      L(I) = POS
      RETURN
      END
      SUBROUTINE SMUMPS_450(IP,LENL,LENH,W,WLEN,A,NVAL,VAL)
      IMPLICIT NONE
      INTEGER WLEN,NVAL
      INTEGER IP(*),LENL(*),LENH(*),W(*)
      REAL A(*),VAL
      INTEGER XX,J,K,II,S,POS
      PARAMETER (XX=10)
      REAL SPLIT(XX),HA
      NVAL = 0 
      DO 10 K = 1,WLEN
        J = W(K)
        DO 15 II = IP(J)+LENL(J),IP(J)+LENH(J)-1
          HA = A(II)
          IF (NVAL.EQ.0) THEN
            SPLIT(1) = HA
            NVAL = 1
          ELSE
            DO 20 S = NVAL,1,-1
              IF (SPLIT(S).EQ.HA) GO TO 15
              IF (SPLIT(S).GT.HA) THEN
                POS = S + 1
                GO TO 21
              ENDIF
  20        CONTINUE
            POS = 1
  21        DO 22 S = NVAL,POS,-1
              SPLIT(S+1) = SPLIT(S)
  22        CONTINUE
            SPLIT(POS) = HA
            NVAL = NVAL + 1
          ENDIF
          IF (NVAL.EQ.XX) GO TO 11
  15    CONTINUE
  10  CONTINUE
  11  IF (NVAL.GT.0) VAL = SPLIT((NVAL+1)/2)
      RETURN
      END
      SUBROUTINE SMUMPS_451(N,NE,IP,IRN,A)
      IMPLICIT NONE
      INTEGER N,NE
      INTEGER IP(N+1),IRN(NE)
      REAL A(NE)
      INTEGER THRESH,TDLEN
      PARAMETER (THRESH=15,TDLEN=50)
      INTEGER J,IPJ,K,LEN,R,S,HI,FIRST,MID,LAST,TD
      REAL HA,KEY
      INTEGER TODO(TDLEN)
      DO 100 J = 1,N
        LEN = IP(J+1) - IP(J)
        IF (LEN.LE.1) GO TO 100
        IPJ = IP(J)
        IF (LEN.LT.THRESH) GO TO 400
        TODO(1) = IPJ
        TODO(2) = IPJ + LEN
        TD = 2
  500   CONTINUE
        FIRST = TODO(TD-1)
        LAST = TODO(TD)
        KEY = A((FIRST+LAST)/2)
        DO 475 K = FIRST,LAST-1
          HA = A(K)
          IF (HA.EQ.KEY) GO TO 475
          IF (HA.GT.KEY) GO TO 470
          KEY = HA
          GO TO 470
  475   CONTINUE
        TD = TD - 2
        GO TO 425
  470   MID = FIRST
        DO 450 K = FIRST,LAST-1
          IF (A(K).LE.KEY) GO TO 450
          HA = A(MID)
          A(MID) = A(K)
          A(K) = HA
          HI = IRN(MID)
          IRN(MID) = IRN(K)
          IRN(K) = HI
          MID = MID + 1
  450   CONTINUE
        IF (MID-FIRST.GE.LAST-MID) THEN
          TODO(TD+2) = LAST
          TODO(TD+1) = MID
          TODO(TD) = MID
        ELSE
          TODO(TD+2) = MID
          TODO(TD+1) = FIRST
          TODO(TD) = LAST
          TODO(TD-1) = MID
        ENDIF
        TD = TD + 2
  425   CONTINUE
        IF (TD.EQ.0) GO TO 400 
        IF (TODO(TD)-TODO(TD-1).GE.THRESH) GO TO 500
        TD = TD - 2
        GO TO 425
  400   DO 200 R = IPJ+1,IPJ+LEN-1
          IF (A(R-1) .LT. A(R)) THEN
            HA = A(R)
            HI = IRN(R)
            A(R) = A(R-1)
            IRN(R) = IRN(R-1)
            DO 300 S = R-1,IPJ+1,-1
              IF (A(S-1) .LT. HA) THEN
                A(S) = A(S-1)
                IRN(S) = IRN(S-1)
              ELSE
                A(S) = HA
                IRN(S) = HI
                GO TO 200 
              END IF
  300       CONTINUE
            A(IPJ) = HA
            IRN(IPJ) = HI
          END IF
  200   CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE SMUMPS_452(M,N,NE,IP,IRN,A,IPERM,NUMX,
     &           W,LEN,LENL,LENH,FC,IW,IW4,RLX,RINF)
      IMPLICIT NONE
      INTEGER M,N,NE,NUMX
      INTEGER IP(N+1),IRN(NE),IPERM(N), 
     &        W(N),LEN(N),LENL(N),LENH(N),FC(N),IW(M),IW4(3*N+M)
      REAL A(NE),RLX,RINF
      INTEGER NUM,NVAL,WLEN,II,I,J,K,L,CNT,MOD,IDUM1,IDUM2,IDUM3
      REAL    BVAL,BMIN,BMAX
      EXTERNAL SMUMPS_450,SMUMPS_453,SMUMPS_455
      DO 20 J = 1,N
        FC(J) = J
        LEN(J) = IP(J+1) - IP(J)
   20 CONTINUE
      DO 21 I = 1,M
        IW(I) = 0
   21 CONTINUE
      CNT = 1
      MOD = 1
      NUMX = 0
      CALL SMUMPS_453(CNT,MOD,M,N,IRN,NE,IP,LEN,FC,IW,NUMX,N,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(2*N+M+1))
      NUM = NUMX
      IF (NUM.NE.N) THEN
        BMAX = RINF
      ELSE
        BMAX = RINF
        DO 30 J = 1,N
          BVAL = 0.0
          DO 25 K = IP(J),IP(J+1)-1
            IF (A(K).GT.BVAL) BVAL = A(K)
   25     CONTINUE
          IF (BVAL.LT.BMAX) BMAX = BVAL
   30   CONTINUE
        BMAX = 1.001 * BMAX
      ENDIF
      BVAL = 0.0
      BMIN = 0.0
      WLEN = 0
      DO 48 J = 1,N
        L = IP(J+1) - IP(J)
        LENH(J) = L
        LEN(J) = L
        DO 45 K = IP(J),IP(J+1)-1
          IF (A(K).LT.BMAX) GO TO 46
   45   CONTINUE
        K = IP(J+1)
   46   LENL(J) = K - IP(J)
        IF (LENL(J).EQ.L) GO TO 48
        WLEN = WLEN + 1
        W(WLEN) = J
   48 CONTINUE
      DO 90 IDUM1 = 1,NE
        IF (NUM.EQ.NUMX) THEN
          DO 50 I = 1,M
            IPERM(I) = IW(I)
   50     CONTINUE
          DO 80 IDUM2 = 1,NE
            BMIN = BVAL
            IF (BMAX-BMIN .LE. RLX) GO TO 1000
            CALL SMUMPS_450(IP,LENL,LEN,W,WLEN,A,NVAL,BVAL)
            IF (NVAL.LE.1) GO TO 1000
            K = 1
            DO 70 IDUM3 = 1,N
              IF (K.GT.WLEN) GO TO 71
              J = W(K)
              DO 55 II = IP(J)+LEN(J)-1,IP(J)+LENL(J),-1
                IF (A(II).GE.BVAL) GO TO 60 
                I = IRN(II)
                IF (IW(I).NE.J) GO TO 55
                IW(I) = 0
                NUM = NUM - 1
                FC(N-NUM) = J
   55         CONTINUE
   60         LENH(J) = LEN(J)
              LEN(J) = II - IP(J) + 1
              IF (LENL(J).EQ.LENH(J)) THEN
                W(K) = W(WLEN)
                WLEN = WLEN - 1
              ELSE
                K = K + 1
              ENDIF
   70       CONTINUE
   71       IF (NUM.LT.NUMX) GO TO 81
   80     CONTINUE
   81     MOD = 1
        ELSE
          BMAX = BVAL
          IF (BMAX-BMIN .LE. RLX) GO TO 1000
          CALL SMUMPS_450(IP,LEN,LENH,W,WLEN,A,NVAL,BVAL)
          IF (NVAL.EQ.0. OR. BVAL.EQ.BMIN) GO TO 1000
          K = 1
          DO 87 IDUM3 = 1,N
            IF (K.GT.WLEN) GO TO 88
            J = W(K)
            DO 85 II = IP(J)+LEN(J),IP(J)+LENH(J)-1
              IF (A(II).LT.BVAL) GO TO 86
   85       CONTINUE
   86       LENL(J) = LEN(J)
            LEN(J) = II - IP(J)
            IF (LENL(J).EQ.LENH(J)) THEN
              W(K) = W(WLEN)
              WLEN = WLEN - 1
            ELSE
              K = K + 1
            ENDIF
   87     CONTINUE
   88     MOD = 0
        ENDIF
        CNT = CNT + 1
        CALL SMUMPS_453(CNT,MOD,M,N,IRN,NE,IP,LEN,FC,IW,NUM,NUMX,
     &              IW4(1),IW4(N+1),IW4(2*N+1),IW4(2*N+M+1))
   90 CONTINUE 
 1000 IF (M.EQ.N .and. NUMX.EQ.N) GO TO 2000
      CALL SMUMPS_455(M,N,IPERM,IW,W)
 2000 RETURN
      END
      SUBROUTINE SMUMPS_453
     &           (ID,MOD,M,N,IRN,LIRN,IP,LENC,FC,IPERM,NUM,NUMX,
     &           PR,ARP,CV,OUT)
      IMPLICIT NONE
      INTEGER ID,MOD,M,N,LIRN,NUM,NUMX
      INTEGER ARP(N),CV(M),IRN(LIRN),IP(N),
     &        FC(N),IPERM(M),LENC(N),OUT(N),PR(N)
      INTEGER I,II,IN1,IN2,J,J1,JORD,K,KK,LAST,NFC,
     &        NUM0,NUM1,NUM2,ID0,ID1
      IF (ID.EQ.1) THEN
        DO 5 I = 1,M
          CV(I) = 0
    5   CONTINUE
        DO 6 J = 1,N
          ARP(J) = 0
    6   CONTINUE
        NUM1 = N
        NUM2 = N
      ELSE
        IF (MOD.EQ.1) THEN
          DO 8 J = 1,N
            ARP(J) = 0
    8     CONTINUE
        ENDIF
        NUM1 = NUMX
        NUM2 = N - NUMX
      ENDIF
      NUM0 = NUM
      NFC = 0
      ID0 = (ID-1)*N 
      DO 100 JORD = NUM0+1,N
        ID1 = ID0 + JORD
        J = FC(JORD-NUM0)
        PR(J) = -1
        DO 70 K = 1,JORD
          IF (ARP(J).GE.LENC(J)) GO TO 30
          IN1 = IP(J) + ARP(J)
          IN2 = IP(J) + LENC(J) - 1
          DO 20 II = IN1,IN2
            I = IRN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
          ARP(J) = LENC(J)
   30     OUT(J) = LENC(J) - 1
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENC(J) - 1
            IN1 = IN2 - IN1
            DO 40 II = IN1,IN2
              I = IRN(II)
              IF (CV(I).EQ.ID1) GO TO 40
              J1 = J
              J = IPERM(I)
              CV(I) = ID1
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
   50       J1 = PR(J)
            IF (J1.EQ.-1) THEN
              NFC = NFC + 1
              FC(NFC) = J
              IF (NFC.GT.NUM2) THEN
                LAST = JORD
                GO TO 101
              ENDIF
              GO TO 100
            ENDIF
            J = J1
   60     CONTINUE
   70   CONTINUE
   80   IPERM(I) = J
        ARP(J) = II - IP(J) + 1
        NUM = NUM + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 95
          II = IP(J) + LENC(J) - OUT(J) - 2
          I = IRN(II)
          IPERM(I) = J
   90   CONTINUE
   95   IF (NUM.EQ.NUM1) THEN
          LAST = JORD
          GO TO 101
        ENDIF
  100 CONTINUE
      LAST = N
  101 DO 110 JORD = LAST+1,N
        NFC = NFC + 1
        FC(NFC) = FC(JORD-NUM0)
  110 CONTINUE
      RETURN
      END
      SUBROUTINE SMUMPS_454(M,N,NE,IP,IRN,A,IPERM,NUM,
     &           JPERM,OUT,PR,Q,L,U,D,RINF) 
      IMPLICIT NONE
      INTEGER M,N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(M),JPERM(N),OUT(N),PR(N),Q(M),L(M)
      REAL A(NE),U(M),D(M),RINF,RINF3
      INTEGER I,I0,II,J,JJ,JORD,Q0,QLEN,JDUM,ISP,JSP,
     &        K,K0,K1,K2,KK,KK1,KK2,UP,LOW
      REAL    CSP,DI,DMIN,DNEW,DQ0,VJ,RLX
      LOGICAL LORD
      REAL    ZERO
      PARAMETER (ZERO=0.0D+0)
      EXTERNAL SMUMPS_445, SMUMPS_446, SMUMPS_447, SMUMPS_455
      RLX = U(1)
      RINF3 = U(2)
      LORD = (JPERM(1).EQ.6)
      NUM = 0
      DO 10 K = 1,N
        JPERM(K) = 0
        PR(K) = IP(K)
        D(K) = RINF
   10 CONTINUE
      DO 15 K = 1,M
        U(K) = RINF3
        IPERM(K) = 0
        L(K) = 0
   15 CONTINUE
      DO 30 J = 1,N
         IF (IP(J+1)-IP(J) .GT. N/10 .AND. N.GT.50) GO TO 30
        DO 20 K = IP(J),IP(J+1)-1
          I = IRN(K)
          IF (A(K).GT.U(I)) GO TO 20
          U(I) = A(K)
          IPERM(I) = J
          L(I) = K
   20   CONTINUE
   30 CONTINUE
      DO 40 I = 1,M
        J = IPERM(I)
        IF (J.EQ.0) GO TO 40
        IF (JPERM(J).EQ.0) THEN
          JPERM(J) = L(I)
          D(J) = U(I)
          NUM = NUM + 1
        ELSEIF (D(J).GT.U(I)) THEN
          K = JPERM(J)
          II = IRN(K)
          IPERM(II) = 0
          JPERM(J) = L(I)
          D(J) = U(I)
        ELSE
          IPERM(I) = 0
        ENDIF
   40 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 45 K = 1,M
        D(K) = ZERO
   45 CONTINUE
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        K1 = IP(J)
        K2 = IP(J+1) - 1
        IF (K1.GT.K2) GO TO 95
        VJ = RINF
        DO 50 K = K1,K2
          I = IRN(K)
          DI = A(K) - U(I)
          IF (DI.GT.VJ) GO TO 50
          IF (DI.LT.VJ .OR. DI.EQ.RINF) GO TO 55
          IF (IPERM(I).NE.0 .OR. IPERM(I0).EQ.0) GO TO 50
   55     VJ = DI
          I0 = I
          K0 = K
   50   CONTINUE
        D(J) = VJ
        K = K0
        I = I0
        IF (IPERM(I).EQ.0) GO TO 90
        DO 60 K = K0,K2
          I = IRN(K)
          IF (A(K)-U(I).GT.VJ) GO TO 60 
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 60
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).GT.0) GO TO 70
            IF (A(KK)-U(II).LE.D(JJ)) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   60   CONTINUE
        GO TO 95
   80   JPERM(JJ) = KK
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = K
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 99 I = 1,M
        D(I) = RINF
        L(I) = 0
   99 CONTINUE
      DO 100 JORD = 1,N
        IF (JPERM(JORD).NE.0) GO TO 100
        DMIN = RINF
        QLEN = 0
        LOW = M + 1
        UP = M + 1
        CSP = RINF
        J = JORD
        PR(J) = -1
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = A(K) - U(I)
          IF (DNEW.GE.CSP) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = K
            JSP = J
          ELSE
            IF (DNEW.LT.DMIN) DMIN = DNEW
            D(I) = DNEW
            QLEN = QLEN + 1
            Q(QLEN) = K
          ENDIF
  115   CONTINUE
        Q0 = QLEN
        QLEN = 0
        DO 120 KK = 1,Q0
          K = Q(KK)
          I = IRN(K)
          IF (CSP.LE.D(I)) THEN
            D(I) = RINF
            GO TO 120
          ENDIF
          IF (D(I).LE.DMIN) THEN
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
          ELSE
            QLEN = QLEN + 1
            L(I) = QLEN
            CALL SMUMPS_445(I,M,Q,D,L,2)
          ENDIF
          JJ = IPERM(I)
          OUT(JJ) = K
          PR(JJ) = J
  120   CONTINUE
        DO 150 JDUM = 1,NUM
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (D(I).LT.RINF) DMIN = D(I)*(1+RLX)
            IF (DMIN.GE.CSP) GO TO 160
  152       CALL SMUMPS_446(QLEN,M,Q,D,L,2)
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
            IF (QLEN.EQ.0) GO TO 153
            I = Q(1)
            IF (D(I).GT.DMIN) GO TO 153
            GO TO 152
          ENDIF
  153     Q0 = Q(UP-1)
          DQ0 = D(Q0)
          IF (DQ0.GE.CSP) GO TO 160
          IF (DMIN.GE.CSP) GO TO 160
          UP = UP - 1
          J = IPERM(Q0)
          VJ = DQ0 - A(JPERM(J)) + U(Q0)
          K1 = IP(J+1)-1
          IF (LORD) THEN
            IF (CSP.NE.RINF) THEN
              DI = CSP - VJ
              IF (A(K1).GE.DI) THEN
                K0 = JPERM(J)
                IF (K0.GE.K1-6) GO TO 178
  177           CONTINUE
                  K = (K0+K1)/2
                  IF (A(K).GE.DI) THEN 
                    K1 = K
                  ELSE 
                    K0 = K
                  ENDIF
                  IF (K0.GE.K1-6) GO TO 178
                GO TO 177
  178           DO 179 K = K0+1,K1          
                  IF (A(K).LT.DI) GO TO 179
                  K1 = K - 1
                  GO TO 181
  179           CONTINUE
              ENDIF
            ENDIF
  181       IF (K1.EQ.JPERM(J)) K1 = K1 - 1
          ENDIF
          K0 = IP(J)
          DI = CSP - VJ
          DO 155 K = K0,K1
            I = IRN(K)
            IF (L(I).GE.LOW) GO TO 155
            DNEW = A(K) - U(I)
            IF (DNEW.GE.DI) GO TO 155
            DNEW = DNEW + VJ
            IF (DNEW.GT.D(I)) GO TO 155
            IF (IPERM(I).EQ.0) THEN
              CSP = DNEW
              ISP = K
              JSP = J
              DI = CSP - VJ
            ELSE
              IF (DNEW.GE.D(I)) GO TO 155
              D(I) = DNEW
              IF (DNEW.LE.DMIN) THEN
                IF (L(I).NE.0) THEN
                  CALL SMUMPS_447(L(I),QLEN,M,Q,D,L,2)
                ENDIF
                LOW = LOW - 1
                Q(LOW) = I
                L(I) = LOW
              ELSE   
                IF (L(I).EQ.0) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL SMUMPS_445(I,M,Q,D,L,2)
              ENDIF
              JJ = IPERM(I)
              OUT(JJ) = K
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE
  160   IF (CSP.EQ.RINF) GO TO 190
        NUM = NUM + 1
        I = IRN(ISP)
        J = JSP
        IPERM(I) = J
        JPERM(J) = ISP
        DO 170 JDUM = 1,NUM
          JJ = PR(J) 
          IF (JJ.EQ.-1) GO TO 180
          K = OUT(J)
          I = IRN(K)
          IPERM(I) = JJ
          JPERM(JJ) = K
          J = JJ
  170   CONTINUE
  180   DO 182 KK = UP,M
          I = Q(KK)
          U(I) = U(I) + D(I) - CSP
  182   CONTINUE 
  190   DO 191 KK = UP,M
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  191   CONTINUE
        DO 192 KK = LOW,UP-1
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  192   CONTINUE 
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  193   CONTINUE
  100 CONTINUE
 1000 CONTINUE
      DO 1200 J = 1,N
        K = JPERM(J)
        IF (K.NE.0) THEN
          D(J) = A(K) - U(IRN(K))
        ELSE
          D(J) = ZERO
        ENDIF
 1200 CONTINUE
      DO 1201 I = 1,M
        IF (IPERM(I).EQ.0) U(I) = ZERO
 1201 CONTINUE
      IF (M.EQ.N .and. NUM.EQ.N) GO TO 2000
      CALL SMUMPS_455(M,N,IPERM,L,JPERM)
 2000 RETURN
      END
      SUBROUTINE SMUMPS_457
     &           (M,N,IRN,LIRN,IP,LENC,IPERM,NUM,PR,ARP,CV,OUT)
      IMPLICIT NONE
      INTEGER LIRN,M,N,NUM
      INTEGER ARP(N),CV(M),IRN(LIRN),IP(N),IPERM(M),LENC(N),OUT(N),PR(N)
      INTEGER I,II,IN1,IN2,J,J1,JORD,K,KK
      EXTERNAL SMUMPS_455
      DO 10 I = 1,M
        CV(I) = 0
        IPERM(I) = 0
   10 CONTINUE
      DO 12 J = 1,N
        ARP(J) = LENC(J) - 1
   12 CONTINUE
      NUM = 0
      DO 1000 JORD = 1,N
        J = JORD
        PR(J) = -1
        DO 70 K = 1,JORD
          IN1 = ARP(J)
          IF (IN1.LT.0) GO TO 30
          IN2 = IP(J) + LENC(J) - 1
          IN1 = IN2 - IN1
          DO 20 II = IN1,IN2
            I = IRN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
          ARP(J) = -1
   30     CONTINUE
          OUT(J) = LENC(J) - 1
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENC(J) - 1
            IN1 = IN2 - IN1
            DO 40 II = IN1,IN2
              I = IRN(II)
              IF (CV(I).EQ.JORD) GO TO 40
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
   50       CONTINUE
            J = PR(J)
            IF (J.EQ.-1) GO TO 1000
   60     CONTINUE
   70   CONTINUE
   80   CONTINUE
        IPERM(I) = J
        ARP(J) = IN2 - II - 1
        NUM = NUM + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 1000
          II = IP(J) + LENC(J) - OUT(J) - 2
          I = IRN(II)
          IPERM(I) = J
   90   CONTINUE
 1000 CONTINUE
      IF (M.EQ.N .and. NUM.EQ.N) GO TO 2000
      CALL SMUMPS_455(M,N,IPERM,CV,ARP)
 2000 RETURN
      END
      SUBROUTINE SMUMPS_455(M,N,IPERM,RW,CW)
      IMPLICIT NONE
      INTEGER M,N
      INTEGER RW(M),CW(N),IPERM(M)
      INTEGER I,J,K
      DO 10 J = 1,N
        CW(J) = 0
   10 CONTINUE
      K = 0
      DO 20 I = 1,M
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          RW(K) = I
        ELSE
          J = IPERM(I)
          CW(J) = I
        ENDIF
   20 CONTINUE
      K = 0
      DO 30 J = 1,N
        IF (CW(J).NE.0) GO TO 30
        K = K + 1
        I = RW(K)
        IPERM(I) = -J
   30 CONTINUE
      DO 40 J = N+1,M
        K = K + 1
        I = RW(K)
        IPERM(I) = -J
   40 CONTINUE
      RETURN
      END
