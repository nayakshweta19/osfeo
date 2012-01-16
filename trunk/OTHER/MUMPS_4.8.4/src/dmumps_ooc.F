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
      MODULE DMUMPS_OOC
      USE MUMPS_OOC_COMMON
      IMPLICIT NONE
      INTEGER NOT_IN_MEM,BEING_READ,NOT_USED,PERMUTED,USED,
     $     USED_NOT_PERMUTED,ALREADY_USED
      PARAMETER (NOT_IN_MEM=0,BEING_READ=-1,NOT_USED=-2,
     $     PERMUTED=-3,USED=-4,USED_NOT_PERMUTED=-5,ALREADY_USED=-6)
      INTEGER OOC_NOT_IN_MEM,OOC_PERMUTED,OOC_NOT_PERMUTED
      PARAMETER (OOC_NOT_IN_MEM=NOT_IN_MEM,
     $     OOC_PERMUTED=NOT_USED,OOC_NOT_PERMUTED=PERMUTED)
      INTEGER, DIMENSION(:,:),POINTER :: SIZE_OF_BLOCK
      INTEGER, DIMENSION(:),POINTER :: TOTAL_NB_OOC_NODES
      INTEGER :: OOC_SOLVE_TYPE_FCT
      INTEGER, DIMENSION(:),ALLOCATABLE :: IO_REQ
      INTEGER, DIMENSION(:),ALLOCATABLE :: LRLUS_SOLVE,LRLU_SOLVE_T,
     $     LRLU_SOLVE_B,POSFAC_SOLVE,IDEB_SOLVE_Z, SIZE_SOLVE_Z
     $     ,PDEB_SOLVE_Z
      INTEGER,SAVE :: FACT_AREA_SIZE,FIRST_ACTIVE_Z,
     $     SIZE_ZONE_SOLVE,SIZE_SOLVE_EMM,TMP_SIZE_FACT,TMP_NB_NODES,
     $     MAX_NB_NODES_FOR_ZONE,MAX_NB_REQ,CURRENT_SOLVE_READ_ZONE,
     $     CUR_POS_SEQUENCE,NB_Z,MIN_SIZE_READ,SOLVE_STEP,
     $     NB_ZONE_REQ,MAX_SIZE_FACTOR_OOC,MTYPE_OOC,NB_ACT
#if defined (NEW_PREF_SCHEME)
      INTEGER,SAVE :: MAX_PREF_SIZE
#endif
     $     ,NB_CALLED,REQ_ACT,NB_CALL
      INTEGER*8, SAVE :: OOC_VADDR_PTR
      INTEGER*8,SAVE :: SIZE_ZONE_REQ
      DOUBLE PRECISION,SAVE :: MAX_OOC_FILE_SIZE
      INTEGER,DIMENSION(:),ALLOCATABLE :: FIRST_POS_IN_READ,
     $     SIZE_OF_READ,READ_DEST,READ_MNG,REQ_TO_ZONE,POS_HOLE_T,
     $     POS_HOLE_B,REQ_ID,OOC_STATE_NODE
      INTEGER DMUMPS_ELEMENTARY_DATA_SIZE,N_OOC
      INTEGER, DIMENSION(:), ALLOCATABLE :: POS_IN_MEM, INODE_TO_POS
      INTEGER, DIMENSION(:), ALLOCATABLE :: CURRENT_POS_T,CURRENT_POS_B
      PUBLIC :: DMUMPS_575,DMUMPS_576,
     $     DMUMPS_577,
     $     DMUMPS_578,
     $     DMUMPS_579,
     $     DMUMPS_581,DMUMPS_582,
     $     DMUMPS_583,DMUMPS_584,
     $     DMUMPS_585,DMUMPS_586
         INTEGER, PARAMETER, PUBLIC :: TYPEF_BOTH_LU = -99976
         PUBLIC DMUMPS_688,
     &        DMUMPS_690,
     &        DMUMPS_OOC_PANEL_SIZE
         PRIVATE DMUMPS_695, 
     &        DMUMPS_697
      CONTAINS
      SUBROUTINE DMUMPS_711( STRAT_IO_ARG,
     & STRAT_IO_ASYNC_ARG, WITH_BUF_ARG, LOW_LEVEL_STRAT_IO_ARG )
      IMPLICIT NONE
      INTEGER, intent(out) :: LOW_LEVEL_STRAT_IO_ARG
      LOGICAL, intent(out) :: STRAT_IO_ASYNC_ARG, WITH_BUF_ARG
      INTEGER, intent(in)  :: STRAT_IO_ARG
      INTEGER TMP
      CALL MUMPS_OOC_IS_ASYNC_AVAIL(TMP)
      STRAT_IO_ASYNC_ARG=.FALSE.
      WITH_BUF_ARG=.FALSE.
      IF(TMP.EQ.1)THEN
         IF((STRAT_IO_ARG.EQ.1).OR.(STRAT_IO_ARG.EQ.2))THEN
            STRAT_IO_ASYNC=.TRUE.
            WITH_BUF=.FALSE.
         ELSEIF((STRAT_IO_ARG.EQ.4).OR.(STRAT_IO_ARG.EQ.5))THEN
            STRAT_IO_ASYNC_ARG=.TRUE.
            WITH_BUF_ARG=.TRUE.
         ELSEIF(STRAT_IO_ARG.EQ.3)THEN
            STRAT_IO_ASYNC_ARG=.FALSE.
            WITH_BUF_ARG=.TRUE.
         ENDIF
         LOW_LEVEL_STRAT_IO_ARG=mod(STRAT_IO_ARG,3)
      ELSE
         LOW_LEVEL_STRAT_IO_ARG=0
         IF(STRAT_IO_ARG.GE.3)THEN
            WITH_BUF_ARG=.TRUE.
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_711
      FUNCTION DMUMPS_579(INODE,ZONE)
      IMPLICIT NONE
      INTEGER INODE,ZONE
      LOGICAL DMUMPS_579
      DMUMPS_579=(LRLUS_SOLVE(ZONE).GE.
     $     SIZE_OF_BLOCK(STEP_OOC(INODE),OOC_FCT_TYPE))
      RETURN
      END FUNCTION DMUMPS_579
      SUBROUTINE DMUMPS_590(LA)
      IMPLICIT NONE
      INTEGER LA
      FACT_AREA_SIZE=LA
      END SUBROUTINE DMUMPS_590
      SUBROUTINE DMUMPS_575(id, MAXS)
      USE DMUMPS_STRUC_DEF
      USE DMUMPS_OOC_BUFFER
      IMPLICIT NONE
      INTEGER TMPDIR_MAX_LENGTH, PREFIX_MAX_LENGTH
      PARAMETER (TMPDIR_MAX_LENGTH=255, PREFIX_MAX_LENGTH=63)
      INTEGER, intent(in) :: MAXS  
      TYPE(DMUMPS_STRUC), TARGET :: id
      INTEGER IERR
      INTEGER allocok
      INTEGER ASYNC
      CHARACTER*1 TMP_DIR(TMPDIR_MAX_LENGTH),
     *            TMP_PREFIX(PREFIX_MAX_LENGTH)
      INTEGER DIM_DIR,DIM_PREFIX
      INTEGER, DIMENSION(:), ALLOCATABLE :: FILE_FLAG_TAB
      INTEGER TMP
      ICNTL1=id%ICNTL(1)
      MAX_SIZE_FACTOR_OOC=0
      OOC_FCT_TYPE=1
      N_OOC=id%N
      IS_INIT_CALLED=.TRUE.
      ASYNC=0
      SOLVE=.FALSE.
      IERR=0
      IF(allocated(IO_REQ))THEN
         DEALLOCATE(IO_REQ)
      ENDIF
      IF(associated(KEEP_OOC))THEN
         NULLIFY(KEEP_OOC)
      ENDIF
      IF(associated(STEP_OOC))THEN
         NULLIFY(STEP_OOC)
      ENDIF
      IF(associated(PROCNODE_OOC))THEN
         NULLIFY(PROCNODE_OOC)
      ENDIF
      IF(associated(OOC_INODE_SEQUENCE))THEN
         NULLIFY(OOC_INODE_SEQUENCE)
      ENDIF
      IF(associated(TOTAL_NB_OOC_NODES))THEN
         NULLIFY(TOTAL_NB_OOC_NODES)
      ENDIF
      IF(associated(SIZE_OF_BLOCK))THEN
         NULLIFY(SIZE_OF_BLOCK)
      ENDIF
      IF(associated(OOC_VADDR))THEN
         NULLIFY(OOC_VADDR)
      ENDIF
      IF(allocated(I_CUR_HBUF_NEXTPOS))THEN
         DEALLOCATE(I_CUR_HBUF_NEXTPOS)
      ENDIF
      CALL DMUMPS_588(id,IERR)
      IF(IERR.LT.0)THEN
         IF (ICNTL1 > 0)
     *   WRITE(ICNTL1,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
         id%INFO(1) = IERR
         id%INFO(2) = 0
         RETURN
      ENDIF
      STEP_OOC=>id%STEP
      PROCNODE_OOC=>id%PROCNODE_STEPS
      MYID_OOC=id%MYID
      SLAVEF_OOC=id%NSLAVES
      KEEP_OOC => id%KEEP
      SIZE_OF_BLOCK=>id%OOC_SIZE_OF_BLOCK
      OOC_VADDR=>id%OOC_VADDR
      IF(id%KEEP(107).GT.0)THEN
         SIZE_SOLVE_EMM=max(id%KEEP(202),int(dble(MAXS*
     $        0.9)*0.2))
         SIZE_ZONE_SOLVE=max(SIZE_SOLVE_EMM,
     $        int((dble(MAXS*0.9)-
     $        dble(SIZE_SOLVE_EMM))/dble(id%KEEP(107))))
         IF(SIZE_ZONE_SOLVE.EQ.SIZE_SOLVE_EMM)THEN
            SIZE_SOLVE_EMM=id%KEEP(202)
            SIZE_ZONE_SOLVE=int((dble(MAXS*0.9)-
     $           dble(SIZE_SOLVE_EMM))/dble(id%KEEP(107)))
         ENDIF
      ELSE
         SIZE_ZONE_SOLVE=int(dble(MAXS*0.9))
         SIZE_SOLVE_EMM=SIZE_ZONE_SOLVE
      ENDIF
      DMUMPS_ELEMENTARY_DATA_SIZE = id%KEEP(35)
      SIZE_OF_BLOCK=0
      ALLOCATE(id%OOC_NB_FILES(OOC_NB_FILE_TYPE), stat=allocok)
      IF (allocok .GT. 0) THEN
         IF (ICNTL1.GT.0)
     *   WRITE(ICNTL1,*) 'PB allocation in DMUMPS_INIT_OOC'
         id%INFO(1) = -13
         id%INFO(2) = OOC_NB_FILE_TYPE
         RETURN
      ENDIF
      id%OOC_NB_FILES=0
      OOC_VADDR_PTR=0
      CALL DMUMPS_711( id%KEEP(99), STRAT_IO_ASYNC,
     &                                WITH_BUF, LOW_LEVEL_STRAT_IO )
      TMP_SIZE_FACT=0
      TMP_NB_NODES=0
      MAX_NB_NODES_FOR_ZONE=0
      OOC_INODE_SEQUENCE=>id%OOC_INODE_SEQUENCE
      ALLOCATE(I_CUR_HBUF_NEXTPOS(OOC_NB_FILE_TYPE),
     $     stat=allocok)
      IF (allocok .GT. 0) THEN
         IF (ICNTL1.GT.0)
     *   WRITE(ICNTL1,*) 'PB allocation in DMUMPS_INIT_OOC'
         id%INFO(1) = -13
         id%INFO(2) = OOC_NB_FILE_TYPE
         RETURN
      ENDIF
      I_CUR_HBUF_NEXTPOS = 1
      IF(WITH_BUF)THEN
         CALL DMUMPS_669(id%INFO(1),id%INFO(2),IERR)
         IF(IERR.LT.0)THEN
            RETURN
         ENDIF
      ENDIF
      ALLOCATE(IO_REQ(id%KEEP(28)),
     $     stat=allocok)
      IF (allocok .GT. 0) THEN
         IF (ICNTL1.GT.0)
     *   WRITE(ICNTL1,*) 'PB allocation in DMUMPS_INIT_OOC'
         id%INFO(1) = -13
         id%INFO(2) = id%KEEP(28)
         RETURN
      ENDIF
      IF(STRAT_IO_ASYNC)THEN
         ASYNC=1
      ENDIF
      DIM_ERR_STR_OOC = ERR_STR_OOC_MAX_LEN
      CALL MUMPS_LOW_LEVEL_INIT_ERR_STR(DIM_ERR_STR_OOC,ERR_STR_OOC)
      DIM_DIR=len(trim(id%OOC_TMPDIR))
      DIM_PREFIX=len(trim(id%OOC_PREFIX))
      CALL DMUMPS_589(TMP_DIR(1),
     $     id%OOC_TMPDIR, TMPDIR_MAX_LENGTH, DIM_DIR )
      CALL DMUMPS_589(TMP_PREFIX(1),
     $     id%OOC_PREFIX, PREFIX_MAX_LENGTH, DIM_PREFIX)
      CALL MUMPS_LOW_LEVEL_INIT_PREFIX(DIM_PREFIX, TMP_PREFIX)
      CALL MUMPS_LOW_LEVEL_INIT_TMPDIR(DIM_DIR, TMP_DIR)
      ALLOCATE(FILE_FLAG_TAB(OOC_NB_FILE_TYPE),
     $     stat=allocok)
      IF (allocok .GT. 0) THEN
         IF (ICNTL1 .GT. 0)
     *   WRITE(ICNTL1,*) 'PB allocation in DMUMPS_INIT_OOC'
         id%INFO(1) = -13
         id%INFO(2) = OOC_NB_FILE_TYPE
         RETURN
      ENDIF
      FILE_FLAG_TAB(1:OOC_NB_FILE_TYPE)=0
      IERR=0
      TMP=int(id%KEEP8(11)/1000000_8,kind=4)+1
      IF((id%KEEP(201).EQ.1).AND.(id%KEEP(50).EQ.0))THEN
         TMP=MAX(1,TMP/2)
      ENDIF
      CALL MUMPS_LOW_LEVEL_INIT_OOC_C(MYID_OOC,TMP,
     $     id%KEEP(35),LOW_LEVEL_STRAT_IO,id%KEEP(211),OOC_NB_FILE_TYPE,
     $     FILE_FLAG_TAB,IERR)
      IF(IERR.LT.0)THEN
         IF (ICNTL1 .GT. 0 ) THEN
           WRITE(ICNTL1,*)MYID_OOC,': PB in MUMPS_LOW_LEVEL_INIT_OOC_C'
           WRITE(ICNTL1,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
         ENDIF
         id%INFO(1) = IERR
         id%INFO(2) = 0
         RETURN
      ENDIF
      CALL MUMPS_GET_MAX_FILE_SIZE_C(MAX_OOC_FILE_SIZE)
      DEALLOCATE(FILE_FLAG_TAB)
      RETURN
      END SUBROUTINE DMUMPS_575
      SUBROUTINE DMUMPS_576(INODE,PTRFAC,KEEP,KEEP8,
     *           A,LA,SIZE,IERR)
      USE DMUMPS_OOC_BUFFER
      IMPLICIT NONE
      INTEGER INODE,KEEP(500),LA
      INTEGER*8 KEEP8(150)
      INTEGER PTRFAC(KEEP(28))
      DOUBLE PRECISION A(LA)
      INTEGER SIZE,IERR,I,FLAG,NODE,ASYNC,REQUEST
      LOGICAL IO_C
      INTEGER ADDR_INT1,ADDR_INT2
      INTEGER TYPE
      TYPE=FCT
      IF(STRAT_IO_ASYNC)THEN
         ASYNC=1
      ELSE
         ASYNC=0
      ENDIF
      IERR=0
      IO_C=.TRUE.
      SIZE_OF_BLOCK(STEP_OOC(INODE),OOC_FCT_TYPE)=SIZE
      MAX_SIZE_FACTOR_OOC=max(MAX_SIZE_FACTOR_OOC,SIZE)
      OOC_VADDR(STEP_OOC(INODE),OOC_FCT_TYPE)=OOC_VADDR_PTR
      OOC_VADDR_PTR=OOC_VADDR_PTR+SIZE
      TMP_SIZE_FACT=TMP_SIZE_FACT+SIZE
      TMP_NB_NODES=TMP_NB_NODES+1
      IF(TMP_SIZE_FACT.GT.SIZE_ZONE_SOLVE)THEN
         MAX_NB_NODES_FOR_ZONE=max(MAX_NB_NODES_FOR_ZONE,TMP_NB_NODES)
         TMP_SIZE_FACT=0
         TMP_NB_NODES=0
      ENDIF
      IF (.NOT. WITH_BUF) THEN
         CALL MUMPS_677(ADDR_INT1,ADDR_INT2,
     $        OOC_VADDR(STEP_OOC(INODE),OOC_FCT_TYPE))
        CALL MUMPS_LOW_LEVEL_WRITE_OOC_C(LOW_LEVEL_STRAT_IO,
     &       A(PTRFAC(STEP_OOC(INODE))),SIZE,
     $       INODE,REQUEST,TYPE,ADDR_INT1,ADDR_INT2,IERR)
        IF(IERR.LT.0)THEN
           IF (ICNTL1.GT.0)
     *     WRITE(ICNTL1,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
           RETURN
        ENDIF
        IF(I_CUR_HBUF_NEXTPOS(OOC_FCT_TYPE).GT.KEEP_OOC(28))THEN
             WRITE(*,*)MYID_OOC,': Internal error (37) in OOC '
           CALL MUMPS_ABORT()
        ENDIF
        OOC_INODE_SEQUENCE(I_CUR_HBUF_NEXTPOS(OOC_FCT_TYPE),
     $       OOC_FCT_TYPE)=INODE
        I_CUR_HBUF_NEXTPOS(OOC_FCT_TYPE)=
     $       I_CUR_HBUF_NEXTPOS(OOC_FCT_TYPE)+1
      ELSE
         IF(SIZE.LE.HBUF_SIZE)THEN
            CALL DMUMPS_678
     *           (A(PTRFAC(STEP_OOC(INODE))),SIZE,IERR)
            OOC_INODE_SEQUENCE(I_CUR_HBUF_NEXTPOS(OOC_FCT_TYPE),
     $           OOC_FCT_TYPE) = INODE
            I_CUR_HBUF_NEXTPOS(OOC_FCT_TYPE) =
     $           I_CUR_HBUF_NEXTPOS(OOC_FCT_TYPE) + 1
#if ! defined (OOC_DEBUG)
            PTRFAC(STEP_OOC(INODE))=-777777
#endif
            RETURN
         ELSE
            CALL DMUMPS_707(OOC_FCT_TYPE,IERR)
            IF(IERR.LT.0)THEN
               RETURN
            ENDIF
            CALL DMUMPS_707(OOC_FCT_TYPE,IERR)
            IF(IERR.LT.0)THEN
               RETURN
            ENDIF
            CALL MUMPS_677(ADDR_INT1,ADDR_INT2,
     $           OOC_VADDR(STEP_OOC(INODE),OOC_FCT_TYPE))
            CALL MUMPS_LOW_LEVEL_WRITE_OOC_C(LOW_LEVEL_STRAT_IO,
     &           A(PTRFAC(STEP_OOC(INODE))),SIZE,
     $           INODE,REQUEST,TYPE,ADDR_INT1,ADDR_INT2,IERR)
            IF(IERR.LT.0)THEN
               IF (ICNTL1.GT.0)
     *         WRITE(*,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
               RETURN
            ENDIF
            IF(I_CUR_HBUF_NEXTPOS(OOC_FCT_TYPE).GT.KEEP_OOC(28))THEN
             WRITE(*,*)MYID_OOC,': Internal error (38) in OOC '
               CALL MUMPS_ABORT()
            ENDIF
            OOC_INODE_SEQUENCE(I_CUR_HBUF_NEXTPOS(OOC_FCT_TYPE),
     $           OOC_FCT_TYPE)=INODE
            I_CUR_HBUF_NEXTPOS(OOC_FCT_TYPE)=
     $           I_CUR_HBUF_NEXTPOS(OOC_FCT_TYPE)+1
            CALL DMUMPS_689(OOC_FCT_TYPE)
         ENDIF
      END IF
      NODE=-9999
#if ! defined (OOC_DEBUG)
      PTRFAC(STEP_OOC(INODE))=-777777
#endif
      IF(STRAT_IO_ASYNC)THEN
         IERR=0
         CALL MUMPS_WAIT_REQUEST(REQUEST,IERR)
         IF(IERR.LT.0)THEN
            IF (ICNTL1 .GT. 0)
     *      WRITE(ICNTL1,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
            RETURN
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_576
      SUBROUTINE DMUMPS_577(DEST,INODE,IERR
     $  )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR,INODE,STATUS
      DOUBLE PRECISION DEST
      INTEGER FLAG,NODE,ASYNC,REQUEST
      LOGICAL IO_C
      INTEGER ADDR_INT1,ADDR_INT2
      INTEGER TYPE
      TYPE=OOC_SOLVE_TYPE_FCT
      IF(STRAT_IO_ASYNC)THEN
        ASYNC=1
      ELSE
        ASYNC=0
      ENDIF
      IERR=0
      IO_C=.TRUE.
#if ! defined(OLD_READ)
      OOC_STATE_NODE(STEP_OOC(INODE))=NOT_USED
      CALL MUMPS_677(ADDR_INT1,ADDR_INT2,
     $     OOC_VADDR(STEP_OOC(INODE),OOC_FCT_TYPE))
      CALL MUMPS_LOW_LEVEL_DIRECT_READ(DEST,
     $     SIZE_OF_BLOCK(STEP_OOC(INODE),OOC_FCT_TYPE),
     $     TYPE,ADDR_INT1,ADDR_INT2,IERR)
      IF(IERR.LT.0)THEN
         IF (ICNTL1.GT.0) THEN
           WRITE(ICNTL1,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
           WRITE(ICNTL1,*)MYID_OOC,
     *     ': Problem in MUMPS_LOW_LEVEL_READ_OOC'
         ENDIF
         RETURN
      ENDIF
#else
      OOC_STATE_NODE(STEP_OOC(INODE))=NOT_USED
      CALL MUMPS_677(ADDR_INT1,ADDR_INT2,
     $     OOC_VADDR(STEP_OOC(INODE),OOC_FCT_TYPE))
      CALL MUMPS_LOW_LEVEL_READ_OOC_C(LOW_LEVEL_STRAT_IO,
     $     DEST,SIZE_OF_BLOCK(STEP_OOC(INODE),OOC_FCT_TYPE),
     $     INODE,REQUEST,TYPE,ADDR_INT1,ADDR_INT2,IERR)
      IF(IERR.LT.0)THEN
         IF (ICNTL1 .GT. 0 ) THEN
           WRITE(ICNTL1,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
           WRITE(ICNTL1,*)MYID_OOC,
     *     ': Problem in MUMPS_LOW_LEVEL_READ_OOC'
         ENDIF
         RETURN
      ENDIF
      IF(STRAT_IO_ASYNC)THEN
         IERR=0
         CALL MUMPS_WAIT_REQUEST(REQUEST,IERR)
         IF(IERR.LT.0)THEN
            IF (ICNTL1 .GT. 0 )
     *      WRITE(ICNTL1,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
            RETURN
         ENDIF
      ENDIF
#endif
      IF(.NOT.DMUMPS_SOLVE_IS_END_REACHED())THEN
         IF(OOC_INODE_SEQUENCE(CUR_POS_SEQUENCE,OOC_FCT_TYPE).EQ.
     $        INODE)THEN
            IF(SOLVE_STEP.EQ.0)THEN
               CUR_POS_SEQUENCE=CUR_POS_SEQUENCE+1
            ELSEIF(SOLVE_STEP.EQ.1)THEN
               CUR_POS_SEQUENCE=CUR_POS_SEQUENCE-1
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_577
      SUBROUTINE DMUMPS_591(IERR)
      USE DMUMPS_OOC_BUFFER
      IMPLICIT NONE
      INTEGER, intent(out):: IERR
      IERR=0
      IF (WITH_BUF) THEN
        CALL DMUMPS_675(IERR)
        IF(IERR.LT.0)THEN
           RETURN
        ENDIF
      END IF
      RETURN
      END SUBROUTINE DMUMPS_591
      SUBROUTINE DMUMPS_592(id,IERR)
      USE DMUMPS_OOC_BUFFER
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC), TARGET :: id
      INTEGER, intent(out) :: IERR
      INTEGER I,SOLVE_OR_FACTO
      IERR=0
      IF(WITH_BUF)THEN
         CALL DMUMPS_659()
      ENDIF
      IF(allocated(IO_REQ))THEN
         DEALLOCATE(IO_REQ)
      ENDIF
      IF(associated(KEEP_OOC))THEN
         NULLIFY(KEEP_OOC)
      ENDIF
      IF(associated(STEP_OOC))THEN
         NULLIFY(STEP_OOC)
      ENDIF
      IF(associated(PROCNODE_OOC))THEN
         NULLIFY(PROCNODE_OOC)
      ENDIF
      IF(associated(OOC_INODE_SEQUENCE))THEN
         NULLIFY(OOC_INODE_SEQUENCE)
      ENDIF
      IF(associated(TOTAL_NB_OOC_NODES))THEN
         NULLIFY(TOTAL_NB_OOC_NODES)
      ENDIF
      IF(associated(SIZE_OF_BLOCK))THEN
         NULLIFY(SIZE_OF_BLOCK)
      ENDIF
      IF(associated(OOC_VADDR))THEN
         NULLIFY(OOC_VADDR)
      ENDIF
      CALL MUMPS_OOC_END_WRITE_C(IERR)
      IF(IERR.LT.0)THEN
         IF (ICNTL1 .GT. 0)
     *   WRITE(ICNTL1,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
         GOTO 500
      ENDIF
      id%OOC_MAX_NB_NODES_FOR_ZONE=max(MAX_NB_NODES_FOR_ZONE,
     $     TMP_NB_NODES)
      IF(allocated(I_CUR_HBUF_NEXTPOS))THEN
         DO I=1,OOC_NB_FILE_TYPE
            id%OOC_TOTAL_NB_NODES(I)=I_CUR_HBUF_NEXTPOS(I)-1
         ENDDO
         DEALLOCATE(I_CUR_HBUF_NEXTPOS)
      ENDIF
      id%KEEP(203)=MAX_SIZE_FACTOR_OOC
      CALL DMUMPS_613(id,IERR)
      IF(IERR.LT.0)THEN
         GOTO 500
      ENDIF      
 500  CONTINUE
      SOLVE_OR_FACTO=0
      CALL MUMPS_CLEAN_IO_DATA_C(MYID_OOC,SOLVE_OR_FACTO,IERR)
      IF(IERR.LT.0)THEN
         IF (ICNTL1.GT.0)
     *   WRITE(ICNTL1,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
         RETURN
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_592
      SUBROUTINE DMUMPS_588(id,IERR)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      EXTERNAL MUMPS_OOC_REMOVE_FILE_C
      TYPE(DMUMPS_STRUC), TARGET :: id
      INTEGER IERR
      INTEGER I,J,I1,K
      CHARACTER*1 TMP_NAME(350)
      IERR=0
      K=1
      IF(associated(id%OOC_FILE_NAMES).AND.
     $     associated(id%OOC_FILE_NAME_LENGTH))THEN
         DO I1=1,OOC_NB_FILE_TYPE
            DO I=1,id%OOC_NB_FILES(I1)
               DO J=1,id%OOC_FILE_NAME_LENGTH(K)
                  TMP_NAME(J)=id%OOC_FILE_NAMES(K,J)
               ENDDO
               CALL MUMPS_OOC_REMOVE_FILE_C(IERR, TMP_NAME(1))
               IF(IERR.LT.0)THEN
                  IF (ICNTL1.GT.0)THEN
                     WRITE(ICNTL1,*)MYID_OOC,': ',
     *                    ERR_STR_OOC(1:DIM_ERR_STR_OOC)
                     RETURN
                  ENDIF
               ENDIF
               K=K+1
            ENDDO
         ENDDO
      ENDIF
      IF(associated(id%OOC_FILE_NAMES))THEN
         DEALLOCATE(id%OOC_FILE_NAMES)
         NULLIFY(id%OOC_FILE_NAMES)
      ENDIF
      IF(associated(id%OOC_FILE_NAME_LENGTH))THEN
         DEALLOCATE(id%OOC_FILE_NAME_LENGTH)
         NULLIFY(id%OOC_FILE_NAME_LENGTH)
      ENDIF      
      IF(associated(id%OOC_NB_FILES))THEN
         DEALLOCATE(id%OOC_NB_FILES)
         NULLIFY(id%OOC_NB_FILES)
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_588
      SUBROUTINE DMUMPS_587(id,IERR)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC), TARGET :: id
      INTEGER IERR
      IERR=0
      CALL DMUMPS_588(id,IERR)
      IF(associated(id%OOC_TOTAL_NB_NODES))THEN
         DEALLOCATE(id%OOC_TOTAL_NB_NODES)
         NULLIFY(id%OOC_TOTAL_NB_NODES)
      ENDIF
      IF(associated(id%OOC_INODE_SEQUENCE))THEN
         DEALLOCATE(id%OOC_INODE_SEQUENCE)
         NULLIFY(id%OOC_INODE_SEQUENCE)
      ENDIF
      IF(associated(id%OOC_SIZE_OF_BLOCK))THEN
         DEALLOCATE(id%OOC_SIZE_OF_BLOCK)
         NULLIFY(id%OOC_SIZE_OF_BLOCK)
      ENDIF
      IF(associated(id%OOC_VADDR))THEN
         DEALLOCATE(id%OOC_VADDR)
         NULLIFY(id%OOC_VADDR)
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_587
      SUBROUTINE DMUMPS_586(id)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      TYPE(DMUMPS_STRUC), TARGET :: id
      INTEGER TMP,TMP_SIZE,I,J
      INTEGER ZONE,INODE,SIZE,allocok,IERR,ASYNC
      IERR=0
      ICNTL1=id%ICNTL(1)
      SOLVE=.TRUE.
      N_OOC=id%N
      IF(allocated(LRLUS_SOLVE))THEN
         DEALLOCATE(LRLUS_SOLVE)
      ENDIF
      IF(allocated(LRLU_SOLVE_T))THEN
         DEALLOCATE(LRLU_SOLVE_T)
      ENDIF
      IF(allocated(LRLU_SOLVE_B))THEN
         DEALLOCATE(LRLU_SOLVE_B)
      ENDIF
      IF(allocated(POSFAC_SOLVE))THEN
         DEALLOCATE(POSFAC_SOLVE)
      ENDIF
      IF(allocated(IDEB_SOLVE_Z))THEN
         DEALLOCATE(IDEB_SOLVE_Z)
      ENDIF
      IF(allocated(PDEB_SOLVE_Z))THEN
         DEALLOCATE(PDEB_SOLVE_Z)
      ENDIF
      IF(allocated(SIZE_SOLVE_Z))THEN
         DEALLOCATE(SIZE_SOLVE_Z)
      ENDIF
      IF(allocated(CURRENT_POS_T))THEN
         DEALLOCATE(CURRENT_POS_T)
      ENDIF
      IF(allocated(CURRENT_POS_B))THEN
         DEALLOCATE(CURRENT_POS_B)
      ENDIF
      IF(allocated(POS_HOLE_T))THEN
         DEALLOCATE(POS_HOLE_T)
      ENDIF
      IF(allocated(POS_HOLE_B))THEN
         DEALLOCATE(POS_HOLE_B)
      ENDIF
      IF(allocated(OOC_STATE_NODE))THEN
         DEALLOCATE(OOC_STATE_NODE)
      ENDIF
      IF(allocated(POS_IN_MEM))THEN
         DEALLOCATE(POS_IN_MEM)
      ENDIF
      IF(allocated(INODE_TO_POS))THEN
         DEALLOCATE(INODE_TO_POS)
      ENDIF
      IF(allocated(SIZE_OF_READ))THEN
         DEALLOCATE(SIZE_OF_READ)
      ENDIF
      IF(allocated(FIRST_POS_IN_READ))THEN
         DEALLOCATE(FIRST_POS_IN_READ)
      ENDIF
      IF(allocated(READ_DEST))THEN
         DEALLOCATE(READ_DEST)
      ENDIF
      IF(allocated(READ_MNG))THEN
         DEALLOCATE(READ_MNG)
      ENDIF
      IF(allocated(REQ_TO_ZONE))THEN
         DEALLOCATE(REQ_TO_ZONE)
      ENDIF
      IF(allocated(REQ_ID))THEN
         DEALLOCATE(REQ_ID)
      ENDIF
      IF(allocated(IO_REQ))THEN
         DEALLOCATE(IO_REQ)
      ENDIF
      IF(associated(KEEP_OOC))THEN
         NULLIFY(KEEP_OOC)
      ENDIF
      IF(associated(STEP_OOC))THEN
         NULLIFY(STEP_OOC)
      ENDIF
      IF(associated(PROCNODE_OOC))THEN
         NULLIFY(PROCNODE_OOC)
      ENDIF
      IF(associated(TOTAL_NB_OOC_NODES))THEN
         NULLIFY(TOTAL_NB_OOC_NODES)
      ENDIF
      IF(associated(SIZE_OF_BLOCK))THEN
         NULLIFY(SIZE_OF_BLOCK)
      ENDIF
      IF(associated(OOC_INODE_SEQUENCE))THEN
         NULLIFY(OOC_VADDR)
      ENDIF
      DIM_ERR_STR_OOC = ERR_STR_OOC_MAX_LEN
      CALL MUMPS_LOW_LEVEL_INIT_ERR_STR(DIM_ERR_STR_OOC,ERR_STR_OOC)
      CALL DMUMPS_614(id)
      IF(id%INFO(1).LT.0)THEN
         RETURN
      ENDIF
      STEP_OOC=>id%STEP
      PROCNODE_OOC=>id%PROCNODE_STEPS
      SLAVEF_OOC=id%NSLAVES
      MYID_OOC=id%MYID
      KEEP_OOC => id%KEEP
      SIZE_OF_BLOCK=>id%OOC_SIZE_OF_BLOCK
      OOC_INODE_SEQUENCE=>id%OOC_INODE_SEQUENCE
      OOC_VADDR=>id%OOC_VADDR
      ALLOCATE(IO_REQ(id%KEEP(28)),
     $     stat=allocok)
      IF (allocok .GT. 0) THEN
         IF (ICNTL1.GT.0)
     *   WRITE(ICNTL1,*) 'PB allocation in DMUMPS_INIT_OOC'
         id%INFO(1) = -13
         id%INFO(2) = id%KEEP(28)
         RETURN
      ENDIF
      DMUMPS_ELEMENTARY_DATA_SIZE = id%KEEP(35)
      MAX_NB_NODES_FOR_ZONE=id%OOC_MAX_NB_NODES_FOR_ZONE
      TOTAL_NB_OOC_NODES=>id%OOC_TOTAL_NB_NODES
      CALL DMUMPS_711( id%KEEP(204), STRAT_IO_ASYNC,
     & WITH_BUF, LOW_LEVEL_STRAT_IO)
      IF(id%KEEP(107).GT.0)THEN
         SIZE_SOLVE_EMM=max(id%KEEP(203),
     $        int(dble(FACT_AREA_SIZE)*0.2))
         SIZE_ZONE_SOLVE=max(SIZE_SOLVE_EMM,
     $        int((dble(FACT_AREA_SIZE)-
     $        dble(SIZE_SOLVE_EMM))/dble(id%KEEP(107))))
         SIZE_ZONE_SOLVE=max(SIZE_ZONE_SOLVE,0)
         IF(SIZE_ZONE_SOLVE.EQ.SIZE_SOLVE_EMM)THEN
            SIZE_SOLVE_EMM=id%KEEP(203)
            SIZE_ZONE_SOLVE=int((dble(FACT_AREA_SIZE)-
     $           dble(SIZE_SOLVE_EMM))/dble(id%KEEP(107)))
            SIZE_ZONE_SOLVE=max(SIZE_ZONE_SOLVE,0)
         ENDIF
      ELSE
         SIZE_ZONE_SOLVE=int(dble(FACT_AREA_SIZE))
         SIZE_SOLVE_EMM=SIZE_ZONE_SOLVE
      ENDIF
      IF(SIZE_SOLVE_EMM.LT.id%KEEP(203))THEN
         IF (ICNTL1.GT.0)
     *   WRITE(ICNTL1,*)MYID_OOC,': More space needed for
     $ solution step in DMUMPS_586'
         id%INFO(1) = -11
         id%INFO(2) = id%KEEP(203)
         RETURN
      ENDIF
      TMP=MAX_NB_NODES_FOR_ZONE
      CALL MPI_ALLREDUCE(TMP,MAX_NB_NODES_FOR_ZONE,1,
     $     MPI_INTEGER,MPI_MAX,id%COMM_NODES, IERR)
      NB_Z=KEEP_OOC(107)+1
      ALLOCATE(POS_IN_MEM(MAX_NB_NODES_FOR_ZONE*NB_Z),
     $     INODE_TO_POS(KEEP_OOC(28)),
     $     stat=allocok)
      IF (allocok .GT. 0) THEN
         IF (ICNTL1.GT.0)
     *   WRITE(ICNTL1,*) 'PB allocation in DMUMPS_586'
         id%INFO(1) = -13
         id%INFO(2) = id%KEEP(28)+(MAX_NB_NODES_FOR_ZONE*NB_Z)
         RETURN
      ENDIF
      ALLOCATE(OOC_STATE_NODE(KEEP_OOC(28)),stat=allocok)
      IF (allocok .GT. 0) THEN
         IF (ICNTL1.GT.0)
     *   WRITE(ICNTL1,*) 'PB allocation in DMUMPS_586'
         id%INFO(1) = -13
         id%INFO(2) = id%KEEP(28)
         RETURN
      ENDIF
      OOC_STATE_NODE(1:KEEP_OOC(28))=0
      INODE_TO_POS=0
      POS_IN_MEM=0
      ALLOCATE(LRLUS_SOLVE(NB_Z), LRLU_SOLVE_T(NB_Z),LRLU_SOLVE_B(NB_Z),
     $     POSFAC_SOLVE(NB_Z),IDEB_SOLVE_Z(NB_Z),
     $     PDEB_SOLVE_Z(NB_Z),SIZE_SOLVE_Z(NB_Z),
     $     CURRENT_POS_T(NB_Z),CURRENT_POS_B(NB_Z),
     $     POS_HOLE_T(NB_Z),POS_HOLE_B(NB_Z),
     $     stat=allocok)
      IF (allocok .GT. 0) THEN
         IF (ICNTL1.GT.0)
     *   WRITE(ICNTL1,*) 'PB allocation in DMUMPS_INIT_OOC_SOLVE'
         id%INFO(1) = -13
         id%INFO(2) = 9*(NB_Z+1)
         RETURN
      ENDIF
      IERR=0
      CALL MUMPS_GET_MAX_NB_REQ_C(MAX_NB_REQ,IERR)
      ALLOCATE(SIZE_OF_READ(MAX_NB_REQ),FIRST_POS_IN_READ(MAX_NB_REQ),
     $     READ_DEST(MAX_NB_REQ),READ_MNG(MAX_NB_REQ),
     $     REQ_TO_ZONE(MAX_NB_REQ),REQ_ID(MAX_NB_REQ),stat=allocok)
      SIZE_OF_READ=-9999
      FIRST_POS_IN_READ=-9999
      READ_DEST=-9999
      READ_MNG=-9999
      REQ_TO_ZONE=-9999
      REQ_ID=-9999
      IF (allocok .GT. 0) THEN
         IF (ICNTL1.GT.0)
     *   WRITE(ICNTL1,*) 'PB allocation in DMUMPS_INIT_OOC_SOLVE'
         id%INFO(1) = -13
         id%INFO(2) = 6*(NB_Z+1)
         RETURN
      ENDIF
      MIN_SIZE_READ=min(max((1024*1024)/id%KEEP(35),SIZE_ZONE_SOLVE/3),
     $     SIZE_ZONE_SOLVE)
      TMP_SIZE=1
      J=1
      DO I=1,NB_Z-1
         IDEB_SOLVE_Z(I)=TMP_SIZE
         POSFAC_SOLVE(I)=TMP_SIZE
         LRLUS_SOLVE(I)=SIZE_ZONE_SOLVE
         LRLU_SOLVE_T(I)=SIZE_ZONE_SOLVE
         LRLU_SOLVE_B(I)=0
         SIZE_SOLVE_Z(I)=SIZE_ZONE_SOLVE
         CURRENT_POS_T(I)=J
         CURRENT_POS_B(I)=J
         PDEB_SOLVE_Z(I)=J
         POS_HOLE_T(I)=J
         POS_HOLE_B(I)=J
         J=J+MAX_NB_NODES_FOR_ZONE
         TMP_SIZE=TMP_SIZE+SIZE_ZONE_SOLVE
      ENDDO
      IDEB_SOLVE_Z(NB_Z)=TMP_SIZE
      PDEB_SOLVE_Z(NB_Z)=J
      POSFAC_SOLVE(NB_Z)=TMP_SIZE
      LRLUS_SOLVE(NB_Z)=SIZE_SOLVE_EMM
      LRLU_SOLVE_T(NB_Z)=SIZE_SOLVE_EMM
      LRLU_SOLVE_B(NB_Z)=0
      SIZE_SOLVE_Z(NB_Z)=SIZE_SOLVE_EMM
      CURRENT_POS_T(NB_Z)=J
      CURRENT_POS_B(NB_Z)=J
      POS_HOLE_T(NB_Z)=J
      POS_HOLE_B(NB_Z)=J
      IO_REQ=-77777
      REQ_ACT=0
      OOC_STATE_NODE(1:KEEP_OOC(28))=NOT_IN_MEM
      FIRST_ACTIVE_Z=1
      NB_ZONE_REQ=0
      SIZE_ZONE_REQ=0
      CURRENT_SOLVE_READ_ZONE=0
      NB_CALLED=0
      NB_CALL=0
      SOLVE_STEP=-9999
#if defined (NEW_PREF_SCHEME)
      MAX_PREF_SIZE=(1024*1024*2)/8
#endif
      RETURN
      END SUBROUTINE DMUMPS_586
      SUBROUTINE DMUMPS_585(A,LA,PTRFAC,NSTEPS,IERR)
      IMPLICIT NONE
      INTEGER LA,NSTEPS,IERR
      DOUBLE PRECISION A(LA)
      INTEGER PTRFAC(NSTEPS)
      INTEGER I
      IERR=0
      IF(NB_Z.GT.1)THEN
         IF(STRAT_IO_ASYNC)THEN
            DO I=1,NB_Z-1
               CALL DMUMPS_594(A,LA,PTRFAC,NSTEPS,IERR)
               IF(IERR.LT.0)THEN
                  RETURN
               ENDIF
            ENDDO
         ELSE
            CALL DMUMPS_594(A,LA,PTRFAC,NSTEPS,IERR)
            IF(IERR.LT.0)THEN
               RETURN
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_585
      SUBROUTINE DMUMPS_594(A,LA,PTRFAC,NSTEPS,IERR)
      IMPLICIT NONE
      INTEGER NSTEPS,LA,IERR
      DOUBLE PRECISION A(LA)
      INTEGER PTRFAC(NSTEPS)
      INTEGER ZONE,SIZE,DEST,FLAG
      CALL DMUMPS_603(ZONE)
      IERR=0
      CALL DMUMPS_611(ZONE,A,LA,PTRFAC,NSTEPS,IERR)
      RETURN
      END SUBROUTINE DMUMPS_594 
      SUBROUTINE DMUMPS_595(DEST,INDICE,SIZE,
     $     ZONE,PTRFAC,NSTEPS,POS_SEQ,NB_NODES,FLAG,IERR)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER ZONE,SIZE,NSTEPS,INDICE,FLAG,POS_SEQ,NB_NODES
      DOUBLE PRECISION DEST
      INTEGER PTRFAC(NSTEPS)
      INTEGER REQUEST,INODE,IERR
      INTEGER ADDR_INT1,ADDR_INT2
      INTEGER TYPE
      TYPE=OOC_SOLVE_TYPE_FCT
      IERR=0
      INODE=OOC_INODE_SEQUENCE(POS_SEQ,OOC_FCT_TYPE)
      CALL MUMPS_677(ADDR_INT1,ADDR_INT2,
     $     OOC_VADDR(STEP_OOC(INODE),OOC_FCT_TYPE))
      CALL MUMPS_LOW_LEVEL_READ_OOC_C(LOW_LEVEL_STRAT_IO,
     $     DEST,SIZE,
     $     INODE,REQUEST,TYPE,ADDR_INT1,ADDR_INT2,IERR)
      IF(IERR.LT.0)THEN
         IF (ICNTL1.GT.0)
     *   WRITE(ICNTL1,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
         RETURN
      ENDIF
      IF(STRAT_IO_ASYNC)THEN
         CALL DMUMPS_597(INODE,SIZE,INDICE,ZONE,
     $        REQUEST,POS_SEQ,NB_NODES,FLAG,PTRFAC,NSTEPS,IERR)
         IF(IERR.LT.0)THEN
            RETURN
         ENDIF
      ELSE
         CALL DMUMPS_597(INODE,SIZE,INDICE,ZONE,
     $        REQUEST,POS_SEQ,NB_NODES,FLAG,PTRFAC,NSTEPS,IERR)
         IF(IERR.LT.0)THEN
            RETURN
         ENDIF
         CALL DMUMPS_596(IO_REQ(STEP_OOC(INODE)),
     $        PTRFAC,NSTEPS)
         REQ_ACT=REQ_ACT-1
      ENDIF
      END SUBROUTINE DMUMPS_595
      SUBROUTINE DMUMPS_596(REQUEST,PTRFAC,
     $     NSTEPS)
      IMPLICIT NONE
      INTEGER NSTEPS,REQUEST
      INTEGER PTRFAC(NSTEPS)
      INTEGER ZONE
      INTEGER POS,J,POS_REQ,I,TMP_NODE,LAST,POS_IN_S,POS_IN_MANAGE
      INTEGER SIZE, STATUS
      LOGICAL DONT_USE
      EXTERNAL MUMPS_330,MUMPS_275
      INTEGER MUMPS_330,MUMPS_275
      POS_REQ=mod(REQUEST,MAX_NB_REQ)+1
      SIZE=SIZE_OF_READ(POS_REQ)
      I=FIRST_POS_IN_READ(POS_REQ)
      POS_IN_S=READ_DEST(POS_REQ)
      POS_IN_MANAGE=READ_MNG(POS_REQ)
      ZONE=REQ_TO_ZONE(POS_REQ)
      DONT_USE=.FALSE.
      J=0
      DO WHILE((J.LT.SIZE).AND.(I.LE.TOTAL_NB_OOC_NODES(OOC_FCT_TYPE)))
         TMP_NODE=OOC_INODE_SEQUENCE(I,OOC_FCT_TYPE)
         LAST=SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE)
         IF((INODE_TO_POS(STEP_OOC(TMP_NODE)).NE.0).AND.
     $        (INODE_TO_POS(STEP_OOC(TMP_NODE)).LT.
     $        -((N_OOC+1)*NB_Z)))THEN
            DONT_USE=
     $           (((MTYPE_OOC.EQ.1).AND.(KEEP_OOC(50).EQ.0).AND.
     $           (SOLVE_STEP.EQ.1).AND.
     $           ((MUMPS_330(STEP_OOC(TMP_NODE),PROCNODE_OOC,
     $           SLAVEF_OOC).EQ.2).AND.(MUMPS_275(
     $           STEP_OOC(TMP_NODE),PROCNODE_OOC,SLAVEF_OOC).NE.
     $           MYID_OOC)))
     $           .OR.
     $           ((MTYPE_OOC.NE.1).AND.(KEEP_OOC(50).EQ.0).AND.
     $           (SOLVE_STEP.EQ.0).AND.
     $           ((MUMPS_330(STEP_OOC(TMP_NODE),PROCNODE_OOC,
     $           SLAVEF_OOC).EQ.2).AND.(MUMPS_275(
     $           STEP_OOC(TMP_NODE),PROCNODE_OOC,SLAVEF_OOC).NE.
     $           MYID_OOC)))).OR.
     $           (OOC_STATE_NODE(STEP_OOC(TMP_NODE)).EQ.ALREADY_USED)
            IF(DONT_USE)THEN
               PTRFAC(STEP_OOC(TMP_NODE))=-POS_IN_S
            ELSE
               PTRFAC(STEP_OOC(TMP_NODE))=POS_IN_S
            ENDIF
            IF(abs(PTRFAC(STEP_OOC(TMP_NODE))).LT.
     $           IDEB_SOLVE_Z(ZONE))THEN
               WRITE(*,*)MYID_OOC,': Inernal error (42) in OOC ', 
     $              PTRFAC(STEP_OOC(TMP_NODE)),IDEB_SOLVE_Z(ZONE)
               CALL MUMPS_ABORT()
            ENDIF
            IF(abs(PTRFAC(STEP_OOC(TMP_NODE))).GT.
     $           (IDEB_SOLVE_Z(ZONE)+SIZE_SOLVE_Z(ZONE)-1))THEN
               WRITE(*,*)MYID_OOC,': Inernal error (43) in OOC '
               CALL MUMPS_ABORT()
            ENDIF
            IF(DONT_USE)THEN
               POS_IN_MEM(POS_IN_MANAGE)=-TMP_NODE
               INODE_TO_POS(STEP_OOC(TMP_NODE))=-POS_IN_MANAGE
               IF(OOC_STATE_NODE(STEP_OOC(TMP_NODE)).NE.
     $              ALREADY_USED)THEN
                  OOC_STATE_NODE(STEP_OOC(TMP_NODE))=USED_NOT_PERMUTED
               ENDIF
               LRLUS_SOLVE(ZONE)=LRLUS_SOLVE(ZONE)+LAST
            ELSE
               POS_IN_MEM(POS_IN_MANAGE)=TMP_NODE
               INODE_TO_POS(STEP_OOC(TMP_NODE))=POS_IN_MANAGE
               OOC_STATE_NODE(STEP_OOC(TMP_NODE))=NOT_USED
            ENDIF
            IO_REQ(STEP_OOC(TMP_NODE))=-7777
         ELSE
            POS_IN_MEM(POS_IN_MANAGE)=0
         ENDIF
         POS_IN_S=POS_IN_S+LAST
         POS_IN_MANAGE=POS_IN_MANAGE+1
         J=J+LAST
         I=I+1
      ENDDO
      SIZE_OF_READ(POS_REQ)=-9999
      FIRST_POS_IN_READ(POS_REQ)=-9999
      READ_DEST(POS_REQ)=-9999
      READ_MNG(POS_REQ)=-9999
      REQ_TO_ZONE(POS_REQ)=-9999
      REQ_ID(POS_REQ)=-9999
      RETURN
      END SUBROUTINE DMUMPS_596
      SUBROUTINE DMUMPS_597(INODE,SIZE,DEST,ZONE,
     $     REQUEST,POS_SEQ,NB_NODES,FLAG,PTRFAC,NSTEPS,IERR)
      IMPLICIT NONE
      INTEGER INODE,SIZE,ZONE,REQUEST,FLAG,POS_SEQ,NB_NODES,NSTEPS
      INTEGER PTRFAC(NSTEPS)
      INTEGER DEST,LOCAL_DEST
      INTEGER I,J,LAST,TMP_NUM_FILE,TMP_NODE,LOC_I,POS_REQ,NB
      INTEGER, intent(out) :: IERR
      IERR=0
      IF(CUR_POS_SEQUENCE.GT.TOTAL_NB_OOC_NODES(OOC_FCT_TYPE))THEN
         RETURN
      ENDIF
      NB=0
      LOCAL_DEST=DEST
      I=POS_SEQ
      POS_REQ=mod(REQUEST,MAX_NB_REQ)+1
      IF(REQ_ID(POS_REQ).NE.-9999)THEN
         CALL MUMPS_WAIT_REQUEST(REQ_ID(POS_REQ),IERR)
         IF(IERR.LT.0)THEN
            IF (ICNTL1.GT.0)
     *      WRITE(ICNTL1,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
            RETURN
         ENDIF
         CALL DMUMPS_596(REQUEST,PTRFAC,NSTEPS)
         REQ_ACT=REQ_ACT-1
      ENDIF
      SIZE_OF_READ(POS_REQ)=SIZE
      FIRST_POS_IN_READ(POS_REQ)=I
      READ_DEST(POS_REQ)=DEST
      IF(FLAG.EQ.0)THEN
         READ_MNG(POS_REQ)=CURRENT_POS_B(ZONE)-NB_NODES+1
      ELSEIF(FLAG.EQ.1)THEN
         READ_MNG(POS_REQ)=CURRENT_POS_T(ZONE)
      ENDIF
      REQ_TO_ZONE(POS_REQ)=ZONE
      REQ_ID(POS_REQ)=REQUEST
      J=0
      IF(FLAG.EQ.0)THEN
         LOC_I=CURRENT_POS_B(ZONE)-NB_NODES+1
      ENDIF
      DO WHILE((J.LT.SIZE).AND.(I.LE.TOTAL_NB_OOC_NODES(OOC_FCT_TYPE)))
         TMP_NODE=OOC_INODE_SEQUENCE(I,OOC_FCT_TYPE)
         LAST=SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE)
         IF((IO_REQ(STEP_OOC(TMP_NODE)).GE.0).OR.
     $        (INODE_TO_POS(STEP_OOC(TMP_NODE)).NE.0))THEN
            IF(FLAG.EQ.1)THEN
               POS_IN_MEM(CURRENT_POS_T(ZONE))=0
            ELSEIF(FLAG.EQ.0)THEN
               POS_IN_MEM(CURRENT_POS_B(ZONE))=0
            ENDIF
         ELSE
            IO_REQ(STEP_OOC(TMP_NODE))=REQUEST
            LRLUS_SOLVE(ZONE)=LRLUS_SOLVE(ZONE)-LAST
            IF(FLAG.EQ.1)THEN
               IF(POSFAC_SOLVE(ZONE).EQ.IDEB_SOLVE_Z(ZONE))THEN
                  POS_HOLE_B(ZONE)=-9999
                  CURRENT_POS_B(ZONE)=-9999
                  LRLU_SOLVE_B(ZONE)=0
               ENDIF
               POSFAC_SOLVE(ZONE)=POSFAC_SOLVE(ZONE)+LAST
               LRLU_SOLVE_T(ZONE)=LRLU_SOLVE_T(ZONE)-LAST
               POS_IN_MEM(CURRENT_POS_T(ZONE))=-TMP_NODE-
     $              ((N_OOC+1)*NB_Z)
               INODE_TO_POS(STEP_OOC(TMP_NODE))=-CURRENT_POS_T(ZONE)-
     $              ((N_OOC+1)*NB_Z)
               OOC_STATE_NODE(STEP_OOC(TMP_NODE))=BEING_READ
               PTRFAC(STEP_OOC(TMP_NODE))=-LOCAL_DEST
               LOCAL_DEST=LOCAL_DEST+SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),
     $              OOC_FCT_TYPE)
            ELSEIF(FLAG.EQ.0)THEN
               LRLU_SOLVE_B(ZONE)=LRLU_SOLVE_B(ZONE)-LAST
               POS_IN_MEM(LOC_I)=-TMP_NODE-((N_OOC+1)*NB_Z)
               IF(LOC_I.EQ.POS_HOLE_T(ZONE))THEN
                  IF(POS_HOLE_T(ZONE).LT.CURRENT_POS_T(ZONE))THEN
                     POS_HOLE_T(ZONE)=POS_HOLE_T(ZONE)+1
                  ENDIF
               ENDIF
               INODE_TO_POS(STEP_OOC(TMP_NODE))=-LOC_I-((N_OOC+1)*NB_Z)
               OOC_STATE_NODE(STEP_OOC(TMP_NODE))=BEING_READ
               PTRFAC(STEP_OOC(TMP_NODE))=-LOCAL_DEST
               LOCAL_DEST=LOCAL_DEST+SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),
     $              OOC_FCT_TYPE)
            ELSE
             WRITE(*,*)MYID_OOC,': Internal error (39) in OOC ',
     &                 ' Invalid Flag Value in ', 
     $                 ' DMUMPS_597',FLAG
               CALL MUMPS_ABORT()
            ENDIF
         ENDIF
         IF(POS_IN_MEM(CURRENT_POS_T(ZONE)).NE.0)THEN
            IF(POS_IN_MEM(CURRENT_POS_T(ZONE)).EQ.
     $           POS_IN_MEM(PDEB_SOLVE_Z(ZONE)))THEN
               IF(CURRENT_POS_T(ZONE).NE.PDEB_SOLVE_Z(ZONE))THEN
             WRITE(*,*)MYID_OOC,': Internal error (40) in OOC ',
     &                      CURRENT_POS_T(ZONE),
     $                 PDEB_SOLVE_Z(ZONE),
     $                 POS_IN_MEM(CURRENT_POS_T(ZONE)),
     $                 POS_IN_MEM(PDEB_SOLVE_Z(ZONE))
                  CALL MUMPS_ABORT()
               ENDIF
            ENDIF
         ENDIF
         J=J+LAST
         IF(LRLUS_SOLVE(ZONE).LT.0)THEN
             WRITE(*,*)MYID_OOC,': Internal error (41) in OOC ',
     &           ' LRLUS_SOLVE must be (1) > 0',
     $           LRLUS_SOLVE(ZONE)
            CALL MUMPS_ABORT()
         ENDIF
         I=I+1
         IF(FLAG.EQ.1)THEN
            CURRENT_POS_T(ZONE)=CURRENT_POS_T(ZONE)+1
            IF(CURRENT_POS_T(ZONE).GT.
     $           MAX_NB_NODES_FOR_ZONE+PDEB_SOLVE_Z(ZONE))THEN
               WRITE(*,*)MYID_OOC,': Internal error (1) in OOC '
               CALL MUMPS_ABORT()
            ENDIF
            POS_HOLE_T(ZONE)=CURRENT_POS_T(ZONE)
         ELSEIF(FLAG.EQ.0)THEN
            IF(POS_HOLE_B(ZONE).LT.PDEB_SOLVE_Z(ZONE))THEN
               WRITE(*,*)MYID_OOC,': Internal error (2) in OOC ', 
     $              POS_HOLE_B(ZONE),LOC_I
               CALL MUMPS_ABORT()
            ENDIF
            CURRENT_POS_B(ZONE)=CURRENT_POS_B(ZONE)-1
            POS_HOLE_B(ZONE)=CURRENT_POS_B(ZONE)
            IF(POS_HOLE_B(ZONE).LT.PDEB_SOLVE_Z(ZONE))THEN
               POS_HOLE_B(ZONE)=-9999
               LRLU_SOLVE_B(ZONE)=0
            ENDIF
         ELSE
            WRITE(*,*)MYID_OOC,': Internal error (3) in OOC ', 
     &       ' Invalid Flag Value in ', 
     $       ' DMUMPS_597',FLAG
            CALL MUMPS_ABORT()
         ENDIF
         IF(FLAG.EQ.0)THEN
            LOC_I=LOC_I+1
         ENDIF
         NB=NB+1
      ENDDO
      IF(NB.NE.NB_NODES)THEN
         WRITE(*,*)MYID_OOC,': Internal error (4) in OOC ',
     $        ' DMUMPS_597 ',NB,NB_NODES
      ENDIF
      IF(SOLVE_STEP.EQ.0)THEN
         CUR_POS_SEQUENCE=I
      ELSE
         CUR_POS_SEQUENCE=POS_SEQ-1
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_597
      SUBROUTINE DMUMPS_598(INODE,PTRFAC,NSTEPS,A,
     $     LA,FLAG,IERR)
      IMPLICIT NONE
      INTEGER LA
      INTEGER, intent(out):: IERR
      DOUBLE PRECISION A(LA)
      INTEGER INODE,NSTEPS
      INTEGER PTRFAC(NSTEPS)
      LOGICAL FLAG
      INTEGER TMP,FREE_SIZE,TMP_NODE,I,ZONE,J,FREE_HOLE_POS,
     $     FREE_HOLE_FLAG
      LOGICAL CONTIGUOUS
      INTEGER WHICH
      IERR = 0
      WHICH=-1
      IF(INODE_TO_POS(STEP_OOC(INODE)).LE.0)THEN
         WRITE(*,*)MYID_OOC,': Internal error (5) in OOC ',
     &            ' Problem in DMUMPS_598',
     &        INODE, STEP_OOC(INODE), INODE_TO_POS(STEP_OOC(INODE))
         CALL MUMPS_ABORT()
      ENDIF
      CALL DMUMPS_600(INODE,ZONE,PTRFAC,NSTEPS)
      TMP=INODE_TO_POS(STEP_OOC(INODE))
      INODE_TO_POS(STEP_OOC(INODE))=-TMP
      POS_IN_MEM(TMP)=-INODE
      PTRFAC(STEP_OOC(INODE))=-PTRFAC(STEP_OOC(INODE))
      IF(OOC_STATE_NODE(STEP_OOC(INODE)).NE.PERMUTED)THEN
         WRITE(*,*)MYID_OOC,': INTERNAL ERROR (53) in OOC',INODE,
     $        OOC_STATE_NODE(STEP_OOC(INODE))
         CALL MUMPS_ABORT()
      ENDIF
      OOC_STATE_NODE(STEP_OOC(INODE))=USED
      LRLUS_SOLVE(ZONE)=LRLUS_SOLVE(ZONE)+
     $        SIZE_OF_BLOCK(STEP_OOC(INODE),OOC_FCT_TYPE)         
      IF(LRLUS_SOLVE(ZONE).LT.0)THEN
         WRITE(*,*)MYID_OOC,': Internal error (6) in OOC ',
     &                      ': LRLUS_SOLVE must be (2) > 0'
         CALL MUMPS_ABORT()
      ENDIF
      FREE_HOLE_FLAG=0
      IF(SOLVE_STEP.EQ.0)THEN
         IF(TMP.GT.POS_HOLE_B(ZONE))THEN
            WHICH=0
         ELSEIF(TMP.LT.POS_HOLE_T(ZONE))THEN
            WHICH=1
         ENDIF
      ELSEIF(SOLVE_STEP.EQ.1)THEN
         IF(TMP.LT.POS_HOLE_T(ZONE))THEN
            WHICH=1
         ELSEIF(TMP.GT.POS_HOLE_B(ZONE))THEN
            WHICH=0
         ENDIF         
      ENDIF
#if ! defined (OLD_OOC_CODE)
      IF(WHICH.EQ.1)THEN
         J=max(PDEB_SOLVE_Z(ZONE),POS_HOLE_T(ZONE))
         J=min(J,PDEB_SOLVE_Z(ZONE)+MAX_NB_NODES_FOR_ZONE-1)
         FREE_SIZE=0
         DO I=J,TMP,-1
            IF((POS_IN_MEM(I).LT.0).AND.(POS_IN_MEM(I).GT.
     $           -(N_OOC+1)*NB_Z))THEN
               TMP_NODE=-POS_IN_MEM(I)
               FREE_SIZE=FREE_SIZE+SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),
     $              OOC_FCT_TYPE)
            ELSEIF(POS_IN_MEM(I).NE.0)THEN
               GOTO 666
            ENDIF
         ENDDO
         POS_HOLE_T(ZONE)=TMP
 666     CONTINUE
      ELSEIF(WHICH.EQ.0)THEN
         J=max(PDEB_SOLVE_Z(ZONE),POS_HOLE_B(ZONE))
         J=min(J,PDEB_SOLVE_Z(ZONE)+MAX_NB_NODES_FOR_ZONE-1)
         FREE_SIZE=0
         DO I=J,TMP
            IF((POS_IN_MEM(I).LT.0).AND.(POS_IN_MEM(I).GT.
     $           -(N_OOC+1)*NB_Z))THEN
               TMP_NODE=-POS_IN_MEM(I)
               FREE_SIZE=FREE_SIZE+SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),
     $              OOC_FCT_TYPE)
            ELSEIF(POS_IN_MEM(I).NE.0)THEN
               IF(J.EQ.PDEB_SOLVE_Z(ZONE))THEN
                  POS_HOLE_B(ZONE)=-9999
                  LRLU_SOLVE_B(ZONE)=0
                  CURRENT_POS_B(ZONE)=-9999
               ENDIF
               GOTO 777
            ENDIF
         ENDDO
         POS_HOLE_B(ZONE)=TMP
 777     CONTINUE
      ENDIF
#else
      IF((TMP.EQ.CURRENT_POS_T(ZONE)-1).OR.
     $     (TMP.EQ.POS_HOLE_T(ZONE)-1))THEN
         J=max(PDEB_SOLVE_Z(ZONE),CURRENT_POS_T(ZONE))
         FREE_SIZE=0
         DO I=CURRENT_POS_T(ZONE)-1,PDEB_SOLVE_Z(ZONE),-1
            IF((POS_IN_MEM(I).LT.0).AND.(POS_IN_MEM(I).GT.
     $           -(N_OOC+1)*NB_Z))THEN
               TMP_NODE=-POS_IN_MEM(I)
               FREE_SIZE=FREE_SIZE+SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),
     $              OOC_FCT_TYPE)
            ELSEIF(POS_IN_MEM(I).NE.0)THEN
               GOTO 666
            ENDIF
            J=J-1
         ENDDO
 666     CONTINUE
         POS_HOLE_T(ZONE)=max(J,PDEB_SOLVE_Z(ZONE))
         IF(POS_HOLE_T(ZONE).EQ.PDEB_SOLVE_Z(ZONE))THEN
            CURRENT_POS_B(ZONE)=-9999
            POS_HOLE_B(ZONE)=-9999
            LRLU_SOLVE_B(ZONE)=0
         ENDIF
      ELSEIF(TMP.GT.POS_HOLE_B(ZONE))THEN
         J=max(POS_HOLE_B(ZONE),PDEB_SOLVE_Z(ZONE))
         FREE_SIZE=0
         DO I=TMP,J,-1
            IF((POS_IN_MEM(I).LT.0).AND.(POS_IN_MEM(I).GT.
     $           -(N_OOC+1)*NB_Z))THEN
               TMP_NODE=-POS_IN_MEM(I)
               FREE_SIZE=FREE_SIZE+SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),
     $              OOC_FCT_TYPE)
            ELSEIF(POS_IN_MEM(I).NE.0)THEN
               IF(J.EQ.PDEB_SOLVE_Z(ZONE))THEN
                  POS_HOLE_B(ZONE)=-9999
                  LRLU_SOLVE_B(ZONE)=0
               ENDIF
               GOTO 777
            ENDIF
         ENDDO
         POS_HOLE_B(ZONE)=TMP
 777     CONTINUE
      ENDIF
#endif
 222  CONTINUE
      IERR=0
      IF((NB_Z.GT.1).AND.FLAG)THEN
         CALL DMUMPS_601(ZONE)
         IF((LRLUS_SOLVE(ZONE).GE.MIN_SIZE_READ).OR.
     $        (LRLUS_SOLVE(ZONE).GE.
     $        int(0.3*dble(SIZE_SOLVE_Z(ZONE)))))THEN
            CALL DMUMPS_594(A,LA,PTRFAC,NSTEPS,IERR)
            IF(IERR.LT.0)THEN
               RETURN
            ENDIF
         ELSE
            CALL DMUMPS_603(ZONE)
         ENDIF
      ENDIF     
      RETURN
      END SUBROUTINE DMUMPS_598
      FUNCTION DMUMPS_SOLVE_IS_INODE_IN_MEM(INODE,PTRFAC,NSTEPS,A,LA,
     $     IERR)
      IMPLICIT NONE
      INTEGER INODE,NSTEPS,LA
      INTEGER, INTENT(out)::IERR
      DOUBLE PRECISION A(LA)
      INTEGER PTRFAC(NSTEPS)
      INTEGER DMUMPS_SOLVE_IS_INODE_IN_MEM
      INTEGER ZONE
      IERR=0
      IF(INODE_TO_POS(STEP_OOC(INODE)).GT.0)THEN
         IF(OOC_STATE_NODE(STEP_OOC(INODE)).EQ.PERMUTED)THEN
            DMUMPS_SOLVE_IS_INODE_IN_MEM=OOC_PERMUTED
         ELSE
            DMUMPS_SOLVE_IS_INODE_IN_MEM=OOC_NOT_PERMUTED
         ENDIF
         IF(.NOT.DMUMPS_SOLVE_IS_END_REACHED())THEN
            IF(OOC_INODE_SEQUENCE(CUR_POS_SEQUENCE,OOC_FCT_TYPE)
     $           .EQ.INODE)THEN
               IF(SOLVE_STEP.EQ.0)THEN
                  CUR_POS_SEQUENCE=CUR_POS_SEQUENCE+1
               ELSEIF(SOLVE_STEP.EQ.1)THEN
                  CUR_POS_SEQUENCE=CUR_POS_SEQUENCE-1
               ENDIF
            ENDIF
         ENDIF
      ELSEIF(INODE_TO_POS(STEP_OOC(INODE)).LT.0)THEN
         IF(INODE_TO_POS(STEP_OOC(INODE)).LT.-((N_OOC+1)*NB_Z))THEN
            CALL MUMPS_WAIT_REQUEST(IO_REQ(STEP_OOC(INODE)),IERR)
            IF(IERR.LT.0)THEN
               IF (ICNTL1.GT.0)
     *         WRITE(ICNTL1,*)MYID_OOC,': Internal error (7) in OOC ',
     &                   ERR_STR_OOC(1:DIM_ERR_STR_OOC)
               RETURN
            ENDIF
            CALL DMUMPS_596(IO_REQ(STEP_OOC(INODE)),
     $           PTRFAC,NSTEPS)
            REQ_ACT=REQ_ACT-1
         ELSE
            CALL DMUMPS_599(INODE,PTRFAC,NSTEPS)
            IF(.NOT.DMUMPS_SOLVE_IS_END_REACHED())THEN
               IF(OOC_INODE_SEQUENCE(CUR_POS_SEQUENCE,OOC_FCT_TYPE).EQ.
     $              INODE)THEN
                  IF(SOLVE_STEP.EQ.0)THEN
                     CUR_POS_SEQUENCE=CUR_POS_SEQUENCE+1
                  ELSEIF(SOLVE_STEP.EQ.1)THEN
                     CUR_POS_SEQUENCE=CUR_POS_SEQUENCE-1
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF(OOC_STATE_NODE(STEP_OOC(INODE)).EQ.PERMUTED)THEN
            DMUMPS_SOLVE_IS_INODE_IN_MEM=OOC_PERMUTED
         ELSE
            DMUMPS_SOLVE_IS_INODE_IN_MEM=OOC_NOT_PERMUTED
         ENDIF
      ELSE
         DMUMPS_SOLVE_IS_INODE_IN_MEM=OOC_NOT_IN_MEM
      ENDIF
      RETURN
      END FUNCTION DMUMPS_SOLVE_IS_INODE_IN_MEM
      SUBROUTINE DMUMPS_682(INODE)
      IMPLICIT NONE
      INTEGER INODE
      IF(OOC_STATE_NODE(STEP_OOC(INODE)).NE.NOT_USED)THEN
         WRITE(*,*)MYID_OOC,': INTERNAL ERROR (51) in OOC',INODE,
     $        OOC_STATE_NODE(STEP_OOC(INODE))
         CALL MUMPS_ABORT()
      ENDIF
      OOC_STATE_NODE(STEP_OOC(INODE))=PERMUTED
      END SUBROUTINE DMUMPS_682
      SUBROUTINE DMUMPS_599(INODE,PTRFAC,NSTEPS)
      IMPLICIT NONE
      INTEGER INODE,NSTEPS
      INTEGER PTRFAC(NSTEPS)
      INTEGER ZONE
      INODE_TO_POS(STEP_OOC(INODE))=-INODE_TO_POS(STEP_OOC(INODE))
      POS_IN_MEM(INODE_TO_POS(STEP_OOC(INODE)))=
     $     -POS_IN_MEM(INODE_TO_POS(STEP_OOC(INODE)))
      PTRFAC(STEP_OOC(INODE))=-PTRFAC(STEP_OOC(INODE))
      IF(OOC_STATE_NODE(STEP_OOC(INODE)).EQ.USED_NOT_PERMUTED)THEN
         OOC_STATE_NODE(STEP_OOC(INODE))=NOT_USED
      ELSEIF(OOC_STATE_NODE(STEP_OOC(INODE)).EQ.USED)THEN
         OOC_STATE_NODE(STEP_OOC(INODE))=PERMUTED
      ELSE
         WRITE(*,*)MYID_OOC,': Internal error (52) in OOC',INODE,
     $        OOC_STATE_NODE(STEP_OOC(INODE)),
     $        INODE_TO_POS(STEP_OOC(INODE))
         CALL MUMPS_ABORT()
      ENDIF
      CALL DMUMPS_610(PTRFAC(STEP_OOC(INODE)),ZONE)
      IF(INODE_TO_POS(STEP_OOC(INODE)).LE.POS_HOLE_B(ZONE))THEN
         IF(INODE_TO_POS(STEP_OOC(INODE)).GT.
     $        PDEB_SOLVE_Z(ZONE))THEN
            POS_HOLE_B(ZONE)=
     $           INODE_TO_POS(STEP_OOC(INODE))-1
         ELSE
            CURRENT_POS_B(ZONE)=-9999
            POS_HOLE_B(ZONE)=-9999
            LRLU_SOLVE_B(ZONE)=0
         ENDIF
      ENDIF
      IF(INODE_TO_POS(STEP_OOC(INODE)).GE.POS_HOLE_T(ZONE))THEN
         IF(INODE_TO_POS(STEP_OOC(INODE)).LT.
     $        CURRENT_POS_T(ZONE)-1)THEN
            POS_HOLE_T(ZONE)=INODE_TO_POS(STEP_OOC(INODE))+1
         ELSE
            POS_HOLE_T(ZONE)=CURRENT_POS_T(ZONE)
         ENDIF
      ENDIF
      CALL DMUMPS_609(INODE,PTRFAC,NSTEPS,1)
      END SUBROUTINE DMUMPS_599
      SUBROUTINE DMUMPS_600(INODE,ZONE,PTRFAC,NSTEPS)
      IMPLICIT NONE
      INTEGER INODE,ZONE,NSTEPS
      INTEGER PTRFAC(NSTEPS)
      ZONE=1
      DO WHILE (ZONE.LE.NB_Z)
         IF(PTRFAC(STEP_OOC(INODE)).LT.IDEB_SOLVE_Z(ZONE))THEN
            ZONE=ZONE-1
            EXIT
         ENDIF
         ZONE=ZONE+1
      ENDDO
      IF(ZONE.EQ.NB_Z+1)THEN
         ZONE=ZONE-1
      ENDIF
      END SUBROUTINE DMUMPS_600
      SUBROUTINE DMUMPS_601(ZONE)
      IMPLICIT NONE
      INTEGER ZONE
      ZONE=mod((CURRENT_SOLVE_READ_ZONE+1),NB_Z-1)+1
      END SUBROUTINE DMUMPS_601
      SUBROUTINE DMUMPS_603(ZONE)
      IMPLICIT NONE
      INTEGER ZONE
      IF(NB_Z.GT.1)THEN
         CURRENT_SOLVE_READ_ZONE=mod((CURRENT_SOLVE_READ_ZONE+1),NB_Z-1)
         ZONE=CURRENT_SOLVE_READ_ZONE+1
      ELSE
         ZONE=NB_Z
      ENDIF
      END SUBROUTINE DMUMPS_603
      SUBROUTINE DMUMPS_578(INODE,PTRFAC,
     *     KEEP,KEEP8,
     $     A,IERR)
      IMPLICIT NONE
      INTEGER INODE,KEEP(500)
      INTEGER, intent(out)::IERR
      INTEGER*8 KEEP8(150)
      INTEGER PTRFAC(KEEP(28))
      DOUBLE PRECISION A(FACT_AREA_SIZE)
      INTEGER REQUESTED_SIZE
      INTEGER ZONE,IFLAG
      LOGICAL FLAG
      IERR=0
      IFLAG=0
      REQUESTED_SIZE=SIZE_OF_BLOCK(STEP_OOC(INODE),OOC_FCT_TYPE)
      ZONE=NB_Z
      IF(CURRENT_POS_T(ZONE).GT.
     $     (PDEB_SOLVE_Z(ZONE)+MAX_NB_NODES_FOR_ZONE-1))THEN
         CALL DMUMPS_608(A,FACT_AREA_SIZE,
     $        REQUESTED_SIZE,PTRFAC,KEEP(28),ZONE,IERR)
         IF(IERR.LT.0)THEN
            RETURN
         ENDIF
      ENDIF
      IF((LRLU_SOLVE_T(ZONE).GT.SIZE_OF_BLOCK(STEP_OOC(INODE),
     $     OOC_FCT_TYPE)).AND.
     $     (CURRENT_POS_T(ZONE).LE.
     $     (PDEB_SOLVE_Z(ZONE)+MAX_NB_NODES_FOR_ZONE-1)))THEN
         CALL DMUMPS_606(INODE,PTRFAC,
     *        KEEP,KEEP8,A,ZONE)
      ELSEIF(LRLU_SOLVE_B(ZONE).GT.SIZE_OF_BLOCK(STEP_OOC(INODE),
     $        OOC_FCT_TYPE).AND.
     $        (CURRENT_POS_B(ZONE).GT.0))THEN
         CALL DMUMPS_607(INODE,PTRFAC,
     *        KEEP,KEEP8,A,ZONE)
      ELSE
         IF(DMUMPS_579(INODE,ZONE))THEN
            IF(SOLVE_STEP.EQ.0)THEN
               CALL DMUMPS_604(A,FACT_AREA_SIZE,
     $              REQUESTED_SIZE,PTRFAC,
     $              KEEP(28),ZONE,IFLAG,IERR)
               IF(IERR.LT.0)THEN
                  RETURN
               ENDIF
               IF(IFLAG.EQ.1)THEN
                  CALL DMUMPS_606(INODE,PTRFAC,
     $                 KEEP,KEEP8,A,ZONE)
               ELSEIF(IFLAG.EQ.0)THEN
                  CALL DMUMPS_605(A,FACT_AREA_SIZE,
     $                 REQUESTED_SIZE,PTRFAC,
     $                 KEEP(28),ZONE,IFLAG,IERR)           
                  IF(IERR.LT.0)THEN
                     RETURN
                  ENDIF
                  IF(IFLAG.EQ.1)THEN
                     CALL DMUMPS_607(INODE,PTRFAC,
     $                    KEEP,KEEP8,A,ZONE)
                  ENDIF
               ENDIF
            ELSE
               CALL DMUMPS_605(A,FACT_AREA_SIZE,
     $              REQUESTED_SIZE,PTRFAC,
     $              KEEP(28),ZONE,IFLAG,IERR)               
               IF(IERR.LT.0)THEN
                  RETURN
               ENDIF
               IF(IFLAG.EQ.1)THEN
                  CALL DMUMPS_607(INODE,PTRFAC,
     $                 KEEP,KEEP8,A,ZONE)
               ELSEIF(IFLAG.EQ.0)THEN
                  CALL DMUMPS_604(A,FACT_AREA_SIZE,
     $                 REQUESTED_SIZE,PTRFAC,
     $                 KEEP(28),ZONE,IFLAG,IERR)
                  IF(IERR.LT.0)THEN
                     RETURN
                  ENDIF
                  IF(IFLAG.EQ.1)THEN
                     CALL DMUMPS_606(INODE,PTRFAC,
     $                    KEEP,KEEP8,A,ZONE)
                  ENDIF
               ENDIF
            ENDIF
            IF(IFLAG.EQ.0)THEN
               CALL DMUMPS_608(A,FACT_AREA_SIZE,
     $              REQUESTED_SIZE,PTRFAC,KEEP(28),ZONE,IERR)
               IF(IERR.LT.0)THEN
                  RETURN
               ENDIF
               CALL DMUMPS_606(INODE,PTRFAC,
     $              KEEP,KEEP8,A,ZONE)
            ENDIF
         ELSE
            WRITE(*,*)MYID_OOC,': Internal error (8) in OOC ',
     &                         ' Not enough space for Solve',INODE,
     $           SIZE_OF_BLOCK(STEP_OOC(INODE),OOC_FCT_TYPE),
     $           LRLUS_SOLVE(ZONE)
            CALL MUMPS_ABORT()
         ENDIF
      ENDIF
      IF(LRLUS_SOLVE(ZONE).LT.0)THEN
         WRITE(*,*)MYID_OOC,': Internal error (9) in OOC ',
     &                      ' LRLUS_SOLVE must be (3) > 0'
         CALL MUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_578
      SUBROUTINE DMUMPS_604(A,LA,REQUESTED_SIZE,PTRFAC,
     $     NSTEPS,ZONE,FLAG,IERR)
      IMPLICIT NONE
      INTEGER REQUESTED_SIZE,LA,NSTEPS,ZONE,FLAG
      INTEGER PTRFAC(NSTEPS)
      DOUBLE PRECISION A(LA)
      INTEGER I,TMP_NODE,FREE_SIZE,HOLE,FREE_HOLE_POS,FREE_HOLE_FLAG,
     $     FREE_HOLE,J
      INTEGER, intent(out)::IERR
      LOGICAL FLAG_HOLE
      IERR=0
      FLAG=0
      IF(LRLU_SOLVE_T(ZONE).EQ.SIZE_SOLVE_Z(ZONE).AND.     
     $     (.NOT.(CURRENT_POS_T(ZONE)
     $     .GT.PDEB_SOLVE_Z(ZONE)+MAX_NB_NODES_FOR_ZONE-1)))THEN
         GOTO 50
      ENDIF
      J=max(POS_HOLE_B(ZONE),PDEB_SOLVE_Z(ZONE))
      J=min(J,PDEB_SOLVE_Z(ZONE)+MAX_NB_NODES_FOR_ZONE-1)
      DO I=POS_HOLE_T(ZONE)-1,J,-1
         IF((POS_IN_MEM(I).LT.0).AND.(POS_IN_MEM(I).GT.
     $        -(N_OOC+1)*NB_Z))THEN
            TMP_NODE=-POS_IN_MEM(I)
         ELSEIF(POS_IN_MEM(I).NE.0)THEN
            EXIT
         ENDIF
      ENDDO
      POS_HOLE_T(ZONE)=I+1
      IF((POS_HOLE_T(ZONE).EQ.PDEB_SOLVE_Z(ZONE)).OR.
     $     (POS_HOLE_T(ZONE).LE.POS_HOLE_B(ZONE)).OR.
     $     (POS_HOLE_T(ZONE).EQ.POS_HOLE_B(ZONE)+1))THEN
         CURRENT_POS_B(ZONE)=-9999
         POS_HOLE_B(ZONE)=-9999
         LRLU_SOLVE_B(ZONE)=0
         POS_HOLE_T(ZONE)=PDEB_SOLVE_Z(ZONE)
      ENDIF
      FREE_HOLE=0
      FREE_SIZE=0
      FREE_HOLE_FLAG=0      
      FREE_HOLE_POS=POSFAC_SOLVE(ZONE)
      DO I=CURRENT_POS_T(ZONE)-1,POS_HOLE_T(ZONE),-1
         IF((POS_IN_MEM(I).LT.0).AND.(POS_IN_MEM(I).GT.
     $        -(N_OOC+1)*NB_Z))THEN
            TMP_NODE=-POS_IN_MEM(I)
            IF(FREE_HOLE_FLAG.EQ.1)THEN
               FREE_HOLE=FREE_HOLE_POS-
     $              (abs(PTRFAC(STEP_OOC(TMP_NODE)))+
     $              SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE))
               FREE_HOLE_FLAG=0
               FREE_SIZE=FREE_SIZE+FREE_HOLE
            ENDIF
            FREE_HOLE_POS=abs(PTRFAC(STEP_OOC(TMP_NODE)))
            PTRFAC(STEP_OOC(TMP_NODE))=-777777
            INODE_TO_POS(STEP_OOC(TMP_NODE))=0
            OOC_STATE_NODE(STEP_OOC(TMP_NODE))=ALREADY_USED
            POS_IN_MEM(I)=0
            FREE_SIZE=FREE_SIZE+SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),
     $           OOC_FCT_TYPE)
         ELSEIF(POS_IN_MEM(I).EQ.0)THEN
            FREE_HOLE_FLAG=1
         ELSEIF(POS_IN_MEM(I).NE.0)THEN
            WRITE(*,*)MYID_OOC,': Internal error (10) in OOC ',
     $            ' DMUMPS_604',
     $           CURRENT_POS_T(ZONE)-1,POS_HOLE_T(ZONE),I
            CALL MUMPS_ABORT()
         ENDIF
      ENDDO
      IF(POS_HOLE_T(ZONE).EQ.PDEB_SOLVE_Z(ZONE))THEN
         IF(FREE_HOLE_FLAG.EQ.0)THEN
            FREE_HOLE_FLAG=1
         ENDIF
      ENDIF
      IF(FREE_HOLE_FLAG.EQ.1)THEN
         IF(POS_HOLE_T(ZONE)-1.GT.PDEB_SOLVE_Z(ZONE))THEN
            I=POS_HOLE_T(ZONE)-1
            TMP_NODE=abs(POS_IN_MEM(I))
            IF(TMP_NODE.GT.(N_OOC+1)*NB_Z)THEN
               TMP_NODE=TMP_NODE-(N_OOC+1)*NB_Z
               CALL MUMPS_WAIT_REQUEST(IO_REQ(STEP_OOC(TMP_NODE)),IERR)
               IF(IERR.LT.0)THEN
                  WRITE(*,*)MYID_OOC,': Internal error (11) in OOC ',
     &                               ERR_STR_OOC(1:DIM_ERR_STR_OOC)
                  CALL MUMPS_ABORT()    
                  RETURN
               ENDIF
               REQ_ACT=REQ_ACT-1
               CALL DMUMPS_596(
     $              IO_REQ(STEP_OOC(TMP_NODE)),PTRFAC,NSTEPS)
               FREE_HOLE=FREE_HOLE_POS-
     $              (abs(PTRFAC(STEP_OOC(TMP_NODE)))+
     $              SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE))
            ELSEIF(TMP_NODE.EQ.0)THEN
               DO J=I,PDEB_SOLVE_Z(ZONE),-1
                  IF(POS_IN_MEM(J).NE.0) EXIT
               ENDDO
               IF(POS_IN_MEM(J).LT.0)THEN
                  WRITE(*,*)MYID_OOC,': Internal error (12) in OOC ',
     $                 ' DMUMPS_604'
                  CALL MUMPS_ABORT()
               ENDIF
               IF(J.GE.PDEB_SOLVE_Z(ZONE))THEN
                  TMP_NODE=POS_IN_MEM(J)
                  FREE_HOLE=FREE_HOLE_POS-
     $                 (abs(PTRFAC(STEP_OOC(TMP_NODE)))+
     $                 SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE))
               ELSE
                  FREE_HOLE=FREE_HOLE_POS-IDEB_SOLVE_Z(ZONE)
               ENDIF
            ELSEIF(TMP_NODE.LT.0)THEN
               WRITE(*,*)MYID_OOC,': Internal error (13) in OOC', 
     $           ' DMUMPS_604'
               CALL MUMPS_ABORT()
            ELSE
               FREE_HOLE=FREE_HOLE_POS-
     $              (abs(PTRFAC(STEP_OOC(TMP_NODE)))+
     $              SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE))
            ENDIF
         ELSE
            FREE_HOLE=FREE_HOLE_POS-IDEB_SOLVE_Z(ZONE)
         ENDIF
         FREE_SIZE=FREE_SIZE+FREE_HOLE
      ENDIF
      CURRENT_POS_T(ZONE)=POS_HOLE_T(ZONE)
      LRLU_SOLVE_T(ZONE)=LRLU_SOLVE_T(ZONE)+FREE_SIZE
      POSFAC_SOLVE(ZONE)=POSFAC_SOLVE(ZONE)-FREE_SIZE
 50   CONTINUE
      IF(REQUESTED_SIZE.LE.LRLU_SOLVE_T(ZONE))THEN
         FLAG=1
      ELSE
         FLAG=0
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_604
      SUBROUTINE DMUMPS_605(A,LA,REQUESTED_SIZE,
     $     PTRFAC,NSTEPS,ZONE,FLAG,IERR)
      IMPLICIT NONE
      INTEGER REQUESTED_SIZE,LA,NSTEPS,ZONE,FLAG
      INTEGER PTRFAC(NSTEPS)
      DOUBLE PRECISION A(LA)
      INTEGER I,J,TMP_NODE,FREE_SIZE,FREE_HOLE_POS,FREE_HOLE_FLAG,
     $     FREE_HOLE,TMP_LRLU
      INTEGER, intent(out) :: IERR
      IERR=0
      FLAG=0
      IF(LRLU_SOLVE_B(ZONE).EQ.SIZE_SOLVE_Z(ZONE))THEN
         GOTO 50
      ENDIF
      IF(POS_HOLE_B(ZONE).EQ.-9999)THEN
         GOTO 50
      ENDIF
      J=max(PDEB_SOLVE_Z(ZONE),POS_HOLE_T(ZONE))
      J=min(J,PDEB_SOLVE_Z(ZONE)+MAX_NB_NODES_FOR_ZONE-1)
      DO I=POS_HOLE_B(ZONE)+1,J
         IF((POS_IN_MEM(I).LT.0).AND.(POS_IN_MEM(I).GT.
     $        -(N_OOC+1)*NB_Z))THEN
            TMP_NODE=-POS_IN_MEM(I)
            FREE_SIZE=FREE_SIZE+SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),
     $           OOC_FCT_TYPE)
         ELSEIF(POS_IN_MEM(I).NE.0)THEN
            EXIT
         ENDIF
      ENDDO
      POS_HOLE_B(ZONE)=I-1
      IF((POS_HOLE_T(ZONE).EQ.PDEB_SOLVE_Z(ZONE)).OR.
     $     (POS_HOLE_T(ZONE).LE.POS_HOLE_B(ZONE)).OR.
     $     (POS_HOLE_T(ZONE).EQ.POS_HOLE_B(ZONE)+1))THEN
         CURRENT_POS_B(ZONE)=-9999
         POS_HOLE_B(ZONE)=-9999
         LRLU_SOLVE_B(ZONE)=0
         POS_HOLE_T(ZONE)=PDEB_SOLVE_Z(ZONE)
      ENDIF
      FREE_HOLE=0
      FREE_SIZE=0
      FREE_HOLE_FLAG=0
      FREE_HOLE_POS=IDEB_SOLVE_Z(ZONE)
      IF(POS_HOLE_B(ZONE).EQ.-9999)THEN
         GOTO 50
      ENDIF
      DO I=PDEB_SOLVE_Z(ZONE),POS_HOLE_B(ZONE)
         IF((POS_IN_MEM(I).LE.0).AND.(POS_IN_MEM(I).GT.
     $        -(N_OOC+1)*NB_Z))THEN
            TMP_NODE=-POS_IN_MEM(I)
            IF(TMP_NODE.NE.0)THEN
               IF(I.EQ.PDEB_SOLVE_Z(ZONE))THEN
                  IF(abs(PTRFAC(STEP_OOC(TMP_NODE))).NE.
     $                 IDEB_SOLVE_Z(ZONE))THEN
                     FREE_SIZE=FREE_SIZE+abs(PTRFAC(STEP_OOC(TMP_NODE)))
     $                    -IDEB_SOLVE_Z(ZONE)
                  ENDIF
               ENDIF
               IF(FREE_HOLE_FLAG.EQ.1)THEN
                  FREE_HOLE=abs(PTRFAC(STEP_OOC(TMP_NODE)))-
     $                 FREE_HOLE_POS
                  FREE_HOLE_FLAG=0
                  FREE_SIZE=FREE_SIZE+FREE_HOLE
               ENDIF
               FREE_HOLE_POS=abs(PTRFAC(STEP_OOC(TMP_NODE)))+
     $              SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE)
               PTRFAC(STEP_OOC(TMP_NODE))=-777777
               INODE_TO_POS(STEP_OOC(TMP_NODE))=0
               OOC_STATE_NODE(STEP_OOC(TMP_NODE))=ALREADY_USED
               FREE_SIZE=FREE_SIZE+SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),
     $              OOC_FCT_TYPE)
            ELSE
               FREE_HOLE_FLAG=1
            ENDIF
            POS_IN_MEM(I)=0
         ELSEIF(POS_IN_MEM(I).NE.0)THEN
            WRITE(*,*)MYID_OOC,': Internal error (14) in OOC ',
     $            ' DMUMPS_605',
     $           CURRENT_POS_T(ZONE)-1,POS_HOLE_B(ZONE),I,POS_IN_MEM(I)
            CALL MUMPS_ABORT()
         ENDIF
      ENDDO
      IF(FREE_HOLE_FLAG.EQ.1)THEN
         IF(POS_HOLE_B(ZONE)+1.LT.CURRENT_POS_T(ZONE)-1)THEN
            I=POS_HOLE_B(ZONE)+1
            TMP_NODE=abs(POS_IN_MEM(I))
            IF(TMP_NODE.GT.(N_OOC+1)*NB_Z)THEN
               TMP_NODE=TMP_NODE-(N_OOC+1)*NB_Z
               CALL MUMPS_WAIT_REQUEST(IO_REQ(STEP_OOC(TMP_NODE)),IERR)
               IF(IERR.LT.0)THEN
                 WRITE(*,*)MYID_OOC,': Internal error (15) in OOC ',
     &                               ERR_STR_OOC(1:DIM_ERR_STR_OOC)
                 CALL MUMPS_ABORT() 
                  RETURN
               ENDIF
               REQ_ACT=REQ_ACT-1
               CALL DMUMPS_596(
     $              IO_REQ(STEP_OOC(TMP_NODE)),PTRFAC,NSTEPS)
               FREE_HOLE=abs(PTRFAC(STEP_OOC(TMP_NODE)))-FREE_HOLE_POS
            ELSEIF(TMP_NODE.EQ.0)THEN
               DO J=I,CURRENT_POS_T(ZONE)-1
                  IF(POS_IN_MEM(J).NE.0) EXIT
               ENDDO
               IF(POS_IN_MEM(J).LT.0)THEN
                  WRITE(*,*)MYID_OOC,': Internal error (16) in OOC ',
     $                  ' DMUMPS_605'
                  CALL MUMPS_ABORT()
               ENDIF
               IF(J.LE.CURRENT_POS_T(ZONE)-1)THEN
                  TMP_NODE=POS_IN_MEM(J)
                  FREE_HOLE=abs(PTRFAC(STEP_OOC(TMP_NODE)))-
     $                 FREE_HOLE_POS
               ELSE
                  FREE_HOLE=POSFAC_SOLVE(ZONE)-FREE_HOLE_POS
               ENDIF
            ELSEIF(TMP_NODE.LT.0)THEN
               WRITE(*,*)MYID_OOC,': Internal error (17) in OOC ',
     $           ' DMUMPS_605'
               CALL MUMPS_ABORT()
            ELSE
               FREE_HOLE=abs(PTRFAC(STEP_OOC(TMP_NODE)))-
     $              FREE_HOLE_POS
            ENDIF
         ELSE
            FREE_HOLE=POSFAC_SOLVE(ZONE)-FREE_HOLE_POS
         ENDIF
         FREE_SIZE=FREE_SIZE+FREE_HOLE
      ENDIF
      LRLU_SOLVE_B(ZONE)=FREE_SIZE
      IF(POS_HOLE_B(ZONE).LT.CURRENT_POS_T(ZONE)-1)THEN
         TMP_NODE=POS_IN_MEM(POS_HOLE_B(ZONE)+1)
         IF(TMP_NODE.LT.-(N_OOC+1)*NB_Z)THEN
            TMP_NODE=abs(TMP_NODE)-(N_OOC+1)*NB_Z
            CALL MUMPS_WAIT_REQUEST(IO_REQ(STEP_OOC(TMP_NODE)),IERR)
            IF(IERR.LT.0)THEN
               WRITE(*,*)MYID_OOC,': Internal error (18) in OOC ',
     &                            ERR_STR_OOC(1:DIM_ERR_STR_OOC)
               CALL MUMPS_ABORT() 
               RETURN
            ENDIF
            REQ_ACT=REQ_ACT-1
            CALL DMUMPS_596(
     $           IO_REQ(STEP_OOC(TMP_NODE)),PTRFAC,NSTEPS)
         ENDIF
         LRLU_SOLVE_B(ZONE)=LRLU_SOLVE_B(ZONE)+
     $        (abs(PTRFAC(STEP_OOC(abs(TMP_NODE))))-IDEB_SOLVE_Z(ZONE)-
     $        LRLU_SOLVE_B(ZONE))
      ENDIF
      CURRENT_POS_B(ZONE)=POS_HOLE_B(ZONE)
 50   CONTINUE
      IF((POS_HOLE_B(ZONE).EQ.-9999).AND.(LRLU_SOLVE_B(ZONE).NE.0))THEN
         WRITE(*,*)MYID_OOC,': Internal error (19) in OOC ',
     &             'DMUMPS_605'
         CALL MUMPS_ABORT()
      ENDIF
      IF((REQUESTED_SIZE.LE.LRLU_SOLVE_B(ZONE)).AND.
     $     (POS_HOLE_B(ZONE).NE.-9999))THEN
         FLAG=1
      ELSE
         FLAG=0
      ENDIF
      END SUBROUTINE DMUMPS_605
      SUBROUTINE DMUMPS_606(INODE,PTRFAC,
     *           KEEP,KEEP8, A,ZONE)
      IMPLICIT NONE
      INTEGER INODE,KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER PTRFAC(KEEP(28))
      DOUBLE PRECISION A(FACT_AREA_SIZE)
      INTEGER REQUESTED_SIZE
      INTEGER ZONE,IFLAG
      LOGICAL FLAG
      LRLU_SOLVE_T(ZONE)=LRLU_SOLVE_T(ZONE)-
     $     SIZE_OF_BLOCK(STEP_OOC(INODE),OOC_FCT_TYPE)
      LRLUS_SOLVE(ZONE)=LRLUS_SOLVE(ZONE)-
     $     SIZE_OF_BLOCK(STEP_OOC(INODE),OOC_FCT_TYPE)
      PTRFAC(STEP_OOC(INODE))=POSFAC_SOLVE(ZONE)
      OOC_STATE_NODE(STEP_OOC(INODE))=NOT_USED
      IF(POSFAC_SOLVE(ZONE).EQ.IDEB_SOLVE_Z(ZONE))THEN
         POS_HOLE_B(ZONE)=-9999
         CURRENT_POS_B(ZONE)=-9999
         LRLU_SOLVE_B(ZONE)=0
      ENDIF
      IF(PTRFAC(STEP_OOC(INODE)).LT.IDEB_SOLVE_Z(ZONE))THEN
         WRITE(*,*)MYID_OOC,': Internal error (20) in OOC ',
     &                      ' Problem avec debut (2)',INODE,
     $              PTRFAC(STEP_OOC(INODE)),IDEB_SOLVE_Z(ZONE),ZONE
         CALL MUMPS_ABORT()
      ENDIF
      INODE_TO_POS(STEP_OOC(INODE))=CURRENT_POS_T(ZONE)
      POS_IN_MEM(CURRENT_POS_T(ZONE))=INODE
      IF(CURRENT_POS_T(ZONE).GT.(PDEB_SOLVE_Z(ZONE)+
     $     MAX_NB_NODES_FOR_ZONE-1))THEN
         WRITE(*,*)MYID_OOC,': Internal error (21) in OOC ',
     &                      ' Problem with CURRENT_POS_T',
     $        CURRENT_POS_T(ZONE),ZONE
         CALL MUMPS_ABORT()
      ENDIF
      CURRENT_POS_T(ZONE)=CURRENT_POS_T(ZONE)+1
      POS_HOLE_T(ZONE)=CURRENT_POS_T(ZONE)
      POS_HOLE_T(ZONE)=CURRENT_POS_T(ZONE)
      POSFAC_SOLVE(ZONE)=POSFAC_SOLVE(ZONE)+
     $     SIZE_OF_BLOCK(STEP_OOC(INODE),OOC_FCT_TYPE)
      END SUBROUTINE DMUMPS_606
      SUBROUTINE DMUMPS_607(INODE,PTRFAC,
     *     KEEP,KEEP8,
     $     A,ZONE)
      IMPLICIT NONE
      INTEGER INODE,KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER PTRFAC(KEEP(28))
      DOUBLE PRECISION A(FACT_AREA_SIZE)
      INTEGER REQUESTED_SIZE
      INTEGER ZONE,IFLAG,TMP_NODE
      LOGICAL FLAG
      IF(POS_HOLE_B(ZONE).EQ.-9999)THEN
         WRITE(*,*)MYID_OOC,': Internal error (22) in OOC ',
     $        ' DMUMPS_607'
         CALL MUMPS_ABORT()
      ENDIF
      LRLUS_SOLVE(ZONE)=LRLUS_SOLVE(ZONE)-
     $     SIZE_OF_BLOCK(STEP_OOC(INODE),OOC_FCT_TYPE)
      LRLU_SOLVE_B(ZONE)=LRLU_SOLVE_B(ZONE)-
     $        SIZE_OF_BLOCK(STEP_OOC(INODE),OOC_FCT_TYPE)
      PTRFAC(STEP_OOC(INODE))=IDEB_SOLVE_Z(ZONE)+
     $     LRLU_SOLVE_B(ZONE)
      OOC_STATE_NODE(STEP_OOC(INODE))=NOT_USED
      IF(PTRFAC(STEP_OOC(INODE)).LT.IDEB_SOLVE_Z(ZONE))THEN
         WRITE(*,*)MYID_OOC,': Internal error (23) in OOC ',
     $              PTRFAC(STEP_OOC(INODE)),IDEB_SOLVE_Z(ZONE)
         CALL MUMPS_ABORT()
      ENDIF
      INODE_TO_POS(STEP_OOC(INODE))=CURRENT_POS_B(ZONE)
      IF(CURRENT_POS_B(ZONE).EQ.0)THEN
         WRITE(*,*)MYID_OOC,': Internal error (23b) in OOC '
         CALL MUMPS_ABORT()
      ENDIF
      POS_IN_MEM(CURRENT_POS_B(ZONE))=INODE
      CURRENT_POS_B(ZONE)=CURRENT_POS_B(ZONE)-1
      POS_HOLE_B(ZONE)=CURRENT_POS_B(ZONE)
      END SUBROUTINE DMUMPS_607
      SUBROUTINE DMUMPS_608(A,LA,REQUESTED_SIZE,PTRFAC,
     $     NSTEPS,ZONE,IERR)
      IMPLICIT NONE
      INTEGER REQUESTED_SIZE,LA,NSTEPS,ZONE
      INTEGER, intent(out) :: IERR
      INTEGER PTRFAC(NSTEPS)
      DOUBLE PRECISION A(LA)
      INTEGER J,I,TMP_NODE,FIRST_FREE,SIZE_HOLE,K,
     $     APOS_FIRST_FREE,NB_FREE,IPOS_FIRST_FREE
      INTEGER AREA_POINTER,FREE_HOLE_POS,FREE_HOLE_FLAG,FREE_HOLE,TMPP
      IERR=0
      IF(LRLU_SOLVE_T(ZONE).EQ.SIZE_SOLVE_Z(ZONE))THEN
         RETURN
      ENDIF
      AREA_POINTER=IDEB_SOLVE_Z(ZONE)
      SIZE_HOLE=0
      DO I=PDEB_SOLVE_Z(ZONE),CURRENT_POS_T(ZONE)-1
         IF((POS_IN_MEM(I).LE.0).AND.
     $        (POS_IN_MEM(I).GT.-((N_OOC+1)*NB_Z))) GOTO 666
         TMP_NODE=abs(POS_IN_MEM(I))
         IF(TMP_NODE.GT.((N_OOC+1)*NB_Z))THEN
            TMP_NODE=TMP_NODE-((N_OOC+1)*NB_Z)
         ENDIF
         AREA_POINTER=AREA_POINTER+
     $        abs(SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE))
      ENDDO
 666  CONTINUE
      IF((I.EQ.CURRENT_POS_T(ZONE)-1).AND.
     $     (PDEB_SOLVE_Z(ZONE).NE.CURRENT_POS_T(ZONE)-1))THEN
         IF((POS_IN_MEM(I).GT.0).OR.
     $        (POS_IN_MEM(I).LT.-((N_OOC+1)*NB_Z)))THEN
            WRITE(*,*)MYID_OOC,': Internal error (25) in OOC ',
     &                      ': There are no free blocks ',
     $         'in DMUMPS_608',PDEB_SOLVE_Z(ZONE),
     $           CURRENT_POS_T(ZONE)
            CALL MUMPS_ABORT()
         ENDIF
      ENDIF
      IF(POS_IN_MEM(I).EQ.0)THEN
         APOS_FIRST_FREE=AREA_POINTER
         FREE_HOLE_POS=AREA_POINTER
      ELSE
         TMP_NODE=abs(POS_IN_MEM(I))
         APOS_FIRST_FREE=abs(PTRFAC(STEP_OOC(TMP_NODE)))
      ENDIF
      IF(POS_IN_MEM(PDEB_SOLVE_Z(ZONE)).NE.0)THEN
         IF(POS_IN_MEM(PDEB_SOLVE_Z(ZONE)).LT.-((N_OOC+1)*NB_Z))THEN
            TMP_NODE=abs(POS_IN_MEM(PDEB_SOLVE_Z(ZONE)))-
     $           ((N_OOC+1)*NB_Z)
            CALL MUMPS_WAIT_REQUEST(IO_REQ(STEP_OOC(TMP_NODE)),IERR)
            IF(IERR.LT.0)THEN
               RETURN
            ENDIF
            REQ_ACT=REQ_ACT-1
            CALL DMUMPS_596(
     $           IO_REQ(STEP_OOC(TMP_NODE)),PTRFAC,NSTEPS)
         ELSE
            TMP_NODE=abs(POS_IN_MEM(PDEB_SOLVE_Z(ZONE)))            
         ENDIF
         IF(abs(PTRFAC(STEP_OOC(TMP_NODE))).NE.IDEB_SOLVE_Z(ZONE))THEN
            IF((POS_IN_MEM(I).NE.0).OR.(I.EQ.CURRENT_POS_T(ZONE)))THEN
               SIZE_HOLE=abs(PTRFAC(STEP_OOC(TMP_NODE)))-
     $              IDEB_SOLVE_Z(ZONE)
            ENDIF
            APOS_FIRST_FREE=IDEB_SOLVE_Z(ZONE)
            IF(POS_IN_MEM(PDEB_SOLVE_Z(ZONE)).GT.0)THEN
               DO J=PDEB_SOLVE_Z(ZONE),I-1
                  TMP_NODE=POS_IN_MEM(J)
                  IF(TMP_NODE.LE.0)THEN
                     IF(TMP_NODE.LT.-((N_OOC+1)*NB_Z))THEN
                        TMP_NODE=abs(POS_IN_MEM(J))-((N_OOC+1)*NB_Z)
                        CALL MUMPS_WAIT_REQUEST(
     $                       IO_REQ(STEP_OOC(TMP_NODE)),IERR)
                        IF(IERR.LT.0)THEN
                           RETURN
                        ENDIF
                        REQ_ACT=REQ_ACT-1
                        CALL DMUMPS_596(
     $                       IO_REQ(STEP_OOC(TMP_NODE)),PTRFAC,NSTEPS)
                        TMP_NODE=POS_IN_MEM(J)
                     ELSE
                    WRITE(*,*)MYID_OOC,': Internal error (26) in OOC ',
     $                      ' DMUMPS_608',TMP_NODE,
     $                       J,I-1,(N_OOC+1)*NB_Z
                        CALL MUMPS_ABORT()
                     ENDIF
                  ENDIF
                  DO K=1,SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE)
                     A(APOS_FIRST_FREE+K-1)=
     $                    A(PTRFAC(STEP_OOC(TMP_NODE))+K-1)
                  ENDDO
                  PTRFAC(STEP_OOC(TMP_NODE))=APOS_FIRST_FREE
                  APOS_FIRST_FREE=APOS_FIRST_FREE+
     $                 SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE)
               ENDDO
            ENDIF
         ENDIF
      ENDIF
      NB_FREE=0
      FREE_HOLE=0
      FREE_HOLE_FLAG=0
      DO J=I,CURRENT_POS_T(ZONE)-1
         TMP_NODE=abs(POS_IN_MEM(J))
         IF(POS_IN_MEM(J).LT.-((N_OOC+1)*NB_Z))THEN               
            TMP_NODE=TMP_NODE-((N_OOC+1)*NB_Z)
            CALL MUMPS_WAIT_REQUEST(IO_REQ(STEP_OOC(TMP_NODE)),IERR)
            IF(IERR.LT.0)THEN
               RETURN
            ENDIF
            REQ_ACT=REQ_ACT-1
            CALL DMUMPS_596(
     $           IO_REQ(STEP_OOC(TMP_NODE)),PTRFAC,NSTEPS)
            TMP_NODE=abs(POS_IN_MEM(J))
         ENDIF
         IF(POS_IN_MEM(J).GT.0)THEN
            DO K=1,SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE)
               A(APOS_FIRST_FREE+K-1)=A(PTRFAC(STEP_OOC(TMP_NODE))+K-1)
            ENDDO
            IF(FREE_HOLE_FLAG.EQ.1)THEN
               FREE_HOLE=abs(PTRFAC(STEP_OOC(TMP_NODE)))-
     $              FREE_HOLE_POS
               FREE_HOLE_FLAG=0
               SIZE_HOLE=SIZE_HOLE+FREE_HOLE               
            ENDIF
            FREE_HOLE_POS=abs(PTRFAC(STEP_OOC(TMP_NODE)))+
     $           SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE)
            PTRFAC(STEP_OOC(TMP_NODE))=APOS_FIRST_FREE
            APOS_FIRST_FREE=APOS_FIRST_FREE+
     $           SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE)
         ELSEIF(POS_IN_MEM(J).EQ.0)THEN
            FREE_HOLE_FLAG=1
            NB_FREE=NB_FREE+1
         ELSE                   
            NB_FREE=NB_FREE+1
            IF(FREE_HOLE_FLAG.EQ.1)THEN
               FREE_HOLE=abs(PTRFAC(STEP_OOC(TMP_NODE)))-
     $              FREE_HOLE_POS
               FREE_HOLE_FLAG=0
               SIZE_HOLE=SIZE_HOLE+FREE_HOLE
            ENDIF
            FREE_HOLE_POS=abs(PTRFAC(STEP_OOC(TMP_NODE)))+
     $           SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE)
            SIZE_HOLE=SIZE_HOLE+SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),
     $           OOC_FCT_TYPE)
            PTRFAC(STEP_OOC(abs(POS_IN_MEM(J))))=-77777
         ENDIF
      ENDDO
      IF(FREE_HOLE_FLAG.EQ.1)THEN
         FREE_HOLE=POSFAC_SOLVE(ZONE)-FREE_HOLE_POS
         FREE_HOLE_FLAG=0
         SIZE_HOLE=SIZE_HOLE+FREE_HOLE
      ENDIF
      IPOS_FIRST_FREE=I
      DO J=I,CURRENT_POS_T(ZONE)-1
         IF(POS_IN_MEM(J).LT.0)THEN
            TMP_NODE=abs(POS_IN_MEM(J))
            INODE_TO_POS(STEP_OOC(TMP_NODE))=0
            POS_IN_MEM(J)=0
            OOC_STATE_NODE(STEP_OOC(TMP_NODE))=ALREADY_USED
          ELSEIF(POS_IN_MEM(J).GT.0)THEN
             TMP_NODE=abs(POS_IN_MEM(J))
             POS_IN_MEM(IPOS_FIRST_FREE)=POS_IN_MEM(J)
             INODE_TO_POS(STEP_OOC(TMP_NODE))=IPOS_FIRST_FREE
             IPOS_FIRST_FREE=IPOS_FIRST_FREE+1
         ENDIF
      ENDDO
 777  CONTINUE
      LRLU_SOLVE_T(ZONE)=LRLU_SOLVE_T(ZONE)+SIZE_HOLE
      POSFAC_SOLVE(ZONE)=POSFAC_SOLVE(ZONE)-SIZE_HOLE
      CURRENT_POS_T(ZONE)=CURRENT_POS_T(ZONE)-NB_FREE
      POS_HOLE_T(ZONE)=CURRENT_POS_T(ZONE)
      LRLU_SOLVE_B(ZONE)=0
      POS_HOLE_B(ZONE)=-9999
      CURRENT_POS_B(ZONE)=-9999
      LRLU_SOLVE_B(ZONE)=0
      IF(LRLU_SOLVE_T(ZONE).NE.LRLUS_SOLVE(ZONE))THEN
         WRITE(*,*)MYID_OOC,': Internal error (27) in OOC ',
     &                 LRLU_SOLVE_T(ZONE),
     $                 LRLUS_SOLVE(ZONE)
         CALL MUMPS_ABORT()
      ENDIF
      LRLU_SOLVE_T(ZONE)=LRLUS_SOLVE(ZONE)
      IF(LRLUS_SOLVE(ZONE).LT.0)THEN
         WRITE(*,*)MYID_OOC,': Internal error (28) in OOC ',
     &                      ' LRLUS_SOLVE must be (4) > 0'
         CALL MUMPS_ABORT()
      ENDIF
      IF(POSFAC_SOLVE(ZONE).LT.IDEB_SOLVE_Z(ZONE))THEN
         WRITE(*,*)MYID_OOC,': Internal error (29) in OOC ',
     $        POSFAC_SOLVE(ZONE),IDEB_SOLVE_Z(ZONE)
         CALL MUMPS_ABORT()
      ENDIF
      IF(POSFAC_SOLVE(ZONE).NE.(IDEB_SOLVE_Z(ZONE)+SIZE_SOLVE_Z(ZONE)-
     $     LRLUS_SOLVE(ZONE)))THEN
         WRITE(*,*)MYID_OOC,': Internal error (30) in OOC ',
     &                      ' Problem avec debut POSFAC_SOLVE',
     $        POSFAC_SOLVE(ZONE),(SIZE_SOLVE_Z(ZONE)-
     $     LRLUS_SOLVE(ZONE))+IDEB_SOLVE_Z(ZONE),LRLUS_SOLVE(ZONE)
         CALL MUMPS_ABORT()
      ENDIF      
      IF(POSFAC_SOLVE(ZONE).GT.
     $     (IDEB_SOLVE_Z(ZONE)+SIZE_SOLVE_Z(ZONE)-1))THEN
         WRITE(*,*)MYID_OOC,': Internal error (31) in OOC ',
     $        POSFAC_SOLVE(ZONE),IDEB_SOLVE_Z(ZONE)+
     $        SIZE_SOLVE_Z(ZONE)-1
         CALL MUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_608
      SUBROUTINE DMUMPS_581()
      IMPLICIT NONE
      IS_INIT_CALLED=.FALSE.
      END SUBROUTINE DMUMPS_581
      SUBROUTINE DMUMPS_609(INODE,PTRFAC,NSTEPS,FLAG)
      IMPLICIT NONE
      INTEGER INODE,NSTEPS,FLAG
      INTEGER PTRFAC(NSTEPS)
      INTEGER ZONE
      IF((FLAG.LT.0).OR.(FLAG.GT.1))THEN
         WRITE(*,*)MYID_OOC,': Internal error (32) in OOC ',
     $        ' DMUMPS_609'
         CALL MUMPS_ABORT()
      ENDIF
      CALL DMUMPS_610(PTRFAC(STEP_OOC(INODE)),ZONE)
      IF(LRLUS_SOLVE(ZONE).LT.0)THEN
         WRITE(*,*)MYID_OOC,': Internal error (33) in OOC ',
     &                      ' LRLUS_SOLVE must be (5) ++ > 0'
         CALL MUMPS_ABORT()
      ENDIF      
      IF(FLAG.EQ.0)THEN
         LRLUS_SOLVE(ZONE)=LRLUS_SOLVE(ZONE)+
     $        SIZE_OF_BLOCK(STEP_OOC(INODE),OOC_FCT_TYPE)
      ELSE
         LRLUS_SOLVE(ZONE)=LRLUS_SOLVE(ZONE)-
     $        SIZE_OF_BLOCK(STEP_OOC(INODE),OOC_FCT_TYPE)         
      ENDIF
      IF(LRLUS_SOLVE(ZONE).LT.0)THEN
         WRITE(*,*)MYID_OOC,': Internal error (34) in OOC ',
     &                      ' LRLUS_SOLVE must be (5) > 0'
         CALL MUMPS_ABORT()
      ENDIF
      END SUBROUTINE DMUMPS_609
      SUBROUTINE DMUMPS_610(ADDR,ZONE)
      IMPLICIT NONE
      INTEGER ADDR,ZONE
      INTEGER I
      I=1
      DO WHILE (I.LE.NB_Z)
         IF(ADDR.LT.IDEB_SOLVE_Z(I))THEN
            EXIT
         ENDIF
         I=I+1
      ENDDO
      ZONE=I-1
      END SUBROUTINE DMUMPS_610
      FUNCTION DMUMPS_SOLVE_IS_END_REACHED()
      IMPLICIT NONE
      LOGICAL DMUMPS_SOLVE_IS_END_REACHED
      DMUMPS_SOLVE_IS_END_REACHED=.FALSE.
      IF(SOLVE_STEP.EQ.0)THEN
         IF(CUR_POS_SEQUENCE.GT.TOTAL_NB_OOC_NODES(OOC_FCT_TYPE))THEN
            DMUMPS_SOLVE_IS_END_REACHED=.TRUE.
         ENDIF
      ELSEIF(SOLVE_STEP.EQ.1)THEN
         IF(CUR_POS_SEQUENCE.LT.1)THEN
            DMUMPS_SOLVE_IS_END_REACHED=.TRUE.
         ENDIF
      ENDIF
      RETURN
      END FUNCTION DMUMPS_SOLVE_IS_END_REACHED
      SUBROUTINE DMUMPS_611(ZONE,A,LA,PTRFAC,NSTEPS,IERR)
      IMPLICIT NONE
      INTEGER NSTEPS,LA,ZONE
      INTEGER, intent(out) :: IERR
      DOUBLE PRECISION A(LA)
      INTEGER PTRFAC(NSTEPS)
      INTEGER FLAG,TMP_FLAG,SIZE,DEST,POS_SEQ,NEEDED_SIZE,TMP_NODE,
     $     NB_NODES
      EXTERNAL MUMPS_330,MUMPS_275
      INTEGER MUMPS_330,MUMPS_275
      IERR=0
      TMP_FLAG=0
      FLAG=0
      IF(DMUMPS_SOLVE_IS_END_REACHED())THEN
         RETURN
      ENDIF      
      IF(SOLVE_STEP.EQ.0)THEN
         IF(CUR_POS_SEQUENCE.LE.TOTAL_NB_OOC_NODES(OOC_FCT_TYPE))THEN
            TMP_NODE=OOC_INODE_SEQUENCE(CUR_POS_SEQUENCE,
     $           OOC_FCT_TYPE)
            DO WHILE(SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE).GT.
     $           SIZE_SOLVE_Z(ZONE))
               CUR_POS_SEQUENCE=CUR_POS_SEQUENCE+1
               IF(DMUMPS_SOLVE_IS_END_REACHED())THEN
                  RETURN
               ENDIF            
               TMP_NODE=OOC_INODE_SEQUENCE(CUR_POS_SEQUENCE,
     $              OOC_FCT_TYPE)
            ENDDO      
            NEEDED_SIZE=max(MIN_SIZE_READ,
     $           SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE))
         ELSE
            NEEDED_SIZE=MIN_SIZE_READ
         ENDIF
      ELSEIF(SOLVE_STEP.EQ.1)THEN
         IF(CUR_POS_SEQUENCE.GE.1)THEN
            TMP_NODE=OOC_INODE_SEQUENCE(CUR_POS_SEQUENCE,
     $           OOC_FCT_TYPE)
            DO WHILE(SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE).GT.
     $           SIZE_SOLVE_Z(ZONE))
               CUR_POS_SEQUENCE=CUR_POS_SEQUENCE-1
               IF(DMUMPS_SOLVE_IS_END_REACHED())THEN
                  RETURN
               ENDIF              
               TMP_NODE=OOC_INODE_SEQUENCE(CUR_POS_SEQUENCE,
     $              OOC_FCT_TYPE)
            ENDDO  
            NEEDED_SIZE=max(MIN_SIZE_READ,
     $           SIZE_OF_BLOCK(STEP_OOC(TMP_NODE),OOC_FCT_TYPE))
         ELSE
            NEEDED_SIZE=MIN_SIZE_READ
         ENDIF
      ENDIF
      IF(LRLUS_SOLVE(ZONE).LT.NEEDED_SIZE)THEN
         RETURN
      ELSEIF((LRLU_SOLVE_T(ZONE).LT.NEEDED_SIZE).AND.
     $        (LRLU_SOLVE_B(ZONE).LT.NEEDED_SIZE).AND.
     $        (LRLUS_SOLVE(ZONE).LT.0.3*SIZE_SOLVE_Z(ZONE)))THEN
         RETURN
      ENDIF
      IF((LRLU_SOLVE_T(ZONE).GT.NEEDED_SIZE).AND.(SOLVE_STEP.EQ.0).AND.
     $     ((CURRENT_POS_T(ZONE)-PDEB_SOLVE_Z(ZONE)+1).LT.
     $     MAX_NB_NODES_FOR_ZONE))THEN
         FLAG=1
      ELSE
         IF(SOLVE_STEP.EQ.0)THEN
            CALL DMUMPS_604(A,FACT_AREA_SIZE,
     $           NEEDED_SIZE,PTRFAC,NSTEPS,ZONE,TMP_FLAG,IERR)
            IF(IERR.LT.0)THEN
               RETURN
            ENDIF
            FLAG=1
            IF(TMP_FLAG.EQ.0)THEN
               CALL DMUMPS_605(A,FACT_AREA_SIZE,
     $              NEEDED_SIZE,PTRFAC,NSTEPS,ZONE,TMP_FLAG,IERR)
               IF(IERR.LT.0)THEN
                  RETURN
               ENDIF
               FLAG=0
            ENDIF
         ELSE
            CALL DMUMPS_605(A,FACT_AREA_SIZE,
     $           NEEDED_SIZE,PTRFAC,NSTEPS,ZONE,TMP_FLAG,IERR)
            IF(IERR.LT.0)THEN
               RETURN
            ENDIF
            FLAG=0
            IF(TMP_FLAG.EQ.0)THEN
               CALL DMUMPS_604(A,FACT_AREA_SIZE,
     $              NEEDED_SIZE,PTRFAC,NSTEPS,ZONE,TMP_FLAG,IERR)
               IF(IERR.LT.0)THEN
                  RETURN
               ENDIF
               FLAG=1
            ENDIF
         ENDIF
         IF(TMP_FLAG.EQ.0)THEN
            CALL DMUMPS_608(A,FACT_AREA_SIZE,
     $           NEEDED_SIZE,PTRFAC,NSTEPS,ZONE,IERR)
            IF(IERR.LT.0)THEN
               RETURN
            ENDIF
            FLAG=1
         ENDIF
      ENDIF
#if defined (NEW_PREF_SCHEME)
      DO WHILE (REQ_ACT.LE.MAX_NB_REQ)
#endif         
      CALL DMUMPS_602(ZONE,SIZE,DEST,POS_SEQ,
     $     NB_NODES,FLAG,PTRFAC,NSTEPS)
      IF(SIZE.EQ.0)THEN
         RETURN
      ENDIF
      NB_ZONE_REQ=NB_ZONE_REQ+1
      SIZE_ZONE_REQ=SIZE_ZONE_REQ+SIZE
      REQ_ACT=REQ_ACT+1
      CALL DMUMPS_595(A(DEST),DEST,SIZE,ZONE,PTRFAC,NSTEPS,
     $     POS_SEQ,NB_NODES,FLAG,IERR)
      IF(IERR.LT.0)THEN
         RETURN
      ENDIF
#if defined (NEW_PREF_SCHEME)
      ENDDO
#endif
      END SUBROUTINE DMUMPS_611
      SUBROUTINE DMUMPS_602(ZONE,SIZE,DEST,POS_SEQ,
     $     NB_NODES,FLAG,PTRFAC,NSTEPS)
      IMPLICIT NONE
      INTEGER ZONE,SIZE,DEST,FLAG,POS_SEQ,NSTEPS
      INTEGER PTRFAC(NSTEPS)
      INTEGER I,J,LAST,TMP_NUM_FILE,START_NODE,MAX_SIZE,K,TMP,MAX_NB,
     $     NB_NODES
      LOGICAL ALREADY
      IF(DMUMPS_SOLVE_IS_END_REACHED())THEN
         SIZE=0
         RETURN
      ENDIF
      IF(FLAG.EQ.0)THEN
         MAX_SIZE=LRLU_SOLVE_B(ZONE)
         MAX_NB=max(0,CURRENT_POS_B(ZONE)-PDEB_SOLVE_Z(ZONE)+1)
      ELSEIF(FLAG.EQ.1)THEN
         MAX_SIZE=LRLU_SOLVE_T(ZONE)
         MAX_NB=MAX_NB_NODES_FOR_ZONE
      ELSE
         WRITE(*,*)MYID_OOC,': Internal error (35) in OOC ',
     &                      ' Unknown Flag value in ',
     $         ' DMUMPS_602',FLAG
         CALL MUMPS_ABORT()
      ENDIF
      I=CUR_POS_SEQUENCE
      START_NODE=OOC_INODE_SEQUENCE(I,OOC_FCT_TYPE)
      ALREADY=.FALSE.
      NB_NODES=0
#if defined (NEW_PREF_SCHEME)
      IF(MAX_SIZE.GT.MAX_PREF_SIZE)THEN
         MAX_SIZE=MIN(MAX(MAX_PREF_SIZE,
     $        SIZE_OF_BLOCK(STEP_OOC(START_NODE),OOC_FCT_TYPE)),
     $        MAX_SIZE)
      ENDIF
#endif
      IF(ZONE.EQ.NB_Z)THEN
         SIZE=SIZE_OF_BLOCK(STEP_OOC(START_NODE),OOC_FCT_TYPE)
      ELSE
         J=0
         IF(FLAG.EQ.0)THEN
            K=0
         ELSEIF(FLAG.EQ.1)THEN
            K=CURRENT_POS_T(ZONE)-PDEB_SOLVE_Z(ZONE)+1
         ENDIF
         IF(SOLVE_STEP.EQ.0)THEN
            I=CUR_POS_SEQUENCE
            DO WHILE(I.LE.TOTAL_NB_OOC_NODES(OOC_FCT_TYPE))
               IF(SIZE_OF_BLOCK(STEP_OOC(OOC_INODE_SEQUENCE(I,
     $              OOC_FCT_TYPE)),
     $              OOC_FCT_TYPE)
     $              .NE.0)THEN
                  EXIT
               ENDIF
               I=I+1
            ENDDO
            CUR_POS_SEQUENCE=min(I,TOTAL_NB_OOC_NODES(OOC_FCT_TYPE))
            I=CUR_POS_SEQUENCE
            DO WHILE((J.LE.MAX_SIZE).AND.
     $           (I.LE.TOTAL_NB_OOC_NODES(OOC_FCT_TYPE)).AND.
     $           (K.LT.MAX_NB))
               IF((INODE_TO_POS(STEP_OOC(OOC_INODE_SEQUENCE(I,
     $              OOC_FCT_TYPE)))
     $              .NE.0).OR.
     $              (IO_REQ(STEP_OOC(OOC_INODE_SEQUENCE(I,
     $              OOC_FCT_TYPE))).GE.
     $              0))THEN
                  IF(.NOT.ALREADY)THEN
                     CUR_POS_SEQUENCE=CUR_POS_SEQUENCE+1
                     I=I+1
                     CYCLE
                  ELSE
                     EXIT
                  ENDIF
               ENDIF
               ALREADY=.TRUE.
               LAST=SIZE_OF_BLOCK(STEP_OOC(OOC_INODE_SEQUENCE(I,
     $              OOC_FCT_TYPE)),
     $              OOC_FCT_TYPE)
               J=J+LAST
               I=I+1
               K=K+1
               NB_NODES=NB_NODES+1
            ENDDO
            IF(J.GT.MAX_SIZE)THEN
               SIZE=J-LAST
               NB_NODES=NB_NODES-1
            ELSE
               SIZE=J
            ENDIF
            DO WHILE (CUR_POS_SEQUENCE+NB_NODES-1.GT.0)
               IF(SIZE_OF_BLOCK(STEP_OOC(
     $              OOC_INODE_SEQUENCE(CUR_POS_SEQUENCE+NB_NODES-1,
     $              OOC_FCT_TYPE)),
     $              OOC_FCT_TYPE)
     $              .NE.0)THEN
                  EXIT
               ENDIF
               NB_NODES=NB_NODES-1
            ENDDO
            POS_SEQ=CUR_POS_SEQUENCE
         ELSEIF(SOLVE_STEP.EQ.1)THEN
            DO WHILE(I.GE.1)
               IF(SIZE_OF_BLOCK(STEP_OOC(OOC_INODE_SEQUENCE(I,
     $              OOC_FCT_TYPE)),
     $              OOC_FCT_TYPE)
     $              .NE.0)THEN
                  EXIT
               ENDIF
               I=I-1
            ENDDO
            CUR_POS_SEQUENCE=max(I,1)
            I=CUR_POS_SEQUENCE
            DO WHILE((J.LE.MAX_SIZE).AND.(I.GE.1).AND.
     $           (K.LT.MAX_NB))
               IF((INODE_TO_POS(STEP_OOC(OOC_INODE_SEQUENCE(I,
     $              OOC_FCT_TYPE)))
     $              .NE.0).OR.
     $              (IO_REQ(STEP_OOC(OOC_INODE_SEQUENCE(I,
     $              OOC_FCT_TYPE))).GE.
     $              0))THEN
                  IF(.NOT.ALREADY)THEN
                     I=I-1
                     CUR_POS_SEQUENCE=CUR_POS_SEQUENCE-1
                     CYCLE
                  ELSE
                     EXIT
                  ENDIF
               ENDIF
               ALREADY=.TRUE.
               LAST=SIZE_OF_BLOCK(STEP_OOC(OOC_INODE_SEQUENCE(I,
     $              OOC_FCT_TYPE)),
     $              OOC_FCT_TYPE)
               J=J+LAST
               I=I-1
               K=K+1
               NB_NODES=NB_NODES+1
            ENDDO
            IF(J.GT.MAX_SIZE)THEN
               SIZE=J-LAST
               NB_NODES=NB_NODES-1
            ELSE
               SIZE=J
            ENDIF
            I=CUR_POS_SEQUENCE-NB_NODES+1
            DO WHILE (I.GT.0)
               IF(SIZE_OF_BLOCK(STEP_OOC(
     $              OOC_INODE_SEQUENCE(I,OOC_FCT_TYPE)),
     $              OOC_FCT_TYPE).NE.0)THEN
                  EXIT
               ENDIF
               I=I-1
               NB_NODES=NB_NODES-1
            ENDDO
            POS_SEQ=CUR_POS_SEQUENCE-NB_NODES+1
         ENDIF
      ENDIF
      IF(FLAG.EQ.0)THEN
         DEST=IDEB_SOLVE_Z(ZONE)+LRLU_SOLVE_B(ZONE)-SIZE
      ELSE
         DEST=POSFAC_SOLVE(ZONE)
      ENDIF
      END SUBROUTINE DMUMPS_602
      SUBROUTINE DMUMPS_582(IERR)
      IMPLICIT NONE
      INTEGER SOLVE_OR_FACTO
      INTEGER, intent(out) :: IERR
      IERR=0
      IF(allocated(LRLUS_SOLVE))THEN
         DEALLOCATE(LRLUS_SOLVE)
      ENDIF
      IF(allocated(LRLU_SOLVE_T))THEN
         DEALLOCATE(LRLU_SOLVE_T)
      ENDIF
      IF(allocated(LRLU_SOLVE_B))THEN
         DEALLOCATE(LRLU_SOLVE_B)
      ENDIF
      IF(allocated(POSFAC_SOLVE))THEN
         DEALLOCATE(POSFAC_SOLVE)
      ENDIF
      IF(allocated(IDEB_SOLVE_Z))THEN
         DEALLOCATE(IDEB_SOLVE_Z)
      ENDIF
      IF(allocated(PDEB_SOLVE_Z))THEN
         DEALLOCATE(PDEB_SOLVE_Z)
      ENDIF
      IF(allocated(SIZE_SOLVE_Z))THEN
         DEALLOCATE(SIZE_SOLVE_Z)
      ENDIF
      IF(allocated(CURRENT_POS_T))THEN
         DEALLOCATE(CURRENT_POS_T)
      ENDIF
      IF(allocated(CURRENT_POS_B))THEN
         DEALLOCATE(CURRENT_POS_B)
      ENDIF
      IF(allocated(POS_HOLE_T))THEN
         DEALLOCATE(POS_HOLE_T)
      ENDIF
      IF(allocated(POS_HOLE_B))THEN
         DEALLOCATE(POS_HOLE_B)
      ENDIF
      IF(allocated(OOC_STATE_NODE))THEN
         DEALLOCATE(OOC_STATE_NODE)
      ENDIF
      IF(allocated(POS_IN_MEM))THEN
         DEALLOCATE(POS_IN_MEM)
      ENDIF
      IF(allocated(INODE_TO_POS))THEN
         DEALLOCATE(INODE_TO_POS)
      ENDIF
      IF(allocated(IO_REQ))THEN
         DEALLOCATE(IO_REQ)
      ENDIF
      IF(allocated(SIZE_OF_READ))THEN
         DEALLOCATE(SIZE_OF_READ)
      ENDIF
      IF(allocated(FIRST_POS_IN_READ))THEN
         DEALLOCATE(FIRST_POS_IN_READ)
      ENDIF
      IF(allocated(READ_DEST))THEN
         DEALLOCATE(READ_DEST)
      ENDIF
      IF(allocated(READ_MNG))THEN
         DEALLOCATE(READ_MNG)
      ENDIF
      IF(allocated(REQ_TO_ZONE))THEN
         DEALLOCATE(REQ_TO_ZONE)
      ENDIF
      IF(allocated(REQ_ID))THEN
         DEALLOCATE(REQ_ID)
      ENDIF
      SOLVE_OR_FACTO=1
      CALL MUMPS_CLEAN_IO_DATA_C(MYID_OOC,SOLVE_OR_FACTO,IERR)
      IF(IERR.LT.0)THEN
         IF (ICNTL1.GT.0)
     *   WRITE(ICNTL1,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
         RETURN
      ENDIF
      END SUBROUTINE DMUMPS_582
      SUBROUTINE DMUMPS_612(PTRFAC,NSTEPS)
      IMPLICIT NONE
      INTEGER NSTEPS
      INTEGER PTRFAC(NSTEPS)
      INTEGER I,TMP,ZONE,SAVE_PTR
      IF(SOLVE_STEP.EQ.0)THEN
         DO I=1,TOTAL_NB_OOC_NODES(OOC_FCT_TYPE)
            TMP=INODE_TO_POS(STEP_OOC(OOC_INODE_SEQUENCE(I,
     $           OOC_FCT_TYPE)))
            IF(TMP.EQ.0)THEN
               EXIT
            ELSEIF(TMP.LT.0)THEN
               IF(TMP.GT.-(N_OOC+1)*NB_Z)THEN
                  SAVE_PTR=PTRFAC(STEP_OOC(OOC_INODE_SEQUENCE(I,
     $                 OOC_FCT_TYPE)))
                  PTRFAC(STEP_OOC(OOC_INODE_SEQUENCE(I,
     $                 OOC_FCT_TYPE)))=abs(SAVE_PTR)
                  CALL DMUMPS_600(OOC_INODE_SEQUENCE(I,
     $                 OOC_FCT_TYPE),
     $                 ZONE,PTRFAC,NSTEPS)
                  PTRFAC(STEP_OOC(OOC_INODE_SEQUENCE(I,OOC_FCT_TYPE)))=
     $                 SAVE_PTR
                  IF(ZONE.EQ.NB_Z)THEN
                     IF(.NOT.((OOC_INODE_SEQUENCE(I,OOC_FCT_TYPE).EQ.
     $                    KEEP_OOC(38)).AND.(KEEP_OOC(38).NE.0)))THEN
                        EXIT
                     ENDIF
                  ENDIF
                  CALL DMUMPS_599(
     $                 OOC_INODE_SEQUENCE(I,OOC_FCT_TYPE),PTRFAC,NSTEPS)
               ENDIF
            ENDIF
         ENDDO
         CUR_POS_SEQUENCE=I
      ELSE
         DO I=TOTAL_NB_OOC_NODES(OOC_FCT_TYPE),1,-1
            TMP=INODE_TO_POS(STEP_OOC(OOC_INODE_SEQUENCE(I,
     $           OOC_FCT_TYPE)))
            IF(TMP.EQ.0)THEN
               EXIT
            ELSEIF(TMP.LT.0)THEN
               IF(TMP.GT.-(N_OOC+1)*NB_Z)THEN
                  SAVE_PTR=PTRFAC(STEP_OOC(OOC_INODE_SEQUENCE(I,
     $                 OOC_FCT_TYPE)))
                  PTRFAC(STEP_OOC(OOC_INODE_SEQUENCE(I,
     $                 OOC_FCT_TYPE)))=abs(SAVE_PTR)
                  CALL DMUMPS_600(OOC_INODE_SEQUENCE(I,
     $                 OOC_FCT_TYPE)
     $                 ,ZONE,
     $                 PTRFAC,NSTEPS)
                  PTRFAC(STEP_OOC(OOC_INODE_SEQUENCE(I,
     $                 OOC_FCT_TYPE)))=SAVE_PTR
                  IF(ZONE.EQ.NB_Z)THEN
                     IF(.NOT.((OOC_INODE_SEQUENCE(I,
     $                    OOC_FCT_TYPE).EQ.KEEP_OOC(38))
     $                    .AND.(KEEP_OOC(38).NE.0)))THEN
                        EXIT
                     ENDIF
                  ENDIF
                  CALL DMUMPS_599(
     $                 OOC_INODE_SEQUENCE(I,OOC_FCT_TYPE),PTRFAC,NSTEPS)
               ENDIF
            ENDIF
         ENDDO
         CUR_POS_SEQUENCE=I
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_612
      SUBROUTINE DMUMPS_583(PTRFAC,NSTEPS,MTYPE,
     *                                    A,LA,DOPREFETCH,IERR)
      IMPLICIT NONE
      INTEGER NSTEPS,MTYPE
      INTEGER, intent(out)::IERR
      INTEGER LA 
      DOUBLE PRECISION A(LA)
      INTEGER PTRFAC(NSTEPS)
      LOGICAL DOPREFETCH
      IERR = 0
      OOC_FCT_TYPE=1
      OOC_SOLVE_TYPE_FCT=FCT
      IF (KEEP_OOC(201).EQ.1) THEN
       OOC_FCT_TYPE=1
       OOC_SOLVE_TYPE_FCT=TYPEF_L-1
      ENDIF
      IF (KEEP_OOC(201).EQ.1) THEN
       IF((MTYPE.NE.1).AND.(KEEP_OOC(50).EQ.0))THEN
         OOC_SOLVE_TYPE_FCT=TYPEF_U-1
         OOC_FCT_TYPE=2
       ENDIF
      ENDIF
      SOLVE_STEP=0
      CUR_POS_SEQUENCE=1
      MTYPE_OOC=MTYPE
      IF ( KEEP_OOC(201).NE.1 
#if ! defined(TEMP_VERSION_TO_FORCE_DATA_REINIT)
     *  .OR. KEEP_OOC(50).NE.0
#endif
     *  ) THEN
        CALL DMUMPS_612(PTRFAC,NSTEPS)
      ELSE 
        CALL DMUMPS_683(KEEP_OOC(28),
     $              KEEP_OOC(38), KEEP_OOC(20) )
      ENDIF
      IF (DOPREFETCH) THEN 
          CALL DMUMPS_585(A,LA,PTRFAC,
     *                                 KEEP_OOC(28),IERR)
      ELSE
          CUR_POS_SEQUENCE = TOTAL_NB_OOC_NODES(OOC_FCT_TYPE)
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_583
      SUBROUTINE DMUMPS_584(PTRFAC,NSTEPS,MTYPE,
     *           I_WORKED_ON_ROOT,IROOT,A,LA,IERR)
      IMPLICIT NONE
      INTEGER NSTEPS,LA
      INTEGER PTRFAC(NSTEPS)
      INTEGER MTYPE
      INTEGER IROOT
      LOGICAL I_WORKED_ON_ROOT
      INTEGER, intent(out):: IERR
      DOUBLE PRECISION A(LA)
      IERR=0
      OOC_FCT_TYPE=1
      OOC_SOLVE_TYPE_FCT=FCT
      IF (KEEP_OOC(201).EQ.1) THEN
       IF(KEEP_OOC(50).EQ.0)THEN
         OOC_SOLVE_TYPE_FCT=TYPEF_U-1
         OOC_FCT_TYPE=2
       ELSE
         OOC_SOLVE_TYPE_FCT=TYPEF_L-1
         OOC_FCT_TYPE=1
       ENDIF
      ENDIF
      IF (KEEP_OOC(201).EQ.1) THEN
       IF((MTYPE.NE.1).AND.(KEEP_OOC(50).EQ.0))THEN
         OOC_FCT_TYPE=1
         OOC_SOLVE_TYPE_FCT=TYPEF_L-1
       ENDIF
      ENDIF
      SOLVE_STEP=1
      CUR_POS_SEQUENCE=TOTAL_NB_OOC_NODES(OOC_FCT_TYPE)
      MTYPE_OOC=MTYPE
      IF ( KEEP_OOC(201).NE.1 
#if ! defined(TEMP_VERSION_TO_FORCE_DATA_REINIT)
     *  .OR. KEEP_OOC(50).NE.0
#endif
     *  ) THEN
        CALL DMUMPS_612(PTRFAC,NSTEPS)
        IF (I_WORKED_ON_ROOT) THEN
          CALL DMUMPS_598( IROOT,
     *         PTRFAC, KEEP_OOC(28), A, LA,.FALSE.,IERR)
          IF (IERR .LT. 0) RETURN
        ENDIF
        IF (NB_Z.GT.1) THEN
          CALL DMUMPS_594(A,LA,PTRFAC,
     *                                  KEEP_OOC(28),IERR)
          IF (IERR .LT. 0) RETURN
        ENDIF
      ELSE 
        CALL DMUMPS_683(KEEP_OOC(28),
     $              KEEP_OOC(38), KEEP_OOC(20) )
        CALL DMUMPS_585(A,LA,PTRFAC,KEEP_OOC(28),IERR)
        IF (IERR .LT. 0 ) RETURN
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_584
      SUBROUTINE DMUMPS_613(id,IERR)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC), TARGET :: id
      INTEGER, intent(out) :: IERR
      INTEGER I,DIM,J,TMP,SIZE,K,I1
      CHARACTER*1 TMP_NAME(350)
      EXTERNAL MUMPS_OOC_GET_NB_FILES_C, MUMPS_OOC_GET_FILE_NAME_C
      IERR=0
      SIZE=0
      DO J=1,OOC_NB_FILE_TYPE
         TMP=J-1
         CALL MUMPS_OOC_GET_NB_FILES_C(TMP,I)
         id%OOC_NB_FILES(J)=I
         SIZE=SIZE+I
      ENDDO
      IF(associated(id%OOC_FILE_NAMES))THEN
         DEALLOCATE(id%OOC_FILE_NAMES)
         NULLIFY(id%OOC_FILE_NAMES)
      ENDIF
      ALLOCATE(id%OOC_FILE_NAMES(SIZE,150),stat=IERR)
      IF (IERR .GT. 0) THEN
         IF (ICNTL1.GT.0)
     *   WRITE(ICNTL1,*) 'PB allocation in DMUMPS_613'
         IERR=-1
         IF(id%INFO(1).GE.0)THEN
            id%INFO(1) = -13
            id%INFO(2) = SIZE*150
            RETURN
         ENDIF
      ENDIF
      IF(associated(id%OOC_FILE_NAME_LENGTH))THEN
         DEALLOCATE(id%OOC_FILE_NAME_LENGTH)
         NULLIFY(id%OOC_FILE_NAME_LENGTH)
      ENDIF
      ALLOCATE(id%OOC_FILE_NAME_LENGTH(SIZE),stat=IERR)
      IF (IERR .GT. 0) THEN
         IERR=-1
         IF(id%INFO(1).GE.0)THEN
            IF (ICNTL1.GT.0)
     &      WRITE(ICNTL1,*)
     &      'PB allocation in DMUMPS_613'
            id%INFO(1) = -13
            id%INFO(2) = SIZE
            RETURN
         ENDIF
      ENDIF
      K=1
      DO I1=1,OOC_NB_FILE_TYPE
         TMP=I1-1
         DO I=1,id%OOC_NB_FILES(I1)
            CALL MUMPS_OOC_GET_FILE_NAME_C(TMP,I,DIM,TMP_NAME(1))
            DO J=1,DIM+1
               id%OOC_FILE_NAMES(K,J)=TMP_NAME(J)
            ENDDO
            id%OOC_FILE_NAME_LENGTH(K)=DIM+1
            K=K+1
         ENDDO
      ENDDO
      END SUBROUTINE DMUMPS_613
      SUBROUTINE DMUMPS_614(id)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC), TARGET :: id
      CHARACTER*1 TMP_NAME(350)
      INTEGER I,I1,TMP,J,K,L,DIM,IERR
      INTEGER, DIMENSION(:),ALLOCATABLE :: NB_FILES
      INTEGER K211
      ALLOCATE(NB_FILES(OOC_NB_FILE_TYPE),stat=IERR)
      IF (IERR .GT. 0) THEN
         IERR=-1
         IF(id%INFO(1).GE.0)THEN
            IF (ICNTL1.GT.0)
     *      WRITE(ICNTL1,*)
     *      'PB allocation in DMUMPS_614'
            id%INFO(1) = -13
            id%INFO(2) = OOC_NB_FILE_TYPE
            RETURN
         ENDIF
      ENDIF
      IERR=0
      NB_FILES=id%OOC_NB_FILES
      I=id%MYID
      K=id%KEEP(35)
      L=mod(id%KEEP(204),3)
      K211=id%KEEP(211)
      CALL MUMPS_OOC_ALLOC_POINTERS_C(OOC_NB_FILE_TYPE,NB_FILES,IERR)
      IF(IERR.LT.0)THEN
         IF (ICNTL1.GT.0)
     *   WRITE(ICNTL1,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
         id%INFO(1)=IERR
         RETURN
      ENDIF
      CALL MUMPS_OOC_INIT_VARS_C(I,K,L,K211,IERR)
      IF(IERR.LT.0)THEN
         IF (ICNTL1.GT.0)
     *   WRITE(ICNTL1,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
         id%INFO(1)=IERR
         RETURN
      ENDIF
      K=1
      DO I1=1,OOC_NB_FILE_TYPE
         DO I=1,NB_FILES(I1)
            DIM=id%OOC_FILE_NAME_LENGTH(K)
            DO J=1,DIM
               TMP_NAME(J)=id%OOC_FILE_NAMES(K,J)
            ENDDO
            TMP=I1-1
            CALL MUMPS_OOC_SET_FILE_NAME_C(TMP,I,DIM,IERR,TMP_NAME(1))
            IF(IERR.LT.0)THEN
               IF (ICNTL1.GT.0)
     *         WRITE(ICNTL1,*)MYID_OOC,': ',
     *         ERR_STR_OOC(1:DIM_ERR_STR_OOC)
               id%INFO(1)=IERR
               RETURN
            ENDIF
            K=K+1
         ENDDO
      ENDDO
      CALL MUMPS_OOC_START_LOW_LEVEL(IERR)
      IF(IERR.LT.0)THEN
         IF (ICNTL1.GT.0)
     *   WRITE(ICNTL1,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
         id%INFO(1)=IERR
         RETURN
      ENDIF
      DEALLOCATE(NB_FILES)
      RETURN
      END SUBROUTINE DMUMPS_614      
      SUBROUTINE DMUMPS_589(DEST,SRC,NB,NB_EFF)
      IMPLICIT NONE
      INTEGER NB, NB_EFF
      CHARACTER(LEN=NB) SRC
      CHARACTER*1 DEST(NB)
      INTEGER I
      DO I=1,NB_EFF
         DEST(I)=SRC(I:I)
      ENDDO
      END SUBROUTINE DMUMPS_589
      SUBROUTINE DMUMPS_580(IERR)
      USE DMUMPS_OOC_BUFFER
      IMPLICIT NONE
      INTEGER, intent(out) :: IERR
      IERR=0
      IF(.NOT.WITH_BUF)THEN
         RETURN
      ENDIF
      CALL DMUMPS_707(OOC_FCT_TYPE,IERR)
      IF (IERR < 0) THEN
        RETURN
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_580
      SUBROUTINE DMUMPS_681(IERR)
      USE DMUMPS_OOC_BUFFER
      IMPLICIT NONE
      INTEGER, intent(out) :: IERR
      IERR=0
      IF(.NOT.WITH_BUF)THEN
         RETURN
      ENDIF
      CALL DMUMPS_707(TYPEF_L_LOC,IERR)
      IF (IERR < 0) RETURN
      IF (KEEP_OOC(50).EQ.0) THEN
         CALL DMUMPS_707(TYPEF_U_LOC,IERR)
         IF (IERR < 0) RETURN
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_681
       SUBROUTINE DMUMPS_683(NSTEPS,
     $     KEEP38, KEEP20)
      IMPLICIT NONE
      INTEGER MASTER
      PARAMETER(MASTER=0)
      INTEGER NSTEPS
      INTEGER I, TMPNODE, TMP, TMP_SIZE, ZONE, J
      INTEGER KEEP38, KEEP20
      INODE_TO_POS = 0
      POS_IN_MEM   = 0
      OOC_STATE_NODE(1:NSTEPS)=0
      FIRST_ACTIVE_Z=1
      TMP_SIZE=1
      J=1
      DO I=1,NB_Z-1
         IDEB_SOLVE_Z(I)=TMP_SIZE
         PDEB_SOLVE_Z(I)=J
         POSFAC_SOLVE(I)=TMP_SIZE
         LRLUS_SOLVE(I) =SIZE_ZONE_SOLVE
         LRLU_SOLVE_T(I)=SIZE_ZONE_SOLVE
         LRLU_SOLVE_B(I)=0
         SIZE_SOLVE_Z(I)=SIZE_ZONE_SOLVE
         CURRENT_POS_T(I)=J
         CURRENT_POS_B(I)=J
         POS_HOLE_T(I)   =J
         POS_HOLE_B(I)   =J
         J = J + MAX_NB_NODES_FOR_ZONE
         TMP_SIZE = TMP_SIZE + SIZE_ZONE_SOLVE
      ENDDO
      IDEB_SOLVE_Z(NB_Z)=TMP_SIZE
      PDEB_SOLVE_Z(NB_Z)=J
      POSFAC_SOLVE(NB_Z)=TMP_SIZE
      LRLUS_SOLVE(NB_Z) =SIZE_SOLVE_EMM
      LRLU_SOLVE_T(NB_Z)=SIZE_SOLVE_EMM
      LRLU_SOLVE_B(NB_Z)=0
      SIZE_SOLVE_Z(NB_Z)=SIZE_SOLVE_EMM
      CURRENT_POS_T(NB_Z)=J
      CURRENT_POS_B(NB_Z)=J
      POS_HOLE_T(NB_Z)   =J
      POS_HOLE_B(NB_Z)   =J
      IO_REQ=-77777
      SIZE_OF_READ=-9999
      FIRST_POS_IN_READ=-9999
      READ_DEST=-9999
      READ_MNG=-9999
      REQ_TO_ZONE=-9999
      REQ_ID=-9999
      END SUBROUTINE DMUMPS_683
      SUBROUTINE DMUMPS_688 
     &     ( STRAT, TYPEFile, 
     &     AFAC, LAFAC, MonBloc,
     &     LNextPiv2beWritten, UNextPiv2beWritten,
     &     IW, LIWFAC,
     &     MYID, FILESIZE, IERR )
      IMPLICIT NONE
      TYPE(IO_BLOCK), INTENT(INOUT):: MonBloc   
      INTEGER,        INTENT(IN)   :: STRAT, LIWFAC, LAFAC, 
     &     MYID, TYPEFile
      INTEGER,   INTENT(INOUT)        :: IW(0:LIWFAC-1) 
      DOUBLE PRECISION, INTENT(IN) :: AFAC(LAFAC)
      INTEGER,   INTENT(INOUT) :: LNextPiv2beWritten, 
     &     UNextPiv2beWritten
      INTEGER*8, INTENT(INOUT) :: FILESIZE
      INTEGER,   INTENT(OUT) :: IERR
      INTEGER L_AddVirtNode(2), L_TailleEcrite
      INTEGER U_AddVirtNode(2), U_TailleEcrite, TempFTYPE
      LOGICAL FIRSTL, FIRSTU, WRITE_L, WRITE_U
      LOGICAL DO_U_FIRST
      INCLUDE 'mumps_headers.h'
      IERR = 0                  
      WRITE_L =  (TYPEFile.EQ.TYPEF_BOTH_LU .OR. TYPEFile.EQ.TYPEF_L)
      WRITE_U =  (TYPEFile.EQ.TYPEF_BOTH_LU .OR. TYPEFile.EQ.TYPEF_U)
      DO_U_FIRST = .FALSE.
      IF ( TYPEFile.EQ.TYPEF_BOTH_LU ) THEN
         IF ( LNextPiv2bewritten .GT. UNextPiv2beWritten ) THEN
            DO_U_FIRST = .TRUE.
         END IF
      END IF
      IF (DO_U_FIRST) GOTO 200
 100  IF (WRITE_L) THEN
         L_TailleEcrite  = IW(XXWRITTENL)
         FIRSTL = (L_TailleEcrite.EQ.0)
         IF (L_TailleEcrite.LT.0) THEN
            WRITE(6,*) ' Internal error in DMUMPS_688 ',
     &           L_TailleEcrite, ' bad default value in header '
            CALL MUMPS_ABORT()
         ENDIF
         IF (.NOT.FIRSTL) THEN
            L_AddVirtNode (1) = max(IW(XXVADDRL),0)
            L_AddVirtNode (2) = max(IW(XXVADDRL+1),0)
         ENDIF
         IF ((MonBloc%Typenode.EQ.2).AND.(.NOT.MonBloc%MASTER)) 
     &        THEN
            IF (L_TailleEcrite.GT.0) THEN 
               LNextPiv2beWritten = L_TailleEcrite / MonBloc%NROW + 1
            ELSE
               LNextPiv2beWritten = 1 
            ENDIF
         ENDIF
         IF (MonBloc%LastPiv.GE.LNextPiv2beWritten
     &       .OR. MonBloc%Last) THEN
           TempFTYPE    = TYPEF_L
           CALL DMUMPS_695( STRAT, FIRSTL,
     &        TempFTYPE, AFAC, LAFAC, MonBloc,
     &        IERR,
     &        LNextPiv2beWritten, 
     &        L_AddVirtNode, L_TailleEcrite) 
           IF (IERR .LT. 0) RETURN
           IF (L_TailleEcrite.GT.0) THEN
             IW(XXVADDRL)   = L_AddVirtNode (1)
             IW(XXVADDRL+1) = L_AddVirtNode (2)
             IW(XXWRITTENL) = L_TailleEcrite
           ENDIF
         ENDIF
         IF (STRAT.EQ.STRAT_WRITE_MAX .AND. MonBloc%Last) THEN
            FILESIZE = FILESIZE + L_TailleEcrite
         ENDIF
         IF (DO_U_FIRST) GOTO 300
      ENDIF
 200  IF (WRITE_U) THEN
         U_TailleEcrite  = IW(XXWRITTENU)
         FIRSTU = (U_TailleEcrite.EQ.0)
         IF (U_TailleEcrite.LT.0) THEN
            WRITE(6,*) ' Internal error in DMUMPS_688 ',
     &           U_TailleEcrite, ' bad default value in header '
            CALL MUMPS_ABORT()
         ENDIF
         IF (.NOT.FIRSTU) THEN
            U_AddVirtNode (1) = max(IW(XXVADDRU),0)
            U_AddVirtNode (2) = max(IW(XXVADDRU+1),0)
         ENDIF
         IF (MonBloc%LastPiv.GE.UNextPiv2beWritten .OR.
     &       MonBloc%Last) THEN
           TempFTYPE  = TYPEF_U
           CALL DMUMPS_695( STRAT, FIRSTU,
     &        TempFTYPE, AFAC, LAFAC, MonBloc,
     &        IERR,
     &        UNextPiv2beWritten, 
     &        U_AddVirtNode, 
     &        U_TailleEcrite) 
            IF (IERR .LT. 0) RETURN
            IF (U_TailleEcrite.GT.0) THEN
              IW(XXVADDRU)   = U_AddVirtNode (1)
              IW(XXVADDRU+1) = U_AddVirtNode (2)
              IW(XXWRITTENU) = U_TailleEcrite
            ENDIF
         ENDIF
         IF ( STRAT.EQ.STRAT_WRITE_MAX .AND. MonBloc%last ) THEN
           FILESIZE = FILESIZE + U_TailleEcrite
         ENDIF
         IF (DO_U_FIRST) GOTO 100
      ENDIF
 300  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_688 
      SUBROUTINE DMUMPS_695( STRAT, FIRST, TYPEF, 
     &     AFAC, LAFAC, MonBloc,
     &     IERR,
     &     LorU_NextPiv2beWritten, 
     &     LorU_AddVirtNode, LorU_TailleEcrite
     &     )
      USE DMUMPS_OOC_BUFFER
      IMPLICIT NONE
      LOGICAL FIRST
      INTEGER, INTENT(IN) :: STRAT
      INTEGER, INTENT(IN) :: TYPEF, LAFAC
      DOUBLE PRECISION, INTENT(IN) :: AFAC(LAFAC)
      INTEGER, INTENT(INOUT) :: LorU_NextPiv2beWritten
      INTEGER, INTENT(INOUT) :: LorU_AddVirtNode(2)
      INTEGER, INTENT(INOUT) :: LorU_TailleEcrite
      TYPE(IO_BLOCK), INTENT(INOUT) :: MonBloc   
      INTEGER, INTENT(OUT)  :: IERR
      INTEGER NNMAX             
      INTEGER TOTSIZE, UnChouia, TailleEcrite, SIZE_PANEL
      INTEGER*8 AddVirtNode, AddVirtCour
      IERR = 0
      IF (TYPEF == TYPEF_L ) THEN
         NNMAX = MonBloc%NROW   
      ELSE
         NNMAX = MonBloc%NCOL
      ENDIF
      SIZE_PANEL = DMUMPS_OOC_PANEL_SIZE(NNMAX) 
      IF ( (.NOT.MonBloc%Last) .AND.
     &     (MonBloc%LastPiv-LorU_NextPiv2beWritten+1.LT.SIZE_PANEL))
     &     THEN
      RETURN
      ENDIF
      CALL MUMPS_676(
     &     LorU_AddVirtNode(1), LorU_AddVirtNode(2), 
     &     AddVirtNode)
      IF (FIRST) THEN
         AddVirtNode = AddVirtLibre(TYPEF)
         CALL MUMPS_677(
     &        LorU_AddVirtNode(1), LorU_AddVirtNode(2), 
     &        AddVirtNode )
         IF (MonBloc%Typenode.eq.3 .AND. TYPEF == TYPEF_L) THEN
            TOTSIZE = MonBloc%NCOL
         ELSE
            TOTSIZE = DMUMPS_OOC_ESTIM_FRONTSIZE_12
     &           (MonBloc%NFS, NNMAX,
     &           MonBloc%MASTER, SIZE_PANEL,
     &           MonBloc%Typenode)
         ENDIF
         AddVirtLibre(TYPEF) = AddVirtLibre(TYPEF) + TOTSIZE
      ELSE
      ENDIF
      AddVirtCour = AddVirtNode + LorU_TailleEcrite
      CALL DMUMPS_697( STRAT, TYPEF, MonBloc,
     &     SIZE_PANEL, 
     &     AFAC, LAFAC,
     &     LorU_NextPiv2beWritten, AddVirtCour, 
     &     TailleEcrite, 
     &     IERR )
      IF ( IERR .LT. 0 ) RETURN
      LorU_TailleEcrite = LorU_TailleEcrite + TailleEcrite
      IF (LorU_TailleEcrite.EQ.0) THEN
         AddVirtLibre(TYPEF) = AddVirtLibre(TYPEF) - TOTSIZE
         LorU_AddVirtNode(1:2) = 0
         RETURN
      ELSE IF (FIRST) THEN
         OOC_VADDR(STEP_OOC(MonBloc%INODE),TYPEF)= AddVirtNode
         SIZE_OF_BLOCK(STEP_OOC(MonBloc%INODE),TYPEF) =
     $        TOTSIZE
         MAX_SIZE_FACTOR_OOC=max(MAX_SIZE_FACTOR_OOC,TOTSIZE)
         TMP_SIZE_FACT=TMP_SIZE_FACT+TOTSIZE
         TMP_NB_NODES=TMP_NB_NODES+1
         IF(TMP_SIZE_FACT.GT.SIZE_ZONE_SOLVE)THEN
            MAX_NB_NODES_FOR_ZONE=max(MAX_NB_NODES_FOR_ZONE,
     $           TMP_NB_NODES)
            TMP_SIZE_FACT=0
            TMP_NB_NODES=0
         ENDIF
         OOC_INODE_SEQUENCE(I_CUR_HBUF_NEXTPOS(TYPEF),
     $        TYPEF) = MonBloc%INODE
         I_CUR_HBUF_NEXTPOS(TYPEF) = I_CUR_HBUF_NEXTPOS(TYPEF) + 1
      ENDIF
      IF (MonBloc%Last) THEN
         KEEP_OOC(228) = max(KEEP_OOC(228),
     &        (MonBloc%LastPiv+SIZE_PANEL-1) / SIZE_PANEL)
      ENDIF
      IF (MonBloc%Last .AND. MonBloc%LastPiv .LT. LorU_NextPiv2beWritten
     &   .AND. MonBloc%Typenode.NE.3.AND.
     &     (
     &     (MonBloc%LastPiv.NE.MonBloc%NFS).OR. 
     &     (KEEP_OOC(50).EQ.2)  
     &     ) 
     &     ) THEN
         TOTSIZE = DMUMPS_OOC_ESTIM_FRONTSIZE_12
     &        (MonBloc%NFS, MonBloc%NROW, 
     &        MonBloc%MASTER, SIZE_PANEL,
     &        MonBloc%Typenode)
         IF (AddVirtLibre(TYPEF).EQ.
     &        (AddVirtNode+TOTSIZE) ) THEN
            IF (TOTSIZE < LorU_TailleEcrite) THEN
               IERR = -5
               write(6,*) ' ... Internal error', IERR, 
     &              ' in DMUMPS_695 : '
               write(6,*) ' ... TOTSIZE, TYPEF,LorU_TailleEcrite ', 
     &              TOTSIZE, TYPEF, LorU_TailleEcrite
               CALL MUMPS_ABORT()
            ENDIF
            AddVirtLibre(TYPEF) = AddVirtNode+LorU_TailleEcrite
            SIZE_OF_BLOCK(STEP_OOC(MonBloc%INODE),TYPEF)
     *           = LorU_TailleEcrite
         ENDIF
      ENDIF
      IF (MonBloc%Last .AND. STRAT.EQ.STRAT_WRITE_MAX) THEN
         MAX_SIZE_FACTOR_OOC=max(MAX_SIZE_FACTOR_OOC,
     *        SIZE_OF_BLOCK(STEP_OOC(MonBloc%INODE),TYPEF))
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_695
      SUBROUTINE DMUMPS_697(
     &     STRAT, TYPEF, MonBloc,
     &     SIZE_PANEL,
     &     AFAC, LAFAC,
     &     NextPiv2beWritten, AddVirtCour, 
     &     TailleEcrite, IERR )
      USE DMUMPS_OOC_BUFFER
      IMPLICIT NONE
      INTEGER,       INTENT(IN) :: STRAT, TYPEF, LAFAC, SIZE_PANEL
      INTEGER*8,     INTENT(IN) ::  AddVirtCour
      DOUBLE PRECISION, INTENT(IN) :: AFAC(LAFAC)
      INTEGER,       INTENT(INOUT) :: NextPiv2beWritten
      TYPE(IO_BLOCK),INTENT(INOUT) :: MonBloc   
      INTEGER,         INTENT(OUT) :: TailleEcrite
      INTEGER, INTENT(OUT)  :: IERR
      INTEGER   :: I, NBeff, LPANELeff, IEND
      INTEGER*8 :: AddVirtDeb
      IERR = 0
      TailleEcrite = 0
      AddVirtDeb   = AddVirtCour
      I = NextPiv2beWritten
      IF ( NextPiv2beWritten .GT. MonBloc%LastPiv ) THEN
        RETURN
      ENDIF
 10   CONTINUE
      NBeff  = min(SIZE_PANEL,MonBloc%LastPiv-I+1 )
      IF ((NBeff.NE.SIZE_PANEL) .AND. (.NOT.MonBloc%Last)) THEN
         GOTO 20
      ENDIF
      IF (TYPEF.EQ.TYPEF_L.AND.MonBloc%MASTER.AND.
     *     KEEP_OOC(50).EQ.2 .AND. MonBloc%Typenode.NE.3) THEN
         IF (MonBloc%INDICES(NBeff+I-1) < 0)
     *        THEN
            NBeff=NBeff+1
         ENDIF
      ENDIF
      IEND   = I + NBeff -1
#if defined(check_vaddr)
      WRITE(*,*) MYID_OOC,": Writing to buffer",AddVirtDeb
#endif
      CALL DMUMPS_653( STRAT, TYPEF, MonBloc,
     &     AFAC, LAFAC,
     &     AddVirtDeb, I, IEND, LPANELeff,
     &     IERR)
      IF ( IERR .LT. 0 ) THEN
        RETURN
      ENDIF
      IF ( IERR .EQ. 1 ) THEN
         IERR=0
         GOTO 20
      ENDIF
      IF (TYPEF .EQ. TYPEF_L) THEN
         MonBloc%LastPanelWritten_L = MonBloc%LastPanelWritten_L+1
      ELSE
         MonBloc%LastPanelWritten_U = MonBloc%LastPanelWritten_U+1
      ENDIF
      AddVirtDeb   = AddVirtDeb + LPANELeff
      TailleEcrite = TailleEcrite + LPANELeff
      I=I+NBeff
      IF ( I .LE. MonBloc%LastPiv ) GOTO 10
 20   CONTINUE
      NextPiv2beWritten = I
      RETURN
      END SUBROUTINE DMUMPS_697
      INTEGER FUNCTION DMUMPS_OOC_ESTIM_FRONTSIZE_12
     &     (NFS, NNMAX, MASTER, SIZE_PANEL, TYPENODE)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NFS, NNMAX, SIZE_PANEL
      LOGICAL, INTENT(IN) :: MASTER
      INTEGER, INTENT(IN) :: TYPENODE
      INTEGER :: I, TOTSIZE
      TOTSIZE = 0
      IF (TYPENODE.NE.1.AND.TYPENODE.NE.2) THEN
         WRITE(*,*)
     &        "Internal error in DMUMPS_OOC_ESTIM_FRONTSIZE_12",
     &        TYPENODE
      ENDIF
      IF ((KEEP_OOC(50).EQ.2).AND.MASTER) THEN
         DO I=1, NFS, SIZE_PANEL
            TOTSIZE = TOTSIZE + 
     &           (NNMAX-I+1) *
     &           ( min(SIZE_PANEL,NFS-I+1 )  + 1   )
         ENDDO
      ELSE IF (.NOT. MASTER) THEN
         TOTSIZE = NFS * NNMAX 
      ELSE
         DO I=1, NFS, SIZE_PANEL
            TOTSIZE = TOTSIZE + 
     &           (NNMAX-I+1) *
     &           min(SIZE_PANEL,NFS-I+1 ) 
         ENDDO
      ENDIF
      DMUMPS_OOC_ESTIM_FRONTSIZE_12 = TOTSIZE
      RETURN
      END FUNCTION DMUMPS_OOC_ESTIM_FRONTSIZE_12
      SUBROUTINE DMUMPS_690(
     &     INODE, TYPEF, IW, LIWFAC, NNMAX, PSIZE, 
     &     TailleEcrite )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  INODE, LIWFAC, NNMAX, TYPEF
      INTEGER, INTENT(OUT) ::  PSIZE, TailleEcrite
      INTEGER, INTENT(IN) ::  IW(0:LIWFAC) 
      INTEGER TOTSIZE
      INCLUDE 'mumps_headers.h'
      PSIZE = DMUMPS_OOC_PANEL_SIZE( NNMAX )
      IF (TYPEF.EQ.TYPEF_L) TailleEcrite = IW(XXWRITTENL)
      IF (TYPEF.EQ.TYPEF_U) TailleEcrite = IW(XXWRITTENU)
      TOTSIZE =  SIZE_OF_BLOCK(STEP_OOC(INODE),TYPEF)
      RETURN
      END SUBROUTINE DMUMPS_690
      INTEGER FUNCTION DMUMPS_OOC_PANEL_SIZE( NNMAX )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NNMAX
      INTEGER DMUMPS_OOC_GET_PANEL_SIZE
      DMUMPS_OOC_PANEL_SIZE=DMUMPS_OOC_GET_PANEL_SIZE(
     *     HBUF_SIZE, NNMAX, KEEP_OOC(227),KEEP_OOC(50))
      RETURN
      END FUNCTION DMUMPS_OOC_PANEL_SIZE
      END MODULE DMUMPS_OOC
