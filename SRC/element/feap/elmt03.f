      subroutine ELMT03(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)
      implicit none

c.... common declarations
      integer         numnp,numel,nummat,nen,neq,ipr
      common /cdata/  numnp,numel,nummat,nen,neq,ipr

      real*8          dm
      integer            n,ma,mct,iel,nel
      common /eldata/ dm,n,ma,mct,iel,nel
      
      integer         ia      ,ir      ,ea      ,er
      common /mdata/  ia(2,50),ir(2,50),ea(2,50),er(2,50)
      
      integer         ior,iow,ilg
      common /iofile/ ior,iow,ilg
      
      integer         nh1,nh2,nh3
      common /hdata/  nh1,nh2,nh3

      real*8  hr
      common  hr(10000)
      
c ... subroutine arguments
      real*8  d(*), xl(ndm,*), ul(ndf,*), s(nst,*), p(nst), tl(*)
      integer ix(*), ndf, ndm, nst, isw, tdof, i, j

c     Input material properties

      if(isw.eq.1) then

C        if(ior.lt.0) write(*,2000)
C        write(iow,2000)
Cc       call inmate(d,tdof,ndf*4,5)
C        call inmate(d,tdof, 24  ,5)
C
Cc       Deactivate dof in element for dof > 6
C
C        do i = 7,ndf
C          ix(i) = 0
C        end do ! i
C
Cc       History for through thickness integrations
C
C        if(nint(d(102)).gt.1) then
C          nh1 = nh1*nint(d(102))
C        endif
C
Cc       Construct rotation parameters: u-x = 1; u-y = 2 (same as defaults)
C
C        ea(1,-iel) = 1
C        ea(2,-iel) = 2
C
Cc       Construct rotation parameters: theta-x = 4; theta-y = 5
C
C        er(1,-iel) = 4
C        er(2,-iel) = 5
C
Cc       Set plot sequence
C
C        call plqud4(iel)
C
Cc       Set maximum number of stress plots
C
CC        istv = 24
C
Cc       Initialize the finite deformation shell
C
C        if(d(18).lt.0.0d0) then
C          call shl3df(d,ul,xl,s,p,ndf,ndm,nst,isw)
C        endif

c     Initialize element strains for activation

      elseif(isw.eq.17) then

        j = nh3
        do i = 1,nel
          hr(j  ) = ul(1,i)
          hr(j+1) = ul(2,i)
          hr(j+2) = ul(3,i)
          hr(j+3) = ul(4,i)
          hr(j+4) = ul(5,i)
          hr(j+5) = ul(6,i)
          j       = j + 6
        end do ! i

c     Initialize element strains for deactivation

      elseif(isw.eq.18) then

        do i = 0,6*nel-1
          hr(nh3+i) = 0.0d0
        end do ! i

c     External node check

      elseif(isw.eq.26) then

        call pcorner2d()

C     Remaining options

      else

Cc       Small deformation
C
C        if(d(18).gt.0.0d0) then
C          if(nint(d(102)).gt.1) then
C            call shl3di(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)
C          else
C            call shl3ds(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)
C          endif
C
Cc       Finite deformation
C
C        else
          call shl3df(d,ul,xl,s,p,ndf,ndm,nst,isw)
C        endif

      endif

c     Format

2000  format(5x,'T h r e e   D i m e n s i o n a l   S h e l l',
     &          '   E l e m e n t')

      end

      
      subroutine pcorner2d()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute corner nodes

c     Inputs:
c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      real*8          dm
      integer            n,ma,mct,iel,nel
      common /eldata/ dm,n,ma,mct,iel,nel
      
      integer          ncorner,icorner
      common /corner/  ncorner,icorner(24,4)

      integer    jcorner3(2,6),jcorner4(2,8),  nn

      data       jcorner3 / 1,3, 1,2,  2,1, 2,3, 3,2, 3,1 /
      data       jcorner4 / 1,4, 1,2,  2,1, 2,3, 3,2, 3,4, 4,3, 4,1 /

      if(nel.eq.3 .or. nel.eq.6) then
        ncorner = 6
        do nn = 1,ncorner
          icorner(nn, 1) = jcorner3(1,nn) 
          icorner(nn, 2) = jcorner3(2,nn) 
        end do ! nn
      else
        ncorner = 8
        do nn = 1,ncorner
          icorner(nn, 1) = jcorner4(1,nn) 
          icorner(nn, 2) = jcorner4(2,nn) 
        end do ! nn
      endif

      end

      
      subroutine shl3df ( d, ul, xl, s, p, ndf, ndm, nst, isw )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c       Description:    SHL3DF is a 4-node bilinear element for
c                       static/transient finite deformation analysis
c                       of inextensible 3-dimensional elastic shell
c                       problems.

c       Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c       Date:           January 1991.

c       Revised:        R.L. Taylor  - - february 1997

c       Modification:   J.C. Simo and N. Tarnow to handle
c                       Energy-Momentum method and shell-intersections.
c                       January 1992.
c                       R.L. Taylor, Adapted to version 6.0 and later.
c                       february 1997.

c       features:       Static and transient inextensible analysis.
c                       Full Isoparametric interpolation, including
c                             inextensible director field.
c                       2-DOF director increments (Delta-T).
c                       Bathe-Dvorkin shear treatment.
c                       Displacement version or
c                       Pian-Summihara membrane/bending treatment.

c       Caution:        This element cannot be used with rotational
c                       increments in excess of 180-degrees.

c       References:     J.C. Simo, D.D. Fox, M.S. Rifai,
c                       "On a Stress Resultant Geometrically Exact
c                       Shell Model. Parts I-VII,"
c                       Comp. Meth. Appl. Mech. & Engng, 1989-1991.
c-----[--.----+----.----+----.-----------------------------------------]
c       Control Information:
c       -------------------------
c       ndm = 3 ....... X, Y, Z Cartesian Coordinates at Nodes.
c       ndf = 5 ....... U-x, U-y, U-z, DT-1, DT-2 DOF at Nodes.
c       nen = 4 ....... # of Element Nodes: Input Counterclockwise.
c       nst = 20 ...... # of Element DOF's.
c-----[--.----+----.----+----.-----------------------------------------]
c       Element Material Set Input:
c       --------------
c       mate, Material-Group-Number
c         shell
c         Elastic Isotropic  E, Nu
c         Thickness shell    h, kappa
c         Density   mass     Rho, rot-mass-factor
c         displacement or mixed
c         finite

c       Input Variables:
c       ----------------
c       E ............. Young's modulus.
c       Nu ............ Poisson's ratio.
c       G-s ........... Transverse shear modulus - including
c                       shear correction factor.
c       thk ........... Shell thickness. If set to zero, shell
c                       thickness is assumed to be average of
c                       global thicknesses defined at nodes.
c                       If non-zero, this value overrides global
c                       input.
c       Rho-trn ....... Translational density.
c       Rho-rot ....... Rotational density.
c                       The reason for inputing two different density
c                       values is to be able to run fictitious
c                       material properties invented in
c                       Simo & Vu-Quoc beam formulation papers,
c                       e.g., pencil toss problem.
c       Opt-mix ....... Mixed treatment option:
c                       0 => Galerkin (displacement formulation),
c                       1 => Hellinger-Reissner (Pian-Summihara).
c-----[--.----+----.----+----.-----------------------------------------]

c       Quasi-Static Element Operation:
c       -------------------------------
c       DT,, Time-Step-Size
c       TOL,, Tolerance-Level           (Dflt: 1.d-10, Sugstd: 1.d-15)
c       PROP,, Proportional-Load-Option (Dflt: 1 - Continuous fctn)
c           Proportional-Load-Input     (Dflt: Linear with time)
c       LOOP,, Number-Of-Time-Steps
c           TIME
c           LOOP,, Number-Of-Newton-Iterations
c               TANG,,1
c           NEXT
c       NEXT

c       Transient Element Operation:
c       ----------------------------
c       DT,, Time-Step-Size
c       TOL,, Tolerance-Level           (Dflt: 1.d-10, Sugstd: 1.d-15)
c       PROP,, Proportional-Load-Option (Dflt: 1 - Continuous fctn)
c           Proportional-Load-Input     (Dflt: Linear with time)
c       BETA, CONS, Beta, Gamma, Alpha  (Dflt: 0.5, 1.0, 0.5)
c       CMAS
c       LOOP,, Number-Of-Time-Steps
c           TIME
c           LOOP,, Number-Of-Newton-Iterations
c               TANG,,1
c           NEXT
c       NEXT

c       ArcLength Operation:
c       --------------------
c       DT,, 1.0                        (Time step must be chosen 1)
c       TOL,, Tolerance-Level           (Dflt: 1.d-10, Sugstd: 1.d-15)
c       PROP,, Proportional-Load-Option (Dflt: 1 - Continuous fctn)
c           Proportional-Load-Input     (Dflt: Linear with time)
c       ARCL,, 2
c       LOOP,, Number-Of-Time-Steps
c           TIME
c           LOOP,, Number-Of-Newton-Iterations
c               TANG,,1
c           NEXT
c       NEXT

c       Note:
c       -----
c       ARCL,, 1 ...... Can be used to switch arclength direction
c                       during execution.

c-----[--.----+----.----+----.-----------------------------------------]
c       FEAP Common Arrays:
c       -------------------
c       hr............. History variable block.
c       mh ............ Global FEAP block common.

c       FEAP Solution Array:
c       --------------------
c       ul ............ Localized displacements and Velocities:
c                       ul(5,Node,1) := Total displacement T0-Tn+a.
c                       ul(5,Node,2) := Step displacement Tn-Tn+a.
c                       ul(5,Node,3) := Displacement increment. Tn+1
c                       ul(5,Node,4) := Velocity at Tn+a.
c                       ul(5,Node,5) := Acceleration at Tn+a.

c       Major Global Variables:
c       -----------------------
c       enflag ........ Energy conservation flag.
c       etotn,etot1 ... Total Energy at t-n and t-n+1.
c       epotn,epot1 ... Potential Energy at t-n and t-n+1.
c       alphar ........ Alpha residual.
c       alphat ........ Alpha tangent.
c       alphak ........ Alpha at previous iteration.
c       dalpha ........ Alpha Increment.
c       xln ........... Localized nodal orthogonal transformation
c                       matrices, xln(3,3,ElementNode,Ntime):
c                       Ntime=1 => Value at time t-n,
c                       Ntime=2 => Value at time t-n+alpha,
c                       Ntime=3 => Value at time t-n+1,
c                       Ntime=4 => Value at time 0.
c       rots .......... Localized nodal rotation vectors,
c                       rots(3,ElementNode,1) = Incremental rotation,
c                       rots(3,ElementNode,2) = Total rotation t-n+1.
c       rvel .......... Localized nodal rot. velocity vectors,
c                       rvel(3,ElementNode,Ntime):
c                       Ntime=1 => Value at time t-n,
c                       Ntime=2 => Value at time t-n+1.
c       racc .......... Localized nodal rot. acceleration vectors,
c                       racc(3,ElementNode,Ntime):
c                       Ntime=1 => Value at time t-n,
c                       Ntime=2 => Value at time t-n+1.
c       thkl .......... Localized nodal thickness.
c-----[--.----+----.----+----.-----------------------------------------]
c       Major Element Variables:
c       ------------------------
c       m ............. Transformation matrix current index:
c                       In stiffness assembly, m=2 => t-n+alpha,
c                       In time-step output,   m=3 => t-n+1.
c       dir ........... Element local directors
c                       dir(3,ElementNode,Ntime):
c                       Ntime=1 => Value at time t-n,
c                       Ntime=2 => Value at time t-n+alpha,
c                       Ntime=3 => Value at time t-n+1,
c                       Ntime=4 => Value at time 0.
c       xi ...........  Nodal transformation matrix at time m
c                       xi(3,3,ElementNode)
c       shp ........... Nodal shape function values.
c       shp1,shp2 ..... Nodal shape function natural derivatives.
c       shx1,shx2 ..... Nodal shape function cartesian derivatives.
c       sg ............ Gauss point coordinates and weights.
c       b ............. Discrete strain-displacement operator.
c       c ............. Resultant stress constitutive matrix.
c       zphG .......... Initial local coordinates
c                       at Gauss points.
c       zph1,zph2 ..... Initial coordinate local derivatives
c                       at Gauss points.
c       zpx1,zpx2 ..... Initial coordinate global derivatives
c                       at Gauss points.
c       zphm .......... Initial coordinate local derivatives
c                       at midside nodes.
c       cphi .......... Current local nodal coordinates.
c       cph1,cph2 ..... Current coordinate local derivatives
c                       at Gauss points.
c       cpx1,cpx2 ..... Current coordinate global derivatives
c                       at Gauss points.
c       cphm .......... Current coordinate local derivatives
c                       at midside nodes.
c       xjinv ......... Local-global coordinate Jacobian matrix.
c       xjs ........... Reference jacobian at Gauss points.
c       xjw ........... Reference jacobian weighted for integration.
c       zdir .......... Initial local directors
c                       at Gauss points.
c       zdx1,zdx2 ..... Initial director global derivatives
c                       at Gauss points.
c       zdrm .......... Initial local directors
c                       at midside nodes.
c       cdir .......... Current local directors
c                       at Gauss points.
c       cdx1,cdx2 ..... Current director global derivatives
c                       at Gauss points.
c       cdrm .......... Current local directors
c                       at midside nodes.
c       ze,ce ......... Initial, current membrane strain measure.
c       zx,cx ......... Initial, current shear strain measure.
c       zr,cr ......... Initial, current bending strain measure.
c       sn ............ Membrane stress.
c       sq ............ Shear stress.
c       sm ............ Bending stress.
c       aa ............ Gauss point accelerations at t-n+a.
c       cmem .......... Membrane constitutive relations coeffiecient
c       cbnd .......... Bending constitutive relations coeffiecient
c       nu ............ Poisson's ratio.
c       rhoa .......... Trans. thickness weighted density multiplyied
c                       by jacobian and weight at Gauss point.
c       rhoi .......... Rot. inertia weighted density multiplyied
c                       by jacobian and weight at Gauss point.
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      character*4     o,head
      common /bdata/  o,head(20)
      
      integer         numnp,numel,nummat,nen,neq,ipr
      common /cdata/  numnp,numel,nummat,nen,neq,ipr
      
      real*8          theta
      integer                  nrk,nrc,nrm,nrt,noi,nt,nrkn,nrcn,nrmn
      common /ddata/  theta(3),nrk,nrc,nrm,nrt,noi,nt,nrkn,nrcn,nrmn
      
      integer         ndof
      logical                 fsetr,frotas
      common /crotas/ ndof(9),fsetr,frotas
      
      real*8          dm
      integer            n,ma,mct,iel,nel
      common /eldata/ dm,n,ma,mct,iel,nel

      real*8          tt
      common /elplot/ tt(200)
      
      real*8          xln
      real*8          rots       ,rvel       ,racc       ,thkl
      integer                                                     rotyp
      common /erotas/ xln(9,9,4),
     &                rots(3,9,2),rvel(3,9,2),racc(3,9,2),thkl(9),rotyp
      
      logical         fl,    pfr
      common /fdata/  fl(12),pfr
      
      integer         nh1,nh2,nh3,nlm
      common /hdata/  nh1,nh2,nh3,nlm
      
      integer         ior,iow,ilg
      common /iofile/ ior,iow,ilg
      
      integer         stype,etype,dtype
      logical                           plasfl,viscfl,hflag,gflag
      common /pmod2d/ stype,etype,dtype,plasfl,viscfl,hflag,gflag
      
      integer         nph,ner
      real*8                  erav,j_int   ,jshft
      common /prstrs/ nph,ner,erav,j_int(3),jshft
      
      integer         istv, iste, istp
      common /strnum/ istv, iste, istp
      
      real*8          ttim,dt,c1,c2,c3,c4,c5, chi
      common /tdata/  ttim,dt,c1,c2,c3,c4,c5, chi

      real*8          hr
      integer               mr
      common          hr(1),mr(1)

      logical   optmix
      logical   rhs       , lhs        , errck

      integer   i         , j          , k          , l
      integer   i5        , i6         , j5         , j6
      integer   ij        , ik         , lint
      integer   ndf       , ndm        , nst        , isw

      real*8    rhoa      , rhoi       , xjdet      , xlgnrm
      real*8    xj1       , xj2        , xj3
      real*8    xj11      , xj12       , xj21       , xj22
      real*8    cmem      , cbnd
      real*8    dthk      , fac1       , facn

      real*8    f1   (3)  , f2  (3)    , y1  (4)    , y2  (4)
      real*8    d    (*)  , xl  (ndm,*), ul  (ndf,nen,*)
      real*8    s   (nst,*), p   (nst)
      real*8    sg   (3,4), shp (4,4)
      real*8    shp1 (4,4), shp2(4,4)  , shx1(4,4)  , shx2(4,4)
      real*8    fphm (3,4), fph1(3,4)  , fph2(3,4)  , fphi(3,4)
      real*8    bphm (3,4),                           bphi(3,4)
      real*8    zphm (3,4), zph1(3,4)  , zph2(3,4)
      real*8    cphm (3,4), cph1(3,4)  , cph2(3,4)  , cphi(3,4)
      real*8    zpx1 (3,4), zpx2(3,4)  , cpx1(3,4)  , cpx2(3,4)
      real*8    bpx1 (3,4), bpx2(3,4)  , fpx1(3,4)  , fpx2(3,4)
      real*8    zdrm (3,4), zdx1(3,4)  , zdx2(3,4)  , zdir(3,4)
      real*8    bdrm (3,4), bdx1(3,4)  , bdx2(3,4)  , bdir(3,4)
      real*8    cdrm (3,4), cdx1(3,4)  , cdx2(3,4)  , cdir(3,4)
      real*8    fdrm (3,4), fdx1(3,4)  , fdx2(3,4)  , fdir(3,4)
      real*8    ze   (3,4), zx  (2,4)  , zr  (3,4)
      real*8    be   (3,4), bx  (2,4)  , br  (3,4)
      real*8    ce   (3,4), cx  (2,4)  , cr  (3,4)
      real*8    fe   (3,4), fx  (2,4)  , fr  (3,4)
      real*8    sn   (3,4), sq  (2,4)  , sm  (3,4)
      real*8    aa   (3)  , xnt  (13)
      real*8    ploc (24) , sloc(24,24), stmp(3,3)
      real*8    c    (8,8), b(8,24,4)  , b1(8,24,4)
      real*8    xjinv(2,2), xjs (4)    , xjw (4)    , xlg (3,3)
      real*8    dir(3,9,4), xi(3,3,9)  ,xj(3,3,9)

      save

c     Set Thickness from Global Input

      if (isw.gt.1) then
        if(d(14).lt.1.d-5) then
          dthk  = 0.25d0 * (thkl(1)+thkl(2)+thkl(3)+thkl(4))
          if(dthk.lt.0.0d0) then
            write(iow,3002) dthk
C            call plstop()
          endif
        else
          dthk  = d(14)
        endif
      endif

c     Unused FEAP Task Options (isw=2,7,9,10,11,12,14):

c     Input Material Properties:

      if(isw.eq.1) then

c       Check input data

        errck = .false.
        if(d(1).le.0.0d0) then
          write(iow,3000)
          errck = .true.
        endif
        if(d(14).le.0.0d0) then
          write(iow,3001)
        endif
        if(errck) stop '  Shell material property errors'

c       Initialize Shape Functions.

        l     = 2
        call int2d      ( l       , lint    , sg )
        do i = 1 , lint
          call sh3fshap ( sg (1,i), shp(1,i), shp1(1,i), shp2(1,i)  )
        end do ! i

c       Set default rotational update type

        rotyp = -1

c     Compute Element Dynamic Tangent and Residual Array:

      elseif(isw.ge.3 .and. isw.le.6  .or. isw.eq.8 .or. isw.eq.13) then

c       Set local element director fields

        call sh3finte (xl, xln, ndof, ndm, nel,  dir, xi, xj)

c       Compute Midsurface Gauss-Point Interpolation

        do i = 1 , lint
          do j = 1 , 3

c           Midsurface

            zph1(j,i) = shp1(1,i)*xl(j,1) + shp1(2,i)*xl(j,2)
     &                + shp1(3,i)*xl(j,3) + shp1(4,i)*xl(j,4)
            zph2(j,i) = shp2(1,i)*xl(j,1) + shp2(2,i)*xl(j,2)
     &                + shp2(3,i)*xl(j,3) + shp2(4,i)*xl(j,4)
          end do ! j

c         Compute Reference Surface Jacobian

          xj1    =  zph1(2,i)*zph2(3,i) - zph1(3,i)*zph2(2,i)
          xj2    = -zph1(1,i)*zph2(3,i) + zph1(3,i)*zph2(1,i)
          xj3    =  zph1(1,i)*zph2(2,i) - zph1(2,i)*zph2(1,i)
          xjs(i) =  sqrt( xj1**2 + xj2**2 + xj3**2 )
          if(isw.eq.4) then
            xjw(i) =  1.d0
          else
            xjw(i) =  sg(3,i) * xjs(i)
          endif

c       End Gauss Loop

        end do ! i

c       Do mass computation

        if(isw.eq.5) then
          call sh3fmasg( d, dthk, shp, xjw, s, lint, nst)
          return
        endif

c       Set Flags from Input

        optmix = d(17).gt.1.d0
        lhs    = isw.eq.3
        rhs    = .true.

c       Compute Current Configuration

        fac1 = 1.d0 / theta(3) - 1.d0
        do i = 1 , 4
          do j = 1 , 3
            cphi(j,i) = xl(j,i)   + ul(j,i,1)
            bphi(j,i) = cphi(j,i) - ul(j,i,2)
            fphi(j,i) = cphi(j,i) + ul(j,i,2)*fac1
          end do ! j
        end do ! i

c       Compute Midside Interpolations

        call sh3fintr ( xl   , dir(1,1,4) , zphm , zdrm )
        call sh3fintr ( bphi , dir(1,1,1) , bphm , bdrm )
        call sh3fintr ( cphi , dir(1,1,2) , cphm , cdrm )
        call sh3fintr ( fphi , dir(1,1,3) , fphm , fdrm )

c       Compute Gauss-Point Interpolations and Derivatives

        fac1 = theta(3)
        facn = 1.d0 - fac1

        do i = 1 , lint
          do j = 1 , 3

c           Midsurface

            cph1(j,i) = shp1(1,i)*cphi(j,1)  + shp1(2,i)*cphi(j,2)
     &                + shp1(3,i)*cphi(j,3)  + shp1(4,i)*cphi(j,4)
            cph2(j,i) = shp2(1,i)*cphi(j,1)  + shp2(2,i)*cphi(j,2)
     &                + shp2(3,i)*cphi(j,3)  + shp2(4,i)*cphi(j,4)
            fph1(j,i) = shp1(1,i)*fphi(j,1)  + shp1(2,i)*fphi(j,2)
     &                + shp1(3,i)*fphi(j,3)  + shp1(4,i)*fphi(j,4)
            fph2(j,i) = shp2(1,i)*fphi(j,1)  + shp2(2,i)*fphi(j,2)
     &                + shp2(3,i)*fphi(j,3)  + shp2(4,i)*fphi(j,4)

c           Directors

            zdir(j,i) = shp (1,i)*dir(j,1,4) + shp (2,i)*dir(j,2,4)
     &                + shp (3,i)*dir(j,3,4) + shp (4,i)*dir(j,4,4)
            bdir(j,i) = shp (1,i)*dir(j,1,1) + shp (2,i)*dir(j,2,1)
     &                + shp (3,i)*dir(j,3,1) + shp (4,i)*dir(j,4,1)
            cdir(j,i) = shp (1,i)*dir(j,1,2) + shp (2,i)*dir(j,2,2)
     &                + shp (3,i)*dir(j,3,2) + shp (4,i)*dir(j,4,2)
            fdir(j,i) = shp (1,i)*dir(j,1,3) + shp (2,i)*dir(j,2,3)
     &                + shp (3,i)*dir(j,3,3) + shp (4,i)*dir(j,4,3)
          end do ! j

c         Compute Surface Normal

          call vecp ( zph1(1,i) , zph2(1,i) , xlg(1,3) )

          xlgnrm   = 1.d0/sqrt(xlg(1,3)**2 + xlg(2,3)**2 + xlg(3,3)**2)
          xlg(1,3) = xlg(1,3) * xlgnrm
          xlg(2,3) = xlg(2,3) * xlgnrm
          xlg(3,3) = xlg(3,3) * xlgnrm

c         Compute Local-Global Jacobian

          call sh3flmda ( xlg(1,3) , xlg  )

          xj11       = xlg(1,1)*zph1(1,i) + xlg(2,1)*zph1(2,i)
     &               + xlg(3,1)*zph1(3,i)
          xj21       = xlg(1,2)*zph1(1,i) + xlg(2,2)*zph1(2,i)
     &               + xlg(3,2)*zph1(3,i)
          xj12       = xlg(1,1)*zph2(1,i) + xlg(2,1)*zph2(2,i)
     &               + xlg(3,1)*zph2(3,i)
          xj22       = xlg(1,2)*zph2(1,i) + xlg(2,2)*zph2(2,i)
     &               + xlg(3,2)*zph2(3,i)

          xjdet      = 1.d0/(xj11*xj22 - xj12*xj21)

          xjinv(1,1) = xj22 * xjdet
          xjinv(1,2) =-xj12 * xjdet
          xjinv(2,1) =-xj21 * xjdet
          xjinv(2,2) = xj11 * xjdet

c         Compute Global Shape Function Derivatives

          call sh3fmdsh ( shp1(1,i) , shp2(1,i) , xjinv     ,
     &                    shx1(1,i) , shx2(1,i) )

c         Interpolate Global Derivatives

          do j = 1 , 3

c           Midsurface

            zpx1(j,i) = shx1(1,i)*xl  (j,1)    + shx1(2,i)*xl  (j,2)
     &                + shx1(3,i)*xl  (j,3)    + shx1(4,i)*xl  (j,4)
            zpx2(j,i) = shx2(1,i)*xl  (j,1)    + shx2(2,i)*xl  (j,2)
     &                + shx2(3,i)*xl  (j,3)    + shx2(4,i)*xl  (j,4)
            cpx1(j,i) = shx1(1,i)*cphi(j,1)    + shx1(2,i)*cphi(j,2)
     &                + shx1(3,i)*cphi(j,3)    + shx1(4,i)*cphi(j,4)
            cpx2(j,i) = shx2(1,i)*cphi(j,1)    + shx2(2,i)*cphi(j,2)
     &                + shx2(3,i)*cphi(j,3)    + shx2(4,i)*cphi(j,4)
            fpx1(j,i) = shx1(1,i)*fphi(j,1)    + shx1(2,i)*fphi(j,2)
     &                + shx1(3,i)*fphi(j,3)    + shx1(4,i)*fphi(j,4)
            fpx2(j,i) = shx2(1,i)*fphi(j,1)    + shx2(2,i)*fphi(j,2)
     &                + shx2(3,i)*fphi(j,3)    + shx2(4,i)*fphi(j,4)
            bpx1(j,i) = shx1(1,i)*bphi(j,1)    + shx1(2,i)*bphi(j,2)
     &                + shx1(3,i)*bphi(j,3)    + shx1(4,i)*bphi(j,4)
            bpx2(j,i) = shx2(1,i)*bphi(j,1)    + shx2(2,i)*bphi(j,2)
     &                + shx2(3,i)*bphi(j,3)    + shx2(4,i)*bphi(j,4)

c           Directors

            zdx1(j,i) = shx1(1,i)*dir(j,1,4) + shx1(2,i)*dir(j,2,4)
     &                + shx1(3,i)*dir(j,3,4) + shx1(4,i)*dir(j,4,4)
            zdx2(j,i) = shx2(1,i)*dir(j,1,4) + shx2(2,i)*dir(j,2,4)
     &                + shx2(3,i)*dir(j,3,4) + shx2(4,i)*dir(j,4,4)
            bdx1(j,i) = shx1(1,i)*dir(j,1,1) + shx1(2,i)*dir(j,2,1)
     &                + shx1(3,i)*dir(j,3,1) + shx1(4,i)*dir(j,4,1)
            bdx2(j,i) = shx2(1,i)*dir(j,1,1) + shx2(2,i)*dir(j,2,1)
     &                + shx2(3,i)*dir(j,3,1) + shx2(4,i)*dir(j,4,1)
            cdx1(j,i) = shx1(1,i)*dir(j,1,2) + shx1(2,i)*dir(j,2,2)
     &                + shx1(3,i)*dir(j,3,2) + shx1(4,i)*dir(j,4,2)
            cdx2(j,i) = shx2(1,i)*dir(j,1,2) + shx2(2,i)*dir(j,2,2)
     &                + shx2(3,i)*dir(j,3,2) + shx2(4,i)*dir(j,4,2)
            fdx1(j,i) = shx1(1,i)*dir(j,1,3) + shx1(2,i)*dir(j,2,3)
     &                + shx1(3,i)*dir(j,3,3) + shx1(4,i)*dir(j,4,3)
            fdx2(j,i) = shx2(1,i)*dir(j,1,3) + shx2(2,i)*dir(j,2,3)
     &                + shx2(3,i)*dir(j,3,3) + shx2(4,i)*dir(j,4,3)
          end do ! j

c         Compute Strain Measures

          call sh3fstrn ( shp1(1,i)     , shp2(1,i)   ,
     &                    zphm          ,
     &                    zpx1(1,i)     , zpx2(1,i)   ,
     &                    zdrm          ,
     &                    zdx1(1,i)     , zdx2(1,i)   ,
     &                    ze  (1,i)     , zx  (1,i)   ,
     &                    zr  (1,i)                   )
          call sh3fstrn ( shp1(1,i)     , shp2(1,i)   ,
     &                    fphm          ,
     &                    fpx1(1,i)     , fpx2(1,i)   ,
     &                    fdrm          ,
     &                    fdx1(1,i)     , fdx2(1,i)   ,
     &                    fe  (1,i)     , fx  (1,i)   ,
     &                    fr  (1,i)                   )

          if (isw.eq.13 .or. .not.fl(9)) then
            do j = 1,3
              ce(j,i) = fe(j,i) - ze(j,i)
              cr(j,i) = fr(j,i) - zr(j,i)
            end do ! j
            do j = 1,2
              cx(j,i) = fx(j,i) - zx(j,i)
            end do ! j
          else
            call sh3fstrn ( shp1(1,i)     , shp2(1,i)   ,
     &                      bphm          ,
     &                      bpx1(1,i)     , bpx2(1,i)   ,
     &                      bdrm          ,
     &                      bdx1(1,i)     , bdx2(1,i)   ,
     &                      be  (1,i)     , bx  (1,i)   ,
     &                      br  (1,i)                   )
            do j = 1 , 3
              ce(j,i) = fac1*fe(j,i) + facn*be(j,i) - ze(j,i)
              cr(j,i) = fac1*fr(j,i) + facn*br(j,i) - zr(j,i)
            end do ! j
            do j = 1 , 2
              cx(j,i) = fac1*fx(j,i) + facn*bx(j,i) - zx(j,i)
            end do ! j
          endif

c         Set Up Strain-Displacement Matrices

          call sh3fbmtx ( shp1(1,i)     , shp2(1,i)   ,
     &                    shx1(1,i)     , shx2(1,i)   ,
     &                    cphm          ,
     &                    cpx1(1,i)     , cpx2(1,i)   ,
     &                    cdrm          ,
     &                    cdx1(1,i)     , cdx2(1,i)   ,
     &                    b   (1,1,i)                 )

          if (lhs) then

            call sh3fbmtx ( shp1(1,i)     , shp2(1,i)   ,
     &                      shx1(1,i)     , shx2(1,i)   ,
     &                      fphm          ,
     &                      fpx1(1,i)     , fpx2(1,i)   ,
     &                      fdrm          ,
     &                      fdx1(1,i)     , fdx2(1,i)   ,
     &                      b1   (1,1,i)                 )
          endif

c       End Gauss Loop

        end do ! i

c       Zero full (local) stiffness and residuals

        do i = 1,24
          do j = 1,24
            sloc(j,i) = 0.0d0
          end do ! j
          ploc(i) = 0.0d0
        end do ! i

c       General Pian-Sumihara Calculations

        if (optmix) then
          cmem = d(1)*dthk
          cbnd = d(1)*dthk**3/12.d0
          call sh3fpian ( xl   , xjw , sg , cmem , cbnd ,
     &                    d(2) , f1  , f2 , y1 , y2   , xnt  )

c         Pian-Sumihara Membrane Treatment

          call sh3fmbrn ( xjw  , sg   , cmem , d(2) ,
     &                    f1   , f2   , y1  , y2   , xnt  ,
     &                    6    , 24   , ce  , b    , sn   ,
     &                    ploc , sloc )

c         Pian-Sumihara Bending Treatment

          call sh3fbend ( xjw  , sg  , cbnd , d(2) ,
     &                    f1   , f2  , y1 , y2   , xnt  ,
     &                    6    , 24  , cr , b    , sm   ,
     &                    ploc , sloc )
        endif

c       Galerkin and Shear Treatment

        do i = 1,lint

c         Compute Constitutive Relations (weighted by jacobian)

          call sh3fcmtx ( d      , dthk      , xjw(i)    ,
     &                    optmix , zph1(1,i) , zph2(1,i) , c )

c         Compute Stresses (weighted by jacobian)

          call sh3fstrs ( c         ,
     &                    ce(1,i)   , cx(1,i)   , cr(1,i)   ,
     &                    optmix    ,
     &                    sn(1,i)   , sq(1,i)   , sm(1,i)   )

c         Store time history plot data for element

          i6 = 8*(l-1)
          do j = 1,3
            tt(j+i6  ) = sn(j,i)
            tt(j+i6+3) = sm(j,i)
          end do ! j
          tt(7+i6) = sq(1,i)
          tt(8+i6) = sq(2,i)

c         Compute Residual

          call sh3fresd ( b (1,1,i) ,
     &                    sn(1,i)   , sq(1,i)   , sm(1,i)   ,
     &                    optmix    , ploc      )

          if (lhs) then

c           Compute Material Stiffness

            call sh3fstif ( b (1,1,i) , b1(1,1,i) , c         ,
     &                      optmix    , sloc                  )

c           Compute Geometric Stiffness

            if(gflag) then
              call sh3fgeom ( shp1(1,i) , shp2(1,i) ,
     &                        shx1(1,i) , shx2(1,i) ,
     &                        cphm      , cpx1(1,i) , cpx2(1,i) ,
     &                        dir(1,1,2),
     &                        sn  (1,i) , sq  (1,i) , sm  (1,i) ,
     &                        sloc                              )
            endif
          endif

c       End Gauss Loop

        end do ! i

c       Compute Transient Terms

        if (fl(9) .and. noi.ne.6) then

          do i = 1,lint

c           Interpolate Accelerations

            do j = 1 , 3

              aa(j) = shp(1,i)*ul(j,1,5) + shp(2,i)*ul(j,2,5)
     &              + shp(3,i)*ul(j,3,5) + shp(4,i)*ul(j,4,5)

            end do ! j

c           Compute Weighted Density/Inertia

            rhoa = d(4) * d(14) * xjw(i)
            rhoi = d(8) * rhoa * d(14)**2/12.d0

c           Compute Transient Residual and Tangent

            call sh3ftran ( dt,  shp(1,i)         , aa               ,
     &                      rots(1,1,1)           , rots(1,1,2)      ,
     &                      rvel(1,1,1)           , rvel(1,1,2)      ,
     &                      racc(1,1,1)           , racc(1,1,2)      ,
     &                      dir(1,1,2)            , dir(1,1,4)       ,
     &                      rhoa         , rhoi   , ploc      , sloc )

c         End Gauss Loop

          end do ! i

        endif

c       Stress and energy processes

        if (isw.eq.4 .or. isw.eq.8)  then
          call sh3fstre ( xl,shp, ce,cr,cx,sn,sq,sm,xjw, lint,ndm,isw )
          rhs = .false.
        elseif (isw.eq.13) then
          call sh3fener ( d,dthk, ce,cr,cx,sn,sq,sm,
     &                    ul,shp,fphi,xjw, ndf,lint )
          rhs = .false.
        elseif(lhs) then

c         Scale tangent stiffness for total increment du_t-n+1

          do j = 1,24
            do i = 1,24
              sloc(i,j) = fac1*sloc(i,j)
            end do ! i
          end do ! j

c         Transform (Reduce) Stiffness/Residual - First Node Loop

          i5 = 0
          i6 = 0

c         First Node Loop

          do i = 1 , 4

            j5 = 0
            j6 = 0

c           Second Node Loop

            do j = 1 , 4

c             Stiffness Translational Terms

              do l = 1 , 3
                do k = 1 , 3
                  s(i5+k,j5+l) = sloc(i6+k,j6+l)
                end do ! k
              end do ! l

c             Stiffness Off-diagonal Terms

              do l = 1 , ndof(j) - 3
                do k = 1 , 3
                  s(i5+k,j5+l+3) = sloc(i6+k,j6+4)*xj(1,l,j)
     &                           + sloc(i6+k,j6+5)*xj(2,l,j)
     &                           + sloc(i6+k,j6+6)*xj(3,l,j)
                end do ! k
              end do ! l
              do l = 1 , ndof(i) - 3
                do k = 1 , 3
                  s(i5+l+3,j5+k) = sloc(i6+4,j6+k)*xi(1,l,i)
     &                           + sloc(i6+5,j6+k)*xi(2,l,i)
     &                           + sloc(i6+6,j6+k)*xi(3,l,i)
                end do ! k
              end do ! l

c             Stiffness Rotational Terms

              do l = 1 , ndof(j) - 3
                do k = 1 , 3
                  stmp(k,l) = sloc(i6+k+3,j6+4)*xj(1,l,j)
     &                      + sloc(i6+k+3,j6+5)*xj(2,l,j)
     &                      + sloc(i6+k+3,j6+6)*xj(3,l,j)
                end do ! k
              end do ! l
              do l = 1 , ndof(j) - 3
                do k = 1 , ndof(i) - 3
                  s(i5+k+3,j5+l+3) = stmp(1,l)*xi(1,k,i)
     &                             + stmp(2,l)*xi(2,k,i)
     &                             + stmp(3,l)*xi(3,k,i)
                end do ! k
              end do ! l

              j5 = j5 + ndf
              j6 = j6 + 6

c           End Second Node Loop

            end do ! j

            i5 = i5 + ndf
            i6 = i6 + 6

c         End First Node Loop

          end do ! i

        endif

c       Node Loop

        if(rhs) then

          i5 = 0
          i6 = 0

          do i = 1 , 4

c           Residual Terms

            do j = 1 , 3
              p(i5+j) = ploc(i6+j)
              ij      = i6 + 3 + j
              do k = 1 , ndof(i) - 3
                ik = i5 + 3 + k
                p(ik) = p(ik) + ploc(ij)*xi(j,k,i)
              end do ! k
            end do ! j

            i5 = i5 + ndf
            i6 = i6 + 6

c         End Node Loop

          end do ! i

c         Add geometric stiffness term for 6 dof nodes

          if(lhs .and. ndf.eq.6 .and. gflag)then
            call sh3fgeo6 (lint , shp1 , shp2 , shx1 , shx2 ,
     &                    cphm , cpx1 , cpx2 ,  dir (1,1,2) ,
     &                    sq  , sm    , s                   )
          end if

c         Modify residual and stiffnes for pressure &body loading

          call sh3fsurf(d,xl,ul,xjw, ndf,ndm,nst, p,s)

        endif
      endif

c     Formats:

3000  format(' *ERROR* Modulus is zero or negative for shell')

3001  format(' *WARNING* Thickness zero: Compute from nodal input.')

3002  format(' *ERROR* Thickness is zero or negative for shell')

      end

      
      subroutine int2d(l,lint,sg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form Gauss points and weights for two dimensions

c      Inputs:
c         l       - Number of points/direction

c      Outputs:
c         lint    - Total number of points
c         sg(3,*) - Array of points and weights
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8          dm
      integer            n,ma,mct,iel,nel
      common /eldata/ dm,n,ma,mct,iel,nel
      
      integer         ior,iow,ilg
      common /iofile/ ior,iow,ilg

      integer   i,j,k,l,lint, lr(9),lz(9),lw(9)
      real*8    g,h, third, sg(3,*),ss(5),ww(5)

      save

      data      lr/-1,1,1,-1,0,1,0,-1,0/,lz/-1,-1,1,1,-1,0,1,0,0/
      data      lw/4*25,4*40,64/
      data      third / 0.3333333333333333d0 /

c     Set number of total points

      lint = l*l

c     5 pt. integration

      if(l.eq.0) then

        lint = 5
        g    = sqrt(0.6d0)
        do i = 1,4
          sg(1,i) = g*lr(i)
          sg(2,i) = g*lz(i)
          sg(3,i) = 5.d0/9.d0
        end do ! i

        sg(1,5) = 0.0d0
        sg(2,5) = 0.0d0
        sg(3,5) = 16.d0/9.d0

c     1x1 integration

      elseif(l.eq.1) then
        sg(1,1) = 0.d0
        sg(2,1) = 0.d0
        if(nel.eq.3) sg(2,1) = -third
        sg(3,1) = 4.d0

c     2x2 integration

      elseif(l.eq.2) then
        g = sqrt(third)
        do i = 1,4
          sg(1,i) = g*lr(i)
          sg(2,i) = g*lz(i)
          sg(3,i) = 1.d0
        end do ! i

c     3x3 integration

      elseif(l.eq.3) then
        g = sqrt(0.6d0)
        h = 1.d0/81.d0
        do i = 1,9
          sg(1,i) = g*lr(i)
          sg(2,i) = g*lz(i)
          sg(3,i) = h*lw(i)
        end do ! i

c     4x4 integration

      elseif(l.eq.4) then
        g     = sqrt(4.8d0)
        h     = third/g
        ss(1) = sqrt((3.d0+g)/7.d0)
        ss(4) = - ss(1)
        ss(2) = sqrt((3.d0-g)/7.d0)
        ss(3) = -ss(2)
        ww(1) = 0.5d0 - h
        ww(2) = 0.5d0 + h
        ww(3) = 0.5d0 + h
        ww(4) = 0.5d0 - h
        i = 0
        do j = 1,4
          do k = 1,4
            i = i + 1
            sg(1,i) = ss(k)
            sg(2,i) = ss(j)
            sg(3,i) = ww(j)*ww(k)
          end do ! k
        end do ! i

c     5x5 integration

      elseif(l.eq.5) then

        g     =  sqrt(1120.d0)
        ss(1) =  sqrt((70.d0 + g)/126.d0)
        ss(2) =  sqrt((70.d0 - g)/126.d0)
        ss(3) =  0.0d0
        ss(4) = -ss(2)
        ss(5) = -ss(1)

        ww(1) =  (21.d0*g + 117.6d0)/(g*(70.d0 + g))
        ww(2) =  (21.d0*g - 117.6d0)/(g*(70.d0 - g))
        ww(3) =  2.d0*(1.d0 - ww(1) - ww(2))
        ww(4) =  ww(2)
        ww(5) =  ww(1)

        i = 0
        do j = 1,5
          do k = 1,5
            i = i + 1
            sg(1,i) = ss(k)
            sg(2,i) = ss(j)
            sg(3,i) = ww(j)*ww(k)
          end do ! k
        end do ! j

c     Error

      else

        write(ilg,2000) l
        write(iow,2000) l
        if(ior.lt.0) then
          write(*,2000) l
        endif
C        call plstop()

      endif

c     Format

2000  format(' *ERROR* INT2D: Illegal quadrature order =',i16)

      end

      
      subroutine sh3fpian ( xl  , xjw  , sg  , eh  , eh12 ,
     &                        xnu , f1   , f2  , y1  , y2  , xnt  )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    SH3FPIAN is a subroutine which computes some
c                      objects needed by Pian-Summihara-type
c                      Hellinger-Reissner mixed element.

c      Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c      Date:           February 1991.

c      Version:         This routine was tested in FEAP version 6.3
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      xl ............ Initial nodal coordinates
c      sg ............ Gauss point coordinates.
c      xjw ........... Reference jacobian weighted for integration.
c      eh ............ Membrane inverse constitutive coefficient.
c      eh12 .......... bending inverse constitutive coefficient.
c      xnu ........... Poisson's ratio.

c      Routine Output:
c      ---------------
c      f1,f2 ......... Tensor transformations in vector format
c                      at center of element.
c      y1,y2, xnt .... Commonly used parts of [H].
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      integer   i
      real*8    xlcnrm , eta   , xi     , etaxi  , xjweh   , xjwE2
      real*8    eh     , eh12  , xnu    , xnu12
      real*8    xl(3,4), xjw(4), zpc1(3), zpc2(3), xjc(2,2), xlc(3,3)
      real*8    sg(3,4), f1(3) , f2(3)  , y1(4)  , y2(4)   , xnt(13)
      real*8    dot

c     Compute central Local Derivatives of Coordinates

      do i = 1,3
        zpc1(i) = 0.25d0 * (-xl(i,1)+xl(i,2)+xl(i,3)-xl(i,4))
        zpc2(i) = 0.25d0 * (-xl(i,1)-xl(i,2)+xl(i,3)+xl(i,4))
      end do ! i

c     Compute central Surface Normal

      call vecp ( zpc1(1) , zpc2(1) , xlc(1,3) )

      xlcnrm   = 1.d0/sqrt ( dot(xlc(1,3),xlc(1,3),3) )
      xlc(1,3) = xlc(1,3) * xlcnrm
      xlc(2,3) = xlc(2,3) * xlcnrm
      xlc(3,3) = xlc(3,3) * xlcnrm

c     Compute Local-Global Jacobian

      call sh3flmda ( xlc(1,3) , xlc  )

      xjc(1,1) = dot ( xlc(1,1) , zpc1(1) ,3 )
      xjc(2,1) = dot ( xlc(1,2) , zpc1(1) ,3 )
      xjc(1,2) = dot ( xlc(1,1) , zpc2(1) ,3 )
      xjc(2,2) = dot ( xlc(1,2) , zpc2(1) ,3 )

c     Compute central Transformation Tensor in Vector Form

      f1(1) = xjc(1,1) * xjc(1,1)
      f1(2) = xjc(2,1) * xjc(2,1)
      f1(3) = xjc(2,1) * xjc(1,1)
      f2(1) = xjc(1,2) * xjc(1,2)
      f2(2) = xjc(2,2) * xjc(2,2)
      f2(3) = xjc(2,2) * xjc(1,2)

c     Compute Element Integrals of 1,xi,eta and Multiples thereof

      do i = 1 , 13
        xnt(i) = 0.0d0
      end do ! i
      do i = 1 , 4
        xnt(1) = xnt(1) + xjw(i)
        xnt(2) = xnt(2) + xjw(i) * sg(2,i)
        xnt(3) = xnt(3) + xjw(i) * sg(1,i)
        xnt(4) = xnt(4) + xjw(i) * sg(2,i) * sg(1,i)
      end do ! i

      xnt(2) = xnt(2) / xnt(1)
      xnt(3) = xnt(3) / xnt(1)
      xnt(4) = xnt(4) / xnt(1)

      do i = 1 , 4
        eta     = sg(2,i)       - xnt(2)
        xi      = sg(1,i)       - xnt(3)
        etaxi   = sg(2,i)*sg(1,i) - xnt(4)
        xjweh   = xjw(i) / eh
        xjwE2   = xjw(i) / eh12
        xnt(5 ) = xnt(5 ) + xjweh * eta ** 2
        xnt(6 ) = xnt(6 ) + xjweh * xi  ** 2
        xnt(7 ) = xnt(7 ) + xjweh * eta *  xi
        xnt(8 ) = xnt(8 ) + xjweh * xi  *  etaxi
        xnt(9 ) = xnt(9 ) + xjweh * eta *  etaxi
        xnt(10) = xnt(10) + xjweh * etaxi ** 2
        xnt(11) = xnt(11) + xjwE2 * eta ** 2
        xnt(12) = xnt(12) + xjwE2 * xi  ** 2
        xnt(13) = xnt(13) + xjwE2 * eta *  xi
      end do ! i

c     Compute Common Parts of [H]

      xnu12 = 2.d0  * (1.d0  + xnu)
      y1(1) = f1(1) - xnu * f1(2)
      y1(2) = f1(2) - xnu * f1(1)
      y1(3) = f1(3) * xnu12
      y2(1) = f2(1) - xnu * f2(2)
      y2(2) = f2(2) - xnu * f2(1)
      y2(3) = f2(3) * xnu12
      y1(4) = (f1(1)+f1(2)) * xnu
      y2(4) = (f2(1)+f2(2)) * xnu

      end

      subroutine sh3fmbrn ( xjw , sg  , eh  , xnu ,
     &                      f1  , f2  , y1  , y2  , xnt ,
     &                      ndf , nst , ce  , b   , sn  ,
     &                      p   , s   )
c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    SH3FMBRN is a subroutine which computes the
c                      stiffness and residual contributions of the
c                      membrane field, using Pian-Sumihara
c                      stress interpolations for Shell.

c      Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c      Date:           February 1991.

c      Version:        This routine was tested in FEAP version 6.3.
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      sg ............ Gauss point coordinates.
c      xjw ........... Reference jacobian weighted for integration.
c      eh ............ Membrane inverse constitutive coefficient.
c      xnu ........... Poisson's ratio.
c      f1,f2 ......... Tensor transformations in vector format
c                      at center of element.
c      y1,y2, xnt .... Commonly used parts of [H].
c      ce ............ Current membrane strain measure.
c      b ............. Discrete strain-displacement operator.
c      ndf ........... Number of degrees of freedom per node.
c      nst ........... Number of degrees of freedom per element.

c      Routine Output:
c      ---------------
c      sn ............ Membrane stress.
c      s ............. Element stiffness.
c      p ............. Element residual.
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      integer   i      , j      , k        , l       , ndf
      integer   k1     , k2     ,  nst     , n
      integer   ni     , nj
      real*8    eta    , xi     , fac0     , fac1    , fac2
      real*8    fach1  , fach2  , eh       , xnu     , hdet
      real*8    xjw(4) , sg(3,4), f1(3)    , f2(3)
      real*8    xnt(13), ce(3,4), sn(3,4)  , htm(2,2), hin(5,5)
      real*8    ee(5)  , be(5)  , b(8,24,4), g(4,5,3), gh(3,5)
      real*8    y1(4)  , y2(4)  , p(nst)   , s(nst,nst)
      real*8    dot

      save

c     Zero Matrices

      do i = 1,5
        ee(i) = 0.0d0
        do j = 1,5
          hin(j,i) = 0.0d0
        end do ! j
        do j = 1,4
          g(j,i,1) = 0.0d0
          g(j,i,2) = 0.0d0
          g(j,i,3) = 0.0d0
        end do ! j
      end do ! i

c     Gauss Point Loop

      do i = 1 , 4

c       Common Factors

        eta   = sg(2,i) - xnt(2)
        xi    = sg(1,i) - xnt(3)
        fac0  = xjw(i)
        fac1  = fac0 * eta
        fac2  = fac0 * xi

c       Modified Strains [ee] = [S][ce]

        ee(1) = ee(1) + fac0 * ce(1,i)
        ee(2) = ee(2) + fac0 * ce(2,i)
        ee(3) = ee(3) + fac0 * ce(3,i)
        ee(4) = ee(4) + fac1 * dot ( f1 , ce(1,i) ,3 )
        ee(5) = ee(5) + fac2 * dot ( f2 , ce(1,i) ,3 )

c       Integrated Strain-Displacement Matrix [G] = [S][B]

        do k = 1 , 4
          do j = 1 , 3
            l = (k-1)*ndf + j
            g(k,1,j) = g(k,1,j) + fac0 * b(1,l,i)
            g(k,2,j) = g(k,2,j) + fac0 * b(2,l,i)
            g(k,3,j) = g(k,3,j) + fac0 * b(3,l,i)
            g(k,4,j) = g(k,4,j) + fac1 * dot(f1,b(1,l,i),3 )
            g(k,5,j) = g(k,5,j) + fac2 * dot(f2,b(1,l,i),3 )
          end do ! j
        end do ! k

      end do ! i

c     Compute Second Block of [H-inv]

      htm(1,1) = xnt(6) * dot ( f1 , y1 ,3 )
      htm(2,2) = xnt(5) * dot ( f2 , y2 ,3 )
      htm(1,2) = xnt(7) * dot ( f1 , y2 ,3 )
      hdet     = 1.d0 / (htm(1,1)*htm(2,2) - htm(1,2)**2)

c     Constitutive Relations

c     Compute [H-inv]

      fach1    = eh / (xnt(1) * (1.d0 - xnu * xnu))
      fach2    = eh / (xnt(1) * (1.d0 + xnu) * 2.d0)
      hin(1,1) = fach1
      hin(1,2) = fach1 * xnu
      hin(2,1) = hin(1,2)
      hin(2,2) = hin(1,1)
      hin(3,3) = fach2
      hin(4,4) = htm(2,2) * hdet
      hin(5,5) = htm(1,1) * hdet
      hin(4,5) =-htm(1,2) * hdet
      hin(5,4) = hin(4,5)

c     Compute Stress Variables

      do i = 1,5
        be(i) = hin(i,1)*ee(1)
     1        + hin(i,2)*ee(2)
     2        + hin(i,3)*ee(3)
     3        + hin(i,4)*ee(4)
     4        + hin(i,5)*ee(5)
      end do ! i

c     Gauss Point Loop

      do i = 1 , 4

c       Compute Stresses

        eta     = sg(2,i) - xnt(2)
        xi      = sg(1,i) - xnt(3)
        sn(1,i) = be(1) + eta * f1(1)*be(4) + xi * f2(1)*be(5)
        sn(2,i) = be(2) + eta * f1(2)*be(4) + xi * f2(2)*be(5)
        sn(3,i) = be(3) + eta * f1(3)*be(4) + xi * f2(3)*be(5)

c       Weight Stresses by Jacobian

        sn(1,i) = sn(1,i) * xjw(i)
        sn(2,i) = sn(2,i) * xjw(i)
        sn(3,i) = sn(3,i) * xjw(i)

      end do ! i

c     Node Loop #1

      do i = 1 , 4
        ni = (i-1)*ndf

c       Compute Residual

        do k = 1 , 3
          p(ni+k) = p(ni+k) - g(i,1,k)*be(1)
     &                      - g(i,2,k)*be(2)
     &                      - g(i,3,k)*be(3)
     &                      - g(i,4,k)*be(4)
     &                      - g(i,5,k)*be(5)
        end do ! i

c       Node Loop #2

        do j = 1 , 4
          nj = (j-1)*ndf

c         Compute Stiffness (m-m) Part

          do n = 1,5
            do k1 = 1,3
              gh(k1,n) = g(i,1,k1)*hin(1,n)
     &                 + g(i,2,k1)*hin(2,n)
     &                 + g(i,3,k1)*hin(3,n)
     &                 + g(i,4,k1)*hin(4,n)
     &                 + g(i,5,k1)*hin(5,n)
            end do ! k1
          end do ! n

          do k1 = 1 , 3
            do k2 = 1 , 3
              s(ni+k1,nj+k2) = s(ni+k1,nj+k2)
     &                       + gh(k1,1) * g(j,1,k2)
     &                       + gh(k1,2) * g(j,2,k2)
     &                       + gh(k1,3) * g(j,3,k2)
     &                       + gh(k1,4) * g(j,4,k2)
     &                       + gh(k1,5) * g(j,5,k2)
            end do ! k2
          end do ! k1

        end do ! j

      end do ! i

      end
  

      subroutine sh3fbend ( xjw  , sg  , eh12, xnu ,
     &                        f1   , f2  , y1  , y2  , xnt ,
     &                        ndf  , nst , cr  , b   , sm  ,
     &                        p    , s   )
c*********************************************************************
c      Description:   SH3FBEND is a subroutine which computes the
c                     stiffness and residual contributions of the
c                     bending field, using Pian-Sumihara
c                     stress interpolations for Shell.

c      Authors:       M.S. Rifai, J.C. Simo, & D.D. Fox.

c      Date:          February 1991.

c      Version:       This routine was tested in FEAP version 6.3.

c-----[--.----+----.----+----.-----------------------------------------]

c      Routine Input:
c      --------------
c      sg ............ Gauss point coordinates.
c      xjw ........... Reference jacobian weighted for integration.
c      eh12 .......... Membrane inverse constitutive coefficient.
c      xnu ........... Poisson's ratio.
c      f1,f2 ......... Tensor transformations in vector format
c                      at center of element.
c      y1,y2, xnt .... Commonly used parts of [H].
c      cr ............ Current bending strain measure.
c      B ............. Discrete strain-displacement operator.
c      ndf ........... Number of degrees of freedom per node.
c      nst ........... Number of degrees of freedom per element.

c      Routine Output:
c      ---------------
c      sm ............ Bending stress.
c      S ............. Element stiffness.
c      P ............. Element residual.

c*********************************************************************

      implicit  none

      integer   i      , j       , k       , l        , ndf
      integer   k1     , k2      , nst     , n
      integer   ni     , nj
      real*8    eta    , xi      , fac0    , fac1     , fac2
      real*8    fach1    , fach2 , hdet    , eh12     , xnu
      real*8    xjw(4) , sg(3,4) , f1(3)   , f2(3)
      real*8    xnt(13), cr(3,4) , sm(3,4) , htm(2,2) , hin(5,5)
      real*8    br(5)  , y1(4)   , y2(4)   , b(8,24,4), g(4,5,6)
      real*8    er(5)  , p(nst)  , s(nst,nst), gh(6,5)
      real*8    dot

      save

c     Zero Matrices

      do i = 1,5
        er(i) = 0.0d0
        do j = 1,5
          hin(j,i) = 0.0d0
        end do ! j
        do j = 1,4
          g(j,i,1) = 0.0d0
          g(j,i,2) = 0.0d0
          g(j,i,3) = 0.0d0
        end do ! j
      end do ! i

c       Gauss Point Loop

      do i = 1 , 4

c       Common Factors

        eta   = sg(2,i) - xnt(2)
        xi    = sg(1,i) - xnt(3)
        fac0  = xjw(i)
        fac1  = fac0 * eta
        fac2  = fac0 * xi

c       Modified Strains [er] = [S][cr]

        er(1) = er(1) + fac0 * cr(1,i)
        er(2) = er(2) + fac0 * cr(2,i)
        er(3) = er(3) + fac0 * cr(3,i)
        er(4) = er(4) + fac1 * dot ( f1 , cr(1,i) ,3 )
        er(5) = er(5) + fac2 * dot ( f2 , cr(1,i) ,3 )

c       Integrated Strain-Displacement Matrix [G] = [S][B]

        do k = 1 , 4
          do j = 1 , 6
            l = (k-1)*ndf + j
            g(k,1,j) = g(k,1,j) + fac0 * b(6,l,i)
            g(k,2,j) = g(k,2,j) + fac0 * b(7,l,i)
            g(k,3,j) = g(k,3,j) + fac0 * b(8,l,i)
            g(k,4,j) = g(k,4,j) + fac1 * dot(f1,b(6,l,i),3 )
            g(k,5,j) = g(k,5,j) + fac2 * dot(f2,b(6,l,i),3 )
          end do ! j
        end do ! k

      end do ! i

c     Compute Second Block of [H-inv]

      htm(1,1) = xnt(12) * dot ( f1 , y1 ,3 )
      htm(2,2) = xnt(11) * dot ( f2 , y2 ,3 )
      htm(1,2) = xnt(13) * dot ( f1 , y2 ,3 )
      hdet = 1.d0 / (htm(1,1)*htm(2,2) - htm(1,2)**2)

c     Constitutive Relations

c     Compute [H-inv]

      fach1    = eh12 / (xnt(1) * (1.d0 - xnu * xnu))
      fach2    = eh12 / (xnt(1) * (1.d0 + xnu) * 2.d0)
      hin(1,1) = fach1
      hin(1,2) = fach1 * xnu
      hin(2,1) = hin(1,2)
      hin(2,2) = hin(1,1)
      hin(3,3) = fach2
      hin(4,4) = htm(2,2) * hdet
      hin(5,5) = htm(1,1) * hdet
      hin(4,5) =-htm(1,2) * hdet
      hin(5,4) = hin(4,5)

c     Compute Stress Variables

      do i = 1,5
        br(i) = hin(i,1)*er(1)
     &        + hin(i,2)*er(2)
     &        + hin(i,3)*er(3)
     &        + hin(i,4)*er(4)
     &        + hin(i,5)*er(5)
      end do ! i

c     Gauss Point Loop

      do i = 1 , 4

c       Compute Stresses

        eta   = sg(2,i) - xnt(2)
        xi    = sg(1,i) - xnt(3)
        sm(1,i) = br(1) + eta * f1(1)*br(4) + xi * f2(1)*br(5)
        sm(2,i) = br(2) + eta * f1(2)*br(4) + xi * f2(2)*br(5)
        sm(3,i) = br(3) + eta * f1(3)*br(4) + xi * f2(3)*br(5)

c       Weight Stresses by Jacobian

        sm(1,i) = sm(1,i) * xjw(i)
        sm(2,i) = sm(2,i) * xjw(i)
        sm(3,i) = sm(3,i) * xjw(i)

      end do ! i

c     Node Loop #1

      do i = 1 , 4
        ni = (i-1)*ndf

c       Compute Residual

        do k = 1 , 6
          p(ni+k) = p(ni+k) - g(i,1,k)*br(1)
     &                      - g(i,2,k)*br(2)
     &                      - g(i,3,k)*br(3)
     &                      - g(i,4,k)*br(4)
     &                      - g(i,5,k)*br(5)
        end do ! k

c       Node Loop #2

        do j = 1 , 4
          nj = (j-1)*ndf

c         Compute Stiffness

          do n = 1,5
            do k1 = 1,6
              gh(k1,n) = g(i,1,k1)*hin(1,n)
     &                 + g(i,2,k1)*hin(2,n)
     &                 + g(i,3,k1)*hin(3,n)
     &                 + g(i,4,k1)*hin(4,n)
     &                 + g(i,5,k1)*hin(5,n)
            end do ! k1
          end do ! n

          do k1 = 1 , 6
            do k2 = 1 , 6
              s(ni+k1,nj+k2) = s(ni+k1,nj+k2)
     &                       + gh(k1,1) * g(j,1,k2)
     &                       + gh(k1,2) * g(j,2,k2)
     &                       + gh(k1,3) * g(j,3,k2)
     &                       + gh(k1,4) * g(j,4,k2)
     &                       + gh(k1,5) * g(j,5,k2)
            end do ! k2
          end do ! k1

        end do ! j

      end do ! i

      end

      subroutine sh3fbmtx ( shp1 , shp2 , shx1 , shx2 ,
     &                        cphm , cpx1 , cpx2 , cdrm ,
     &                        cdx1 , cdx2 , b           )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]

c      Description:    SH3FBMTX is the subroutine which computes the
c                      discrete strain-displacement operator (matrix)
c                      for the general shell inextensible element.
c                      The membrane and bending strains are assumed
c                      to be defined in a cartesian reference frame.
c                      The shear strains are computed in the natural
c                      frame using Bathe-Dvorkin interpolations.

c      Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.
c      Date:           January, 1991.
c      Revised:        R.L. Taylor  - - February 1997

c      Version:        This routine was tested in FEAP version 6.3
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      shp1,shp2 ..... Nodal shape function natural derivatives.
c      shx1,shx2 ..... Nodal shape function cartesian derivatives.
c      cpx1,cpx2 ..... Current coordinate global derivatives
c                      at the Gauss points.
c      cphm .......... Current coordinate local derivatives
c                      at the midside nodes.
c      cdx1,cdx2 ..... Current director global derivatives
c                      at the Gauss points.
c      cdrm .......... Current local directors
c                      at the midside nodes.
c      xln ........... Localized nodal orthogonal transformation
c                      matrices, xln(3,3,ElementNode).

c      Routine Output:
c      ---------------
c      b ............. Discrete strain-displacement operator.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i            , j        ,  nm      , nb
      integer   nsb1 (4)     , nsb2 (4) , nsm1 (4) , nsm2 (4)
      real*8    shp1 (4)     , shp2 (4) , shx1 (4) , shx2 (4)
      real*8    cphm (3,4)   , cpx1 (3) , cpx2 (3)
      real*8    cdrm (3,4)   , cdx1 (3) , cdx2 (3) , b(8,24)

      save

c     Shear Scatter Data:

      data      nsb1 / 2 , 2 , 3 , 3 /
      data      nsb2 / 4 , 3 , 3 , 4 /
      data      nsm1 / 2 , 2 , 4 , 4 /
      data      nsm2 / 1 , 3 , 3 , 1 /

c     [Bmm] Part (Membrane), [Bbb] Part (Bending), and
c     [Bbm] Part (Bending-Membrane Coupling Terms)

c     [Bsm] Part (Shear - Displacement)
c     [Bsb] Part (Shear - Rotation)

      nm = 0
      nb = 3
      do j = 1,4
        do i = 1,3
          b(1,nm+i) = shx1(j)*cpx1(i)
          b(2,nm+i) = shx2(j)*cpx2(i)
          b(3,nm+i) = shx2(j)*cpx1(i) + shx1(j)*cpx2(i)

          b(6,nm+i) = shx1(j)*cdx1(i)
          b(7,nm+i) = shx2(j)*cdx2(i)
          b(8,nm+i) = shx1(j)*cdx2(i) + shx2(j)*cdx1(i)

          b(6,nb+i) = b(1,nm+i)
          b(7,nb+i) = b(2,nm+i)
          b(8,nb+i) = b(3,nm+i)

          b(4,nm+i) = shp1(j) * cdrm(i,nsm1(j))
          b(5,nm+i) = shp2(j) * cdrm(i,nsm2(j))

          b(4,nb+i) = shp1(nsb1(j)) * cphm(i,nsm1(j))
          b(5,nb+i) = shp2(nsb2(j)) * cphm(i,nsm2(j))

        end do ! i
        nm = nm + 6
        nb = nb + 6
      end do ! j

      end

     
      subroutine sh3fstrn ( shp1 , shp2 , cphm , cpx1 , cpx2 , cdrm ,
     &                        cdx1 , cdx2 , ce   , cx   , cr          )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    SH3FSTRN is the subroutine which computes the
c                      strain measures for the inextensible shell
c                      element. The membrane and bending strains
c                      are defined in a cartesian reference frame.
c                      The shear strains are computed in the natural
c                      frame using Bathe-Dvorkin interpolations.

c      Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c      Date:           January 1991.

c      Version:        This routine was tested in FEAP version 6.3
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      cpx1,cpx2 ..... Current coordinate global derivatives
c                      at the Gauss points.
c      cphm .......... Current coordinate local derivatives
c                      at the midside nodes.
c      cdx1,cdx2 ..... Current director global derivatives
c                      at the Gauss points.
c      cdrm .......... Current local directors
c                      at the midside nodes.

c      Routine Output:
c      ---------------
c      ce............. Current membrane strain measure.
c      cx............. Current shear strain measure.
c      cr............. Current bending strain measure.
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      real*8    shp1 (4)   , shp2 (4)
      real*8    cphm (3,4) , cpx1 (3) , cpx2 (3)
      real*8    cdrm (3,4) , cdx1 (3) , cdx2 (3)
      real*8    ce   (3)   , cx   (2) , cr   (3)

      save

c     Compute Membrane Strains:

      ce(1) = (cpx1(1)*cpx1(1)+cpx1(2)*cpx1(2)+cpx1(3)*cpx1(3))*0.5d0
      ce(2) = (cpx2(1)*cpx2(1)+cpx2(2)*cpx2(2)+cpx2(3)*cpx2(3))*0.5d0
      ce(3) =  cpx1(1)*cpx2(1)+cpx1(2)*cpx2(2)+cpx1(3)*cpx2(3)

c     Compute Shear Strains:

      cx(1) = 2.d0 * ( shp1(2) * (cphm(1,2)*cdrm(1,2)
     &                         +  cphm(2,2)*cdrm(2,2)
     &                         +  cphm(3,2)*cdrm(3,2))
     &               + shp1(3) * (cphm(1,4)*cdrm(1,4)
     &                         +  cphm(2,4)*cdrm(2,4)
     &                         +  cphm(3,4)*cdrm(3,4)) )

      cx(2) = 2.d0 * ( shp2(4) * (cphm(1,1)*cdrm(1,1)
     &                         +  cphm(2,1)*cdrm(2,1)
     &                         +  cphm(3,1)*cdrm(3,1))
     &               + shp2(3) * (cphm(1,3)*cdrm(1,3)
     &                         +  cphm(2,3)*cdrm(2,3)
     &                         +  cphm(3,3)*cdrm(3,3)) )

c     Compute Bending Strains:

      cr(1) = cpx1(1)*cdx1(1) + cpx1(2)*cdx1(2) + cpx1(3)*cdx1(3)
      cr(2) = cpx2(1)*cdx2(1) + cpx2(2)*cdx2(2) + cpx2(3)*cdx2(3)
      cr(3) = cpx1(1)*cdx2(1) + cpx1(2)*cdx2(2) + cpx1(3)*cdx2(3)
     &      + cpx2(1)*cdx1(1) + cpx2(2)*cdx1(2) + cpx2(3)*cdx1(3)

      end

      subroutine sh3fcmtx ( d      , dthk , xjw  ,
     &                      optmix , zph1 , zph2 , c   )
c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    SH3FCMTX is the subroutine which sets-up the
c                      linear resultant constitutive relations.

c      Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c      Date:           January, 1991.

c      Version:        This routine was tested in FEAP version 6.01
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      d(*) .......... Material parameters from inmate
c      dthk .......... Element thickness
c      xjw ........... Jacobian times quadrature weight
c      optmix ........ Mixed treatment option:
c                      0 => Galerkin (displacement formulation),
c                      1 => Hellinger-Reissner (Pian-Summihara).
c      zph1,zph2 ..... Initial coordinate local derivatives.

c      Routine Output:
c      ---------------
c      c ............. Resultant stress constitutive matrix.
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      logical   optmix
      real*8    a11     , a22     , a12     , adet
      real*8    dthk    , dthk3   , xjw
      real*8    d(*)    , zph1(3) , zph2(3) , c(8,8) , ain(2,2)

c     Calculate Metric

      a11 = zph1(1)*zph1(1) + zph1(2)*zph1(2) + zph1(3)*zph1(3)
      a22 = zph2(1)*zph2(1) + zph2(2)*zph2(2) + zph2(3)*zph2(3)
      a12 = zph1(1)*zph2(1) + zph1(2)*zph2(2) + zph1(3)*zph2(3)

c     Calculate Inverse of Metric

      adet     = 1.d0/(a11*a22-a12**2)
      ain(1,1) = a22*adet
      ain(1,2) =-a12*adet
      ain(2,1) = ain(1,2)
      ain(2,2) = a11*adet

      if (.not.optmix) then

c       [Cm] Part (Membrane):

        c(1,1) = d(21)*dthk*xjw
        c(1,2) = d(24)*dthk*xjw
        c(1,3) = 0.d0
        c(2,2) = d(22)*dthk*xjw
        c(2,3) = 0.d0
        c(3,3) = d(27)*dthk*xjw
        c(2,1) = c(1,2)
        c(3,1) = c(1,3)
        c(3,2) = c(2,3)

c       [Cb] Part (Bending):

        dthk3  = dthk**3/12.d0*xjw
        c(6,6) = d(21)*dthk3
        c(6,7) = d(24)*dthk3
        c(6,8) = 0.d0
        c(7,7) = d(22)*dthk3
        c(7,8) = 0.d0
        c(8,8) = d(27)*dthk3
        c(7,6) = c(6,7)
        c(8,6) = c(6,8)
        c(8,7) = c(7,8)
      endif

c     [Cs] Part (Shear):

      dthk3  = d(37)*d(28)*dthk*xjw
      c(4,4) = dthk3 * ain(1,1)
      c(4,5) = dthk3 * ain(1,2)
      c(5,5) = dthk3 * ain(2,2)
      c(5,4) = c(4,5)

      end

      subroutine sh3fstrs ( c      , ce     , cx     , cr ,
     &                        optmix , sn     , sq     , sm )
c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    SH3FSTRS is a subroutine which computes the
c                      stress resultants, given the strain measures
c                      and a linear constitutive matrix.

c      Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c      Date:           January, 1991.

c      Version:        This routine was tested in FEAP version 6.3
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      c ............. Resultant stress constitutive matrix.
c      ce............. Current membrane strain measure.
c      cx............. Current shear strain measure.
c      cr............. Current bending strain measure.
c      optmix ........ Mixed treatment option:
c                      0 => Galerkin (displacement formulation),
c                      1 => Hellinger-Reissner (Pian-Summihara).

c      Routine Output:
c      ---------------
c      sn ............ Membrane stress.
c      sq ............ Shear stress.
c      sm ............ Bending stress.
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      logical   optmix
      real*8    c  (8,8) ,
     1            ce (3)   , cx (2) , cr (3) ,
     2            sn (3)   , sq (2) , sm (3)

      if (.not.optmix) then

c       Membrane Stresses:

        sn(1) = c(1,1)*ce(1) + c(1,2)*ce(2) + c(1,3)*ce(3)
        sn(2) = c(2,1)*ce(1) + c(2,2)*ce(2) + c(2,3)*ce(3)
        sn(3) = c(3,1)*ce(1) + c(3,2)*ce(2) + c(3,3)*ce(3)

c       Calculate Bending Stresses:

        sm(1) = c(6,6)*cr(1) + c(6,7)*cr(2) + c(6,8)*cr(3)
        sm(2) = c(7,6)*cr(1) + c(7,7)*cr(2) + c(7,8)*cr(3)
        sm(3) = c(8,6)*cr(1) + c(8,7)*cr(2) + c(8,8)*cr(3)
      endif

c     Calculate Shear Stresses:

      sq(1) = c(4,4)*cx(1) + c(4,5)*cx(2)
      sq(2) = c(5,4)*cx(1) + c(5,5)*cx(2)

      end

      subroutine sh3fener ( d, dthk, ce, cr, cx, sn, sq, sm,
     &                      ul, shp, fphi, xjw, ndf, lint)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Compute momenta and energy

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer         numnp,numel,nummat,nen,neq,ipr
      common /cdata/  numnp,numel,nummat,nen,neq,ipr

      real*8          xln
      real*8          rots       ,rvel       ,racc       ,thkl
      integer                                                     rotyp
      common /erotas/ xln(9,9,4),
     &                rots(3,9,2),rvel(3,9,2),racc(3,9,2),thkl(9),rotyp

      real*8          dm
      integer            n,ma,mct,iel,nel
      common /eldata/ dm,n,ma,mct,iel,nel

      logical         fl,    pfr
      common /fdata/  fl(12),pfr

      integer         ior,iow,ilg
      common /iofile/ ior,iow,ilg

      real*8          epl
      integer                  iepl,       neplts
      common /ptdat6/ epl(200),iepl(2,200),neplts

      real*8          ttim,dt,c1,c2,c3,c4,c5, chi
      common /tdata/  ttim,dt,c1,c2,c3,c4,c5, chi

      integer   i         , j          , lint       , ndf
      real*8    rhoa      , rhoi       , dthk       , dot
      real*8    d    (*)  , ul  (ndf,nen,*)         , shp (4,4)
      real*8    ce   (3,4), cx  (2,4)  , cr  (3,4)  , fphi(3,4)
      real*8    sn   (3,4), sq  (2,4)  , sm  (3,4)  , aa   (3)
      real*8    pg   (3)  , vd  (3)    , w1  (3)    , xjw (4)

      save

c     Galerkin and Shear Treatment

      rhoa = d( 4)*dthk
      rhoi = d(8)*rhoa*dthk**2/12.d0

      do i = 1,lint

c       Interpolate Velocities and Accelerations

        do j = 1 , 3
          pg(j)   = shp(1,i)*fphi(j, 1)  + shp(2,i)*fphi(j, 2)
     &            + shp(3,i)*fphi(j, 3)  + shp(4,i)*fphi(j, 4)

          vd(j)   = shp(1,i)*ul  (j,1,4) + shp(2,i)*ul  (j,2,4)
     &            + shp(3,i)*ul  (j,3,4) + shp(4,i)*ul  (j,4,4)

          w1(j)   = shp(1,i)*rvel(j,1,2) + shp(2,i)*rvel(j,2,2)
     &            + shp(3,i)*rvel(j,3,2) + shp(4,i)*rvel(j,4,2)
        end do ! j

        call vecp ( pg , vd , aa )

c       Integrate Momenta

        do j = 1 , 3
          epl(j  ) = epl(j  ) + aa(j) * rhoa * xjw(i)
          epl(j+3) = epl(j+3) + w1(j) * rhoi * xjw(i)
        end do ! j

c       Integrate Energy

        epl(7) = epl(7) + 0.5d0*(dot(vd,vd,3) * rhoa
     &                         + dot(w1,w1,3) * rhoi) * xjw(i)
        epl(8) = epl(8) + 0.5d0*(dot(sn(1,i),ce(1,i),3)
     &                         + dot(sq(1,i),cx(1,i),2)
     &                         + dot(sm(1,i),cr(1,i),3)) * xjw(i)
      end do ! i

      end

     
      subroutine sh3fgeo6 ( lint , shp1 , shp2 , shx1 , shx2 ,
     &                      cphm , cpx1 , cpx2 , dir  ,
     &                      sq   , sm   , s    )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FGEO6 adds additional terms
c                        static geometric tangent stiffness,
c                        in case of a 6 dof formulation.

c        Author:         N. Tarnow and J.C. Simo

c        Date:           February 1992

c       Version:        This routine was tested in FEAP version 6.01
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        shp1,shp2 ..... Nodal shape function natural derivatives.
c        shx1,shx2 ..... Nodal shape function cartesian derivatives.
c        cpx1,cpx2 ..... Current coordinate global derivatives
c                        at Gauss points.
c        cphm .......... Current coordinate local derivatives
c                        at midside nodes.
c        cdrm .......... Current local directors
c                        at midside nodes.
c        sn ............ Membrane stress.
c        sq ............ Shear stress.
c        sm ............ Bending stress.

c        Routine Output:
c        ---------------
c        s ............. Element stiffness.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer         ndof
      logical                 fsetr,frotas
      common /crotas/ ndof(9),fsetr,frotas

      integer   ni,k,l,i,lint
      real*8    shp1(4,4),shp2(4,4),shx1(4,4),shx2(4,4)
      real*8    cphm(3,4),cpx1(3,4),cpx2(3,4),dir(3,9)
      real*8    sq(2,4),sm(3,4),s(24,24),fac(4),cfac(4,3)
      real*8    xm(3,4),xn1(4),xn2(4),xn3(4),xn4(4), fact

      save

c     Evaluate Shape Functions:

      data      xn1  / 0.5d0, 0.0d0, 0.0d0, 0.5d0 /
      data      xn2  / 0.5d0, 0.5d0, 0.0d0, 0.0d0 /
      data      xn3  / 0.0d0, 0.5d0, 0.5d0, 0.0d0 /
      data      xn4  / 0.0d0, 0.0d0, 0.5d0, 0.5d0 /

c     Loop over Gauss Points

      do k = 1,4
        fac(k) = 0.0d0
      end do ! k
      do l = 1,lint
        fac(1) = fac(1) + 2.d0 * shp1(2,l) * sq(1,l)
        fac(2) = fac(2) + 2.d0 * shp1(3,l) * sq(1,l)
        fac(3) = fac(3) + 2.d0 * shp2(4,l) * sq(2,l)
        fac(4) = fac(4) + 2.d0 * shp2(3,l) * sq(2,l)
      end do ! l

      do k = 1,3
        cfac(1,k) = fac(1)*cphm(k,2)
        cfac(2,k) = fac(2)*cphm(k,4)
        cfac(3,k) = fac(3)*cphm(k,1)
        cfac(4,k) = fac(4)*cphm(k,3)
      end do ! k

c     Loop over Nodes

      ni = 0
      do i = 1,4
        if(ndof(i).eq.6)then
          do k = 1,3
            xm(k,i) = cfac(1,k)*xn2(i) + cfac(2,k)*xn4(i)
     &              + cfac(3,k)*xn1(i) + cfac(4,k)*xn3(i)
          end do ! k

c         Loop over Gauss Points

          do l = 1,lint
            do k = 1,3
              xm(k,i) = xm(k,i) + sm(1,l)*cpx1(k,l)*shx1(i,l)
     &                          + sm(2,l)*cpx2(k,l)*shx2(i,l)
     &                          + sm(3,l)*cpx1(k,l)*shx2(i,l)
     &                          + sm(3,l)*cpx2(k,l)*shx1(i,l)
            end do ! k
          end do ! l

c         Compute Projection

          fact = xm(1,i)*dir(1,i) + xm(2,i)*dir(2,i) + xm(3,i)*dir(3,i)
          do k = 1,3
            xm(k,i) = xm(k,i) - fact*dir(k,i)
          end do ! k

c         Add Contributions to Element Stiffness

          do k = 1,3
            do l = 1,3
              s(ni+3+k,ni+3+l) = s(ni+3+k,ni+3+l)
     &                         + 0.5d0*(xm(k,i)*dir(l,i)
     &                                + xm(l,i)*dir(k,i))
            end do ! l
          end do ! k
        endif
        ni = ni + 6
      end do ! i

      end

     
      subroutine sh3fgeom ( shp1 , shp2 , shx1 , shx2 ,
     &                        cphm , cpx1 , cpx2 , dir  ,
     &                        sn   , sq   , sm   , s           )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    SH3FGEOM is the subroutine which constructs
c                      the static geometric tangent stiffness.

c      Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c      Date:           January 1991.

c      Version:        This routine was tested in FEAP version 6.3
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      shp1,shp2 ..... Nodal shape function natural derivatives.
c      shx1,shx2 ..... Nodal shape function cartesian derivatives.
c      cpx1,cpx2 ..... Current coordinate global derivatives
c                      at the Gauss points.
c      cphm .......... Current coordinate local derivatives
c                      at the midside nodes.
c                      at the midside nodes.
c      sn ............ Membrane stress.
c      sq ............ Shear stress.
c      sm ............ Bending stress.

c      Routine Output:
c      ---------------
c      s ............. Element stiffness.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i  , j  , k  , ni , nj
      integer   ns1  (4), ns2  (4),  ns3  (4), ns4  (4)
      real*8    termb   , termm   , termd   , sq1h    , sq2h
      real*8    sm1i    , sm2i    , sn1i    , sn2i
      real*8    shp1 (4), shp2 (4), shx1 (4), shx2 (4), dir(3,9)
      real*8    cphm (3,4)   , cpx1 (3)     , cpx2 (3)
      real*8    sn   (3)     , sq   (2)     , sm   (3)
      real*8    s    (24,24) , qt   (4,4)   , qd   (4)

c     Scatter Data:

      data      ns1  / 2 , 2 , 4 , 4 /
      data      ns2  / 1 , 3 , 3 , 1 /
      data      ns3  / 1 , 1 , 4 , 4 /
      data      ns4  / 1 , 2 , 2 , 1 /

c     Set Up Shear Terms:

      sq1h    = 0.5d0 * sq(1)
      sq2h    = 0.5d0 * sq(2)

      qt(1,1) = shp1(1)*sq1h + shp2(1)*sq2h
      qt(1,2) = shp1(1)*sq1h
      qt(1,3) = 0.0d0
      qt(1,4) = shp2(1)*sq2h

      qt(2,1) = shp1(2)*sq1h
      qt(2,2) = shp1(2)*sq1h + shp2(2)*sq2h
      qt(2,3) = shp2(2)*sq2h
      qt(2,4) = 0.0d0

      qt(3,1) = 0.0d0
      qt(3,2) = shp2(3)*sq2h
      qt(3,3) = shp1(3)*sq1h + shp2(3)*sq2h
      qt(3,4) = shp1(3)*sq1h

      qt(4,1) = shp2(4)*sq2h
      qt(4,2) = 0.0d0
      qt(4,3) = shp1(4)*sq1h
      qt(4,4) = shp1(4)*sq1h + shp2(4)*sq2h

      do i = 1,4
        qd(i)   = sq(1)*shp1(ns3(i)) * (cphm(1,ns1(i))*dir(1,i)
     &                               +  cphm(2,ns1(i))*dir(2,i)
     &                               +  cphm(3,ns1(i))*dir(3,i))
     &          + sq(2)*shp2(ns4(i)) * (cphm(1,ns2(i))*dir(1,i)
     &                               +  cphm(2,ns2(i))*dir(2,i)
     &                               +  cphm(3,ns2(i))*dir(3,i))
      end do ! i

      ni = 0
      do i = 1,4
        sm1i= sm(1)*shx1(i) + sm(3)*shx2(i)
        sm2i= sm(2)*shx2(i) + sm(3)*shx1(i)

        sn1i= sn(1)*shx1(i) + sn(3)*shx2(i)
        sn2i= sn(2)*shx2(i) + sn(3)*shx1(i)

        nj = 0
        do j = 1,4

c         Bending Term

          termb = sm1i*shx1(j) + sm2i*shx2(j)

c         Membrane Term

          termm = sn1i*shx1(j) + sn2i*shx2(j)

c         Assemble Membrane, Bending and Shear Parts

          do k = 1,3
            s(ni+k,nj+k  ) = s(ni+k,nj+k  ) + termm
            s(ni+k,nj+3+k) = s(ni+k,nj+3+k) + termb + qt(i,j)
            s(ni+3+k,nj+k) = s(ni+3+k,nj+k) + termb + qt(j,i)
          end do ! k
          nj = nj + 6
        end do ! j

c       Diagonal Bending Term

        termd = (sm1i*cpx1(1) + sm2i*cpx2(1))*dir(1,i)
     &        + (sm1i*cpx1(2) + sm2i*cpx2(2))*dir(2,i)
     &        + (sm1i*cpx1(3) + sm2i*cpx2(3))*dir(3,i)

c       Assemble Diagonal Bending and Shear Parts

        do k = 1,3
          s(ni+3+k,ni+3+k) = s(ni+3+k,ni+3+k) + qd(i) - termd
        end do ! k

        ni = ni + 6
      end do ! l

      end

     

      subroutine sh3finte (xl,xln,ndof,ndm,nel, dir, xi, xj)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c        Description:   SH3FINTE is subroutine which initializes
c                       appropriately local element director field
c                       according on number of degrees of freedom
c                       of nodes in element; as defined by the
c                       integer array ndof(nel). If a node posseses 6
c                       dof's if indicates presence of a shell
c                       intersection; otherwise, node is associated
c                       with a smooth shell mid-surface.

c        Author:        M.S. Rifai J.C. Simo D.D. Fox

c        Date:          January 1992.
c        Revised:       February 1997: reduce to 1-call

c        Version:       This routine was tested in FEAP version 6.3
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        xl(ndm,*) ..... Local nodal coordinates of element
c        xln(3,3,9,4) .. local nodal rotation matrices
c                        xln(3,3,9,1): time t_n
c                        xln(3,3,9,2): time t_n+a
c                        xln(3,3,9,3): time t_n+1
c                        xln(3,3,9,4): time t_0
c        ndof(nel) ..... Array containing number of DOF/node
c        ndm ........... Spacial dimension (ndm=3)
c        nel ........... Number of nodes in element

c        Routine Output:
c        ---------------
c        dir(3,9,4) .... Local element director field
c                        dir(3,9,1): time t_n
c                        dir(3,9,2): time t_n+a
c                        dir(3,9,3): time t_n+1
c                        dir(3,9,4): time t_0
c        xi(3,3,9) ..... Transformation matrix for nodal DOF at t_n+a
c        xj(3,3,9) ..... Transformation matrix for nodal DOF at t_n+1
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8          theta
      integer                  nrk,nrc,nrm,nrt,noi,nt,nrkn,nrcn,nrmn
      common /ddata/  theta(3),nrk,nrc,nrm,nrt,noi,nt,nrkn,nrcn,nrmn

      logical         fl,    pfr
      common /fdata/  fl(12),pfr

      integer   i, j, k, l, m, node, ndm, nel, ntime
      integer   ii(4), jj(4), kk(4), ndof(nel)
      real*8    xnorm, facn, fac1,tt1(3), tt2(3), tt3(3), dir(3,9,4)
      real*8    xl(ndm,*), xln(3,3,9,4), xi(3,3,9), xj(3,3,9)

      save

c     Scatter Data (nel = 4)

      data      ii / 1 , 2 , 3 , 4 /
      data      jj / 2 , 3 , 4 , 1 /
      data      kk / 3 , 4 , 1 , 2 /

c     Define local element directors

      fac1 = theta(3)
      facn = 1.d0 - fac1

      do j = 1,4
        node = jj(j)

        if    (ndof(node) .eq. 5) then
          do ntime = 1,4
            do l = 1,3
              dir(l,node,ntime) = xln(l,3,node,ntime)
            end do ! l
          end do ! ntime

c         Interpolate director to t_n+alpha
          if(fl(9))then
            do l=1,3
              dir(l,node,2) = facn*dir(l,node,1) + fac1*dir(l,node,3)
            end do ! l
          endif

        elseif(ndof(node) .eq. 6) then

          i = ii(j)
          k = kk(j)

c         Compute reference directors

          do l = 1,3
            tt1(l) = xl(l,k) - xl(l,node)
            tt2(l) = xl(l,i) - xl(l,node)
          end do ! l
          call vecp( tt1, tt2, tt3 )
          xnorm = 1.d0/sqrt( tt3(1)**2 + tt3(2)**2 + tt3(3)**2 )
          do l = 1,3
            dir(l,node,4) = tt3(l) * xnorm
          end do ! l

c         Compute current directors

          do l = 1,3

            tt3(l) = xln(1,l,node,4)*dir(1,node,4)
     &             + xln(2,l,node,4)*dir(2,node,4)
     &             + xln(3,l,node,4)*dir(3,node,4)

            dir(l,node,1) = xln(l,1,node,1)*tt3(1)
     &                    + xln(l,2,node,1)*tt3(2)
     &                    + xln(l,3,node,1)*tt3(3)

            dir(l,node,3) = xln(l,1,node,3)*tt3(1)
     &                    + xln(l,2,node,3)*tt3(2)
     &                    + xln(l,3,node,3)*tt3(3)

          end do ! l
          if(fl(9))then
            do l = 1,3
              dir(l,node,2) = facn*dir(l,node,1) + fac1*dir(l,node,3)
            end do ! l
          else
            do l = 1,3
              dir(l,node,2) = xln(l,1,node,2)*tt3(1)
     &                      + xln(l,2,node,2)*tt3(2)
     &                      + xln(l,3,node,2)*tt3(3)
            end do ! l
          endif
        endif
      end do ! j

c     Compute transformation matrix

      do node = 1,4

        if (ndof(node) .eq. 5) then
          do m = 1,3
            do l = 1,3
              xi(l,m,node) = xln(l,m,node,2)
              xj(l,m,node) = xln(l,m,node,3)
            end do ! l
          end do ! m

        elseif (ndof(node) .eq. 6) then
          xi(1,1,node) =  0.d0
          xi(2,2,node) =  0.d0
          xi(3,3,node) =  0.d0
          xi(1,2,node) =  dir(3,node,2)
          xi(2,3,node) =  dir(1,node,2)
          xi(3,1,node) =  dir(2,node,2)
          xi(2,1,node) = -dir(3,node,2)
          xi(3,2,node) = -dir(1,node,2)
          xi(1,3,node) = -dir(2,node,2)

          xj(1,1,node) =  0.d0
          xj(2,2,node) =  0.d0
          xj(3,3,node) =  0.d0
          xj(1,2,node) =  dir(3,node,3)
          xj(2,3,node) =  dir(1,node,3)
          xj(3,1,node) =  dir(2,node,3)
          xj(2,1,node) = -dir(3,node,3)
          xj(3,2,node) = -dir(1,node,3)
          xj(1,3,node) = -dir(2,node,3)
        endif

      end do ! node

      end

      
      
        subroutine sh3fintr ( cphi , dir , cphm , cdrm )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FINTR is the subroutine which interpolates
c                        the midsurface derivatives and directors
c                        at the midside points.

c        Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c        Date:           January 1991.

c        Version:        This routine was tested in FEAP version 6.01
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        cphi .......... Current local nodal coordinates.
c        dir ........... Localized nodal director
c                        vectors, dir(3,ElementNode).

c        Routine Output:
c        ---------------
c        cphm .......... Current coordinate local derivatives
c                        at the midside nodes.
c        cdrm .......... Initial local directors
c                        at the midside nodes.
c-----[--.----+----.----+----.-----------------------------------------]

        implicit  none

        integer   i

        real*8    cphi (3,*) , cphm (3,4) , dir(3,9) , cdrm (3,4)

        save

c       Interpolate Midsurface Positions

        do i = 1,3
          cphm(i,1) = ( cphi(i,4) - cphi(i,1) ) * 0.5d0
          cphm(i,2) = ( cphi(i,2) - cphi(i,1) ) * 0.5d0
          cphm(i,3) = ( cphi(i,3) - cphi(i,2) ) * 0.5d0
          cphm(i,4) = ( cphi(i,3) - cphi(i,4) ) * 0.5d0
        end do ! i

c       Interpolate Directors to Midsurface

        do i = 1,3
          cdrm(i,1) = ( dir(i,4) + dir(i,1) ) * 0.5d0
          cdrm(i,2) = ( dir(i,1) + dir(i,2) ) * 0.5d0
          cdrm(i,3) = ( dir(i,2) + dir(i,3) ) * 0.5d0
          cdrm(i,4) = ( dir(i,3) + dir(i,4) ) * 0.5d0
        end do ! i

        end

        subroutine sh3fshap ( xi , shp , shp1 , shp2 )
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FSHAP is the subroutine which computes the
c                        element shape functions and local derivatives
c                        for the 4-noded quadrilateral element.
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        xi ............ Local coordinates of the current point.

c        Routine Output:
c        ---------------
c        shp ........... Nodal shape function values.
c        shp1,shp2 ..... Nodal shape function natural derivatives.
c-----[--.----+----.----+----.-----------------------------------------]

        implicit  none

        real*8    xi(2) , shp(4) , shp1(4) , shp2(4)

c       Evaluate (xi) Derivative:

        shp1(2) =  0.25d0 * (1.d0-xi(2))
        shp1(3) =  0.25d0 * (1.d0+xi(2))
        shp1(4) = -shp1(3)
        shp1(1) = -shp1(2)

c       Evaluate (xi) Derivative:

        shp2(3) =  0.25d0 * (1.d0+xi(1))
        shp2(4) =  0.25d0 * (1.d0-xi(1))
        shp2(1) = -shp2(4)
        shp2(2) = -shp2(3)

c       Evaluate Shape Function:

        shp (1) =  shp2(4) * (1.d0-xi(2))
        shp (2) =  shp2(3) * (1.d0-xi(2))
        shp (3) =  shp2(3) * (1.d0+xi(2))
        shp (4) =  shp2(4) * (1.d0+xi(2))

c       Exit

        end

        subroutine sh3fmdsh ( shp1  , shp2 , xjinv , shx1  , shx2 )
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FMDSH is a subroutine which transforms local
c                        shape function gradients into global
c                        (Cartesian) gradients.
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        shp1,shp2 ..... Nodal shape function natural derivatives.
c        xjinv ......... Local-global coordinate Jacobian matrix.

c        Routine Output:
c        ---------------
c        shx1,shx2 ..... Nodal shape function cartesian derivatives.
c-----[--.----+----.----+----.-----------------------------------------]

        implicit  none

        integer   j

        real*8    shp1  (4)   , shp2 (4) ,
     &            shx1  (4)   , shx2(4)  ,
     &            xjinv (2,2)

c       Transform Shape Function Derivatives:

        do j = 1 , 4
           shx1(j) = xjinv(1,1)*shp1(j) + xjinv(2,1)*shp2(j)
           shx2(j) = xjinv(1,2)*shp1(j) + xjinv(2,2)*shp2(j)
        end do ! j

c       Exit Subroutine.

        end

        subroutine sh3fresd ( b      , sn     , sq     , sm ,
     2                        optmix , p                    )
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FRESD is the routine which computes the
c                        static residual.
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        B ............. Discrete strain-displacement operator.
c        sn ............ Membrane stress.
c        sq ............ Shear stress.
c        sm ............ Bending stress.
c        optmix ........ Mixed treatment option:
c                        0 => Galerkin (displacement formulation),
c                        1 => Hellinger-Reissner (Pian-Summihara).

c        Routine Output:
c        ---------------
c        P ............. Element residual.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      logical    optmix
      integer    i
      real*8     b(8,24) , sn(3) , sq(2) , sm(3) , p(24)

      if (.not.optmix) then

c       Calculate Membrane, Shear and Bending Parts

        do i = 1,24
          p(i) = p(i) - b(1,i)*sn(1) - b(2,i)*sn(2) - b(3,i)*sn(3)
     &                - b(4,i)*sq(1) - b(5,i)*sq(2)
     &                - b(6,i)*sm(1) - b(7,i)*sm(2) - b(8,i)*sm(3)
        end do ! i

      else

c       Calculate Shear Part

        do i = 1,24
          p(i) = p(i) - b(4,i)*sq(1) - b(5,i)*sq(2)
        end do ! i

      endif


      end

      subroutine sh3fstif ( b      , b1     , c      ,
     &                       optmix , s      )
c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    SH3FSTIF is the routine which computes the
c                      material tangent, i.e., [B-t][C][B].
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      B ............. Discrete strain-displacement operator.
c      C ............. Resultant stress constitutive matrix.
c      OptMix ........ Mixed treatment option:
c                      0 => Galerkin (displacement formulation),
c                      1 => Hellinger-Reissner (Pian-Summihara).

c      Routine Output:
c      ---------------
c      s ............. Element stiffness.
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      logical   optmix
      integer   i       , j

      real*8    b(8,24) , b1(8,24) ,c(8,8) , s(24,24) , temp(24,8)

      if (.not.optmix) then

c       Calculate Membrane [B-transpose] . [C]

        do j = 1,3
          do i = 1,24
            temp(i,j) = b(1,i)*c(1,j) + b(2,i)*c(2,j) + b(3,i)*c(3,j)
          end do ! i
        end do ! j

c       Calculate Bending [B-transpose] . [C]

        do j = 6,8
          do i = 1,24
            temp(i,j) = b(6,i)*c(6,j) + b(7,i)*c(7,j) + b(8,i)*c(8,j)
          end do ! i
        end do ! j
      endif

c     Calculate Shear [B-transpose] . [C]

      do j = 4,5
        do i = 1,24
          temp(i,j) = b(4,i)*c(4,j) + b(5,i)*c(5,j)
        end do ! i
      end do ! j

c     Calculate [B-transpose] . [C] . [B]

      if (.not.optmix) then

        do j = 1,24
          do i = 1,24
            s(i,j) = s(i,j) + temp(i,1)*b1(1,j) + temp(i,2)*b1(2,j)
     &                      + temp(i,3)*b1(3,j) + temp(i,4)*b1(4,j)
     &                      + temp(i,5)*b1(5,j) + temp(i,6)*b1(6,j)
     &                      + temp(i,7)*b1(7,j) + temp(i,8)*b1(8,j)
          end do ! i
        end do ! j

      else

c       Calculate Shear [B-transpose] . [C] . [B]

        do j = 1,24
          do i = 1,24
            s(i,j) = s(i,j) + temp(i,4)*b1(4,j) + temp(i,5)*b1(5,j)
          end do ! i
        end do ! j

      endif

      end

      subroutine sh3fplst ( sn , sq , sm  ,
     &                        dt , st , shp , xjw     )
c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    SH3FPLST is the stress projection subroutine
c                      that computes the element contributions to the
c                      global equation system of the least-squares
c                      problem. The `mass' matrix is lumped with
c                      row-sum.
c                      The stresses are arranged as:
c                      1-3 := Membrane resultants.
c                      4-5 := Shear resultants.
c                      6-8 := Bending resultants.

c      Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c      Date:           January, 1991.
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      sn ............ Membrane stress.
c      sq ............ Shear stress.
c      sm ............ Bending stress.
c      shp ........... Nodal shape function values.
c      xjw ........... Reference jacobian weighted for integration.

c      Routine Output:
c      ---------------
c      dt ............ R.H.S. of stress projection equation system.
c      st ............ L.H.S. of stress projection equation system,
c                      lumped by row-sum.
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      integer         numnp,numel,nummat,nen,neq,ipr
      common /cdata/  numnp,numel,nummat,nen,neq,ipr

      integer         istv, iste, istp
      common /strnum/ istv, iste, istp

      integer   i
      real*8    xjw        , xjshp

      real*8    sn (3)     , sq(2)       , sm(3) ,
     &          dt (*)     , st(nen,*)   , shp(4)

c     Loop over Nodes:

      do i = 1,4
        xjshp      = xjw*shp(i)
        dt(i)   = dt(i) + xjshp

c       Membrane:

        st(i,1) = st(i,1) + sn(1)*xjshp
        st(i,2) = st(i,2) + sn(2)*xjshp
        st(i,3) = st(i,3) + sn(3)*xjshp

c       Shear:

        st(i,4) = st(i,4) + sq(1)*xjshp
        st(i,5) = st(i,5) + sq(2)*xjshp

c       Bending:

        st(i,6) = st(i,6) + sm(1)*xjshp
        st(i,7) = st(i,7) + sm(2)*xjshp
        st(i,8) = st(i,8) + sm(3)*xjshp
      end do ! i

      iste = 8

      end


      subroutine sh3flmda ( v , xlm )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Description:    SH3FLMDA is a subroutine which computes a unique
c                      orthogonal transformation matrix which rotates
c                      the 3-rd canonical basis vector, E3 into the
c                      vector v --- without drill.
c                      The singularity of the formula at (or near) the
c                      case: v = -E3 is avoided by computing the matrix
c                      which transforms -E3 into v when v.E3 < 0, and
c                      then applying a 180-degree rotation to it.

c      Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c      Date:           January 1991.

c      Version:        This routine was tested in FEAP version 6.3
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      v ............. Arbitrary vector in R-3.

c      Routine Output:
c      ---------------
c      xlm ........... Orthogonal transformation matrix, which
c                      rotates E3 into v without drill.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    fac, v(3) , xlm(3,2)

c     Compute [XLm], when v.E3 > 0

      if (v(3).gt.0.d0) then
         fac      =   1.d0 / ( 1.d0 + v(3) )
         xlm(1,1) =   v(3) + fac * v(2) * v(2)
         xlm(1,2) =        - fac * v(2) * v(1)
         xlm(2,1) =        - fac * v(1) * v(2)
         xlm(2,2) =   v(3) + fac * v(1) * v(1)
         xlm(3,1) = - v(1)
         xlm(3,2) = - v(2)

c     Compute [XLm], when v.E3 < 0

      else
         fac      =   1.d0 / ( 1.d0 - v(3) )
         xlm(1,1) = - v(3) + fac * v(2) * v(2)
         xlm(1,2) =          fac * v(2) * v(1)
         xlm(2,1) =        - fac * v(1) * v(2)
         xlm(2,2) =   v(3) - fac * v(1) * v(1)
         xlm(3,1) =   v(1)
         xlm(3,2) = - v(2)
      endif

      end

      
      subroutine sh3fmasg( d, dthk, shp, xjw, s, lint, nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FMASG is the subroutine which computes the
c                        mass matrix for rigid body and eigen analyses.
c        Authors:        R. Taylor
c        Date:           March, 1995.
c        Version:        This routine tested in FEAP version 6.01
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i , j , k , l , j1 , k1 , lint , nst
      real*8    rhoa  , rhoi  , dthk    , fac  , fac1
      real*8    d(*)  , s(nst,*) , shp(4,4), xjw(4)

c     Compute Mass Array:

      do i = 1,lint

c       Compute Weighted Density/Inertia

        rhoa = d(4)*dthk * xjw(i)
        rhoi = d(8)*rhoa*dthk**2/12.d0

c       Compute Mass

        j1 = 0
        do j = 1 , 4

c         Translational Mass

          fac = rhoa*shp(j,i)
          k1 = 0
          do k = 1 , 4
            fac1 = fac*shp(k,i)
            do l = 1 , 3
              s(j1+l,k1+l) = s(j1+l,k1+l) + fac1
            end do ! l
            k1   = k1 + 6
          end do ! k

c         Rotational Mass

          fac = rhoi*shp(j,i)
          do l = 4 , 6
            s(j1+l,j1+l) = s(j1+l,j1+l) + fac
          end do ! l

          j1 = j1 + 6
        end do ! j

      end do ! i

      end

      
      
      subroutine sh3fstre ( xl, shp, ce, cr, cx, sn, sq, sm, xjw,
     &                      lint, ndm, isw )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Stress output routine

c      Inputs:
c         xl(3,*)  - Nodal coordinates
c         shp(4,*) - Element shape functions
c         ce(3,*)  - Membrane strains
c         cr(3,*)  - Bending  strains
c         cx(2,*)  - Shear    strains
c         sn(3,*)  - Membrane stresses
c         sm(3,*)  - Bending  stresses
c         sq(3,*)  - Shear    stresses
c         xjw(*)   - Jacobian weight
c         lint     - Number of output points
c         ndm      - Mesh spatial dimension
c         isw      - Output switch: 4 = print; 8 = project.

c      Outputs:
c         Prints or projected values (in hr array)
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      integer         numnp,numel,nummat,nen,neq,ipr
      common /cdata/  numnp,numel,nummat,nen,neq,ipr

      real*8          dm
      integer            n,ma,mct,iel,nel
      common /eldata/ dm,n,ma,mct,iel,nel

      integer         ior,iow,ilg
      common /iofile/ ior,iow,ilg

      integer          num_nps
      parameter       (num_nps = 400)

      integer          np
      common /pointer/ np(num_nps)

      integer         nph,ner
      real*8                  erav,j_int   ,jshft
      common /prstrs/ nph,ner,erav,j_int(3),jshft

      integer         istv, iste, istp
      common /strnum/ istv, iste, istp

      real*8          hr
      integer               mr
      common          hr(1),mr(1)

      integer   i, j, lint, ndm, isw

      real*8    xl(ndm,*), shp(4,4), zphg(3,4)
      real*8    ce(3,4)  , cx(2,4) , cr(3,4)
      real*8    sn(3,4)  , sq(2,4) , sm(3,4) , xjw(4)

      save

c     Output Element Stresses:

      if(isw.eq.4) then
        mct = mct - 1
        if (mct.le.0) then
          if (ior.lt.0) write (*,4000)
          write (iow,4000)
          mct = 12
        endif

c       Gauss Loop

        do i = 1, lint

c         Interpolate Position

          do j = 1, 3
            zphg(j,i) = shp(1,i)*xl(j,1) + shp(2,i)*xl(j,2)
     &                + shp(3,i)*xl(j,3) + shp(4,i)*xl(j,4)
          end do ! i

c         Output to Screen

          if (ior.lt.0) then
            write (*,4100)  n, i, (zphg(j,i),j=1,3),
     &                     (ce(j,i),j=1,3), (sn(j,i),j=1,3),
     &                     (cr(j,i),j=1,3), (sm(j,i),j=1,3),
     &                     (cx(j,i),j=1,2), (sq(j,i),j=1,2)
          endif

c         Output to File

          write (iow,4100)  n, i, (zphg(j,i),j=1,3),
     &                     (ce(j,i),j=1,3), (sn(j,i),j=1,3),
     &                     (cr(j,i),j=1,3), (sm(j,i),j=1,3),
     &                     (cx(j,i),j=1,2), (sq(j,i),j=1,2)

        end do ! i

c       Blank record between elements

        if (ior.lt.0) write (*,*) ' '
        write (iow,*) ' '


c     Compute Nodal Stresses:

      elseif(isw.eq.8) then

c       Set Number of Printed/Plotted Stresses to 8.

        istv = 8

c       Project Stresses

        do i = 1,lint
          call sh3fplst ( sn(1,i)   , sq(1,i)   , sm(1,i),
     &                    hr(np(35)), hr(np(36)), shp(1,i), xjw(i) )
        end do ! i

      endif

c     Element Stress Format:

4000  format(/
     & ' ----------------------------------------',
     & '------------------------------------- '/,
     & '  Dynamic Finite-Deformation Elastic ',
     & '4-Node Shell Element -- Gauss Point Info  '/,
     & ' ----------------------------------------',
     & '------------------------------------- '//,
     & '  Elment Num   GaussPnt #   Mid-surf-X ',
     & '  Mid-surf-Y   Mid-surf-Z'/,
     & '  MemStrn-xx   MemStrn-yy   MemStrn-zz ',
     & '  MemStrs-xx   MemStrs-yy   MemStrs-zz  '/,
     & '  BndStrn-xx   BndStrn-yy   BndStrn-zz ',
     & '  BndStrs-xx   BndStrs-yy   BndStrs-zz  '/,
     & '  ShrStrn-1    ShrStrn-2    ShrStrs-1  ',
     & '  ShrStrs-2                             '/,
     & ' ------------ ------------ ------------',
     & ' ------------ ------------ ------------ '/)

4100  format(i7,12x,i1,6x,3e13.6/(6e13.6))

      end

     

      subroutine sh3fsurf(d,xl,ul,xjw, ndf,ndm,nst, p,s)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Nodal force and tangent array for pressure loading

c      INPUT variables
c        d(10)      Value of constant pressure on face
c        xl(4,*)    Nodal coordinates
c        ul(4,*)    Nodal displacements
c        xjw(*)     Jacobian times weight
c        ndf        Number of DOF / node
c        ndm        Space dimension
c        nst        Dimension of residual vector

c      OUTPUT variables
c        p(nst)     Contribution to residual
c        s(nst,nst) Contribution to striffness matrix

c      PARAMATER set up:

c         ma   =  0  (for 3-d analysis in reference coordinates)
c              =  1  (for 3-d analysis in current coordinates;
c                     an unsymmetric tangent matrix also computed.)

c         nel  =  4  (for 4-node surface of bi-linear element)

c                        4                3
c                        o----------------o
c                        |       A        |
c                        |       |eta     |
c                        |       |        |
c                        |       +---> xi |
c                        |                |
c                        |                |
c                        |                |
c         Nodes are:     o----------------o
c                        1                2

c       Author:         J.C. Simo

c       Revised:        R.L. Taylor  - - February 1997

c       Version:        This routine was tested in FEAP version 6.01

c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      real*8          theta
      integer                  nrk,nrc,nrm,nrt,noi,nt,nrkn,nrcn,nrmn
      common /ddata/  theta(3),nrk,nrc,nrm,nrt,noi,nt,nrkn,nrcn,nrmn

      real*8          dm
      integer            n,ma,mct,iel,nel
      common /eldata/ dm,n,ma,mct,iel,nel

      real*8          bpr   ,ctan   ,psil
      common /eltran/ bpr(3),ctan(3),psil

      integer   l, lint, ndf, ndm, nst
      integer   ii, jj, i1, j1
      real*8    pn, pp, fact

      real*8    sg(3,4), shp(3,4), xl(ndm,*), xu(3,4), xjw(*)
      real*8    d(*), ul(ndf,*), p(ndf,*), s(nst,*), dx(3,2)
      real*8    xii(4), eti(4), bg(3),bl(3)

      save

      data      xii / -0.5d0, 0.5d0, 0.5d0,-0.5d0 /
      data      eti / -0.5d0,-0.5d0, 0.5d0, 0.5d0 /
      data      bl  /  3*0.0d0 /

c     Set body loading factors

      call sbodyf(d, bg)

      if(int(d(76)).gt.0) then
        bl(3) = d(10)
      else
        bl(3) = d(10)*dm
      endif

c     Compute nodal coordinates in correct reference frame

      if(d(68).eq.0.0d0) then
        fact = 0.d0
      else
        fact = theta(3)
      endif
      do ii = 1,4
        xu(1,ii) = xl(1,ii) + fact*ul(1,ii)
        xu(2,ii) = xl(2,ii) + fact*ul(2,ii)
        xu(3,ii) = xl(3,ii) + fact*ul(3,ii)
      end do ! ii

c     Get quadrature information

      l  = 2
      call int2d (l, lint, sg)

c     First loop over quadrature points

      do l = 1,lint

c       Compute shape functioins and geometric factors

        do ii = 1,4
          shp(1,ii) = xii(ii)*(0.5d0+eti(ii)*sg(2,l))
          shp(2,ii) = eti(ii)*(0.5d0+xii(ii)*sg(1,l))
          shp(3,ii) = (0.5d0+xii(ii)*sg(1,l))*(0.5d0 + eti(ii)*sg(2,l))
        end do ! ii

        do ii = 1,3
          dx(ii,1) = shp(1,1)*xu(ii,1) + shp(1,2)*xu(ii,2)
     &             + shp(1,3)*xu(ii,3) + shp(1,4)*xu(ii,4)
          dx(ii,2) = shp(2,1)*xu(ii,1) + shp(2,2)*xu(ii,2)
     &             + shp(2,3)*xu(ii,3) + shp(2,4)*xu(ii,4)
        end do ! ii

c       Compute nodal loads for pressures

        pn = bl(3)*sg(3,l)
        do ii = 1,4
          pp     = shp(3,ii)*pn
          p(1,ii) = p(1,ii) + pp*(dx(2,1)*dx(3,2) - dx(3,1)*dx(2,2))
     &                      + bg(1)*shp(3,ii)*xjw(l)
          p(2,ii) = p(2,ii) + pp*(dx(3,1)*dx(1,2) - dx(1,1)*dx(3,2))
     &                      + bg(2)*shp(3,ii)*xjw(l)
          p(3,ii) = p(3,ii) + pp*(dx(1,1)*dx(2,2) - dx(2,1)*dx(1,2))
     &                      + bg(3)*shp(3,ii)*xjw(l)
        end do ! ii

c       Compute a tangent if necessary

        if(d(68).gt.0.0d0) then
          i1 = 0
          do ii = 1,4
            pp = shp(3,ii)*pn*ctan(1)
            j1 = 0
            do jj = 1,4
              s(i1+1,j1+2) = s(i1+1,j1+2)
     &                     - pp*(shp(1,jj)*dx(3,2) - dx(3,1)*shp(2,jj))
              s(i1+2,j1+3) = s(i1+2,j1+3)
     &                     - pp*(shp(1,jj)*dx(1,2) - dx(1,1)*shp(2,jj))
              s(i1+3,j1+1) = s(i1+3,j1+1)
     &                     - pp*(shp(1,jj)*dx(2,2) - dx(2,1)*shp(2,jj))
              s(i1+1,j1+3) = s(i1+1,j1+3)
     &                     - pp*(shp(2,jj)*dx(2,1) - dx(2,2)*shp(1,jj))
              s(i1+2,j1+1) = s(i1+2,j1+1)
     &                     - pp*(shp(2,jj)*dx(3,1) - dx(3,2)*shp(1,jj))
              s(i1+3,j1+2) = s(i1+3,j1+2)
     &                     - pp*(shp(2,jj)*dx(1,1) - dx(1,2)*shp(1,jj))
              j1 = j1 + ndf
            end do ! jj
            i1 = i1 + ndf
          end do ! ii
        endif

      end do ! l

      end

      
      subroutine sh3ftran ( dt   , shp  , aa   ,
     &                      dr   , r1   , vrn  , vr1  , arn  , ar1,
     &                      dira , dir0 , rhoa , rhoi , p    , s  )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FTRAN is the subroutine which computes the
c                        transient (inertial) terms of the residual
c                        vector and tangent stiffness matrix
c                        for the Energy-Momentum method

c        Authors:        N.Tarnow & J.C. Simo

c        Date:           March 1993.

c        Version:        This routine was tested in FEAP version 6.01

c        Caution:        Must use the conservative time-stepping
c                        algorithms in conjunction with this routine,
c                        i.e., beta,cons, etc... (nop=5).
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        dt ............ Time step.
c        shp ........... Nodal shape function values.
c        aa ............ Gauss point acceleration at T-n+a.
c        dr ............ Localized nodal incremental rotation vector.
c        r1 ............ Localized nodal total rotation vector and T-n+1.
c        vrn,vr1 ....... Gauss point rot. velocities at T-n and T-n+1.
c        arn,ar1 ....... Gauss point rot accelerations at T-n and T-n+1.
c-----[--.----+----.----+----.-----------------------------------------]
c       Remark: In case of shell intersections dr,r1, vrn, vr1, arn and
c               ar1 are used to store d/dt Lambda at tn and tn+1
c-----[--.----+----.----+----.-----------------------------------------]
c        dira .......... Localized nodal director
c                        vectors, dir(3,ElementNode) at time T-n+a.
c        dir0 .......... Localized nodal director
c                        vectors, dir(3,ElementNode) at time T-0.
c        rhoa .......... Trans. thickness weighted density multiplyied
c                        by the jacobian and weight at the Gauss point.
c        rhoi .......... Rot. inertia weighted density multiplyied
c                        by the jacobian and weight at the Gauss point.

c        Routine Output:
c        ---------------
c        s ............. Element stiffness.
c        p ............. Element residual.
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      integer         ndof
      logical                 fsetr,frotas
      common /crotas/ ndof(9),fsetr,frotas

      integer   j,   k, l, j1, k1    , m
      real*8    dt    , dtr  , dtr4
      real*8    rhoa  , rhoi , fac   , fac1 , fac2

      real*8    w1(3,3)      , w0(3,3)      , shp  (4)     , aa   (3)
      real*8    dr   (3,9)   , r1   (3,9)   , vrn  (3,9)   , vr1  (3,9)
      real*8    arn  (3,9)   , ar1  (3,9)   , wd   (3)
      real*8    dira (3,9)   , dir0(3,9)
      real*8    s    (24,24) , p    (24)

      save

c     Compute director velocities at shell intersections

      dtr  = 1.d0/dt
      dtr4 = 4.d0*dtr*dtr

      do j=1,4
        if(ndof(j).eq.6)then
          do k=1,3
            w1(k,1) = dr(k,j)
            w1(k,2) = vr1(k,j)
            w1(k,3) = ar1(k,j)
            w0(k,1) = r1(k,j)
            w0(k,2) = vrn(k,j)
            w0(k,3) = arn(k,j)
          end do ! k
          do k=1,3
            vr1(k,j) = 0.0d0
            vrn(k,j) = 0.0d0
            do l=1,3
               vr1(k,j) = vr1(k,j) + w1(k,l)*dir0(l,j)
               vrn(k,j) = vrn(k,j) + w0(k,l)*dir0(l,j)
            end do ! l
          end do ! k
        endif
      end do ! j

c     Node Loop 1

      j1 = 0
      do j = 1 , 4

c       Multiply Delta.w * t-n+a

        do k  = 1 , 3
          wd(k) = (vr1(k,j)-vrn(k,j))*dtr
        end do ! k

c       Set Up Residual

        do k = 1 , 3
          p(j1+k  ) = p(j1+k  ) - shp(j) * rhoa * aa(k)
          p(j1+k+3) = p(j1+k+3) - shp(j) * rhoi * wd(k)
        end do ! k

        j1 = j1 + 6

c     End Node Loop 1

      end do ! j

c     Node Loop 2

      j1 = 0
      do j = 1 , 4
        fac = rhoi*shp(j)*dtr
        do l = 1 , 3
          wd(l) = fac*(vr1(l,j) - vrn(l,j))
        end do ! l
        fac1 = wd(1)*dira(1,j) + wd(2)*dira(2,j) + wd(3)*dira(3,j)

c       Final Displacement Stiffness

        fac = rhoa*shp(j) * dtr4
        k1 = 0
        do k = 1 , 4
          fac2 = fac*shp(k) - fac1
          do l = 1 , 3
            s(j1+l,k1+l) = s(j1+l,k1+l) + fac2
          end do ! l
          k1 = k1 + 6
        end do ! k

c       Final Rotational Stiffness

        fac = rhoi*shp(j) * dtr4
        do m = 4 , 6
          s(j1+m,j1+m) = s(j1+m,j1+m) + fac
        end do ! m

        j1 = j1 + 6

c     End Node Loop 2

      end do ! j

      end
