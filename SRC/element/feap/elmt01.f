      subroutine elmt01(d,ul,xl,ix,tl,s,r,ndf,ndm,nst,isw)

c     double *h;            // the elements local h array
c     double *u;            // to hold u^(k-1)_(n+1) as nodes don't save this
c     double *d; 	          // data stored for each element
c     int ndf;              // number of dof at nodes  nst = nen x ndf
c     int nen;              // number of element nodes
c     int ndm;	          // dimension of mesh
c     int nh1, nh3;         // the size of nh1(nh2) and nh3 data in h
c     double *s;            // element tangent (nst x nst)
c     double *r;            // element residual (ndf x nen) == (nst)
c     double *ul;           // nodal responses (ndf x nen x 5) == (nst X 5)
c     double *xl;           // nodal coordinates (ndf x nen) == (nst)
c     double *tl;           // nodal temp (nen)
c     int    *ix;           // nodal tags (nen)

c-----[--.----+----.----+----.-----------------------------------------]
c     Two dimensional frame element

c     Control Information:

c       ndm  - Spatial dimension of problem       = 2
c       ndf  - Number degree-of-freedoms at node >= 3
c              ( 1 = u_1 ; 2 = u_2 ; 3 = theta )
c       nen  - Two node element (nel = 2)        >= 2
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

c.... common declarations
      character*4 o,head
      common /bdata/ o,head(20)
      
      integer         numnp,numel,nummat,nen,neq,ipr
      common /cdata/  numnp,numel,nummat,nen,neq,ipr

      real*8          dm
      integer            n,ma,mct,iel,nel
      common /eldata/ dm,n,ma,mct,iel,nel

      real *8         bpr, ctan
      common /eltran/ bpr(3), ctan(3)
      
      integer         nh1,nh2,nh3
      common /hdata/  nh1,nh2,nh3
      
      integer         ior,iow
      common /iofile/ ior,iow
      
      real*8  hr
      common  hr(10000)
      
c ... subroutine arguments
      integer ix(*), ndf, ndm, nst, isw
      real*8  d(*),ul(ndf,nen,*),xl(ndm,*),tl(*),s(nst,*),r(*)
      
c ... local variables
      integer  i,j
      real*8  cs, sn, L, dx, dy, k, m, A, E, rho, tran(4), eps, force
      
      save

c     Set Element type

c     CHECK ELEMENT FOR ERRORS
      if(isw.eq.2) then
c     check of mesh if desired (chec)

c     Body force computation
      elseif(isw.eq.15 .or. isw.eq.23) then

        call fbody2d(d,xl,ul, r,s, ndm,ndf,nst, isw)

c     Initialize element strains for activation

      elseif(isw.eq.17) then

        j = nh3
        do i = 1,2
          hr(j  ) = ul(1,i,1)
          hr(j+1) = ul(2,i,1)
          hr(j+2) = ul(3,i,1)
          j       = j + 3
        end do ! i

c     Initialize element strains for deactivation

      elseif(isw.eq.18) then

        do i = 0,5
          hr(nh3+i) = 0.0d0
        end do ! i

c     External nodes

      elseif(isw.eq.26) then

c     Compute element arrays

      else

c       INPUT MATERIAL PARAMETERS

        if(isw.eq.1) then
          if(ior.lt.0) write(*,2000)
          write(iow,2000)
c          call inmate(d,tdof,3*2,3)
          nh1 = 0
          nh3 = 0
c         Set plot sequence
c          call pltln2(iel)

        endif

c       Small deformaion

        call frams2e(d,ul,xl,s,r,ndf,ndm,nst,isw)
        
        if(d(18).gt.0.0d0 ) then

c         Shear deformable

          if(d(79).eq.0.0d0) then

c           2-node: cubic/quadratic enhanced interpolations)

            if(nint(d(17)).eq.3) then

c              call frams2e(d,ul,xl,s,r,ndf,ndm,nst,isw)

c           2-node: linear interpolations

            else

c              call frams2d(d,ul,xl,s,r,ndf,ndm,nst,isw)

            endif

c         Euler-Bernoulli  (2-node: cubic interpolations)

          else

c            call franf2d(d,ul,xl,s,r,ndf,ndm,nst,isw)

          endif

c       Finite deformation (2-node: linear interpolations)

        else

c         Shear deformable (2-node: linear, finite displacements)

          if(d(79).eq.0.0d0) then

c            call framf2d(d,ul,xl,s,r,ndf,ndm,nst,isw)

c         No shear case    (2-node: cubic, 2-nd order displacements)

          else

c            call franf2d(d,ul,xl,s,r,ndf,ndm,nst,isw)

          endif
        endif
      endif

c     Formats

2000  format(5x,'T w o    D i m e n s i o n a l    F r a m e'/)

4000  format(' *ERROR* Element',i7,' has unspecified node: ix =',2i7)
4001  format(' *ERROR* Element',i7,' has zero length')

      end
c-----[--.----+----.----+----.-----------------------------------------]
      
      
c-----[--.----+----.----+----.-----------------------------------------]
      subroutine frams2e(d,ul,xl,s,r,ndf,ndm,nst,isw)

c     double *h;            // the elements local h array
c     double *u;            // to hold u^(k-1)_(n+1) as nodes don't save this
c     double *d; 	          // data stored for each element
c     int ndf;              // number of dof at nodes  nst = nen x ndf
c     int nen;              // number of element nodes
c     int ndm;	          // dimension of mesh
c     int nh1, nh3;         // the size of nh1(nh2) and nh3 data in h
c     double *s;            // element tangent (nst x nst)
c     double *r;            // element residual (ndf x nen) == (nst)
c     double *ul;           // nodal responses (ndf x nen x 5) == (nst X 5)
c     double *xl;           // nodal coordinates (ndf x nen) == (nst)
c     double *tl;           // nodal temp (nen)
c     int    *ix;           // nodal tags (nen)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Two dimensional small displacement frame element
c              Enhanced formulation
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

c.... common declarations
      character*4 o,head
      common /bdata/ o,head(20)
      
      real*8           defa     ,strs     ,sl
      common /bm2com/ defa(3,2),strs(3,2),sl(2,6)
      
      real*8         siglr    ,epslr
      integer                            nout
      common /bm2str/ siglr(50),epslr(50),nout
      
      integer         numnp,numel,nummat,nen,neq,ipr
      common /cdata/  numnp,numel,nummat,nen,neq,ipr
      
      integer        ndd,nie,nud
      common /cdat1/ ndd,nie,nud
      
      real*8          dm
      integer            n,ma,mct,iel,nel
      common /eldata/ dm,n,ma,mct,iel,nel
      
      real *8         bpr, ctan
      common /eltran/ bpr(3), ctan(3)
      
      integer         imtyp,mf,mq
      common /evdata/ imtyp,mf,mq
      
      integer         nh1,nh2,nh3
      common /hdata/  nh1,nh2,nh3
      
      integer         ior,iow
      common /iofile/ ior,iow
      
      integer         npart,ndfp    ,ndfo    ,ndfg
      common /part0/  npart,ndfp(20),ndfo(50),ndfg(20)
      
      integer         nph,ner
      real*8                  erav,j_int   ,jshft
      common /prstrs/ nph,ner,erav,j_int(3),jshft
      
      real*8          epl
      integer                  iepl,       neplts
      common /ptdat6/ epl(200),iepl(2,200),neplts
      
      real*8          tol,rnmax,shift
      logical                         linear,shflg
      common /rdata/  tol,rnmax,shift,linear,shflg
      
      real*8          ttim,dt,c1,c2,c3,c4,c5, chi
      common /tdata/  ttim,dt,c1,c2,c3,c4,c5, chi
      
      real*8  hr
      common  hr(10000)

      logical   enflg
      integer   ndf,ndm,nst,isw, nlay,nlob
      integer   i,ii,istrt, j,jj, ll,lint, nn,nhv,niter
      integer   is(12)
      real*8    cs,sn,a1,a2,a3,le,lr,dx,dva,dvi,cfac,dfac,xjac
      real*8    energy, rhoa, rhoi, ctan1, ctan3, shpen,two3,hdet, itol
      real*8    sg(2,4),xx(2),dma(3), shp(2,3)
      real*8    xl(ndm,*),ul(ndf,nen,*)
      real*8    d(*),r(ndf,*),s(nst,nst),pp(3,3),alp(2),dalp(2)
      real*8    eps0(3,2),eps(3,2),bc(2,3),bb(2),hh(2,2),gg(6,2),ghi(2)
      real*8    cc(3,3,2),aa(3,3),ae(2,3),sig(3,2),forc(3),mass(6,6)

      save

      data      two3 / 0.6666666666666667d0 /
      data      itol / 1.d-8 /

c     Compute element length and direction cosines

      if(isw.ge.2) then
        cs = xl(1,nel) - xl(1,1)
        sn = xl(2,nel) - xl(2,1)
        le = sqrt(cs*cs + sn*sn)
        lr = 1.d0/le
        cs = cs*lr
        sn = sn*lr

        nlay = nint(d(101))
        if(nlay.eq.0) then
          nout = 0
          nhv  = 0
          rhoa = d(4)*d(32)
          rhoi = d(4)*d(33)
        else
          nlob = nint(d(102))
          nout = (nlay - 1)*(nlob - 1) + 1
          nhv  = nint(d(15))*nout
          call int1dl(nlob,sl)
          call bm2rho(nlay,d,d(103), rhoa,rhoi)
        endif
      endif

c     Read data
      if(isw.eq.0) then
c     output element type
         if(iow.lt.0) then
            write(*,*) '   Elmt01: 2d Frame Element.'
         else
            write(iow,*) '   Elmt01: 2d Frame Element.'
         endif

      else if(isw.eq.1) then

c       Increment history storage if necessary

        nh1 = nh1 + 2  ! Add for enhanced parameters

        do i = 1,3
          is(i  ) = i
          is(i+3) = i + ndf
        end do ! i

c     Compute mass or geometric stiffness array

      elseif(isw.eq.5) then

c       Mass matrix

        if(imtyp.eq.1) then

c         Lumped mass factors

          dma(1) = 0.5d0*rhoa*le
          dma(2) = dma(1)
          dma(3) = 0.5d0*rhoi*d(8)*le

          do i = 1,3
            r(i,1) = dma(i)
            r(i,2) = dma(i)
          end do !i

c         Set mass factors         

          cfac      = d(7)
          dfac      = 1.d0 - d(7)
          mass(1,1) =  rhoa*le*(dfac*0.5d0 + cfac/3.d0)
          mass(1,4) =  mass(1,1)*0.5d0*cfac
          mass(2,2) =  mass(1,1)
          mass(2,3) =  rhoa*le*le/24.d0*cfac
          mass(2,5) =  mass(1,4)
          mass(2,6) = -mass(2,3)
          mass(3,2) =  mass(2,3)
          mass(3,3) = (rhoa*le*le*le/120.d0 + rhoi*le/3.0d0)*cfac
          mass(3,5) =  mass(2,3)
          mass(3,6) = -mass(3,3)            - rhoi*le/6.0d0*cfac
          mass(3,3) =  mass(3,3)            + rhoi*le*0.5d0*dfac
          mass(4,1) =  mass(1,4)
          mass(4,4) =  mass(1,1)
          mass(5,2) =  mass(2,5)
          mass(5,3) =  mass(3,5)
          mass(5,5) =  mass(2,2)
          mass(5,6) = -mass(2,3)
          mass(6,2) =  mass(2,6)
          mass(6,3) =  mass(3,6)
          mass(6,5) =  mass(5,6)
          mass(6,6) =  mass(3,3)

          do j = 1,6
            do i = 1,6
              s(is(i),is(j)) = s(is(i),is(j)) + mass(i,j)
            end do ! i
          end do ! j

c       Geometric stiffness

        elseif(imtyp.eq.2) then

        endif

c      elseif(isw.eq.12) then

c        lint = nel
c        call int1dg(lint,sg)
c        call bm2trn (ul(1,1,1),cs,sn,ndf*nel,ndf,2)
c        call bm2trn (ul(1,1,2),cs,sn,ndf*nel,ndf,2)
c        call bm2trn (xl       ,cs,sn,ndm*nel,ndm,2)
c        nn = 2
c        do ll = 1,lint
c          call b2mode(d,sig,eps,cc,energy,h(nh1+nn),h(nh2+nn),
c     &                  ndf,nen,isw,ll)
c          nn = nn + nhv
c        end do ! ll

c     Compute energy

      elseif(isw.eq.13) then

        lint = nel - 1
        call int1dg(lint,sg)
        call bm2trn (ul(1,1,1),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,2),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,4),cs,sn,ndf*nel,ndf,2)
        call bm2trn (xl       ,cs,sn,ndm*nel,ndm,2)
        dva = 0.5d0*rhoa*le
        dvi  =0.5d0*rhoi*d(8)*le

c       Compute internal energy

        nn = 2
        do ll = 1,lint

c         Compute energy density from stress and deformation

          call shp1d(sg(1,ll),xl,shp,ndm,nel,xjac)
          dx = sg(2,ll)*xjac
          call b2mode(d,sig,eps,cc,energy,hr(nh1+nn),hr(nh2+nn),
     &                  ndf,nen,isw,ll)

c         Accumulate energy

          epl(8) = epl(8) + 0.5d0*energy*dx

          nn = nn + nhv
        end do ! ll

c       Compute kinetic energy for lumped mass

        epl(7) = epl(7) + 0.5d0*dva*(ul(1,1,4)**2 + ul(1,2,4)**2
     &                             + ul(2,1,4)**2 + ul(2,2,4)**2)
     &                  + 0.5d0*dvi*(ul(3,1,4)**2 + ul(3,2,4)**2)

c     Initialize history variables

      elseif(isw.eq.14) then

        istrt  = nint(d(84))
        call modl1d(d,le,xx,hr(nh1),hr(nh2),nint(d(15)),1,istrt,
     &              forc,aa,isw)

c     Residual, tangent, stress, projection options

      elseif(isw.eq.3 .or. isw.eq.4 .or. isw.eq.6 .or. isw.eq.8) then

c       Transform nodal parameters to local frame

        call bm2trn (ul(1,1,1),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,2),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,5),cs,sn,ndf*nel,ndf,2)
        call bm2trn (xl       ,cs,sn,ndm*nel,ndm,2)

c       Extract enhanced parameters

        alp(1) = hr(nh2)
        alp(2) = hr(nh2+1)

c       Set quadrature

        lint   = nel
        call int1dg(lint,sg)

c       Compute enhanced terms

        eps0(1,1) = (ul(1,2,1) - ul(1,1,1))*lr     ! Axial
        eps0(2,1) = (ul(2,2,1) - ul(2,1,1))*lr     ! Shear
     &            - (ul(3,2,1) + ul(3,1,1))*0.5d0
        eps0(3,1) = (ul(3,2,1) - ul(3,1,1))*lr     ! Bending

        niter = 0
        enflg = .true.
        do while (enflg)

          niter   = niter + 1

          hh(1,1) = 0.0d0
          hh(1,2) = 0.0d0
          hh(2,1) = 0.0d0
          hh(2,2) = 0.0d0
          bb(1)   = 0.0d0
          bb(2)   = 0.0d0

          do j = 1,3
            do i = 1,3
              aa(i,j) = 0.0d0
            end do ! i
            ae(1,j) = 0.0d0
            ae(2,j) = 0.0d0
            forc(j) = 0.0d0
          end do ! j
            

          nn = 2
          do ll = 1,lint
            
            dx     = sg(2,ll)*le*0.5d0
            shpen  = 4.d0*sg(1,ll)*lr
            eps(1,1) = eps0(1,1) + shpen*alp(1)
            eps(2,1) = eps0(2,1) + two3 *alp(2)
            eps(3,1) = eps0(3,1) + shpen*alp(2)

            call b2mode(d,sig,eps,cc,energy,hr(nh1+nn),hr(nh2+nn),
     &                  ndf,nen,isw,ll)

c           Multiply moduli by solution parameter: ctan(1)

c           ctan1 = ctan(1) + d(78)*ctan(2)
            ctan1 = ctan(1)*dx
            do jj = 1,3
              do ii = 1,3
                cc(ii,jj,1) = cc(ii,jj,1)*ctan1
              end do ! ii
              sig(jj,1) = sig(jj,1)*dx
            end do ! jj

c           Compute enhanced residual and tangent

            bb(1) = bb(1) - shpen*sig(1,1) 
            bb(2) = bb(2) -  two3*sig(2,1) - shpen*sig(3,1)
            do i = 1,3
              bc(1,i) = shpen*cc(1,i,1) 
              bc(2,i) = two3*cc(2,i,1) + shpen*cc(3,i,1)
            end do ! i
            hh(1,1) = hh(1,1) + bc(1,1)*shpen
            hh(1,2) = hh(1,2) + bc(1,2)*two3 + bc(1,3)*shpen
            hh(2,1) = hh(2,1) + bc(2,1)*shpen
            hh(2,2) = hh(2,2) + bc(2,2)*two3 + bc(2,3)*shpen

c           Accumulate enhanced integrals for stiffness

            do j = 1,3
              do i = 1,3
                aa(i,j) = aa(i,j) + cc(i,j,1)
              end do ! i
              ae(1,j) = ae(1,j) + bc(1,j)
              ae(2,j) = ae(2,j) + bc(2,j)
              forc(j) = forc(j) + sig(j,1)
            end do ! j

c           Increment history counter

            nn = nn + nhv

          end do ! ll

c         Invert enhanced stiffness and compute incremental alpha

          hdet    =  1.d0/(hh(1,1)*hh(2,2) - hh(1,2)*hh(1,2))
          dx      =  hh(1,1)*hdet
          hh(1,1) =  hh(2,2)*hdet
          hh(1,2) = -hh(1,2)*hdet
          hh(2,1) = -hh(2,1)*hdet
          hh(2,2) =  dx
          dalp(1) =  hh(1,1)*bb(1) + hh(1,2)*bb(2)
          dalp(2) =  hh(2,1)*bb(1) + hh(2,2)*bb(2)
          alp(1)  =  alp(1) + dalp(1)
          alp(2)  =  alp(2) + dalp(2)

c         Check convergence

          if(max(abs(dalp(1)),abs(dalp(2))) .le.
     &       max(abs( alp(1)),abs( alp(2)))*itol) then
            enflg = .false.
          elseif(niter.gt.3) then
            enflg = .false.
          endif
        end do ! while

c       Save history variables

        hr(nh2  ) = alp(1) 
        hr(nh2+1) = alp(2) 
          
c       Residual
          
        do i = 1,3
          r(i,1) = r(i,1) + forc(i)*lr
          r(i,2) = r(i,2) - forc(i)*lr
        end do ! i
        r(3,1) = r(3,1) + 0.5d0*forc(2)
        r(3,2) = r(3,2) + 0.5d0*forc(2)

c       Stiffness for main terms

        if(isw.eq.3) then
          a1 = lr*lr
          a2 = lr*0.5d0
          a3 = 0.25d0*aa(2,2)
          do j = 1,3
            do i = 1,3
              s(i    ,j    ) = s(i    ,j    ) + a1*aa(i,j)
              s(i+ndf,j    ) = s(i+ndf,j    ) - a1*aa(i,j)
              s(i    ,j+ndf) = s(i    ,j+ndf) - a1*aa(i,j)
              s(i+ndf,j+ndf) = s(i+ndf,j+ndf) + a1*aa(i,j)
            end do ! i

            s(j    ,3    ) = s(j    ,3    ) +a2*aa(j,2)
            s(j    ,3+ndf) = s(j    ,3+ndf) +a2*aa(j,2)
            s(j+ndf,3    ) = s(j+ndf,3    ) -a2*aa(j,2)
            s(j+ndf,3+ndf) = s(j+ndf,3+ndf) -a2*aa(j,2)

            s(3    ,j    ) = s(3    ,j    ) +a2*aa(2,j)
            s(3+ndf,j    ) = s(3+ndf,j    ) +a2*aa(2,j)
            s(3    ,j+ndf) = s(3    ,j+ndf) -a2*aa(2,j)
            s(3+ndf,j+ndf) = s(3+ndf,j+ndf) -a2*aa(2,j)
          end do ! j
          s(3    ,3    ) = s(3    ,3    ) + a3
          s(3    ,3+ndf) = s(3    ,3+ndf) + a3
          s(3+ndf,3    ) = s(3+ndf,3    ) + a3
          s(3+ndf,3+ndf) = s(3+ndf,3+ndf) + a3

c         Coupling and static condensation

          do i = 1,3
            gg(i  ,1) = -lr*ae(i,1)
            gg(i  ,2) = -lr*ae(i,2)
            gg(i+3,1) =  lr*ae(i,1)
            gg(i+3,2) =  lr*ae(i,2)
          end do ! i
          do i = 1,2
            gg(3,i) = gg(3,i) - 0.5d0*ae(2,i)
            gg(6,i) = gg(6,i) - 0.5d0*ae(2,i)
          end do ! i

          do i = 1,6
            ghi(1) = gg(i,1)*hh(1,1) + gg(i,2)*hh(2,1)
            ghi(2) = gg(i,1)*hh(1,2) + gg(i,2)*hh(2,2)
            do j = 1,6
              s(is(i),is(j)) = s(is(i),is(j)) - ghi(1)*gg(j,1)
     &                                        - ghi(2)*gg(j,2)
            end do ! j
          end do ! i
        endif

c       Lumped and consistent inertia contributions
       
        if((ndfo(1).gt.0 .or. shflg)  .and.
     &          isw.eq.3 .or. isw.eq.6) then
          cfac = d(7)
          dfac = 1.d0 - d(7)

          dma(1) = 0.5d0*rhoa*le
          dma(2) = dma(1)
          dma(3) = 0.5d0*rhoi*d(8)*le

c         Set mass factors         

          mass(1,1) =  rhoa*le*(dfac*0.5d0 + cfac/3.d0)
          mass(1,4) =  mass(1,1)*0.5d0*cfac
          mass(2,2) =  mass(1,1)
          mass(2,3) =  rhoa*le*le/24.d0*cfac
          mass(2,5) =  mass(1,4)
          mass(2,6) = -mass(2,3)
          mass(3,2) =  mass(2,3)
          mass(3,3) = (rhoa*le*le*le/120.d0 + rhoi*le/3.0d0)*cfac
          mass(3,5) =  mass(2,3)
          mass(3,6) = -mass(3,3)            - rhoi*le/6.0d0*cfac
          mass(3,3) =  mass(3,3)            + rhoi*le*0.5d0*dfac
          mass(4,1) =  mass(1,4)
          mass(4,4) =  mass(1,1)
          mass(5,2) =  mass(2,5)
          mass(5,3) =  mass(3,5)
          mass(5,5) =  mass(2,2)
          mass(5,6) = -mass(2,3)
          mass(6,2) =  mass(2,6)
          mass(6,3) =  mass(3,6)
          mass(6,5) =  mass(5,6)
          mass(6,6) =  mass(3,3)

          do i = 1,3
            do j = 1,3
              r(i,1) = r(i,1)- mass(i  ,j  )*(ul(i,1,5)+d(77)*ul(i,1,4))
     &                       - mass(i  ,j+3)*(ul(i,2,5)+d(77)*ul(i,2,4))
              r(i,2) = r(i,2)- mass(i+3,j  )*(ul(i,1,5)+d(77)*ul(i,1,4))
     &                       - mass(i+3,j+3)*(ul(i,2,5)+d(77)*ul(i,2,4))
            end do ! j
          end do !i

          if(isw.eq.3) then
            ctan3  = ctan(3) + d(77)*ctan(2)
            do j = 1,6
              do i = 1,6
                s(is(i),is(j)) = s(is(i),is(j)) + mass(i,j)*ctan3
              end do ! i
            end do ! j
          endif
        endif

c       Transform stiffness and residual to global coordinates

        if(isw.eq.3) then
          call bm2trn (s,cs,sn,nst,ndf,1)
        endif
        if(isw.eq.3 .or. isw.eq.6) then
          call bm2trn ( r,cs,-sn,nst,ndf,2)

c         Set body loads

          call fbody2d(d,xl,ul, r,s, ndm,ndf,nst, isw)

c       Output forces

        elseif(isw.eq.4 .or. isw.eq.8) then

          call bm2trn (xl       ,cs,-sn,ndm*nel,ndm,2)

          do i = 1,3
            pp(i,1) =  r(i,1)
            pp(i,2) = -r(i,2)
            r(i,1)  =  0.0d0
            r(i,2)  =  0.0d0
          end do ! i

          if(isw.eq.4) then
            mct = mct - 3
            if (mct.le.0) then
              write(iow,2001) o,head,ttim
              if(ior.lt.0) write(*,2001) o,head,ttim
              mct = 50
            endif
            write(iow,2002) n,ma,((xl(i,j),i=1,2),(pp(i,j),i=1,3),j=1,2)
            if(ior.lt.0) then
              write(*,2002) n,ma,((xl(i,j),i=1,2),(pp(i,j),i=1,3),j=1,2)
            endif

c         Stress projections save

          else

c            call frcn2d(pp,r,s)

          endif

        endif

      endif

c     Formats

2001  format(a1,20a4/5x,'Time',1p,e12.4,6x,' Element Forces '//
     & 43x,'*************  FORCE  *************'/
     &  2x,'Element Material',
     &  5x,'1-Coord',5x,'2-Coord',5x,'n-dir',7x,'s-dir',7x,'m-dir'/)

2002  format(2i9,0p,2e12.3,1p,3e12.3/18x,0p,2e12.3,1p,3e12.3)

      end

      subroutine b2mode(d,sig,eps,cc,engy,hn,h1,ndf,nen,isw,ll)

      implicit none

      real*8           defa     ,strs     ,sl
      common /bm2com/ defa(3,2),strs(3,2),sl(2,6)
      
      real*8          theta
      integer                  nrk,nrc,nrm,nrt,noi,nt,nrkn,nrcn,nrmn
      common /ddata/  theta(3),nrk,nrc,nrm,nrt,noi,nt,nrkn,nrcn,nrmn
      
      real*8          dm
      integer            n,ma,mct,iel,nel
      common /eldata/ dm,n,ma,mct,iel,nel
      
      real*8          ttim,dt,c1,c2,c3,c4,c5, chi
      common /tdata/  ttim,dt,c1,c2,c3,c4,c5, chi

      integer  ndf,nen, ii,isw, ll, nh
      real*8   d(*),cc(3,3,2), eps(3,2),sig(3,2), engy,hn(*),h1(*)

      save

c     Compute constitutive model type

      nh   = nint(d(15))

c     Compute forces

      call bm2con (d,hn(1),h1(1),nh,cc,sig,eps,eps,isw)

c     Save forces and deformation for outputs

      ii       = 6*(ll-1)
c      tt(ii+1) = sig(1,1)
c      tt(ii+2) = eps(1,1)
c      tt(ii+3) = sig(2,1)
c      tt(ii+4) = eps(2,1)
c      tt(ii+5) = sig(3,1)
c      tt(ii+6) = eps(3,1)

c     Compute stored energy density

      if(isw.eq.13) then

        engy = sig(1,1)*eps(1,1)+sig(2,1)*eps(2,1)+sig(3,1)*eps(3,1)

      endif

      end

      

      subroutine bm2con(d,hn,h1,nh, cc,strs,def,def1,isw)

      implicit  none

      integer         nh1,nh2,nh3
      common /hdata/  nh1,nh2,nh3
      
      real*8  h
      common  h(1000)

      integer   isw, i, nh,nlay
      real*8    d(*),hn(*),h1(*), cc(3,3,2),strs(3,2),def(3,2),def1(3)
      real*8    dl(3)

      save

c     Resultant model

      nlay = int(d(101))

      if(nlay.eq.0) then
        dl(1) = d(1)*d(32)
        dl(2) = d(1)*d(32)*d(37)*0.5d0/(1.d0+d(2))
        dl(3) = d(1)*d(33)

        do i = 1,9
          cc(i,1,1) = 0.0d0
          cc(i,1,2) = 0.0d0
        end do ! i

c       Elastic resultant model only

        do i = 1,3
          cc(i,i,1) = dl(i)
          cc(i,i,2) = dl(i)
          strs(i,1) = dl(i)*def(i,1)
          strs(i,2) = dl(i)*def(i,2)
        end do ! i

        if(nint(d(160)).eq.2) then ! Add augmented value
          strs(1,1) = strs(1,1) + h(nh2)
          if(isw.eq.10) then ! Update for augmentation
            h(nh2) = strs(1,1)
          endif
        endif

c     Thickness quadrature model

      else

        if(isw.eq.12) then
          call bm2ups(nlay,d,d(103),hn,h1,nh,def1, isw)
        else
          call bm2res(nlay,d,d(103),hn,h1,nh,def, strs,cc, isw)
        endif

      endif

      end

      subroutine bm2trn(s,cs,sn,nst,ndf,itype)

c     Itype:

c        1  Transform matrix s(nst,nst)
c        2  Transform vector s(nst,1)

      implicit  none

      integer   nst,ndf,itype, i,j,n
      real*8    cs,sn,t,tol, s(nst,*)

      save

      data      tol/ 1.d-12 /

      if(cs.gt.1.0d0-tol) return

      if(itype.eq.1) then

        do i = 1,nst,ndf
          j = i + 1
          do n = 1,nst
            t      = s(n,i)*cs - s(n,j)*sn
            s(n,j) = s(n,i)*sn + s(n,j)*cs
            s(n,i) = t
          end do ! n
        end do ! i
        do i = 1,nst,ndf
          j = i + 1
          do n = 1,nst
            t      = s(i,n)*cs - s(j,n)*sn
            s(j,n) = s(i,n)*sn + s(j,n)*cs
            s(i,n) = t
          end do ! n
        end do ! i

      else

        do i=1,nst,ndf
          j = i + 1
          t      =  s(i,1)*cs + s(j,1)*sn
          s(j,1) = -s(i,1)*sn + s(j,1)*cs
          s(i,1) =  t
        end do ! i

      endif

      end
      
      
      
      subroutine bm2res(nlay,d,bz,hn,h1,nh,def, forc,aa, isw)

c     Inputs:
c        nlay      - Number of z-levels
c        d(*)      - Material parameters
c        bz(2,*)   - Z-coordinate and B-width
c        hn(*)     - History terms at t_n
c        h1(*)     - History terms at t_n+1
c        nh        - Number of history terms per level
c        def(3,2)  - Axial, shear, and bending strains

c     Outputs:
c        forc(3,2) - Force resultants: Axial, shear, bending
c        aa(3,3,2) - Modulus array

      implicit  none

      real*8           defa     ,strs     ,sl
      common  /bm2com/ defa(3,2),strs(3,2),sl(2,6)

            real*8         siglr    ,epslr
      integer                            nout
      common/bm2str/ siglr(50),epslr(50),nout

      logical   rayl
      integer   i,ii,istrt,n,nlay,nn,nh,isw
      real*8    ta,eps(2),sig(2),dd(2),kg,ks(2), bb,zz,df,dz
      real*8    d(*),bz(2,*),hn(*),h1(*),def(3,2),forc(3,2),aa(3,3,2)

      save

      data     ta /0.0d0/

c     Set flag for Rayleigh damping

      rayl = (d(77).ne.0.0d0  .or. d(78).ne.0.0d0)

c     Set shear parameters to elastic

      kg    = d(1)*d(37)*0.5d0/(1.d0+d(2))
      ks(1) = kg*def(2,1)
      ks(2) = kg*def(2,2)

c     Initialize arrays

      do i = 1,3
        forc(i,1) = 0.0d0
        forc(i,2) = 0.0d0
        aa(i,1,1) = 0.0d0
        aa(i,2,1) = 0.0d0
        aa(i,3,1) = 0.0d0
        aa(i,1,2) = 0.0d0
        aa(i,2,2) = 0.0d0
        aa(i,3,2) = 0.0d0
      end do ! i

c     Compute constitution using Gauss-Lobbato quadrature in layers

      nn     = 1
      ii     = 1
      zz     = bz(1,1)
      eps(1) = def(1,1) - zz*def(3,1)
      eps(2) = def(1,2) - zz*def(3,2)
      istrt  = nint(d(84))
      call modl1d(d,ta,eps,hn(nn),h1(nn),nh,ii,istrt, dd,sig, isw)
      siglr(1) = sig(1)
      epslr(1) = eps(1)
c      tt(7)    = sig(1)
c      tt(8)    = eps(1)

      do n = 1,nlay-1

c       Layer halfthickness

        dz        = 0.5d0*(bz(1,n+1) - bz(1,n))

c       First point for each layer

        bb        = bz(2,n)*dz*sl(2,1)

        df        = bb*sig(1)
        forc(1,1) = forc(1,1) + df
        forc(2,1) = forc(2,1) + ks(1)*bb
        forc(3,1) = forc(3,1) - df*zz

        df        = bb*dd(1)
        aa(1,1,1) = aa(1,1,1) + df
        aa(1,3,1) = aa(1,3,1) - df*zz
        aa(3,3,1) = aa(3,3,1) + df*zz*zz
        aa(2,2,1) = aa(2,2,1) + kg*bb

        if(rayl) then
          df        = bb*sig(2)
          forc(1,2) = forc(1,2) + df
          forc(2,2) = forc(2,2) + ks(2)*bb
          forc(3,2) = forc(3,2) - df*zz

          df        = bb*dd(2)
          aa(1,1,2) = aa(1,1,2) + df
          aa(1,3,2) = aa(1,3,2) - df*zz
          aa(3,3,2) = aa(3,3,2) + df*zz*zz
          aa(2,2,2) = aa(2,2,2) + kg*bb

        endif ! rayl

c       Add remaining points

        do i = 2,int(d(102))

          zz = 0.5d0*((1.d0 - sl(1,i))*bz(1,n)
     &              + (1.d0 + sl(1,i))*bz(1,n+1))
          bb = 0.5d0*((1.d0 - sl(1,i))*bz(2,n)
     &              + (1.d0 + sl(1,i))*bz(2,n+1))*dz*sl(2,i)
          eps(1) = def(1,1) - zz*def(3,1)
          eps(2) = def(1,2) - zz*def(3,2)

c         Increment counters

          nn  = nn + nh
          ii  = ii + 1

c         Get constitution: Stress (sig) and moduli (dd)

          istrt  = nint(d(84))
          call modl1d(d,ta,eps,hn(nn),h1(nn),nh,ii,istrt, dd,sig, isw)

          siglr(ii) = sig(1)
          epslr(ii) = eps(1)
c          tt(2*ii+5)= sig(1)
c          tt(2*ii+6)= eps(1)

          df          = bb*sig(1)
          forc(1,1)   = forc(1,1) + df
          forc(2,1)   = forc(2,1) + ks(1)*bb
          forc(3,1)   = forc(3,1) - df*zz

          df          = bb*dd(1)
          aa(1,1,1)   = aa(1,1,1) + df
          aa(1,3,1)   = aa(1,3,1) - df*zz
          aa(3,3,1)   = aa(3,3,1) + df*zz*zz
          aa(2,2,1)   = aa(2,2,1) + kg*bb

          if(rayl) then
            df          = bb*sig(2)
            forc(1,2)   = forc(1,2) + df
            forc(2,2)   = forc(2,2) + ks(2)*bb
            forc(3,2)   = forc(3,2) - df*zz

            df          = bb*dd(1)
            aa(1,1,2)   = aa(1,1,2) + df
            aa(1,3,2)   = aa(1,3,2) - df*zz
            aa(3,3,2)   = aa(3,3,2) + df*zz*zz
            aa(2,2,2)   = aa(2,2,2) + kg*bb
          endif ! rayl

        end do ! i
      end do ! n

      aa(3,1,1) = aa(1,3,1)
      aa(3,1,2) = aa(1,3,2)

c      tt(1)   = forc(1,1)
c      tt(2)   = def (1,1)
c      tt(3)   = forc(2,1)
c      tt(4)   = def (2,1)
c      tt(5)   = forc(3,1)
c      tt(6)   = def (3,1)

      end

      subroutine bm2ups(nlay,d,bz,hn,h1,nh,def,isw)

c     Inputs:
c        nlay    - Number of z-levels
c        d(*)    - Material parameters
c        bz(2,*) - Z-coordinate and B-width
c        hn(*)   - History terms at t_n
c        h1(*)   - History terms at t_n+1
c        nh      - Number of history terms per level
c        def(3)  - Axial, shear, and bending strains at t_n+alpha
c        isw     - Element control parameter

      implicit  none

      real*8           defa     ,strs     ,sl
      common  /bm2com/ defa(3,2),strs(3,2),sl(2,6)

      integer   i,istrt, n,nlay,nn,nh,isw
      real*8    ta,eps,zz,d(*),bz(2,*),hn(*),h1(*),def(3)

      save

      data      ta /0.0d0/

c     Compute constitution updates using Gauss-Lobbato quadrature

      nn    = 1
      eps   = def(1) - bz(1,1)*def(3)
      istrt = nint(d(84))
      call modl1u(d,ta,eps,hn(nn),h1(nn),nh,istrt, isw)
      nn       = nn + nh
      do n = 1,nlay-1
        do i = 2,int(d(102))
          zz  = 0.5d0*((1.d0 - sl(1,i))*bz(1,n)
     &               + (1.d0 + sl(1,i))*bz(1,n+1))
          eps = def(1) - zz*def(3)
          call modl1u(d,ta,eps,hn(nn),h1(nn),nh,istrt, isw)
          nn  = nn + nh
        end do ! i
      end do ! n

      end

      subroutine bm2rho(nlay,d,bz, rhoa,rhoi)

c     Inputs:
c        nlay    - Number of z-levels
c        d(*)    - Material parameters
c        bz(2,*) - Z-coordinate and B-width

c     Outputs:
c        rhoa    - integral * rho * b(z) * dz
c        rhoi    - integral * rho * b(z) * z*z * dz

      implicit  none

      real*8           defa     ,strs     ,sl
      common  /bm2com/ defa(3,2),strs(3,2),sl(2,6)

      integer   i,n,nlay
      real*8    rhoa,rhoi, zz,df,dz, d(*),bz(2,*),brho(50)

      save

c     Compute density times width at each level

      do n = 1,nlay
        brho(n) = bz(2,n)*d(4)
      end do ! n

c     Initialize values

      rhoa = 0.0d0
      rhoi = 0.0d0

c     Integrate over each layer

      do n = 1,nlay-1
        dz = 0.5d0*(bz(1,n+1) - bz(1,n))
        if(dz .gt. 0.0d0) then
          do i = 1,int(d(102))
            zz   = 0.5d0*((1.d0 - sl(1,i))*bz(1,n)
     &                  + (1.d0 + sl(1,i))*bz(1,n+1))
            df   = 0.5d0*((1.d0 - sl(1,i))*brho(n)
     &                  + (1.d0 + sl(1,i))*brho(n+1))*dz*sl(2,i)

            rhoa =  rhoa + df
            rhoi =  rhoi + df*zz*zz
          end do ! i
        endif
      end do ! n

      end