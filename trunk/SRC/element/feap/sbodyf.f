c$Id: sbodyf.f,v 1.1 2004/07/09 15:49:36 rlt Exp $
      subroutine sbodyf(d, body)

c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Compute body force values

c     Inputs:
c       d(*)    - Material parameters

c     Outputs:
c       body(*) - Body force intensities
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      real*8          dm
      integer            n,ma,mct,iel,nel
      common /eldata/ dm,n,ma,mct,iel,nel
      
      real*8         prldv    ,propo
      common /prld1/ prldv(50),propo


      integer    ii
      real*8     d(*), body(3)

c     Set body load levels

      do ii = 1,3
        if(int(d(73+ii)).gt.0) then
          body(ii) = d(10+ii) + prldv(int(d(73+ii)))*d(70+ii)
        else
          body(ii) = d(10+ii)*dm
        endif
      end do ! ii

      end
