c$Id: hvisc.f,v 1.1 2004/01/11 19:11:18 rlt Exp $
      function hvisc(x,expx)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute integration factor for viscoelastic material
c              H_visc(x) = [ 1 - exp(-x)]/x

c     Inputs:
c       x        - Current relaxation time parameter
c       expx     - Exponential of x

c     Outputs:
c       hvisc    - Viscoelastic function
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    hvisc, x, expx

      if(x.lt.1.d-04) then
        hvisc = 1.d0 - 0.5d0*x*(1.d0 - x/3.d0*(1.d0
     &               - 0.25d0*x*(1.d0 - 0.2d0*x)))
      else
        hvisc = (1.d0 - expx)/x
      endif

      end
