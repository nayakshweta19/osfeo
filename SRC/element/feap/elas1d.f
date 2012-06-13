c$Id: elas1d.f,v 1.1 2004/01/11 19:11:18 rlt Exp $
      subroutine elas1d(d,ta, eps, sig,dd)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2004: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  1-D elastic model

c      Inputs:
c        d(*)  - Material property parameters
c        ta    - Thermal strain
c        eps   - Strain

c      Outputs:
c        sig   - Stress
c        dd(2) - Moduli
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    ta,eps, sig,dd(2), d(*)

c     Compute stress

      sig   = sig + d(1)*(eps - d(3)*ta)

c     Set modulus

      dd(1) = d(1)
      dd(2) = d(1)

      end
