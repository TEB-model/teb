! MC2 (Mesoscale Compressible Community) model.
! Copyright (c) 2001-2020 Environment Canada (EC-RPN COMM Group).
! Licensed under the MIT license.

      SUBROUTINE VSLOG(PA,PLOG,N)

*   COMPUTES THE LOGARITHM

#include "impnone.cdk"

      INTEGER N
      REAL PA(N), PLOG(N)

      INTEGER J

      DO J=1,N
        PLOG(J) = LOG(PA(J))
      END DO

      END
