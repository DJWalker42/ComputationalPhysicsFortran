MODULE EULER2_MODULE

IMPLICIT NONE

CONTAINS

!***************************************************

SUBROUTINE EULER2( DF, X, Y, A, B, N, FNAME)

DOUBLE PRECISION DF, X, Y(2), A, B, K0(2), H

INTEGER I, J, N

EXTERNAL DF

CHARACTER(*) FNAME

IF(A.GE.B)STOP 'User error: B must be .gt. A'

OPEN(UNIT=1, FILE=FNAME)

X = A

H = (B - A)/N

!Write the initial values to file

WRITE(1,100) X, Y(1), Y(2)

DO I = 1,N-1

CALL DERIV( X, Y, K0)

DO J = 1,2

Y(J) = Y(J) + H * K0(J)

END DO

X = X + H

!Write each integration step values to file

WRITE(1,100) X, Y(1), Y(2)

END DO

100 FORMAT(3(2X,F10.6))

CLOSE(1)

END SUBROUTINE EULER2

!***************************************************

SUBROUTINE DERIV( X, Y, K )

DOUBLE PRECISION X, Y(2), K(2), DF

EXTERNAL DF

K(1) = Y(2)

K(2) = DF( X, Y )

RETURN

END SUBROUTINE DERIV

END MODULE

!END OF FILE******************************************