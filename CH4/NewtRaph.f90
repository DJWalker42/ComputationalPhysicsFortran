!*******************************************************

! Program to compute the roots of a function using the

! Newton-Raphson method

!*******************************************************

PROGRAM NRROOTS

IMPLICIT NONE

DOUBLE PRECISION X, F, DF

EXTERNAL F, DF

X = 0.8D0

CALL NEWTRAPH( X, F, DF )

PRINT *, 'Root of f(x) = cos(x)-x is ', X

END PROGRAM NRROOTS

!*******************************************************

DOUBLE PRECISION FUNCTION F(X)

DOUBLE PRECISION X

F = COS(X) - X

END FUNCTION F

!*******************************************************

DOUBLE PRECISION FUNCTION DF(X)

DOUBLE PRECISION X

DF = -SIN(X) - 1.D0

END FUNCTION DF

!*******************************************************

SUBROUTINE NEWTRAPH( X, F, DF )

DOUBLE PRECISION X, F, DF, DELTA, ERROR

DOUBLE PRECISION, PARAMETER :: TOL = 1.D-08

INTEGER I

DO WHILE (ERROR.GT.TOL)

I = I + 1

DELTA = -F(X)/DF(X)

X = X + DELTA

ERROR = ABS( DELTA )

END DO

RETURN

END SUBROUTINE

! END OF FILE ******************************************