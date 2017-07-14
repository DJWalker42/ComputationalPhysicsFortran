MODULE RK4_MODULE

IMPLICIT NONE

CONTAINS

!*****************************************************

DOUBLE PRECISION FUNCTION DF(X,Y)

DOUBLE PRECISION X, Y

DF = -X*Y

END FUNCTION DF

!*****************************************************

DOUBLE PRECISION FUNCTION SOL(X)

DOUBLE PRECISION X

SOL = EXP(-X*X/2)

END FUNCTION SOL

!*****************************************************

SUBROUTINE RK4( DF, X, Y, A, B, N, FNAME )

DOUBLE PRECISION X, Y, A, B

DOUBLE PRECISION DF, H

DOUBLE PRECISION K0, K1, K2, K3

INTEGER I, N

!Open file to save data

CHARACTER(*) FNAME

OPEN(UNIT=1, FILE=FNAME)

IF(A.GE.B)STOP 'User error: B must be .gt. A'

!Initialise X and compute step size H

X = A

H = (B - A)/N

!Write the initial values to file

WRITE(1,100) X, Y, SOL(X)

!Perform the RK4 integration

DO I = 1, N-1

!Compute Ks

K0 = DF( X, Y )

K1 = DF( X + H/2, Y + K0*H/2 )

K2 = DF( X + H/2, Y + K1*H/2 )

K3 = DF( X + H, Y + K2*H )

!Advance one step

Y = Y + (H/6)*(K0 + 2*K1 + 2*K2 + K3)

X = X + H

!Write each integration step values to file

WRITE(1,100) X, Y, SOL(X)

END DO

100 FORMAT(3(2X,F10.6))

CLOSE(1)

END SUBROUTINE RK4

!****************************************************

SUBROUTINE RK4A( DF, X, Y, A, B, TOL, FNAME )

DOUBLE PRECISION X, Y, A, B, TOL

DOUBLE PRECISION DF, H, Y0, Y1, X1, YT, XT

DOUBLE PRECISION K0, K1, K2, K3, ERR

INTEGER I

CHARACTER(*) FNAME

OPEN(UNIT=1, FILE=FNAME)

IF(A.GE.B)STOP 'User error: B must be .gt. A'

!First step size equal to whole interval

X = A

H = B - A

!Save initial values to file

WRITE(1,101) X, Y, SOL(X)

!Set the condition to exit the while loop

DO WHILE ( (X.LT.B) )

!Store values at start of step

XT = X

YT = Y

!Compute Ks for y0

K0 = DF( X, Y )

K1 = DF( X + H/2, Y + K0*H/2 )

K2 = DF( X + H/2, Y + K1*H/2 )

K3 = DF( X + H, Y + K2*H )

!Compute y0

Y0 = Y + (H/6)*(K0 + 2*K1 + 2*K2 + K3)

!Halve the step size

H = H/2

!Perform same integration using h/2

DO I = 1,2

K0 = DF( X, Y )

K1 = DF( X + H/2, Y + K0*H/2 )

K2 = DF( X + H/2, Y + K1*H/2 )

K3 = DF( X + H, Y + K2*H )

Y1 = Y + (H/6)*(K0 + 2*K1 + 2*K2 + K3)

X1 = X + H

Y = Y1

X = X1

END DO

!Estimate the relative error in Y1

ERR = ABS(Y1 - Y0)/Y1

IF(ERR.GT.TOL)THEN

!Go back to start of the step with h/2

Y = YT

X = XT

ELSE

!Accept the step using extrapolation

Y = (16*Y - Y0)/15

!Double the starting step size (currently h/2)

H = H*4

!X has already been advanced.

WRITE(1,101) X, Y, SOL(X)

END IF

END DO

101 FORMAT(3(2X,F10.6))

CLOSE(1)

END SUBROUTINE RK4A

!****************************************************

END MODULE

!END OF FILE***************************************