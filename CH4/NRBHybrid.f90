!****************************************************

! Program that implements a hybrid Newton-Raphson

! Bisection method to find the roots of a function.

! The root is known to lie within an interval [A,B]

! and the result is returned to R. If the NR step

! lies with the bounds then it is accepted, else a

! bisection step is taken.

!****************************************************

PROGRAM NRBHYBRID

IMPLICIT NONE

DOUBLE PRECISION A, B, R, F, DF

EXTERNAL F, DF

A = 0.3D0

B = 0.7D0

CALL NRBISEC( A, B, R, F, DF )

PRINT *, 'Root at xr =',R

PRINT *, 'f(xr) =',F(R)

END PROGRAM NRBHYBRID

!****************************************************

DOUBLE PRECISION FUNCTION F(X)

DOUBLE PRECISION X, X8, X6, X4, X2

DOUBLE PRECISION, PARAMETER :: C8 = 6435, C6 = -12012

DOUBLE PRECISION, PARAMETER :: C4 = 6930, C2 = -1260

DOUBLE PRECISION, PARAMETER :: C0 = 35, D = 128

X8 = X*X*X*X*X*X*X*X

X6 = X*X*X*X*X*X

X4 = X*X*X*X

X2 = X*X

F = (C8*X8 + C6*X6 + C4*X4 + C2*X2 + C0)/D

END FUNCTION F

!****************************************************

DOUBLE PRECISION FUNCTION DF(X)

DOUBLE PRECISION X, X7, X5, X3, X1

DOUBLE PRECISION, PARAMETER :: C8 = 6435, C6 = -12012

DOUBLE PRECISION, PARAMETER :: C4 = 6930, C2 = -1260

DOUBLE PRECISION, PARAMETER :: D = 128

X7 = X*X*X*X*X*X*X

X5 = X*X*X*X*X

X3 = X*X*X

X1 = X

DF = (8*C8*X7 + 6*C6*X5 + 4*C4*X3 + 2*C2*X1)/D

END FUNCTION DF

!****************************************************

SUBROUTINE NRBISEC( A, B, R, F, DF )

DOUBLE PRECISION A, B, R, F, DF

DOUBLE PRECISION FA, FB, FR, DFR, DELTA, ERROR

INTEGER I

DOUBLE PRECISION, PARAMETER :: TOL = 5.D-03

! initialise variables

FA = F(A)

FB = F(B)

ERROR = 1.D0

I = 0

IF(FA * FB .GT. 0)STOP 'Root not bracketed'

! Select a "best" starting point - the limit with the

! smallest function value is likely closest to the root

IF( ABS(FA) .LE. ABS(FB) ) THEN

R = A

FR = FA

ELSE

R = B

FR = FB

END IF

DFR = DF(R)

! Set up the while loop to exit when the tolerance met

DO WHILE( ERROR .GT. TOL )

! if iteration counter goes above 30 stop execution

I = I + 1

IF(I.GT.30)STOP 'Root not converged after 30 iters'

! Decide if using NR or Bisection for this step

IF( (DFR * (R - A) - FR)&

*(DFR * (R - B) - FR).LE.0 )THEN ! NR

PRINT *, 'NR'

DELTA = -FR/DFR

R = R + DELTA

ELSE ! bisection

PRINT *, 'BI'

DELTA = (B - A)/2.D0

R = (A + B)/2.D0

END IF

FR = F(R)

DFR = DF(R)

IF(FA * FR .LE. 0)THEN

B = R

FB = FR

ELSE

A = R

FA = FR

END IF

ERROR = ABS(DELTA/R)

END DO

PRINT*, 'No of iters:', I

RETURN

END SUBROUTINE NRBISEC

! END OF FILE ****************************************