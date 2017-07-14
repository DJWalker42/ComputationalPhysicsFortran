!**************************************************************************

MODULE UNITS

DOUBLE PRECISION, PARAMETER :: H  = 6.62606957D0

DOUBLE PRECISION, PARAMETER :: ME = 9.10938215D0

DOUBLE PRECISION, PARAMETER :: EV = 1.602176565D0

DOUBLE PRECISION, PARAMETER :: V0 = 1.D1, A = 5

DOUBLE PRECISION  PI, HBAR2

END MODULE UNITS

!************************************************************

PROGRAM FINITEWELL1

USE UNITS

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3) :: AE, BE, AO, BO

DOUBLE PRECISION, DIMENSION(6) :: E, ALPHA, BETA, C, D, PSI

DOUBLE PRECISION, EXTERNAL     :: EVEN, ODD

DOUBLE PRECISION, PARAMETER :: XMIN = -1.D1, XMAX = 1.D1

DOUBLE PRECISION, PARAMETER :: XSTEP = 1.D-1

DOUBLE PRECISION X, V

INTEGER I

OPEN(UNIT=1,FILE='rootfunc.txt')

! set the root brackets from plot of f(E)

AE(1) = 1.D-1; BE(1) = 5.D-1

AE(2) = 2.5D0; BE(2) = 3.D0

AE(3) = 7.D0 ; BE(3) = 7.3D0

AO(1) = 1.D0 ; BO(1) = 1.4D0

AO(2) = 4.5D0; BO(2) = 4.8D0

AO(3) = 9.6D0; BO(3) = 9.9D0

PI = 4.D0*ATAN(1.D0)

HBAR2 = H * H * 1.D2 / 4 / PI/ PI/ EV / ME

! We set ground state, which is even, in E(1) hence

! even I = odd function, odd I = even function

DO I = 1, 3

CALL BISECANT ( AE(I), BE(I), E(2*I-1), EVEN )

CALL BISECANT ( AO(I), BO(I), E(2*I), ODD )

END DO

DO I = 1,6

ALPHA(I) = SQRT( 2 * E(I)/HBAR2 )

BETA(I)  = SQRT( 2 * (V0 - E(I))/HBAR2 )

D(I) = SQRT( BETA(I)/(1+BETA(I)*A) )

IF(MOD(I,2).EQ.0)THEN ! odd function

C(I) = -D(I) * SIN(ALPHA(I)*A) * EXP(BETA(I)*A)

ELSE ! even function

C(I) = D(I) * COS(ALPHA(I)*A) * EXP(BETA(I)*A)

END IF

END DO

X = XMIN

DO WHILE ( X .LT. XMAX )

DO I = 1, 6

IF(X .LT. -A) THEN ! REGION I

V = 1.D1

PSI(I) = C(I)* EXP( BETA(I)*X )

ELSE IF ( X .LE. A) THEN ! REGION II

V = 0.D0

IF(MOD(I,2).EQ.0) THEN ! odd function

PSI(I) = D(I)*SIN(ALPHA(I)*X)

ELSE ! even function

PSI(I) = D(I)*COS(ALPHA(I)*X)

END IF

ELSE ! REGION III

V = 1.D1

IF(MOD(I,2).EQ.0)THEN !odd function

PSI(I) = -C(I)*EXP( -BETA(I)*X )

ELSE ! even function

PSI(I) = C(I)*EXP( -BETA(I)*X )

END IF

END IF

END DO

WRITE(1,'(8F12.6)') X, V, (PSI(I)*PSI(I)+E(I), I=1,6)

X = X + XSTEP

END DO

CLOSE(1)

END PROGRAM FINITEWELL1

!************************************************************

DOUBLE PRECISION FUNCTION EVEN(E)

USE UNITS

DOUBLE PRECISION ALPHA, BETA, E

! Mass = electron mass = 1

ALPHA = SQRT( 2 * E/HBAR2 )

BETA  = SQRT( 2 * (V0 - E)/HBAR2 )

EVEN  = BETA * COS(ALPHA * A) - ALPHA * SIN(ALPHA * A)

END FUNCTION EVEN

!**********************************************************

DOUBLE PRECISION FUNCTION ODD(E)

USE UNITS

DOUBLE PRECISION ALPHA, BETA, E

! Mass = electron mass = 1

ALPHA = SQRT( 2 * E/HBAR2 )

BETA  = SQRT( 2 * (V0 - E)/HBAR2 )

ODD  = ALPHA * COS(ALPHA * A) + BETA * SIN(ALPHA * A)

END FUNCTION ODD

!****************************************************

SUBROUTINE BISECANT( A, B, R, F )

DOUBLE PRECISION A, B, R, F

DOUBLE PRECISION FA, FB, FR, DFR, DELTA, ERROR

INTEGER I

DOUBLE PRECISION, PARAMETER :: TOL = 1.D-08

! initialise variables

FA = F(A)

FB = F(B)

ERROR = 1.D0

I = 0

IF(FA * FB .GT. 0)STOP 'Root not bracketed'

! Select a "best" starting point - arbitrary

IF( ABS(FA) .LE. ABS(FB) ) THEN

R = A

FR = FA

ELSE

R = B

FR = FB

END IF

DFR = (FB - FA)/(B - A)

! Set up the while loop to exit when the tolerance met

DO WHILE( ERROR .GT. TOL )

! if iteration counter goes above 30 stop execution

I = I + 1

IF(I.GT.30)STOP 'Root not converged after 30 iters'

! Decide if using Secant or Bisection for this step

IF( (DFR * (R - A) - FR)&

*(DFR * (R - B) - FR).LE.0 )THEN ! Secant

DELTA = -FR/DFR

R = R + DELTA

ELSE ! Bisection

DELTA = (B - A)/2.D0

R = (A + B)/2.D0

END IF

! Evaluate function at the new value

FR = F(R)

! Decide if root in left or right "half" of subinterval

! and adjust limits accordingly

IF(FA * FR .LE. 0)THEN ! In left

B = R

FB = FR

ELSE ! In right

A = R

FA = FR

END IF

! Update the estimate of the derivative using the new value

DFR = (FB - FA)/(B - A)

ERROR = ABS(DELTA/R)

END DO

RETURN

END SUBROUTINE

! END OF FILE ****************************************