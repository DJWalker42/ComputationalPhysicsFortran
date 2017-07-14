!***************************************************************

! Program to perform an adaptive strip quadrature using the

! trapezoidal method.

! If T0 defines the numerical value of the trapezium rule

! integration with one strip, then T1 defines it with 2 strips

! (Tp, n = 2^p). Using the fact that the trapezium rule is of

! order h^2 we can obtain an estimate of the error using the

! formula (T0 - T1)/3. By comparing this estimate to the

! overall tolerance required we can either accept the numerical

! integration value, or halve the strip width and repeat

! the process until the desired accuracy is achieved.

!***************************************************************

PROGRAM ADAPTIVEQUAD

IMPLICIT NONE

! Declare variables

DOUBLE PRECISION, PARAMETER :: EPS = 1.D-3

DOUBLE PRECISION A, B, S, F

EXTERNAL F

OPEN(UNIT=1, FILE='adaptive_data.txt', ACTION='WRITE')

! Initialise integration limits

A = 6.D0

B = 1.4D1

! Initialise integration sum

S = 0.D0

CALL ADAPTIVE(F,A,B,S,EPS,0)

PRINT *, 'Numerical integration =',S

CLOSE(1)

END PROGRAM ADAPTIVEQUAD

!**************************************************************

DOUBLE PRECISION FUNCTION F(X)

DOUBLE PRECISION X, I0, G, L0

I0 = 1

L0 = 10

G  = 1

F  = I0/(1 + 4*(X - L0)*(X - L0)/G/G)

END FUNCTION F

!**************************************************************

RECURSIVE SUBROUTINE ADAPTIVE(F,A,B,S,EPS,CT)

DOUBLE PRECISION A,B,S,EPS

DOUBLE PRECISION F,H,C,T0,T1

INTEGER CT

H = B - A ! interval width for T0

T0 = (F(B) + F(A)) * H/2

H = H/2 ! halve the interval for T1

C = (B + A)/2 ! midpoint

T1 = H/2 * (2*F(C) + F(A) + F(B))

CT = CT + 1

IF (ABS( (T1-T0)/3 ).LE.EPS .AND. CT.GT.1 )THEN ! Accept integration

S = S + T1

WRITE(1,100) A, F(A)! save left hand values for plotting

100 FORMAT(2(2X,F10.6))

ELSE ! Call subroutine with adjusted limits & eps

CALL ADAPTIVE(F,A,C,S,EPS/2,CT) ! left hand strip

CALL ADAPTIVE(F,C,B,S,EPS/2,CT) ! right hand strip

END IF

!Save upper limit and function evaluation for plotting

WRITE(1,100) B, F(B)

END SUBROUTINE ADAPTIVE

!END OF FILE***************************************************