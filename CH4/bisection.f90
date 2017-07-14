!********************************************************

! Program performs root finding by the bisection method

!********************************************************

PROGRAM BISECTION

IMPLICIT NONE

! Declarations

DOUBLE PRECISION A, B, R, F

! A -- initial interval value

! B -- final interval value

! R -- computed root

! F -- function to be evaluated

! External means F defined in separate subprogram unit

EXTERNAL F

! Initialise variables

A = 0.D0

B = 1.57D0

! Call subroutine to perform bisection

CALL BISECT( A, B, R, F )

! print results to screen or to file

PRINT *, 'Root of f(x) = cos(x) - x is ', R

END PROGRAM BISECTION

!********************************************************

DOUBLE PRECISION FUNCTION F(X)

DOUBLE PRECISION X

F = COS(X) - X

END FUNCTION F

!*********************************************************

SUBROUTINE BISECT( LEFT, RIGHT, MID, F)

! Declare variables

DOUBLE PRECISION LEFT, RIGHT, MID, F

DOUBLE PRECISION FL, FR, FM, ERROR

DOUBLE PRECISION, PARAMETER :: TOL = 1.D-8

INTEGER I

! Initialise variables

I = 0

ERROR = 1.D0

FL = F(LEFT)

FR = F(RIGHT)

! Check that interval brackets root

IF(FL*FR.GT.0)THEN ! i.e. either both +ve or both -ve

MID = -9999 ! unlikely the root will be -9999

RETURN

END IF

DO WHILE (ERROR.GT.TOL)

I = I + 1

MID = (LEFT + RIGHT)/2

FM = F(MID)

IF(FL*FM.LT.0)THEN ! Root bracketed in left interval

RIGHT = MID

FR = FM

ELSE ! Root bracketed in the right interval

LEFT = MID

FL = FM

END IF

ERROR = ABS( RIGHT-LEFT )

END DO

PRINT *, 'No of iters:', I

RETURN

END SUBROUTINE BISECT

! END OF FILE **********************************************