PROGRAM RICHARDSON_EXTRAP

IMPLICIT NONE

DOUBLE PRECISION, PARAMETER :: TOL = 1.D-4, X = 1.D0

INTEGER,          PARAMETER :: ROW_MAX = 4

DOUBLE PRECISION FUNC, Y, H, A(ROW_MAX, ROW_MAX)

INTEGER M, L, I, J

LOGICAL SOLUTION

EXTERNAL FUNC

H = 4.D-1

A = 0.D0

SOLUTION = .FALSE.

CALL CENTRALDIFF(X, Y, H, FUNC)

A(1,1) = Y

DO M = 1, ROW_MAX - 1

H = H/2

CALL CENTRALDIFF(X, Y, H, FUNC)

A(M + 1, 1) = Y

DO L = 1,M

A(M + 1, L + 1) = (2**(2*L) * A(M + 1, L) - A(M, L))/(2**(2*L) - 1)

END DO

IF(ABS(A( M+1, M+1 ) - A( M, M )) .LT. TOL) THEN

PRINT *, 'Solution found within tolerance'

DO I = 1,M+1

WRITE (*,"(10F14.10)") ( A(I,J), J = 1, M + 1)

SOLUTION = .TRUE.

END DO

EXIT

END IF

END DO

IF(SOLUTION .EQV. .FALSE.) THEN

PRINT *, "Extrapolation not found within tolerance specified"

DO I = 1,M+1

WRITE (*,"(10F14.10)") ( A(I,J), J = 1, M + 1)

END DO

END IF

END PROGRAM RICHARDSON_EXTRAP

!******************************************************************

SUBROUTINE CENTRALDIFF (X, Y, H, FUNC)

DOUBLE PRECISION X, Y, H, FUNC

EXTERNAL FUNC

Y = (FUNC( X + H ) - FUNC ( X - H ))/2/H

END SUBROUTINE CENTRALDIFF

!*****************************************************************

DOUBLE PRECISION FUNCTION FUNC(X)

DOUBLE PRECISION X

FUNC = SIN(X)

END FUNCTION FUNC

!END OF FILE ***************************************************