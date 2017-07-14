MODULE MOD_MC

IMPLICIT NONE

CONTAINS

!******************************************************************************

SUBROUTINE MONTE_CARLO(A, B, N, F,SUM,SIG)

IMPLICIT NONE

INTEGER            :: I, N

DOUBLE PRECISION   :: F, A, B, R, SUM, SUM2, SIG, S, S2

EXTERNAL F

IF(N.LE.0)STOP "User error: N must be .gt. zero"

IF(A.GT.B)STOP "User error: A must be .lt. B"

SUM = 0.D0

SUM2 = 0.D0

SIG = 0.D0

CALL INIT_RANDOM_SEED()

DO I = 1,N

CALL RANDOM_NUMBER(R)

R = A + (B - A)* R

SUM = SUM + F(R)

SUM2 = SUM2 + F(R)*F(R)

S = SUM/DBLE(I)

S2 = SUM2/DBLE(I)

IF(I.NE.1)SIG = SQRT( (S2 - (S*S))/(DBLE(I)-1.D0) )

END DO

SUM = (B - A)*SUM/N

SIG = (B - A)*SIG

END SUBROUTINE MONTE_CARLO

!****************************************************************************

SUBROUTINE INIT_RANDOM_SEED()

INTEGER :: I, N

REAL    :: T

INTEGER, DIMENSION(:), ALLOCATABLE :: SEED

CALL RANDOM_SEED(SIZE = N)

ALLOCATE( SEED(N) )

T = SECNDS(0.0)

SEED = INT(T) + 37 * (/ (I - 1, I = 1, N) /)

DO I =1,N

IF (MOD(SEED(I),2).EQ.0) SEED(I) = SEED(I) -1

END DO

CALL RANDOM_SEED(PUT = SEED)

DEALLOCATE(SEED)

END SUBROUTINE INIT_RANDOM_SEED

END MODULE MOD_MC

!END OF FILE *******************************************************************
