PROGRAM MONTE_DECAY

IMPLICIT NONE

!------------------------------------------------------------------------

INTEGER,          PARAMETER :: MAX = 1000, T_MAX = 500

DOUBLE PRECISION, PARAMETER :: LAMBDA = 0.01

INTEGER            :: ATOM, T, N, NCOUNT

DOUBLE PRECISION   :: DECAY

OPEN(UNIT=1, FILE='decay_data.txt')

10 FORMAT(2(I5))

! Initialise variables

N = MAX

NCOUNT = MAX

CALL INIT_RANDOM_SEED()

DO T = 1,T_MAX

DO ATOM = 1,N

CALL RANDOM_NUMBER(DECAY)

IF(DECAY.LT.LAMBDA) NCOUNT = NCOUNT - 1

END DO

N = NCOUNT

WRITE(1,10) T, N

END DO

CLOSE(1)

STOP

END PROGRAM MONTE_DECAY

!END OF FILE *******************************************************************