PROGRAM PI_MONTE

IMPLICIT NONE

!------------------------------------------------------------------------

! A program to calculate pi using the dart theory. A dart is randomly

! thrown into the unit square, x [0,1], y [0,1].

! From the comparison of the areas of the square and

! the quarter circle, the chance of the dart landing within a unit circle

! centred on the origin is pi/4.

!-------------------------------------------------------------------------

!-------------------------------------------------------------------------

INTEGER,          PARAMETER :: N = 1000, C = 10

DOUBLE PRECISION, PARAMETER :: ONE = 1.0D+0

INTEGER            :: I, M

DOUBLE PRECISION   :: R(2), NORM, PIE, PI

! n - total number of darts to throw

! c - print out pie every c darts thrown

! i - iteration counter; number of darts thrown.

! m - number of darts successfully landed within quarter circle.

! r - position vector representing the coordinates of where the dart lands

!     in the xy plane

! norm - size of the position vector from the origin. norm.le.1.0 dart lands

!        within quarter circle, else lands outside.

! pie - approximation of pi.

! pi - actual pi

OPEN(UNIT=1, FILE='PIdata.txt') ! file for pi approximation

OPEN(UNIT=2, FILE='RanXY.txt') ! file to check "randomness" of throws

PI = 4.D0*ATAN(1.D0)

10 FORMAT(2F12.6)

11 FORMAT(I10, 5X, 2F11.8)

! Initialise variables

M = 0

CALL INIT_RANDOM_SEED()

DO I = 1,N

! Generate random position vector

CALL RANDOM_NUMBER(R) ! populates R with 2 random numbers.

IF(MOD(I,C).EQ.0) THEN

WRITE (2,10) R(1), R(2) ! print random numbers to file

END IF

NORM = R(1)*R(1) + R(2)*R(2) ! compute distance from origin

IF(NORM.LE.ONE) THEN ! hit the circle add one to counter

M = M + 1

END IF

IF(MOD(I,C).EQ.0) THEN ! Print out every C dart thrown.

PIE = 4.0D0*M/I ! ensure double precision

WRITE (1,11) I, PIE, PI

END IF

END DO

CLOSE(1)

CLOSE(2)

END PROGRAM PI_MONTE

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

!END OF FILE *******************************************************************