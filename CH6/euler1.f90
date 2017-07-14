PROGRAM EULER1A_DRIVER

IMPLICIT NONE

DOUBLE PRECISION,  PARAMETER :: A = 0.D0

DOUBLE PRECISION,  PARAMETER :: B = 2.D0

CHARACTER(LEN=13), PARAMETER :: FHEAD = 'euler1a_data_'

CHARACTER(LEN=19) FNAME

CHARACTER(LEN=2) NUM

DOUBLE PRECISION X, Y

DOUBLE PRECISION DF, SOL, ERR

INTEGER N

EXTERNAL DF, SOL

DO N = 5,20,5

WRITE(NUM,101) N

FNAME = FHEAD//NUM//'.txt'

Y = 1.D0

CALL EULER1A( DF, X, Y, A, B, N, FNAME )

ERR = ABS(SOL(X) - Y)

PRINT*, 'H =',(B - A)/N, 'ERR =',ERR

END DO

101 FORMAT(I2.2)

END PROGRAM EULER1A_DRIVER

!*************************************************

DOUBLE PRECISION FUNCTION DF(X,Y)

DOUBLE PRECISION X, Y

DF = -X * Y

END FUNCTION

!*************************************************

DOUBLE PRECISION FUNCTION SOL(X)

DOUBLE PRECISION X

SOL = EXP(-X*X/2)

END FUNCTION

!*************************************************

SUBROUTINE EULER1A( DF, X, Y, A, B, N , FNAME )

DOUBLE PRECISION X, Y, DF, H, A, B, SOL

CHARACTER(*) FNAME

INTEGER I, N

IF(A.GT.B)STOP 'User error: B must be .gt. A'

OPEN(UNIT=1, FILE=FNAME)

X = A

H = (B - A)/N

WRITE(1,100) X, Y, SOL(X)

DO I = 1,N-1

Y = Y + H * DF(X,Y)

X = X + H

WRITE(1,100) X, Y, SOL(X)

END DO

100 FORMAT(3(2X,F10.6))

CLOSE(1)

RETURN

END SUBROUTINE EULER1A

!END OF FILE ********** *************************