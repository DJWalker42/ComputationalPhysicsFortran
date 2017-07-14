PROGRAM PI_POLYGON

IMPLICIT NONE

INTEGER, PARAMETER ::  N = 4

INTEGER, PARAMETER :: NRHS = 1

INTEGER  I, J

INTEGER  INFO, IPIV(N)

DOUBLE PRECISION A(N,N), P

DOUBLE PRECISION B(N)

EXTERNAL DGESV

DO J=1,N

DO I=1,N

P = (I+2)*(J-1)

A(I,J) = 1/(2**P)

END DO

END DO

PRINT *, 'Coefficient matrix A:'

!

DO I=1,N

WRITE(*,100) A(I,1), A(I,2), A(I,3), A(I,4)

100  FORMAT(4('  ', f10.6))

END DO

!

B(1) = 3.061467

B(2) = 3.121445

B(3) = 3.136548

B(4) = 3.140331

!

CALL DGESV(N, NRHS, A, N, IPIV, B, N, INFO)

PRINT *, 'Display the solution x'

DO I = 1,N

WRITE(*,101) B(I)

101 FORMAT(1('  ', f10.6))

END DO

STOP

END PROGRAM PI_POLYGON

!END OF FILE*************************************************