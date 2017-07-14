!**********************************************************

! Program to test the performance of the naive approach to

! matrix multiplication.

!**********************************************************

PROGRAM MATRIXMUL1

USE OMP_LIB

IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), C(:,:)

DOUBLE PRECISION T1, T2

INTEGER N, I, J, K

INTEGER, PARAMETER :: N_MAX = 200

OPEN(UNIT=1, FILE='matrix1_times.txt', ACTION='WRITE')

100 FORMAT(I4,2X,ES12.6)

DO N = N_MAX,N_MAX

ALLOCATE ( A(N,N), B(N,N), C(N,N) )

DO J = 1,N

DO I = 1,N

A(I,J) = I+J

B(I,J) = J

C(I,J) = 0.0

END DO

END DO

T1 = OMP_GET_WTIME()

CALL MATMUL(A,B,C,N)

T2 = OMP_GET_WTIME()

WRITE(1,100), N, T2-T1

PRINT *, T2-T1

DEALLOCATE( A, B, C )

PRINT *, '% done:',N*100/N_MAX

END DO

CLOSE(1)

END PROGRAM MATRIXMUL1

!*****************************************************

SUBROUTINE MATMUL(A,B,C,N)

INTEGER N, I, J, K

DOUBLE PRECISION A(N,N), B(N,N), C(N,N)

DO I = 1,N

DO J = 1,N

DO K = 1,N

C(I,J) = C(I,J) + A(I,K)*B(K,J)

END DO

END DO

END DO

RETURN

END SUBROUTINE MATMUL

!END OF FILE *****************************************