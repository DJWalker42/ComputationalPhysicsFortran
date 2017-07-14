!**********************************************************

! Program to test the performance of the naive approach to

! matrix multiplication.

!**********************************************************

PROGRAM MATRIXMUL2

USE OMP_LIB

IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), C(:,:)

DOUBLE PRECISION T1, T2

INTEGER I, J, K, NB, M

INTEGER IMIN, IMAX, JMIN ,JMAX, KMIN, KMAX

INTEGER, PARAMETER :: N = 2**9

OPEN(UNIT=1, FILE='matrix2_times.txt', ACTION='WRITE')

100 FORMAT(I4,2X,F10.6)

ALLOCATE ( A(N,N), B(N,N), C(N,N) )

DO J = 1,N

DO I = 1,N

A(I,J) = I+J

B(I,J) = J

C(I,J) = 0.0

END DO

END DO

M = 2**3

NB = N/M

T1 = OMP_GET_WTIME()

DO K = 1,M

KMIN = (K-1)*NB + 1

KMAX = K*NB

DO I = 1,M

IMIN = (I-1)*NB + 1

IMAX = I*NB

DO J = 1,M

JMIN = (J-1)*NB + 1

JMAX = J*NB

CALL MATMUL( A(IMIN:IMAX, KMIN:KMAX), &

B(KMIN:KMAX, JMIN:JMAX), &

C(IMIN:IMAX, JMIN:JMAX), NB )

END DO

END DO

END DO

T2 = OMP_GET_WTIME()

CALL DGEMM('N', 'N', N, N, N, 1.D0, A, N, B, N, 0.D0, C, N)

T3 = OMP_GET_WTIME()

WRITE(1,100), NB, T2-T1

PRINT *,'Block time:', T2-T1

PRINT *,'LAPACK time:', T3-T2

DEALLOCATE( A, B, C )

CLOSE(1)

END PROGRAM MATRIXMUL2

!*****************************************************

SUBROUTINE MATMUL(AS,BS,CS,N)

INTEGER N, I, J, K

DOUBLE PRECISION AS(N,N), BS(N,N), CS(N,N)

DO I = 1,N

DO J = 1,N

DO K = 1,N

CS(I,J) = CS(I,J) + AS(I,K)*BS(K,J)

END DO

END DO

END DO

RETURN

END SUBROUTINE MATMUL

!END OF FILE *****************************************