!!$ Program to test the performance of array addition when

!!$ row and column indices are swapped.

PROGRAM MATRIXADD

USE OMP_LIB

IMPLICIT NONE

INTEGER,          PARAMETER    :: N = 1000

DOUBLE PRECISION, ALLOCATABLE  :: A(:,:),B(:,:),C(:,:)

DOUBLE PRECISION T1,T2,T3

INTEGER I,J

ALLOCATE( A(N,N), B(N,N), C(N,N) )

A = 1.D0

B = 2.D0

T1 = OMP_GET_WTIME()

DO J = 1,N

DO I = 1,N

C(I,J) = A(I,J) + B(I,J)

END DO

END DO

T2 = OMP_GET_WTIME()

DO I = 1,N

DO J = 1,N

C(I,J) = A(I,J) + B(I,J)

END DO

END DO

T3 = OMP_GET_WTIME()

PRINT *, 'Time for JI loop:',T2-T1

PRINT *, 'Time for IJ loop:',T3-T2

DEALLOCATE( A, B, C )

END PROGRAM MATRIXADD

!END OF FILE*******************************************************