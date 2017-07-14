!!$ Program to test the performance loop unrolling

PROGRAM LOOPROLL

USE OMP_LIB

IMPLICIT NONE

INTEGER,          PARAMETER       :: N = 250, MMAX = 1000

DOUBLE PRECISION, DIMENSION (N,N) :: A,B,C

DOUBLE PRECISION T1,T2,T3,T4,TA,TB,TC

INTEGER I,J,M

A = 1.D0

B = 2.D0

M = 1

TA = 0.D0

TB = 0.D0

TC = 0.D0

DO WHILE (M .LE. MMAX)

T1 = OMP_GET_WTIME()

DO J = 1,N

DO I = 1,N

C(I,J) = A(I,J) + B(I,J)

END DO

END DO

T2 = OMP_GET_WTIME()

DO J = 1,N

DO I = 1,N,2

C(I,J)   = A(I,J)   + B(I,J)

C(I+1,J) = A(I+1,J) + B(I+1,J)

END DO

END DO

T3 = OMP_GET_WTIME()

DO J = 1,N

DO I = 1,N,4

C(I,J)   = A(I,J)   + B(I,J)

C(I+1,J) = A(I+1,J) + B(I+1,J)

C(I+2,J) = A(I+2,J) + B(I+2,J)

C(I+3,J) = A(I+3,J) + B(I+3,J)

END DO

END DO

T4 = OMP_GET_WTIME()

TA = TA + T2 - T1

TB = TB + T3 - T2

TC = TC + T4 - T3

M = M + 1

END DO

PRINT *, 'Time loop 1:',TA

PRINT *, 'Time loop 2:',TB

PRINT *, 'Time loop 3:',TC

END PROGRAM LOOPROLL

!END OF FILE*****************************************************

OMP_HelloWorld.f90

PROGRAM HELLO

USE OMP_LIB

IMPLICIT NONE

INTEGER NTHREADS, TID

!Fork a team of threads giving them their own copies of variables

!$OMP PARALLEL PRIVATE(NTHREADS, TID)

!Obtain thread number

TID = OMP_GET_THREAD_NUM()

PRINT *, 'Hello World from thread = ', TID

!Only master thread does this

IF (TID .EQ. 0) THEN

NTHREADS = OMP_GET_NUM_THREADS()

PRINT *, 'Number of threads = ', NTHREADS

END IF

!All threads join master thread and disband

!$OMP END PARALLEL

END PROGRAM HELLO

!END OF FILE**********************************************

omp_param.f90

PROGRAM OMPPARAM

USE OMP_LIB

IMPLICIT NONE

PRINT *, 'Max threads', OMP_GET_MAX_THREADS()

PRINT *, 'Num of cores', OMP_GET_NUM_PROCS()

END PROGRAM OMPPARAM

!END OF FILE*********************************************