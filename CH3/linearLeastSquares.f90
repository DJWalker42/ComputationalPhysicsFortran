!********************************************************************

!  Program to compute the least linear squares data fit to the data

!  stored in ‘llsdata.txt’

!********************************************************************

PROGRAM LLS

IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:), XI(:), YI(:), T(:)

INTEGER I,J,K,N,M,INFO

INTEGER, ALLOCATABLE :: IPIV(:)

! M - order of the approximating polynomial

!     i.e. y(x) ~ a0 + a1*x + ... + am*x^m

! N - total number of data point pairs

! INFO and IPIV are required for the LAPACK routine

M = 2

ALLOCATE ( A(M+1,M+1), B(M+1), IPIV(M+1) )

A = 0.0D0

B = 0.0D0

OPEN(UNIT=1, FILE="llsdata.txt")

READ(1,*) N ! N saved in first line of text file

ALLOCATE ( XI(N), YI(N), T(N) )

T = 1.0D0

I = 1

DO WHILE(1.EQ.1)

READ(1,*, END=10) XI(I), YI(I)

! x & y data saved in two separated columns

I = I + 1

END DO

10 CONTINUE

CLOSE(1)

! Error check total number of data point pairs match

IF((I-1).NE.N) THEN

PRINT *, 'Number of data points does not match header'

STOP

END IF

! Set up matrix A and RHS vector B

DO J = 1,M+1

DO K = 1,N

A(1,J) = A(1,J) + T(K)

A(2,J) = A(2,J) + T(K)*XI(K)

A(3,J) = A(3,J) + T(K)*XI(K)*XI(K)

B(J) = B(J) + YI(K)*T(K)

END DO

DO I = 1,N

T(I) = T(I)*XI(I)

END DO

END DO

! Solve the linear system of equations using LAPACK's DGESV routine

CALL DGESV(M+1, 1, A, M+1, IPIV, B, M+1,INFO)

WRITE(6,*) 'X = '

WRITE(6,101) (B(I), I=1,M+1)

101 FORMAT (F10.6)

DEALLOCATE (A,B,XI,YI, T, IPIV)

END PROGRAM LLS

!END OF FILE ***********************************************************