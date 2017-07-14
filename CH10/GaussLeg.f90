! Program to integrate f(x) = exp(x) on the interval [-1,4]

! using Gauss-Legendre quadrature.

! Here we set quadruple precision

PROGRAM GAUSS

IMPLICIT NONE

INTEGER, PARAMETER :: p = 16 ! quadruple precision

INTEGER            :: N = 10, K

REAL(kind=p), ALLOCATABLE :: R(:,:)

REAL(kind=p)       :: Z, A, B, EXACT

! Perform N point quadrature

DO N = 1,20

A = 1.0_p

B = 4.0_p

! R contains abscissas in 1st row and weights in 2nd row

R = GAUSSQUAD(N)

! Compute the integral the dot_product(j,k) = sum over n (jk)

Z = (B-A)/2._p*DOT_PRODUCT(R(2,:),EXP((B+A)/2._p+R(1,:)*(B-A)/2._p))

! Calculate the exact value

EXACT = EXP(B)-EXP(A)

! Print to screen the results for each quadrature.

PRINT "(I0,1X,G0,1X,G10.2)",N, Z, Z-EXACT

END DO

CONTAINS

! Function computes Legendre polynomials, their roots

! and corresponding weights

FUNCTION GAUSSQUAD(N) RESULT(R)

INTEGER                 :: N

REAL(kind=p), PARAMETER :: PI = 4._p*ATAN(1._p)

REAL(kind=p)            :: R(2, N), X, F, DF, DX

INTEGER                 :: I,  ITER

REAL(kind = p), ALLOCATABLE :: P0(:), P1(:), TMP(:)

! P0, P1, and TMP contain the polynomial coefficients

! they expand with K up to N. The [] concatenates the arrays.

P0 = [1._p]

P1 = [1._p, 0._p]

! compute polynomial coefficients

DO K = 2, N

TMP = ((2*K-1)*[P1,0._p]-(K-1)*[0._p, 0._p,P0])/K

P0 = P1

P1 = TMP

END DO

DO I = 1, N

! Compute initial guess

X = COS(PI*(I-0.25_p)/(N+0.5_p))

DX = 1._p

ITER = 1

! perform Newton Raphson root search

! here we compute F and DF for the given X

! EXIT once root reaches sufficient precision

! OR the iteration count gt 10.

DO WHILE (ABS(DX).GT.10.0*EPSILON(DX) .AND. ITER .LE. 10)

F = P1(1)

DF = 0._p

!as P1 can expand we set the upper limit as the

!size of P1.

DO K = 2, SIZE(P1)

DF = F + X*DF

F  = P1(K) + X * F

END DO

DX =  F / DF

X = X - DX

ITER = ITER + 1

END DO

!once root is found save to array R and compute weight

R(1,I) = X

R(2,I) = 2.0_p/((1.0_p-X*X)*DF*DF)

END DO

END FUNCTION GAUSSQUAD

END PROGRAM GAUSS

!END OF FILE ***************************************************************