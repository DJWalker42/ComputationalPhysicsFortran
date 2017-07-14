MODULE FFT_MOD

IMPLICIT NONE

CONTAINS

!*****************************************************************

SUBROUTINE FFT (A,M,INV)

! An example of the fast Fourier transform subroutine with N = 2**M.

! A(N) is the complex data in the input and corresponding Fourier

! coefficients in the output. INV is flag to compute inverse(1)

! transform or not(0).

IMPLICIT NONE

INTEGER N,N2,M,I,J,K,L,L1,L2,IP,INV

DOUBLE PRECISION PI, Q, WR, WI

COMPLEX*16 :: A(*), U, W, T

PI = 4.D0*ATAN(1.D0)

N = 2**M

N2 = N/2

! Rearrange the data to the bit reversed order

! This interweaves the even and odd function evaluations

J = 1

DO I = 1, N-1

IF (I.LT.J) THEN

T    = A(J)

A(J) = A(I)

A(I) = T

END IF

K   = N2

DO WHILE (K.LT.J)

J = J - K

K = K/2

END DO

J = J + K

END DO

! Perform additions at all levels with reordered data

! i.e. perform the FFT

L2 = 1

DO L = 1, M

L1 =  L2

L2 =  2*L2

Q = PI/L1

U = ( 1.D0, 0.D0 )

WR = COS(Q)

WI = SIN(Q)

W = DCMPLX( WR, -WI )

IF(INV .EQ. 1) W = DCONJG(W) ! complex conjugate

DO J = 1, L1

DO I = J, N, L2

IP    = I + L1

T     = A(IP)*U

A(IP) = A(I)-T

A(I)  = A(I)+T

END DO

U = U * W

END DO

END DO

! Apply the factor of 1/N if not taking the inverse transform

IF(INV .NE. 1) THEN

DO I = 1, N

A(I) = A(I)/N

END DO

END IF

END SUBROUTINE FFT

END MODULE FFT_MOD

!END OF FILE ************************************************