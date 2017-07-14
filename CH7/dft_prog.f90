!*****************************************************************

!Program to perform the discrete Fourier transform of the

!example function. Here we separate out the real and imaginary

!parts of the function and the transform. The inverse transform

!is also performed for comparison.

!*****************************************************************

PROGRAM DFT_EXAMPLE

IMPLICIT NONE

INTEGER, PARAMETER :: N=64,M=8

DOUBLE PRECISION, DIMENSION (N) :: FR,FI,GR,GI

DOUBLE PRECISION, EXTERNAL :: FUNC

INTEGER I

DOUBLE PRECISION F0,H,X

CHARACTER(LEN=15), PARAMETER :: FNAME = 'dft_example.txt'

OPEN(UNIT=1,FILE=FNAME)

! Initialise the variables.

! Here we compute a single factor of 1/N for convenience

F0 = 1.0/N

H  = 1.0/(N-1)

DO I = 1, N

X = H*(I-1)

FR(I) = FUNC(X) ! Re(fm)

FI(I) = 0.0     ! Im(fm)

END DO

! Perform the Fourier transform

CALL DFT (FR,FI,GR,GI,N)

! Remember to include the factor 1/N

! We negate the imaginary part of g in preparation for

! the inverse transform.

DO I = 1, N

GR(I) = F0*GR(I)

GI(I) = -F0*GI(I)

END DO

! Perform the inverse transform

CALL DFT (GR,GI,FR,FI,N)

! The factor here is one as we have used 1/N previously

! Again negate the imaginary part to obtain the correct sign.

! Done for completeness.

DO I = 1, N

FI(I) = -FI(I)

END DO

! Output the transformed data and function values

! to file in steps of M

WRITE (1,"(3F16.8)") (H*(I-1),FR(I),FUNC(H*(I-1)),I=1,N,M)

! Capture the last data point

WRITE (1,"(3F16.8)") H*(N-1),FR(N),FUNC(H*(N-1))

END PROGRAM DFT_EXAMPLE

!************************************************************

SUBROUTINE DFT (FR,FI,GR,GI,N)

! Subroutine to perform the discrete Fourier transform with

! FR and FI as the real and imaginary parts of the signal and

! GR and GI as the corresponding parts of the transform.

IMPLICIT NONE

INTEGER I,J,N

DOUBLE PRECISION PI,X,Q

DOUBLE PRECISION, DIMENSION (N) :: FR,FI

DOUBLE PRECISION, DIMENSION (N) :: GR,GI

PI = 4.D0*ATAN(1.D0)

X  = 2.D0*PI/N

DO I = 1, N

GR(I) = 0.0

GI(I) = 0.0

DO J = 1, N

Q = X*(J-1)*(I-1)

GR(I) = GR(I)+FR(J)*COS(Q)+FI(J)*SIN(Q)

GI(I) = GI(I)+FI(J)*COS(Q)-FR(J)*SIN(Q)

END DO

END DO

END SUBROUTINE DFT

!**************************************************************

DOUBLE PRECISION FUNCTION FUNC(X)

DOUBLE PRECISION X, PI, SI, MU

! Gaussian distribution function

MU = 5.D-1 ! Mean

SI = 1.D-1 ! sigma

PI = 4.D0*ATAN(1.D0)

FUNC = (1/SI/SQRT(2.D0*PI))*EXP(-(X-MU)*(X-MU)/2.D0/SI/SI)

END FUNCTION FUNC

!END OF FILE***************************************************