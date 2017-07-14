!****************************************************************

! Program to calculate the Lagrange interpolation of the function

! f(x) = sinc(x^2) using 10 points as simulated measured data.

! Polynomials of order n = 1 to 7 are calculated.

! Data is saved to a text file for visuals elsewhere,

! e.g. gnuplot, Excel, MATLAB,...

!******************************************************************

PROGRAM LAGRANGEINTERP

IMPLICIT NONE

! Line spaces represented by arrays defined by following parameters

INTEGER,          PARAMETER :: PTS = 100

INTEGER,          PARAMETER :: DPTS = 10

DOUBLE PRECISION, PARAMETER :: INTERVAL = 5.0

INTEGER,          PARAMETER :: N_MAX = 9

! Declare the variables you are going to use.

DOUBLE PRECISION, DIMENSION(PTS)          :: X

DOUBLE PRECISION, DIMENSION(PTS, N_MAX)   :: G

DOUBLE PRECISION, DIMENSION(DPTS, 2)      :: D

INTEGER M                 ! Multiplier to obtain equidistant points

INTEGER I,J,K,N           ! Loop index integers

DOUBLE PRECISION ALPHA, P, F

EXTERNAL F

ALPHA = INTERVAL/PTS

M = (PTS - 1)/(DPTS -1)

DO I = 1,PTS

P = I

X(I) = P*ALPHA ! Set up the x values which are floats.

END DO

! Assign the data points we are taking as "measurements"

! Here we use 10 equidistant points

DO  K = 1, DPTS

D(K,1) = X((K-1)*M +1)

D(K,2) = F( D(K,1) )

END DO

! Initialise the approximating functions to zero

DO J = 1, N_MAX

DO I = 1, PTS

G(I,J) = 0.0

END DO

END DO

! Loop over polynomial orders n=1 to n=N_MAX=7

DO N = 1, N_MAX

CALL LAGRANGE(X, D, G, M, PTS, DPTS, N)

END DO

! Write the data to a text file.

! The output will be space delimited columns.

OPEN(UNIT=1,FILE="lagrange_interp_data.txt")

WRITE(1,100)( X(I), F(X(I)), G(I,2), G(I,4), &

G(I,6), G(I,8), G(I,9), I = 1, PTS)

100 FORMAT(7(2X,F16.10))

CLOSE(1) ! File must be closed prior to exiting

END PROGRAM

!**********************************************************************

SUBROUTINE LAGRANGE(X, D, G, M, PTS, DPTS, N)

! Declare the variables you are going to use.

INTEGER PTS, DPTS, N, M

DOUBLE PRECISION, DIMENSION(PTS)      :: X

DOUBLE PRECISION, DIMENSION(PTS, N)   :: G

DOUBLE PRECISION, DIMENSION(DPTS, 2)  :: D

DOUBLE PRECISION, DIMENSION(DPTS, PTS):: LAM

INTEGER I,J,K,L

! Initialise the Lagrange coefficients to one

DO J = 1,PTS

DO I = 1, DPTS

LAM(I,J) = 1.0

END DO

END DO

!Loop for the Lagrange interpolation

DO J = 1, DPTS-N, N              ! intervals; stride = n

DO I = (J-1)*M+2, M*(J-1+N)+1 ! values between intervals

DO K = J, J+N              ! data points defining interval

DO L = J, J+N           ! loop to determine coefficients

! ensure we miss the case when k == l

IF(K.NE.L)THEN

LAM(K,I) = LAM(K,I)*(X(I)-D(L,1))/(D(K,1)-D(L,1))

END IF

END DO ! loop l

! Calculate the approximation.

G(I,N) = G(I,N) + LAM(K,I)*D(K,2)

END DO ! loop k

END DO ! loop i

END DO ! loop j

RETURN

END SUBROUTINE

!********************************************************************

DOUBLE PRECISION FUNCTION F(X)

DOUBLE PRECISION X

F = SIN(X*X)/(X*X)

END FUNCTION

!END OF FILE*********************************************************