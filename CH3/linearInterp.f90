!***************************************************************

! Program to calculate the linear interpolation of the function

! f(x) = sinc(x^2) using 10 points as simulated measured data.

! Data is saved to a text file for graphical use elsewhere --

! gnuplot, Excel, MATLAB,...

!***************************************************************

PROGRAM LINEARINTERP

IMPLICIT NONE

! Line spaces represented by arrays, pts defines the number of

! points to use.

INTEGER,          PARAMETER :: PTS = 100

INTEGER,          PARAMETER :: DPTS = 10

DOUBLE PRECISION, PARAMETER :: INTERVAL = 5.D0

! Declare the variables you are going to use.

! single dimension arrays of size pts

DOUBLE PRECISION, DIMENSION(PTS)     :: X,G1, G2

DOUBLE PRECISION, DIMENSION(DPTS, 2) :: D

DOUBLE PRECISION F, ALPHA, P

INTEGER I, K, M

EXTERNAL F

ALPHA = INTERVAL/PTS

M = (PTS - 1)/(DPTS -1)

DO I = 1, PTS

P = I            ! Assignment converts int into real

X(I) =  P*ALPHA  ! Set up the x values which are floats.

END DO

! Assign the data points we are taking as "measurements"

! Here we use 10 equidistant points

DO  K = 1, DPTS

D(K,1) = X((K-1)*M +1)

D(K,2) = F( D(K,1) )

END DO

CALL LINEAR(X, D, G1, G2, M, PTS, DPTS)

! Write the data to a text file.

! fw.d float total-width . decimal-range

! The output will be space delimited using 2X

OPEN(UNIT=1, FILE="linear_interp_data.txt")

WRITE(1,100) (X(I), F(X(I)), G1(I), G2(I), I=1,PTS)

100 FORMAT(4(f16.10,2X))

CLOSE(1) ! File must be closed prior to exiting

END PROGRAM

!***************************************************************

DOUBLE PRECISION FUNCTION F(X)

DOUBLE PRECISION X

F = SIN(X*X)/(X*X)

END FUNCTION

!***************************************************************

SUBROUTINE LINEAR(X, D, G1, G2, M, PTS, DPTS)

! Declare the variables you are going to use.

! single dimension arrays of size pts

INTEGER I, J, K, M, PTS, DPTS

DOUBLE PRECISION, DIMENSION(PTS)    :: X, G1, G2

DOUBLE PRECISION, DIMENSION(DPTS,2) :: D

DOUBLE PRECISION C1, C2

! Assign the data points we are taking as "measurements"

! to our approximation function g(x)

DO K = 1, DPTS

G1((K-1)*M +1) = D(K, 2)

G2((K-1)*M +1) = D(K, 2)

END DO

! Calculate the linear interpolation. Max J is reduced to dpts-1 as

! index I would go past last "measured" data point otherwise.

! The values for I select the x and f(x) between "measured" points

DO J = 1, DPTS-1

DO I = (J-1)*M+2, M*J

! Nested loop. Here j remains constant while i iterates

! d[j,1], d[j+1,1], d[j,2], and d[j+1,2] are our

! x and f(x) "measurement" values, respectively

C1 = ( X(I) - D(J+1,1) )/( D(J,1)- D(J+1,1) )

C2 = ( X(I) - D(J,1) ) /( D(J+1,1) - D(J,1) )

G1(I) = D(J,2) + C2*(D(J+1,2) - D(J,2)) ! non-symmetrical form

G2(I) = C1*D(J,2) + C2*D(J+1,2)         ! symmetrical form

END DO

END DO

RETURN

END SUBROUTINE

!END OF FILE ***************************************************