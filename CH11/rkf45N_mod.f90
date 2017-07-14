!************************************************************

!Module containing Runge-Kutta-Fehlberg algorithms to deal

!with both first (RK4F1) and second ordered (RK4FN) ODEs.

!The 2nd ordered ODEs are reduced to coupled 1st order ODEs.

!The R-K methods used are fourth/fifth order.

!The subroutine DERIV calculates the derivatives of all

!components of the derivative vector

!************************************************************

MODULE RKF45_MODULE

IMPLICIT NONE

CONTAINS

!********************************************************

SUBROUTINE RKF451( DF, X, Y, A, B, TOL )

!******************************************************

!RK4F1 uses a fourth/fifth order Runge-Kutta-Fehlberg

!algorithm to solve first order ODEs.

!******************************************************

!Declare subroutine arguments

DOUBLE PRECISION DF, X, Y, A, B, TOL

EXTERNAL DF

!Adjustable parameters *******************************

!Define the conservative factor for new step size

DOUBLE PRECISION, PARAMETER :: ALPHA = 9.D-1

!FRACMAX: Fraction of interval defining HMAX (B-A/FRACMAX)

!FRACMIN: Fraction of interval defining HMIN (B-A/FRACMIN)

DOUBLE PRECISION, PARAMETER :: FRACMAX = 1.D1

DOUBLE PRECISION, PARAMETER :: FRACMIN = 1.D5

!INC: Max tolerated increase factor of h

!DEC: Max tolerated decrease factor of h

DOUBLE PRECISION, PARAMETER :: INC = 4.D0

DOUBLE PRECISION, PARAMETER :: DEC = 1.D-1

!DO NOT CHANGE the coefficients below ***************

! Am define the coefficients used to calculate x

DOUBLE PRECISION, PARAMETER :: A1 = 2.5D-1

DOUBLE PRECISION, PARAMETER :: A2 = 3.75D-1

DOUBLE PRECISION, PARAMETER :: A3 = 1.2D1/1.3D1

DOUBLE PRECISION, PARAMETER :: A4 = 1.D0

DOUBLE PRECISION, PARAMETER :: A5 = 5.D-1

!Bmn define the coefficients used to calculate y

DOUBLE PRECISION, PARAMETER :: B10 = 2.5D-1

DOUBLE PRECISION, PARAMETER :: B20 = 3.D0/3.2D1

DOUBLE PRECISION, PARAMETER :: B21 = 9.D0/3.2D1

DOUBLE PRECISION, PARAMETER :: B30 = 1.932D3/2.197D3

DOUBLE PRECISION, PARAMETER :: B31 = -7.2D3/2.197D3

DOUBLE PRECISION, PARAMETER :: B32 = 7.296D3/2.2197D3

DOUBLE PRECISION, PARAMETER :: B40 = 4.39D2/2.16D2

DOUBLE PRECISION, PARAMETER :: B41 = -8.D0

DOUBLE PRECISION, PARAMETER :: B42 = 3.68D3/5.13D2

DOUBLE PRECISION, PARAMETER :: B43 = -8.45D2/4.104D3

DOUBLE PRECISION, PARAMETER :: B50 = -8.D0/2.7D1

DOUBLE PRECISION, PARAMETER :: B51 = 2.D0

DOUBLE PRECISION, PARAMETER :: B52 = -3.544D3/2.565D3

DOUBLE PRECISION, PARAMETER :: B53 = 1.859D3/4.104D3

DOUBLE PRECISION, PARAMETER :: B54 = -1.1D1/4.D1

!Cm are used to evaluate YHAT

DOUBLE PRECISION, PARAMETER :: C0 = 1.6D1/1.35D2

DOUBLE PRECISION, PARAMETER :: C2 = 6.656D3/1.2825D4

DOUBLE PRECISION, PARAMETER :: C3 = 2.8561D4/5.643D4

DOUBLE PRECISION, PARAMETER :: C4 = -9.D0/5.D1

DOUBLE PRECISION, PARAMETER :: C5 = 2.D0/5.5D1

!Dm are used to calculate the error

DOUBLE PRECISION, PARAMETER :: D0 = 1.D0/3.6D2

DOUBLE PRECISION, PARAMETER :: D2 = -1.28D2/4.275D3

DOUBLE PRECISION, PARAMETER :: D3 = -2.197D3/7.524D4

DOUBLE PRECISION, PARAMETER :: D4 = 1.D0/5.D1

DOUBLE PRECISION, PARAMETER :: D5 = 2.D0/5.5D1

!Declare local variables

DOUBLE PRECISION H, HMAX, HNEW, HMIN, X0, Y0, YHAT

DOUBLE PRECISION K0, K1, K2, K3, K4, K5, ERR

IF(A.GE.B)STOP 'User error: B must be .gt. A'

IF(TOL.LT.EPSILON(0.D0))THEN

PRINT*, 'User error: Tol must be .gt.',EPSILON(0.D0)

STOP

END IF

!Open a file to save data

OPEN(UNIT=1, FILE='rk4f_data.txt')

!Set maximum step size 1/10 of full range

HMAX = (B - A)/FRACMAX

!Set minimum step size 1/10,000 of full range

HMIN = (B - A)/FRACMIN

!Initialise the integration;

!initial y defined in calling program.

H = HMAX

X0 = A

Y0 = Y

!Perform the integration

DO WHILE (X0.LE.B)

K0 = DF(X0,Y0)

X = X0 + A1*H

Y = Y0 + B10*H*K0

K1 = DF(X,Y)

X = X0 + A2*H

Y = Y0 + H*(B20*K0 + B21*K1)

K2 = DF(X,Y)

X = X0 + A3*H

Y = Y0 + H*(B30*K0 + B31*K1 + B32*K2)

K3 = DF(X,Y)

X = X0 + A4*H

Y = Y0 + H*(B40*K0 + B41*K1 + B42*K2 + B43*K3)

K4 = DF(X,Y)

X = X0 + A5*H

Y = Y0 + H*(B50*K0 + B51*K1 + B52*K2 + B53*K3 &

+ B54*K4)

K5 = DF(X,Y)

YHAT = Y0 + H*(C0*K0 + C2*K2 + C3*K3 + C4*K4 &

+ C5*K5)

ERR = H * ABS(D0*K0 + D2*K2 + D3*K3 + D4*K4 &

+ D5*K5)

HNEW = ALPHA * H * SQRT( SQRT( H*TOL/ERR ) )

!Check prediction of step size falls within

!acceptable limits and adjust if necessary

IF( HNEW .GT. INC * H )HNEW = INC * H

IF( HNEW .LT. DEC * H )HNEW = DEC * H

IF( HNEW .GT. HMAX)HNEW = HMAX

IF( HNEW .LT. HMIN)THEN

PRINT *, H,'->', HNEW

PRINT *,'Step size too small; exiting integration'

WRITE (1,100) X, YHAT

STOP 'Possible problem with RK4F subroutine'

END IF

!Check error vs. max_error = h*eps

!if error too large repeat step using hnew

!else accept propagation

IF( ERR .GT. (H*TOL) )THEN

H = HNEW

ELSE

X0 = X0 + H

Y0 = YHAT

H = HNEW

!Save the data to file

WRITE(1,100) X0, Y0

END IF

END DO

100 FORMAT(2(2X,F10.6))

!Assign calculated values back to subroutine arguments

X = X0

Y = Y0

CLOSE(1)

RETURN

END SUBROUTINE RKF451

!************************************************************

SUBROUTINE RKF45N( N, DF, X, Y, A, B, TOL, FNAME, FREQ )

!******************************************************

!RK4FN uses a fourth/fifth order Runge-Kutta-Fehlberg

!algorithm to solve second order ODEs.

!2nd order ODE -> pair of coupled 1st order ODEs

!for any number of dependent functions

!backwards integration is allowed i.e. B .lt. A

!******************************************************

!Declare subroutine arguments

INTEGER FREQ, N

DOUBLE PRECISION DF, X, Y(N), A, B, TOL

CHARACTER(*) FNAME

EXTERNAL DF

!Define the conservative factor for new step size

DOUBLE PRECISION, PARAMETER :: ALPHA = 9.D-1

! Am define the coefficients used to calculate x

DOUBLE PRECISION, PARAMETER :: A1 = 2.5D-1

DOUBLE PRECISION, PARAMETER :: A2 = 3.75D-1

DOUBLE PRECISION, PARAMETER :: A3 = 1.2D1/1.3D1

DOUBLE PRECISION, PARAMETER :: A4 = 1.D0

DOUBLE PRECISION, PARAMETER :: A5 = 5.D-1

!Bmn define the coefficients used to calculate y

DOUBLE PRECISION, PARAMETER :: B10 = 2.5D-1

DOUBLE PRECISION, PARAMETER :: B20 = 3.D0/3.2D1

DOUBLE PRECISION, PARAMETER :: B21 = 9.D0/3.2D1

DOUBLE PRECISION, PARAMETER :: B30 = 1.932D3/2.197D3

DOUBLE PRECISION, PARAMETER :: B31 = -7.2D3/2.197D3

DOUBLE PRECISION, PARAMETER :: B32 = 7.296D3/2.2197D3

DOUBLE PRECISION, PARAMETER :: B40 = 4.39D2/2.16D2

DOUBLE PRECISION, PARAMETER :: B41 = -8.D0

DOUBLE PRECISION, PARAMETER :: B42 = 3.68D3/5.13D2

DOUBLE PRECISION, PARAMETER :: B43 = -8.45D2/4.104D3

DOUBLE PRECISION, PARAMETER :: B50 = -8.D0/2.7D1

DOUBLE PRECISION, PARAMETER :: B51 = 2.D0

DOUBLE PRECISION, PARAMETER :: B52 = -3.544D3/2.565D3

DOUBLE PRECISION, PARAMETER :: B53 = 1.859D3/4.104D3

DOUBLE PRECISION, PARAMETER :: B54 = -1.1D1/4.D1

!Cm are used to evaluate YHAT

DOUBLE PRECISION, PARAMETER :: C0 = 1.6D1/1.35D2

DOUBLE PRECISION, PARAMETER :: C2 = 6.656D3/1.2825D4

DOUBLE PRECISION, PARAMETER :: C3 = 2.8561D4/5.643D4

DOUBLE PRECISION, PARAMETER :: C4 = -9.D0/5.D1

DOUBLE PRECISION, PARAMETER :: C5 = 2.D0/5.5D1

!Dm are used to calculate the error

DOUBLE PRECISION, PARAMETER :: D0 = 1.D0/3.6D2

DOUBLE PRECISION, PARAMETER :: D2 = -1.28D2/4.275D3

DOUBLE PRECISION, PARAMETER :: D3 = -2.197D3/7.524D4

DOUBLE PRECISION, PARAMETER :: D4 = 1.D0/5.D1

DOUBLE PRECISION, PARAMETER :: D5 = 2.D0/5.5D1

!Declare local variables

DOUBLE PRECISION H, HMAX, HNEW, HMIN, X0, ERR, BIGERR

DOUBLE PRECISION, DIMENSION(N) :: Y0, YHAT, K0, K1

DOUBLE PRECISION, DIMENSION(N) :: K2, K3, K4, K5

INTEGER I, COUNT

IF(MOD(N,2).NE.0) STOP 'User error: N must be even'

IF(TOL.LT.EPSILON(0.D0))THEN

PRINT*, 'User error: Tol must be .gt.',EPSILON(0.D0)

STOP

END IF

IF(FREQ.GT.1000)THEN

STOP 'User error: only a MAX of 1000 data points allowed'

END IF

!Open a file to save data

OPEN(UNIT=1, FILE=FNAME)

!Initialise print count

COUNT = 1

!Set maximum step size as 1/freq of full range

!i.e. freq is the number of points to save

HMAX = (B - A)/FREQ

!Set minimum step size to a reasonable (+ve) value say,

HMIN = 1.D1*EPSILON(0.D0)

!Initialise the integration;

!initial Y values defined in calling program.

H = HMAX

X0 = A

DO I = 1,N

Y0(I) = Y(I)

END DO

!Save the initial conditions

WRITE(1,100) X0, (Y0(I), I=1,N)

!Perform the integration

DO WHILE ( ABS(X0-A).LT.ABS(B-A) )

CALL DERIV( N, X0, Y0, K0 )

X = X0 + A1*H

DO I = 1,N

Y(I) = Y0(I) + B10*H*K0(I)

END DO

CALL DERIV( N, X, Y, K1 )

X = X0 + A2*H

DO I = 1,N

Y(I) = Y0(I) + H*(B20*K0(I) + B21*K1(I))

END DO

CALL DERIV( N, X, Y, K2 )

X = X0 + A3*H

DO I = 1,N

Y(I) = Y0(I) + H*(B30*K0(I) + B31*K1(I) &

+ B32*K2(I))

END DO

CALL DERIV( N, X, Y, K3 )

X = X0 + A4*H

DO I = 1,N

Y(I) = Y0(I) + H*(B40*K0(I) + B41*K1(I) &

+ B42*K2(I) + B43*K3(I))

END DO

CALL DERIV( N, X, Y, K4 )

X = X0 + A5*H

DO I = 1,N

Y(I) = Y0(I) + H*(B50*K0(I) + B51*K1(I) &

+ B52*K2(I) + B53*K3(I)  + B54*K4(I))

END DO

CALL DERIV( N, X, Y, K5 )

BIGERR = 0.D0

DO I = 1,N

YHAT(I) = Y0(I) + H*(C0*K0(I) + C2*K2(I) &

+ C3*K3(I) + C4*K4(I)  + C5*K5(I))

ERR = H * ABS( D0*K0(I) + D2*K2(I) + D3*K3(I) &

+ D4*K4(I) + D5*K5(I) )

IF( ABS(ERR) .GT. ABS(BIGERR) )BIGERR = ERR

END DO

HNEW = ALPHA * H * SQRT( SQRT( H*TOL/BIGERR ) )

!Check prediction of step size falls within

!acceptable limits and adjust if necessary

IF( ABS(HNEW) .GT. 4.D0 * ABS(H) ) HNEW = 4.D0 * H

IF( ABS(HNEW) .LT. .1D0 * ABS(H) ) HNEW = .1D0 * H

IF( ABS(HNEW) .GT. ABS(HMAX) )     HNEW = HMAX

! No ABS required for HMIN as it is a +ve value

IF( ABS(HNEW) .LT. HMIN)THEN

PRINT *, H,'->', HNEW

WRITE (1,100) X, (YHAT(I),I=1,N)

STOP 'Step size too small; exiting integration'

END IF

!Check error vs. Max error = h*eps

!if error too large repeat step using hnew

!else accept propagation

IF( ABS(BIGERR) .GT. ABS(H)*TOL )THEN

!check new step doesn't take us past B

IF(ABS(X0 - A + HNEW).GT. ABS(B - A)) HNEW = B - X0

H = HNEW

ELSE

X0 = X0 + H

DO I = 1,N

Y0(I) = YHAT(I)

END DO

!check new step doesn't take us past B

IF(ABS(X0 - A + HNEW).GT. ABS(B - A)) HNEW = B - X0

H = HNEW

!Save the data to file

IF( ABS(X0 - A)/COUNT .GE. ABS(HMAX) )THEN

WRITE(1,100) X0, (Y0(I),I=1,N)

COUNT = COUNT + 1

END IF

END IF

END DO

100 FORMAT(999(2X,ES12.5))

!Assign calculated values back to subroutine arguments

X = X0

DO I = 1,N

Y(I) = Y0(I)

END DO

!write final values to file

WRITE(1,100) X0, (Y0(I),I=1,N)

!Add a blank line at end of data file before closing

WRITE(1,*)

CLOSE(1)

RETURN

END SUBROUTINE RKF45N

!*******************************************************************

SUBROUTINE DERIV( N, X, Y, K )

INTEGER N,I,N2

DOUBLE PRECISION X, Y(N), K(N), DF

EXTERNAL DF

N2 = N/2

DO I = 1,N2

K(I) = Y(I+N2)

K(I+N2) = DF( X, Y, I )

END DO

RETURN

END SUBROUTINE DERIV

!********************************************************************

END MODULE RKF45_MODULE

!END OF FILE ******************************************************

rkf45FFT_mod.f90

!************************************************************

!Module containing Runge-Kutta-Fehlberg algorithm to produce

!data for the FFT subroutine.

!************************************************************

MODULE RKF45FFT_MODULE

IMPLICIT NONE

CONTAINS

!************************************************************

SUBROUTINE RKF45FFT( DF, X, Y, A, B, TOL, DATA, DELTA )

!******************************************************

!RK4F2 uses a fourth/fifth order Runge-Kutta-Fehlberg *

!algorithm to solve second order ODEs.                *

!2nd order ODE -> pair of coupled 1st order ODEs      *

!******************************************************

!Declare subroutine arguments

DOUBLE PRECISION DF, X, Y(2), A, B, TOL, DELTA

EXTERNAL DF

!Declare variables for FFT sampling

DOUBLE PRECISION, PARAMETER :: GOAL = 100

INTEGER,          PARAMETER :: NFFT = 512

COMPLEX*16 :: DATA(NFFT*2)

INTEGER IFFT

IFFT = 1

!Define the conservative factor for new step size

DOUBLE PRECISION, PARAMETER :: ALPHA = 9.D-1

!FRACMAX: Fraction of interval defining HMAX (B-A/FRACMAX)

!FRACMIN: Fraction of interval defining HMIN (B-A/FRACMIN)

DOUBLE PRECISION, PARAMETER :: FRACMAX = 1.D1

DOUBLE PRECISION, PARAMETER :: FRACMIN = 1.D5

!INC: Max tolerated increase factor of h

!DEC: Max tolerated decrease factor of h

DOUBLE PRECISION, PARAMETER :: INC = 4.D0

DOUBLE PRECISION, PARAMETER :: DEC = 1.D-1

! Am define the coefficients used to calculate x

DOUBLE PRECISION, PARAMETER :: A1 = 2.5D-1

DOUBLE PRECISION, PARAMETER :: A2 = 3.75D-1

DOUBLE PRECISION, PARAMETER :: A3 = 1.2D1/1.3D1

DOUBLE PRECISION, PARAMETER :: A4 = 1.D0

DOUBLE PRECISION, PARAMETER :: A5 = 5.D-1

!Bmn define the coefficients used to calculate y

DOUBLE PRECISION, PARAMETER :: B10 = 2.5D-1

DOUBLE PRECISION, PARAMETER :: B20 = 3.D0/3.2D1

DOUBLE PRECISION, PARAMETER :: B21 = 9.D0/3.2D1

DOUBLE PRECISION, PARAMETER :: B30 = 1.932D3/2.197D3

DOUBLE PRECISION, PARAMETER :: B31 = -7.2D3/2.197D3

DOUBLE PRECISION, PARAMETER :: B32 = 7.296D3/2.197D3

DOUBLE PRECISION, PARAMETER :: B40 = 4.39D2/2.16D2

DOUBLE PRECISION, PARAMETER :: B41 = -8.D0

DOUBLE PRECISION, PARAMETER :: B42 = 3.68D3/5.13D2

DOUBLE PRECISION, PARAMETER :: B43 = -8.45D2/4.104D3

DOUBLE PRECISION, PARAMETER :: B50 = -8.D0/2.7D1

DOUBLE PRECISION, PARAMETER :: B51 = 2.D0

DOUBLE PRECISION, PARAMETER :: B52 = -3.544D3/2.565D3

DOUBLE PRECISION, PARAMETER :: B53 = 1.859D3/4.104D3

DOUBLE PRECISION, PARAMETER :: B54 = -1.1D1/4.D1

!Cm are used to evaluate YHAT

DOUBLE PRECISION, PARAMETER :: C0 = 1.6D1/1.35D2

DOUBLE PRECISION, PARAMETER :: C2 = 6.656D3/1.2825D4

DOUBLE PRECISION, PARAMETER :: C3 = 2.8561D4/5.643D4

DOUBLE PRECISION, PARAMETER :: C4 = -9.D0/5.D1

DOUBLE PRECISION, PARAMETER :: C5 = 2.D0/5.5D1

!Dm are used to calculate the error

DOUBLE PRECISION, PARAMETER :: D0 = 1.D0/3.6D2

DOUBLE PRECISION, PARAMETER :: D2 = -1.28D2/4.275D3

DOUBLE PRECISION, PARAMETER :: D3 = -2.197D3/7.524D4

DOUBLE PRECISION, PARAMETER :: D4 = 1.D0/5.D1

DOUBLE PRECISION, PARAMETER :: D5 = 2.D0/5.5D1

!Declare local variables

DOUBLE PRECISION H, HMAX, HNEW, HMIN, X0, ERR, BIGERR

DOUBLE PRECISION, DIMENSION(2) :: Y0, YHAT, K0, K1

DOUBLE PRECISION, DIMENSION(2) :: K2, K3, K4, K5

INTEGER I

IF(A.GE.B)STOP 'User error: B must be .gt. A'

IF(TOL.LT.EPSILON(0.D0))THEN

PRINT*, 'User error: Tol must be .gt.',EPSILON(0.D0)

STOP

END IF

!Set maximum step size

HMAX = (B - A)/FRACMAX

!Set minimum step size

HMIN = (B - A)/FRACMIN

!Initialise the integration;

!initial Y values defined in calling program.

H = HMAX

X0 = A

DO I = 1,2

Y0(I) = Y(I)

END DO

! save initial value to data array

DATA(IFFT) = Y(1)

!Perform the integration; exit after collecting required data

DO WHILE (IFFT.LT.NFFT)

CALL DERIV( X0, Y0, K0 )

X = X0 + A1*H

DO I = 1,2

Y(I) = Y0(I) + B10*H*K0(I)

END DO

CALL DERIV( X, Y, K1 )

X = X0 + A2*H

DO I = 1,2

Y(I) = Y0(I) + H*(B20*K0(I) + B21*K1(I))

END DO

CALL DERIV( X, Y, K2 )

X = X0 + A3*H

DO I = 1,2

Y(I) = Y0(I) + H*(B30*K0(I) + B31*K1(I) &

+ B32*K2(I))

END DO

CALL DERIV( X, Y, K3 )

X = X0 + A4*H

DO I = 1,2

Y(I) = Y0(I) + H*(B40*K0(I) + B41*K1(I) &

+ B42*K2(I) + B43*K3(I))

END DO

CALL DERIV( X, Y, K4 )

X = X0 + A5*H

DO I = 1,2

Y(I) = Y0(I) + H*(B50*K0(I) + B51*K1(I) &

+ B52*K2(I) + B53*K3(I)  + B54*K4(I))

END DO

CALL DERIV( X, Y, K5 )

BIGERR = 0.D0

DO I = 1,2

YHAT(I) = Y0(I) + H*(C0*K0(I) + C2*K2(I) &

+ C3*K3(I) + C4*K4(I)  + C5*K5(I))

ERR = H * ABS(D0*K0(I) + D2*K2(I) + D3*K3(I) &

+ D4*K4(I) + D5*K5(I))

IF(ERR .GT. BIGERR)BIGERR = ERR

END DO

HNEW = ALPHA * H * SQRT( SQRT( H*TOL/BIGERR ) )

!Check prediction of step size falls within

!acceptable limits and adjust if necessary

IF( HNEW .GT. INC * H ) HNEW = INC * H

IF( HNEW .LT. DEC * H ) HNEW = DEC * H

IF( HNEW .GT. HMAX) HNEW = HMAX

IF( HNEW .LT. HMIN)THEN

PRINT *, H,'->', HNEW

WRITE (1,100) X, YHAT(1), YHAT(2), YHAT(3), YHAT(4)

STOP 'Step size too small; exiting integration'

END IF

!Check error vs. max_error = h*eps

!if error too large repeat step using hnew

!else accept propagation

IF( BIGERR .GT. (H*TOL) )THEN

H = HNEW

ELSE

X0 = X0 + H

IF( ABS((X0-GOAL)/X0) .LT. 1.D-5)THEN

IFFT = IFFT + 1

DATA(IFFT) = YHAT(1)

GOAL = GOAL + DELTA

END IF

IF(H .GT. GOAL-X0) H = GOAL-X0

DO I = 1,2

Y0(I) = YHAT(I)

END DO

H = HNEW

END IF

END DO

100 FORMAT(5(2X,F10.6))

!Assign calculated values back to subroutine arguments

X = X0

DO I = 1,2

Y(I) = Y0(I)

END DO

!Add a blank line at end of data file before closing

WRITE(1,*)

CLOSE(1)

RETURN

END SUBROUTINE RKF45FFT

!*******************************************************************

SUBROUTINE DERIV( X, Y, K )

DOUBLE PRECISION X, Y(2), K(2), DF

EXTERNAL DF

K(1) = Y(2)

K(2) = DF( X, Y )

RETURN

END SUBROUTINE DERIV

!********************************************************************

END MODULE RKF45FFT_MODULE

!END OF FILE*******************************************************