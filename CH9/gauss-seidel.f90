SUBROUTINE GAUSS_SIEDEL( F, N, A, B, C, D, X1, XN, TOL, FNAME)

IMPLICIT NONE

! F(N) on entry contains initial guess, on exit contains iterated value

! N dimension of the problem; H = (XN - X1)/(N-1)

! A thru D are the coefficients of the general second order ODE

!          aF'' + bF' + cF = dX

! X1, XN start and end of computational domain on x

! TOL tolerance we want in the iteration will exit once met or itcount > itmax

! FNAME file name to which to store data

INTEGER :: N

DOUBLE PRECISION :: F(N), A, B, C, D, X1, XN, TOL

CHARACTER(*) :: FNAME

LOGICAL DONE

DOUBLE PRECISION H, FF, THE, PSI, PHI, X

INTEGER, PARAMETER :: ITMAX = 100

INTEGER I, J, ITCOUNT

DOUBLE PRECISION Z(ITMAX+1,N)

OPEN(UNIT=1,FILE=FNAME)

100 FORMAT(9999(2X,F10.6))

! calculate parameters and initialise variables

H = (XN - X1)/(N-1)

THE = C - 2*A/H/H

PHI = (A/H/H) + (B/2/H)

PSI = (A/H/H) - (B/2/H)

ITCOUNT = 0

DONE = .FALSE.

! perform the iteration loop

DO WHILE (DONE .EQV. .FALSE. .AND. ITCOUNT .LE. ITMAX)

DONE = .TRUE.

ITCOUNT = ITCOUNT + 1

DO I = 2, N -1

X = X1 + H * (I - 1)

FF = -( PHI * F(I+1) + PSI * F(I-1) -D*X )/THE

IF ( ABS( (FF - F(I))/FF ) .GT. TOL ) DONE = .FALSE.

F(I) = FF

! use z to print data to file

Z(ITCOUNT, I) = F(I)

END DO

END DO

! print message to screen if tolerance not reached

IF(ITCOUNT .GT. ITMAX) THEN

PRINT *, 'Iteration count gone above iteration maximum'

END IF

! add x1 an xn values to z

DO I = 1,ITCOUNT

Z(I, 1) = F(1)

Z(I, N) = F(N)

END DO

! print results to file

DO I = 1,ITCOUNT

WRITE(1,100) (Z(I,J), J=1,N)

END DO

! print to screen iteration count

PRINT *, 'Number of iterations =',ITCOUNT

CLOSE(1)

END SUBROUTINE GAUSS_SIEDEL

!END OF FILE**********************************************************