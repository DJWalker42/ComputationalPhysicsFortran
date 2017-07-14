PROGRAM HEAT_EQUATION

IMPLICIT NONE

INTEGER, PARAMETER :: M = 6, N = 101

DOUBLE PRECISION :: AD(M-2), AS(M-3), B(M-2,M), U(M,N), US(M-2)

DOUBLE PRECISION :: AD2(M-2), AS2(M-3)

DOUBLE PRECISION H, K, R, W, L, T_END, D1, D2

INTEGER I, J, INFO

W = 5.D-1 ! Crank-Nicolson

L = 1.D0  ! length of conducting rod

D1 = 1.D2 ! Dirichlet boundary condtion at x = 0

D2 = 1.D2 ! Dirichlet boundary condtion at x = L

T_END = 1.D0 ! time to run system

H = L/(M-1) ! H = 0.2

K = T_END/(N-1) ! K = 0.01

R = K/H/H ! R = 0.01/0.2**2 = 1/4

! Initialise our variables *****************************

!Set up the diagonal of A.

!Store to another array for later use.

DO I = 1, M-2

AD(I) = 1 + 2 * R * W

AD2(I) = AD(I)

END DO

!set up the sub/super diagonal of A

!Store to another array for later use.

DO I = 1, M-3

AS(I) = -R * W

AS2(I) = AS(I)

END DO

!Initialise B to zero

B = 0.D0

!Set up non-zero values of B

DO I = 1, M -2

B(I,I)     = R * (1 - W)

B(I,I + 1) = 1 - 2 * R *(1 - W)

B(I,I + 2) = R * (1 - W)

END DO

!Adjust first and last entries for boundary values

B(1,1)   = R

B(M-2,M) = R

!Initial conditions: U(x,0) = 0.d0 Celsius

! set the entire system to the IC, advanced time step

! values will be overwritten.

U = 0.D0

!Account for boundary conditions

DO J = 1,N

U(1,J) = D1

U(M,J) = D2

END DO

! Main loop *******************************************

DO J = 1, N-1

!Compute the known right hand sides

DO I = 1, M-2

US(I) = B(I,I)*U(I,J) + B(I,I+1)*U(I+1,J) + &

B(I,I+2)*U(I+2,J)

END DO

!Solve the tridiagonal system using LAPACK routine DPTSV

CALL DPTSV(M-2, 1, AD, AS, US, M-2, INFO)

!Write the solution back to U at the advanced time j+1

!Remember to write to interior grid points i=2, m-1.

DO I = 1, M-2

U(I+1,J+1) = US(I)

END DO

!As the Lapack subroutine overwrites the arrays AD and AS

!with their factors we need to reset these.

!NOTE: THIS IS NOT OPTIMAL. How can it be improved? Hint:

!Take a look at the DPTSV source code located at

!http://www.netlib.org/lapack/double/dptsv.f

DO I = 1, M-2

AD(I) = AD2(I)

END DO

DO I = 1, M-3

AS(I) = AS2(I)

END DO

END DO

! Extract data to file *********************************

OPEN(UNIT=1, FILE='heat_eqn_data.txt')

DO I = 1,M

WRITE (1,'(1001f10.4)') ( U(I,J), J=1,N )

END DO

CLOSE(1)

END PROGRAM HEAT_EQUATION

!END OF FILE *****************************************************