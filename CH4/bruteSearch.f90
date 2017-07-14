!****************************************************

SUBROUTINE SEARCH ( A, B, STEP, MAXVAL, N, M, F )

! Searches for the roots of the function F on the

! interval [A,B] using STEP size. There should be

! N roots in this interval. M found roots returned.

DOUBLE PRECISION, DIMENSION (N) :: A, B

DOUBLE PRECISION STEP, FA, FB, MAXVAL

INTEGER I, RCOUNT, M

! Min value store in A(1)

! Evaluate function at minval for initial check

FA = F(A(1))

M = N

RCOUNT = 0

I = 1

! Terminate loop if the search reaches the end of the interval

DO WHILE ( B(I).LT.MAXVAL )

! Move the search forward one step

B(I) = A(I) + STEP

! Evaluate the function at new position

FB = F( B(I) )

! Check for existence of root

IF (FA * FB .LT. 0)THEN

! Add 1 to root count

RCOUNT = RCOUNT + 1

! Exit if reached expected no of roots

IF(RCOUNT .EQ. N)RETURN

! Update function value

FA = FB

! Move to next root

I = I + 1

! Assign current search pos. to new search pos.

A(I) = B(I - 1)

ELSE

! Update function value

FA = FB

! Store current position

A(I) = B(I)

END IF

END DO

! Arrive here if root count not met within defined interval

PRINT *, 'Warning: Fewer roots found than predicted in interval'

PRINT *, 'No of roots found:', RCOUNT

M = RCOUNT

RETURN

END SUBROUTINE SEARCH

!END OF FILE *************************************************************