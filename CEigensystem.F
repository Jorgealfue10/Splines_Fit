* CEigensystem.F
* diagonalization of a complex n-by-n matrix using the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 24 Aug 15 th

#include "diag-f.h"


************************************************************************
** CEigensystem diagonalizes a general complex n-by-n matrix.
** Input:	n, A = n-by-n matrix
** Output:	d = vector of eigenvalues,
**		U = transformation matrix,
** these fulfill
**	d = U A U^-1,  A = U d U^-1,  U A = d U  (UCOLS=0),
**	d = U^-1 A U,  A = U^-1 d U,  A U = U d  (UCOLS=1).

	subroutine CEigensystem(n, A,ldA, d, U,ldU, sort)
	implicit none
	integer n, ldA, ldU, sort
	ComplexType A(ldA,*), U(ldU,*), d(*)

	integer p, q, j
	RealType red, off, thresh, norm
	RealType sqp,sqq
	ComplexType delta, t, s, invc, sx, sy, tx, ty
	ComplexType x, y
	ComplexType ev(2,MAXDIM)

	integer sweep
	common /nsweeps/ sweep

	if( n .gt. MAXDIM ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, n
	  ev(1,p) = 0
	  ev(2,p) = A(p,p)
	  d(p) = ev(2,p)
	enddo

	do p = 1, n
	  do q = 1, n
	    U(q,p) = 0
	  enddo
	  U(p,p) = 1
	enddo

	red = .01D0/n**4

	do sweep = 1, 50
	  off = 0
	  do q = 2, n
	    do p = 1, q - 1
			sqq = Sq(A(p,q))
			sqp = Sq(A(q,p))
	      off = off + sqq+sqp
	    enddo
	  enddo
	  if( .not. off .gt. EPS ) goto 1

	  thresh = 0
	  if( sweep .lt. 4 ) thresh = off*red

	  do q = 2, n
	    do p = 1, q - 1
		  sqq = Sq(A(p,q))
		  sqp = Sq(A(q,p))
	      off = sqq+sqp
		  sqq = Sq(ev(2,q))
		  sqp = Sq(ev(2,p))
	      if( sweep .gt. 4 .and. off .lt. EPS*(sqq+sqp) ) then
c     &              EPS*(Sq(ev(2,p)) + Sq(ev(2,q))) ) then
	        A(p,q) = 0
	        A(q,p) = 0
	      else if( off .gt. thresh ) then
	        delta = A(p,q)*A(q,p)
	        x = .5D0*(ev(2,p) - ev(2,q))
	        y = sqrt(x**2 + delta)
	        t = x - y
	        s = x + y
	        if( Sq(t) .lt. Sq(s) ) t = s

	        t = 1/t
	        delta = delta*t
	        ev(1,p) = ev(1,p) + delta
	        ev(2,p) = d(p) + ev(1,p)
	        ev(1,q) = ev(1,q) - delta
	        ev(2,q) = d(q) + ev(1,q)

	        invc = sqrt(delta*t + 1)
	        s = t/invc
	        t = t/(invc + 1)
	        sx = s*A(p,q)
	        ty = t*A(p,q)
	        sy = s*A(q,p)
	        tx = t*A(q,p)

	        do j = 1, n
	          x = A(j,p)
	          y = A(j,q)
	          A(j,p) = x + sy*(y - ty*x)
	          A(j,q) = y - sx*(x + tx*y)
	          x = A(p,j)
	          y = A(q,j)
	          A(p,j) = x + sx*(y - tx*x)
	          A(q,j) = y - sy*(x + ty*y)
	        enddo

	        A(p,q) = 0
	        A(q,p) = 0

	        do j = 1, n
	          x = UL(p,j)
	          y = UL(q,j)
#if UCOLS
	          UL(p,j) = x + sy*(y - ty*x)
	          UL(q,j) = y - sx*(x + tx*y)
#else
	          UL(p,j) = x + sx*(y - tx*x)
	          UL(q,j) = y - sy*(x + ty*y)
#endif
	        enddo
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    d(p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in CEigensystem"

1	continue

* normalize the eigenvectors
	do p = 1, n
	  norm = 0
	  do q = 1, n
	    norm = norm + Sq(UL(p,q))
	  enddo
	  norm = 1/sqrt(norm)
	  do q = 1, n
	    UL(p,q) = UL(p,q)*norm
	  enddo
	enddo

	if( sort .eq. 0 ) return

* sort the eigenvalues by their real part
	do p = 1, n - 1
	  j = p
	  t = d(p)
	  do q = p + 1, n
	    if( sort*(Re(t) - Re(d(q))) .gt. 0 ) then
	      j = q
	      t = d(q)
	    endif
	  enddo

	  if( j .ne. p ) then
	    d(j) = d(p)
	    d(p) = t
	    do q = 1, n
	      x = UL(p,q)
	      UL(p,q) = UL(j,q)
	      UL(j,q) = x
	    enddo
	  endif
	enddo
	end
