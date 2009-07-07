***
*     Author:         Peter Dunn
*     Creation date:  15 February 1999
*     Last revision:  11 January 2005
***
*
*** 25 NOV 2005:  removed the IF statements in the integration loops
***               where we stopped if flag==1.
***               Now we continue anyway.  We did this because often we
***               stopped after 2, 3 or 4 iterations, without really knowing why.
*
* The entry points are  pdf  and  cdf  and employ the asymptotic Sidi method,
* or the exact zeros method
*
* IN R, print statements for debugging as follows:
*     call dblepr("The value of y is ",-1, y, 1)
*     call intpr("The value of exact is ",-1, exact, 1)
******************************************************************

      subroutine pdf(p, phi, y, mu, exact, verbose,
     &               funvalue, exitstatus, relerr, its )

***
*
*     Calculates the density of the log-likelihood
*     function of a Poisson-gamma distribution by inverting the
*     moment generating function.
*
***
*     IN:   p, phi, y, mu, exact
*     OUT:  funvalue, exitstatus, relerr, its
***
*** NOTE:  WE SHOULD ALWAYS GET mu=1 FROM R

      double precision  p, phi, y, funvalue, mu, savemu,
     &                  lambda, calclambda, pi, area,
     &                  result, relerr, aimrerr
      integer  ier, maxit, m, iteratn, exitstatus,
     &         its, exact, verbose

***
* VARIABLES:
*    y        : the value at which the function is to be evaluated
*    x        : the range over which the integral is to be integrated;
*               an internal variable and NOT the value at which
*               the function is to be evaluated
*    c        : the constant mapping to the distribution with mean=1
*    lambda   : for 1<p<2, P(Y=0) = exp( -lambda )
*    p        : the index (ie variance function is V(mu) = mu^p)
*    phi      : the dispersion parameter
*    funvalue : the value of the function at the given value of  x
*    bound    : The bound using Chebyshev theorm.
*    exitstatus:  1  if relative error is smaller than wished (aimrerr)
*                -1  if not, but the absolute error is less than aimrerr
*               -10  if neither rel or abs error any good
*    exact    : 1 if the exact zero acceleration algorithms is used;
*               0 if the approx zeros acceleration algorithm is used.
*    verbose  : 1 to print lots of diagnostic information
*             : 0 to keep quiet
***

*     Set defaults
      exitstatus = 1
      relerr = 0.0d00
      its = 0
      
*     SPECIAL CASES {
      if ( p .EQ. 1.0d00 ) then
         funvalue = -10.d00
         return
      endif

      if ( ( p .LT. 2.0d00 ) .AND. 
     &     ( p .GT. 1.0d00 ) ) then

         if ( y .LT. 0.0d00 ) then

            funvalue = 0.0d00
            return

         elseif ( y .EQ. 0.0d00 ) then

*           If  y=0 in P-G case, then we don't need to integrate, just
*           evaluate; otherwise, integrate

            lambda = calclambda(p, phi, mu)
            funvalue = exp( -1.0d00 * lambda )

            return

         endif
    
      elseif ( p .GE. 2.0d00 ) then

         if ( y .LE. 0.0d00 ) then

            funvalue = 0.0d00
            return

         endif

      endif

*     } END SPECIAL CASES

***

*     SET ACCURACY REQUIREMENTS
*     maximum number of iterations in calculating errors
      maxit = 100
      aimrerr = 1.0d-10

*     set other parameters
      m = -1
      pi = acos( -1.0d00 )
*      pi = atan( 1.0d00 ) * 4.0d00
*     pi = 3.14159 26535 89793 23846 26433d00
      area = 0.0d00
      iteratn = 0
      relerr = 1.0d00
      ier = 0
      savemu = mu

*     For Tweedie distributions,
*        f(y; \mu, \phi) = cf(cy; c\mu, c^{2-p}\phi)
*     When mu\ne 1, we can set c = 1/\mu and then we
*     evaluate  f( y/\mu, 1, (\mu)^{2-p}\phi)
*     and the actual density should be c times this.
        
     
*     DETERMINE a(y,phi) (except in normal and gamma cases, where
*     the distribution is determined from the closed form result).


      if ( exact .EQ. 1 ) then
*        Use exact zero acceleration algorithm
         if ( (p .GT. 1.0d00).AND.(p .LE. 2.0d00) ) then
            call smallp( p, phi, y, mu, aimrerr, result,
     &                   maxit, ier, exitstatus, relerr, 
     &                   its, verbose )
*           This has found P(y|y>0); must convert back to P(Y)

         elseif ( p .GT. 2.0d00 ) then

            call bigp( p, phi, y, mu, aimrerr, result,
     &                 maxit, ier, exitstatus, relerr, 
     &                 its, verbose )

         endif

      else

*        Use approx zero acceleration algorithm
         call evlint( p, phi, y, mu, aimrerr, result,
     &                maxit, ier, exitstatus, relerr, 
     &                its, verbose )
*         result = -1000.0d00

      endif

      funvalue = result

*     Some tidy-ups
      if ( funvalue .LT. 0.0d00 ) funvalue = 0.0d00

      return
      end

******************************************************************
******************************************************************
******************************************************************

      subroutine cdf(p, phi, y, mu, exact,
     &              funvalue, exitstatus, relerr, its )

***
*
*     Calculates the cdf of the log-likelihood
*     function of a Poisson-gamma distribution by inverting the
*     moment generating function.
*
***
*     IN:   p, phi, y, mu, exact, exact
*     OUT:  funvalue, exitstatus, relerr, its
***

      double precision  p, phi, y, funvalue, mu, 
     &                  lambda, calclambda, resulta,
     &                  result, relerr, aimrerr,
     &                  result0
      integer  ier, maxit, iteratn, exitstatus,
     &         its, exact, verbose

***
* VARIABLES:
*    y        : the value at which the function is to be evaluated
*    x        : the range over which the integral is to be integrated;
*               an internal variable and NOT the value at which
*               the function is to be evaluated
*    c        : the constant mapping to the distribution with mean=1
*    lambda   : for 1<p<2, P(Y=0) = exp( -lambda )
*    p        : the index (ie variance function is V(mu) = mu^p)
*    phi      : the dispersion parameter
*    funvalue : the value of the function at the given value of  x
*    bound    : The bound using Chebyshev theorm.
*    exitstatus:  1  if relative error is smaller than wished (aimrerr)
*                -1  if not, but the absolute error is less than aimrerr
*               -10  if neither rel or abs error any good
*    exact    : 1 if the exact zero acceleration algorithms is used;
*               0 if the approx zeros acceleration algorithm is used.
***

*     Set defaults
      verbose = 0
      exitstatus = 1
      relerr = 0.0d00
      its = 0


*     SPECIAL CASES {

      if ( ( p .LT. 2.0d00 ) .AND. 
     &     ( p .GT. 1.0d00 ) ) then

         if ( y .LT. 0.0d00 ) then

            funvalue = 0.0d00
            return

         elseif ( y .EQ. 0.0d00 ) then

*           If  y=0 in P-G case, then we don't need to integrate, just
*           evaluate; otherwise, integrate

            lambda = calclambda(p, phi, mu)
            funvalue = exp( -lambda )

            return

         endif

      elseif ( p .GE. 2.0d00 ) then

         if ( y .LE. 0.0d00 ) then

            funvalue = 0.0d00
            return

         endif

      endif

*     } SPECIAL CASES
****
****


*     SET ACCURACY REQUIREMENTS
      maxit = 100
      aimrerr = 1.0d-10

*     set other parameters
      iteratn = 0
      ier = 0

      resulta = 0.0d00
      result = 0.0d00
      result0 = 0.5d00

      if ( exact .EQ. 0 ) then
*         Then use approx zero acceleration algorithm

         call evlintc( p, phi, y, mu, aimrerr, result,
     &                 maxit, ier, exitstatus, relerr, 
     &                 its )


         if ( (p.GT. 1.0d00 ) .AND. 
     &        (p .LT. 2.0d00) ) then

            lambda = calclambda(p, phi, mu)

*            Now convert from conditional density to ordinary density
            funvalue = ( ( result + result0 )
     &               * ( 1.0d00 - exp( -lambda ) ) )
     &                 + exp( -lambda )


         elseif ( p .GT. 2.0d00 ) then

            funvalue = result + result0

          endif

      else
*         Then use exact zeros acceleration algorithm

         if ( (p .GT. 1.0d00).AND.(p .LT. 2.0d00) ) then

*           Poisson-gamma case
            call cumsmallp( p, phi, y, mu, aimrerr, 
     &      resulta, ier, relerr, its, verbose )

            result0 = 0.5d00

            lambda = calclambda(p, phi, mu)

*           Now convert from conditional density to ordinary density
            funvalue = ( ( resulta + result0 )
     &              * ( 1.0d00 - exp( -lambda ) ) )
     &            + exp( -lambda )

         elseif ( p .GT. 2.0d00) then

             call cumbigp( p, phi, y, mu, aimrerr, 
     &               resulta, maxit, ier, exitstatus, 
     &               relerr, its, verbose )
             result0 = 0.5d00

             funvalue = resulta + result0

         endif

*        Some tidy-ups
         if ( funvalue .GT. 1.0d00 ) funvalue = 1.0d00

         if ( p .GT. 2.0d00 ) then

            if ( funvalue .LT. 0.0d00 ) then
               funvalue = 0.0d00
            endif

         else

            if ( funvalue .LT. exp(-lambda) ) then

               funvalue = exp(-lambda)

            endif

         endif


      endif

      return
      end

******************************************************************
****************************************************************

****************************************************************
****************************************************************

      subroutine evlint( p, phi, y, mu, aimrerr, result,
     &              maxit, ier, exitstatus, relerr, its, verbose )

***
*     Calculates the density in the case of distributions.
***

      double precision  p, phi, y,  pi, area, aimrerr,
     &                  relerr, result, zero1, zero2,
     &                  mmatrix(2, 200), nmatrix(2, 200),
     &                  xvec(200), w, wold(3), area0, 
     &                  mu, f, tmax, kmax, sumarea
      integer  its, ier, maxit, flag, exitstatus, itsidi, 
     &         mmax, verbose

      external  f, f2

***
*     VARIABLES:
*     exitstatus:    1  if relative error is smaller than wished (aimrerr)
*                   -1  if not, but the absolute error is less than aimrerr
*                  -10  if neither rel or abs error any good
***
*     SET OTHER PARAMETERS
      pi = acos( -1.0d00 )
*      pi = atan( 1.0d00 ) * 4.0d00
*      pi = 3.14159 26535 89793 23846 26433d00
      area = 0.0d00
      area0 = 0.0d00
      its = 0
      relerr = 1.0d00
      flag = 0
      itsidi = 0

      wold(1) = 0.0d00
      wold(2) = 0.0d00
      wold(3) = 0.0d00


      if ( p .GT. 2.0d00 ) then

         kmax = 0.0d00
         tmax = 0.0d00
         mmax = -1

         if ( y .LT. 1.0d00 ) then
*           Find kmax
            call fndkmax(p, phi, y, kmax, tmax, mmax, ier)
         endif


         zero2 = 0.0d00

*        Now, while t < tmax, do not use the Sidi acceleration.
*        Instead, use ...

 400     if ( zero2 .LE. tmax ) then

*           get next zeros: jump the zeros by pi/(y)
            zero1 = zero2
            zero2 = zero2 + ( pi / y )

            call gaussq( f, result, zero1, zero2, 
     &                   p, phi, y, mu )
            area0 = area0 + result
            its = its + 1

            goto 400
            
         endif
         
         xvec(1) = zero2

      else

*        The initial region is done separately
*        since it isn't actually between zeros.
         zero1 = 0.0d00
         zero2 = pi / ( 2.0d00 * y )
         xvec(1) = zero2

         call gaussq( f2, area0, zero1, zero2, 
     &                p, phi, y, mu )
         its = its + 1

      endif

******Works OK to here

*     Now for Sidi acceleration!

*     Now do some more integrations and use Sidi acceleration
  500 if (    ( ( itsidi .LT. 10 ) .AND.
     &          ( flag .NE. 1 )
     &        )
     &     .OR.
     &        ( ( itsidi .LT. maxit ) .AND.
     &          ( flag .NE. 1 ) .AND.
     &          ( abs(relerr) .GT. aimrerr )
     &        )
     &   ) then

*        get next zeros: jump the zeros by pi/(y)
         zero1 = zero2
         zero2 = zero2 + ( pi / y )

*        integrate between zeros
         if ( p .GT. 2.0d00 ) then
            call gaussq( f,  result, zero1, zero2, 
     &                   p, phi, y, mu )
         else
            call gaussq( f2, result, zero1, zero2, 
     &                   p, phi, y, mu )
         endif
         
*        Update interation count         
         its = its + 1
         itsidi = itsidi + 1

*        accelerate convergence of infinite sequence
         xvec( itsidi+1 ) = zero2

         call sidiacc( area, result, xvec, mmatrix,
     &        nmatrix, w, itsidi, relerr, wold, 
     &        sumarea, flag, verbose )

         relerr = ( abs( w-wold(1) ) + 
     &               abs( ( w-wold(2) ) ) ) / w 

         area = area + result

         go to 500

      endif

      result = ( area0 + w ) / pi
      if ( result .LT. 0.0d00 ) then
         result = 0.0d00
      endif
*     occasionally, a result may be very small, but negative

***
**
**     We have integrated to find  \int_0^{\infty}.  The integral
**     required actually goes from -infty to +infty, but is
**     symmetric about the y-axis, so the integral is _twice_ the result
**     obtained above.
**
***

      if ( flag .EQ. 1 ) then
         ier = -10
      endif

*     Determine the error
*     (Keep in this order so the most important aspect is returned)
*     Note also that if we don't use Sidi, and w=0 as set initially,
*     we have good relative error, so that line should be OK.
      if ( abs(w) .LT. aimrerr ) then
         exitstatus = -1
      else
         exitstatus = -10
      endif
      if ( abs(relerr) .LT. aimrerr ) exitstatus = 1

      return
      end

****************************************************************
****************************************************************

      double precision function f(p, phi, y, x)

***
*     A function to be numerically integrated
***

      double precision  x, p, phi, y, rl, im

***
*
*   Tweedie EDM's are characterised by variance functions of the
*   form
*      V( mu ) = phi * mu^p
*   for some p, where  p  is _outside_ the range (0, 1).
*
*   The function to be integrated is a Tweedie EDM of the form
*
*      f = b(y, phi) * exp{ (1/phi) * (y*(mu^[1-p]-1)/(1-p) -
*                                      (mu^[2-p])) /(2-p)}         (1)
*
*   which is of the form
*
*      f = a(y,phi) exp{ (1/phi) * (y * theta - kappa(theta) ) }  (2)
*
*   The density cannot be written in closed form when 1<p<2
*   due to the function outside the exponential
*   (ie, b(y, phi) in (1), or a(y,phi) in (2)),
*   and so there are two options:  (i) Evaluate this function by
*   an infinite summation; or (ii) evaluate by an infinitie integral, as
*   here.  This approach is prefered, since even when the density
*   can be written in closed form, the numerical evaluation
*   is still valid.
*
*   Given the cgf, K(.) for a distribution, the distribution can
*   be reconstructed:
*
*      f_Y(y) = 1/2*pi  \int_{\infty}^{\infty} exp( K(ix) - ixy ) dx
*
*   This function then defines the integrand  exp( K(ix) - ixy )
*   in the above expression.
*
***
*
* VARIABLES:
*   p          : the index in the variance function, V(mu) = phi * mu^p
*   phi        : the dispersion parameter
*   y          : the point at which the density function is to be evaluated
*   mu         : the mean value
*   x          : the internal variable over which to integrate in this function
*   rl         : the real part of the cgf
*   im         : the imaginary part of the cgf
*
***

      call calccgf(p, phi, y, x, rl, im)
      f = dexp( rl ) * cos( im )

      return
      end

*******************************************************************
*******************************************************************

      double precision function f2( p, phi, y, mu, x )

***
*     A function to be numerically integrated in density, 1<p<2
*     using the conditional density.  Note that
*        lim (n->infty) exp( Re k) = exp( -lambda).
***
*     IN:  p, phi, y, mu, x
*     OUT: f2
***
*
* MAJOR VARIABLES:
*   p          : the index in the variance function, V(mu) = phi * mu^p
*   phi        : the dispersion parameter
*   y          : the point at which the density function is to be evaluated
*   mu         : the mean value
*   x          : the internal variable over which to integrate in this function
*   rl         : the real part of the cgf
*   im         : the imaginary part of the cgf
*
***
      double precision  x, p, phi, y, rl, im, mu,
     &                  calclambda, lambda

***

      lambda = calclambda( p, phi, mu )

      if ( x .EQ. 0.0d00 ) then

         f2 = 1.0d00
*        True in the limit

      else

         call calccgf(p, phi, y, x, rl, im)
         f2 = exp( rl ) * cos( im ) - 
     &         exp( -lambda ) * cos( x*y )

      endif

      return
      end

*******************************************************************
*******************************************************************

      subroutine calccgf(p, phi, y, x, rl, im)

***
*     Determines the cgf
***

      double precision  x, p, phi, y, pi, rl, im,
     &                  psi, alpha, front, denom

***
*
*   Given the cgf, K(.) for a distribution, the distribution can
*   be reconstructed thus:
*
*      f_Y(y) = 1/(2\pi)  \int_{\infty}^{\infty} exp( K(ix) - ixy ) dx
*
*   This function determines  K(ix)  in the above expression
*   in terms of the  kappa  and  theta.
*
***
*
* VARIABLES:
*   rl         : the real part of the cgf
*   im         : the imagniary part of the cgf
*   eps        : small value for comparing doubles
*   cgf        : K(s) = (kappa(theta+s*phi) - kappa(theta)) / phi,
*                the cumulant generating function
*   p          : the index in the variance function, V(mu) = phi * mu^p
*   phi        : the dispersion parameter
*   mu         : the mean value
*   lambda     : P(Y=0) = exp( -lambda )
*   x          : the internal variable over which to integrate in this function
*   y          : the point at which the density function is to be evaluated
*   theta      : <see above>
*   calclambda : a function for calculating  lambda
*   kappa      : <see above>
*   t_itphi    : the expression  (theta + i*t*phi)
*   eye        : i (sqrt(-1))
*   one        : 1 + 0i
*
***

*     Set-up some parameters
      pi = acos( -1.0d00 )
*      pi = atan( 1.0d00 ) * 4.0d00
*      pi = 3.14159 26535 89793 32846 26433d+00

*     Now, in the case 1<p<2, we need to alter the mgf as we have to
*     use the conditional distribution of y given y>0.
*     We know that, if  M  is the mgf of the distribution itself,
*     then the mgf of the conditional distribution, M*, is
*     given by
*
*                 M(t) - exp( -lambda )
*        M*(t) = ----------------------
*                   1 - exp( -lambda )
*
*     In terms of cgf's, then, we write as
*
*        K*(t) = log M* = log( exp( cgf ) - exp( -lambda ) )   -
*                         log( 1.0 - exp( -lambda ) )
*

      psi = atan( ( 1.0d00 - p ) * x * phi )
  
      alpha = ( 2.0d00 - p ) / ( 1.0d00 - p )
      front = 1.0d00 / ( phi * ( 2.0d00 - p ) )
      denom = cos( psi ) ** alpha
*        NOTE:  cos(psi) > 0
*      denom = exp( alpha * log(cos(psi)) )

      im = front * sin ( psi * alpha ) / denom - x*y
      rl = front * cos ( psi * alpha ) / denom - front

      return
      end

*******************************************************************
*******************************************************************

      double precision function calclambda(p, phi, mu )

***
*     Calculates lambda, such that  P(Y=0) = exp( -lambda )
*     in the case  1<p<2.
*     NOTE:  IF p>2, this code is not used
***

      double precision  p,  phi, mu

***

      calclambda = mu**( 2.0d00-p ) / 
     &               ( phi * (2.0d00 - p) )

      return
      end

*******************************************************************
*******************************************************************

      subroutine sidiacc( FF, psi, xvec, mmatrix, nmatrix, w,
     &                    znum, relerr, wold, sumarea, 
     &                    flag, verbose )

***
*     Accelerates the series using Sidi's (1988) method.  Works
*     faster than Shanks' method for the series with which
*     we are concerned.
***

      double precision  mmatrix(2, 200), nmatrix(2, 200),
     &          xvec(200), FF, psi, w, relerr, wold(3),
     &          denom, sumarea, abserr, largest, wsave
      integer  i, znum, ell, flag, verbose

***
*     VARIABLES
*     FF      : The cdf up to x_l
*     psi     : the integral between x_l and x_{l+1}
*     xvec    : contains the x-values used up to  x_{l+1}
*     znum    : the iteration number for the table
*     sumarea : the sum of the areas; rarely, psi = 0.0, which
*               causes things to die.  In this case, use sumarea.
*     flag    : if FLAG=1,  then the limits of the machine are being
*               reached (largest integer).
*     wsave   : the old value of w in case something goes awry
*     w       : the value of the series summation, since the first
*               value of x used.
***

      ell = znum - 1

      flag = 0
      largest = 1.0d30

      if ( abs(psi) .LT. 1.0d-31 ) then

         w = FF
         relerr = 0.0d00

         return

      else

         mmatrix(2, 1) = FF / psi
         nmatrix(2, 1) = 1.0d00 / psi
         sumarea = sumarea + psi
*         call intpr("In sidiacc, verbose = ", -1, verbose, 1)
         
         if ( verbose .EQ. 1 ) then
            call dblepr("    w(x) = ", -1, psi, 1)
            call dblepr("    F(x) = ", -1, FF, 1)
            call dblepr("    M-matrix (2,1) = ", -1, 
     &                  mmatrix(2,1), 1)
            call dblepr("    N-matrix (2,1) = ", -1, 
     &                  nmatrix(2,1), 1)
         endif

*        Add the new information
         flag = 0
         do i = 2, znum
   
            if ( verbose .EQ. 1 ) then
               call intpr("    Adding new info at element ", 
     &                     -1, i, 1)
            endif
            
            denom = 1.0d00 / xvec(znum+1-i) - 
     &                        1.0d00 / xvec(znum)
   
            mmatrix(2, i) = ( mmatrix(1, i-1) - 
     &                        mmatrix(2, i-1) ) 
     &                      / denom
            nmatrix(2, i) = ( nmatrix(1, i-1) - 
     &                        nmatrix(2, i-1) ) 
     &                      / denom

            
            if ( verbose .EQ. 1 ) then
               call dblepr("    demoninator = ", -1, denom, 1)
               call dblepr("    New M-matrix entry = ", -1, 
     &                     mmatrix(2,i), 1)
               call dblepr("    New N-matrix entry = ", -1, 
     &                     nmatrix(2,i), 1)
            endif
            
            if ( (abs(mmatrix(2, i)) .GT. largest) .OR.
     &           (abs(nmatrix(2, i)) .GT. largest) ) then

               flag = 1

            endif

         enddo

         
         
         
         if ( (abs(mmatrix(2,znum)) .GT. largest) .OR.
     &        (abs(nmatrix(2,znum)) .GT. largest) ) then

            flag = 1

            wsave = w

         else

            if ( znum .GT. 1 ) then
               w = mmatrix(2, znum) / nmatrix(2, znum)
               if ( verbose .EQ. 1 ) then
                  call dblepr("    New W value = ", -1, w, 1)
               endif
            endif

            wold(1) = wold(2)
            wold(2) = wold(3)
            wold(3)= w

         endif


*        This is the error used in DQEXT:
*         relerr = abs( w - wold(1) ) + abs( ( w - wold(2) ) )
*        This causes problems here though, and failure to work?????
*        (eg, try it with p=1.5, mu=phi=1, y=0(100)10)?????????????

         if ( ell .GT. 3 ) then
            relerr = abs( w - wold(1) ) + 
     &               abs( ( w - wold(2) ) ) / w
            abserr = abs( wold(3) - wold(2) )
            if ( verbose .EQ. 1 ) then
               call dblepr("    Rel. error estimate = ", 
     &                     -1, relerr, 1)
            endif
         else
            relerr = 1.0d00
         endif

         if ( w .NE. 0.0d00 ) then
            relerr = relerr
         endif

*        Now drop first row
         do  i = 1, znum

            mmatrix(1, i) = mmatrix(2, i)
            nmatrix(1, i) = nmatrix(2, i)

         enddo

      endif

      return
      end
  
*******************************************************************
****************************************************************
  
      subroutine fndkmax(p, phi, y, kmax, tmax, mmax, ier)
***
*     Finds  k_max and t_max
*     We solve k'(t)=0 using Newton's method.
***

      double precision  kmax, tmax, imgdcgf, imgddcgf, 
     &            imgcgf, p, largest, dpmmax, zhi, zlo, 
     &            phi, y, pi, psi, front, inner, z,flo, 
     &            fhi, dz1, dz2, z1, z2, sfzro, dk, dz
      integer  mmax, ier, allok, lrgint

      external  dk, imgddcgf, imgdcgf

***

      pi = acos( -1.0d00 )
*      pi = atan( 1.0d00 ) * 4.0d00
*      pi = 3.14159 26535 89793 23846 26433d00
      ier = 0
      allok = 1
      largest = 1.0d30
      lrgint = 100000000

*     A good starting point is sometime crucial, so we spend a 
*     little time finding a decent one.

*     First, the point of inflection, which will always
*     work, but may be very slow.

      psi = ( pi / 2.0d00 ) * ( 1.0d00 - p ) / 
     &      ( 2.0d00*p - 1.0d00 )

*     psi is point of inflection, so start Newton's method at this point
      z2 = 1.0d00 / ( phi * (1.0d00 - p) ) * tan( psi )
      
      dz2 = imgdcgf(p, phi, z2) - y

*     But an alternative, closer to the real value for small y 
*     can be found provided p>2:
      if ( p .GT. 2.0d00) then
         front = -1.0d00 / ( phi * (1.0d00 - p ) )
         inner = ( 1.0d00 / y ) * 
     &             cos( -pi / (2.0d00*(1.0d00-p)) )
*        NOTE:   inner  will always be positive when p>2

         z1 = front * ( inner ) ** ( p - 1.0d00 )
      
         dz1 = imgdcgf(p, phi, z1) - y

*     Now choose the starting point that is the larger, but is still
*     to the left of the k_max; this can be done by examining the
*     derivative at each of the two points.  
*     In effect, we choose the largest starting value for which
*     the derivative is positive.

         if ( dz1 .GT. 0.0d00 ) then
*        If both starting points are to the left of kmax,
*        choose the larger
         
            if ( z1 .GT. z2 ) then
               z = z1
               dz = dz1
            else
               z = z2
               dz = dz2
            endif
         
         else
      
            if ( dz2 .LT. 0.0d00 ) then
*           That is, both to the right.
*           So, use the smaller
            
               if (z1 .GT. z2 ) then
                  z = z2
                  dz = dz2
               else
                  z = z1
                  dz = dz1
               endif
            
            else

*           Use inflection starting point if the alternative is to 
*           the right of kmax
               z = z2
               dz = dz2
         
            endif
               
         endif
      else
*        The case 1 < p < 2
         z = z2
         dz = dz2
            
      endif

*     Now solve for kmax
*     Need bounds on either side
      
      if ( dz .GT. 0.0d00 ) then
*        z  is chosen to be on the left of  kmax
         zlo = z
         zhi = z + 10.d00
         
         flo = dk(p, phi, y, zlo )
         fhi = dk(p, phi, y, zhi )
         allok = 1
      
 565     if ( ( allok .EQ. 1 ) .AND.
     &        ( fhi .GT. 0.d00 ) ) then

            zlo = zhi
            zhi = ( 1.1d00 * zhi ) + 1.0d00

            flo = fhi
            fhi = dk(p, phi, y, zhi )

            if ( zhi .GT. largest/10.d00 ) allok = 0
*              This may not be a problem; it may be that
*              it is so far out no-one really cares!  So
*              let it go, and if problems emerge, they will
*              (hopefully!) be picked up elsewhere.

            goto 565
   
         endif

      else
         zlo = z / 2.0d00
*        As z is chosen to be on the right of  z

         zhi = z
         
         flo = dk(p, phi, y, zlo )
         fhi = dk(p, phi, y, zhi )
         allok = 1
      
 566     if ( ( allok .EQ. 1 ) .AND.
     &        ( flo .LT. 0.d00 ) ) then

            zhi = zlo
            zlo = zlo / 2.0d00

            fhi = flo
            flo = dk(p, phi, y, zlo )


            goto 566
   
         endif

      endif

*     We have the bounds of  zlo  and  zhi
*     To find a useful starting point now, 
*     use linear interpolation:
      z = zlo - flo * (zhi - zlo ) / ( fhi - flo )

*     At last, find the value of tmax (ie. z) and co
      z = sfzro( p, phi, y, zlo, zhi, z,
     &           dk, imgddcgf, ier )

      tmax = z
      kmax = imgcgf( p, phi, y, tmax )
      
      if ( kmax .LT. 0 ) then

         kmax = abs(kmax)
         mmax = lrgint
         allok = 0

      else
      
         dpmmax = ( kmax / pi ) - 0.5d00
         
         if ( dpmmax .GT. dble( lrgint ) ) then
            mmax = lrgint
            allok = 0
         else
            mmax = int( dpmmax )
         endif

      endif

*     If no convergence for tmax, we can probably take it if it is 
*     large enough.  In general, say, the maximum  m  that is used is 5000,
*     when an error is issued.  So if we can make do with m<5000,
*     take the value as being good enough and continue.

      ier = 0

      return
      end

****************************************************************
*******************************************************************

      subroutine gaussq( f1, sum, a, b, p, phi, y, mu )

***
*     Integrates using the Gauss quadrature rules.
*     Currently, the 512-pt rule is implemented.
*
*     The integral requires access to the externally defined
*     function  f  that is of the form
*
*        int_a^b   f(x)  dx
*
*     The integral is transformed to the region [-1 +1].
*
***
*
* INPUTS:
*    f1     :The externally declared, double precision function
*            to be integerated
*    a, b   :The lower and upper integration limits
*
* OUTPUTS:
*    sum    :The value of the integral
*
***

      double precision  weights(256), absc(256),
     &       sum, f1, a, b, xl, xu, p, phi, y, mu
      integer  i, npoints

      external  f1

***
*
* VARIABLES:
*
*    f1       :must be a double precision function declared externally
*    weights  :The weights
*    absc     :The abscissae
*    sum      :The value of the integral
*    a, b     :The lower and upper integration limits
*    x        :The abscissae in the original variable scale
*    npoints  :Half the number of points in the Gaussian integration
*
***
*     The weights and abscicca are  generated using  GRULE
*     The numbers are stoted in  bp.txt  and wf.txt.
***
*NOTE:  There are a limit on the number of continuation lines
*      allowed in FORTRAN; using the data statement, this could easily
*      be exceeded.  Reading in a file of the data also is dangerous;
*      S-Plus with the WATCOM compiler seems to baulk at this.
*      Hence, this is the best (only?) option.


      absc(  1) =   0.003064962185159399599837515282
      absc(  2) =   0.009194771386432905660446301965
      absc(  3) =   0.015324235084898182521206955187
      absc(  4) =   0.021453122959774875710969865850
      absc(  5) =   0.027581204711919792005314633343
      absc(  6) =   0.033708250072480593073631638390
      absc(  7) =   0.039834028811548446991075422829
      absc(  8) =   0.045958310746809061253514983036
      absc(  9) =   0.052080865752192069539905361353
      absc( 10) =   0.058201463766518225784185602834
      absc( 11) =   0.064319874802144239023249383536
      absc( 12) =   0.070435868953604666153900382142
      absc( 13) =   0.076549216406251049948927800415
      absc( 14) =   0.082659687444887164353701791697
      absc( 15) =   0.088767052462401033197103572547
      absc( 16) =   0.094871081968392542704826553290
      absc( 17) =   0.100971546597796779654032661711
      absc( 18) =   0.107068217119502664957941817647
      absc( 19) =   0.113160864444966535735659363127
      absc( 20) =   0.119249259636820398311485291742
      absc( 21) =   0.125333173917474477443434466295
      absc( 22) =   0.131412378677713714836272629327
      absc( 23) =   0.137486645485288105916765744041
      absc( 24) =   0.143555746093496028326086388915
      absc( 25) =   0.149619452449761269896555404557
      absc( 26) =   0.155677536704201868733576930026
      absc( 27) =   0.161729771218192097670396378817
      absc( 28) =   0.167775928572916122050173726166
      absc( 29) =   0.173815781577913441857674570201
      absc( 30) =   0.179849103279615923911549657532
      absc( 31) =   0.185875666969875702472236866925
      absc( 32) =   0.191895246194484031532212497950
      absc( 33) =   0.197907614761680478165928320777
      absc( 34) =   0.203912546750652373672707540209
      absc( 35) =   0.209909816520023939645511745766
      absc( 36) =   0.215899198716335033454427616562
      absc( 37) =   0.221880468282509041300087915261
      absc( 38) =   0.227853400466309585770119383596
      absc( 39) =   0.233817770828785853609588230029
      absc( 40) =   0.239773355252706182882960206371
      absc( 41) =   0.245719929950979243393760498293
      absc( 42) =   0.251657271475063337717870126653
      absc( 43) =   0.257585156723362629360707387605
      absc( 44) =   0.263503362949610298038294331491
      absc( 45) =   0.269411667771238594326632664888
      absc( 46) =   0.275309849177735044278847453825
      absc( 47) =   0.281197685538984665232220550024
      absc( 48) =   0.287074955613597970760508815147
      absc( 49) =   0.292941438557224431704639755480
      absc( 50) =   0.298796913930850727147969791986
      absc( 51) =   0.304641161709084229425315015760
      absc( 52) =   0.310473962288420446409276109989
      absc( 53) =   0.316295096495494865163067288449
      absc( 54) =   0.322104345595318808381790631756
      absc( 55) =   0.327901491299498415443736121233
      absc( 56) =   0.333686315774437136649765989205
      absc( 57) =   0.339458601649521019005817379366
      absc( 58) =   0.345218132025286672526220854706
      absc( 59) =   0.350964690481571417457473671675
      absc( 60) =   0.356698061085645612422467820579
      absc( 61) =   0.362418028400326441840206825873
      absc( 62) =   0.368124377492073107109860075070
      absc( 63) =   0.373816893939063366048003445030
      absc( 64) =   0.379495363839250532400626525487
      absc( 65) =   0.385159573818401101963360133595
      absc( 66) =   0.390809311038112505709563038181
      absc( 67) =   0.396444363203810545837058043617
      absc( 68) =   0.402064518572726958822727283405
      absc( 69) =   0.407669565961855551172732248233
      absc( 70) =   0.413259294755887629513324554864
      absc( 71) =   0.418833494915126280933037605791
      absc( 72) =   0.424391956983378670908990670796
      absc( 73) =   0.429934472095826580861910315434
      absc( 74) =   0.435460831986874741250659326397
      absc( 75) =   0.440970828997976571628214514931
      absc( 76) =   0.446464256085437494192547092098
      absc( 77) =   0.451940906828194099986717446882
      absc( 78) =   0.457400575435571277171931114935
      absc( 79) =   0.462843056755014803371750531369
      absc( 80) =   0.468268146279799957198974880157
      absc( 81) =   0.473675640156716648565549121486
      absc( 82) =   0.479065335193728458751394327919
      absc( 83) =   0.484437028867608643345477048570
      absc( 84) =   0.489790519331549878412346288314
      absc( 85) =   0.495125605422748638062557802186
      absc( 86) =   0.500442086669964369960439398710
      absc( 87) =   0.505739763301052192012718933256
      absc( 88) =   0.511018436250469942905283460277
      absc( 89) =   0.516277907166757477064322756632
      absc( 90) =   0.521517978419990813065965085116
      absc( 91) =   0.526738453109207749314180091460
      absc( 92) =   0.531939135069806612321485772554
      absc( 93) =   0.537119828880917804525552128325
      absc( 94) =   0.542280339872746153240257172001
      absc( 95) =   0.547420474133886614254151936620
      absc( 96) =   0.552540038518610221451865527342
      absc( 97) =   0.557638840654121947792987157300
      absc( 98) =   0.562716688947789034358493154286
      absc( 99) =   0.567773392594340675643138638407
      absc(100) =   0.572808761583037395759276932949
      absc(101) =   0.577822606704811114752828871133
      absc(102) =   0.582814739559374461741469986009
      absc(103) =   0.587784972562300778164967596240
      absc(104) =   0.592733118952072146612408687361
      absc(105) =   0.597658992797097665672367838852
      absc(106) =   0.602562409002699417293058559153
      absc(107) =   0.607443183318068347098517278937
      absc(108) =   0.612301132343186949036351052200
      absc(109) =   0.617136073535721196847703140520
      absc(110) =   0.621947825217879390891084767645
      absc(111) =   0.626736206583239252587702594610
      absc(112) =   0.631501037703541601153744977637
      absc(113) =   0.636242139535451722842651633982
      absc(114) =   0.640959333927286434295922390447
      absc(115) =   0.645652443625708949426211802347
      absc(116) =   0.650321292282389107342055467598
      absc(117) =   0.654965704460630293581857586105
      absc(118) =   0.659585505641960279099578201567
      absc(119) =   0.664180522232690528916521088831
      absc(120) =   0.668750581570438429324099161022
      absc(121) =   0.673295511930615209195138959331
      absc(122) =   0.677815142532878778247606987861
      absc(123) =   0.682309303547550927149245580949
      absc(124) =   0.686777826101999111507723227987
      absc(125) =   0.691220542286981598500972268084
      absc(126) =   0.695637285162956975348436117201
      absc(127) =   0.700027888766357242467108790152
      absc(128) =   0.704392188115823825178551942372
      absc(129) =   0.708730019218407059078401744046
      absc(130) =   0.713041219075728482934550811478
      absc(131) =   0.717325625690105272980190420640
      absc(132) =   0.721583078070637928824737628020
      absc(133) =   0.725813416239259323603505436040
      absc(134) =   0.730016481236746561656048015720
      absc(135) =   0.734192115128692979197921886225
      absc(136) =   0.738340161011444173766449239338
      absc(137) =   0.742460463017992289280755358050
      absc(138) =   0.746552866323834107831203255046
      absc(139) =   0.750617217152788063216917180398
      absc(140) =   0.754653362782772507699746711296
      absc(141) =   0.758661151551544898907764036267
      absc(142) =   0.762640432862400241553757496149
      absc(143) =   0.766591057189829894191746006982
      absc(144) =   0.770512876085140518966909439769
      absc(145) =   0.774405742182031731069002944423
      absc(146) =   0.778269509202133780156884768076
      absc(147) =   0.782104031960504153531132942589
      absc(148) =   0.785909166371082990032448378770
      absc(149) =   0.789684769452107193643541904748
      absc(150) =   0.793430699331483024749900323513
      absc(151) =   0.797146815252117502126338877133
      absc(152) =   0.800832977577207061337105642451
      absc(153) =   0.804489047795484579772562483413
      absc(154) =   0.808114888526424324233232709958
      absc(155) =   0.811710363525404265949703130900
      absc(156) =   0.815275337688824874859960800677
      absc(157) =   0.818809677059186835634818635299
      absc(158) =   0.822313248830123577626238784433
      absc(159) =   0.825785921351392504519139947661
      absc(160) =   0.829227564133821259950707371900
      absc(161) =   0.832638047854211360565557242808
      absc(162) =   0.836017244360197420149916069931
      absc(163) =   0.839365026675062742000932303199
      absc(164) =   0.842681269002510502375002943154
      absc(165) =   0.845965846731390636037417607440
      absc(166) =   0.849218636440382645957924978575
      absc(167) =   0.852439515902632671817684695270
      absc(168) =   0.855628364090346482662141625042
      absc(169) =   0.858785061179337283476797892945
      absc(170) =   0.861909488553529001819697441533
      absc(171) =   0.865001528809411501796944321541
      absc(172) =   0.868061065760453942630192614160
      absc(173) =   0.871087984441469842522565159015
      absc(174) =   0.874082171112937289514377425803
      absc(175) =   0.877043513265272300927222204336
      absc(176) =   0.879971899623057107753254513227
      absc(177) =   0.882867220149221032521325014386
      absc(178) =   0.885729366049175403929893946042
      absc(179) =   0.888558229774901842112910799187
      absc(180) =   0.891353705028992693293332649773
      absc(181) =   0.894115686768646500404145172070
      absc(182) =   0.896844071209613846740182907524
      absc(183) =   0.899538755830097902510544827237
      absc(184) =   0.902199639374606787711741162639
      absc(185) =   0.904826621857757973366176429408
      absc(186) =   0.907419604568035498282085882238
      absc(187) =   0.909978490071499224178808162833
      absc(188) =   0.912503182215446018155091678636
      absc(189) =   0.914993586132022862500434712274
      absc(190) =   0.917449608241791114693342024111
      absc(191) =   0.919871156257243582921034885658
      absc(192) =   0.922258139186271863607657905959
      absc(193) =   0.924610467335585606285519588710
      absc(194) =   0.926928052314082817630946919962
      absc(195) =   0.929210807036171093642451523920
      absc(196) =   0.931458645725040335072719699383
      absc(197) =   0.933671483915885391802191861643
      absc(198) =   0.935849238459080523533373252576
      absc(199) =   0.937991827523303123292919281084
      absc(200) =   0.940099170598609368276754594262
      absc(201) =   0.942171188499458911458361853875
      absc(202) =   0.944207803367690501339382080914
      absc(203) =   0.946208938675447974731014255667
      absc(204) =   0.948174519228055068253979698056
      absc(205) =   0.950104471166841935136915253679
      absc(206) =   0.951998721971919814599516485032
      absc(207) =   0.953857200464905963244177655724
      absc(208) =   0.955679836811598848456128507678
      absc(209) =   0.957466562524601938477530893579
      absc(210) =   0.959217310465897199378559889738
      absc(211) =   0.960932014849367743813957076782
      absc(212) =   0.962610611243270297698870763270
      absc(213) =   0.964253036572656041514051139529
      absc(214) =   0.965859229121740714418820061837
      absc(215) =   0.967429128536223759127210541919
      absc(216) =   0.968962675825556618569578404276
      absc(217) =   0.970459813365158741049754098640
      absc(218) =   0.971920484898583625366086380382
      absc(219) =   0.973344635539632463405723683536
      absc(220) =   0.974732211774417045546670124168
      absc(221) =   0.976083161463370263533079196350
      absc(222) =   0.977397433843205876158322098490
      absc(223) =   0.978674979528826316510503602331
      absc(224) =   0.979915750515178207713518077071
      absc(225) =   0.981119700179057141475880143844
      absc(226) =   0.982286783280859610023583172733
      absc(227) =   0.983416955966283978796127485111
      absc(228) =   0.984510175767978390481971473491
      absc(229) =   0.985566401607137931861757351726
      absc(230) =   0.986585593795049176080169672787
      absc(231) =   0.987567714034582877502543851733
      absc(232) =   0.988512725421635041200829618901
      absc(233) =   0.989420592446515700935094628221
      absc(234) =   0.990291280995286848920500233362
      absc(235) =   0.991124758351048074089817419008
      absc(236) =   0.991920993195171463163717362477
      absc(237) =   0.992679955608486541684953863296
      absc(238) =   0.993401617072414810927227790671
      absc(239) =   0.994085950470055879080177874130
      absc(240) =   0.994732930087228184312664325262
      absc(241) =   0.995342531613465753004277303262
      absc(242) =   0.995914732142977210394008125149
      absc(243) =   0.996449510175577368720212234621
      absc(244) =   0.996946845617603827349739731289
      absc(245) =   0.997406719782849782163225427212
      absc(246) =   0.997829115393562893210344100225
      absc(247) =   0.998214016581612795242506308568
      absc(248) =   0.998561408890039747809908021736
      absc(249) =   0.998871279275449386325647083140
      absc(250) =   0.999143616112378230020851788140
      absc(251) =   0.999378409202599238270181558619
      absc(252) =   0.999575649798310816862567662611
      absc(253) =   0.999735330671042699002271092468
      absc(254) =   0.999857446369979419031892575731
      absc(255) =   0.999941994606845629967040167685
      absc(256) =   0.999988990984381875826159102871


      weights(  1) =   0.006129905175405764294893629085
      weights(  2) =   0.006129674838036492517945319491
      weights(  3) =   0.006129214171953068987508395082
      weights(  4) =   0.006128523194465529920493818139
      weights(  5) =   0.006127601931538031489188345091
      weights(  6) =   0.006126450417787949499770494555
      weights(  7) =   0.006125068696484561869830542946
      weights(  8) =   0.006123456819547496918221263229
      weights(  9) =   0.006121614847544605726714639360
      weights( 10) =   0.006119542849689838838467270676
      weights( 11) =   0.006117240903840640703359454733
      weights( 12) =   0.006114709096494903503571372028
      weights( 13) =   0.006111947522787882815242799239
      weights( 14) =   0.006108956286488514790533610466
      weights( 15) =   0.006105735499995448845034218266
      weights( 16) =   0.006102285284333060395856040969
      weights( 17) =   0.006098605769146656953305640769
      weights( 18) =   0.006094697092697685079920599804
      weights( 19) =   0.006090559401858643313876218173
      weights( 20) =   0.006086192852107506767733724473
      weights( 21) =   0.006081597607521639116401335201
      weights( 22) =   0.006076773840772089693706980995
      weights( 23) =   0.006071721733116766488158599913
      weights( 24) =   0.006066441474393663782493923975
      weights( 25) =   0.006060933263013820911091489307
      weights( 26) =   0.006055197305953895908769979428
      weights( 27) =   0.006049233818748141547350094527
      weights( 28) =   0.006043043025480822859341056841
      weights( 29) =   0.006036625158776993613218841972
      weights( 30) =   0.006029980459794644954973907858
      weights( 31) =   0.006023109178214984885113558732
      weights( 32) =   0.006016011572233289327049643447
      weights( 33) =   0.006008687908549399311897154519
      weights( 34) =   0.006001138462357170390293337192
      weights( 35) =   0.005993363517334814559445188564
      weights( 36) =   0.005985363365633674000154673678
      weights( 37) =   0.005977138307867538649653660343
      weights( 38) =   0.005968688653101249068366751516
      weights( 39) =   0.005960014718839098078750904364
      weights( 40) =   0.005951116831012847295523382485
      weights( 41) =   0.005941995323969674266950669050
      weights( 42) =   0.005932650540459486442068648415
      weights( 43) =   0.005923082831621807528565959444
      weights( 44) =   0.005913292556972981305063452595
      weights( 45) =   0.005903280084392509806379134574
      weights( 46) =   0.005893045790109107881504790782
      weights( 47) =   0.005882590058686672750132284904
      weights( 48) =   0.005871913283009922226995946914
      weights( 49) =   0.005861015864269402374231443531
      weights( 50) =   0.005849898211946678167061364206
      weights( 51) =   0.005838560743798747003363569519
      weights( 52) =   0.005827003885842343793022291010
      weights( 53) =   0.005815228072338094258975083051
      weights( 54) =   0.005803233745774116596194414086
      weights( 55) =   0.005791021356849213735928927349
      weights( 56) =   0.005778591364456383931702543322
      weights( 57) =   0.005765944235664991271428370112
      weights( 58) =   0.005753080445703687324787711788
      weights( 59) =   0.005740000477942362212824267687
      weights( 60) =   0.005726704823874017614981912772
      weights( 61) =   0.005713193983096204360550007806
      weights( 62) =   0.005699468463292455683300019587
      weights( 63) =   0.005685528780212967606133567244
      weights( 64) =   0.005671375457655504839782345528
      weights( 65) =   0.005657009027445287531465911712
      weights( 66) =   0.005642430029415474758425208535
      weights( 67) =   0.005627639011386637545031330632
      weights( 68) =   0.005612636529146218002106483169
      weights( 69) =   0.005597423146427572132610706035
      weights( 70) =   0.005581999434888915492813943331
      weights( 71) =   0.005566365974091757977404437696
      weights( 72) =   0.005550523351479155591270409076
      weights( 73) =   0.005534472162353648236332581689
      weights( 74) =   0.005518213009854874839810179310
      weights( 75) =   0.005501746504936822455833489443
      weights( 76) =   0.005485073266345072764971213530
      weights( 77) =   0.005468193920593387644113470003
      weights( 78) =   0.005451109101940144682774125329
      weights( 79) =   0.005433819452364725861859273692
      weights( 80) =   0.005416325621543047544315108155
      weights( 81) =   0.005398628266823572718902113365
      weights( 82) =   0.005380728053202113274344764449
      weights( 83) =   0.005362625653297344377468114374
      weights( 84) =   0.005344321747325165260222856745
      weights( 85) =   0.005325817023073333225657854939
      weights( 86) =   0.005307112175875508715272577120
      weights( 87) =   0.005288207908585203231854876549
      weights( 88) =   0.005269104931549262356427210108
      weights( 89) =   0.005249803962581399939535398147
      weights( 90) =   0.005230305726934937789185386947
      weights( 91) =   0.005210610957275768270746674204
      weights( 92) =   0.005190720393654706284192190680
      weights( 93) =   0.005170634783479783128101736622
      weights( 94) =   0.005150354881487986119514843608
      weights( 95) =   0.005129881449717141328470404460
      weights( 96) =   0.005109215257477111096773292331
      weights( 97) =   0.005088357081320884003905469228
      weights( 98) =   0.005067307705015408093862649963
      weights( 99) =   0.005046067919512328692199787383
      weights(100) =   0.005024638522917936923894988155
      weights(101) =   0.005003020320463469512717313847
      weights(102) =   0.004981214124474673925202505842
      weights(103) =   0.004959220754341320605562692947
      weights(104) =   0.004937041036486572963271068915
      weights(105) =   0.004914675804335549846868502755
      weights(106) =   0.004892125898284490834178050989
      weights(107) =   0.004869392165668924923882521227
      weights(108) =   0.004846475460731644938072726347
      weights(109) =   0.004823376644591045349363955808
      weights(110) =   0.004800096585208414069756432951
      weights(111) =   0.004776636157355448886185911306
      weights(112) =   0.004752996242581400063165197878
      weights(113) =   0.004729177729179922379243450337
      weights(114) =   0.004705181512155669557029291639
      weights(115) =   0.004681008493190663179162047669
      weights(116) =   0.004656659580610528897937072657
      weights(117) =   0.004632135689350181870227451952
      weights(118) =   0.004607437740919578979259529916
      weights(119) =   0.004582566663369058192201155322
      weights(120) =   0.004557523391254390821014652602
      weights(121) =   0.004532308865601840722203696998
      weights(122) =   0.004506924033872596394023624100
      weights(123) =   0.004481369849927314963355939881
      weights(124) =   0.004455647273990272390353784004
      weights(125) =   0.004429757272613176269371315641
      weights(126) =   0.004403700818638965619467029455
      weights(127) =   0.004377478891165111941907728266
      weights(128) =   0.004351092475507059922912311833
      weights(129) =   0.004324542563160944756706083325
      weights(130) =   0.004297830151766558401393858446
      weights(131) =   0.004270956245069621078080945864
      weights(132) =   0.004243921852884336397282449838
      weights(133) =   0.004216727991055227442451780462
      weights(134) =   0.004189375681419113366110718033
      weights(135) =   0.004161865951766540415446282708
      weights(136) =   0.004134199835803465360173358789
      weights(137) =   0.004106378373112003904443767510
      weights(138) =   0.004078402609111729873458962459
      weights(139) =   0.004050273595020140865452518142
      weights(140) =   0.004021992387813410133046154726
      weights(141) =   0.003993560050186265031335608455
      weights(142) =   0.003964977650512617468603338011
      weights(143) =   0.003936246262804875099827750518
      weights(144) =   0.003907366966673973297796695903
      weights(145) =   0.003878340847288519813162999128
      weights(146) =   0.003849168995334298088578650621
      weights(147) =   0.003819852506973011631308256852
      weights(148) =   0.003790392483801296400619529336
      weights(149) =   0.003760790032809268722963080833
      weights(150) =   0.003731046266338848560462082560
      weights(151) =   0.003701162302042068034252375597
      weights(152) =   0.003671139262839005247551771305
      weights(153) =   0.003640978276875699460451984990
      weights(154) =   0.003610680477481611767159863646
      weights(155) =   0.003580247003126992098170910950
      weights(156) =   0.003549678997380481711154676105
      weights(157) =   0.003518977608865717348479718041
      weights(158) =   0.003488143991218293615136358810
      weights(159) =   0.003457179303042459076605874557
      weights(160) =   0.003426084707867667507319442421
      weights(161) =   0.003394861374104690254077665301
      weights(162) =   0.003363510475001769365471782081
      weights(163) =   0.003332033188600574697552092474
      weights(164) =   0.003300430697691895606804557417
      weights(165) =   0.003268704189771169145439788650
      weights(166) =   0.003236854856993954046573414018
      weights(167) =   0.003204883896131119781075513586
      weights(168) =   0.003172792508523682164511825476
      weights(169) =   0.003140581900037948612225413569
      weights(170) =   0.003108253281020026750902651713
      weights(171) =   0.003075807866250348642650491726
      weights(172) =   0.003043246874898155977795521920
      weights(173) =   0.003010571530475542479515782546
      weights(174) =   0.002977783060791477226514345489
      weights(175) =   0.002944882697905871343779793392
      weights(176) =   0.002911871678083000243575373389
      weights(177) =   0.002878751241745274112165953184
      weights(178) =   0.002845522633426492749991743025
      weights(179) =   0.002812187101725080462522043945
      weights(180) =   0.002778745899257360849748943465
      weights(181) =   0.002745200282610227044549633391
      weights(182) =   0.002711551512294096445698787790
      weights(183) =   0.002677800852695407917564152100
      weights(184) =   0.002643949572029331059747070398
      weights(185) =   0.002609998942291816281802141475
      weights(186) =   0.002575950239212147514084039202
      weights(187) =   0.002541804742204615413792012646
      weights(188) =   0.002507563734320777774911004343
      weights(189) =   0.002473228502201043829678006603
      weights(190) =   0.002438800336026468017908142016
      weights(191) =   0.002404280529470125687963033556
      weights(192) =   0.002369670379648584207510353394
      weights(193) =   0.002334971187073220984936616773
      weights(194) =   0.002300184255601201484958684418
      weights(195) =   0.002265310892386644160689801453
      weights(196) =   0.002230352407831378593050519754
      weights(197) =   0.002195310115535779524331694290
      weights(198) =   0.002160185332249389255493410289
      weights(199) =   0.002124979377821471521886609324
      weights(200) =   0.002089693575151347661178480308
      weights(201) =   0.002054329250138731913916112504
      weights(202) =   0.002018887731633921007318166474
      weights(203) =   0.001983370351387804628867650436
      weights(204) =   0.001947778444001947023567211659
      weights(205) =   0.001912113346878267219550173728
      weights(206) =   0.001876376400168962114978210565
      weights(207) =   0.001840568946726033397465194241
      weights(208) =   0.001804692332050859783151852689
      weights(209) =   0.001768747904243635031551473702
      weights(210) =   0.001732737013952766495436530469
      weights(211) =   0.001696661014324110389878130789
      weights(212) =   0.001660521260950091485333879326
      weights(213) =   0.001624319111818701803773290493
      weights(214) =   0.001588055927262728948823333752
      weights(215) =   0.001551733069908424120231238419
      weights(216) =   0.001515351904624347286962282588
      weights(217) =   0.001478913798470224069681044909
      weights(218) =   0.001442420120645365689757144700
      weights(219) =   0.001405872242437510889756513421
      weights(220) =   0.001369271537171097780430373270
      weights(221) =   0.001332619380155811508736896087
      weights(222) =   0.001295917148634917904354013629
      weights(223) =   0.001259166221733550722339245453
      weights(224) =   0.001222367980406950313532199459
      weights(225) =   0.001185523807388666516285380403
      weights(226) =   0.001148635087138642285956025013
      weights(227) =   0.001111703205791432849669497784
      weights(228) =   0.001074729551104118682389176875
      weights(229) =   0.001037715512404510350211173098
      weights(230) =   0.001000662480539097265105907830
      weights(231) =   0.000963571847821184873997268916
      weights(232) =   0.000926445007979156964772471383
      weights(233) =   0.000889283356104514192963517161
      weights(234) =   0.000852088288600481814569209682
      weights(235) =   0.000814861203130772785048485662
      weights(236) =   0.000777603498568675576864406285
      weights(237) =   0.000740316574946985831578993853
      weights(238) =   0.000703001833408759078253291719
      weights(239) =   0.000665660676159933917435396200
      weights(240) =   0.000628294506424452014331505367
      weights(241) =   0.000590904728403224407604077406
      weights(242) =   0.000553492747240409460780796724
      weights(243) =   0.000516059969000759975743530816
      weights(244) =   0.000478607800667961002377692736
      weights(245) =   0.000441137650179552207388433693
      weights(246) =   0.000403650926533314074063502064
      weights(247) =   0.000366149040035628262328842863
      weights(248) =   0.000328633402852310262977353350
      weights(249) =   0.000291105430251488786486113725
      weights(250) =   0.000253566543570602365847976856
      weights(251) =   0.000216018177976967721805670597
      weights(252) =   0.000178461805545972196102716412
      weights(253) =   0.000140899017388190165439576518
      weights(254) =   0.000103331903496931828490348892
      weights(255) =   0.000065765731659236768705603660
      weights(256) =   0.000028252637373961186168999649
      

      sum = 0.0d00
      npoints = 256

      do i = 1, npoints

*        adjust abscissae
         xl = ( b - a ) / 2.0d00 * 
     &         absc(i) + ( b + a ) / 2.0d00
         xu = ( a - b ) / 2.0d00 * 
     &         absc(i) + ( b + a ) / 2.0d00

*        evaluate
         sum = sum + weights(i) * 
     &            ( f1( p, phi, y, mu, xl ) +
     &              f1( p, phi, y, mu, xu )   )

      enddo

      sum = sum * ( b - a ) / 2.0d00
    
      return
      end

*****************************************************
*******************************************************************

      double precision function imgdcgf(p, phi, x)

***
*     Calculates the imag part of the derivative of cgf
***

      double precision  x, top, bottom, alpha, p, phi, 
     &                  psi, logb

***

      psi = atan( (1.0d00 - p) * x * phi )

      alpha = 1.0d00 /( 1.0d00 - p )
      top = cos( psi * alpha )

      logb = alpha  * log( cos( psi ) )
      bottom = exp( logb )
*     Appear to need thi fix to get it to work...?      
*        NOTE:  cos(psi) > 0, so it should work without logs

      imgdcgf = top/bottom

      return
      end

*******************************************************************
*******************************************************************

      double precision function dk(p, phi, y, x)

***
*     Evaluates the derivative of the  k  function
***

      double precision  x, p, phi, y, rl

***

      call calcdcgf( p, phi, y, x, rl, dk )

      return
      end

*******************************************************************
*******************************************************************

      double precision function imgddcgf(p, phi, x)

***
*     Calculates the imag part of the 2nd derivative of cgf
***

      double precision  x, top, bottom, alpha, p, phi, psi

***

      psi = atan( (1.0d00 - p) * x * phi )

      alpha = p / ( 1.0d00 - p )
      top = sin( psi * alpha)
*      bottom = ( cos( psi ) ) ** alpha
      bottom = exp( alpha * log( abs(cos(psi)) ) )


*        NOTE:  cos(psi) > 0

      imgddcgf = -phi * top/bottom

      return
      end

*******************************************************************
*******************************************************************

      double precision function sfzro( p, phi, y, x1, x2, x0,
     &                                 fun, dfun, ier )

***
*     Uses a modified Newton's Method to find a root between
*     x1 and x2 to find kmax
***

      double precision  df, dx, dxold, f, fh, fl, temp, 
     &                  xh, xl, x1, x2, fun, dfun, y, p, 
     &                  phi, x0
      integer  j, maxit, ier

      external  fun, dfun

***

*     SET PARAMETERS
      ier = 0
      maxit = 100
      
      fl = fun(p, phi, y, x1)
      fh = fun(p, phi, y, x2)

      if ( fl .EQ. 0.0d00 ) then
         sfzro = x1
         return
      elseif ( fh .EQ. 0.0d00 ) then
         sfzro = x2
         return
      elseif ( fl .LT. 0.0d00 ) then
         xl = x1
         xh = x2
      else
         xl = x2
         xh = x1
      endif
         
      sfzro = x0
      dxold = abs( x2-x1 )
      dx = dxold

      f = fun( p, phi, y, sfzro )
      df = dfun( p, phi, y, sfzro )

      do j = 1, maxit

         if (( (sfzro-xh)*df-f)*((sfzro-xl)*df-f).GT.0
     &                .OR.
     &        abs(2.0d00*f) .GT. abs( dxold*df)) then

            dxold = dx
            dx = 0.5d00 * ( xh - xl )
            sfzro = xl + dx      
            if ( xl .EQ. sfzro ) return

         else

            dxold = dx
            if ( df .EQ. 0.0d00 ) return
               
            dx = f/df
            temp = sfzro
            sfzro = sfzro - dx
            if ( temp .EQ. sfzro ) return

         endif

!           if ( abs( dx ) .LT. 1.0d-13 ) return

         f = fun( p, phi, y, sfzro )
         df = dfun( p, phi, y, sfzro )

         if ( f .LT. 0.0d00 ) then
            xl = sfzro
         else
            xh = sfzro
         endif

      enddo

      ier = -20

      return
      end
      
      
*******************************************************************
*******************************************************************

      double precision function imgcgf(p, phi, y, x)

***
*     Calculates the imag part of the cgf
***

      double precision  x, p, phi, y, rl

***

      call calccgf( p, phi, y, x, rl, imgcgf )

      return
      end

*******************************************************************
******************************************************************

      subroutine calcdcgf(p, phi, y, x, rl, im)

***
*     Calculates the derivative of the cgf
***
*     IN:  p, phi, y, x
*     OUT: rl, im
***

      double precision  p, phi, y, x, rl, im,
     &                  psi, alpha, denom

***
* MAJOR VARIABLES:
*   rl         : the real part of the derivative of the cgf
*   im         : the imaginary part of the derivative of the cgf
***

      psi = atan( ( 1.0d00 - p ) * x * phi )

      alpha = 1.0d00 / ( 1.0d00 - p )
      denom = cos( psi ) ** alpha
*        NOTE:  cos(psi) > 0

      rl = -( sin( psi * alpha ) / denom )
      im = cos ( psi * alpha ) / denom - y

      return
      end

****************************************************************
*****************************************************************

 
      subroutine smallp( p, phi, y, mu, aimrerr, result,
     &                   maxit, ier, exitstatus, relerr, its, 
     &                   verbose )

***
*     Calculates the density in the case of distributions with 
*     1 < p < 2
***

      double precision  p, phi, y,  pi, area, aimrerr,
     &         relerr, result, zero1, zero2,
     &         f, g, intim, flo, fhi, t0, dk,
     &         mmatrix(2, 200), nmatrix(2, 200),
     &         xvec(500), w, wold(3), area0, 
     &         sumarea, sfzro2, mu, area1, 
     &         sbuffer, sfzro, zerofn, zerodfn, 
     &         lower, upper, tstep,     
     &         zarea0, z1lo, z1hi, zdelta, resultp
      integer  m, iteratn, ier, maxit, flag, numzr, tier,
     &         exitstatus, its, i, go, totalits,
     &         verbose

      external  f, g, intim, dk, sfzro2, zerofn,
     &          zerodfn, sfzro, f2

***
*     VARIABLES:
*     exitstatus:    1  if relative error is smaller than wished (aimrerr)
*                   -1  if not, but the absolute error is less than aimrerr
*                  -10  if neither rel or abs error any good
***   

*     SET OTHER PARAMETERS
      m = -1
      pi = acos( -1.0d00 )
*      pi = atan( 1.0d00 ) * 4.0d00
*      pi = 3.14159 26535 89793 23846 26433d00
      area = 0.0d00
      area0 = 0.0d00
      area1 = 0.0d00
      iteratn = 0
      relerr = 1.0d00
      sbuffer = 0
      flag = 0
      totalits = 1
*     totalits = 1 (not 0) as the first region is done separately; 
*     this is that one region     
      ier = 0
      tier = 0
      resultp = 0.0d00

      wold(1) = 0.0d00
      wold(2) = 0.0d00
      wold(3) = 0.0d00

*     FIND FIRST ZERO
      zero1 = 0.0d00
      
*     Find bounds on the other zero
      call findsp( p, mu, phi, y, lower, upper, flo, fhi )

*     This is linear interpolation between lower and upper:
      t0 = upper - fhi * ( upper - lower ) / ( fhi - flo )

      zero2 = sfzro( p, phi, y, lower, upper, t0, 
     &               zerofn, zerodfn, ier )
      tier = tier + ier

*     FIND FIRST AREA
*     The first region can be very strange.  For care, we use a
*     very high order Gaussian method, and break the region into
*     small pieces and operate on each separately.  Any funny business
*     should then hopefully be cornered.
      numzr = 20
* WAS 20:  Changed 07/Dec/2005
      
      if ( verbose .EQ. 1 ) then
         call dblepr(" Integrating between ",-1,zero1,1)
         call dblepr(" and ",-1,zero2,1)
         call intpr(" using this many sub-regions: ",-1,numzr,1)
      endif
      
      
      zdelta = zero2 / dble( numzr )
      z1lo = 0.0d00
      z1hi = 0.0d00
      do i = 1, numzr
         zarea0 = 0.0d00
         z1lo = z1hi
         z1hi = z1hi + zdelta
         call gaussq(f2, zarea0, z1lo, z1hi, p, phi, y, mu )
         area0 = area0 + zarea0
      enddo
*     So that is one iteration (between t=0 and t=<first zero>: 
*     hence totalits = 1 up to here

      zero1 = zero2
      tstep = zero2 / 2.0d00

*     NOW, DO A FEW MORE ITS FOR SAFETY'S SAKE
*     as the regions can be dodgy when 1<p<2
      go = 1
      flag = 0
      area1 = 0.0d00

550   if ( go .EQ. 1 ) then

         totalits = totalits + 1

*        FIND THE ZERO
*        We don't jump too far ahead of ourselves--especially
*        early on, when some intervals can be pretty dodgy.
         lower = zero1 + tstep*0.05d00
         upper = zero1 + 0.3d00*tstep

         flo = zerofn( p, phi, y, lower )
         fhi = zerofn( p, phi, y, upper )

*        Try harder to bound the zero        
 650     if ( ( flo * fhi ) .GT. 0.0 ) then

            lower = upper
            upper = upper + 0.5d00*tstep

            flo = zerofn(p, phi, y, lower)
            fhi = zerofn(p, phi, y, upper)

            goto 650
         
         endif


         zero2 = sfzro( p, phi, y, lower, upper, t0,
     &                   zerofn, zerodfn, ier )
         tier = tier + ier

*        Keep track of last result, too
         resultp = result

*        INTEGRATE
         call gaussq( f2, result, zero1, zero2, 
     &                p, phi, y, mu )
         if ( verbose .EQ. 1 ) then
            call dblepr(" Integrating between ",-1,zero1,1)
            call dblepr(" and ",-1,zero2,1)
         endif

*        SUM AREA
         area1 = area1 + result

*        PREPARE FOR NEXT ITERATION
         tstep = zero2 - zero1
         zero1 = zero2 
         t0 = zero2 + ( 0.8d00 * tstep )

*        See if we can stop now
         if ( sbuffer .GE. 3 ) then
            go = 0
         endif

         sbuffer = sbuffer + 1

*        Also stop if the areas are getting really small
*         if ( ( abs(result) .LT. aimrerr/10000.d00 ) .AND.
*     &        (abs(resultp) .LT. aimrerr/10000.d00 ) .AND.
*     &        ( its .GT. 10 ) ) then
*            go = 0
*         endif

*         if ( result .EQ. 0.0d00 ) then
*            go=0
*         endif
* Those lines above commented out 15 Sep 2005

         goto 550

      endif

*     NOW, INTEGRATE WITH ACCELERATION
*     We only need to do this if the regions have any significant area
*      if ( ( abs(result) .LT. aimrerr/10000.d00 ) .AND.
*     &     (abs(resultp) .LT. aimrerr/10000.d00 ) .AND.
*     &     ( its .GT. 10 ) ) then
*         go = 0
*         exitflag = 1
*      else
*         go = 1
*      endif
      go = 1
      flag = 0
      its = 0
      area = 0.0d00
      xvec(1) = zero2

 1550 if ( go .EQ. 1 ) then

         its = its + 1
         totalits = totalits + 1
         
*        FIND THE ZERO
*        Our jumps here can be a little more bold, since most of the
*        initial antics should have been sorted out.
         lower = zero1 + 0.05d00*tstep
         upper = zero1 + 0.8d00 *tstep
 
*        Note:  zero1 has been set above

         flo = zerofn( p, phi, y, lower )
         fhi = zerofn( p, phi, y, upper)

 1650    if ( ( flo * fhi ) .GT. 0.0 ) then

            lower = upper
            upper = upper + 0.5d00*tstep

            flo = zerofn(p, phi, y, lower)
            fhi = zerofn(p, phi, y, upper)

            goto 1650

         endif

*        This is linear interpolation:
         t0 = lower - flo * (upper -  lower ) / 
     &            ( fhi - flo )  
         zero2 = sfzro( p, phi, y, lower, upper, t0,
     &                   zerofn, zerodfn, ier )
          tier = tier + ier

            call gaussq( f2, result, zero1, zero2,
     &                     p, phi, y, mu )

*           ACCELERATE CONVERGENCE
            xvec( its + 1 ) = zero2
            call sidiacc( area, result, xvec, mmatrix,
     &           nmatrix, w, its, relerr, wold, sumarea, 
     &           flag, verbose )
            
            if ( its .GE. 3 ) then
               
               relerr =  ( abs(   w-wold(1) ) + 
     &                     abs( ( w-wold(2) ) ) )
     &               / (area0 + area1 + w)
            endif
         
         if ( flag .EQ.1 ) then
*            print *,'Machine limits being reached...'
         endif

*        SUM AREA
         area = area + result

*        PREPARE FOR NEXT ITERATION
         tstep = zero2 - zero1
         zero1 = zero2 
         t0 = zero2 + ( 0.8d00 * tstep )

*        NOTE IF FLAG=1, the limits of the machine are being reached.
*        We stop and report that required accuracy may not be achieved.

         if (     ( its .LT. 3 ) 
     &         .OR.
     &            ( ( its .LT. maxit ) .AND.
     &              ( abs(relerr) .GT. aimrerr )
     &            )
     &      ) then

*           THEN keep going...
            go = 1

         else

*           THEN stop
            go = 0

         endif

*        Now sometimes we get a very small w, whose relative
*        error isn't that small.  Since w is small, we can make
*        do with it anyway.  So we check this case too.

         if ( flag .EQ. 1 ) then

            ier = -90
            tier = tier + ier

         endif

         goto 1550

      endif



      result = area0 + area1 + w
      result = result / pi
      ier = tier
      
*     Now report the total number of iterations.
*     This is stored as  totalits; we used  its  above
*     as it was needed for Sidi acceleration.  Now we
*     revert to  its
      its = totalits

*     Recall that the function we integrate (which one?) has the division by
*     (1.0d00 - exp( -lambda ) ) already built into it.
            
*     Determine the error
*     In this order, the most important aspect is returned.
*     Note also that if we don't use Sidi, and w=0 as set initially,
*     we have good relative error, so that line should be OK.
      if ( ( abs(w-wold(1))+abs(w-wold(2)) ) .LT. aimrerr ) then
*        Absolute error isn't too bad
         exitstatus = -1
      else
*        Even absolute error isn't too good
         exitstatus = -10
      endif
      if ( abs(relerr) .LT. aimrerr ) then
*        Relative error is inside the required accuracy
         exitstatus = 1
      endif

      return

      end

*****************************************************************
*******************************************************************

      subroutine findsp( p, mu, phi, y, lowerb, upperb, 
     &                   flo, fhi )

***
*     Determines a lower and upper bound for the first zero
*     in the density, conditional trick.
***
*     IN:  p, mu, phi, y
*     OUT: lowerb, upperb, flo, fhi
***

      double precision p, mu, phi, y, lowerb, upperb, pi,
     &          calclambda, lambda, t, told, f1, f2,
     &          zerofn, t3, f3, othzero, rl, im,
     &          wt1, wt2, tstep, flo, fhi

      external calclambda

***

      pi = acos( -1.0d00 )
*      pi = atan( 1.0d00 ) * 4.0d00
*      pi = 3.14159 26535 89793 23846 26433d00

      
***SURELY can improve when y is small, when t gets large.  This
*  routine takes a long time to find the zero.

*     We look at two points of interest.  We have the integrand as
*        exp( rek ) * cos( imk ) - exp( -lambda ) * cos(t*y)  
*
*     For possible zeros, we look at the first zero of each of the
*     two parts of the function.

*     FIRST, the first zero of the second part of the function, because
*     that's easy:   exp( -lambda ) * cos(t*y)

***NECESSARY??? (see later also)
      told = 1.0d-02

      t = pi / y
      f2 = zerofn( p, phi, y, t )
*      print *,'First option: ',t, f2

*     SECOND, we now examine the first part of the function for it's
*     first zero.  This is much trickier.

*     We look at the sign of the zeroing function at the turning
*     points and the zeros of the  sin(t*y)  part.

      f1 = zerofn( p, phi, y, told )
      f2 = zerofn( p, phi, y, t )
*      print *,'f1, zerodfn(0) are: ',f1, zerodfn(p,phi,y,mu,0.0d00)
*AAAAAAAAAAAAAAA
***IT looks to me like told is here just to establish the sign of the
*  derivative as we start.  So can we just evaluate  zerodfn() at x=0???

*Well, it appears that  f1  and  zerodfn(..., 0) gives the same sign.
*But I'm not sure it reads any easier.
      call calccgf(p,phi,y,t,rl,im)
      wt1 = exp( rl )

      t3 = othzero( p, phi, y )

*     So our two candidate values so far are   t  and  t3.

      f3 = zerofn( p, phi, y, t3 )

      lambda = calclambda(p, phi, mu)
      wt2 = exp( -lambda )

*     Now we pick the smaller.  It may not be brilliant, but it the only
*     way to guarantee we won't miss the first zero.  The larger one may
*     give a smaller function value, but may be waaaayyy past the first zero.

*      print *,'CANDIDATES:  ', t, t3
      t = min( t, t3 )
      f2 = zerofn( p, phi, y, t )
      tstep = 0.2d00 * t

 100  if ( ( f1 * f2 ) .GT. 0.0d00 ) then

         told = t
         t = told + tstep

         f1 = f2
         f2 = zerofn( p, phi, y, t )
*      print *,'f1*f2 = ',f1*f2

          goto 100

      endif

      lowerb = told
      upperb = t

      flo = f1
      fhi = f2

      return
      end

*******************************************************************
*******************************************************************

      double precision function zerofn( p, phi, y, x )

***
*     This evaluates the function we are trying to find the zeros
*     of in finding the conditional density.
***
*     IN:  p, phi, y, x
*     OUT: zerofn
***

      double precision  p, mu, phi, y, rl, im, calclambda,
     &                  lambda, x

***

*     When calculating the density, we always have mu=1
      mu = 1.0d00
      call calccgf(p, phi, y, x, rl, im)
      
      lambda = calclambda( p, phi, mu )

       zerofn = exp( rl ) * cos( im )  -
     &         exp( -lambda ) * cos( x*y )

*      print *,'p, phi, y, mu, lambda, im, rl, fn1, fn2'
*      print *,p, phi, y, mu, lambda, im, rl, exp( rl ) * sin( im ),
*     &        exp( -lambda ) * sin( x*y ), zerofn

      return
      end

*****************************************************************
*****************************************************************

      double precision function zerodfn( p, phi, y, x )

***
*     This evaluates the derivative of the function we are trying to
*     find the zeros of in finding the conditional density.
***
*     IN:  p, phi, y, mu, x
*     OUT: zerodfn
***

      double precision  p, mu, phi, y, rl, im, calclambda,
     &                  lambda, x, drl, dim

***

      mu = 1.0d00
      call calccgf(p, phi, y, x, rl, im)
      call calcdcgf(p, phi, y, x, drl, dim)

      lambda = calclambda( p, phi, mu )

      zerodfn = exp( rl ) * ( -dim * sin( im ) ) +
     &          exp( rl )* drl * cos( im )+
     &          exp( -lambda ) * y * sin( x*y )
*      print *,'ZERODFN: rl, im, drl, dim',rl,im,drl,dim

      return
      end



*****************************************************************
****************************************************************

      double precision function g(p, phi, y, x)

***
*     A function to be numerically integrated
*     This function is for use with the Gauss-cos quadrature
*     method, y>1.
***

      double precision  x, p, phi, y, rl, im, imkdash,
     &                  imgdcgf


***

      call calccgf(p, phi, y, x, rl, im)
      imkdash = imgdcgf(p, phi, x)
      
! *     Compute theta
!       theta = ( mu ** (1.0d00-p) - 1.0d00 ) / 
!      &         ( 1.0d00 - p ) )
! 
! *     Compute kappa(theta)
!       if ( abs( 2.0d00 - p ) .LT. 1.0d-06 ) 
! *        Use one more term in Taylor series expansion
!          kappa = log( mu ) + 
!      &      (2.0d00 - p) / 2.0d00 * ( ( log(mu) ) ^ 2.0d00 )
!       
!       else
!          mu ** (2.0d00-p) - 1.0d00 ) / ( 2.0d00 - p )
!       endif
! 
!       D = y * theta - kappa

      g = exp( rl ) / imkdash
      g = exp( rl ) / ( imkdash - y )
***Changed 03 Dec 1999...look at the definition of imgdcgf.
*  Of couse, if thsi work, CHANGE THE FUNCTION IMGCGF AND ALL CALLS TO IT

      return
      end

*******************************************************************
*******************************************************************

      double precision function intim( p, phi, y, x, m )

***
*     Computes
*       Im( k(t) ) - pi/2 - m*pi
*     for finding zeros of the imginary part of he integrand for given  m
***

      double precision  pi, x, p, phi, y, im, rl
      integer  m

***
      pi = acos( -1.0d00 )
*      pi = atan( 1.0d00 ) * 4.0d00
*      pi = 3.14159 26535 89793 23846 26433d00

      call calccgf(  p, phi, y, x, rl, im )

      intim = -pi/2.0d00 - dble(m)*pi + im

      return
      end

******************************************************************
*******************************************************************

      double precision function sfzro2( p, phi, y, x1, x2,
     &                       x0, fun, dfun, m, ier )

***
*     Uses a modified Newton's Method to find a root between
*     x1 and x2 to find kmax
***

      double precision  df, dx, dxold, f, fh, fl, temp, 
     &         xh, xl, x1, x2, fun, dfun, y,
     &         p, phi, x0
      integer  j, maxit, m, ier

      external  fun, dfun

***

*     SET PARAMETERS
      maxit = 100
      ier = 0

*         call dblepr("x1:",-1, x1, 1)
*         call dblepr("x2:",-1, x2, 1)
         
      fl = fun(p, phi, y, x1, m)
      fh = fun(p, phi, y, x2, m)

      if ( ( fl .GT. 0.0d00 .AND. fh .GT. 0.0d00 ) 
     &          .OR.
     &     ( fl .LT. 0.0d00 .AND. fh .LT. 0.0d00 ) ) then

*         print *,'Error: Root must be bounded.'
*         print *,'xs are : ',x1,x2
*         print *,'giving :', fl, fh
*         print *,'for m = ',m
         return

      endif

      if ( fl .EQ. 0.0d00 ) then
      
         sfzro2 = x1
*      print *,'POI Return: sfzro2=',sfzro2
         return
      elseif ( fh .EQ. 0.0d00 ) then
         sfzro2 = x2
*      print *,'ASP Return: sfzro2=',sfzro2
         return
      elseif ( fl .LT. 0.0d00 ) then
         xl = x1
         xh = x2
      else
         xl = x2
         xh = x1
      endif
         
      if ( ( x0 .GT. xl ) .AND. ( x0 .LT. xh ) ) then
         sfzro2 = x0
      else
         sfzro2 = ( xl + xh ) / 2.0d00
      endif
      dxold = abs( x2-x1 )
      dx = dxold

      f = fun( p, phi, y, sfzro2, m )
      df = dfun( p, phi, y, sfzro2 )

      do j = 1, maxit

         if ( ( ( sfzro2-xh*df-f)*(sfzro2-xl)*df-f) 
     &          .GT. 0.0d00
     &                .OR.
     &        abs(2.0d00*f) .GT. abs( dxold*df) ) then
*           Then use bisection

            dxold = dx
            dx = 0.5d00 * ( xh - xl )
            sfzro2 = xl + dx
            if ( xl .EQ. sfzro2 ) then
               return
            endif

         else
*           Then use Newton's method

            dxold = dx
            dx = f/df
            temp = sfzro2
            sfzro2 = sfzro2 - dx
            if ( temp .EQ. sfzro2 ) then
               return

            endif

         endif
         
*         call dblepr("sfzro2:",-1, sfzro2, 1)
*         call intpr("Iteration:",-1, j, 1)

         if ( abs( dx ) .LT. 1.0d-11 ) then
            return
         endif

         f = fun( p, phi, y, sfzro2, m )
         df = dfun( p, phi, y, sfzro2 )

         if ( f .LT. 0.0d00 ) then
            xl = sfzro2
         else
            xh = sfzro2
         endif

      enddo


      ier = -30

      return
      end


*******************************************************************
*****************************************************************

      double precision function othzero(p, phi, y )

***
*     Finds the `other zero' in the density integrand
*     in finding a starting point.  (This `other root' is the zero
*     of sin( Im k).)
***
*     IN:  p, phi, y, mu
*     OUT: othzero
***

      double precision  p, phi, y, pi, psi, inflec, tlo,
     &                  thi, t0, intim, dk, tmax, kmax,
     &                  sfzro2, smallest, largest,
     &                  flo, fhi, zstep
      integer  ier, mmax, maxit, m

      external  intim, dk, sfzro2

***

      largest = 1.0d30
      smallest= 1.0d-30

      ier = 0
      maxit = 100
      pi = acos( -1.0d00 )
*      pi = atan( 1.0d00 ) * 4.0d00
*      pi = 3.14159 26535 89793 23846 26433d00

*     This is the point of inflection: a starting point for the zeroing:
      psi = ( pi / 2.0d00 ) * ( 1.0d00-p ) / 
     &          ( 2.0d00*p - 1.0d00 )
      inflec = atan( psi ) / ( ( 1.0d00-p) * phi )

*     First, we establish whether Im k is going up or down
*     initially.  This established what value of  m  we are
*     using to get the zero.  If y>mu, we are heading up first.

      if ( y .GE. 1.0d00 ) then
*        In this case, we initially head down.

         m = -1
         kmax = 0.0d00
         tmax = 0.0d00
         tlo = 1.0d-05
         thi = inflec

      else

*        In this case, we have to find k_max.  To do so, we
*        solve Im k'(t) = 0.
         call fndkmax( p, phi, y, 
     &                 kmax, tmax, mmax, ier )

         if ( kmax .GE. pi/2.0d00 ) then

            m = 0
            tlo = smallest
            thi = tmax

         else

            m = -1 
            tlo = min(tmax, inflec)
            thi = max(tmax, inflec)

         endif

      endif

      flo = intim(p, phi, y, tlo, m)
      fhi = intim(p, phi, y, thi, m)
      zstep = thi - tlo

 565  if ( ( flo * fhi ) .GT. 0.0d00 ) then

           tlo = thi
           thi = thi + 0.2d00*zstep

           flo = intim(p, phi, y, tlo, m)
           fhi = intim(p, phi, y, thi, m)

           goto 565

        endif

*     Having established  m, this means we can find the root
*     we are after a bit more easily.

*     Linear interpolation for the estimate:
      t0 = tlo - flo * ( thi - tlo ) / ( fhi - flo )

      othzero =  sfzro2( p, phi, y, tlo, thi, t0,
     &                intim, dk, m, ier)

      return
      end

*****************************************************************
****************************************************************

      subroutine bigp( p, phi, y, mu, aimrerr, result,
     &                    maxit, ier, exitstatus, relerr, 
     &                    its, verbose )


***
*     Calculates the density in the case of distributions with
*     p>2.
***

      double precision  p, phi, y, mu, pi, area, aimrerr,
     &          relerr, result, zero1, zero2, zero,
     &          f, g, zlo, zhi, intim, flo, fhi, kmax,
     &          tmax, mmatrix(2, 200), nmatrix(2, 200),
     &          xvec(200), w, wold(3), area0, sumarea,
     &          dk, sfzro2, diff, largest, smallest
      integer  m, its, mmax, firstm, ier, maxit, flag,
     &         allok, kmaxok, tier, exitstatus, verbose

      external  f, g, intim, dk, sfzro2

***
*     VARIABLES:
*     exitstatus:    1  if relative error is smaller than wished (aimrerr)
*                   -1  if not, but the absolute error is less than aimrerr
*                  -10  if neither rel or abs error any good
***

      if ( verbose .EQ. 1 ) then
         call dblepr("Using p>2 code since p = ",-1, p, 1)
      endif


*     SET ACCURACY REQUIREMENTS
      largest = 1.0d30
      smallest = 1.0d-30

*     SET OTHER PARAMETERS
      m = -1
      pi = acos( -1.0d00 )
*      pi = atan( 1.0d00 ) * 4.0d00
*      pi = 3.14159 26535 89793 23846 26433d00
      area = 0.0d00
      area0 = 0.0d00
      its = 0
      relerr = 1.0d00
      flag = 0
      tier = 0
      allok = 1

      wold(1) = 0.0d00
      wold(2) = 0.0d00
      wold(3) = 0.0d00



      if ( y .GE. 1.0) then
      
         if ( verbose .EQ. 1 ) then
            call dblepr(" Using y .GE. 1 since y = ",-1, y, 1)
         endif
*        In this case, Im(k) heads down straight away.

*        FIND ZEROS
         zero1 = 0.0d00

*        An approximation to the first zero:
         zero = pi / ( 2.0d00 * y )

*        Bracket first zero
         zlo = 0.9d00 * pi / (2.0d00 * y )
         
         if ( y .GT. 1.0d00 ) then
            zhi = pi / (2.0d00 * ( y - 1.0d00 ) )
            fhi = intim( p, phi, y, zhi, m )
         else
            zhi = zero * 2.0d00  
            fhi = intim( p, phi, y, zhi, m )
         endif
         flo = intim( p, phi, y, zlo, m )
         
         allok = 1
 
 565     if ( ( allok .EQ. 1 ) .AND.
     &        (fhi * flo ) .GT. 0.0d00 ) then

            zlo = zhi 
            zhi = zhi * 1.5d00

            flo = intim( p, phi, y, zlo, m )
            fhi = intim( p, phi, y, zhi, m )

            if ( zhi .GT. largest/10.0d00 ) allok = 0

            goto 565

         endif

         if ( zhi .GT. largest/10.0d00 ) allok = 0

         if ( allok .EQ. 0 ) then
            ier = -50
            tier = tier + ier
            result = 0.0d00
            exitstatus = -10
            return
         endif

         if ( verbose .EQ. 1 ) then
            call dblepr("  First zero bound by ",-1, zlo, 1)
            call dblepr("  and ",-1, zhi, 1)
         endif
         zero2 = sfzro2( p, phi, y, zlo, zhi, zero, 
     &                   intim, dk, m, ier )
         if ( verbose .EQ. 1 ) then
            call dblepr("  Found first zero = ",-1, zero2, 1)
         endif

         xvec( 1 ) = zero2

*        special case: between 0 and first zero
         if ( verbose .EQ. 1 ) then
            call dblepr("   Integrating between y = ",-1, zero1, 1)
            call dblepr("                   and y = ",-1, zero2, 1)
         endif
         call gaussq( f, area0, zero1, zero2, 
     &                p, phi, y, mu )
         if ( verbose .EQ. 1 ) then
            call dblepr("   Area is ",-1, area0, 1)
         endif

*        NOTE:  We don't update the iteration count here, since
*        we keep it all to work with Sidi acceleration; we instead
*        add one later (after accelerating)

*        Now do some more integrations and use Sidi acceleration
!   500    if (    ( ( its .LT. 4 ) .AND.
!      &             ( flag .NE. 1 )
!      &           )
!      &        .OR.
!      &           ( ( its .LT. maxit ) .AND.
!      &             ( flag .NE. 1 ) .AND.
!      &             ( abs(relerr) .GT. aimrerr )
!      &           ) 
!      &      ) then
  500    if (    ( its .LT. 4 )
     &        .OR.
     &           ( ( its .LT. maxit ) .AND.
     &             ( abs(relerr) .GT. aimrerr )
     &           ) 
     &      ) then
     
            if ( verbose .EQ. 1 ) then
               call intpr("   Iterating; iteration ",-1, its, 1)
            endif
     

*           get next zeros
            m = m - 1
            zero1 = zero2
            zero = zero2
            zlo = zero2
            zhi = zero2 * 1.5d00

            if ( zhi .GT. largest/10.0d00 ) then

               allok = 0
               flo = 0.0d00
               fhi = 0.0d00

            else

               allok = 1
               flo = intim( p, phi, y, zlo, m )
               fhi = intim( p, phi, y, zhi, m )
 
            endif

 765        if ( ( allok .EQ. 1 ) .AND.
     &           (fhi * flo ) .GT. 0.0d00 ) then

               zlo = zhi
               zhi = zhi * 1.5d00

               flo = intim( p, phi, y, zlo, m )
               fhi = intim( p, phi, y, zhi, m )

*              Linearly interpolate for zero:
               zero = zhi - fhi * ( zhi - zlo ) / 
     &                  ( fhi - flo )

               if ( zhi .GT. largest/10.0d00 ) allok = 0

               goto 765

            endif

            if ( zhi .GT. largest/10.0d00 ) allok = 0

            if ( allok .EQ. 0 ) then
               ier = -50
               tier = tier + ier
               result = 0.0d00
               exitstatus = -10
               return
            endif
     
            zero2 = sfzro2( p, phi, y, zlo, zhi, zero, 
     &                      intim, dk, m, ier )

            if ( ier .NE. 0 ) then
               tier = tier + ier
            endif


*           integrate between zeros
            if ( verbose .EQ. 1 ) then
               call dblepr("   Integrating between y = ",-1, zero1, 1)
               call dblepr("                   and y = ",-1, zero2, 1)
            endif

            call gaussq( f, result, zero1, zero2, 
     &                   p, phi, y, mu )
            if ( verbose .EQ. 1 ) then
               call dblepr("   Giving area = ",-1, result, 1)
            endif

*           Update iteration count
            its = its + 1

*           accelerate convergence of infinite sequence
            xvec( its+1 ) = zero2

            call sidiacc( area, result, xvec, mmatrix,
     &           nmatrix, w, its, relerr, wold, sumarea, 
     &           flag, verbose )
            if ( verbose .EQ. 1 ) then
               call dblepr("   Accelerating; w = ",-1, w, 1)
            endif

            relerr = (abs( w-wold(1))+abs(( w-wold(2))))
     &               / (area0 + w)       

            area = area + result

            if ( verbose .EQ. 1 ) then
               call dblepr("  Area = ",-1, area, 1)
            endif
            
            go to 500

         endif

         if ( ( its .GE.  maxit ) .AND.
     &        ( abs(relerr) .GT. aimrerr) ) then
            ier = -40
            tier = tier + ier
         endif
         if ( flag .EQ. 1 ) then
            ier = -70
            tier = tier + ier
         endif

         result = area0 + w

*        Now, the very first integration has not been counted
*        (since that can foul up the Sidi acceleration iteration count)
*        so update now
         its = its + 1

      else
*        that is:    if ( y .LT. 1) then
*        In this case, Im(k) may head up before going to  -infinity
         if ( verbose .EQ. 1 ) then
            call dblepr(" Using y .LT. 1 since y = ",-1, y, 1)
         endif


*        FIND k_max AND t_max
         kmaxok = 1
         call fndkmax(p, phi, y, kmax, tmax, mmax, ier)
         if ( verbose .EQ. 1 ) then
            call dblepr(" Found kmax = ",-1, kmax, 1)
         endif

         if ( ier .NE. 0 ) then
            tier = tier + ier
            kmaxok = 0
         endif
         

         if ( kmax .LT. pi/2.0d00 ) then

            if ( verbose .EQ. 1 ) then
               call dblepr(" Using kmax .LT. pi/2; kmax = ",
     &                     -1, kmax, 1)
            endif
            m = -1

*           FIND ZEROS
            zero1 = 0.0d00
            zero = tmax + pi/( 2.0d00*y )

*           BOUNDS ON `OTHER' ZERO:
            zlo = tmax
            zhi = zero*2.0d00

            if ( zhi .GT. largest/10.0d00 ) then

               allok = 0
               flo = 0.0d00
               fhi = 0.0d00

            else

               allok = 1
               flo = intim( p, phi, y, zlo, m )
               fhi = intim( p, phi, y, zhi, m )

            endif

 1565       if ( (allok .EQ. 1 ) .AND.
     &            ( (fhi * flo ) .GT. 0.0d00 ) ) then

               zlo = zhi
               zhi = zhi * 1.5d00

               flo = intim( p, phi, y, zlo, m )
               fhi = intim( p, phi, y, zhi, m )

               if ( zhi .GT. largest/10.0d00 ) allok = 0

               goto 1565

            endif

            if ( zhi .GT. largest/10.0d00 ) allok = 0

            if ( allok .EQ. 0 ) then
               ier = -50
               tier = tier + ier
               result = 0.0d00
               exitstatus = -10
               return
            endif

            zero2 = sfzro2( p, phi, y, zlo, zhi, zero, 
     &                      intim, dk, m, ier )

            if ( ier .NE. 0 ) then
               tier = tier + ier
            endif

            xvec(1) = zero2

*           integrate between zeros
            if ( verbose .EQ. 1 ) then
               call dblepr("  Integrating between y = ",
     &                     -1, zero1, 1)
               call dblepr("                  and y = ",
     &                     -1, zero2, 1)
            endif
            
            call gaussq( f, area0, zero1, zero2, 
     &                   p, phi, y, mu )
            
            if ( verbose .EQ. 1 ) then
               call dblepr("  Giving area = ",-1, area0, 1)
            endif

*        NOTE:  We don't update the iteration count here, since
*        we keep it all to work with Sidi acceleration; we instead
*        add one later (after accelerating)


  600       if (    ( its .LT. 4 )
     &            .OR.
     &               ( ( its .LT. maxit ) .AND.
     &                 ( abs(relerr) .GT. aimrerr )
     &               ) ) then

            if ( verbose .EQ. 1 ) then
               call intpr("  Iteration ",-1, its, 1)
            endif

*              get next zeros
               m = m - 1
               diff = zero2 - zero1
               zero1 = zero2

               zlo = zero2 - 0.01d00*diff
               zhi = zero2 + 2.0d00*diff

               if ( zhi .GT. largest/10.d00 ) then

                  allok = 0
                  flo = 0.0d00
                  fhi = 0.0d00

               else

                  allok = 1
                  flo = intim( p, phi, y, zlo, m )
                  fhi = intim( p, phi, y, zhi, m )

               endif

 1665          if ( ( allok .EQ. 1 ) .AND.
     &              (fhi * flo ) .GT. 0.0d00 ) then

                  zlo = zhi
                  zhi = zhi * 1.5d00

                  flo = intim( p, phi, y, zlo, m )
                  fhi = intim( p, phi, y, zhi, m )

                  if ( zhi .GT. largest/10.0d000 ) then
                     allok = 0
                  endif

                  goto 1665
   
               endif

               if ( zhi .GT. largest/10.d00 ) allok = 0

               if ( allok .EQ. 0 ) then
                  ier = -50
               tier = tier + ier
                  result = 0.0d00
               exitstatus = -10
   
                  return

               endif

*              Approximate zero with linear interpolation
               zero = zlo - flo * ( zhi - zlo ) / 
     &                     ( fhi - flo )

               zero2 = sfzro2( p, phi, y, zlo, zhi, zero,
     &                         intim, dk, m, ier )

               if ( ier .NE. 0 ) then
                  tier = tier + ier
               endif


*              integrate between zeros
               if ( verbose .EQ. 1 ) then
                  call dblepr("  Integrating between y = ",-1, zero1, 1)
                  call dblepr("                  and y = ",-1, zero2, 1)
               endif
               call gaussq( f, result, zero1, zero2, 
     &                      p, phi, y, mu )
               if ( verbose .EQ. 1 ) then
                  call dblepr("  Giving area ",-1, result, 1)
               endif
     
*              Update interation count
               its = its + 1

*              accelerate convergence of infinite sequence
               xvec( its+1 ) = zero2

               call sidiacc( area, result, xvec, mmatrix, 
     &              nmatrix, w, its, relerr, wold, 
     &              sumarea, flag, verbose )
               if ( verbose .EQ. 1 ) then
                  call dblepr("  Accelerating; w = ",-1, w, 1)
               endif
               relerr = ( abs( w-wold(1) ) + 
     &                 abs( ( w-wold(2) ) ) )
     &                  / (area0 + w)       

               area = area + result
               
               if ( verbose .EQ. 1 ) then
                  call dblepr("  Area = ",-1, area, 1)
               endif

               go to 600

            endif

            if ( ( its .GE.  maxit ) .AND.
     &           ( abs(relerr) .GT. aimrerr) ) then
               ier = -40
               tier = tier + ier
            endif

            if ( flag .EQ. 1 ) then
               ier = -70
               tier = tier + ier
            endif
      
            result = area0 + w

*           Now, the very first integration has not been counted
*           (since that can foul up the Sidi acceleration iteration count)
*           so update now
            its = its + 1

         else
         
*           that is:   case where kmax >= pi/2
*           In this case, the upward trend is goes above  pi/2, and so
*           the first zero will be at  m=0.

            if ( verbose .EQ. 1 ) then
               call dblepr(" Using kmax .GE. pi/2; kmax = ",
     &                     -1, kmax, 1)
            endif

*           Now, kmax may not have been found accurately.  IF, however,
*           the corresponding max > maxit, it won't matter a great deal
*           and we can proceed.  If not, accuracy cannot be ensured
*           unless the maximum  m  used is less than maxit.

            if ( ier .EQ. -80 ) then
               tier = tier + ier
               result = 0.0d00
               exitstatus = -10
               return
            endif

            wold(1) = 0.0d00
            wold(2) = 0.0d00
            wold(3) = 0.0d00

*           FIND ZEROS
            zero1 = 0.0d00
            zero = pi / ( 2.0d00*( 1.0d00 - y ) )

            m = 0
            firstm = 1

            zlo = smallest
            zhi = tmax

            if ( zhi .GT. largest/10.d00 ) then

              allok = 0
              flo = 0.0d00
              fhi = 0.0d00

            else

               allok = 1
               flo = intim( p, phi, y, zlo, m )
               fhi = intim( p, phi, y, zhi, m )

            endif

            diff = zhi - zlo

 2565       if ( ( allok .EQ. 1 ) .AND. 
     &           (fhi * flo ) .GT. 0.0d00 ) then

               zlo = zhi
               zhi = zhi + 0.1d00*diff

               flo = intim( p, phi, y, zlo, m )
               fhi = intim( p, phi, y, zhi, m )

               if ( zhi .GT. largest / 10.0d00 ) then
                  allok = 0
               endif

               goto 2565

            endif

            if ( zhi .GT. largest/10.d00 ) allok = 0

            if ( allok .EQ. 0 ) then
               ier = -50
               tier = tier + ier
               result = 0.0d00
               exitstatus = -10

               return

            endif

            zero2 = sfzro2( p, phi, y, zlo, zhi, zero, 
     &                      intim, dk, m, ier )

            if ( ier .NE. 0 ) then
               tier = tier + ier
            endif

            xvec( 1 ) = zero2

            if ( verbose .EQ. 1 ) then
               call dblepr("  Integrating between y = ",-1, zero1, 1)
               call dblepr("                  and y = ",-1, zero2, 1)
            endif
            call gaussq( f, area0, zero1, zero2, 
     &                   p, phi, y, mu )
            if ( verbose .EQ. 1 ) then
               call dblepr("  Giving area = ",-1, area0, 1)
            endif
*           NOTE:  We don't update the iteration count here, since
*           we keep it all to work with Sidi acceleration; we instead
*           add one later (after accelerating)


            diff = zero2 - zero1

  700       if (    ( its .LT. 4 )
     &           .OR.
     &              ( ( its .LT. maxit ) .AND.
     &                ( abs(relerr) .GT. aimrerr )
     &              ) ) then

               if ( verbose .EQ. 1 ) then
                  call intpr("  Iteration ", 1, its, 1)
               endif

*              get next zeros

               zlo = zero2 - 1.0d-05*diff
               zhi = zero2 + 2.0d00*diff

               zero1 = zero2

*              FIND THE NEXT VALUE OF m
               call nextm( tmax, mmax, zero2, m, firstm, 
     &                     zlo, zhi, zero )

               if ( zhi .GT. largest/10.d00 ) then

                  allok = 0
                  flo = 0.0d00
                  fhi = 0.0d00

               else

                  allok = 1
                  flo = intim( p, phi, y, zlo, m )
                  fhi = intim( p, phi, y, zhi, m )

               endif

 2665          if ( ( allok .EQ. 1 ) .AND.
     &              (fhi * flo ) .GT. 0.0d00 ) then

                  zlo = zhi
                  zhi = zhi * 1.5d00

                  flo = intim( p, phi, y, zlo, m )
                  fhi = intim( p, phi, y, zhi, m )

                  if ( zhi .GT. largest/10.0d00 ) then
                     allok = 0
                  endif

                  goto 2665

               endif

               if ( zhi .GT. largest/10.0d00 ) allok = 0

               if ( allok .EQ. 0 ) then
                  ier = -50
               tier = tier + ier
                  result = 0.0d00
                  exitstatus = -10
   
                  return
               endif

               zero2 = sfzro2( p, phi, y, zlo, zhi, zero,
     &                         intim, dk, m, ier )

               if ( ier .NE. 0 ) then
               tier = tier + ier
               endif

               if ( verbose .EQ. 1 ) then
                  call dblepr("  Iteragrating between y = ", 
     &                        1, zero1, 1)
                  call dblepr("                   and y = ", 
     &                        1, zero2, 1)
               endif
               call gaussq( f, result,  zero1, zero2,
     &                      p, phi, y, mu )
               if ( verbose .EQ. 1 ) then
                  call dblepr("  giving area = ", 1, 
     &                        result, 1)
               endif

*              Update iteration count
               its = its + 1
               
*              accelerate convergence of infinite sequence
               xvec( its+1 ) = zero2
               call sidiacc( area, result, xvec, mmatrix, 
     &            nmatrix, w, its, relerr, wold, sumarea, 
     &            flag, verbose )
               if ( verbose .EQ. 1 ) then
                  call dblepr("  Accelerating; w = ", 
     &                        1, w, 1)
               endif

               area = area + result
               
               if ( verbose .EQ. 1 ) then
                  call dblepr("  Area = ", 1, area, 1)
               endif

               go to 700

            endif

            if ( ( its .GE.  maxit ) .AND.
     &           ( abs(relerr) .GT. aimrerr) ) then
               ier = -40
               tier = tier + ier
            endif
            if ( flag .EQ. 1 ) then
               ier = -70
               tier = tier + ier
            endif

            if ( kmaxok .EQ. 0 ) then
*              IF finding kmax didn't converge...

               if ( m .LT. mmax-1 )  kmaxok = 1
*              All should be OK if greatest value of m used
*              is less than the `turning' m value.

            endif

            if ( kmaxok .EQ. 0 )  then

               ier = -60
               tier = tier + ier

            endif

            result = area0 + w
            
*           Now, the very first integration has not been counted
*           (since that can foul up the Sidi acceleration iteration count)
*           so update now
            its = its + 1

         endif

      endif
         
***
**
**     We have integrated to find  \int_0^{\infty}.  The integral
**     required actually goes from -infty to +infty, but is
**     symmetric about the y-axis, so the integral is _twice_ the result
**     obtained above.
**
***


      result = abs( result / pi )
*     occasionally, a result may be very small, but negative

      if ( flag .EQ. 1 ) then
         ier = -10
      endif

      if ( ier .NE. 0 ) then
               tier = tier + ier
      endif


*     Determine the error
*     (Keep in this order so the most important aspect is returned)
*     Note also that if we don't use Sidi, and w=0 as set initially,
*     We have good relative error, so that line should be OK.
      if ( abs(w) .LT. aimrerr ) then
         exitstatus = -1
      else
         exitstatus = -10
      endif
      if ( abs(relerr) .LT. aimrerr ) exitstatus = 1

      if ( verbose .EQ. 1 ) then
         call dblepr("Final result = ", 1, result, 1)
      endif
      
      return
      end

*****************************************************************
******************************************************************

      subroutine  nextm( tmax, mmax, zero2, m, firstm, 
     &                   zlo, zhi, zero )

***
*     This subroutine advances to the next value of  m  in the
*     case  p>2, and  y<mu  and  kmax>pi.  It also returns,
*     when appropriate, adjusted bounds on the next zero.
***
*     IN:    tmax, mmax, zero2, m, firstm, zlo
*     OUT:   m, firstm, zlo, zhi, zero
***

      double precision  tmax, zero2, zlo, zhi, zero
      integer  mmax, m, firstm

***

*     This can get tricky, since we need to be careful when  Im(k)  
*     turns back down.


      if ( m .LT. mmax ) then

         if ( firstm .EQ. 1 ) then
            m = m + 1
            zhi = tmax
         else
            m = m - 1
            zlo = max( zlo, tmax )
         endif

      elseif ( m .EQ. mmax ) then

         if ( firstm .EQ. 1 ) then
            firstm = firstm + 1
            zero = tmax + (tmax - zero2)
            zlo = tmax
         else
            m = m - 1
*            zlo = tmax
         endif

      endif

      return
      end

*******************************************************************

*** The uniquely  cdf  stuff

*****************************************************************

      subroutine evlintc( p, phi, y, mu, aimrerr, result,
     &                       maxit, ier, exitstatus, relerr, its )

***
*     Calculates the density in the case of distributions when
*     p>2
***

      double precision  p, phi, y,  pi, area, aimrerr,
     &          relerr, result, zero1, zero2, cumf,
     &          cumintim, imgddcgf, cumsfzro,
     &          mmatrix(2, 200), nmatrix(2, 200),
     &          xvec(200), w, wold(3), area0, sumarea,
     &          cumdk, cumsfzro2, mu,
     &          largest,smallest,
     &          area1, zerofn, zerodfn,
     &          kmax, tmax
      integer  its, ier, maxit, flag, 
     &         mmax, exitstatus, verbose,
     &         itsidi

      external  cumdk, cumsfzro2, cumintim, imgddcgf,
     &          cumsfzro, cumf,zerofn, zerodfn,
     &          cumf2

***

*     SET ACCURACY REQUIREMENTS
      largest = 1.0d30
      smallest = 1.0d-30

*SET M and N matrices, and other sidi things, to 0?


*     SET OTHER PARAMETERS
      pi = acos( -1.0d00 )
*      pi = atan( 1.0d00 ) * 4.0d00
*      pi = 3.14159 26535 89793 23846 26433d00
      area = 0.0d00
      area0 = 0.0d00
      area1 = 0.0d00
      result = 0.0d00
      its = 0
      itsidi = 0
      relerr = 1.0d00
      flag = 0


      wold(1) = 0.0d00
      wold(2) = 0.0d00
      wold(3) = 0.0d00

      if ( p .GT. 2.0d00 ) then

*        Find kmax
         call cumfndkmax(p, phi, y, mu, 
     &            kmax, tmax, mmax, ier)

*        While t < tmax, do not use Sidi acceleration

         zero2 = 0.0d00

 400     if ( zero2 .LE. tmax ) then

*           get next zero: jump the zero by pi/y
            zero1 = zero2
            zero2 = zero2 + (pi/y)

            call gaussq( cumf, result, zero1, zero2,
     &                   p, phi, y, mu )
            area0 = area0 + result

            its = its + 1

            goto 400

         endif

      else

         zero1 = 0.0d00
         zero2 = pi / y
         xvec(1) = zero2


         call gaussq( cumf2, area0, zero1, zero2,
     &                p, phi, y, mu )
         its = its + 1

      endif

*     Now for Sidi acceleration
*     itsidi = itsidi + 1
      xvec(1) = zero2

*     Now do some more integrations and use sidi acceleration
  500 if (    ( ( itsidi .LT. 4 ) .AND.
     &          ( flag .NE. 1 )
     &        )
     &     .OR.
     &        ( ( itsidi .LT. maxit ) .AND.
     &          ( flag .NE. 1 ) .AND.
     &          ( abs(relerr) .GT. aimrerr )
     &        )
     &   ) then

*        get next zeros: jump by pi/y
         zero1 = zero2
         zero2 = zero2 + ( pi / y)

*        integrate between zeros
         if ( p .GT. 2.0d00 ) then
           call gaussq( cumf, result, zero1, zero2,
     &                  p, phi, y, mu )
         else
           call gaussq( cumf2, result, zero1, zero2,
     &                  p, phi, y, mu )
         endif

*        Update interation count
         its = its + 1
         itsidi = itsidi + 1

*        accelerate convergence of infinite sequence
         xvec( itsidi+1 ) = zero2

         call sidiacc( area, result, xvec, mmatrix,
     &        nmatrix, w, itsidi, relerr, wold, sumarea, 
     &        flag, verbose )

         relerr = ( abs( w-wold(1) ) + 
     &              abs( ( w-wold(2) ) ) )
     &            / (area0 + w )

         area = area + result

         go to 500

      endif

      result = -(area+area0) / pi

***
**
**     We have integrated to find  \int_0^{\infty}.  The integral
**     required actually goes from -infty to +infty, but is
**     symmetric about the y-axis, so the integral is _twice_ the result
**     obtained above.
**
***



*     Determine the error
*     (Keep in this order so the most important aspect is returned)
*     Note also that if we don't use Sidi, and w=0 as set initially,
*     We have good relative error, so that line should be OK.
      if ( abs(w) .LT. aimrerr ) then
         exitstatus = -1
      else
         exitstatus = -10
      endif
      if ( abs(relerr) .LT. aimrerr ) exitstatus = 1


*************SIDI-T STUFF IF WE FIND WE NEED IT************************
**     FIND THE UPPER LIMIT OF t so that after this t, use Sidi.
*      sidit = fndsidit( p, mu, phi, y )
*      print *,'SIDI-T = ',sidit
*************SIDI-T STUFF IF WE FIND WE NEED IT************************


      return

      end

****************************************************************
****************************************************************

      subroutine cumfndkmax(p, phi, y, mu,
     &                      kmax, tmax, mmax, ier)

***
*     Finds  k_max and t_max
*     We solve Im{k'(t)}=0 using Newton's method.
***
*     IN:  p, phi, y, maxit, mu
*     OUT: kmax, tmax, mmax, ier
***

      double precision  kmax, tmax, imgddcgf,
     &          delta, z, p, largest, fhi, flo,
     &          phi, y, pi, front, mu, rl,
     &          cumsfzro, zhi, zlo, cumdk, dpmmax
      integer  mmax, ier, allok, lrgint

      external  cumdk, imgddcgf

***
* MAJOR VARIABLES:
*    kmax      : the maximum value of k {we have sin( Im k(t) ) }
*    tmax      : the t-value where the maximum occurs
*    mmax      : the maximum value of  m  we will have
*                (where sin( Im k(t) ) = m\pi )
***

      pi = acos( -1.0d00 )
*      pi = atan( 1.0d00 ) * 4.0d00
*      pi = 3.14159 26535 89793 23846 26433d00
      delta = 1.0d00
      ier = 0
      allok = 1
      largest = 1.0d30
      lrgint = 100000000

*     A good starting point is sometime crucial, so we spend a
*     little time finding a decent one.
*     We try the point of inflection and the approximation for
*     t  large (which can happen, for example, when y is very small).
*     We then choose the point for which the derivative is closest
*     to zero and still positive (ie on the left of the maximum).

      front = mu ** ( 1.0d00 - p ) / 
     &         ( phi * (1.0d00 - p ) )

      z = front * tan( pi/2.0d00 * (1.0d00-p) /
     &                 ( 2.0d00*p - 1.0d00 ) )
      z = abs( z )

      zlo = 0.0d00
      zhi = z

      flo = cumdk(p, phi, y, mu, zlo )
      fhi = cumdk(p, phi, y, mu, zhi )

*     Now we ensure the zero is bounded; if not, we venture further out
 565  if ( ( allok .EQ. 1 ) .AND.
     &     ( fhi*flo .GT. 0.d00 ) ) then

         zlo = zhi
         zhi = 1.1d00 * zhi + 1.0d-05

         flo = cumdk(p, phi, y, mu, zlo )
         fhi = cumdk(p, phi, y, mu, zhi )

         if ( zhi .GT. largest/10.d00 ) allok = 0

         goto 565

      endif

      if ( allok .EQ. 0 ) then

         ier = -80
         z = 0.0d00
         kmax = 0.0d00
         mmax = 0.0d00
         tmax = 0.0d00
         return

      else

         flo = cumdk(p, phi, y, mu, zlo )
         fhi = cumdk(p, phi, y, mu, zhi )

*        This is linear interpolation:
         z = zlo - flo * ( zhi - zlo ) / ( fhi - flo )
         z = cumsfzro( p, phi, y, zlo, zhi, z, mu,
     &                 cumdk, imgddcgf, ier )
      endif

*     We proceed even if ier .NE. 0 since the calling function
*     may be happy about it.  Eg, if the corresponding  mmax
*     is larger than  maxit, it won't matter if  kmax  isn't
*     found to any accuracy, provided  mmax>maxit.

      tmax = z

*     Now we find the value of  kmax  at  tmax:
      call cumcalccgf( p, phi, y, mu, tmax, rl, kmax )

*     Now, we can have real troubles here if  p  gets large
*     since them psi= datan( (1.0d00 - p) * x * phi ) gets
*     really close to -pi/2.  I'm not sure how to best deal with
*     this, but a side-effect is sometimes that kmax<0, which it
*     shouldn't be.  In any case, kmax<0 means PROBLEMS.
*
*     What bothers me is how to detect problems when there *are* problems
*     but  kmax  doesn't become negative...?

      if ( kmax .LT. 0 ) then

         kmax = abs(kmax)
         mmax = lrgint
         allok = 0

      else

         dpmmax = kmax / pi

         if ( dpmmax .GT. dble( lrgint ) ) then
            mmax = lrgint
            allok = 0
         else
            mmax = int( dpmmax )
         endif


      endif

      return
      end

****************************************************************
****************************************************************

      double precision function cumf(p, phi, y, mu, x)

***
*     A function to be numerically integrated in cum density.
*     Note that
*        lim (n->infty) exp( Re k) = exp( -lambda).
*     cumf2  is used for the conditional density; cumf  for p>2.
***
*     IN:  p, phi, y, mu, x
*     OUT: cumf
***

      double precision  x, p, phi, y, mu, lambda,
     &                  calclambda, rl, im

***
*
* MAJOR VARIABLES:
*   p          : the index in the variance function, V(mu) = phi * mu^p
*   phi        : the dispersion parameter
*   y          : the point at which the density function is to be evaluated
*   mu         : the mean value
*   x          : the internal variable over which to integrate in this function
*   rl         : the real part of the cgf
*   im         : the imaginary part of the cgf
*   lambda     : P(Y=0) = exp( -lambda ) when 1<p<2
*
***

      lambda = calclambda( p, phi, mu )

      if ( x .EQ. 0.0d00 ) then

         cumf = mu - y

      else

         call cumcalccgf(p, phi, y, mu, x, rl, im)
         cumf = dexp( rl ) * dsin( im ) / x

      endif

      return
      end

*******************************************************************
****************************************************************

      double precision function cumf2(p, phi, y, mu, x)

***
*     A function to be numerically integrated in cum density
*     using the conditional density.  Note that
*        lim (n->infty) exp( Re k) = exp( -lambda).
*     cumf  is used for p>2.
***
*     IN:  p, phi, y, mu, x
*     OUT: cumf2
***

      double precision  x, p, phi, y, rl, im, mu,
     &                  calclambda, lambda

***
* MAJOR VARIABLES:
*   p          : the index in the variance function, V(mu) = phi * mu^p
*   phi        : the dispersion parameter
*   y          : the point at which the density function is to be evaluated
*   mu         : the mean value
*   x          : the internal variable over which to integrate in this function
*   rl         : the real part of the cgf
*   im         : the imaginary part of the cgf
*
***

      lambda = calclambda( p, phi, mu )

      if ( x .EQ. 0.0d00 ) then

         cumf2 = ( mu - y ) + y * exp( -lambda )

      else

         call cumcalccgf(p, phi, y, mu, x, rl, im)

         cumf2 = exp( rl ) * sin( im ) + 
     &           exp( -lambda ) * sin( x*y )
         cumf2 = cumf2 / x

      endif

*     Now, we use the conditional density so we apply the tranform:
      cumf2 = cumf2 / ( 1.0d00 - exp( -lambda ) )

      return
      end

*******************************************************************
*******************************************************************

      double precision function cumdk(p, phi, y, mu, x)

***
*     Evaluates the derivative of the imaginary part of the  k  function
***
*     IN:  p, phi, y, mu, x
*     OUT: cumdk
***

      double precision  p, phi, y, mu, x, rl

***

      call cumcalcdcgf( p, phi, y, mu, x, rl, cumdk )

      return
      end

*******************************************************************
****************************************************************

      subroutine cumcalccgf(p, phi, y, mu, x, rl, im)

***
*     Calculates the cgf for the cumulative distribution
***
*     IN:  p, phi, y, mu, x
*     OUT: rl, im
***

      double precision  p, phi, y, mu, x, rl, im,
     &                  omega, alpha, front, denom

***
* MAJOR VARIABLES:
*   rl         : the real part of the cgf
*   im         : the imaginary part of the cgf
***

      omega = atan( ( 1.0d00 - p ) * x * phi /
     &              ( mu ** (1.0d00-p) ) )

      alpha = ( 2.0d00 - p ) / ( 1.0d00 - p )
      front = ( mu ** ( 2.0d00 - p ) ) / 
     &         ( phi * ( 2.0d00 - p ) )
      denom = cos( omega ) ** alpha

      rl = front * cos ( omega * alpha ) / denom - front
      im = ( front * sin ( omega * alpha ) / denom ) - 
     &         x*y

      return
      end

****************************************************************
****************************************************************

      subroutine cumcalcdcgf(p, phi, y, mu, x, rl, im)

***
*     Calculates the derivative of the cgf for the cumulative distribution
***
*     IN:  p, phi, y, mu, x
*     OUT: rl, im
***

      double precision  p, phi, y, mu, x, rl, im,
     &                  omega, alpha, denom

***
* MAJOR VARIABLES:
*   rl         : the real part of the derivative of the cgf
*   im         : the imaginary part of the derivative of the cgf
***

      omega = atan( ( 1.0d00 - p ) * x * phi /
     &              ( mu ** (1.0d00-p) ) )

      alpha = 1.0d00 / ( 1.0d00 - p )
      denom = cos( omega ) ** alpha

      rl = -mu * ( sin( omega * alpha ) / denom )
      im = mu * cos ( omega * alpha ) / denom - y

      return
      end

****************************************************************
****************************************************************

      double precision function cumsfzro( p, phi, y, x1, 
     &            x2, x0, mu, fun, dfun, ier )

******************************************************************
*     Uses a modified Newton's Method to find a root between
*     x1 and x2.  Used when functions called doesn't need
*     an  m  input; otherwise, use cumsfzro2.
***
*     IN:  p, phi, y, x1, x2, x0, mu, acc, fun, dfun
*     OUT: cumsfzro
***

      double precision  df, dx, dxold, f, fh, fl, temp, 
     &            xh, xl, x1, x2, fun, dfun, y, mu,
     &            p, phi, x0
      integer  j, maxit, ier

      external  fun, dfun

***

*     SET PARAMETERS
      maxit = 100
      ier = 0

      fl = fun(p, phi, y, mu, x1)
      fh = fun(p, phi, y, mu, x2)

      if ( fl .EQ. 0.0d00 ) then
         cumsfzro = x1
         return
      elseif ( fh .EQ. 0.0d00 ) then
         cumsfzro = x2
         return
      elseif ( fl .LT. 0.0d00 ) then
         xl = x1
         xh = x2
      else
         xl = x2
         xh = x1
      endif


      if ( ( x0 .GT. xl ) .AND. ( x0 .LT. xh ) ) then
         cumsfzro = x0
      else
         cumsfzro = ( xl + xh ) / 2.0d00
      endif
      dxold = abs( x2-x1 )
      dx = dxold

      f = fun( p, phi, y, mu, cumsfzro )
      df = dfun( p, phi, y, mu, cumsfzro )

      do j = 1, maxit

         if ( ( (cumsfzro-xh)*df-f)*((cumsfzro-xl)*df-f)
     &            .GT.0.0d00
     &                .OR.
     &        abs(2.0d00*f) .GT. abs( dxold*df) ) then

            dxold = dx
            dx = 0.5d00 * ( xh - xl )
            cumsfzro = xl + dx

            if ( xl .EQ. cumsfzro ) then
               return
            endif

         else

            dxold = dx
            dx = f/df
            temp = cumsfzro
            cumsfzro = cumsfzro - dx
            if ( temp .EQ. cumsfzro ) then
               return

            endif

         endif

         if ( abs( dx ) .LT. 1.0d-11 ) then
            return
         endif

         f = fun( p, phi, y, mu, cumsfzro )
         df = dfun( p, phi, y, mu, cumsfzro )


         if ( f .LT. 0.0d00 ) then
            xl = cumsfzro
         else
            xh = cumsfzro
         endif

      enddo

*     If we get this far, we haven't returned and so there is an error:
      ier = -20

      return
      end

****************************************************************
******************************************************************

      subroutine cumsmallp( p, phi, y, mu, aimrerr, 
     &               resulta, ier, relerr, its, verbose )

***
*     Calculates the density in the case of distributions with
*     1 < p < 2.
***

      double precision  p, phi, y, pi, area, aimrerr, 
     &         lambda, relerr, result, sidit, cumsfzro,
     &         cumf2,  cumintim, imgddcgf, flo, fhi,
     &         mmatrix(2, 200), nmatrix(2, 200),
     &         xvec(200), w, wold(3), area0, sumarea,
     &         cumdk, cumsfzro2, mu, calclambda, area1,
     &         largest,smallest, resulta, t0, z0, cumf,
     &         delta, oldft, deriv1, deriv2, omega, fndst,
     &         tinfl, ft, zold, z1, dt, olddt, tmp,
     &         alpha, a2, oldt, front, top, bot, rek, imk,
     &         cumddk, rl, workvec1, workvec2, lobnd,
     &         hibnd
      integer  m, ier, n, go, its, ceil, myfloor,
     &         tpt, mdirn, boost, flag, its1, verbose

      external  cumf2, cumdk, cumsfzro2, cumintim, 
     &          imgddcgf, cumsfzro, ceil, myfloor, cumddk, 
     &          cumf, fndst

***

*     SET ACCURACY REQUIREMENTS
      largest = 1.0d30
      smallest = 1.0d-30

*     SET OTHER PARAMETERS
      pi = acos( -1.0d00 )
      wold(1) = 0.0d00
      wold(2) = 0.0d00
      wold(3) = 0.0d00

      m = 0
*     The zero to solve for
      n = 0
*     The turning point we solving for
      deriv2 = mu - y
*     The initial derivative at t=0
      tpt = 1
*     A flag to indicate if turning points found
      tinfl = 0
*     Inflection points at these t
      oldft = 0.0d00
*     Last value of Im k(t)
      ft = 0.0d00
*     Current value of Im k(t)
      area = 0.0d00
      area0 = 0.0d00
      area1 = 0.0d00
*     The integrals values
      zold = 0.0d00
*     The previous zero
      z1 = 0.0d00
*     The current zero
      dt = 0.0d00
*     The current location of Im k'(t) = 0
      olddt  = 0.0d00
*     The previous location of Im k'(t) = 0
      its1 = 0
*     The number of iterations until sidi is used

      alpha = ( 2.0d00 - p ) / ( 1.0d00 - p )
      lambda = calclambda(p, phi, mu )
      a2 = p / ( 1.0d00 - p )
*     Find the value of t after which we can use Sidi acceleration
      sidit = fndst( p, phi, mu)

      if ( sidit .GT. 0.0d00 ) then
         tpt = 1
      else
         tpt = 0
      endif

 4000 if ( tpt .EQ. 1 ) then

*        Update
         n = n + 1
         oldt = tinfl
         deriv1 = deriv2

*        Find inflection points
         front = mu ** (1.0d00 - p ) / 
     &         ( phi * ( 1.0d00 - p ) )
         tinfl = front * tan( dble(n)*pi*(1.0d00-p)/p )

*        Find sign of derivative
         omega = atan( (1.0d00-p) * phi * tinfl * 
     &           mu**(p-1.0d00) )
         top =  omega * alpha
         bot = cos(omega)**alpha
         Rek = lambda * ( cos(top) / bot - 1.0d00 )
         Imk = lambda * sin(top)/bot - tinfl * y
         deriv2 = mu * cos( top ) / bot - y

         tpt = 0
*        test for turning point

         if ( ( deriv1 * deriv2 ) .LT. 0.0d00 ) then

            tpt = 1

*           Now find the zero, a turning point of Im k(t)
            olddt = dt
            t0 = ( oldt + tinfl ) / 2.0d00
            dt = cumsfzro(p, phi, y, oldt, tinfl,
     &                    t0, mu, cumdk, cumddk, ier )

*           Determine function values here
            oldft = ft
            call cumcalccgf(p, phi, y, mu, dt, rl, ft)
*           The ft is the result we needed

*           Determine m values
            workvec1 = oldft / pi
            workvec2 = ft / pi

*           Find direction is which m's should go
            tmp = sign( ft-oldft, 1.0d00)
            if ( tmp .GT. 0.0d00 ) then
               mdirn = 1
            else
               mdirn = -1
            endif

*           Determine the m's
            if ( mdirn .GT. 0 ) then
               workvec1 = dble( ceil(workvec1) )
               workvec2 = dble( myfloor(workvec2) )
            else
               workvec1 = dble( myfloor(workvec1) )
               workvec2 = dble( ceil(workvec2) )
            endif

*           Set up bounds
            lobnd = olddt
            hibnd = dt

*           Now, work at each m
            do m = int(workvec1), int(workvec2), mdirn

*              Find the zeros
               zold = z1
               z0 = ( lobnd + hibnd ) / 2.0d00
               z1 = cumsfzro2( p, phi, y, lobnd, hibnd,
     &                        z0, mu,
     &                        cumintim, cumdk, m, ier)

*              Integrate
               call gaussq( cumf, result, zold, z1,
     &                      p, phi, y, mu )

               area0 = area0 + result
               its1 = its1 + 1

*              Update
               lobnd = z1

            enddo

         endif

         go to 4000

      endif

*     Now do the cases where m is easy to find!
      its = 0
      go = 1

*     First, we evaluate up to sidit without acceleration
*     If sidit=0, we do one extra iteration anyway, since
*     sidi's algorithm cannot have a  xvec  being 0.
*
*     So we do this step for at least one iteration whatever.

      area1 = 0.0d00
      zold = z1

 4300 if ( go .EQ. 1 ) then

         its = its + 1

*        Determine m
         m = m - 1

*        Find zero
         lobnd = zold
         hibnd = zold + pi / y
         flo = cumintim( p, phi, y, lobnd, mu, m)
         fhi = cumintim( p, phi, y, hibnd, mu, m)

         boost = 1
*        boost enables us to get a region a bit quicker
 4550    if ( flo * fhi .GT. 0.0d00 ) then

            boost = boost + 1

*           Try harder for bounds
            delta = hibnd - lobnd
            lobnd = hibnd
            hibnd = hibnd + dble(boost) * delta / 2.0d00

            flo = cumintim( p, phi, y, lobnd, mu, m)
            fhi = cumintim( p, phi, y, hibnd, mu, m)

            go  to 4550

         endif

         z0 = ( lobnd + hibnd ) / 2.0d00

         z1 = cumsfzro2( p, phi, y, lobnd, hibnd,
     &                  z0, mu,
     &                  cumintim, cumdk, m, ier)

*        Integrate
         call gaussq( cumf, result, zold, z1,
     &                p, phi, y, mu )
         area1 = area1 + result

         its1 = its1 + 1

*        Update
         zold = z1
         go = 1

         if ( z1 .GT. sidit ) then
            go = 0
*           Time to start accelerating!
         endif

         go to 4300

      endif

*     area1 is now this area


*     Now start accelerating
      go = 1
      xvec(1) = z1
      its = 0
      w = 0

 4500 if ( go .EQ. 1 ) then

         its = its + 1

*        Determine m
         m = m - 1

*        Find zero
         lobnd = zold + 1.0d-05
         hibnd = zold + pi / y
         flo = cumintim( p, phi, y, lobnd, mu, m)
         fhi = cumintim( p, phi, y, hibnd, mu, m)

         boost = 1
*        boost enables us to get a region a bit quicker
 4350    if ( flo * fhi .GT. 0.0d00 ) then

            boost = boost + 1

*           Try harder for bounds
            delta = hibnd - lobnd
            lobnd = hibnd
            hibnd = hibnd + dble(boost) * delta / 2.0d00

            flo = cumintim( p, phi, y, lobnd, mu, m)
            fhi = cumintim( p, phi, y, hibnd, mu, m)

            go  to 4350

         endif

         z0 = ( lobnd + hibnd ) / 2.0d00

         z1 = cumsfzro2( p, phi, y, lobnd, hibnd,
     &                  z0, mu, cumintim, cumdk, m, ier)

*        Integrate
         call gaussq( cumf, result, zold, z1,
     &                p, phi, y, mu )
         area = area + result

*        Accelerate
         xvec( its+1 ) = z1
         call sidiacc( area, result, xvec, mmatrix,
     &          nmatrix,w, its, relerr, wold, sumarea,
     &          flag, verbose )

*        Update
         zold = z1
         go = 1
         if ( abs(relerr) .LT. aimrerr ) then
            go = 0
         endif

         go to 4500

      endif

*     Combine all the areas
      area = w + area1 + area0

      lambda = mu**(2.0d00-p) / ( phi * (2.0d00-p) )

      resulta = -1.0d00/(pi * ( 1.0d00 - exp(-lambda) ))
     &               * area -
     &          exp(-lambda) / 
     &          (2.0d00 * (1.0d00 - exp(-lambda)))

      its = its + its1

      return

      end

****************************************************************
****************************************************************

      subroutine cumbigp( p, phi, y, mu, aimrerr, 
     &         resulta, maxit, ier, exitstatus, relerr, 
     &         its, verbose )

***
*     Calculates the density in the case of distributions when
*     p>2
***

      double precision  p, phi, y,  pi, area, aimrerr,
     &          relerr, result, zero1, zero2, cumf,
     &          cumintim, imgddcgf, flo, fhi,
     &          mmatrix(2, 200), nmatrix(2, 200),
     &          xvec(200), w, wold(3), area0, sumarea,
     &          cumdk, cumsfzro2, mu, cumsfzro,
     &          largest,smallest, resulta,
     &          area1, zero, zerofn, zerodfn, zstep2,
     &          kmax, tmax, zlo, zhi, zstep
      integer  m, its, ier, maxit, flag, allok, kmaxok, 
     &         mmax, firstm, exitstatus, verbose

      external  cumdk, cumsfzro2, cumintim, imgddcgf,
     &          cumsfzro, cumf,zerofn, zerodfn

***
*     SET ACCURACY REQUIREMENTS
      largest = 1.0d30
      smallest = 1.0d-30

*SET M and N matrices, and other sidi things, to 0?


*     SET OTHER PARAMETERS
      verbose = 0
      pi = acos( -1.0d00 )
*      pi = 3.14159 26535 89793 23846 26433d00
      area = 0.0d00
      area0 = 0.0d00
      area1 = 0.0d00
      result = 0.0d00
      its = 0
      relerr = 1.0d00
      flag = 0


      wold(1) = 0.0d00
      wold(2) = 0.0d00
      wold(3) = 0.0d00

      if ( y .GE. mu ) then
*        In this case, Im(k) heads down straight away.

         kmax = 0.0d00
         mmax = -1
         tmax = 0.0d00
         m = -1

*        FIND ZEROS
         zero1 = 0.0d00
*        An approximation to the first zero:
         zero = pi / y

*        Now try to find this first zero:
         zlo = 0.0d00
         zhi = zero * 2.0d00


         if ( zhi .GT. largest/10.0d00 ) then

            allok = 0
            flo = 0.0d00
            fhi = 0.0d00

         else

            allok = 1
            flo = cumintim( p, phi, y, zlo, mu, m )
            fhi = cumintim( p, phi, y, zhi, mu, m )

         endif

 565     if ( ( allok .EQ. 1 ) .AND.
     &        (fhi * flo ) .GT. 0.0d00 ) then

            zlo = zhi
            zhi = zhi * 1.5d00

            flo = cumintim( p, phi, y, zlo, mu, m )
            fhi = cumintim( p, phi, y, zhi, mu, m )

            if ( zhi .GT. largest/10.0d00 ) allok = 0

            goto 565

         endif

         if ( zhi .GT. largest/10.0d00 ) allok = 0

         if ( allok .EQ. 0 ) then
            ier = -50
            result = 0.0d00

            return
         endif

         zero2 = cumsfzro2( p, phi, y, zlo, zhi, zero, 
     &               mu, cumintim, cumdk, m, ier )

         xvec( 1 ) = zero2

*        special case: between 0 and first zero
         call gaussq( cumf, area0, zero1, zero2, 
     &                p, phi, y, mu )

*        Now do some more integrations and use sidi acceleration
  500    if (    ( ( its .LT. 4 ) .AND.
     &             ( flag .NE. 1 )
     &           )
     &        .OR.
     &           ( ( its .LT. maxit ) .AND.
     &             ( flag .NE. 1 ) .AND.
     &             ( abs(relerr) .GT. aimrerr )
     &           )
     &      ) then

*           get next zeros
            m = m - 1
            zstep = zero2 - zero1
            zero1 = zero2
            zlo = zero2 - 0.01d00 * zstep
            zhi = zero2 + 2.0d00 * zstep

            if ( zhi .GT. largest/10.0d00 ) then

               allok = 0
               flo = 0.0d00
               fhi = 0.0d00

            else

               allok = 1
               flo = cumintim( p, phi, y, zlo, mu, m )
               fhi = cumintim( p, phi, y, zhi, mu, m )

            endif

 765        if ( ( allok .EQ. 1 ) .AND.
     &           (fhi * flo ) .GT. 0.0d00 ) then

               zlo = zhi
               zhi = zhi * 1.5d00

               flo = cumintim( p, phi, y, zlo, mu, m )
               fhi = cumintim( p, phi, y, zhi, mu, m )

               if ( zhi .GT. largest/10.0d00 ) allok = 0

               goto 765

            endif

            if ( zhi .GT. largest/10.0d00 ) allok = 0

            if ( allok .EQ. 0 ) then

               ier = -50
               result = 0.0d00

               return

            endif

            zero2 = cumsfzro2( p, phi, y, zlo, zhi, zero, 
     &                  mu, cumintim, cumdk, m, ier )

*           integrate between zeros
            call gaussq( cumf, result, zero1, zero2,
     &                     p, phi, y, mu )
     
*           Update iteration count
            its = its + 1

*           accelerate convergence of infinite sequence
            xvec( its+1 ) = zero2

            call sidiacc( area, result, xvec, mmatrix,
     &           nmatrix, w, its, relerr, wold, sumarea, 
     &           flag, verbose )

            relerr = ( abs( w-wold(1) ) + 
     &                  abs( ( w-wold(2) ) ) )
     &               / (area0 + w)

            area = area + result

            go to 500

         endif

         if ( ( its .GE.  maxit ) .AND.
     &        ( abs(relerr) .GT. aimrerr) ) then
            ier = -40
         endif

         if ( flag .EQ. 1 ) then
            ier = -70
         endif

         result = area0 + w

*        Now, the very first integration has not been counted
*        (since that can foul up the Sidi acceleration iteration count)
*        so update now
         its = its + 1
         
      else

*        that is:    if ( y .LT. mu) then
*        In this case, Im(k) may head up before going to  -infinity

*        FIND k_max AND t_max
         kmaxok = 1
         call cumfndkmax(p, phi, y, mu, 
     &            kmax, tmax, mmax, ier)

         if ( ier .NE. 0 ) then
             kmaxok = 0
         endif

*        IF KMAX<pi/2 AND NO CONVERGENCE?

         if ( kmax .LT. pi ) then
*           In this case, the upward trend doesn't reach pi/2.  Thus,
*           the first zero is where Im(k) first crosses the t/x-axis
*           again; that is, at m=0.

            if ( ier .EQ. -80 ) then
               result = 0.0d00
               return
            endif

            m = 0

*           FIND ZEROS
            zero1 = 0.0d00
            zero = tmax + ( pi / y )

*           BOUNDS ON OTHER ZERO:
            zlo = tmax
            zhi = zero*1.2d00

            if ( zhi .GT. largest/10.0d00 ) then

               allok = 0
               flo = 0.0d00
               fhi = 0.0d00

            else

               allok = 1
               flo = cumintim( p, phi, y, zlo, mu, m )
               fhi = cumintim( p, phi, y, zhi, mu, m )

            endif

 1565       if ( (allok .EQ. 1 ) .AND.
     &            ( (fhi * flo ) .GT. 0.0d00 ) ) then

               zlo = zhi
               zhi = zhi * 1.5d00

               flo = cumintim( p, phi, y, zlo, mu, m )
               fhi = cumintim( p, phi, y, zhi, mu, m )

               if ( zhi .GT. largest/10.0d000 ) allok = 0

               goto 1565

            endif

            if ( zhi .GT. largest/10.0d00 ) allok = 0

            if ( allok .EQ. 0 ) then
               ier = -50
               result = 0.0d00

               return
            endif

            zero2 = cumsfzro2( p, phi, y, zlo, zhi, zero,
     &                  mu, cumintim, cumdk, m, ier )


            xvec(1) = zero2

*           integrate between zeros
            call gaussq( cumf, area0, zero1, zero2, 
     &                   p, phi, y, mu )

  600       if (    ( ( its .LT. 4 ) .AND.
     &                 ( flag .NE. 1 )
     &               )
     &            .OR.
     &               ( ( its .LT. maxit ) .AND.
     &                 ( flag .NE. 1 ) .AND.
     &                 ( abs(relerr) .GT. aimrerr )
     &               ) ) then

*              get next zeros
               m = m - 1

               zstep = ( zero2 - zero1 )
               zero1 = zero2
               zlo = zero2 - 0.01d00 * zstep
               zhi = zero2 + 2.0d00 * zstep

               if ( zhi .GT. largest/10.d00 ) then

                  allok = 0
                  flo = 0.0d00
                  fhi = 0.0d00

               else

                  allok = 1
                  flo = cumintim( p, phi, y, zlo, mu, m )
                  fhi = cumintim( p, phi, y, zhi, mu, m )

               endif

 1665          if ( ( allok .EQ. 1 ) .AND.
     &              (fhi * flo ) .GT. 0.0d00 ) then

                   zlo = zhi
                   zhi = zhi + zstep

                  flo = cumintim( p, phi, y, zlo, mu, m )
                  fhi = cumintim( p, phi, y, zhi, mu, m )

                  if ( zhi .GT. largest/10.0d000 ) then
                     allok = 0
                  endif

                  goto 1665

               endif

               if ( zhi .GT. largest/10.d00 ) allok = 0

               if ( allok .EQ. 0 ) then
                  ier = -50
                  result = 0.0d00

                  return

               endif

*              Approximate with a linear interpolation
               zero = zlo - flo * ( zhi - zlo ) / 
     &                  ( fhi - flo )

               zero2 = cumsfzro2( p, phi, y, zlo, zhi, 
     &                         zero, mu, cumintim, cumdk, 
     &                         m, ier )

*              integrate between zeros
               call gaussq( cumf, result, zero1, zero2,
     &                        p, phi, y, mu )

*              Update iteration count
               its = its + 1
               
*              accelerate convergence of infinite sequence
               xvec( its+1 ) = zero2

               call sidiacc( area, result, xvec, mmatrix, 
     &            nmatrix, w, its, relerr, wold, 
     &            sumarea, flag, verbose )

               area = area + result

               go to 600

            endif

            if ( ( its .GE.  maxit ) .AND.
     &           ( abs(relerr) .GT. aimrerr) ) then
               ier = -40
            endif

            if ( flag .EQ. 1 ) then
               ier = -70
            endif

            result = area0 + w
            
*           Now, the very first integration has not been counted
*           (since that can foul up the Sidi acceleration iteration count)
*           so update now
            its = its + 1

         else
*           that is:   case where kmax >= pi
*           In this case, the upward trend is goes above  pi, and so
*           the first zero will be at  m=1.

*           Now, kmax may not have been found accurately.  IF, however,
*           the corresponding max > maxit, it won't matter a great deal
*           and we can proceed.  If not, accuracy cannot be ensured
*           unless the maximum  m  used is less than maxit.  So we
*           test  ier.

            if ( ier .EQ. -80 ) then
               result = 0.0d00
               return
            endif

            wold(1) = 0.0d00
            wold(2) = 0.0d00
            wold(3) = 0.0d00

*           FIND ZEROS
            zero1 = 0.0d00
            zero = pi / ( mu - y )
            m = 1
            firstm = 1

            zlo = smallest
            zhi = tmax

            if ( zhi .GT. largest/10.d00 ) then

              allok = 0
              flo = 0.0d00
              fhi = 0.0d00

            else

               allok = 1
               flo = cumintim( p, phi, y, zlo, mu, m )
               fhi = cumintim( p, phi, y, zhi, mu, m )

            endif

            zstep = zhi - zlo

 2565       if ( ( allok .EQ. 1 ) .AND.
     &           (fhi * flo ) .GT. 0.0d00 ) then


               zlo = zhi
               zhi = zhi + 0.1d00 * zstep

               flo = cumintim( p, phi, y, zlo, mu, m )
               fhi = cumintim( p, phi, y, zhi, mu, m )

               if ( zhi .GT. largest / 10.0d00 ) then
                  allok = 0
               endif

               goto 2565

            endif

            if ( zhi .GT. largest/10.d00 ) allok = 0

            if ( allok .EQ. 0 ) then
               ier = -50
               result = 0.0d00

               return
            endif

            zero2 = cumsfzro2( p, phi, y, zlo, zhi, zero, 
     &                  mu, cumintim, cumdk, m, ier )

            xvec( 1 ) = zero2

            call gaussq( cumf, area0, zero1, zero2, 
     &                   p, phi, y, mu )

****NOW INTEGRATE UP TO  tmax  WITHOUT ACCELERATION???
*   This seems to make sense to ensure that Sidi's method converges, but
*   there will be cases where we will take a loooonnnggg time to reach
*   tmax, and we will `never' use Sidi's method.

            zstep2 = zero2-zero1

  700       if (    ( ( its .LT. 4 ) .AND.
     &                ( flag .NE. 1 )
     &              )
     &           .OR.
     &              ( ( its .LT. maxit ) .AND.
     &                ( flag .NE. 1 ) .AND.
     &                ( abs(relerr) .GT. aimrerr )
     &              ) ) then

*              get next zeros

               zlo = zero2 - 1.0d-05 * zstep2
               zhi = zero2 + 2.0d00 * zstep2

               zero1 = zero2

*              FIND THE NEXT VALUE OF  m
               call nextm( tmax, mmax, zero2, m, firstm,
     &                     zlo, zhi, zero )

               if ( zhi .GT. largest/10.d00 ) then

                  allok = 0
                  flo = 0.0d00
                  fhi = 0.0d00

               else

                  allok = 1
                  flo = cumintim( p, phi, y, zlo, mu, m )
                  fhi = cumintim( p, phi, y, zhi, mu, m )

               endif

 2665          if ( ( allok .EQ. 1 ) .AND.
     &              (fhi * flo ) .GT. 0.0d00 ) then

                  zstep = ( zhi - zlo )
                  zlo = zhi
                  zhi = zhi + 0.2d00*zstep2

                  flo = cumintim( p, phi, y, zlo, mu, m )
                  fhi = cumintim( p, phi, y, zhi, mu, m )

                  if ( zhi .GT. largest/10.0d00 ) then
                     allok = 0
                  endif

                  goto 2665

               endif

               if ( zhi .GT. largest/10.0d00 ) allok = 0

               if ( allok .EQ. 0 ) then
                  ier = -50
                  result = 0.0d00

                  return
               endif

               zero2 = cumsfzro2( p, phi, y, zlo, zhi,
     &                         zero, mu, cumintim, 
     &                        cumdk, m, ier )

               call gaussq( cumf, result, zero1, zero2,
     &                        p, phi, y, mu )
     
*              Update iteration count
               its = its + 1
               
*              accelerate convergence of infinite sequence
               xvec( its+1 ) = zero2
               call sidiacc( area, result, xvec, mmatrix,
     &                   nmatrix, w, its, relerr, wold, 
     &                   sumarea, flag, verbose )

               area = area + result

               go to 700

            endif

            if ( ( its .GE.  maxit ) .AND.
     &           ( abs(relerr) .GT. aimrerr) ) then
               ier = -40
            endif

            if ( flag .EQ. 1 ) then
               ier = -70
            endif

            if ( kmaxok .EQ. 0 ) then
*              IF finding kmax didn't converge...

               if ( mmax .LT. maxit ) kmaxok = 1
*              All should be OK if mmax is less than maximum
*              that mmax can reach

               if ( m .LT. mmax-1 )  kmaxok = 1
*              All should be OK if greatest value of m used
*              is less than the `turning' m value.

            endif

            if ( kmaxok .EQ. 0 )  then

               ier = -60

            endif

            result = area0 + w
            
*           Now, the very first integration has not been counted
*           (since that can foul up the Sidi acceleration iteration count)
*           so update now
            its = its + 1

         endif

      endif


***
**
**     We have integrated to find  \int_0^{\infty}.  The integral
**     required actually goes from -infty to +infty, but is
**     symmetric about the y-axis, so the integral is _twice_ the result
**     obtained above.
**
***

      resulta = -result / pi

*     Determine the error
*     (Keep in this order so the most important aspect is returned)
*     Note also that if we don't use Sidi, and w=0 as set initially,
*     We have good relative error, so that line should be OK.
      if ( abs(w) .LT. aimrerr ) then
         exitstatus = -1
      else
         exitstatus = -10
      endif
      if ( abs(relerr) .LT. aimrerr ) exitstatus = 1


*************SIDI-T STUFF IF WE FIND WE NEED IT************************
**     FIND THE UPPER LIMIT OF t so that after this t, use Sidi.
*      sidit = fndsidit( p, mu, phi, y )
*      print *,'SIDI-T = ',sidit
*************SIDI-T STUFF IF WE FIND WE NEED IT************************


      return

      end

******************************************************************
*******************************************************************

      double precision function cumddk(p, phi, mu, x)

***
*     Evaluates the second derivative of the imaginary
*     part of the  k  function
***
*     IN:  p, phi, y, mu, x
*     OUT: cumdk
***

      double precision  p, phi, mu, x, rl

***

      call cumcalcddcgf( p, phi, mu, x, rl, cumddk )

      return
      end

*******************************************************************
*****************************************************************

      integer function ceil( x )

***
* This function rounds  x  towards +infinity
***

      double precision  x

***

      if ( x .GT. 0 ) then
         ceil = int( x ) + 1
      else
         ceil = int( x )
      endif

      return
      end

*****************************************************************
*****************************************************************

      integer function myfloor( x )

***
* This function rounds  x  towards -infinity
***

      double precision  x

***

      if ( x .GT. 0 ) then
         myfloor = int( x )
      else
         myfloor = int( x ) - 1
      endif

      return
      end

*****************************************************************
*******************************************************************

      double precision function fndst( p, phi, mu)

***

      double precision  p, phi, mu, front, pi
      integer  mmax

      external  calclambda

***

      pi = acos( -1.0d00 )
*      pi = atan(1.0d00)*4.0d00
      if ( abs( p - 1.5) .LT. 1.0d-02 ) then
         fndst = 0.0d00
      else
*        Find the largest m for which there are turning points in exp( Re k)
         mmax = int( 1.0d00 / ( 2.0d00 * (1.0d00-p) ) )

*        The t this correpsonds to is...
         front = mu**(1.0d00-p) / ( phi * ( 1.0d00-p ) )
         fndst = front * tan( dble(mmax) * 
     &            pi * (1.0d00-p) )
         fndst = abs( fndst )

      endif

      return
      end

*******************************************************************
*******************************************************************

      double precision function cumintim( p, phi, y, 
     &                                    x, mu, m )

***
*     The imaginary part for zeroing the integrand.
***
*     IN:  p, phi, y, x, mu, m
*     OUT: cumintim
***

      double precision  x, p, phi, y, pi, im, rl, mu
      integer  m

***
* MAJOR VARIABLES:
*   rl         : the real part of the cgf
*   im         : the imaginary part of the cgf
***

      pi = acos( -1.0d00 )
*      pi = 3.14159 26535 89793 23846 26433d00

      if ( x .EQ. 0.0d00 ) then

         cumintim = -dble(m) * pi

      else

         call cumcalccgf( p, phi, y, mu, x, rl, im )
         cumintim = -dble(m)*pi + im

      endif

      return
      end

******************************************************************
****************************************************************

      double precision function cumsfzro2( p, phi, y, x1, 
     &            x2, x0, mu, fun, dfun, m, ier )

***
*     Uses a modified Newton's Method to find a root between
*     x1 and x2 to find kmax.  Used when functions called need
*     an  m  input; otherwise, use cumsfzro.
***
*     IN:  p, phi, y, x1, x2, x0, mu, acc, fun, dfun, m
*     OUT: cumsfzro2
***

      double precision  df, dx, dxold, f, fh, fl, temp, 
     &                  xh, xl, x1, x2, fun, dfun, y,
     &                  mu, p, phi, x0
      integer  j, maxit, ier, m

      external  fun, dfun

***

*     SET PARAMETERS
      maxit = 100
      ier = 0

      fl = fun(p, phi, y, x1, mu, m)
      fh = fun(p, phi, y, x2, mu, m)

      if ( fl .EQ. 0.0d00 ) then
         cumsfzro2 = x1
         return
      elseif ( fh .EQ. 0.0d00 ) then
         cumsfzro2 = x2
         return
      elseif ( fl .LT. 0.0d00 ) then
         xl = x1
         xh = x2
      else
         xl = x2
         xh = x1
      endif


      if ( ( x0 .GT. xl ) .AND. ( x0 .LT. xh ) ) then
         cumsfzro2 = x0
      else
         cumsfzro2 = ( xl + xh ) / 2.0d00
      endif

      dxold = abs( x2-x1 )
      dx = dxold

      f = fun( p, phi, y, cumsfzro2, mu, m )
      df = dfun( p, phi, y, mu, cumsfzro2 )

      do j = 1, maxit

         if ( ( (cumsfzro2-xh)*df-f)*
     &            ((cumsfzro2-xl)*df-f) .GT. 0.0d00
     &                .OR.
     &        abs(2.0d00*f) .GT. abs( dxold*df) ) then

            dxold = dx
            dx = 0.5d00 * ( xh - xl )
            cumsfzro2 = xl + dx

            if ( xl .EQ. cumsfzro2 ) then
               return
            endif

         else

            dxold = dx
            dx = f/df
            temp = cumsfzro2
            cumsfzro2 = cumsfzro2 - dx
            if ( temp .EQ. cumsfzro2 ) then
               return

            endif

         endif

         if ( abs( dx ) .LT. 1.0d-11 ) then
            return
         endif

         f = fun( p, phi, y, cumsfzro2, mu, m )
         df = dfun( p, phi, y, mu, cumsfzro2 )

         if ( f .LT. 0.0d00 ) then
            xl = cumsfzro2
         else
            xh = cumsfzro2
         endif

      enddo

*     If we get this far, we haven't returned and so there is an error:
      ier = -20

      return
      end

*******************************************************************
*****************************************************************

      subroutine cumcalcddcgf( p, phi, mu, x, ddrl, ddim)

***
*     The second derivative of the cgf.
***
*     IN:   p, phi, mu, x
*     OUT:  ddrl, ddim
***

      double precision  p, phi, mu, x, ddrl, ddim, 
     &                  denom, alpha, omega

***
* MAJOR VARIABLES:
*     ddrl      : the real part of the second derivative
*     ddim      : the imaginary part of the second derivative
***

      omega = atan( ( 1.0d00 - p ) * x * phi /
     &              ( mu ** (1.0d00-p) ) )

      alpha = p / ( 1.0d00 - p )
      denom = cos( omega ) ** alpha

      ddim = -phi * mu**p  * ( cos( omega * alpha ) 
     &            / denom )
      ddrl = -phi * mu**p *   sin ( omega * alpha ) 
     &            / denom

      return
      end

*****************************************************************
