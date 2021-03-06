20 January 2020 ; tweedie 2.3.3
  - Added outputs gamma.mean and gamma.phi to tweedie.convert()
  - Minor fixes

15 November 2017: tweedie 2.3.1
   - Updated pwteedie.series to fix a bug (reported by Lu Yang),
     where incorrect answers could sometimes be returned.
   - Other minor fixes
   
06 November 2017: tweedie 2.3.0
	- Fixed an issue with AICtweedie, where the incorrect AIC
	  was given when prior weights used (reported by David Scollnik)
	- Fixed a compilation error (in the subroutine smallp, where 
	  variables were declared as initialised (many thanks to 
	  Iñaki Úcar <i.ucar86@gmail.com>)
	- Other minor edits.

23 August 2017: tweedie 2.2.10
   - Kept it even more quiet


22 August 2017: tweedie 2.2.9

   - Fixed tweedie.f to keep it quiet more often (sometimes, diagnostic reports 
     meant for internal monitoring, were printed)
   - Fixed a problem reported by Gustavo Lacerda, where  ptweedie()  returned NaN
   - As a result, the series is now used in far more cases when  1<xi<2  for  ptweedie()
   - Changed  rtweedie()  algorithm for the case 1<xi<2 (thanks to Carlos J. Gil Bellosta)
   - New function  tweedie.convert()  added

19 December 2016: tweedie 2.2.5

   - Added CITATION file
   - Fixed an issue where the Tweedie cdf could return a value greater than one (reported by Jeremie Juste)
   - Minor tidy of FORTRAN code
   - Minor tidies in R code

June 2014: tweedie 2.1.9 
  - Some administrative fixes for CRAN

10 July 2013: tweedie 2.1.8

	- Fixed an issue where ptweedie() would fail for very small y; set this to 0 when y<1.0e-300
		(based on a report by Johann Cuenin)
		
15 January 2013: tweedie 2.1.7 
	- Admin release (e.g. .First.lib()  removed)
	- Some minor edits to manual
	- Added the  control  input to  tweedie.profile() 
	  (thanks to Giri Khageswor, DPI Victoria)
	- Minor fixes in the manual
	
01 November 2012: tweedie 2.1.5 
	- Changed the example in  tweedie-package  to execute faster (CRAN requirement)

31 October 2012: tweedie 2.1.4 
	- Fixed an error in dtweedie() that reported NA in the case when power=1 and phi!=1 
	  (thanks to Dina Farkas)
	- [dqpr]stable now in package  stabledist  rather then fBasics; fixed
	- Fixed some typos in the help for tweedie-package (thanks to Peng Yu)
	- rtweedie  reported an error if power=1; fixed (thanks for Peng Yu)
   - Edits to conform with new first argument of .Fortran (i.e. .NAME rather than name)
   - Added NAMESPACE
   
09 June 2011: tweedie 2.1.0 
	- Changed tweedie.profile to ignore values of p/xi outside (1,2)
	  rather than report an error.
	- In some unusual cases, when p/xi=0 was used (with add0=TRUE), the
	  mle of p was between 0 and 1 (which is impossible).  We report a
	  warning message to check the data and the call to tweedie.profile(),
	  but then set the mle to the value of p/xi giving the larger value
	  of the likelihood
	- Made the functions usually called by users able to accept xi or power
	- A few minor edits to FORTRAN code; some variables not declared

08 June 2011: tweedie 2.0.8 
	- Minor fixes to documentation
 	- Fixed the  add0  input
 	- Some minor changes to the code to tidy up adding the zero
	- If values of xi/p are given between 0 and 1, or less than 0, 
     they are now omitted (with a message) rather than creating an error.  
     If there are no values left after omitting the problem value,
     an error message is given.
 	
30 September 2010: tweedie 2.0.7 
   - Ensured tweedie.profile does not use power=1.  This case 
    (power=1 and phi not equal to 1) is too hard for me to deal with at present.
   - Fixed an error introduced in version 2.0.5, where the value of xi/vec/p.vec
     was set to 1.2 (y>=0) or 1.5 (y>0) when not explictly specified
   - Fixed an error that reported the wrong mle of phi when the mle 
     occured at an endpoint of the given  xi  values.

26 August 2010: tweedie 2.0.5 
   - Change to dtweedie.inversion to ensure density=0 is returned
     when y<0 (1<p<2) or p<-0 (when p>2)
   - Change to tweedie.profile to fix a problem that p=1 returned an error
   - Changed so that tweedie.profile works withg p/xi = 0 when add0=TRUE (default is FALSE)
   - Fixed some minor outputting messages (when verbose==2)  
   - Location of CITATION file moved to correct location
    
12 July 2010: tweedie 2.0.4 
   - Minor edits
   
18 December 2009: tweedie 2.0.3 
	- Changed default p.vec: There were too many values
	- Added the facility to refer to  p  as  xi  in line with GLMs text

17 November 2009: tweedie 2.0.1 
  - Changed the default p.vec when 1<p<2 to seq(1.2, 1.8, by=0.05)
    (it was  seq(1.2, 1.8, by=0.1)  )
  - Slightly changed default output (added sep="" to some paste commands)
  - Added  AICtweedie  to compute AIC for Tweedie glms

10 August 2009: tweedie 2.0.0 
  - Made  method="inversion" the default (was "series")
  - Slightly changed default output (added sep="" to some paste commands)

tweedie 1.6.9 
  - An error introduced earlier (unsure when exactly; prob v 1.6.1)
  - In trying to identify and fix the error, tidied some of the FORTRAN code
  - Made  do.smooth=TRUE  the default (was FALSE)
  - If  p.vec  is not supplied,  tweedie.profile  makes a sensible guess
  - Made  verbose=FALSE the default  (was TRUE)
  - Made minor changes to output when  verbose=FALSE

tweedie 1.6.8
  - Correct error in dtweedie.saddle y=0 for 1<p,2 
    (Thanks to Glenn Meyers <GMeyers@iso.com>)
  
tweedie 1.6.7
  - Correct error in rtweedie when power=1 and phi\ne 1 
    (Thanks to Frederic Gosselin <frederic.gosselin@cemagref.fr>)
  
tweedie 1.6.6
  - Corrected specification of GPL; no changes to functionality, code or documentation
  
tweedie 1.6.5
  - In dtweedie, added power==3 (inverse Gaussian) as a special case
  - Minor change in documentation

tweedie 1.6.4
  - Added  data, weights, and  offset  input arguments
  - Fixed a bug where "=>" was used rather than ">="
  
tweedie 1.6.3
  - Fixed a bug in tweedie.profile: exact zeros still allowed p=2 (but not p>2)
  
tweedie 1.6.2
  - Fixed a bug in  ptweedie.inversion  where variable  cdf  was undefined
  - Fixed a bug in  qtweedie  where any(power)<1  should have been  any(power<1)

tweedie 1.6.1
  - Fixed error in tweedie.f  (twice) where a float was used in a do-loop

tweedie 1.6
  - Fixed a bug in ptweedie: It now returns 0 when y<0 (rather than reporting an error)
  - Fixed a bug in ptweedie.inversion: y.len was undefined
  - Fixed a bug in dwteedie.inversion: if( any(phi<0) ) written as if (any(phi) <0 )

tweedie 1.5.3:
  - Built for R version 2.6.0
  - Changed the tweedie.plot to produce no default x and y labels
  
tweedie 1.5.2: 
  - Built for R version 2.5.1
  - Removed the inversion computation from the example for  dtweedie;
    it was very slow and hence caused problems

tweedie 1.5.1: 
  - Built for R version 2.5.0
  - Fixed citation details
  - Added note in dtweedie help page indicating where the methods
    are defined (Dunn and Smyth, 2007).
  - Added tweedie-package.Rd
  - Fixed a small bug in  dtweedie.inversion (replaced  p  with   power)

02 May 2007: tweedie  1.5
  - Built for R version 2.5
  - Made change to use of pmatch in function  tweedie
  - Swapped the order of methods 1 and 2 in the call of
    tweedie.inversion  to make consistent with a change
    in the associated paper Dunn and Smyth.
  - Fixed an uninitialized variable error in the FORTRAN code
    which sometimes caused the inversion procedure to fail.
  - Other minor changes (fixing spelling errors, updated
    Dunn and Smyth reference; etc.)

04 December 2006: tweedie 1.4
  - Some very minor changes (some cosmetic) in the code
  - Made the code more robust, exiting gracefully when Inf and NA appear
    in the computed log-likelihood
  - Slightly changed the output plot in tweedie.profile, to reflect
    that some value of  p.vec  may not have a corresponding likelihood
    computed accurately (or, indeed, at all)
  - Added new function   tweedie.plot  to make it easier to plot Tweedie
    densities

tweedie 1.02  
  - Built for R version 2.0.0
  - Major fix: The mle estimate of phi was incorrect in tweedie.profile 
    (thanks Sarah Lennox, QDPI&F) 
  - Some better warnings in  tweedie.profile
  - tweedie.profile now allows the actual computed points to
    be added to the plot using the  do.points  option (default 
	 is TRUE)
  - The tweedie.profile function now smooth and computes 
    confidence intervals even when the computed likelihood 
	 values contain Infs or NAs.

tweedie 1.01  
  - Built for R version 1.9
  - fixed a bug in rtweedie: it failed for p=2 (thanks 
    Sarah Lennox, QDPI&F)
  - fixed a bug in qtweedie: it failed in some cases
   (thanks Gordon Smyth, WEHI)
  - added a proper  .First.lib  function (!) (thanks to
    James Wettenhall and Gordon Smyth, WEHI)
  - other minor bugs fixed
  - the  tweedie  family function removed, and package
    made dependent on the  statmod  package
  - updated for R version 1.9

tweedie 1.0  
  - Built for R version 1.7
