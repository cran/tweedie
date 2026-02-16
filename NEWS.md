tweedie 3.0.14 (Release date: 2026-02-16)
==============

Changes:

* Improved the vignette.
* Some internal renaming.
* Fix some xi = 1 cases (thanks Milan Bouchet-Valat).
* Relocate some messages in tweedie_profile().
* Add poison example to vignette.


tweedie 3.0.12 (Release date: 2026-02-07)
==============

Changes:

* Trying to fix bugs that pop up (seemingly at random) with rhub etc. checks.



tweedie 3.0.5 (Release date: 2026-01-30)
==============

Changes:

* FORTRAN code restructured to make the similar flow in the three zones 
  (initial; pre-acceleration; acceleration) clearer
* Some fixes to documentation to pass tests.
* Some minor fixes to R code.



tweedie 3.0.4 (Release date: 2026-01-20)
==============

Changes:

* Some fixes to implementation of IGexact


tweedie 3.0.3 (Release date: 2025-11-29)
==============

Changes:

* Add IGexact for [dp]tweedie_inversion: whether to use exact values or inversion when p = 3.
* Fixed some comments



tweedie 3.0.2 (Release date: 2025-11-29)
==============

Changes:

* All code moved from FORTRAN77 to FORTRAN90. Almost no FORTRAN code remains from version < 3.
* PDF and CDF computations consolidated and code shared where possible, substantially reducing the amount of FORTRAN code.
* Separated FORTRAN code into different files for easier debugging.
* Improvements to the acceleration algorithm and root-finding algorithms, so should work better for more cases.
* Added verbose (shows what's happening behind the scenes) and details (reports on the fitting) as options for many user-facing R functions.
* Added  ptweedie_inversion()  to the man page for  dtweedie.
* Tidied the man pages; added examples.
* Removed the almost-never used dtweedie.stable() function.
* Changed function names (e.g,  tweedie.convert()  to  tweedie_convert()).
* Separated R functions into separate files depending on purpose (e.g., dtweedie.R  and  ptweedie.R).
* Moved the  tweedie_Extra  files into the main package.
* dtweedie.igrand()  (now tweedie_igrand()) to plot the integrand for the DF also.


tweedie 2.3.5 (Release date: 2022-08-17)
==============

Changes:

* Added outputs gamma.mean and gamma.phi to tweedie.convert()
* Added more error checks to  tweedie.convert()  to prevent a model being provided
* Minor edits

tweedie 2.3.1 (Release date: 2017-11-15)
==============

Changes:

* Updated ptweedie.series() to fix a bug (reported by Lu Yang), where incorrect answers could sometimes be returned.
* Other minor fixes

tweedie 2.3.0 (Release date: 2017-11-06)
==============

Changes:

* Fixed an issue with AICtweedie(), where the incorrect AIC was given when prior weights used (reported by David Scollnik)
* Fixed a compilation error (in the subroutine smallp(), where variables were declared as initialised (thanks to Iñaki Úcar <i.ucar86@gmail.com>)
* Other minor edits.

tweedie 2.2.10 (Release date: 2017-08-23)
==============

Changes:

* Kept it even more quiet

tweedie 2.2.6 (Release date: 2017-08-22)
==============

Changes:

* Fixed tweedie.f() to keep it quiet more often (sometimes, diagnostic reports meant for internal monitoring, were printed)
* Fixed a problem reported by Gustavo Lacerda, where  ptweedie()  returned NaN
* As a result, the series is now used in far more cases when  1<xi<2  for  ptweedie()
* Changed  rtweedie()  algorithm for the case 1 < xi < 2 (thanks to Carlos J. Gil Bellosta)
* New function  tweedie.convert()  added

tweedie 2.2.5 (Release date: 2016-12-19)
==============

Changes:

* Added CITATION file
* Fixed an issue where the Tweedie cdf could return a value greater than one (reported by Jeremie Juste)
* Minor tidy of FORTRAN code
* Minor tidies in R code

tweedie 2.1.9 (Release date: 2014-06-06)
==============

Changes:

* Some administrative fixes for CRAN

tweedie 2.1.8 (Release date: 2013-09-10)
==============

Changes:

* Fixed an issue where ptweedie() would fail for very small y; set this to 0 when y<1.0e-300 (based on a report by Johann Cuenin)

tweedie 2.1.7 (Release date: 2013-01-15)
==============

Changes:

* Admin release (e.g. .First.lib()  removed)
* Some minor edits to manual
* Added the  control  input to  tweedie.profile() (thanks to Giri Khageswor, DPI Victoria)
* Minor fixes in the manual

tweedie 2.1.5 (Release date: 2012-11-01)
==============

Changes:

* Changed the example in  tweedie-package  to execute faster (CRAN requirement)

tweedie 2.1.4 (Release date: 2012-10-31)
==============

Changes:

* Fixed an error in dtweedie() that reported NA in the case when power=1 and phi != 1 (thanks to Dina Farkas)
* [dqpr]stable now in package  stabledist  rather then fBasics; fixed
* Fixed some typos in the help for tweedie-package (thanks to Peng Yu)
* rtweedie()  reported an error if power = 1; fixed (thanks for Peng Yu)
* Edits to conform with new first argument of .Fortran (i.e. .NAME rather than name)
* Added NAMESPACE

tweedie 2.1.0 (Release date: 2011-06-09)
==============

Changes:

* Changed tweedie.profile() to ignore values of p/xi outside (1, 2) rather than report an error.
* In some unusual cases, when p/xi = 0 was used (with add0 = TRUE), the mle of p was between 0 and 1 (which is impossible).  We report a warning message to check the data and the call to tweedie.profile(), but then set the mle to the value of p/xi giving the larger value of the likelihood
* Made the functions usually called by users able to accept xi or power
* A few minor edits to FORTRAN code; some variables not declared

tweedie 2.0.8 (Release date: 2011-06-08)
==============

Changes:

* Minor fixes to documentation
* Fixed the  add0  input
* Some minor changes to the code to tidy up adding the zero
* If values of xi/p are given between 0 and 1, or less than 0, they are now omitted (with a message) rather than creating an error. If there are no values left after omitting the problem value, an error message is given.

tweedie 2.0.7 (Release date: 2010-09-30)
==============

Changes:

* Ensured tweedie.profile() does not use power = 1.  This case (power=1 and phi not equal to 1) is too hard for me to deal with at present.
* Fixed an error introduced in version 2.0.5, where the value of xi.vec/p.vec was set to 1.2 (y >= 0) or 1.5 (y > 0) when not explicitly specified
* Fixed an error that reported the wrong mle of phi when the mle occurred at an endpoint of the given  xi  values.

tweedie 2.0.5 (Release date: 2010-08-27)
==============

Changes:

* Change to dtweedie.inversion() to ensure density = 0 is returned when y < 0 (1 < p < 2) or p <- 0 (when p > 2)
* Change to tweedie.profile() to fix a problem that p = 1 returned an error
* Changed so that tweedie.profile() works with p/xi = 0 when add0=TRUE (default is FALSE)
* Fixed some minor outputting messages (when verbose == 2)  
* Location of CITATION file moved to correct location

tweedie 2.0.4 (Release date: 2010-07-12)
==============

Changes:

* Minor edits

tweedie 2.0.3 (Release date: 2009-12-18)
==============

Changes:

* Changed default p.vec: There were too many values
* Added the facility to refer to  p  as  xi  in line with GLMs text

tweedie 2.0.1 (Release date: 2009-11-17)
==============

Changes:

* Changed the default p.vec when 1 < p < 2 to seq(1.2, 1.8, by = 0.05) (it was  seq(1.2, 1.8, by = 0.1)  )
* Slightly changed default output (added sep = "" to some paste commands)
* Added  AICtweedie()  to compute AIC for Tweedie glms

tweedie 2.0.0 (Release date: 2009-08-10)
==============

Changes:

* Made  method = "inversion" the default (was "series")
* Slightly changed default output (added sep = "" to some paste commands)
* An error introduced earlier (unsure when exactly; prob v 1.6.1)
* In trying to identify and fix the error, tidied some of the FORTRAN code
* Made  do.smooth = TRUE  the default (was FALSE)
* If  p.vec  is not supplied,  tweedie.profile()  makes a sensible guess
* Made  verbose = FALSE the default  (was TRUE)
* Made minor changes to output when  verbose = FALSE

tweedie 1.6.8 (Release date: 2009-07-18)
==============

Changes:

* Correct error in dtweedie.saddle() when y = 0 for 1 < p < 2 (Thanks to Glenn Meyers <GMeyers@iso.com>)

tweedie 1.6.7 (Release date: 2009--06-30)
==============

Changes:

* Correct error in rtweedie() when power = 1 and phi \ne 1 (Thanks to Frederic Gosselin <frederic.gosselin@cemagref.fr>)

tweedie 1.6.6 (Release date: 2009-09-19)
==============

Changes:

* Corrected specification of GPL; no changes to functionality, code or documentation
* In dtweedie(), added power == 3 (inverse Gaussian) as a special case
* Minor change in documentation
* Added  data, weights, and  offset  input arguments
* Fixed a bug where "=>" was used rather than ">="
* Fixed a bug in tweedie.profile(): exact zeros still allowed p = 2 (but not p > 2)

tweedie 1.6.2 (Release date: 2009-02-16)
==============

Changes:

* Fixed a bug in  ptweedie.inversion()  where variable  cdf  was undefined
* Fixed a bug in  qtweedie()  where any(power) < 1  should have been  any(power < 1)

tweedie 1.6.1 (Release date: 2009-02-06)
==============

Changes:

* Fixed error in tweedie.f()  (twice) where a float was used in a do-loop
* Fixed a bug in ptweedie(): It now returns 0 when y < 0 (rather than reporting an error)
* Fixed a bug in ptweedie.inversion(): y.len was undefined
* Fixed a bug in dtweedie.inversion(): if( any(phi < 0) ) written as if (any(phi) < 0 )
* Built for R version 2.6.0
* Changed the tweedie.plot() to produce no default x and y labels

tweedie 1.5.2 (Release date: 2007-09-01)
==============

Changes:

* Built for R version 2.5.1
* Removed the inversion computation from the example for  dtweedie(); it was very slow and hence caused problems

tweedie 1.5.1 (Release date: 2007-05-03)
==============

* Built for R version 2.5.0
* Fixed citation details
* Added note in dtweedie() help page indicating where the methods are defined (Dunn and Smyth, 2007).
* Added tweedie-package.Rd
* Fixed a small bug in  dtweedie.inversion() (replaced  p  with   power)
* Built for R version 2.5
* Made change to use of pmatch() in function  tweedie()
* Swapped the order of methods 1 and 2 in the call of tweedie.inversion  to make consistent with a change in the associated paper Dunn and Smyth.
* Fixed an uninitialized variable error in the FORTRAN code which sometimes caused the inversion procedure to fail.
* Other minor changes (fixing spelling errors, updated Dunn and Smyth reference; etc.)
* Some very minor changes (some cosmetic) in the code
* Made the code more robust, exiting gracefully when Inf and NA appear in the computed log-likelihood
* Slightly changed the output plot in tweedie.profile(), to reflect that some value of  p.vec  may not have a corresponding likelihood computed accurately (or, indeed, at all)
* Added new function   tweedie.plot()  to make it easier to plot Tweedie densities

tweedie 1.2 (Release date: 2006-09-05)
==============

Changes:

* Built for R version 2.0.0
* Major fix: The mle estimate of phi was incorrect in tweedie.profile() (thanks Sarah Lennox, QDPI&F) * Some better warnings in  tweedie.profile()
* tweedie.profile() now allows the actual computed points to be added to the plot using the  do.points  option (default is TRUE)
* The tweedie.profile() function now smooth and computes confidence intervals even when the computed likelihood values contain Infs or NAs.

tweedie 1.01 (Release date: 2004-09-30)
==============

Changes:

* Built for R version 1.9
* Fixed a bug in rtweedie(): it failed for p = 2 (thanks Sarah Lennox, QDPI&F)
* Fixed a bug in qtweedie(): it failed in some cases (thanks Gordon Smyth, WEHI)
* Added a proper  .First.lib  function (!) (thanks to James Wettenhall and Gordon Smyth, WEHI)
* Other minor bugs fixed
* The  tweedie()  family function removed, and package made dependent on the  statmod  package
* Updated for R version 1.9
