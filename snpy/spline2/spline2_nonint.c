/* Hacked version of spine2 that puts the algorithm in its own function that can be called from 
 * other programs.  Meant to be run non-interactively (I think) and memory allocation done by
 * calling program.*/
#define  VERSION "6.1"
/*
 *   Spline2/Supersp
 *
 *      Public version (5.02) of 14 may 1997
 *      Update to version (5.03), 2 october 2001
 *               Changed REJLEV into a command-line option -L for fine tuning.
 *      Update to version (5.1), 16 october 2001
 *               Introduced -K option to specify assumed correlation length.
 *      Update to version (5.11), 26 october 2001
 *               Introduced -e option to allow "equal-information" splines,
 *               i.e. splines with non-optimized breakpoints. By default,
 *               these are now no longer allowed.
 *      Update to version (6.0), 22 january 2002
 *               Automatic search for autocorrelation in noise: -s option.
 *               New output file: splacf.
 *               New output file: splksi.
 *               Option -e inverted: forbid "equal-information" splines.
 *               Option -F now redundant. Full output is now standard.
 *               File splstat is nonsense for -s.
 *               Function lgamma (from Marcus Karolewski) added to code.
 *               Output file splsum contains a five-lines summary of results.
 *      Update to version (6.0a), 25 january 2002
 *       Ksi steps over a finer mesh as long as it is < dxav.
 *      Update to version (6.1), 14 october 2004
 *               No functional change. Mainly: fixed array dimensions and
 *               array overflow checking, and improved help text.
 *
 ****   Info:
 *
 *      For information about splines and the algorithm-kernels used, see
 *   (Ref. 1) Carl de Boor, "A Practical Guide to Splines",
 *       (Springer, 1978).
 *
 *   The curve-fitting algorithms and strategies employed are described in:
 *      (Ref. 2) Barend J. Thijsse, Mark A. Hollanders, and J. Hendrikse,
 *               "A practical algorithm for least-squares spline approximation
 *               of data containing noise", Computers In Physics,
 *               Jul/Aug (1998).
 *      A paper explaining the autocorrelation search method is in preparation.
 *
 *   Authors and owners of this program are:
 *      Barend J. Thijsse and Mark A. Hollanders
 *              Virtual Materials Laboratory
 *      Department of Materials Science and Engineering
 *      Delft University of Technology
 *      Rotterdamseweg 137, 2628 AL  Delft
 *      Netherlands
 *         Phone: +31 15 278 2221
 *         Fax:   +31 15 278 6730
 *         E-mail: B.J.Thijsse@tnw.tudelft.nl
 *         URL: http://dutsm183.stm.tudelft.nl/fcm
 *
 *   Additional Acknowledgments to:
 *      - Mathsoft Engineering & Education, Inc.
 *      - Jan Hendrikse
 *      - Marcus Karolewski
 *
 *****************************************************************************
 *
 *  THE DOCUMENTATION SECTION THAT FOLLOWS IS NOW OBSOLETE (22 Jan 2002)
 *
 *  PLEASE REFER TO THE SEPARATE "SPLINE2 V6.0 TUTORIAL AND USER MANUAL"
 *
 *
 *****************************************************************************
 *
 *  (nevertheless, please leave the documentation in this file, so that
 *   detailed reference to earlier versions remains possible)
 *
 * ---------------------------------------------------------------------------
 *  Full manual pages start after this short introductory section
 * ---------------------------------------------------------------------------
 *
 *   This file "spline2.c" contains the C source-code of "spline2", starting
 *   with a few Unix-style user manual pages. The program is coded in
 *   traditional C (not in ANSI-C), so essentially nothing is ever checked.
 *   However, the authors are reasonably confident that the program
 *   works correctly.
 *   There is also a second file "evalsp.c", containing the program
 *   "evalsp", which   enables you to actually use the results of "supersp".
 *   This file starts with manual pages as well.
 *   Three more files, "example1", "example2", and "tiheat" contain noisy
 *  data to play around with. These data are those  used in Ref. 2, below.
 *
 *
 ****   To be up and running in a few seconds:
 *
 *      Use your favorite C-compiler to compile the programs,
 *
 *      cc -o spline2 spline2.c -lm
 *      cc -o evalsp evalsp.c -lm
 *
 *   or by commands that look like this, and include your favorite
 *      optimization flags. Next, try out:
 *
 *      spline2 example1
 *
 *   and answer 's' to the question about choosing the number of spline-
 *   intervals l. A good fit should have a "dws"-value of about 2.
 *   ("dws" means Durbin-Watson Statistic; it is explained later on.)
 *   Other values for l can also be tried. If l is taken too low,
 *   you'll end up with systematic fit-errors (dws < 2). If l is
 *   taken too high, you'll be overfitting your data, but STATISTICALLY
 *   the fit is still good: you simply have allowed your fitting curve to
 *   be 'too flexible' to be realistic.
 *   The program exits when you answer 'y' to the "Done (y/n)?"-question.
 *
 *      The normal output of "spline2" is only a file named "splrep", which
 *      contains a bunch of numbers representing the fitted spline.
 *      By themselves these numbers will do you no good; you have to do
 *      some post-processing. So, for example, type
 *
 *      evalsp splrep -e -a -v 10000
 *
 *   to obtain the minima and maxima of the spline (-e for extremes), the
 *   area under the spline (-a for area), and the x-value(s) where the
 *      spline has a value of 10000 (-v for value). These results are written
 *      on your terminal (stdout).
 *
 *   The datapoints x, y, and spline values S(x) in the datapoints
 *   are found in the columns 1, 2, and 3 of the file "splres" (spline
 *   results). This file is only produced if you add -F to the command
 *   line, i.e. if you type
 *
 *      spline2 example1 -F
 *
 *   Value pairs x, S(x) in 2001 equidistant points (2000 x-intervals)
 *   are found in the columns 1 and 2 of the file "evlres" (evaluation
 *   results). This file is produced each time you run "evalsp" after
 *   running "spline2", even in its default form
 *
 *      evalsp splrep
 *
 *   After playing around with this for a while, you'll know enough to
 *   try out the other options in the manual section below. For example,
 *   compare the results of
 *
 *      spline2 example2
 *
 *   with those of
 *
 *      spline2 example2 -n 3
 *
 *   Good luck.
 *
 *----------------------------------------------------------------------------
 *            MANUAL PAGES
 *----------------------------------------------------------------------------
 *
 *      spline2 - A program that tries to find the best least-squares spline-
 *        approximation to a set of datapoints in a more or less
 *        user-friendly way.
 *
 *                Note: the highly unimaginative name of this program is
 *        supposed to be a combination of "spline" and "least squares".
 *        Back in Delft we call the program "superspline" or even
 *                "hyperspline", but to the outside world we wanted to appear
 *        a bit more modest.
 *
 ****   Syntax:
 *
 *   spline2 [filename] [-F] [-k order] [-x xbegin xend] [-X Xbegin Xend]
 *      [-a value] [-r value] [-n spacing] [-o lrev] [-q] [-l]
 *              [-L rejlev] [-K ksi] [-e] [-s]
 *
 *   The input file should be a "2-column" or a "3-column" file,
 *   i.e. each line should contain two data:
 *      x y
 *      or three data:
 *      x y sigma
 *   where sigma is (an estimate of) the standard deviation of the
 *   random error in y. The data on a line should be separated by
 *   "white space" - one or more spaces or tabs.
 *
 *   In case of a 2-column file, all sigmas are set internally equal
 *   to 1 (unless another value is prescribed via the options -a or -r,
 *   see below).
 *
 *   The x-values should be listed in increasing order.
 *   The file should not contain more than 10000 datapoints (or else
 *      change the line "#define NDATMAX 10000" in the source code).
 *
 ****   Program summary:
 *
 *      See next section.
 *
 ****   Command line options:
 *
 *   Note: a short summary of the command line options can be viewed
 *      by just typing "spline2" at the prompt.
 *
 *   [filename]
 *   If no filename is given, the program expects input from "stdin",
 *   i.e. data arriving through a pipe or via the keyboard.
 *
 *   [-s]
 *   Searches automatically for the best spline AND for any possible
 *   autocorrelation in the noise. The -s option is meant to produce
 *   the "optimal" results of the program. If it is used, the -n and
 *   -K options are no longer necessary.
 *   For the autocorrelation test the autocorrelation length ksi runs
 *   from 0 to a certain value (ksiend) to find the optimum. For each 
 *   trial-ksi, the "data spacing index" n covers the range from 1 to 
 *   a certain ksi-dependent upper bound (nmax). Ksiend and nmax are
 *  given reasonable values in the code.
 *   NOTE: Without -s, -n, and -K, the program assumes that there is
 *   no autocorrelation in the noise, i.e., that the noise is "white". 
 *
 *   [-n spacing]
 *  Flag to restrict the test on autocorrelation to data points
 *   separated by a fixed difference in index. The effect is that 
 *   the dws parameter, by which the quality of the spline fit is
 *   measured, is calculated for further-than-successive
 *   weighted fit residuals, i.e. using terms like (d[i+3] - d[i])**2
 *   (if spacing=3) rather than (d[i+1] - d[i])**2 (if spacing=1). The
 *   default value for "spacing" is 1.
 *   If the -K option is not used, the -n option implies the assumption
 *  ksi = 0 (ksi is the assumed autocorrelation length).
 *
 *   [-K ksi]
 *   The -K option allows the user to specify an assumed (or known)
 *  autocorrelation length. The default value is ksi = 0 (i.e.,
 *   uncorrrelated -"white"- noise). As for the -s option, the -K
 *   option implies that n runs from 1 to a certain ksi-dependent
 *   upper bound nmax. See code for how nmax depends on ksi.
 *
 *   [-x xbegin xend]
 *   [-X Xbegin Xend]
 *   Normally, the spline-approximation is calculated from the point
 *   with the lowest x-value to the point with the highest x-value.
 *      Extrapolation outside the datarange is impossible.
 *   If an approximation of only a part of the dataset is wanted,
 *   use the option "-x xbegin xend"; if a part of the dataset
 *   should be cut out, use "-X Xbegin Xend". The options
 *   -x and -X can be used simultaneously: this is convenient
 *   for instance for background determination. As an example,
 *      spline2 example1 -x 320 1950 -X 400 1800
 *      fits a spline only to the data between x=320 and x=400
 *   and those between x=1800 and x=1950. Of course this spline is
 *      defined in the cut-out range 400 < x < 1800 as well. This method
 *      can be fruitfully used to determine and subtract background
 *      trends.
 *
 *   [-a value]
 *   If the sigmas are exactly known as absolute values, and/or you want
 *   them to be considered strictly true, use this "absolute freeze" option.
 *   There are 3 possibilities:
 *      value = 0    sigmas are used as specified in the 3rd column
 *              of the input file.
 *      value > 0    sigma=value is used for all y.
 *      value = -n   input y-values are assumed to have n significant
 *              digits. (data values equal to zero get the smallest
 *                   sigma found; this was implemented in version 4.2)
 *   The -a option has the effect that the Chi-squared test
 *   is used instead of the Durbin-Watson test (see section "Program
 *      summary", below).
 *
 *   [-r value]
 *   If the sigmas are exactly known as relative values, and/or you want
 *   them to be considered strictly true, use this "relative freeze" option.
 *   There are 2 possibilities:
 *      value > 0    sigma=value*|y| is used for all y.
 *      value < 0    sigma=|value|*sqrt(|y|) is used for all y.
 *   In both cases the program crashes if a value y=0 is found, because
 *   that would make the weight factor for that point infinitely large.
 *   The -r option has the effect that the Chi-squared test
 *   is used instead of the Durbin-Watson test (see section "Program
 *   Summary", below).
 *
 *   [-e]
 *   Forbid splines with non-optimized breakpoints.
 *
 *   [-L rejlev]
 *   Specifying another rejection level than the default value (0.05) may
 *   be useful for finetuning. A smaller value means: be less critical.
 *
 *   [-q]
 *   The -q option ("quick" or "quiet") is described in a separate
 *   paragraph further down.
 *
 *   [-l]
 *   Sometimes a transformation of the x-values to log10(x) is convenient.
 *   This can be achieved with the -l option.
 *
 *   [-k order]
 *   The default spline order (= polynomial degree + 1) is 4 (cubic).
 *   Using the -k option, any order between 1 and 23 can be selected.
 *
 *   [-o lrev]
 *   The -o option can be used to specify the number of intervals l at which
 *   the "search-direction" reverses and knot-optimization starts ("lrev"
 *   stands for l-reversal). Normally this is done automatically. This
 *   option allows you to interfere with the search-strategy (see below).
 *
 *   [-f]
 *   This option is exactly equal to the -a option. We found -a easier to
 *   remember, that's all.
 *
 *   [-F]
 *   The -F option is used to generate a full collection of
 *   output-files, instead of just "splrep" (see section "Output", below).
 *  As from Version 6.0, the -F option is no longer necessary. Full
 *  output is always generated.
 *
 *
 ****   Program summary:
 *
 *   INTRODUCTION
 *   The important point to realize when one tries to fit a spline
 *   funtion to datapoints is that one can obtain a fit as close to
 *   the data as one wants, simply by adding more and more breakpoints
 *      to the spline funcion. However, what one really wants is a
 *      smooth curve, flexible enough to capture the (unknown) functional
 *      relationship underlying the data, yet smooth enough NOT to follow
 *      the noise component in the data due to e.g. measurement errors.
 *   In other words, one wants an 'adequate' fit, not a precise fit.
 *
 *   The problem of separating the noise from the underlying trend
 *   boils down to two issues:
 *         - Which statistic should be used to assess the 'adequacy'
 *              of the fit, and which statistical test should be
 *              performed on it?
 *         - Which search path should be followed (in 'knot-position space')
 *      so that the number and the distribution of the knots
 *      converge to their 'optimum' values?
 *      The algorithm described below provides solutions to these questions.
 *   After many years of experience with widely varying types of data we
 *      find that this algorithm works reliably and robustly. Almost always
 *      it ends up with the 'correct' spline without any user decisions.
 *   Only in exceptional cases 'manual' tuning of some parameters is
 *   necessary; some tips are given below in the section "User interference
 *      with the search strategy".
 *
 *   TEST STATISTICS USED
 *      The search strategy, to be outlined below, involves the computation
 *   of a considerable number of spline-approximations, each
 *   with a different set of breakpoints. The user does not see this.
 *
 *   For each spline-fit two test statistics, "rms" and "dws", are
 *      determined:
 *       rms = sqrt(D**2/(N-m))
 *   and
 *       dws = G**2/D**2,
 *      where
 *              D**2 = SUM(i=1..N){ (d[i])**2 },
 *              G**2 = SUM(i=1..N-1){ (d[i+1]-d[i])**2 },
 *   in which
 *              d[i] = (y[i]-S(x[i])/sigma[i].
 *   N is the number of datapoints participating in the fit, and N-m is
 *   the number of degrees of freedom. The program uses m = l + [k-1],
 *      where l is the number of intervals between the breakpoints of the
 *   spline and k-1 the polynomial degree of the spline. S(x) is the spline
 *   function itself and sigma[i] is the RMS estimate of the error component
 *   of y[i]. If the sigmas are not given to the program (in the 3rd
 *      data column or via the -a or -r options), they all get the value 1.
 *
 *   If the sigmas are frozen (via the -a or -r options), the statistical
 *   acceptance of a fit is evaluated by the well-known Chi-squared
 *   test on D**2. This is not a common situation. Normally one does not
 *      know the experimental errors so precisely that the Chi-squared
 *   test will perform well. This test is so powerful that even a
 *   modest under- or overestimation of the noise leads to unwanted
 *   over- or underfitting of the data.
 *
 *   If the sigmas are not frozen, which is the normal situation, the
 *   statistical acceptance of a fit is decided upon by the so-called
 *   Durbin-Watson test on dws, using the beta-function approximations
 *   for the two extreme theoretical distributions of dws. A statistically
 *   acceptable fit has a dws-value of around 2. Smaller values indicate
 *   the presence of systematic fit errors, larger values are just lucky
 *      (or indicate that the data have been tampered with).
 *   More on the Durbin-Watson statistic dws and on the associated
 *      theoretical distributions can be found in J. Durbin and G.S. Watson,
 *   Biometrika 37 (1950) 409, Biometrika 38 (1951) 159, and Biometrika 58
 *      (1971) 1, in Chapter IX of Mark A. Hollanders, Thesis TU-Delft
 *      (1990), and in Ref. 2, below.
 *
 *   The crucial property of dws is that any unknown common multiplicative
 *   factor in the sigmas cancels out, because the sigmas appear in the
 *   numerator and in the denominator. Experience shows that even if the
 *   sigmas, or an unknown part of them, vary from datapoint to datapoint,
 *   the value of dws is not too sensitive to this variation (if, at
 *   least, it is not varying too wildly) to make the program find
 *   unsatisfactory splines. This makes dws a far better statistic than
 *      chi**2. In some ways, dws can be compared to the number of sign
 *      changes in the sequence of fit-residuals.
 *
 *   In all cases a significance level of 5% is used, but another value
 *      may be used by using the -L option.
 *
 *   SERIAL CORRELATION
 *   The dws test is crucially dependent on the assumption that the
 *   noise values of successive datapoints are statistically independent.
 *   This is in real situations often not the case, since many datasets
 *   contain some sort of serial correlation, for example as a result of
 *   smoothing, filtering, or limited experimental resolution along
 *   the x-axis (these three are in fact the same). Jan Hendrikse has
 *   shown that there is an elegant way around this, by letting dws
 *   'skip over' the correlation 'length' of the data. The generalized
 *   form of the dws parameter therefore becomes
 *       dws = A*(G**2/D**2),
 *      with
 *              D**2 = SUM(i=1..N) { (d[i])**2 },
 *              G**2 = SUM(i=1..N-spacing) { (d[i+spacing]-d[i])**2 },
 *         A = (N-1)/(N-spacing),
 *   where "spacing" can be specified using the -n option. The default
 *   value is 1, the classical Durbin-Watson value.
 *
 *   Note: the sum G**2 is not cyclic over the data (as probably would be
 *   better), so it is scaled up by (N-1)/(N-spacing) to account for the
 *   loss of terms in the numerator. This is a very crude way to be able
 *   to keep using the same statistical percentages.
 *
 *   Note: In version 5.1 a possibly better way to handle correlations
 *      has been implemented. It is based on the following idea.
 *   When correlations are present, the autocorrelation function of
 *      the fit residuals should no longer be zero. We assume here that
 *   the residuals are correlated according to
 *     <d(x)d(x+r)> = <d(x)^2> exp(-r/ksi),
 *      where r and ksi are distances along the x-axis. Such a model
 *   (approximately) arises, for example, when uncorrelated data
 *   are averaged over a certain range before spline2 gets them.
 *   From this equation it is easy to derive that
 *     <(d[i+n]-d[i])^2> / <d[i]^2> = 2 - 2<exp(-(x[i+n]-x[i])/ksi)>.
 *      Since (N-1)/N times the left hand side is equal to the dws statistic
 *      in unmodified form, adding 2(N-1)/N <exp(-(x[i+n]-x[i])/ksi)> to it
 *   has the effect of retaining the expectation value of dws (2(N-1)/N)
 *      also in the presence of correlation.
 *      In the current implementation this is the only modification
 *   applied, i.e. the statistical tests to which dws is subjected
 *   are not modified (it is not known how to do this).
 *      Usage note: When the correct ksi is chosen, using the -n option
 *      with various spacings not lead to significantly different splines.
 *
 *   INTERNAL SEARCH STRATEGY
 *   The search strategy for the "best" spline-approximation consists of
 *   the following steps (again, l denotes the number of intervals,
 *      i.e. the number of separate polynomial pieces of the spline):
 *     1) First, spline fits are calculated while the number of intervals
 *        increases from l = 1 up to the lowest number of intervals where
 *        the fit is statistically accepted, or (with option -o) up to
 *        "lrev", which stands for l-reversal.
 *        If the Durbin-Watson test is on, there is a certain indecisive
 *        acceptance range, because there are in fact two different
 *        theoretical distribtuions involved. Fits inside the indecisive
 *        range are accepted in this stage.
 *           In this first stage the knot-positions are uniformly distributed,
 *        but uniformly with respect to the datapoints, not uniformly
 *        over the x-axis. Each spline interval is made to contain the same
 *        number of datapoints (or as close to it as possible).
 *        This approach is based on experimental efficiency: one should
 *        measure more datapoints in the ranges where the data vary strongly
 *           than in the ranges where the data follow a slow trend.
 *           (On the other hand, in many practical cases the experimentalist
 *        takes his data at constant increments of 'x'. So maybe we should
 *        have an option in the program to switch from 'equal information
 *        content' in the spline intervals to 'equal x-increment'.)
 *     2) Next, the program calculates a new series of splines, while it
 *        decreases the number of intervals and simultaneously optimizes
 *        the knot-positions according to the "newnot"-algorithm of De Boor.
 *        (To prevent singularities, all intervals are made to contain
 *        at least one datapoint at all times.)
 *     3) The "spline wizard" then suggests one (or more) splines
 *        as "best" if it meets the following criteria:
 *      a) No singularities were found in solving the fit-equations.
 *      b) The fit has passed the acceptance test; fits within the
 *         indecisive area of the Durbin-Watson test are initially
 *         rejected.
 *      c) The spline with the next lower value of l does NOT meet
 *         the criteria a) and b). This is to make sure
 *         that one always gets the fit with the lowest possible
 *         number of intervals (to avoid "overfitting").
 *        4) When no spline-fit lies in the indecisive Durbin-Watson region,
 *           "lrev" is increased, a new equidistant spline is calculated
 *           for l = lrev, and the "newnot"-sequence is repeated. (This
 *        step is not invoked when the -o option is used to fix "lrev".)
 *     5) If no acceptable fit is found, a warning appears and the fit with
 *        the lowest rms is suggested.
 *     6) The number of intervals, together with the values of rms and
 *        dws for the suggested fits are tabulated and you may choose one of
 *        them (or in fact any other if you like). Typing "s" at this
 *        point selects the best "suggested" spline. You'll get as many
 *        opportunities to try an l value as you need.
 *        Note that when you request a spline-fit with l < lrev,
 *        you always get "optimized" breakpoint locations. If l >= lrev,
 *        you get regular "equispaced" breakpoint locations. Hence, by
 *        using the option -o 1 on the command line, you can make sure
 *         that every spline is "equispaced"; by using -o 200 (or so), you
 *        can make sure that every spline is (in a way) "optimized".
 *     7) Finally, if you think you're satisfied, you can inform the spline
 *        wizard that you are done with her assistance, and the program
 *        exits. See under "Using evalsp", below, what to do next.
 *
 ****   User interference with the search strategy:
 *
 *   It sometimes happens that optimization of the knot positions starts
 *   too early (especially when there are only a few data points) or
 *   too late (when there are trends in the deviations which you consider
 *   to be noise, for instance when a 'slow' noise component is
 *   superimposed on a 'fast' one); the best thing to do then, is to use
 *   option "-o lrev" to force the start of knot-number reversal and
 *      optimization at a desired point l=lrev, or to freeze the sigmas
 *   at a desired level (-a or -r options).
 *
 ****   Output:
 *
 *   Always:
 *      File "splrep" contains the following information about
 *      the spline S(x) that was finally accepted (i.e., the last one):
 *      a) the spline order, b) the number of intervals l,
 *      c) the l+1 breakpoints and d) the l+k-1 beta-spline coefficients.
 *      These data form a complete representation of S(x), the
 *      so-called Beta-Spline (BS) representation. These data
 *      are intended to be used by the program "evalsp", with which you
 *      can use, analyze and manipulate the spline.
 *         (v5.1: a new file "splcall" is generated, containing info on
 *         how the program was called. Not meant for the normal user.)
 *
 *   Only if the -F option is used:
 *      File "splstat" contains 6 (or 7) data columns giving statistical
 *      information on all spline approximations in "plottable" order of
 *      calculation (including the ones you tried yourself):
 *      1) # of intervals, initially increasing, then, during knot
 *         optimization, decreasing,
 *      2) # of degrees of freedom,
 *      3) value of rms,
 *          3a) percentage point of rms (only when -a or -r option is used),
 *      4) value of dws,
 *      5) percentage point of lower distribution of dws,
 *      6) percentage point of upper distribution of dws.
 *      (Some data for the ones you tried yourself appear as "-999" in the
 *      file, because we have forgotten where to get the proper
 *      values from).
 *
 *   Only if the -F option is used:
 *      File "splres" contains 7 columns giving data & fit information
 *      for each data point:
 *      1) x (or log(x), when the -l option is used),
 *      2) y,
 *      3) value of final spline approximation in point x,
 *      4) weighted error,
 *      5) weight factor 1/sigma**2,
 *      6) first derivative of the final spline approximation,
 *      7) second derivative of the final spline apprximation.
 *
 ****   Using evalsp
 *
 *   Type "evalsp splrep [options]" (after running "spline2") in order to:
 *      1) Evaluate S(x) and all its derivatives in any set of
 *         x-values [no options, or -x xbegin xend xstep, or
 *         -f xfilename]; results go to file "evlres".
 *      2) Determine the maxima, minima, and inflection points of S(x)
 *         [options -e and -i].
 *      3) Determine the x-values in which S(x) has a certain value,
 *         e.g. for computing level crossing points [option -v value].
 *      4) Calculate the integral of S(x)dx from the first to the
 *         last data point [option -a].
 *      5) Calculate the Fourier transform of S(x) [option -F sbegin
 *         send sstep]; results go to file "evlfour".
 *      6) Estimate the x-values in which a single-peaked spline
 *         reaches half its peak height [option -h, with or without
 *         option -b backgroundvalue].
 *   For detailed information see the "evalsp" manual page.
 *
 ****   The -q option
 *
 *   The -q option is intended for those cases where one trusts
 *   spline2's suggestion for the best approximation, and wants the
 *   program to run silently and quickly. There is no output at all,
 *   except for the rms-, dws-, and l-values of the best spline (to stderr),
 *   and the BS-representation of the spline to stdout (it goes to the file
 *   "splrep" as well). In this way "spline2" can be simply used in pipes,
 *   for example:
 *
 *      generate_data | spline2 -q [options] | evalsp [options]
 *
 *   where "generate_data" denotes any process that writes a two-
 *   or three-column data set to stdout.
 *
 ****   Bugs and version information:
 *
 *   If the x-values of the data are not in increasing order, everything
 *   goes wrong!
 *
 *      Program versions lower than 5.02 sometimes followed different
 *   search paths through knot-space on different platforms.
 *   The reason was that up- or downscaling the number of knots
 *   to an 'exact' integer number (for instance, lnew = 0.95*40) led
 *   on some machines to a double-to-int truncation to the correct
 *   integer (here: 38), while on others to the next lower integer
 *   (here: 37). This 'bug' has been 'fixed' in version 5.02: all
 *      values are still truncated, unless they are within UPMARGIN
 *      of the 'exact' integer.
 *
 *   Don't expect miracles. The program is good, but not perfect.
 *   Be wise and check your spline by plotting it together with the data.
 *   The spline values are most conveniently found in columns 1 and 2 of
 *   the file "evlres" that is created by running "evalsp splrep".
 *
 *----------------------------------------------------------------------------
 *         END OF MANUAL PAGES
 *----------------------------------------------------------------------------
 */

/* ---INCLUDES--- */
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

/* ---DEFINES (NUMBERS)--- */
/* v6.1 Array dimensioning fixed */
/* --------------------------------------------------------------------------*/
#define NDATMAX 100000   /* Maximum number of data points (ntau=1..NDATMAX) */
#define LMAX 20000   /* Maximum number of intervals (l) */
#define KMAX 24      /* Maximum spline order (k) */
#define RMAX 4000   /* Maximum number of intermediate results 
            (Overflow of 'iresul' and 'resul' NOT TESTED) */
#define KSIMAX  200     /* Maximum number of ksi values tested */
#define ACFMAX 24   /* Maximum no. of autocorr.fct. points (njan=1..ACFMAX) */
#define ROOTMAX 100
/* --------------------------------------------------------------------------*/

#define NMAX (KMAX+LMAX-1)   /* Maximum spline dimension (n) */
#define NTMAX (2*NMAX)      /* Maximum number of knots */
#define NBRMAX (LMAX+1)      /* Maximum number of breakpoints */

#define UPFACTOR  1.1   /* Increment factor for l */
#define DOWNFACTOR 0.95 /* Decrement factor for l */
#define MAXHYPFRAC 0.5  /* Max l value (fraction of N) for hyperspline */
#define UPMARGIN 1e-7   /* Calculated (double) values of l this close
            under an integer are rounded upward instead
            of truncated downward (v. 5.02) */
#define ACFRANGE1 3.0   /* The autocorr.fct. of the residuals is tested for.. */
#define ACFRANGE2 3     /* ..njan = 1 to ACFRANGE1*(ksi/dxav)+ACFRANGE2 */

/* ---DEFINES (MACROS)--- */
#define TRUNC_SPEC(i,d) (((double)((i)+1)-(d)) < UPMARGIN ? ((i)+1) : (i))

/* ---ANNOUNCE FUNCTIONS---  */
int   l2sub(int *, int, int);
int newknt(double ara[], double dara[][NMAX+1],int *iptr, int *jptr, double arb[], int *kptr, 
      double darb[][LMAX+1]);
int l2knts(double *ara,int *iptr,int *jptr,double *arb,int *kptr);
int l2appr(double ara[],int *iptr,int *jptr,double dara[][NMAX+1],double arb[], double arc[]);
int bchfac(double dara[][NMAX+1],int *iptr,int *jptr,double ara[]);
int bchslv(double dara[][NMAX+1],int *iptr,int *jptr,double ara[]);
int bsplpp(double ara[],double arb[],int *iptr,int *jptr,double dara[][NMAX+1],double arc[],
      double darb[][NMAX+1],int *kptr);
int bsplvb(double ara[],int jhigh,int index, double *xptr,int *iptr,double arb[]);
int l2err(int prfun,int grfun);
int interv(double ara[],int *iptr,double *xptr,int *jptr,int *kptr);
int dwtest(double x,double *xptr,double *yptr);
int gser(double *xptr, double a, double x,double *yptr);
int gcf(double *xptr,double a,double x,double *yptr);
double ppvalu(double ara[],double dara[][NMAX+1],int *iptr,int *jptr,double *xptr,
      int jderiv);
double chitest(double chi,int fr);
double betai(double a,double b,double x);
double betacf(double a,double b,double x);
double acffunc(int index, double deltax,double ksi);
double lgamma(double xx);
double ppigr(double ara[],double dara[][NMAX+1],int *iptr,int *jptr,double *xptr,double *yptr);
double zbrent(double x1,double x2,int ind,double value);

/*----GLOBAL VARIABLES----*/
int    km1,      /* degree of polynomials */
       k,      /* order of polynomials */
       ntau,      /* number of datapoints */
       l, lnew,      /* number of intervals */
       ibeg,      /* starting index of knot optimization */
       n,      /* order of spline-approximation = k+l-1 */
       freed,      /* degrees of freedom left */
       nsing,      /* number of singularities */
       ihi, ilo = 1,   /* flags for function interv() */
       nl2sub,       /* index of highest-l non-optimized approximation */
       il2sub,      /* counter of and number of approximations */
       isug,      /* index of suggested approximation */
       iksisug,      /* ksi-index of suggested approximation */
       jbsp,      /* counter for function bsplvb() */
       lend, lbeg,   /* endpoints of interval for sigma determination */
       lfin,      /* final (?) choice */
       left,      /* counter for function newknt() */
       wind,      /* flag for indecisive are warning status */
       windbp,      /* flag for non-optimized breakpoints warning status */
       sigma,      /* if =0: sigmas specified in input-file, otherwise =1*/
       fullon,      /* flag, =1 if full output is required */
       fix,      /* flag, =1 if sigma is fixed */
       done,            /* flag, =1 if Done? has been answered by yes */
       njan,      /* durbin-watson distance index */
       nmin, nmax,   /* minimum and maximum durbin-watson distance index */
       acfindex,        /* determines type of assumed autocorr. function */
       iksi,            /* counter for ksi trials */
       nksi,      /* number of ksi trials */
       iresul[RMAX+1][3],   /* l and freed for all approximations */
       iresulksi[KSIMAX+1][7];   /* array to store results for ksi-trial */
double ppvalue,      /* return value of program ppvalu() */
       rmsmin, rmsmax,   /* minimum and maximum rms-value */
       taumin, taumax,   /* minimum and maximum y-value of datapoints */
       dermax,      /* maximum absolute value of the spline derivative */
       dxav,      /* average point-to-point spacing along x-axis */
       totalw,      /* sum of weight factors */
       rms,      /* root mean square error */
       dws,      /* durbin-watson statistic(s), in CIP paper called Q */
       ctest, ltest,   /* percentage points of various distributions */
       utest, test,
       fixval,      /* (if positive) value at which sigma is fixed */
       nonsig,      /* value of first nonsignificant digit range */
       ksi,             /* v5.1 assumed correlation length */
       acfdevmin,   /* minimum deviation from assumed acf */
       acfdevmax,   /* maximum deviation from assumed acf */
       acffit,      /* overall deviation from assumed acf */
       tau[NDATMAX+1],    /* x-values of input data */
       gtau[NDATMAX+1],   /* y-values of input data */
       weight[NDATMAX+1],   /* weight factors of input data */
       maxweight,       /* largest weight (only used for negative fixval) */
       ftau[NDATMAX+1],   /* spline approximation data S(x) at input x-values */
       q[NDATMAX+1],      /* normalized errors */
       deriv[NDATMAX+1],   /* first derivative of spline approximation */
       brek[NBRMAX+1],   /* breakpoints */
       bcoef[NMAX+1],   /* b-coefficients */
       coef[KMAX+1][NMAX+1],   /* pp-representation */
       resul[RMAX+1][8],   /* test statistics for all splines */
       resulksi[KSIMAX+1][6], /* array to store results for ksi-trial */
       acfa[ACFMAX+1],  /* values of actual autocorr.fct. of residuals */
       acfv[ACFMAX+1],  /* values of assumed autocorr.fct. of residuals */
       deltal[KMAX+1],   /* used in function bsplvb() */
       deltar[KMAX+1];   /* used in function bsplvb() */
char   Accept,      /* indicator if final result is satisfying */
       line[200],   /* line read from input */
       rmsversion[80];  /* scaling mode of "rms" variable */

/* Extras from evalsp */
int    hhon, nextremes, nomaxmax=1, firsttime=1;
double halfheight, extrsum, maxmax, zerolevel;


int spline2(double *x, double *y, double *w, int npoints,   /* input data */
       int order,           /* order of the spline */
       int acfsearch,       /* flag, =1 if autocorrelation is automatically searched */
       int acf_ind,          /* determines type of assumed autocorr. function */
       int xflag, double xbegin, double xend,   /* data limits for fitting */
       int Xflag, double Xbegin, double Xend,   /* exclusion limits */
       int logtrans,   /* flag for log10-transformation of x-scale */
       int nset, int ksiset,   /* flags to indicate if n or ksi are specified */
       int n_min, int n_max,   /* minimum and maximum durbin-watson distance index */
       double ksibegin,        /* start of automatic ksi-range */
       double ksiend,          /* ksi end */
       double rejlev,          /* Rejection level of statistical tests */
       int rel,      /* flag, =1 if sigma is relative to y */
       int fixed_sigma, double fixval,  /* flag : fixed sigmas and value of fixed sigmas */
       int lind, int lopt,     /* flag:  if lind, then use lopt initial knot points */
       int allownonopt,     /* flag for allowing non-optimized breakpoints */
       int interactive,   /* flag for the -q option */
       int verbose,  /* verbose output?*/
       /*** RETURN VALUES ***/
       double **tt,  /* The knot locatitons */
       int *n_knots,  /* The number of knots */
       double **cc,  /* The spline coefficients */
       int *kk,   /* The order of the spline */
       double *ret_rms,  /* the returned RMS */
       double *ret_dws,  /* the redturned DWS */
       int *ret_lfin,  /* final value of l */
       double *ret_ksi, /* Returned correlation length */
       double *ret_acffit /* Returned deviation from acf */
       ) {
       
   /* ---DECLARATIONS (note: auto-arrays can't be initialized)--- */
   register int i, j, ii;
   int          nn, ln, lmax;
   double       a, rmmin, dummy;
   char          dum;
   double       ksistep = 1.0;         /* ksi step */
   double       *dptr;  /* pointer to double array */

   acfindex = acf_ind;
   nmin = n_min;  nmax = n_max;
   k = order + 1;
   /* fullon = full_output; */
   fix = fixed_sigma;

   /* ---INITIALIZATION---*/
   done = 0;
   sigma = njan = 1;
   ksi = 0.0;
   iksi = 0;

   i = 1;
   ii = 0;
   for (ii = 0 ; ii < npoints ; ++ii) {   /* run through the input data */
      tau[i] = x[ii];
      gtau[i] = y[ii];
      if (i>NDATMAX) {
         fprintf(stderr,"Line %d: too many data. Increase NDATMAX\n", i);
         return(1);
      }
      /* weight[i] = 1.0/dy[ii]/dy[ii]; */
      weight[i] = w[ii]*w[ii];
      if (fix==1) {
         /* sigmas freezed: override assignments */
         if (rel==0) {
            /* absolute sigma freezed */
            if (fixval>0.0) weight[i] = 1.0/fixval/fixval;
            /* if (fixval==0.0) weight was already assigned */
            if (fixval<0.0) {
               if (gtau[i]==0.0) {
                 /* temporary value: */
                 weight[i]= -999.999;
               } else {
                 nonsig = log10(fabs(gtau[i]));
                 /* repaired 26-2-1996 */
                 /*nonsig = floor(nonsig)+fixval;*/
                 nonsig = floor(nonsig)+fixval+1.0;
                 nonsig = exp(2.302585*nonsig);
                 nonsig /= 3.464102;
                 weight[i] = 1.0/nonsig/nonsig;
                 if (weight[i]>maxweight)
                     maxweight=weight[i];
               }
            }
         }
         if (rel==1) {
            /* relative sigma freezed */
            if (fixval>0.0) {
               weight[i] = fixval*fabs(gtau[i]);
               weight[i] = 1.0/weight[i]/weight[i];
            }
            if (fixval<0.0) {
               weight[i] = -fixval*sqrt(fabs(gtau[i]));
               weight[i] = 1.0/weight[i]/weight[i];
            }
         }
      }
      /* weight assigned; accept point for fit? */
      a = tau[i];
      if ((xflag==0 && Xflag==0) ||
          (xflag==1 && Xflag==0 && a>=xbegin && a<=xend) ||
          (xflag==0 && Xflag==1 && (a>Xend || a<Xbegin)) ||
          (xflag==1 && Xflag==1 &&
          ((a>=xbegin && a<Xbegin) || (a>Xend && a<=xend))))
          i++;
   } /* end of read loop */

   ntau = i-1;

   /* If data have certain number of significant digits, assign
      largest weight to the data value(s) zero */
   if (fix==1 && rel==0 && fixval<0.0) {
      for (i=1; i<=ntau; i++) {
         if (weight[i]== -999.999) weight[i]=maxweight;
      }
   }

   if (lopt >= .8*ntau || lind == 0) lopt = .8*ntau-1;
   km1 = k-1;
   lmax = ntau-km1-1;
   if (lmax < lopt) lopt = lmax;
   for (i = 1; i <= ntau; i++) totalw += weight[i];
   taumax = taumin = gtau[1];
   for (i = 2; i <= ntau; i++) {
      if (gtau[i] < taumin) taumin = gtau[i];
      if (gtau[i] > taumax) taumax = gtau[i];
   }
   if (logtrans != 0) for(i = 1; i <= ntau; i++) {
      if (tau[i]>0.0) {
         tau[i] = log10(tau[i]);
      } else {
        fprintf(stderr,"X-datapoint %d: cannot take log of %g\n", i, tau[i]);
        return(1);
      }
   }

   /* reading done */

   /* ---COMPUTE AVERAGE DATA SPACING--- */   
   dxav = (tau[ntau]-tau[1])/(double)(ntau-1);

   if (verbose) {
      /* ---MESSAGES--- */
      fprintf(stderr,"............. spline2 version %s ...............\n", VERSION);
      if (logtrans==1) fprintf(stderr,"x -> log10(x)\n");
      if (xflag==1) fprintf(stderr,"X-range limited to (%g, %g)\n", xbegin, xend);
      if (Xflag==1) fprintf(stderr,"X-range (%g, %g) cut out\n", Xbegin, Xend);
      fprintf(stderr,"Number of datapoints participating in spline-fit = %d\n", ntau);
      if (lind == 1) fprintf(stderr,"Search-reversal point set at l=%d\n", lopt);
      fprintf(stderr,"Degree of spline = %d\n", k-1);
      if (fix==1) {
         fprintf(stderr,"Uncertainties in y frozen at ");
         if (rel==0) {
            if (fixval==0.0) fprintf(stderr,"values specified in input\n");
            if (fixval>0.0) fprintf(stderr,"%g\n", fixval);
            if (fixval<0.0) fprintf(stderr,"%g-th digit\n", -fixval+1);
         }
         if (rel==1) {
            if (fixval>0.0) fprintf(stderr,"%g*|y|\n", fixval);
            if (fixval<0.0) fprintf(stderr,"%g*sqrt(|y|)\n", -fixval);
         }
      }
      fprintf(stderr,"Spline approximations are tested according to the ");
      if (fix==0) {
         fprintf(stderr,"Durbin-Watson test\n");
      } else {
         fprintf(stderr,"Chi-square test\n");
      }
      if (acfsearch) fprintf(stderr,"Automatic search for autocorrelation in residuals\n");
      if (ksiset) fprintf(stderr,"Assumed autocorrelation length ksi = %g (fixed)\n", ksibegin);
      if (nset) {
         fprintf(stderr,"Autocorrelation tested only for datapoint index spacing = %d\n", nmin);
         if(!ksiset) fprintf(stderr,"Assumed autocorrelation length ksi = %g (fixed)\n", ksibegin);
      }
      if (!nset && !ksiset && !acfsearch) {
         fprintf(stderr,"Autocorrelation tested only for datapoint index spacing = %d\n", nmin);
         fprintf(stderr,"Assumed autocorrelation length ksi = %g (fixed)\n", ksibegin);
      }
      fprintf(stderr,"Assumed autocorrelation function of residuals: ");
      if (acfindex==1) fprintf(stderr,"exponential: exp(-x/ksi)\n");
      if (acfindex==2) fprintf(stderr,"gaussian: exp(-x^2/(2ksi^2))\n");
      if (acfindex==3) fprintf(stderr,"linear: 1-x/(2ksi)\n");
      if (acfindex==4) fprintf(stderr,"sinc: sin(2x/ksi)/(2x/ksi)\n");
      fprintf(stderr,"Average datapoint spacing <delta x> = %g\n", dxav);
      if (allownonopt == 1) {
         fprintf(stderr,"Splines with non-optimized breakpoints are allowed\n");
       } else {
         fprintf(stderr,"Only splines with optimized breakpoints are allowed\n");
       }
   }
         
   /* ---DETERMINE AUTOMATIC KSI-RANGE--- */
   /* v6.0 */
   if (acfsearch) {
      ksibegin = 0.0;
      ksistep = dxav/5.0;
      ksiend = 3.0*dxav;
      /* v6.1 */
      iksi = 0;
      /* array overflow will be tested on the fly */
   }

   /* ---MAIN ALGORITHM--- */

   ksi = ksibegin;
   if (!interactive && verbose) fprintf(stderr,"ksi-sweep: ");

   
   /* ----------------------------------------------------- */
   /* begin step 1: equal-information splines, increasing l */
   /* ----------------------------------------------------- */

   /* v6.0: return point for autocorrelation search */
   lbegin2:
   
   if (!interactive && verbose ) fprintf(stderr,"+");
   il2sub = 0;

   /* ---COMPUTE MAXIMUM N FOR ACF TESTING--- */
   /* v6.0 Combined check on correlation-corrected dws for n=1..nmax */
   if ( acfsearch || (ksiset && !nset) ) {
      nmin = 1;
      nmax = (int)(ACFRANGE1*ksi/dxav+0.5)+ACFRANGE2;
   }
      
   if (nmax > ACFMAX) {
      nmax = ACFMAX;
      if (verbose) printf("Preventing ACF size overflow by limiting nmax to %d\n", nmax);
   }

   lnew = 1;  lfin = 0;

   /* v5.0: return point for hyperspline extension */
   lbegin1:

   test = 0.0;
   while (test < rejlev && lnew <= lopt && lnew <= LMAX) {
      l2sub(&lnew,0,1);
      ln = dummy = UPFACTOR*lnew;
      ln = TRUNC_SPEC(ln,dummy);
      if (ln > lopt && lnew < lopt) ln = lopt;
      if (ln == lnew) {
         lnew += 1; 
      } else {
         lnew = ln;
      }
      if (lind == 1) {
         test = 0.0;
      } else {
         test = ((fix == 0) ? ltest : ctest);
      }
   }
   nl2sub = il2sub;
   lbeg = iresul[nl2sub][1];


   /* -------------------------------------------------- */
   /* begin step 2: knot-optimized splines, decreasing l */
   /* -------------------------------------------------- */

   lnew = dummy = DOWNFACTOR*lbeg;
   lnew = TRUNC_SPEC(lnew,dummy);
   while (lnew > 1 && lnew <= LMAX) {
      l2sub(&lnew,1,1);
      lnew = dummy = DOWNFACTOR*lnew;
      lnew = TRUNC_SPEC(lnew,dummy);
   }
   if (lbeg > 1) {
      /* not essential: just for generating the final lines
         in the "search-path-plot" */
      il2sub++;
      for (j = 1; j <= 2; j++) iresul[il2sub][j] = iresul[1][j];
      for (j = 1; j <= 7; j++) resul[il2sub][j] = resul[1][j];
   }

   /* ---------- */
   /* end step 2 */
   /* ---------- */

   /*  BEtween here */
   Accept = 'n';
   while (Accept == 'n') {
      /* ---------------------------------------------------------- */
      /* step 4: search for best spline in this decreasing-l series */
      /* ---------------------------------------------------------- */

      /* Print some messages */
      if (interactive && !done) {
         printf("'rms' means: (Fit-estimated y-uncertainty)");
         if (fix==0) {
            if (sigma==1) sprintf(rmsversion,".");
            if (sigma==0)
               sprintf(rmsversion," / (uncertainty given in input file).");
         }
         if (fix==1) {
            if (rel==0 && fixval>0.0) {
                sprintf(rmsversion," / %g.", fixval);
            } else {
               if (rel==0 && fixval==0.0) {
                  sprintf(rmsversion," / (uncertainty given in input file.)");
               } else {
                  sprintf(rmsversion," / (uncertainty given on command line.)");
               }
            }
         }
         printf("%s\n", rmsversion);
         printf("...........................................................\n");
         printf("Suggestion(s) for a good spline-fit:\n");
      }

      /* Determine best fit for the current decreasing-l optimization
         series */
      wind = windbp = nn = 0;
      if (fix == 0) {
         /* If errors are not fixed: dws test */
         /* check results so far, first from the conservative
            viewpoint (j=5), and if no acceptable spline
            can be found, also from the liberal viewpoint
            (j=4).
         */
         for (j = 5; j >= 4 && nn == 0; j--) {
            if (j == 4) wind = 1;
            /* only the most recent decreasing-l series */
            for (i = nl2sub; i < il2sub; i++) {
               /* suggest spline i when
                  no singularities were found AND
                  spline i is acceptable AND
                  spline i+1 is not acceptable
               */
               if (resul[i][1] > 0.0 &&
                   resul[i][j] > rejlev &&
                   resul[i+1][j] < rejlev) {
                  if (interactive && !done)
                      printf("rms = %g, dws = %g: l=%d\n",
                         resul[i][1], resul[i][3],
                         iresul[i][1]);
                  if (windbp==1) windbp = 0;
                  if (iresul[i][1]==lbeg) windbp = 1;
                  isug = i;
                  nn++;
               }
            }
            if (nl2sub == 1 && resul[i][j] > rejlev) {
                if (interactive && !done)
                   printf("rms = %g, dws = %g: l=1\n",
                          resul[il2sub][1], resul[il2sub][3]);
                if (iresul[i][1]==lbeg) windbp = 1;
                isug = il2sub;
                nn++;
             }
         }
      } else {   /* If errors are fixed: chi-squared test */
         for (i = nl2sub; i < il2sub; i++) {
            if (resul[i][1] > 0.0 && resul[i][2] > rejlev && resul[i+1][2] < rejlev) {
               if (interactive && !done)
                   printf("rms = %g, dws = %g: l=%d\n",
               resul[i][1], resul[i][3], iresul[i][1]);
               if (windbp==1) windbp = 0;
               if (iresul[i][1]==lbeg) windbp = 1;
               isug = i;
               nn++;
            }
         }
         if (nl2sub == 1) {
            if (interactive && !done)
               printf("rms = %g, dws = %g: l=1\n",
            resul[il2sub][1], resul[il2sub][3]);
            if (iresul[i][1]==lbeg) windbp = 1;
            isug = il2sub;
            nn++;
         }
      }

      /* What if no acceptable fit can be found? */
      if (nn == 0) {
         if (interactive && !done) {
             printf("No really good spline can be found...");
             printf(" You'll have to live with this one:\n");
         }
         rmmin = resul[nl2sub][1];
         j = nl2sub;
         for (i = nl2sub; i <= il2sub; i++) {
             if (resul[i][1] < rmmin) {
                rmmin = resul[i][1];
                j = i;
             }
         }
         if (interactive && !done)
             printf("rms = %g, dws = %g: l=%d\n", resul[j][1],
                    resul[j][3], iresul[j][1]);
         if (windbp==1) windbp = 0;
         if (iresul[i][1]==lbeg) windbp = 1;
         isug = j;
      }

      /* v5.0: "Hyperspline" extension starts here */

      /* Check warning status. Does lrev have to be increased? */
      if ( (wind == 1 || (!allownonopt && windbp == 1) ) && !done) {
         if (interactive && wind == 1)
            printf("Mild warning: all suggested fits are in the indecisive D-W area\n");
         if (interactive && windbp == 1)
            printf("Mild warning: suggested fit has non-optimized breakpoints \n");
         if (lind != 1) {
            lopt = dummy = UPFACTOR*lbeg;
            lopt = TRUNC_SPEC(lopt,dummy);
            if (lopt==lbeg) lopt++;
            if (lopt > MAXHYPFRAC*ntau) {
               if (interactive)
                  fprintf(stderr,"End of lrev increments:  will not increase lrev to %d\n", lopt);
            } else {
               lnew = lopt;
               if (interactive)
                  printf("Now increasing lrev from %d to %d\n", lbeg, lopt);
               goto lbegin1;
            }
         }
      }

      /* Here the best possible spline has been identified for the
         current value of ksi. */

      /* Now store parameters necessary for final re-calculation:
         for spline:
           lbeg      = starting value of l
           lfin      = final value of l
         for Q (dw-statistic):
           nmin      = lowest data index spacing for acf
           nmax      = highest data index spacing for acf
           ksi       = autocorrelation distance
           acfdevmin = minimum deviation from assumed acf
           acfdevmax = maximum deviation from assumed acf
         for diagnostics message:
           nn        = number of good fits found
           wind      = indicator for indecisive fit
      */

      iksi++;
      iresulksi[iksi][1] = lbeg;
      iresulksi[iksi][2] = iresul[isug][1];  /* lfin */
      iresulksi[iksi][3] = nn;
      iresulksi[iksi][4] = wind;
      iresulksi[iksi][5] = nmin;
      iresulksi[iksi][6] = nmax;
      resulksi[iksi][1]  = ksi;
      resulksi[iksi][2]  = resul[isug][6];   /* acfdevmin */
      resulksi[iksi][3]  = resul[isug][7];   /* acfdevmax */
      resulksi[iksi][4]  = resul[isug][1];   /* rms */
      resulksi[iksi][5]  = resul[isug][3];   /* dws */

      /* Next ksi */
      if (iksi > KSIMAX) {
         fprintf(stderr,
            "Warning: ksi testrange abnormally ended because of array overflow\n");
         fprintf(stderr,"Last ksi = %g\n", ksi);
      } else {
         if (ksi < ksiend) {
            /* V6.0a Try out a finer step for ksi < dxav */
            ksi += (ksi < dxav ? ksistep/2.0 : ksistep);
            if (interactive)
               printf("Now increasing ksi to %g\n", ksi);
            goto lbegin2;
         }
      }
      nksi = iksi;

      if (!interactive && verbose) fprintf(stderr," end sweep\n");

      /* All ksi's have been tested. Now find the one that
         yielded a spline fit with the residual acf that agrees
         best with the assumed acf */
         

      acffit = fabs(resulksi[1][2])+fabs(resulksi[1][3]);
      iksisug = 1;
      for(iksi=1; iksi<=nksi; iksi++) {
         dummy = fabs(resulksi[iksi][2])+fabs(resulksi[iksi][3]);
         if (dummy < acffit) {
            iksisug = iksi;
            acffit = dummy;
         }
      }
         
      /*( for(iksi=1; iksi<=nksi; iksi++) {
         fprintf(fpk, "%g %g\n", resulksi[iksi][1], resulksi[iksi][2]);
         }
      for(iksi=1; iksi<=nksi; iksi++) {
         fprintf(fpk, "%g %g\n", resulksi[iksi][1], resulksi[iksi][3]);
         }
      for(iksi=1; iksi<=nksi; iksi++) {
         dummy = fabs(resulksi[iksi][2])+fabs(resulksi[iksi][3]);
         fprintf(fpk, "%g %g\n", resulksi[iksi][1], dummy);
         } */


      /* The best ksi and the corresponding best spline are known,
         so that the final outcome can be presented */

      lbeg = iresulksi[iksisug][1];
      lfin = iresulksi[iksisug][2];
      nn = iresulksi[iksisug][3];
      wind = iresulksi[iksisug][4];
      nmin = iresulksi[iksisug][5];
      nmax = iresulksi[iksisug][6];
      ksi = resulksi[iksisug][1];
      acfdevmin = resulksi[iksisug][2];
      acfdevmax = resulksi[iksisug][3];
      rms = resulksi[iksisug][4];   
      dws = resulksi[iksisug][5];
/* Debug */   
      if (interactive) {
         iksi = iksisug;   
         fprintf(stdout,
            "<> ksi= %g small= %g big= %g rms= %g Q= %g ",
            resulksi[iksi][1],
            resulksi[iksi][2],
            resulksi[iksi][3],
            resulksi[iksi][4],
            resulksi[iksi][5]);
         fprintf(stdout,
            "lbeg= %d lfin= %d nn= %d wind= %d nmin= %d nmax= %d\n",
            iresulksi[iksi][1],
            iresulksi[iksi][2],
            iresulksi[iksi][3],
            iresulksi[iksi][4],
            iresulksi[iksi][5],
            iresulksi[iksi][6]);
      };


      /* In interactive mode: allow user a possibly different choice
         for l (this has become somewhat minimal, because ksi can
         not be hand-picked anymore; add this option? */
      if (interactive) {
         if (!done)
            printf("Non-optimized breakpoints will be used for l=%d and above\n", lbeg);
         printf("Choose the number of intervals for the spline ('s' for suggested)? l=");
         if (scanf("%d", &lfin)==0) lfin=iresulksi[iksisug][2];
      } else { /* non-interactive */
         /* old code
         lfin = iresul[isug][1];
         fprintf(stderr,"Spline-fit: rms=%g, dws=%g, l=%d (%s) ",
         resul[isug][1], resul[isug][3], iresul[isug][1],
         (iresul[isug][1]>=lbeg) ? "eqd" : "opt");
         if (nn == 0) fprintf(stderr,"(no good fit) ");
         if (wind == 1) fprintf(stderr,"(DW indecisive) ");
         fprintf(stderr,"\n");
         */
         /* new code */
         if (verbose) {
            fprintf(stderr,"Spline-fit: rms=%g, dws=%g, l=%d (%s) ksi=%g acffit=%g ",
               rms, dws, lfin, (lfin>=lbeg) ? "eqi" : "opt", ksi, acffit);
            if (nn == 0) fprintf(stderr,"(no good fit) ");
            if (wind == 1) fprintf(stderr,"(DW indecisive) ");
            fprintf(stderr,"\n");
            fprintf(stderr,"rms    %g\n", rms);
            fprintf(stderr,"dws    %g ", dws);
            if (nn == 0) fprintf(stderr,"(no good fit) ");
            if (wind == 1) fprintf(stderr,"(DW indecisive) ");
            fprintf(stderr,"\n");
         }
         /* INCLUDE IN FULL_OUTPUT? */
         /* fprintf(fpu,"l      %d (%s)\n", lfin,
                     (lfin>=lbeg) ? "eqi" : "opt");
         fprintf(fpu,"ksi    %g\n", ksi);
         fprintf(fpu,"acffit %g\n", acffit);*/
      }

      /* ---------------------------------------------- */
      /* begin step 5: re-calculate the selected spline */
      /* ---------------------------------------------- */


      if (lfin > lbeg && lfin <= LMAX) {
         l2sub(&lfin,0,0);
      } else {
         if (lbeg <= LMAX) {
            l2sub(&lbeg,0,0);
            lnew = dummy = DOWNFACTOR*lbeg;
            lnew = TRUNC_SPEC(lnew,dummy);
            while (lnew >= lfin && lnew <= LMAX) {
               l2sub(&lnew,1,0);
               lnew = dummy = DOWNFACTOR*lnew;
               lnew = TRUNC_SPEC(lnew,dummy);
            }
            if (l > lfin && lfin <= LMAX) l2sub(&lfin,1,0);
         }
      }

      if (interactive) {
         printf("rms = %g", rms);
         printf(", dws = %g.", dws);
         if (fix == 0) {
            test = ltest;
         } else {
            test = ctest;
         }
         if (test < rejlev) {
            printf(" <-- WARNING: not ok. ");
         } else {
            if (fix == 0 && utest < rejlev) {
               printf(" <-- Mild warning: indecisive. ");
            }
         }
         printf(" Done (y/n)? ");
         scanf("%c%c", &dum, &Accept);
         scanf("%c", &dum);
         done = 1;
      } else {
         Accept = 'y';
      }
   }
   /* And Here */

   /* Still needs to write to splsum if mode is interactive (added 20/2/02) */
   /* STUFF FOR FULL_OUTPUT?*/
   /* if (interactive) {
      fprintf(fpu,"rms    %g\n", rms);
      fprintf(fpu,"dws    %g ", dws);
      if (nn == 0) fprintf(fpu,"(no good fit) ");
      if (wind == 1) fprintf(fpu,"(DW indecisive) ");
      fprintf(fpu,"\n");
      fprintf(fpu,"l      %d (%s)\n", lfin,
              (lfin>=lbeg) ? "eqi" : "opt");
      fprintf(fpu,"ksi    %g\n", ksi);
      fprintf(fpu,"acffit %g\n", acffit);
      }   */
   *ret_rms = rms;
   *ret_dws = dws;
   *ret_lfin = lfin;
   *ret_ksi = ksi;
   *ret_acffit = acffit;

   /*---WRITE TO OUTPUT FILE(S)---*/
   l2err(1,0);
   *kk = k;
   *n_knots = l+1;
   /* fprintf(fps, "%d %d\n", k, l); 
   for (i = 1; i <= l+1; i++) fprintf(fps, "%g\n", brek[i]);
   for (i = 1; i <= l+km1; i++) fprintf(fps, "%g\n", bcoef[i]); */
   if ( ! (*tt = (double *) malloc((l+1)*sizeof(double)))) {
      fprintf(stderr, "ERROR:  bad malloc for *tt in SPLINE2\n");
      return(1);
   }
   if ( ! (*cc = (double *) malloc((l+km1)*sizeof(double)))) {
      fprintf(stderr, "ERROR:  bad malloc for *cc in SPLINE2\n");
      free(tt);
      return(1);
   }

   dptr = *tt;
   for (i = 1 ; i <= l+1 ; i++) dptr[i-1] = brek[i];
   dptr = *cc;
   for (i = 1 ; i <= l+km1 ; i++) dptr[i-1] = bcoef[i];

   /* ---HELPFUL MESSAGES--- */

   return(0);

} /* --end main-- */

int evalsp(double *in_x, int num_points, int in_k, int in_l, double *in_t, double *in_c, int deriv, 
      double **out_y) {
   /* Evaluate the spline.  *x is the input num_points x-values.  If values=1, then treat x-values as
    * y-values and do a reverse-interpolation.  If deriv>0, return the deriv-th derivative of the
    * interpolation (up to k-1) */
   /* ---DECLARATIONS (note: auto-arrays can't be initialized)--- */
   register int i;
   int          xflag, ntau, fileind, val, ext, inf, are, fou;
   int        background;
   double       xbegin, xend,t[LMAX+KMAX], b, scrtch[KMAX][NMAX+1], xv[NMAX];
   double     *dptr;

   /* ---INITIALIZATION---*/
     xflag = fileind = val = ext = inf = are = fou = 0;
     background = 0;

   /* ---PROCEED---*/

   k = in_k;
   l = in_l;
   for (i = 1; i <= l+1; i++) brek[i] = in_t[i-1];
   km1 = k-1;
   for (i = 1; i <= l+km1; i++) bcoef[i] = in_c[i-1];
   for (i = 1 ; i <= num_points ; ++i) xv[i] = in_x[i-1];
   ntau = num_points;
   xbegin = xv[1];
   xend = xv[ntau];

   l2knts(brek,&l,&k,t,&n);
   bsplpp(t,bcoef,&n,&k,scrtch,brek,coef,&l);

   if( ! (*out_y = (double *) malloc(num_points*sizeof(double)))) {
      fprintf(stderr, "Bad MALLOC for out_y in evalsp\n");
      return(1);
   }

   dptr = *out_y;
   for (i = 1; i <= ntau; i++)
   {
      b = xv[i];
      dptr[i-1] = ppvalu(brek,coef,&l,&k,&b,deriv);
   }

   return(0);
   } /* --end main-- */
            
int eval_extrema(int in_k, int in_l, double *in_t, double *in_c, double **xextr, double **yextr,
      int **signs, int *n_extr){
   register int i;
   int          xflag, ntau, fileind, val, ext, inf, are, fou;
   int        background;
   double       xbegin, xend, xstep, x, t[LMAX+KMAX], b, bb=0, a[KMAX],
           scrtch[KMAX][NMAX+1], a0=0, a1=0, a2=0, ymax;
   double xmaxs[ROOTMAX], ymaxs[ROOTMAX];
   int maxflags[ROOTMAX];
   double *dptr1, *dptr2;
   int *iptr;

   /* ---INITIALIZATION---*/
     xflag = fileind = val = ext = inf = are = fou = 0;
     background = 0;

   /* ---PROCEED---*/

   k = in_k;
   l = in_l;
   if (k <= 2)
   {
      fprintf(stderr, " spline order too low for finding extrema\n");
      return(1);
   }
   for (i = 1; i <= l+1; i++) brek[i] = in_t[i-1];
   km1 = k-1;
   for (i = 1; i <= l+km1; i++) bcoef[i] = in_c[i-1];
   ntau = 2001;
   xbegin = brek[1];
   xend = brek[l+1];
   xstep = (xend - xbegin)/2000;

   l2knts(brek,&l,&k,t,&n);
   bsplpp(t,bcoef,&n,&k,scrtch,brek,coef,&l);
   nextremes = 0;
   for (i = 1; i <= ntau; i++)
   {
      b = xbegin + (i-1)*xstep;
      a[0] = ppvalu(brek,coef,&l,&k,&b,0);
      a[1] = ppvalu(brek,coef,&l,&k,&b,1);
      a[2] = ppvalu(brek,coef,&l,&k,&b,2);

      if (i > 1 && a[1]*a1 <= 0. && a1 != 0 && a[2] != 0)
      {
         if (a[1] == 0.) x = b;
         else x = zbrent(bb,b,1,0);
         ymax = ppvalu(brek,coef,&l,&k,&x,0);
         if (a1 > 0.)
         {
         xmaxs[nextremes] = x;
         ymaxs[nextremes] = ymax;
         maxflags[nextremes] = -1;
         nextremes++;
         }
         else
         {
         xmaxs[nextremes] = x;
         ymaxs[nextremes] = ymax;
         maxflags[nextremes] = 1;
         nextremes++;
         }
      }
      a0 = a[0];
      a1 = a[1];
      a2 = a[2];
      bb = b;
   }
   *n_extr = nextremes;
   if (! (*xextr = (double *)malloc(*n_extr*sizeof(double)))) {
      fprintf(stderr, "Bad MALLOC for out_y in eval_extrema\n");
      return(1);
   }
   if (! (*yextr = (double *)malloc(*n_extr*sizeof(double)))) {
      fprintf(stderr, "Bad MALLOC for out_y in eval_extrema\n");
      return(1);
   }
   if (! (*signs = (int *)malloc(*n_extr*sizeof(int)))) {
      fprintf(stderr, "Bad MALLOC for out_y in eval_extrema\n");
      return(1);
   }
   dptr1 = *xextr;
   dptr2 = *yextr;
   iptr = *signs;
   for (i = 0 ; i < nextremes ; ++i) {
      dptr1[i] = xmaxs[i];
      dptr2[i] = ymaxs[i];
      iptr[i] = maxflags[i];
   }
   return(0);
}

int eval_inflect(int in_k, int in_l, double *in_t, double *in_c, double **out_x, double **out_y, 
      double **out_dy, int *n_inflect){
   register int i;
   int          xflag, ntau, fileind, val, ext, inf, are, fou;
   int        background;
   double       xbegin, xend, xstep, x, t[LMAX+KMAX], b, bb=0, a[KMAX],
           scrtch[KMAX][NMAX+1], a0=0, a1=0, a2=0, yder, ymax;
   double xmaxs[ROOTMAX], ymaxs[ROOTMAX], yders[ROOTMAX];
   double *dptr1,*dptr2, *dptr3;

   /* ---INITIALIZATION---*/
     xflag = fileind = val = ext = inf = are = fou = 0;
     background = 0;

   /* ---PROCEED---*/

   k = in_k;
   l = in_l;
   if (k <= 3)
   {
      fprintf(stderr, " spline order too low for inflection evaluation\n");
      return(1);
   }
   for (i = 1; i <= l+1; i++) brek[i] = in_t[i-1];
   km1 = k-1;
   for (i = 1; i <= l+km1; i++) bcoef[i] = in_c[i-1];
   ntau = 2001;
   xbegin = brek[1];
   xend = brek[l+1];
   xstep = (xend - xbegin)/2000.0;

   l2knts(brek,&l,&k,t,&n);
   bsplpp(t,bcoef,&n,&k,scrtch,brek,coef,&l);
   *n_inflect = 0;
   for (i = 1; i <= ntau; i++)
   {
      b = xbegin + (i-1)*xstep;
      a[0] = ppvalu(brek,coef,&l,&k,&b,0);
      a[1] = ppvalu(brek,coef,&l,&k,&b,1);
      a[2] = ppvalu(brek,coef,&l,&k,&b,2);
      a[3] = ppvalu(brek,coef,&l,&k,&b,3);

      if (i > 1 && a[2]*a2 <= 0. && a2 != 0 && a[3] != 0)
      {
         if (a[2] == 0.) x = b;
         else x = zbrent(bb,b,2,0);
         ymax = ppvalu(brek,coef,&l,&k,&x,0);
         yder = ppvalu(brek,coef,&l,&k,&x,1);
         if(yder>0.0){
            xmaxs[*n_inflect] = x;
            ymaxs[*n_inflect] = ymax;
            yders[*n_inflect] = yder;
            *n_inflect += 1;
         }
         if(yder<0.0){
            xmaxs[*n_inflect] = x;
            ymaxs[*n_inflect] = ymax;
            yders[*n_inflect] = yder;
            *n_inflect += 1;
         }
         if(yder==0.0){
            xmaxs[*n_inflect] = x;
            ymaxs[*n_inflect] = ymax;
            yders[*n_inflect] = yder;
            *n_inflect += 1;
         }
      }
      a0 = a[0];
      a1 = a[1];
      a2 = a[2];
      bb = b;
   }
   if (! (*out_x = (double *)malloc(*n_inflect*sizeof(double)))) {
      fprintf(stderr, "Bad MALLOC for out_x in eval_inflect\n");
      return(1);
   }
   if (! (*out_y = (double *)malloc(*n_inflect*sizeof(double)))) {
      fprintf(stderr, "Bad MALLOC for out_y in eval_inflect\n");
      return(1);
   }
   if (! (*out_dy = (double *)malloc(*n_inflect*sizeof(double)))) {
      fprintf(stderr, "Bad MALLOC for out_dy in eval_inflect\n");
      return(1);
   }
   dptr3 = *out_x;
   dptr1 = *out_y;
   dptr2 = *out_dy;
   for (i = 0 ; i < *n_inflect ; ++i ) {
      dptr3[i] = xmaxs[i];
      dptr1[i] = ymaxs[i];
      dptr2[i] = yders[i];
   }
   return(0);
}


double eval_integ(double x0, double x1, int in_k, int in_l, double *in_t, double *in_c) {
   register int i;
   int          xflag, fileind, val, ext, inf, are, fou;
   int        background;
   double       xbegin, xend, y, t[LMAX+KMAX],
           scrtch[KMAX][NMAX+1];

   /* ---INITIALIZATION---*/
     xflag = fileind = val = ext = inf = are = fou = 0;
     background = 0;

   /* ---PROCEED---*/

   k = in_k;
   l = in_l;
   for (i = 1; i <= l+1; i++) brek[i] = in_t[i-1];
   km1 = k-1;
   for (i = 1; i <= l+km1; i++) bcoef[i] = in_c[i-1];
   xbegin = brek[1];
   xend = brek[l+1];
   if (x0 < xbegin || x1 > xend) {
      fprintf(stderr, "Error:  integration limits beyond spline definition\n");
      return(1);
   }
   l2knts(brek,&l,&k,t,&n);
   bsplpp(t,bcoef,&n,&k,scrtch,brek,coef,&l);

   y = ppigr(brek,coef,&l,&k,&x0,&x1);
   return(y);
}


int eval_x(double value, int in_k, int in_l, double *in_t, double *in_c,double **ret_roots, int *n_roots) {
   register int i;
   int          ntau;
   double  xbegin, xend, xstep, x, t[LMAX+KMAX], b, bb=0, a[KMAX],
           scrtch[KMAX][NMAX+1], a0=0;
   double roots[ROOTMAX];
   double *dptr;

   /* ---PROCEED---*/

   k = in_k;
   l = in_l;
   if (k <= 1)
   {
      fprintf(stderr, "spline order too low for reverse evaluation\n");
      return(1);
   }
   for (i = 1; i <= l+1; i++) brek[i] = in_t[i-1];
   km1 = k-1;
   for (i = 1; i <= l+km1; i++) bcoef[i] = in_c[i-1];
   ntau = 2001;
   xbegin = brek[1];
   xend = brek[l+1];
   xstep = (xend - xbegin)/2000;

   l2knts(brek,&l,&k,t,&n);
   bsplpp(t,bcoef,&n,&k,scrtch,brek,coef,&l);
   *n_roots = 0;
   for (i = 1; i <= ntau; i++)
   {
      b = xbegin + (i-1)*xstep;
      a[0] = ppvalu(brek,coef,&l,&k,&b,0);

      if (i > 1 && (a[0]-value)*(a0-value) <= 0. && a0 != 0)
      {
         if (a[0]-value == 0.) x = b;
         else x = zbrent(bb,b,0,value);
         roots[*n_roots] = x;
         *n_roots += 1;
      }
      a0 = a[0];
      bb = b;
   }
   if (! (*ret_roots = (double *)malloc(*n_roots*sizeof(double)))) {
      fprintf(stderr, "Bad MALLOC for ret_roots in eval_x\n");
      return(1);
   }
   dptr = *ret_roots;
   for(i = 0 ; i < *n_roots ; ++i) {
      dptr[i] = roots[i];
   }
   return(0);
}



/*------------------------------------------------------------------*/
double ppigr(double ara[],double dara[][NMAX+1],int *iptr,int *jptr,double *xptr,double *yptr)

   /*  calculates integral from *xptr to *yptr  */
{
   int     i, j, right, ndummy;
   double  h, aa, bb, ppintgr;
   ppintgr = 0.;
   interv(ara,iptr,xptr,&left,&ndummy);
   h = *xptr-ara[left];
   bb = 0.;
   for (j = *jptr; j >= 1; j--) bb = (bb+dara[j][left])*h/j;
   interv(ara,iptr,yptr,&right,&ndummy);
   for (i = left; i < right; i++)
   {
      aa = 0.;
      h = ara[i+1]-ara[i];
      for (j = *jptr; j >= 1; j--)
      {
         aa = (aa+dara[j][i])*h/j;
      }
      ppintgr += aa;
   }
   h = *yptr-ara[right];
   aa = 0.;
   for (j = *jptr; j >= 1; j--) aa = (aa+dara[j][right])*h/j;
   ppintgr = ppintgr+aa-bb;
   return(ppintgr);
}


/*---------------------------------------------------------------------------*/
/*            Functions                 */
/*---------------------------------------------------------------------------*/
int l2sub(int *iptr,int iopt, int wr)

   /*   Calculates spline approximation to the data.
    *
    *   input:
    *      *iptr   = number of intervals;
    *       iopt   = 0: intervals will contain equal numbers of
    *              data points,
    *         = 1: optimize breakpoint locations with respect
    *         to the previous locations;
    *       wr     = 0: nothing is written to iresul[] and resul[],
    *              = 1: results are written to iresul[] and resul[].
    *
    *   output:
    *      iok = 1 if global Powell test is ok,
    *          = 0 otherwise,
    *      ipowt = sum of local Powell numbers,
    *      rms = weighted rms error for fit, should be equal to
    *            about 1.0 for well estimated data uncertainty,
    *      nsing = number of singularities found while solving the
    *            normal equations,
    *      ipo   = number of intervals where local powell test fails,
    *      brek = the new breakpoint sequence,
    *      coef = the matrix of the (i-1)th derivatives (to the right)
    *             in the jth breakpoint,
    *      ftau = the spline at the data points,
    *      q = normalized errors at the data points.
    */
{
   int   ip;
   double  coefg[3][LMAX+1], p[KMAX+1][NMAX+1], scrtch[NMAX+1],
      t[NTMAX+1], chi;
   ip = *iptr;
   if (iopt == 0) ip = -ip;
   newknt(brek,coef,&l,&k,scrtch,&ip,coefg);
   l2knts(scrtch,&ip,&k,t,&n);
   l2appr(t,&n,&k,p,scrtch,bcoef);
   bsplpp(t,bcoef,&n,&k,p,brek,coef,&l);
   l2err(0,1);
   dwtest(dws,&ltest,&utest);
   if (fix == 1)
   {
      chi = rms*rms*freed;
      ctest = chitest(chi,freed);
   }
   if (wr == 1)
   {
      il2sub++;
      iresul[il2sub][1] = ip;
      iresul[il2sub][2] = freed;
      if (il2sub == 1)
      {
         if (fix == 0) rmsmin = rms;
         else
         {
            rmsmin = 1.;
            rmsmax = 5.;
         }
      }
      if (rms < rmsmin) rmsmin = rms;
      if (rms > rmsmax) rmsmax = rms;
      if (nsing > 0) rms *= -1;
      resul[il2sub][1] = rms;
      if (fix == 1) resul[il2sub][2] = ctest;
      else resul[il2sub][2] = 0;
      resul[il2sub][3] = dws;
      resul[il2sub][4] = ltest;
      resul[il2sub][5] = utest;
      resul[il2sub][6] = acfdevmin;
      resul[il2sub][7] = acfdevmax;
   }
   return(0);
}

/*---------------------------------------------------------------------------*/
int newknt(double ara[], double dara[][NMAX+1],int *iptr, int *jptr, double arb[], int *kptr, 
      double darb[][LMAX+1])

   /*  calculates new ("best") set of breakpoints, no intervals with
    *  fewer then two data-points will be made
    */
{
   int   ip, jp, kp, i, ii, j, lleft, mflag, iend, idif;
   double  a, oneovk, dif, difprv, step, stepi;
   ip = *iptr;
   jp = *jptr;
   kp = *kptr;

   arb[1] = tau[1];
   arb[abs(kp)+1] = tau[ntau];
   if(kp < 0)
   {
      kp = -1*kp;
      /*$dir no_recurrence*/
      for (i = 2; i <= kp; i++)
      {
         a = ((double)(i-1))*(ntau-1)/kp+1;
         ii = a;
         arb[i] = tau[ii]+(a-ii)*(tau[ii+1]-tau[ii]);
      }
      *kptr = kp;
   }
   else
   {
      oneovk = 1./jp;
      darb[1][1] = 0.;
      difprv = fabs((dara[jp][2]-dara[jp][1])/(ara[3]-ara[1]));
      for (i = 2; i <= ip; i++)
      {
   dif = fabs((dara[jp][i]-dara[jp][i-1])/(ara[i+1]-ara[i-1]));
         darb[2][i-1] = pow((dif+difprv),oneovk);
      darb[1][i] = darb[1][i-1]+darb[2][i-1]*(ara[i]-ara[i-1]);
         difprv = dif;
      }
      darb[2][ip] = pow((2*difprv),oneovk);
step = (darb[1][ip]+darb[2][ip]*(ara[ip+1]-ara[ip]))/kp;
      if (step > 0.)
      {
         j = 1;
         for (i = 2; i <= kp; i++)
         {
            stepi = (i-1)*step;
            while (j != ip && stepi > darb[1][j+1]) j++;
            if (darb[2][j] != 0)
            {
            arb[i] = ara[j]+(stepi-darb[1][j])/darb[2][j];
            }
            else arb[i] = (ara[j]+ara[j+1])/2;
         }
      }
      else
      {
         step = (tau[ntau]-tau[1])/kp;
         /*$dir no_recurrence*/
         for (i = 2; i <= kp; i++)
         {
            arb[i] = tau[i]+(i-1)*step;
         }
      }
      lleft = 1;
      iend = idif = 0;
      for (i = 2; i <= kp+1; i++)
      {
         interv(tau,&ntau,&arb[i],&left,&mflag);
         if (left-lleft > 0) lleft = left;
         else
         {
            lleft += 1;
            if (lleft < ntau)
            {
               arb[i] = (tau[lleft]+tau[lleft+1])/2.;
            }
            else
            {
               iend++;
               arb[i] = tau[ntau];
            }
         }
      }
      if (iend != 0)
      {
         i = 1;
         lleft = ntau-1;
         while (idif < 1)
         {
         interv(tau,&ntau,&arb[kp+1-i],&left,&mflag);
            idif = lleft-left;
            if (idif < 1)
            {
               lleft -= 1;
            arb[kp+1-i] = (tau[lleft]+tau[lleft+1])/2.;
               i++;
            }
         }
      }
   }
   return(0);
}

/*---------------------------------------------------------------------------*/
int l2knts(double *ara,int *iptr,int *jptr,double *arb,int *kptr)

   /*  breakpoints to knots  */
{
   int   ip, jp, kp, i, kk;
   kk = km1;
   /*$dir no_recurrence*/
   for (i = 1; i <= kk; i++) arb[i] = ara[1];
   ip = *iptr;
   jp = *jptr;
   /*$dir no_recurrence*/
   for (i = 1; i <= ip; i++) arb[km1+i] = ara[i];
   kp = km1+ip;
   /*$dir no_recurrence*/
   for (i = 1; i <= jp; i++) arb[kp+i] = ara[ip+1];
   *kptr = kp;
   return(0);
}

/*---------------------------------------------------------------------------*/
int l2appr(double ara[],int *iptr,int *jptr,double dara[][NMAX+1],double arb[], double arc[])

   /*  calculates spline approximation to given data (tau,gtau)  */
{
   int   ip, jp, j, i, ll, mm, jj;
   double  dw, biatx[KMAX+1];
   ip = *iptr;
   jp = *jptr;
   ll = jp;
   for (j = 1; j <= ip; j++)
   {
      arc[j] = 0.;
      for (i = 1; i <= ll; i++) dara[i][j] = 0.;
   }
   left = jp;
   for (ll = 1; ll <= ntau; ll++)
   {
      while (left != ip && tau[ll] >= ara[left+1]) left++;
      bsplvb(ara,jp,1,&tau[ll],&left,biatx);
      for(mm = 1; mm <= jp; mm++)
      {
         dw = biatx[mm]*weight[ll];
         j = left-jp+mm;
         arc[j] += dw*gtau[ll];
         i = 1;
         for (jj = mm; jj <= jp; jj++)
         {
            dara[i][j] += biatx[jj]*dw;
            i++;
         }
      }
   }
   bchfac(dara,&jp,&ip,arb);
   bchslv(dara,&jp,&ip,arc);
   return(0);
}

/*---------------------------------------------------------------------------*/
int bchfac(double dara[][NMAX+1],int *iptr,int *jptr,double ara[])

   /*  constructs cholesky factorization  */
{
   int   ip, jp, i, j, imax, jmax, ii;
   double  ratio;
   ip = *iptr;
   jp = *jptr;
   nsing = 0;
   if (jp <= 1 && dara[1][1] > 0) dara[1][1] = 1/dara[1][1];
   else
   {
      /*$dir no_recurrence*/
      for (ii = 1; ii <= jp; ii++) ara[ii] = dara[1][ii];
      for (ii = 1; ii <= jp; ii++)
      {
         if (dara[1][ii]+ara[ii] <= ara[ii])
         {
            nsing++;
            for(j = 1; j <= ip; j++) dara[j][ii] = 0.;
         }
         else
         {
            dara[1][ii] = 1/dara[1][ii];
            if (ip-1 < jp-ii) imax = ip-1;
            else imax = jp-ii;
            if (imax > 0)
            {
               jmax = imax;
               for(i = 1; i <= imax; i++)
               {
               ratio = dara[i+1][ii]*dara[1][ii];
                  /*$dir no_recurrence*/
                  for(j = 1; j <= jmax; j++)
                  {
               dara[j][ii+i] -= dara[j+i][ii]*ratio;
                  }
                  jmax--;
                  dara[i+1][ii] = ratio;
               }
            }
         }
      }
   }
   return(0);
}

/*---------------------------------------------------------------------------*/
int bchslv(double dara[][NMAX+1],int *iptr,int *jptr,double ara[])

   /* solves a banded positive definite set of equations */
{
   int   ip, jp, kk, jmax, j, ii;
   double  nbndm1;
   ip = *iptr;
   jp = *jptr;
   if (jp <= 1) ara[1] *= dara[1][1];
   else
   {
      nbndm1 = ip-1;
      kk = jp;
      for(ii = 1; ii <= kk; ii++)
      {
         if (nbndm1 < jp-ii) jmax = nbndm1;
         else jmax = jp-ii;
         /*$dir no_recurrence*/
         if (jmax > 0) for (j = 1; j <= jmax; j++)
         {
            ara[j+ii] -= dara[j+1][ii]*ara[ii];
         }
      }
      for (ii = jp; ii >= 1; ii--)
      {
         ara[ii] *= dara[1][ii];
         if (nbndm1 < jp-ii) jmax = nbndm1;
         else jmax = jp-ii;
         if (jmax > 0) for (j = 1; j <= jmax; j++)
         {
            ara[ii] -= dara[j+1][ii]*ara[j+ii];
         }
      }
   }
   return(0);
}

/*---------------------------------------------------------------------------*/
int bsplpp(double ara[],double arb[],int *iptr,int *jptr,double dara[][NMAX+1],double arc[],
      double darb[][NMAX+1],int *kptr)

   /*  converts spline to piecewise polynomial representation  */
{
   int   ip, jp, lsofar, j, i, jp1, kmj;
   double  diff, sum, biatx[KMAX+1];
   ip = *iptr;
   jp = *jptr;
   arc[1] = ara[jp];
   lsofar = 0;
   for (left = jp; left <= ip; left++)
   {
      if(ara[left+1] != ara[left])
      {
         lsofar++;
         arc[lsofar+1] = ara[left+1];
         if (jp <= 1) darb[1][lsofar] = arb[left];
         else
         {
            /*$dir no_recurrence*/
            for (i = 1; i <= jp; i++)
            {
               dara[i][1] = arb[left-jp+i];
            }
            for (jp1 = 2; jp1 <= jp; jp1++)
            {
               j = jp1-1;
               kmj = k-j;
               /*$dir no_recurrence*/
               for(i = 1; i <= kmj; i++)
               {
               diff = ara[left+i]-ara[left+i-kmj];
                  if (diff > 0.)
                  {
         dara[i][jp1] = ((dara[i+1][j]-dara[i][j])/diff)*kmj;
                  }
               }
            }
            bsplvb(ara,1,1,&ara[left],&left,biatx);
            darb[jp][lsofar] = dara[1][jp];
            for(jp1 = 2; jp1 <= jp; jp1++)
            {
            bsplvb(ara,jp1,2,&ara[left],&left,biatx);
               kmj = k+1-jp1;
               sum = 0.;
               for(i = 1; i <=jp1; i++)
               {
                  sum += biatx[i]*dara[i][kmj];
                  darb[kmj][lsofar] = sum;
               }
            }
         }
      }
   }
   *kptr = lsofar;
   return(0);
}

/*---------------------------------------------------------------------------*/
int bsplvb(double ara[],int jhigh,int index, double *xptr,int *iptr,double arb[])

   /*  calculates all nonzero beta-splines at *xptr  */
{
   int   ip, jp1, i;
   double  xp, saved, term;
   ip = *iptr;
   xp = *xptr;
   if (index == 1)
   {
      jbsp = 1;
      arb[1] = 1.;
   }
   while (jbsp < jhigh)
   {
      jp1 = jbsp+1;
      deltar[jbsp] = ara[ip+jbsp]-xp;
      deltal[jbsp] = xp-ara[ip+1-jbsp];
      saved = 0.;
      for (i = 1; i <= jbsp; i++)
      {
         term = arb[i]/(deltar[i]+deltal[jp1-i]);
         arb[i] = saved+deltar[i]*term;
         saved = deltal[jp1-i]*term;
      }
      arb[jp1] = saved;
      jbsp++;
   }
   return(0);
}

/*---------------------------------------------------------------------------*/
int l2err(int prfun,int grfun)

   /*  calculates and prints details of spline fit;
    *
    *  input:
    *     prfun = 0: no final results written to splres,
    *          1: final results written to splres;
    *     grfun = 0: no first derivative calculated for plot,
    *           = 1: first derivative calculated for plot.
    */
{
   int   ll;
   double   err, errl2=0, a, dtau, ddtau; 
   double  errm, dipn2av, di2av, d2av, dipndiav, expav, sum1, sum2;
   double  dev, small, big;

   small = 1e30; big = -1e30;

   /* v6.0 Combined check on correlation-corrected dws for n=nmin..nmax
           This code applies to a single ksi-value */

   sum1 = sum2 =0.0;
   acfa[0] = acfv[0] = 1.0;
   
   if (nmax >= ntau) nmax = ntau-1; /* restricted sums may not be empty */
   
   for (njan=nmin; njan<=nmax; njan++) {
      dipn2av = di2av = d2av = dipndiav = expav = 0.0;
     
      for (ll = 1; ll <= ntau; ll++) {
         ftau[ll] = ppvalu(brek,coef,&l,&k,&tau[ll],0);
         /* All residuals are taken in normalized form : */
         q[ll] = (ftau[ll]-gtau[ll])*sqrt(weight[ll]);
         if (grfun == 1 && k != 1) {
            deriv[ll] = ppvalu(brek,coef,&l,&k,&tau[ll],1);
            if (fabs(deriv[ll]) > dermax) dermax = fabs(deriv[ll]);
         }
         err = q[ll];
         d2av += err*err;
         if (ll > njan) {
            errm = q[ll-njan];
            dipn2av  += err*err;
            di2av    += errm*errm;
            dipndiav += err*errm;
            expav    += acffunc(acfindex,tau[ll]-tau[ll-njan],ksi);
         }
      }
       
      errl2 = d2av;
      d2av /= (double)ntau;
      dipn2av /= (double)(ntau-njan);
      di2av /= (double)(ntau-njan);
      dipndiav /= (double)(ntau-njan);
      expav /= (double)(ntau-njan);
      sum1 += ((dipn2av + di2av)/d2av);
      dev = (dipndiav/d2av - expav);
      sum2 += dev;
      if (dev < small) small = dev;
      if (dev > big)   big   = dev;
      acfa[njan] = dipndiav/d2av;
      acfv[njan] = expav;
   }
     
   dws  = (sum1 - 2.0*sum2)/(double)(nmax-nmin+1);
   dws *= (1.0 - 1.0/(double)ntau);

   freed = ntau-km1-l;
   rms = sqrt((errl2/freed));
   acfdevmin = small;
   acfdevmax = big;
/* debug printf("ksi=%g l=%d freed=%d rms=%g nmin=%d nmax=%d dws=%g small=%g, big=%g\n",
ksi,l,freed,rms,nmin,nmax,dws,acfdevmin,acfdevmax); */
   if (prfun == 1) {
      for (ll = 1; ll <= ntau; ll++) {
         a = ntau*weight[ll]/totalw; /* why this variable? */
         dtau = ppvalu(brek,coef,&l,&k,&tau[ll],1);
         ddtau = ppvalu(brek,coef,&l,&k,&tau[ll],2);
   /*      if (fullon) {
    fprintf(fpp, " %g %g %g %g %g %g %g\n",
      tau[ll], gtau[ll], ftau[ll], q[ll], a, dtau, ddtau); */
   /* 22-6-1996: */
   /* fprintf(fpp, " %g %g %g %g %g %g %g\n", tau[ll], gtau[ll], ftau[ll], q[ll], weight[ll], dtau, ddtau); */
   /* } */
      }
   }
   return(0);
}

/*---------------------------------------------------------------------------*/
double ppvalu(double ara[],double dara[][NMAX+1],int *iptr,int *jptr,double *xptr,
      int jderiv)

   /*  evaluates the jderiv-th derivative of a pp-function  */
{
   int   ip, jp, fmmjdr, i, ndummy, m;
   double  xp, h;
   ip = *iptr;
   jp = *jptr;
   xp = *xptr;
   ppvalue = 0.;
   fmmjdr = jp-jderiv;
   if (fmmjdr > 0)
   {
      interv(ara,&ip,&xp,&i,&ndummy);
      h = xp-ara[i];
      for (m = jp; m >= jderiv+1; m--)
      {
         ppvalue = (ppvalue/fmmjdr)*h+dara[m][i];
         fmmjdr--;
      }
   }
   return(ppvalue);
}

/*---------------------------------------------------------------------------*/
int interv(double ara[],int *iptr,double *xptr,int *jptr,int *kptr)

   /*  locates a point within an increasing sequence of points  */
{
   int   ip, jp, kp, istep, middle, ilos=0;
   double   xp;
   ip = *iptr;
   xp = *xptr;
   kp = 10;
   ihi = ilo+1;
   if (ihi >= ip)
   {
      if (xp >= ara[ip]) kp = 1;
      else
      {
         if (ip <= 1) kp = -1;
         else
         {
            ilo = ip-1;
            ihi = ip;
         }
      }
   }
   if (kp == 10)
   {
      if (xp < ara[ihi])
      {
         if (xp >= ara[ilo]) kp = 0;
         else
         {
            istep = 1;
            while (ilo > 1 && xp < ara[ilo])
            {
               ihi = ilo;
               ilo = ihi-istep;
               istep *= 2;
            }
            if (ilo <= 1)
            {
               ilo = 1;
               if (xp < ara[1]) kp = -1;
            }
         }
      }
      else
      {
         istep = 1;
         while (ihi < ip && xp > ara[ihi])
         {
            ilo = ihi;
            ihi = ilo+istep;
            istep *= 2;
         }
         if (ihi >= ip)
         {
            ihi = ip;
            if (xp > ara[ip]) kp = 1;
         }
      }
      if (kp == 10)
      {
         do
         {
            middle = (ilo+ihi)/2;
            if (xp >= ara[middle])
            {
               ilos = ilo;
               ilo = middle;
            }
            else ihi = middle;
         }
         while (middle != ilos);
      }
   }
   if (kp == -1) jp = 1;
   else
   {
      if (kp == 1) jp = ip;
      else
      {
         kp = 0;
         jp = ilo;
      }
   }
   *jptr = jp;
   *kptr = kp;
   return(0);
}

/*---------------------------------------------------------------------------*/
int dwtest(double x,double *xptr,double *yptr)

   /* calculates percentages of lower (*xptr) and upper (*yptr)
    * percentages of x
    */
{
   int   i;
   double   y, pin, ci, sc, scc, aa, pl, pu, xp, yp;
   pin = 4*atan(1.)/ntau;
   sc = scc = 0;
   for (i = 1; i < ntau-freed; i++)
   {
      ci = cos(pin*i);
      sc += ci;
      scc += ci*ci;
   }
   sc = sc/freed;
   aa = (double) freed*(freed+2)*(1-sc*sc)/(ntau-2-2*scc-2*freed*sc*sc)-1;
   pl = aa*(1-sc)/2;
   pu = aa*(1+sc)/2;
   y = x/4;
   xp = betai(pl,pu,y);
   yp = betai(pu,pl,y);
   *xptr = xp;
   *yptr = yp;
   return(0);
}

/*---------------------------------------------------------------------------*/
double betai(double a,double b,double x)

   /* calculates Ix(a,b)   */
{
   double   bt, betai;
   if (x == 0 || x == 1) bt = 0;
   else bt = exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*log(x)+b*log(1-x));
   if (x < (a+1)/(a+b+2)) betai = bt*betacf(a,b,x)/a;
   else betai = 1-bt*betacf(b,a,1-x)/b;
   return(betai);
}

/*---------------------------------------------------------------------------*/
double betacf(double a,double b,double x)

   /* continued fraction for Ix(a,b) */
{
   int   em, tem;
   double   aold, am, bm, az, qab, qap, qam, bz, d, ap, bp, app, bpp;
   aold = 0;
   em = 0;
   am = bm = az = 1;
   qab = a+b;
   qap = a+1;
   qam = a-1;
   bz = 1-qab*x/qap;
   while (fabs(az-aold) > 1e-7*fabs(az))
   {
      em++;
      tem = em+em;
      d = em*(b-em)*x/(qam+tem)/(a+tem);
      ap = az+d*am;
      bp = bz+d*bm;
      d = -(a+em)*(qab+em)*x/(a+tem)/(qap+tem);
      app = ap+d*az;
      bpp = bp+d*bz;
      aold = az;
      am = ap/bpp;
      bm = bp/bpp;
      az = app/bpp;
      bz = 1;
   }
   return(az);
}

/*---------------------------------------------------------------------------*/
double chitest(double chi,int fr)

   /* calculates upper percentage points of chisquare of x with a degrees
      of freedom
    */
{
   double   a, x, gamser, gln, gammcf;
   a = (double) fr/2.;
   x = chi/2;
   if (x < a+1)
   {
      gser(&gamser,a,x,&gln);
      gammcf = 1-gamser;
   }
   else gcf(&gammcf,a,x,&gln);
   return(gammcf);
}

/*---------------------------------------------------------------------------*/
int gser(double *xptr, double a, double x,double *yptr)

   /* returns the incomplete gamma function p(a,x) evaluated by its
      series representation as *xptr; also returns ln(gamma(a)) as
      *yptr
    */
{
   double   xp, yp, ap, sum, del;
   yp = lgamma(a);
   ap = a;
   sum = 1/a;
   del = sum;
   while (fabs(del) > fabs(sum)*1e-7)
   {
      ap++;
      del *= x/ap;
      sum += del;
   }
   xp = sum*exp(-x+a*log(x)-yp);
   *xptr = xp;
   *yptr = yp;
   return(0);
}

/*---------------------------------------------------------------------------*/
int gcf(double *xptr,double a,double x,double *yptr)

   /* returns the incomplete gamma function q(a,x) evaluated by its
      continued fraction representation as *xptr; also returns ln(gamma(a))
      as *yptr
    */
{
   int   n;
   double   xp, yp, gold, b0, g, a0, b1, fac, a1, ana, anf;
   yp = lgamma(a);
   gold = b0 = 0.;
   a0 = b1 = fac = g = 1.;
   a1 = x;
   n = 1;
   while (fabs((g-gold)/g) > 1e-7)
   {
      if (a1 != 0 && n != 1) gold = g;
      ana = n-a;
      a0 = (a1+a0*ana)*fac;
      b0 = (b1+b0*ana)*fac;
      anf = n*fac;
      a1 = x*a0+anf*a1;
      b1 = x*b0+anf*b1;
      if (a1 != 0)
      {
         fac = 1/a1;
         g = b1*fac;
      }
      n++;
   }
   xp = exp(-x+a*log(x)-yp)*g;
   *xptr = xp;
   *yptr = yp;
   return(0);
}

/*---------------------------------------------------------------------------*/
double acffunc(int index, double deltax,double ksi)

   /* returns the value of the assumed autocorrelation function of the
      residuals in points spaced deltax apart. Normalized.
    */
{
   double a;
   if (ksi == 0.0) return(0.0);
   if (index==1) return(exp(-deltax/ksi));
   if (index==2) return(exp(-deltax*deltax/ksi/ksi/2.0));
   if (index==3) {
      a = 1.0-deltax/(2.0*ksi);
      return(a<0.0 ? 0.0 : a);
      }
   if (index==4) {
      a = 2.0*deltax/ksi;
      return(a==0.0 ? 1.0 : sin(a)/a);
      }
   return(0.0);
}

/*---------------------------------------------------------------------------*/
double lgamma(double xx)
{

  /* returns the log of the gamma function. Code supplied by Marcus Karolewski */
  
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
  24.01409824083091,-1.231739572450155, 0.1208650973866179e-2,0.5395239384953e-5};
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

/*------------------------------------------------------------------*/
double zbrent(double x1,double x2,int ind,double value)

   /* using Brent's method, find the root of a function known to lie
    * between x1 and x2; the root returned as x3, will be refined until
    * its accuracy is smaller than 1e-8 times f(x1) or f(x2), depending
    * which one is largest;
    */
{
   double   y1, y2, y, a, b, fa, fb, x3, fc, c=0, d=0, e=0, tol1, xm, s, p,
      q, r, eps;
   eps = 1e-8;
   a = x1;
   b = x2;
   x3 = 0;
   fa = ppvalu(brek,coef,&l,&k,&a,ind)-value;
   fb = ppvalu(brek,coef,&l,&k,&b,ind)-value;
   fc = fb;
   y1 = fabs(a);
   y2 = fabs(b);
   if (y1 > y2) y = y1;
   else y = y2;
   tol1 = 2*eps*y;
   while (x3 == 0)
   {
      if (fb*fc > 0)
      {
         c = a;
         fc = fa;
         d = b-a;
         e = d;
      }
      if (fabs(fc) < fabs(fb))
      {
         a = b;
         b = c;
         c = a;
         fa = fb;
         fb = fc;
         fc = fa;
      }
      xm = .5*(c-b);
      if (fabs(xm) <= tol1 || fb == 0) x3 = b;
      else
      {
         if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
         {
            s = fb/fa;
            if (a == c)
            {
               p = 2*xm*s;
               q = 1-s;
            }
            else
            {
               q = fa/fc;
               r = fb/fc;
               p = s*(2*xm*q*(q-r)-(b-a)*(r-1));
               q = (q-1)*(r-1)*(s-1);
            }
            if (p > 0) q *= -1;
            p = fabs(p);
            y1 = 3*xm*q-fabs(tol1*q);
            y2 = fabs(e*q);
            if (y1 < y2) y = y1;
            else y = y2;
            if (2*p < y)
            {
               e = d;
               d = p/q;
            }
            else
            {
               d = xm;
               e = d;
            }
         }
         else
         {
            d = xm;
            e = d;
         }
         a = b;
         fa = fb;
         if (fabs(d) > tol1) b += d;
         else
         {
            if (xm >= 0) y = tol1;
            else y = -1*tol1;
            b += y;
         }
         fb = ppvalu(brek,coef,&l,&k,&b,ind)-value;
      }
   }
   return(x3);
}
