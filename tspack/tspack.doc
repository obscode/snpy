                TSPACK:  Tension Spline Curve Fitting Package

                               Robert J. Renka
                                  05/27/91


        I.  INTRODUCTION


             The primary purpose of TSPACK is to construct a smooth
        function which interpolates a discrete set of data points.
        The function may be required to have either one or two con-
        tinuous derivatives, and, in the C-2 case, several options
        are provided for selecting end conditions.  If the accuracy
        of the data does not warrant interpolation, a smoothing func-
        tion (which does not pass through the data points) may be
        constructed instead.  The fitting method is designed to avoid
        extraneous inflection points (associated with rapidly varying
        data values) and preserve local shape properties of the data
        (monotonicity and convexity), or to satisfy the more general
        constraints of bounds on function values or first derivatives.
        The package also provides a parametric representation for con-
        structing general planar curves and space curves.

             The fitting function h(x) (or each component h(t) in the
        case of a parametric curve) is defined locally, on each
        interval associated with a pair of adjacent abscissae (knots),
        by its values and first derivatives at the endpoints of the
        interval, along with a nonnegative tension factor SIGMA
        associated with the interval (h is a Hermite interpolatory
        tension spline).  With SIGMA = 0, h is the cubic function
        defined by the endpoint values and derivatives, and, as SIGMA
        increases, h approaches the linear interpolant of the endpoint
        values.  Since the linear interpolant preserves positivity,
        monotonicity, and convexity of the data, h can be forced to
        preserve these properties by choosing SIGMA sufficiently
        large.  Also, since SIGMA varies with intervals, no more
        tension than necessary is used in each interval, resulting in
        a better fit and greater efficiency than is achieved with a
        single constant tension factor.



        II.  USAGE


             TSPACK must be linked to a driver program which re-
        serves storage, reads a data set, and calls the appropriate
        subprograms selected from those described below in section
        III.B.  Header comments in the software modules provide
        details regarding the specification of input parameters and
        the work space requirements.  It is recommended that curves
        be plotted in order to assess their appropriateness for the
        application.  This requires a user-supplied graphics package.



        III.  SOFTWARE


        A)  Code

             The code is written in 1977 ANSI Standard Fortran.  All
        variable and array names conform to the default typing con-
        vention:  I-N for type INTEGER and A-H or O-Z for type REAL.
        (There are no DOUBLE PRECISION variables.)  There are 32
        modules, and they are ordered alphabetically in the package.
        Each consists of the following sections:

            1)  the module name and parameter list with spaces sepa-
                rating the parameters into one to three subsets:
                input parameters, I/O parameters, and output parame-
                ters (in that order);
            2)  type statements in which all parameters are typed
                and arrays are dimensioned;
            3)  a heading with the name of the package, identifica-
                tion of the author, and date of the most recent
                modification to the module;
            4)  a description of the module's purpose and other rel-
                evant information for the user;
            5)  input parameter descriptions and output parameter
                descriptions in the same order as the parameter
                list;
            6)  a list of other modules required (called either
                directly or indirectly);
            7)  a list of intrinsic functions called, if any; and
            8)  the code, including comments.

             Note that it is assumed that floating point underflow
        results in assignment of the value zero.  If not the default,
        this may be specified as either a compiler option or an
        operating system option.  Also, overflow is avoided by re-
        stricting arguments to the exponential function EXP to have
        value at most SBIG=85.  SBIG, which appears in DATA statements
        in the evaluation functions, HVAL, HPVAL, HPPVAL, and TSINTL,
        must be decreased if it is necessary to accomodate a floating
        point number system with fewer than 8 bits in the exponent.
        No other system dependencies are present in the code.

             The modules that solve nonlinear equations, SIGS, SIGBP,
        SIG0, SIG1, SIG2, and SMCRV, include diagnostic print capabi-
        lity which allows the iteration to be traced.  This can be
        enabled by altering logical unit number LUN in a DATA state-
        ment in the relevant module.


        B)  Module Descriptions

             The software modules are divided into three categories,
        referred to as level 1, level 2, and level 3, corresponding
        to the hierarchy of calling sequences:  level 1 modules call
        level 2 modules which call level 3 modules.  For most ap-
        plications, the driver need only call two level 1 modules --
        one from each of groups (a) and (b).  However, additional
        control over various options can be obtained by directly
        calling level 2 modules.  Also, additional fitting methods,
        such as parametric smoothing, can be obtained by calling
        level 2 modules.  Note that, in the case of smoothing or C-2
        interpolation with automatically selected tension, the use
        of level 2 modules requires that an iteration be placed around
        the computation of knot derivatives and tension factors.



        1) Level 1 modules

             These are divided into two groups.

        a)  The following modules return knots (in the parametric
            case), knot derivatives, tension factors, and, in the
            case of smoothing, knot function values, which define
            the fitting function (or functions in the parametric
            case).  The naming convention should be evident from the
            descriptions.


        TSPSI   Subroutine which constructs a shape-preserving or
                  unconstrained interpolatory function.  Refer to
                  TSVAL1.

        TSPSS   Subroutine which constructs a shape-preserving or
                  unconstrained smoothing spline.  Refer to TSVAL1.

        TSPSP   Subroutine which constructs a shape-preserving or
                  unconstrained planar curve or space curve.  Refer
                  to TSVAL2 and TSVAL3.

        TSPBI   Subroutine which constructs a bounds-constrained
                  interpolatory function.  The constraints are defined
                  by a user-supplied array containing upper and lower
                  bounds on function values and first derivatives,
                  along with required signs for the second derivative,
                  for each interval.  Refer to TSVAL1.

        TSPBP   Subroutine which constructs a bounds-constrained
                  parametric planar curve.  The constraints are de-
                  fined by user-supplied arrays containing upper and
                  lower bounds on the signed perpendicular distance
                  between the smooth curve segment and the polygonal
                  line segment associated with each knot interval.
                  The bounds might, for example, be chosen to avoid
                  intersections between smooth contour curves.  Refer
                  to TSVAL2.


        b)  The following modules return values, derivatives, or
            integrals of the fitting function(s).


        TSVAL1  Subroutine which evaluates a Hermite interpolatory
                  tension spline or its first or second derivative at
                  a user-specified set of points.  Note that a smooth-
                  ing curve constructed by TSPSS is the interpolant of
                  the computed knot function values.  The evaluation
                  points need not lie in the interval defined by the
                  knots, but care must be exercised in assessing the
                  accuracy of extrapolation.

        TSVAL2  Subroutine which returns values or derivatives of a
                  pair of Hermite interpolatory tension splines which
                  form the components of a parametric planar curve.
                  The output values may be used to construct unit tan-
                  gent vectors, curvature vectors, etc.

        TSVAL3  Subroutine which returns values or derivatives of
                  three Hermite interpolatory tension splines which
                  form the components of a parametric space curve.

        TSINTL  Function which returns the integral over a specified
                  domain of a Hermite interpolatory tension spline.
                  This provides an effective means of quadrature for
                  a function defined only by a discrete set of values.



        2) Level 2 modules

             These are divided into four groups.

        a)  The following modules are called by TSPSP and TSPBP to
            obtain a sequence of knots (parameter values) associated
            with a parametric curve.  For some data sets, it might
            be advantageous to replace these with routines that
            implement an alternative method of parameterization.


        ARCL2D  Subroutine which computes the sequence of cumulative
                  arc lengths associated with a sequence of points in
                  the plane.

        ARCL3D  Subroutine which computes the sequence of cumulative
                  arc lengths associated with a sequence of points in
                  3-space.


        b)  The following modules are called by the level 1, group (a)
            modules to obtain knot derivatives (and values in the case
            of SMCRV).


	YPC1	Subroutine which employs a monotonicity-constrained
		  quadratic interpolation method to compute locally
                  defined derivative estimates, resulting in a C-1
                  fit.

        YPC1P   Subroutine similar to YPC1 for the case of periodic
                  end conditions.  In the case of a parametric curve
                  fit, periodic end conditions are necessary to ob-
                  tain a closed curve.

        YPC2    Subroutine which determines a set of knot-derivative
		  estimates which result in a tension spline with two
		  continuous derivatives and satisfying specified
		  end conditions.

	YPC2P	Subroutine similar to YPC2 for the case of periodic
		  end conditions.

	SMCRV	Subroutine which, given a sequence of abscissae with
		  associated data values and tension factors, returns
		  a set of function values and first derivative values
		  defining a twice-continuously differentiable tension
                  spline which smoothes the data and satisfies either
                  natural or periodic end conditions.


        c)  The following modules are called by the level 1, group (a)
            modules to obtain tension factors associated with knot
            intervals.


        SIGS    Subroutine which, given a sequence of abscissae,
		  function values, and first derivative values,
		  determines the set of minimum tension factors for
		  which the Hermite interpolatory tension spline
		  preserves local shape properties (monotonicity
                  and convexity) of the data.  SIGS is called by
                  TSPSI, TSPSS, and TSPSP.

        SIGBI   Subroutine which, given a sequence of abscissae,
		  function values, and first derivative values,
		  determines the set of minimum tension factors for
		  which the Hermite interpolatory tension spline
                  satisfies specified bounds constraints.  SIGBI is
                  called by TSPBI.

        SIGBP   Subroutine which, given a sequence of points in the
                  plane with associated derivative vectors, determines
                  the set of minimum tension factors for which the
                  parametric planar tension spline curve defined by
                  the data satisfies specified bounds on the signed
                  orthogonal distance between the parametric curve and
                  the polygonal curve defined by the points.  SIGBP is
                  called by TSPBP.


        d)  The following functions are called by the level 1, group
            (b) modules to obtain values and derivatives.  These pro-
            vide a more convenient alternative to the level 1 routines
            when a single value is needed.


	HVAL	Function which evaluates a Hermite interpolatory ten-
		  sion spline at a specified point.

	HPVAL	Function which evaluates the first derivative of a
		  Hermite interpolatory tension spline at a specified
		  point.

	HPPVAL	Function which evaluates the second derivative of a
		  Hermite interpolatory tension spline at a specified
		  point.



        3) Level 3 modules

             These are divided into three groups.

        a)  The following functions are called by SIGBI to compute
            tension factors, and are convenient for obtaining an
            optimal tension factor associated with a single interval.


	SIG0	Function which, given a pair of abscissae, along with
                  associated endpoint values and derivatives, deter-
                  mines the smallest tension factor for which the
                  corresponding Hermite interpolatory tension spline
                  satisfies a specified bound on function values in
                  the interval.

	SIG1	Function which, given a pair of abscissae, along with
		  associated endpoint data, determines the smallest
		  tension factor for which the corresponding Hermite
		  interpolatory tension spline satisfies a specified
		  bound on first derivative values in the interval.

	SIG2	Function which, given a pair of abscissae, along with
		  associated endpoint data, determines the smallest
		  tension factor for which the corresponding Hermite
		  interpolatory tension spline preserves convexity
		  (or concavity) of the data.


        b)  The following modules are of general utility.


	INTRVL	Function which, given an increasing sequence of ab-
		  scissae, returns the index of an interval containing
                  a specified point.  INTRVL is called by the evalua-
                  tion functions TSINTL, HVAL, HPVAL, and HPPVAL.

	SNHCSH	Subroutine called by several modules to compute
		  accurate approximations to the modified hyperbolic
		  functions which form a basis for exponential ten-
		  sion splines.

        STORE   Function used by SIGBP, SIGS, SIG0, SIG1, and SIG2 in
                  computing the machine precision.  STORE forces a
                  value to be stored in main memory so that the pre-
                  cision of floating point numbers in memory locations
                  rather than registers is computed.


        c)  The remaining modules are listed below.


        B2TRI   Subroutine called by SMCRV to solve the symmetric
		  positive-definite block tridiagonal linear system
		  associated with the gradient of the quadratic
		  functional whose minimum corresponds to a smooth-
		  ing curve with nonperiodic end conditions.

	B2TRIP	Subroutine similar to B2TRI for periodic end
		  conditions.

	ENDSLP	Function which returns the derivative at X1 of a
		  tension spline h(x) which interpolates three
		  specified data points and has third derivative
		  equal to zero at X1.	ENDSLP is called by YPC2
		  when this choice of end conditions is selected
		  by an input parameter.

	YPCOEF	Subroutine called by SMCRV, YPC2, and YPC2P to com-
		  pute coefficients defining the linear system.



        IV.  REFERENCE


        For the theoretical background, consult the following:

          RENKA, R. J.  Interpolatory tension splines with automatic
          selection of tension factors. SIAM J. Sci. Stat. Comput. 8
          (1987), pp. 393-415.
