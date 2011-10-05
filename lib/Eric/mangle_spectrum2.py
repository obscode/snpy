## Automatically adapted for numpy.oldnumeric Feb 04, 2009 by ipython

import numpy.oldnumeric as num
from filters import filters
from filters import filter
from filters import spectra
import Multipack
import deredden

filts = {}

#NAME
# colour_min2
#PURPOSE
# Thing to be minimized to get colours right.  Returns the fluxes of
# the spectrum scaled by scale_factors.  
#COMMON BLOCKS
# filter_info            Not modified.  Zero points must be set. 
#                         See calc_zps
# templates_colour_fit   Internal common block from mangle_spectrum.
#                         Not modified
#PROCEDURES CALLED
# apply_mangle
# filter_integ
#AUTHOR
# Mark Sullivan, August 2004
#-
def colour_min2(scale_factors, fjac=None, scale_factors_wave, tobe_mangled_wave,
      tobe_mangled_flux, method, rv, allfilters, z, wanted):

   # This differs from colour_min in 3 respects: First, it allows the
   # AB photometric system to be used.  Second, startnorm/finalnorm are not
   # normalized by any zeropoint in this version.  Finally, it allows
   # for an arbitrary redshift -- aconley 2005/06/16

   global filts

   #generate mangled spectrum
   spectrum_mangle_flux = apply_mangle(tobe_mangled_wave,tobe_mangled_flux,
                                    scale_factors_wave,scale_factors,
                                    method, Rv=rv)

   # calculate the new fluxes in each filter
   if(useccm):
      temp_flux = zeros((len(allfilters),), num.Float32)
      for i in range(len(allfilters):
         temp_flux[i]=filts[allfilters[i]].response(tobe_mangled_wave,
            spectrum_mangle_flux, z)
         temp_flux[i] /= num.power(10.0, 0.4*filts[allfilters[i]].zp)
   else:
      temp_flux = zeros((len(allfilters-2),), num.Float32)
      for i in range(len(allfilters)-2):
         temp_flux[i]=filts[allfilters[i+1]].response(tobe_mangled_wave,
            spectrum_mangle_flux, z)
         temp_flux[i] /= num.power(10.0, 0.4*filts[allfilters[i+1]].zp)
      
   return (0, temp_flux - wanted)
   
#NAME
# apply_mangle
#PURPOSE
# Apply the mangling to a spectrum given pretabulated coefficients,
#  which are determined by mangle_spectrum2
#INPUTS
# tobe_mangled_wave      Input wavelength vector
# tobe_mangled_flux      Input flux vector
#OPTIONAL INPUTS
# sw                     Wavelengthts of knots (not used in CCM)
# sf                     Scale factors of spline (E(B-V) for CCM,
#                        first element only)
# Rv                     Rv to use (CCM only; def: 3.1)
# ebmv                   alternative way to pass ebmv for CCM
#KEYWORDS
# spline              Use spline interpolation (DEFAULT)
# linear              Use linear interpolation (uses linterp)
# poly                Use polynomial interpolation (Not yet implemented)
# ccm                 Use the CCM dust law
#RETURNS
# Mangled spectrum
#COMMON BLOCKS
#  None
#AUTHOR
# Mark Sullivan
#-
def apply_mangle(tobe_mangled_wave,tobe_mangled_flux,sw,sf,method='spline',
      verbose=0,Rv=3.1):

   ccm = linear = spline = poly = 0
   if method=='spline':  spline = 1
   if method=='linear':  linear = 1
   if method=='poly':  poly = 1
   if method=='ccm':  ccm = 1
   
   if(verbose):
      if(spline): print 'Using spline interpolation'
      if(linear): print 'Using linear interpolation'
      if(poly): print 'Using polynomial interpolation'
      if(ccm): print 'Using CCM dust law'
   
   if(not ccm):
      nscale_factors=len(sf)
      nwavelengths=len(tobe_mangled_wave)
      
      # extrapolate the scalefactors in the ends of the spectrum using
      # linear extrapolation
      start_scale=sf[0]
      end_scale=sf[-1]
      
      manglefactor=num.zeros((nwavelengths,), num.Float32)
      ispline_start = num.nonzero(num.greater(tobe_mangled_wave-sw[0], 0))[0]
      ispline_end = num.nonzero(num.greater(tobe_mangled_wave-sw[-1], 0))[0]

      if ispline_start < 0:  ispline_start = 0
      
      if (ispline_start GT 0) : manglefactor[0:ispline_start-1]=start_scale
      if (ispline_end LT (nwavelengths-1)): 
         manglefactor[ispline_end+1:-1]=end_scale
      
      # do the mangling to the passed spectrum
      
      if(spline):
         tck = Multipack.splrep(sw, sf, k=3)
         manglefactor[ispline_start:ispline_end] =  Multipack.splev(tck, 
               tobe_mangled_wave[ispline_start:ispline_end])
      elif (linear):
         tck = Multipack.splrep(sw, sf, k=1)
         manglefactor[ispline_start:ispline_end] =  Multipack.splev(tck, 
               tobe_mangled_wave[ispline_start:ispline_end])
      elif (poly) :
         raise RuntimeError,"Polynomial mangling not supported"
     
      spectrum_mangle_flux=tobe_mangled_flux*manglefactor
   else:
      #;Using CCM law
      ebmv=sf[0]
      scalefac = sf[1]
      spectrum_mangle_flux,a,b = deredden.unred(tobe_mangled_wave,
            tobe_mangled_flux,-ebmv, R_V=rv 
      spectrum_mangle_flux *= scalefac
   
   return(spectrum_mangle_flux)

#NAME
# mangle_spectrum2
#PURPOSE
# Mangles an input spectrum so that it has the desired colors
# This version avoids some of the common blocks
#USAGE
# mang = mangle_spectrum( spectrum_wave, spectrum_flux, passedfilters,
#                         wanted_colours, vega_struct )
#INPUTS
# spectrum_wave       Wavelength vector for input spectrum
# spectrum_flux       Flux vector for input spectrum
# passedfilters       Filter numbers.
# wanted_colours      Desired colors.  Corresponding to [0]-[1],
#                     [1]-[2], [2]-[3] where [] refers to the entries
#                     of passedfilters
# vega_struct         Structure holding info about vega.  See read_vega
#OPTIONAL INPUTS
# redshift            What redshift to  this at.  Default: 0.0.
#                      Note that the returned SED is still in the rest
#                      frame, as the filters are shifted instead of
#                      the SED.
# normfilternum       Index into passedfilters for the filter used to
#                      keep the normalization of the spectrum.  After
#                      mangling, the flux through this filter should
#                      be the same as before mangling.  Default: 0 --
#                      that is, the first filter in passedfilter is
#                      used to normalize
# fixedfilters        Numbers of anchor filters.  if not provided, the
#                      code will choose its own anchor filters
# photsys             Photometric system colors are defined in.  1 for
#                      Vega, 2 for AB.  Default: 1
# truncate_filters     This adjusts the filter response to neglect all
#                     parts of the filter where the throughput is
#                     <truncate_pc%. Some filters have very wide definitions
#                     even where there is essentially no throughput. Not
#                     recommed for light-curve fitting, but great
#                     for mangling real-life spectra! Default: NO
# truncate_pc         Default 1: percentage throughput at whihc to
#                     truncate the filter responses
# anchorwidth         Width of the automatic anchorfilters. Default: 100A
# Rv                  Rv to use (CCM only; def: 3.1)
# noublesample      Array of bytes specifying whether to not use
#                      ublesampling in filter_integ.  The default is
#                      all 0b's -- i.e., to use ublesampling.
#OPTIONAL OUTPUTS
# factors_w           The wavelengths of the spline knots.  Note that
#                      these will be in the rest frame of the SED you
#                      pass in, so if you set the redshift parameter
#                      these may be divided by a factor of 1+z
#                      you n't expect. Not used for CCM.
# factors_s           The scale factors of the applied function. for
#                     CCM the first element is the best-fit E(B-V)
# ebmv                Alternative way of getting best E(B-V); CCM only
#KEYWORDS
# spline              Use spline interpolation (DEFAULT)
# linear              Use linear interpolation (uses linterp)
# poly                Use polynomial interpolation (Not yet implemented)
# ccm                 Use the CCM dust law
# quiet               Run MPFIT quietly. Default:Yes. You must specify
#                     QUIET=0 to override this. (extreme debug).
# verbose             Verbose output
# usemeanwave         Use the mean wavelengths of the filters for the
#                      knot points rather than the effective
#                      wavelength
# gradient            Use the gradient to set the anchor boundary
#                      conditions instead of fixed levels
#RETURNS
# The output spectrum, normalized so that the flux in the first filter
# is the same before and after mangling
#NOTES
# if factors_w and factors_s are known ahead of time, it is far faster
#  to use apply_mangle.
#COMMON BLOCKS
# filter_info         Note this is modified internally by adding the
#                      anchor filters to the .  These are removed
#                      at the  of this code
# templates_color_fit Created by this routine for passing to
#                      colour_min2, the minimizer function
#PROCEDURES CALLED
# filter_integ
# colour_min2
# effective_wavelength
# calc_ab_zp
# calc_vega_zp
# linterp
#AUTHOR
# Mark Sullivan, August 2004
#-

#;Little function to calculate colour of spectrum through two filters
#;Kind of like kcorr, but both filters are integrated at the object
#;redshift.  Returns the color in magnitudes

def mangle_spectrum2_getcolour(filt1,filt2,wave,flux,redshift):
   global filts

   f1 = filts[filt1]
   f2 = filts[filt2]

   flux1 = f1.response( wave, flux, redshift)
   flux2 = f2.response( wave, flux, redshift)
   
   return(-2.5*log10(flux1/flux2) - f1.zp + f2.zp)


def mangle_spectrum2(s_wave,s_flux,passedfilters, wanted_colours,
      fixedfilters=None, normfilter=None,z=0,verbose=0,
      factors_s=None, factors_w=None, gradient=0, usemeanwave=1, 
      truncate_filters=0, truncate_pc=1.0, quiet=1, anchorwidth=100,
      method='spline', Rv=3.1)

   global filts

   usespline=usepoly=uselinear=useccm = 0
   if method == 'spline':  usespline = 1
   if method == 'linear':  uselinear = 1
   if method == 'poly':  usepoly = 1
   if method == 'ccm':  useccm = 1
   if verbose:
      print ''
      print 'MANGLE_SPECTRUM2 starting.'
      if(usespline): print 'Using spline interpolation'
      if(uselinear): print 'Using linear interpolation'
      if(usepoly): print 'Using polynomial interpolation'
      if(useccm): print 'Using CCM dust law'
   
   # get into percentage
   truncate_pc=truncate_pc/100.
   
   obj_redshift = z
   opz = 1.0 + obj_redshift
   
   #;Set which filter will be used to normalize
   if normfilter is None:
      normfilter = passedfilters[0]
   
   # check spectrum is OK
   if len(s_wave) != len(s_flux):
      raise RuntimeError, 'spectrum wave and flux must have same number' + \
                            ' of elements in mangle_spectrum !'
   if num.maximum.reduce(s_flux) <= 0.0:
      raise RuntimeError, 'ERROR in mangle_spectrum: Your spectrum must ' + \
                          ' contain some positive (non-zero) fluxes.'
   
   if len(passedfilters) != len(wanted_colours)+1:
      raise RuntimeError, 'Number of colours passed must be equal to number '+ \
         'filters - 1 in mangle_spectrum !'
   
   # truncate_filters assumes filter response are normalised to a peak
   # throughput of 1 !
   if(truncate_filters): 
      if(verbose): print ''
      if(verbose): print 'Truncated passed filters:'
      for i in range(len(passedfilters):
   
                                   # find peak throughput
         filts[passedfilters[i]] = filters[passedfilters[i]].copy()
         f = filts[passedfilters[i]]
         f.resp = num.where(num.greater(f.resp/num.maximum.reduce(f.resp), 
            truncate_pc), f.resp, 0.0)
         ids = num.equal(f.resp, 0.0)
         f.resp = f.resp[ids[0]-1:ids[-1]+1]
         f.wave = f.wave[ids[0]-1:ids[-1]+1]
         f.min = num.minimum.reduce(f.wave)
         f.max = num.maximum.reduce(f.wave)
         
   # check that the (possibly truncated) passedfilters and spectrum completely 
   #  overlap
   low_waves = []
   high_waves = []
   for i in range(len(passedfilters)):
       f = filts[passedfilters[i]]
       low_wave.append(f.min)
       high_wave.append(f.max)
   
   lowest_filter_wavelength = min(low_wave)
   highest_filter_wavelength = max(high_wave)
   
   bluest_spectrum_wavelength=s_wave[0]
   reddest_spectrum_wavelength=s_wave[-1]
   
   if(verbose): 
      print ''
      print 'Filter definitions cover  ',lowest_filter_wavelength,'A to ', \
            highest_filter_wavelength,'A'
      print 'SED covers ',bluest_spectrum_wavelength,'A to ',\
            reddest_spectrum_wavelength,'A'
   
   if(lowest_filter_wavelength < bluest_spectrum_wavelength or 
      highest_filter_wavelength > reddest_spectrum_wavelength): 
      print 'Problem in mangle_spectrum: SED es not cover filter definitions'
      if(not verbose): 
         print 'Filter definitions cover  ',lowest_filter_wavelength,\
               'A to ',highest_filter_wavelength,'A'
         print 'SED covers ',bluest_spectrum_wavelength,'A to ', \
               reddest_spectrum_wavelength,'A'
   
   # copy the passed arrays to the common block
   #tobe_mangled_flux=spectrum_flux
   #tobe_mangled_wave=spectrum_wave
   
   nwavelengths=len(s_wave)
   
   # fixed filters contain extra spline knots - these are the "anchor" filters
   # if these are not specified, we must work out what they should be
   
   # make an array containg all the filters
   internalfilters=0
   allfilters=passedfilters[:]
   
   #;Choose the fixedfilters if not set by the caller
   if(not fixedfilters and  not useccm): 
       internalfilters=1
       if(verbose) : print \
         'No fixed filters passed - defining internally..'
       nfilters=len(passedfilters)
       mean_wave=num.zeros((nfilters,), num.Float64)
       for i in range(len(passedfilters)):
           pf = filts[passedfilters[i]]
           if ( usemeanwave ) :  
              mean_wave[i]=pf.ave_wave / opz
           else: 
              mean_wave[i]=pf.eff_wave(s_wave,s_flux,z) / opz
       mean_wave = num.array(mean_wave)
       #dummy=MIN(mean_wave,bluest_filter,SUBSCRIPT_MAX=reddest_filter)
       reddest_filter = num.argmax(mean_wave)
       bluest_filter = num.argmin(mean_wave)
   
       #extract_filter,passedfilters[bluest_filter],filter_wave,dummy
       bluest_filter_wave=filts[bluest_filter].wave[0]
       #extract_filter,passedfilters[reddest_filter],filter_wave,dummy
       reddest_filter_wave=filts[reddest_filter].wave[-1]
   
       # Create two new fake filters
       filts['blue'] = filter('blue_anchor')
       filts['red'] = filter('red_anchor')
   
       filts['blue'].wave = num.array([bluest_filter_wave-anchorwidth-2.,
             bluest_filter_wave-anchorwidth-1., bluest_filter_wave-anchorwidth, 
             bluest_filter_wave-3., bluest_filter_wave-2.,bluest_filter_wave-1.])
       filts['blue'].resp = num.array([0.,0.,1.,1.,0.,0.])
   
       #;red anchor filter
       filts['red'].wave = num.array([reddest_filter_wave+1.,
            reddest_filter_wave+2., reddest_filter_wave+3.,
            reddest_filter_wave+anchorwidth, reddest_filter_wave+anchorwidth+1.,
            reddest_filter_wave+anchorwidth+2.])
       filts['red'].response = num.array([0.,0.0,1.,1.,0.0,0.])
   
       filts['red'].zp = filts['red'].compute_zpt(spectra['vega_Fukugita'], 0.0)
       filts['blue'].zp = filts['blue'].compute_zpt(spectra['vega_Fukugita'], 0.0)
   
       filts['red'].ave_wave = num.sum(filts['red'].wave)/len(filts['red'].wave)
       filts['blue'].ave_wave = num.sum(filts['blue'].wave)/len(filts['blue'].wave)
   
       lowest_filter_wavelength=bluest_filter_wave-anchorwidth-2.
       highest_filter_wavelength=reddest_filter_wave+anchorwidth+2.
   
       if(lowest_filter_wavelength < s_wave[0] or \
          highest_filter_wavelength > s_wave[1]): 
          print 'Problem in mangle_spectrum: SED es not cover anchor filter '+\
                'definitions'
          if(not verbose): 
             print 'Anchor filter definitions cover  ',lowest_filter_wavelength,\
                   'A to ',highest_filter_wavelength,'A'
             print 'SED covers ',bluest_spectrum_wavelength,'A to ',\
                   reddest_spectrum_wavelength,'A'
             print "Try reducing ANCHORWIDTH. But you're probably screwed."
       if(verbose): print 'Anchor filter definitions cover  ',\
             lowest_filter_wavelength,'A to ',highest_filter_wavelength,'A'
   
   if(not useccm): 
       #; add the fixedfilters to the allfilters array if we're not using CCM
       allfilters=['blue'] + allfilters + ['red']
   
   nfilters=len(allfilters)
   
   # calculate effective wavelengths of spectrum for the passed filters
   # these are the positions of the spline knots
   mean_wave=num.zeros((nfilters), num.Float32)
   for i in range(nfilters)
       if ( usemeanwave ) : 
           mean_wave[i] = filts[allfilters[i]].ave_wave / opz
       else: 
           mean_wave[i] = filts[allfilters[i]].eff_wave( s_wave,s_flux,z ) / opz
   if(not useccm): 
      passed_mean_wave=mean_wave[1:-2]
   
      if (mean_wave[0] >= min(passed_mean_wave)) : 
         raise RuntimeError,"Your lower fixed filter must be bluer than the " +\
               "bluest filter passsed!"
      if (mean_wave[nfilters-1] <= max(passed_mean_wave)) : 
         raise RuntimeError,"Your upper fixed filter must be redder than the " + \ 
               "reddest filter passsed!"
   else: 
      passed_mean_wave=mean_wave
   
   
   # B-band flux of tobe_mangled_wave/tobe_mangled_flux
   # Note that we don't bother with the zeropoint here
   startnorm = filts[normfilter].response( s_wave, s_flux, z)
   if(verbose): 
      #; print wanted colours
      print 'The colours you want are:'
      for i in range(len(wanted_colours)):
         print "%s (%10.3f) minus %s (%10.3f) = %8.4f" % (passedfilters[i],\
               passed_mean_wave[i],passedfilters[i+1],\
               passed_mean_wave[i+1],wanted_colours[i]
   
   if(verbose): 
   # print colours of spectrum passed
      spectrum_colour=num.zeros((len(wanted_colours)), num.Float32)
      print 'The colours of the spectrum you passed are:'
      for i in range(len(wanted_colours)):
         #;Kcorr no longer works here because we are generalizing
         #; to other redshifts, and we have to integrate both through
         #; the same filter
         spectrum_colour[i] = \
            mangle_spectrum2_getcolour( passedfilters[i],passedfilters[i+1],\
                                        s_wave, s_flux, obj_redshift)
   
         print "%s (%10.3f) minus %s (%10.3f) = %8.4f" % (passedfilters[i],
               passed_mean_wave[i], passedfilters[i+1],
               passed_mean_wave[i+1],spectrum_colour[i])
   
   if(not useccm): 
      # setup initial scale factors array same size as number of passed filters
      # set to 1.0 at start of fit
      scale_factors_or=num.ones((nfilters,), num.Float32)
      nscale_factors=nfilters
   
      # setup an array of "wanted" fluxes
      wantedflux_at_mean_wave=num.zeros((nfilters-2),num.Float32)
   else: 
      scale_factors_or=num.array([ 0.0, 1.0 ]) #;ebmv, fluxscale
      nscale_factors=2
      wantedflux_at_mean_wave=num.zeros((nfilters),num.Float32)
   
   
   wantedflux_at_mean_wave[0]=filts[passedfilters[0]].response(s_wave,
                                            s_flux,z)
   wantedflux_at_mean_wave[0] /= num.power(10, 0.4*filts[passedfilters[0]].zp)
   
   # calculate relative flux in other filters to first passed filter
   # using passed colours
   for i in range(len(wanted_colours)):
      wantedflux_at_mean_wave[i+1]=  num.power(10.0,0.4*wanted_colours[i]) * \
                                     wantedflux_at_mean_wave[i]
   
   #;Get normfilter to have the same flux before and after
   desired_flux_in_normfilter = startnorm / \
         num.power(10.0, 0.4*filts[inormfilter].zp)
   
   wantedflux_at_mean_wave *= desired_flux_in_normfilter / \
     wantedflux_at_mean_wave[ passedfilters.index(normfilternum) ]
   
   # setup our fitting weights including the two  values used
   # for constraining the spline
   # sigma is not used but I think must be specified
   # sigma = num.ones((len(mean_wave)),num.Float32)*wantedflux_at_mean_wave[0]/1000.0
   # weights = num.ones((len(mean_wave)),num.Float32)*1000.0
   
   # sort the wavelengths of the passed filters
   # these are set as the wavelengths of the scale factors
   sort_index=num.argsort(mean_wave)
   scale_factors_wave=num.take(mean_wave, sort_index)
   #sigma = num.take(sigma, sort_index)
   #weights = num.take(weights, sort_index)
   
   if(useccm): 
      pi = [{'parname':'E(B-V)', 'fixed':0, 'limited':[0,0], 
                  'limits':[0.0,0.0]},
                 {'parname':'scale', 'fixed':0, 'limited':[0,0], 
                  'limits':[0.0,0.0]}]
      pi[0] = 0.0
      pi[1] = 1.0
   else: 
      #; Set the boundaries etc. of the fitted scales
      #; nothing is limited except fluxscale cannot be negative
      pi = [{'fixed':0, 'limited':[1,0], 'limits':[0.0,0.0]}]*nscale_factors
   
      if(gradient): 
         #; use the gradient of the two filters nearest the edge to set the
         #; fluxscale in the anchor filters
         #; gradient = dy/dx i.e. d(fluxscale)/d(wavelength)
         #; red anchor
         wavediff=( mean_wave[-1] - mean_wave[-2] ) / \
                  ( mean_wave[-2] - mean_wave[-3] )
         pi[nfilters-1]['tied']='P[%d] + %8.4f*(P[%d]-P[%d])' % (nfilters-2,
               wavediff, nfilters-2, nfilters-3)
         #blue anchor
         wavediff=(mean_wave[0]-mean_wave[1])/(mean_wave[1]-mean_wave[2])
         pi[0]['tied']='P[1]+%8.4f*(P[1]-P[2])' % (wavediff)
      else: 
         # tie the scales at the fixed filters to the scale in the filter next or
         pi[0]['tied']='P[1]'
         pi[nfilters-1]['tied']='P[%d]' % (nfilters-2)
   
      for i in range(nscale_factors):
         pi[i]['value'] = 1.0
   
   # FIT! Change XTOL to get better fits
   status=0
   if(verbose): print 'Calling MPFITFUN...'
   result = mpfit.mpfit(colour_min2, parinfo=pi, quite=quiet, maxiter=200,
         ftol=1.e-10, gtol=1.0e-10, xtol=1.0e-6, functkw=\
         {'scale_factors_wave':scale_factors_wave,
          'tobe_mangled_wave':s_wave,
          'tobe_mangled_flux':s_flux,
          'method':method, 'rv':Rv, 'allfilters':allfilters, 'z':z,
          'wanted':wantedflux_at_mean_wave})
   if (result.status == 5) : print \
     'Maximum number of iterations exceeded in mangle_spectrum'
   
   scale_factors = num.array(result.params)
   
   # generate final spectrum
   final_mangled_flux = apply_mangle(s_wave,s_flux,
                                       scale_factors_wave,scale_factors,
                                       method, Rv=Rv)
   
   # scale spectrum so it has the same B-flux as that passed
   # No zeropoint here because it wasn't used when startnorm was set
   final_norm=filts[normfilter].response(s_wave, final_mangled_flux,z)
   if final_norm > 0.d0: 
      factor=startnorm/final_norm 
   else: 
      factor=1.d0
   final_mangled_flux *= factor
   
   #  LEFT OFF HERE!!!
   if (verbose): 
      if(not useccm): 
          #;We have two extra filters -- the 'guard' ones
          temp_flux=DBLARR(nfilters-2)
          init_flux = DBLARR(nfilters-2)
      else: 
          temp_flux=DBLARR(nfilters)
          init_flux = DBLARR( nfilters )
      print "Filter:","InitFlux","WantedFlux:","FinalFlux:",\
            forMAT='(A-8,3X,A12,3X,A12,3X,A12)'
      for i=0,N_ELEMENTS(passedfilters)-1  
          cfiltnum = passedfilters[i]
          init_flux[i]=filter_integ( cfiltnum,tobe_mangled_wave,\
                                     tobe_mangled_flux,obj_redshift,\
                                     NOUBLESAMPLE=noublesample[i])
          temp_flux[i]=filter_integ( cfiltnum,tobe_mangled_wave,\
                                     final_mangled_flux,obj_redshift,\
                                     NOUBLESAMPLE=noublesample[i])
          CASE photo_system OF
              1 : 
                  temp_flux[i] /= master_filter[cfiltnum].vegazp
                  init_flux[i] /= master_filter[cfiltnum].vegazp
              
              2 : 
                  temp_flux[i] /= master_filter[cfiltnum].abzp
                  init_flux[i] /= master_filter[cfiltnum].abzp
              
          CASE
          print cfiltnum,init_flux[i],wantedflux_at_mean_wave[i],\
                temp_flux[i],forMAT='(I-8,3X,F12.5,3X,F12.5,3X,F12.5)'
      for
   if
   
   if(verbose): 
   # and check we did it right
      print 'The final colours of the mangled spectrum are:'
      for i=0,N_ELEMENTS(wanted_colours)-1  
         spectrum_colour[i] = \
            mangle_spectrum2_getcolour( passedfilters[i],passedfilters[i+1],\
                                        tobe_mangled_wave, final_mangled_flux,\
                                        obj_redshift, photo_system,\
                                        [noublesample[i],noublesample[i+1]])
   
         if(not useccm): 
            print master_filter[passedfilters[i]].shortid,\
                  mean_wave[i+1],master_filter[passedfilters[i+1]].shortid,\
                  mean_wave[i+2],spectrum_colour[i],\
                  forMAT='(%"%s (%10.3f) minus %s (%10.3f) = %8.4f")'
         else: 
            print master_filter[passedfilters[i]].shortid,\
                  mean_wave[i],master_filter[passedfilters[i+1]].shortid,\
                  mean_wave[i],spectrum_colour[i],\
                  forMAT='(%"%s (%10.3f) minus %s (%10.3f) = %8.4f")'
         ELSE
      for
   if
   
   if(useccm): ebmv=scale_factors[0]
   
   return(final_mangled_flux, scale_factors_wave, scale_factors)
   
