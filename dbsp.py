"""
dbsp.py is a pyraf-based reduction pipeline for spectra taken with the
Palomar 200-inch Double Spectrograph.
"""

# Authors:  Eric Bellm <ebellm@caltech.edu>,
#           Branimir Sesar <bsesar@astro.caltech.edu>
# License: BSD Style.

from pyraf import iraf
import numpy as np
import os
import shutil
import inspect
import pyfits
from scipy.optimize.minpack import leastsq
import copy
from glob import glob
import cosmics

# directory where the reduction code is stored
BASE_DIR = os.path.dirname(os.path.abspath(inspect.getfile(
                inspect.currentframe())))

def is_new_red_camera():
    """Utility for determining red camera version from image size."""
    ids = range(10)
    for id in ids:
        name = 'red{:04d}.fits'.format(id)
        if os.path.exists(name):
            hdr = pyfits.getheader(name)
            if hdr['NAXIS1'] == 4141 or hdr['NAXIS1'] == 4114:
                return True
            elif hdr['NAXIS1'] == 1024 or hdr['NAXIS1'] == 1124:
                return False
            else:
                raise ValueError('Unexpected image size')
        else:
            continue

    #raise ValueError('Could not locate red side files')
    print 'Could not locate red side files--defaulting to new camera'
    return True
    
NEW_RED_SIDE = is_new_red_camera()

# load IRAF packages
iraf.noao(_doprint=0)
iraf.imred(_doprint=0)
iraf.ccdred(_doprint=0)
iraf.twodspec(_doprint=0)
iraf.longslit(_doprint=0)
iraf.kpnoslit(_doprint=0)
iraf.astutil(_doprint=0)
iraf.onedspec(_doprint=0)

# defaults
# (blue trace is usually around 250-260)
det_pars = {'blue':{'gain':0.8,'readnoise':2.7,'trace':253,
                    'crval':4345, 'cdelt':-1.072, 'arc':'FeAr_0.5.fits',
                    'fwhm_arc':2.8}}

if NEW_RED_SIDE:
    det_pars['red'] = {'gain':2.8,'readnoise':8.5,'trace':166,
                       'crval':7502, 'cdelt':1.530, 'arc':'HeNeAr_0.5.fits',
                       'biassec':'[4128:4141,1:440]', 'fwhm_arc':2.4}
else:
    # old CCD
    det_pars['red'] = {'gain':2.05,'readnoise':7.8,'trace':130,
                       'crval':6600, 'cdelt':2.46, 'arc':'HeNeAr_0.5.fits',
                       'biassec':'[1080:1124,1:300]', 'fwhm_arc':1.6}
                    # crval is in Angstrom, cdelt is Angstrom/pixel

def mark_bad(side,numbers):
    """Utility for excluding specific files from further analysis.

    Saturated or mis-configured exposures are suffixed .bad so file searches
    do not find them.
    
    Parameters
    ----------
    side : {'blue' (default), 'red'}
        'blue' or 'red' to indicate the arm of the spectrograph
    numbers : list of int or int
        image id(s) to be marked as bad.

    """

    assert (side in ['blue','red'])
    try:
        for num in numbers:
            name = '{:s}{:04d}.fits'.format(side,num)
            os.rename(name, name+'.bad')
    except TypeError:
        # single number
        num = numbers
        name = '{:s}{:04d}.fits'.format(side,num)
        os.rename(name, name+'.bad')

def create_arc_dome(side='blue',trace=None,arcslit=0.5,overwrite=True):
    """Convenience function which subtracts bias, creates dome flats, and
    creates arc frames.
    
    Parameters
    ----------
    side : {'blue' (default), 'red'}
        'blue' or 'red' to indicate the arm of the spectrograph
    trace : int
        row or column of the spectral trace, if different from default.

    """

    assert ((side == 'blue') or (side == 'red'))
    if trace is None:
        trace = det_pars[side]['trace']
    
    bias_subtract(side=side,trace=trace)

    if side == 'blue':
        fix_bad_column_blue()

    make_flats(side=side,overwrite=overwrite)

    if side == 'blue':
        make_arcs_blue(slit=arcslit,overwrite=overwrite)
    else:
        make_arcs_red(slit=arcslit,overwrite=overwrite)

def bias_subtract(side='blue',trace=None):
    """Use iraf.ccdproc to subtract bias from all frames.

    Parameters
    ----------
    side : {'blue' (default), 'red'}
        'blue' or 'red' to indicate the arm of the spectrograph
    trace : int
        row or column of the spectral trace, if different from default.
    """

    # update the headers
    iraf.asthedit('%s????.fits' % side, '/home/bsesar/opt/python/DBSP.hdr')
    if side == 'blue':
        iraf.hedit('blue*.fits', 'DISPAXIS', 2, update="yes", 
            verify="no", add="yes", show="no")
    else:
        iraf.hedit('red*.fits', 'DISPAXIS', 1, update="yes", 
            verify="no", add="yes", show="no")

    # bias subtraction using the overscan
    filenames = glob("%s????.fits" % side)
    hdr = pyfits.getheader(filenames[0])
    iraf.unlearn('ccdproc')
    iraf.ccdproc.zerocor = "no"
    iraf.ccdproc.flatcor = "no"
    iraf.ccdproc.fixpix = "no"
#    iraf.ccdproc.fixfile = "../bluebpm"
    if side == 'blue':
        iraf.ccdproc.biassec = hdr['BSEC1']
        iraf.ccdproc.trimsec = "[%d:%d,*]" % (trace-100, trace+100)
        iraf.ccdproc.function = "spline3"
        iraf.ccdproc.order = 3
    else:
        # trim the specified region
        iraf.ccdproc.biassec = det_pars['red']['biassec'] # this may not work for the old camera...
        tsec_x = hdr['TSEC1'].split(',')[0]
        iraf.ccdproc.trimsec = tsec_x + ",%d:%d]" % (trace-100, trace+100)
        iraf.ccdproc.function = "legendre"
        iraf.ccdproc.order = 1
    iraf.ccdproc.darkcor = "no"
    iraf.ccdproc.ccdtype = ""
    iraf.ccdproc.niterate = 3
    iraf.ccdproc('%s????.fits' % side)

def fix_bad_column_blue():
    """Uses a science exposure to find the bad column in the blue CCD. 
    Automatically corrects all blue exposures with iraf.fixpix.
    """

    science = iraf.hselect('blue????.fits', '$I', 'TURRET == "APERTURE" & LAMPS == "0000000"', Stdout=1)
    if len(science) > 0:
        f = open('science_dump', 'w')
        for fn in science:
            f.write('%s\n' % fn)
        f.close()
    science = iraf.hselect('@science_dump', '$I', 'TURRET == "APERTURE" & LAMPS == "0000000" & AIRMASS > 1.01', Stdout=1)
    os.unlink('science_dump')
    f = pyfits.open(science[0])
    bad_column = f[0].data[1608,:].argmin() + 1
    f.close()
    f = open('bluebpm', 'w')
    f.write('%d %d 1 2835\n' % (bad_column, bad_column))
    f.close()
    iraf.fixpix('blue????.fits', "bluebpm")

def find_flats(aperture, side='blue'):
    """Finds flat images for the selected aperture and side using 
    header keywords.  
    
    Uses dome flats if present, otherwise falls back on internal flats.

    Parameters
    ----------
    aperture : {'0.5','1.0', '1.5', '2.0'}
        string indicating the slit width
    side : {'blue' (default), 'red'}
        'blue' or 'red' to indicate the arm of the spectrograph
    """

    # find dome flat images
    domeflats = iraf.hselect('%s????.fits' % side, '$I', 'TURRET == "APERTURE" & APERTURE == "%s" & LAMPS == "0000000"' % aperture, Stdout=1)
    # safest to check that IMGTYPE = flat 
    # otherwise zenith pointings can get into the flats
    domeflats = iraf.hselect(','.join(domeflats), '$I', 'TURRET == "APERTURE" & APERTURE == "%s" & LAMPS == "0000000" & AIRMASS < 1.01 & IMGTYPE == "flat"' % aperture, Stdout=1)
    # find internal flat (incandescent lamp) images
    intflats = iraf.hselect('%s????.fits' % side, '$I', 'TURRET == "LAMPS" & APERTURE == "%s" & LAMPS == "0000001"' % aperture, Stdout=1)
    intflats = iraf.hselect(','.join(intflats), '$I', 'TURRET == "LAMPS" & APERTURE == "%s" & LAMPS == "0000001" & AIRMASS < 1.01' % aperture, Stdout=1)
    # dome flats are prefered over internal flats
    flats = []
    if (len(intflats) > 0) & (len(domeflats) == 0):
        flats = intflats
        print "Using %d internal flats for the %s arcsec slit." % (len(intflats), aperture)
    if len(domeflats) > 3:
        flats = domeflats
        print "Using %d dome flats for the %s arcsec slit." % (len(domeflats), aperture)

    return flats

def make_flats(side='blue',overwrite=False):
    """Creates dome flat images using iraf.flatcombine.

    Creates flats for all slit widths in ['0.5','1.0', '1.5', '2.0'] if present.
    Uses dome flats if present, otherwise falls back on internal flats.
    Stores flats as flat_{side}_{aper}.fits.

    Parameters
    ----------
    side : {'blue' (default), 'red'}
        'blue' or 'red' to indicate the arm of the spectrograph
    overwrite : boolean
        If True, overwrites existing flats; if False, skips processing for
        all slit widths with existing flats.
    """

    iraf.unlearn('flatcombine')
    iraf.flatcombine.ccdtype = ""
    iraf.flatcombine.process = "no"
    iraf.flatcombine.subsets = "no"
    iraf.flatcombine.rdnoise = "RON"
    iraf.flatcombine.gain = "GAIN"
    for aperture in ['0.5','1.0', '1.5', '2.0']:
        flats = find_flats(aperture, side=side)
        if len(flats) > 0:
            if overwrite:
                iraf.delete('flat_%s_%s.fits' % (side, aperture), verify='no')
                if len(flats) < 3:
                    iraf.flatcombine(','.join(flats), output='temp', reject='average')
                if len(flats) >= 3:
                    iraf.flatcombine(','.join(flats), output='temp', reject='avsigclip')

                # normalize the flat
                iraf.unlearn('response')
                iraf.response.function = "spline3"
                iraf.response.order = 100
                iraf.response.niterate = 3
                iraf.response.low_rej = 3
                iraf.response.high_rej = 3
                if side == 'blue':
                    iraf.twodspec.longslit.dispaxis=2
                else:
                    iraf.twodspec.longslit.dispaxis=1
                iraf.response('temp', 'temp', 
                    'flat_%s_%s.fits' % (side, aperture), interactive="no")
                os.rename('temp.fits', 'raw_flat_%s_%s.fits' % (side, aperture))

                # measure flat-field error from sigma images
                iraf.unlearn('imcombine')
                iraf.imcombine.reject = 'avsigclip'
                iraf.imcombine(','.join(flats), output='flat', sigma='sigma', scale='mode')
                iraf.imarith('sigma', '/', 'flat', 'frac')
                s = iraf.imstat('frac.fits', fields="mean", nclip=20, Stdout=1, format="no")
                print 'Flat field error: ', np.float(s[0])
                iraf.delete('flat.fits', verify="no")
                iraf.delete('sigma.fits', verify="no")
                iraf.delete('frac.fits', verify="no")
        else:
            print "No dome or internal flats for the %s arcsec slit." % aperture


def make_arcs_blue(slit=0.5, overwrite=False):
    """Creates the master FeAr arc with iraf.imcombine.

    Stores arc as FeAr_{slit}.fits.

    Parameters
    ----------
    slit : {'0.5','1.0', '1.5', '2.0'}
        string indicating the slit width
    overwrite : boolean
        If True, overwrites existing arc; if False, skips processing.
    """

    aperture = "{:3.1f}".format(slit)

    iraf.unlearn('imcombine')
    iraf.imcombine.rdnoise = det_pars['blue']['readnoise']
    iraf.imcombine.gain = det_pars['blue']['gain']
    arcs = iraf.hselect('blue????.fits', '$I', 'TURRET == "LAMPS" & APERTURE == "{aperture}" & LAMPS == "0100000"'.format(aperture=aperture), Stdout=1)
    try:
        arcs = iraf.hselect(','.join(arcs), '$I', 'TURRET == "LAMPS" & APERTURE == "{aperture}" & LAMPS == "0100000" & AIRMASS < 1.01'.format(aperture=aperture), Stdout=1)
    except:
        pass
    if overwrite:
        iraf.delete('FeAr_{aperture}.fits'.format(aperture=aperture), 
            verify='no')
    iraf.imcombine(','.join(arcs), 'FeAr_{}'.format(aperture), reject="none")

def make_arcs_red(slit=0.5, overwrite=False):
    """Creates the master HeNeAr arc with iraf.imcombine.

    Stores arc as HeNeAr_{slit}.fits.

    Parameters
    ----------
    slit : {'0.5','1.0', '1.5', '2.0'}
        string indicating the slit width
    overwrite : boolean
        If True, overwrites existing arc; if False, skips processing.
    """

    aperture = "{:3.1f}".format(slit)

    iraf.unlearn('imcombine')
    iraf.imcombine.rdnoise = det_pars['red']['readnoise']
    iraf.imcombine.gain = det_pars['red']['gain']
    arcs = iraf.hselect('red????.fits', '$I', 'TURRET == "LAMPS" & APERTURE == "{aperture}" & LAMPS == "0001110"'.format(aperture=aperture), Stdout=1)
    try:
        arcs = iraf.hselect(','.join(arcs), '$I', 'TURRET == "LAMPS" & APERTURE == "{aperture}" & LAMPS == "0001110" & AIRMASS < 1.01'.format(aperture=aperture), Stdout=1)
    except:
        pass
    if overwrite:
        iraf.delete('HeNeAr_{aperture}.fits'.format(aperture=aperture), 
            verify='no')
    iraf.imcombine(','.join(arcs), 'HeNeAr_{}'.format(aperture), reject="none")

def preprocess_image(filename, side='blue', flatcor = 'yes', 
    remove_cosmics=True, trace=None):
    """Remove instrumental signatures from a CCD image.
    
    Performs bias subtraction and flat correction, 
    adds header info if needed, and removes cosmic rays.

    Parameters
    ----------
    filename : string
        Name of CCD image to process
    side : {'blue' (default), 'red'}
        'blue' or 'red' to indicate the arm of the spectrograph
    flatcor : {'yes' (default), 'no'}
        Divide the image by the flat field?
    remove_cosmics : boolean
        Apply cosmic ray rejection?
    trace : int
        Row or column of the spectral trace, if different from default.
    """

    assert(flatcor in ['yes','no'])

    if trace is None:
        trace = det_pars[side]['trace']

    # bias subtraction using the overscan
    hdr = pyfits.getheader(filename)
    iraf.unlearn('ccdproc')
    iraf.ccdproc.zerocor = "no"
    iraf.ccdproc.flatcor = flatcor
    iraf.ccdproc.fixpix = "no"
    iraf.hedit(filename, 'GAIN', det_pars[side]['gain'], 
        update="yes", verify="no", show="no")
    iraf.hedit(filename, 'RON', det_pars[side]['readnoise'], 
        update="yes", verify="no", show="no")
    if side == 'blue':
        # update the header
        iraf.hedit(filename, 'DISPAXIS', 2, update="yes", verify="no", add="yes", show="no")
        # trim the specified region
        iraf.ccdproc.biassec = hdr['BSEC1']
        iraf.ccdproc.trimsec = "[%d:%d,*]" % (trace-100, trace+100)
        iraf.ccdproc.function = "spline3"
        iraf.ccdproc.order = 3
    else:
        # update the header
        iraf.hedit(filename, 'DISPAXIS', 1, update="yes", verify="no", add="yes", show='no')
        # trim the specified region
        iraf.ccdproc.biassec = det_pars['red']['biassec']
        tsec_x = hdr['TSEC1'].split(',')[0]
        iraf.ccdproc.trimsec = tsec_x + ",%d:%d]" % (trace-100, trace+100)
        iraf.ccdproc.function = "legendre"
        iraf.ccdproc.order = 1
    iraf.ccdproc.ccdtype = ""
    iraf.ccdproc.darkcor = "no"
    iraf.ccdproc.niterate = 3
    iraf.ccdproc(filename,
        flat="flat_%s_%s" % (side,hdr['APERTURE']))

    if (side == 'blue') and ('FIXPIX' not in hdr):
        iraf.fixpix('blue????.fits', "bluebpm")

    if 'OBSERVAT' not in hdr:
        # update the headers
        iraf.asthedit(filename, '/home/bsesar/opt/python/DBSP.hdr')

    # remove cosmic rays with LA Cosmic
    if remove_cosmics and ('COSMIC' not in hdr) and (hdr['EXPTIME'] > 60) and \
            (hdr['TURRET'] == 'APERTURE'):
        array, header = pyfits.getdata(filename, header=True)
        c = cosmics.cosmicsimage(array, gain=det_pars[side]['gain'], 
        readnoise=det_pars[side]['readnoise'], 
        sigclip = 4.5, sigfrac = 0.5, objlim = 2.0, satlevel=60000)
        c.run(maxiter = 3)
        header.update('COSMIC', 1, '1 if we ran LA Cosmic')
        pyfits.writeto(filename, c.cleanarray, header, clobber=True)

def store_standards(imgID_list, side='blue', trace=None, 
    arc=None, splot='no', redo='no', resize='yes', 
    crval=None, cdelt=None, extract=True, telluric_cal_id=None):
    """Extract spectra for spectroscopic standard stars and determine 
    corrections needed for fluxing.
    
    Wraps iraf.standard and iraf.sensfunc.  
    
    After extracting spectra for each standard star, the code will query you 
    for the name of the standards and prompt you to edit the bandpasses.  
    The calibrated data are stored in std-{side}.

    Next, the software fits a sensitivity function to the calibrated intensity
    bandpasses, prompting you to edit the fit as needed.  The fit sensitivity
    is stored in sens-{side}.

    See the README for step-by-step instructions for interacting with iraf.

    Parameters
    ----------
    imgID_list : list of ints
        List of file numbers for images of spectroscopic standards,
        e.g., red0011.fits and red0015.fits -> [11,15]
    side : {'blue' (default), 'red'}
        'blue' or 'red' to indicate the arm of the spectrograph
    trace : int
        Row or column of the spectral trace, if different from default.
    arc :

    splot : {'yes', 'no' (default)}
        Plot extracted spectra with iraf.splot?
    redo : {'yes', 'no' (default)}
        Redo the spectral extraction from scratch?  Passed to iraf.doslit.
        **Warning**--discards previous wavelength solution!  
        Use extract=False if you simply want to redefine calibration bandpasses
        and refit the sensitivity function.
    resize : {'yes' (default), 'no'}
        Resize the extraction aperture?  Passed to iraf.doslit.    
    crval :

    cdelt :

    extract : boolean (default True)
        Extract spectra?  If False, use existing extractions.
        Useful for redefining the calibration bandpasses and refitting
        the sensitivity function.
    telluric_cal_id : int or None (default)
        If defined, use specified image id to perform telluric correction.
    """

    # first extract all the standard spectra
    if extract:
        for i, imgID in enumerate(imgID_list):
            
            # only redo the first one so we don't have to keep re-defining the
            # dispersion solution
            if i == 0:
                redo = redo
            else:
                redo = 'no'
            
            extract1D(imgID, side=side, trace=trace, arc=arc, splot=splot,
                redo=redo, resize=resize, flux=False, crval=crval, cdelt=cdelt,
                telluric_cal_id = telluric_cal_id)

    iraf.delete('std-{}'.format(side),verify='no')
    iraf.delete('sens-{}'.format(side),verify='no')

    iraf.unlearn('standard')
    iraf.standard.caldir = "onedstds$iidscal/"
    iraf.standard.output = 'std-{}'.format(side)
    iraf.standard.bandwidth = 50
    iraf.standard.bandsep = 100
    # try these one at a time
    for imgID in imgID_list:
        # use the extracted spectrum!
        iraf.standard('%s%04d.spec.fits' % (side,imgID))

    iraf.unlearn('sensfunc')
    iraf.sensfunc.standards = 'std-{}'.format(side)
    iraf.sensfunc.sensitivity = 'sens-{}'.format(side)
    if side == 'blue':
        iraf.sensfunc.order = 3
    else:
        iraf.sensfunc.order = 6
    # varun says to use ignoreaps, but it's causing me problems downstream
    #iraf.sensfunc.ignoreaps = 'yes'
    iraf.sensfunc()

def estimateFWHM(imgID, side='blue'):
    """
    Use IRAF's imexam to measure the FWHM of the trace.
    """
    iraf.unlearn('imexam')
    iraf.rimexam.fittype = "gaussian"
    iraf.delete('trace.xy', verify="no")
    iraf.delete('fwhm.log', verify="no")
    # extract the position of the trace
    f = open('database/ap%s%04d' % (side, imgID), 'r')
    dat = f.read()
    xy = dat.split('\n')[5].split()[1:3]
    f.close()
    f = open('trace.xy', 'w')
    f.write('%s %s\n' % (xy[0], xy[1]))
    f.close()
    # run imexam
    if side == 'blue':
        defkey='j'
    else:
        defkey='k'
    iraf.imexam('%s%04d' % (side, imgID), '1', logfile='fwhm.log', keeplog="yes", defkey=defkey, imagecur='trace.xy', use_display="no", autoredraw="no")
    # load values
    f = open('fwhm.log', 'r')
    dat = f.read()
    fwhm = float(dat.split('\n')[1].split('=')[4].split()[0])
    f.close()
    # cleanup
    os.unlink("fwhm.log")
    os.unlink("trace.xy")

    # update the header
    f = pyfits.open('%s%04d.spec.fits' % (side, imgID))
    f[0].header.update('FWHM', np.round(fwhm, 2), 'FWHM estimate of the trace [pix]')
    f.writeto('%s%04d.spec.fits' % (side, imgID), clobber=True)
    f.close()
    if os.access('%s%04d_flux.spec.fits' % (side, imgID), os.F_OK):
        f = pyfits.open('%s%04d_flux.spec.fits' % (side, imgID))
        f[0].header.update('FWHM', np.round(fwhm, 2), 'FWHM estimate of the trace [pix]')
        f.writeto('%s%04d_flux.spec.fits' % (side, imgID), clobber=True)
        f.close()

def extract1D(imgID, side='blue', trace=None, arc=None, splot='no', redo='no', 
        resize='yes', flux=False, telluric_cal_id=None, reextract=False, 
        crval=None, cdelt=None):

    assert (side in ['blue','red'])
    assert (splot in ['yes','no'])
    assert (redo in ['yes','no'])
    assert (resize in ['yes','no'])

    rootname = '%s%04d' % (side,imgID)

    if trace is None:
        trace = det_pars[side]['trace']

    if arc is None:
        arc = det_pars[side]['arc']

    if crval is None:
        crval = det_pars[side]['crval']

    if cdelt is None:
        cdelt = det_pars[side]['cdelt']

    if reextract:
        apfile = 'database/ap'+rootname
        if os.path.exists(apfile):
            os.remove(apfile)
            

    # preprocess the science image
    preprocess_image(rootname+'.fits', side=side, trace=trace)

    #preprocess the arc image
    preprocess_image(arc, side=side, trace=trace, flatcor='no', 
        remove_cosmics=False)

    # set up doslit
    fwhm = 5.6
    iraf.unlearn('doslit')
    iraf.unlearn('apslitproc')
    iraf.unlearn('aidpars')
    iraf.doslit.readnoise = "RON"
    iraf.doslit.gain = "GAIN"
    iraf.doslit.width = 1.4*fwhm
    iraf.doslit.crval = crval
    iraf.doslit.cdelt = cdelt
    iraf.doslit.nsum = 50
    iraf.doslit.clean = "yes"
    if side == 'blue':
        iraf.doslit.dispaxis = 2
    else:
        iraf.doslit.dispaxis = 1
    iraf.doslit.extras = "yes"
    iraf.doslit.lower = -1*int(round(iraf.doslit.width/2.))
    iraf.doslit.upper = int(round(iraf.doslit.width/2.))
    iraf.doslit.t_function = "legendre"
    iraf.doslit.t_niter = 3
    iraf.doslit.t_order = 4
    iraf.doslit.t_high = 2
    iraf.doslit.t_low = 2
    iraf.doslit.weights = "variance"
    iraf.doslit.b_order = 3
    iraf.doslit.b_niterate = 1
    iraf.doslit.select = "average"
    anullus_start = fwhm*2
    xL = np.floor(np.linspace(-80,-1*anullus_start,10))
    xR = np.floor(np.linspace(anullus_start,80,10))
    background_range = ''
    for i in np.arange(xL.size-1):
        background_range += '%d:%d,' % (np.int(xL[i]), np.int(xL[i+1]-1))
    for i in np.arange(xR.size-1):
        background_range += '%d:%d,' % (np.int(xR[i]+1), np.int(xR[i+1]))
    iraf.doslit.b_sample = background_range[:-1]
    iraf.doslit.i_function = "legendre"
    iraf.doslit.i_order = 4
    if side == 'blue':
        iraf.doslit.coordlist = BASE_DIR + '/cal/FeAr_dbsp.dat'
        fwhm_arc = det_pars['blue']['fwhm_arc'] # input FWHM of arc lines here (in pixels)
    else:
        iraf.doslit.coordlist = BASE_DIR + '/cal/HeNeAr_dbsp.dat'
        fwhm_arc = det_pars['red']['fwhm_arc'] # input FWHM of arc lines here (in pixels)
    iraf.doslit.fwidth = fwhm_arc
    iraf.doslit.match = 10. # positive number is angstrom, negative is pix
    iraf.doslit.i_niterate = 5
    iraf.doslit.addfeatures = 'no'
    iraf.doslit.linearize = "yes"

    # extract 1D spectrum
    #print arc, iraf.doslit.crval, iraf.doslit.cdelt
    #iraf.epar('doslit')
    iraf.doslit(rootname+'.fits', arcs=arc, splot=splot, redo=redo, resize=resize)

    # extract the trace from the flat-field image
    hdr = pyfits.getheader(rootname + '.fits')
    iraf.apall('raw_flat_%s_%s.fits' % (side, hdr['APERTURE']), interactive='no', extras='no', edit="no", nsum=10, recenter="no", trace="no", background="none", output='trace_flat', find="no", reference=rootname)
    # normalize the response with mean
    m = np.float(iraf.imstat('trace_flat.fits', fields='mean', Stdout=1, format="no")[0])
    iraf.imarith('trace_flat.fits', '/', m, 'trace_flat.fits')
    # transform from pixel to wavelength coordinates
    iraf.hedit('trace_flat.fits', 'REFSPEC1', '%s%s.ms' % (rootname, os.path.splitext(arc)[0]), add="yes", verify="no", show="no")
    iraf.dispcor('trace_flat.fits', 'd_trace_flat.fits', table=rootname + '.ms.fits')

    # correct tellurics, if requested
    if telluric_cal_id is not None and side == 'red':
        tell_rootname = '%s%04d' % (side,telluric_cal_id)
        if not os.path.exists('norm_' + tell_rootname + '.fits'):
            normalize_to_continuum(telluric_cal_id,side=side)
        iraf.unlearn('telluric')
        iraf.telluric.input = rootname + '.ms.fits'
        iraf.telluric.output = ""
        iraf.telluric.sample = "7584:7678,9252:9842"
        iraf.telluric.interactive = "no"
        iraf.telluric.cal = 'norm_%s.fits' % tell_rootname
        iraf.telluric.ignoreaps = 'yes'
        iraf.telluric.xcorr = 'yes'
        iraf.telluric.tweakrms = 'yes'
        iraf.telluric()

    # measure shift with sky lines *before* fluxing to avoid floating point errors
    # measure the position and width of sky lines (do this only for exposures longer than 3 min)
    hdr = pyfits.getheader(rootname + '.fits')
    midpoint_loc = {'blue':4750,'red':7400}
    if hdr['EXPTIME'] > 120:
        iraf.unlearn('scopy')
        # background band
        iraf.delete(rootname + '.2001.fits', verify='no')
        iraf.scopy(rootname + '.ms.fits', rootname, band=3, format="onedspec")
        iraf.unlearn('fitprofs')
        iraf.fitprofs.gfwhm = fwhm_arc
        iraf.fitprofs.nerrsample = 100
        iraf.fitprofs.sigma0 = det_pars[side]['readnoise']
        iraf.fitprofs.invgain = 1./det_pars[side]['gain']
        iraf.delete('skyfit*.dat', verify='no')

        sky_lines = {'blue':
            {'wavelength':[4046.565, 4358.335, 5460.750, 5577.340],
            'regs':['4040 4054', '4350 4365', '5455 5469', '5572 5583']},
            'red':
            {'wavelength':[6923.21,7340.885,7821.51,8430.174,8885.83],
            'regs':['6917 6930', '7333 7348', '7813 7830','8422 8439','8875 8895']}}
        
        offsets = []
        for i in range(len(sky_lines[side]['wavelength'])):
            iraf.fitprofs(rootname + '.2001.fits',
                reg=sky_lines[side]['regs'][i], 
                logfile='skyfit_{:s}_{:1d}.dat'.format(side,i), 
                pos=BASE_DIR + '/cal/skyline_{:s}_{:1d}.dat'.format(side,i), 
                verbose='no')

        # dump useful data from skyfit?.dat (center width err_center err_width)
            os.system('fgrep -v "#" skyfit_{side:s}_{num:1d}.dat |perl -pe "s/0\.\n/0\./g;s/^ +//;s/\(/ /g;s/\)/ /g;s/ +/ /g;" |cut -d" " -f1,6,8,13 > wavelength_offset_{side:s}_{num:1d}.dat'.format(side=side,num=i))
            try:
                dat = np.genfromtxt('wavelength_offset_{side:s}_{num:1d}.dat'.format(side=side,num=i) , usecols=(0,2), names="center, error")
            except:
                # keep going if there's a bad profile fit
                print "Warning: bad fit for wavelength_offset_{side:s}_{num:1d}.dat".format(side=side,num=i)
                continue
            assert (dat['center'].size == 1)
            offsets.append(dat['center'] - sky_lines[side]['wavelength'][i])
        print offsets

        offset_final = np.mean(offsets)
        error_at_mid = np.std(offsets, ddof=1)/np.sqrt(len(offsets))/ \
            midpoint_loc[side]*299792.458 # uncertainty in km/s at 4750 A

    # add wavelength shifts/ uncertainty in km/s to headers
    # (CRVAL1 doesn't seem to apply the shift correctly?)
    f = pyfits.open(rootname + '.ms.fits')
    hdr = f[0].header
    if hdr['EXPTIME'] > 120:
        f[0].header.update('WOFF', '%.2f' % offset_final, 'Wavelength offset from sky lines in A at {} A'.format(midpoint_loc[side]))
        f[0].header.update('VERR', '%.2f' % error_at_mid, 'Uncertainty in km/s at {} A'.format(midpoint_loc[side]))
    else:
        f[0].header.update('WOFF', '-99.99', 'Wavelength offset from sky lines in A at {} A'.format(midpoint_loc[side]))
        f[0].header.update('VERR', '-99.99', 'Uncertainty in km/s at {} A'.format(midpoint_loc[side]))
    f.writeto(rootname + '.ms.fits', clobber=True)
    f.close()

    # extract counts and uncertainties
    iraf.unlearn('scopy')
    iraf.delete(rootname + '.0001.fits', verify='no')
    iraf.delete(rootname + '.2001.fits', verify='no')
    iraf.delete(rootname + '.3001.fits', verify='no')
    iraf.scopy(rootname + '.ms.fits', rootname, band=1, format="onedspec")
    iraf.scopy(rootname + '.ms.fits', rootname, band=4, format="onedspec")

    # correct wavelength calibration using sky lines
    if hdr['EXPTIME'] > 120:
        iraf.specshift(rootname + '.0001', '%.3f' % (-offset_final))
        iraf.specshift(rootname + '.3001', '%.3f' % (-offset_final))

    # normalize the spectrum and the uncertainty using the reponse function
    iraf.imarith(rootname + '.0001', '/', 'd_trace_flat.fits', rootname + '.0001')
    iraf.imarith(rootname + '.3001', '/', 'd_trace_flat.fits', rootname + '.3001')
    iraf.delete('*trace_flat.fits', verify="no")

    # flux, if requested
    if flux:
        # TODO: run setairmass; ensure extinction is set up correctly
        iraf.unlearn('calibrate')
        iraf.calibrate.input = rootname+'.0001,'+rootname+'.3001'
        iraf.calibrate.output = rootname+'_flux.0001,'+rootname+'_flux.3001'
        # I'm not sure yet why this gets moved to .0001...
        iraf.calibrate.sensitivity = sensitivity='sens-{}.0001'.format(side)
        iraf.calibrate.ignoreaps = 'yes'
        iraf.calibrate()

    # output to text files
    def text_output(rootname,hdr_arc,flux=False):
        if flux:
            suffix = '_flux'
        else:
            suffix = ''
        iraf.delete(rootname + '%s.spec.txt' % (suffix), verify='no')
        iraf.delete(rootname + '%s.err.txt' % (suffix), verify='no')
        iraf.delete(rootname + '%s.spec.fits' % (suffix), verify='no')
        iraf.delete(rootname + '%s.err.fits' % (suffix), verify='no')
        iraf.unlearn('dispcor')
        iraf.dispcor(rootname + '%s.0001' % (suffix), 
                rootname + '%s.spec' % (suffix), w1=hdr_arc['CRVAL1'], 
                dw=hdr_arc['CDELT1'], nw=hdr_arc['NAXIS1'])
        iraf.wspectext(rootname + '%s.spec.fits' % (suffix), 
                rootname + '%s.spec.txt' % (suffix), header="no")
        iraf.dispcor(rootname + '%s.3001' % (suffix), 
                rootname + '%s.err' % (suffix), w1=hdr_arc['CRVAL1'], 
                dw=hdr_arc['CDELT1'], nw=hdr_arc['NAXIS1'], blank=1.0)
        iraf.wspectext(rootname + '%s.err.fits' % (suffix), 
                rootname + '%s.err.txt' % (suffix), header="no")

    hdr_arc = pyfits.getheader('%s.ms.fits' % os.path.splitext(arc)[0])
    text_output(rootname,hdr_arc,flux=False)
    if flux:
        text_output(rootname,hdr_arc,flux=True)

    # calculate SNR
    iraf.delete(rootname + '.snr.fits', verify='no')
    iraf.imarith(rootname + '.spec.fits', '/', rootname + '.err.fits', 
        rootname + '.snr.fits')

    # cleanup
    iraf.delete('skyfit*.dat', verify='no')
    iraf.delete('wavelength_offset*.dat', verify='no')
    iraf.delete(rootname + '.ms.fits', verify="no")
    iraf.delete(rootname + '.0001.fits', verify="no")
    iraf.delete(rootname + '.3001.fits', verify="no")
    if flux:
        iraf.delete(rootname + '_flux.0001.fits', verify="no")
        iraf.delete(rootname + '_flux.3001.fits', verify="no")


    # statistics
    hdr = pyfits.getheader(rootname + '.spec.fits')
    if hdr['EXPTIME'] > 120:
        print "Wavelengths are offset by %.3f A, zero-point uncertainty is %.2f km/s at %.0f A." % (offset_final, error_at_mid,midpoint_loc[side])
    snr_loc = {'blue':4000,'red':7000}
    wave1 = np.int(np.floor((snr_loc[side]-10 - hdr_arc['CRVAL1'])/hdr_arc['CDELT1']))
    wave2 = np.int(np.floor((snr_loc[side]+10 - hdr_arc['CRVAL1'])/hdr_arc['CDELT1']))
    try:
        s = iraf.imstat(rootname + '.snr.fits[%d:%d]' % (wave1, wave2), 
                fields='mean', nclip=20, Stdout=1, format="no")
        print "SNR = %.1f at %d A" % (np.float(s[0]),snr_loc[side])
    except iraf.IrafError:
        print "Warning: could not imstat SNR"

    # measure FWHM of the trace
    estimateFWHM(imgID, side=side)

def combine_sides(imgID_list_blue, imgID_list_red, output=None, splot='yes'):
    """imgID_lists are lists of numbers of extracted spectra:
    eg, [41] for red0041_flux.spec.fits"""

            
    blue_files = ['blue{:04d}_flux.spec.fits'.format(imgID) for imgID in imgID_list_blue]
    red_files = ['red{:04d}_flux.spec.fits'.format(imgID) for imgID in imgID_list_red]
    blue_err_files = ['blue{:04d}_flux.err.fits'.format(imgID) for imgID in imgID_list_blue]
    red_err_files = ['red{:04d}_flux.err.fits'.format(imgID) for imgID in imgID_list_red]


    if output is None:
        hdr = pyfits.getheader(blue_files[0])
        obj = hdr['OBJECT'].replace(' ','_')
        # create a unique name based on the input

        def ids_to_string(idlist):
            if len(idlist) == 1:
                return "{:d}".format(idlist[0])
            else:
                return "-".join(["{:d}".format(id) for id in idlist])

        output = obj + '_' + \
            ids_to_string(imgID_list_blue) + '+' + \
            ids_to_string(imgID_list_red) 

    # clobber the old output files if they exist
    iraf.delete(output+'.*.fits',verify='no')
    iraf.delete(output+'.*.txt',verify='no')
    
    # determine dispersion: downsample to lower-resolution spectrum
    hdr = pyfits.getheader(blue_files[0])
    dw_blue = hdr['CDELT1']
    hdr = pyfits.getheader(red_files[0])
    dw_red = hdr['CDELT1']
    dw = np.max([dw_blue,dw_red])

    # find wavelength ranges
    def wavelength_range(fits_list):
        mins = []
        maxes = []
        for fname in fits_list:
            spec = np.genfromtxt(fname.replace('fits','txt'), names='wave, flux', 
                    dtype='f4, f4')
            mins.append(spec['wave'].min())
            maxes.append(spec['wave'].max())
        return [np.array(mins).min(),np.array(maxes).max()]

    blue_range = wavelength_range(blue_files)
    red_range = wavelength_range(red_files)

    # find overlap region
    if red_range[0] >= blue_range[1]:
        raise ValueError('No overlap in wavelength solution between sides!')

    # specify total spectral range
    w1 = blue_range[0]
    w2 = red_range[1]

    # re-disperse to common wavelength solution
    def redisperse_list(files,dw,w1,w2,key='spec'):
        input_list = ','.join(files)
        disp_files = [f.replace(key,key+'-disp') for f in files]
        output_disp_list = ','.join(disp_files)
        iraf.unlearn('dispcor')
        iraf.dispcor.input = input_list
        iraf.dispcor.output = output_disp_list
        # keep existing wavelength endpoints
        iraf.dispcor.dw = dw
        iraf.dispcor.w1 = w1
        iraf.dispcor.w2 = w2
        iraf.dispcor.flux = 'no'
        iraf.dispcor()
        # write text files
        for output in disp_files:
            iraf.wspectext(output, output.replace('fits','txt'), header="no")

        return disp_files

    # delete any lingering files
    iraf.delete('*-disp.fits',verify='no')
    iraf.delete('*-disp.txt',verify='no')

    blue_files_redisp = redisperse_list(blue_files,dw,w1,w2)
    red_files_redisp = redisperse_list(red_files,dw,w1,w2)
    blue_err_files_redisp = redisperse_list(blue_err_files,dw,w1,w2,key='err')
    red_err_files_redisp = redisperse_list(red_err_files,dw,w1,w2,key='err')

    # combine individual sides
    coadd_spectra(blue_files_redisp,'tmp-blue')
    coadd_spectra(red_files_redisp,'tmp-red')

    # find optimum weighting between sides

    # combine sides, weighted by uncertainties
    coadd_spectra(['tmp-blue.spec.fits','tmp-red.spec.fits'],output,
        one_side=False)


    # clean up
    iraf.delete('*-disp.fits',verify='no')
    iraf.delete('*-disp.txt',verify='no')
    iraf.delete('tmp-*.fits',verify='no')
    iraf.delete('tmp-*.txt',verify='no')
    
    if splot == 'yes':
        iraf.splot(output+'.spec')

def match_spectra_leastsq(y, yref, yerr, yreferr):
    errfunc = lambda p, y, yref, yerr, yreferr: \
        (yref - p[0]*y)/np.sqrt(yerr**2. + yreferr**2.)
    p0 = [1.]
    # cast arguments to double to avoid fjac being zero:
    # http://stackoverflow.com/questions/12473406/scipy-optimize-leastsq-returns-best-guess-parameters-not-new-best-fit
    p, cov_x, infodict, mesg, ier = leastsq(errfunc,
            p0, args=(y.astype(np.float64), yref.astype(np.float64), 
            yerr.astype(np.float64), yreferr.astype(np.float64)),
            full_output = True)
    if ier >= 1:
        print 'Best-fit ratio: ', p[0]
        return p[0]
    else:
        raise ValueError('Matching did not converge: {}'.format(mesg))

def coadd_spectra(spec_list_fits, out_name, 
    use_ratios=False, ratio_range=[4200,4300], scale_spectra=True,
    one_side=True):
    """Scales input 1D spectra onto the same scale multiplicatively
       and then combines spectra using a weighted mean.
    """

    spec_list_txt = [f.replace('fits','txt') for f in spec_list_fits]

    # first spectrum in the list is always the reference spectrum
    hdr = pyfits.getheader(spec_list_fits[0])
    mjd = hdr['MJD']
    date_obs = hdr['DATE-OBS']
    epoch = hdr['EPOCH']
    observat = hdr['OBSERVAT']
    exptime = hdr['EXPTIME']
    seeing = hdr['FWHM']
    # save some keywords
    keys = ['OBJECT', 'OBSERVER', 'DICHROIC', 'APERTURE', 'LAMPS', 'UTSHUT', 'OBSLST', 'RA', 'DEC', 'HOURANG', 'HA', 'TELFOCUS', 'CASSPA', 'PARALLAC', 'CCDTEMP', 'ANGLE', 'GRATING', 'AIRMASS']
    mjd_blue = hdr['MJD']
    exptime_blue = hdr['EXPTIME']
    hdr_save = {}
    for key in keys:
        hdr_save[key] = hdr[key]
    verr = np.float(hdr['VERR'])**2
    spec_ref = np.genfromtxt(spec_list_txt[0], names='wave, flux', 
            dtype='f4, f4')
    err_ref = np.genfromtxt(spec_list_txt[0].replace('spec','err'), 
            names='wave, flux', dtype='f4, f4')
    wave = spec_ref['wave']
    spec_ref = spec_ref['flux'].view(np.ma.masked_array)
    err_ref = err_ref['flux'].view(np.ma.masked_array)


    #err_ref['flux'] = np.where(err_ref['flux'] <= 0, 1, err_ref['flux']) # reset bad error values to 1
    # boolean array: mask out invalid regions so average excludes zeros
    bad_err = err_ref <= 0
    spec_ref[bad_err] = np.ma.masked
    err_ref[bad_err] = np.ma.masked


    # spectra and their errors will be stored here
    spectra = np.ma.zeros((spec_ref.size, len(spec_list_fits)), dtype='f4')
    spectra_err = np.ma.zeros((spec_ref.size, len(spec_list_fits)), dtype='f4')

    spectra[:, 0] = spec_ref
    spectra_err[:, 0] = err_ref

    ratio = [1]

    for i, fname in enumerate(spec_list_fits[1:]):
        fname_txt = spec_list_txt[i+1]
        hdr = pyfits.getheader(fname)
        exptime += hdr['EXPTIME']
        seeing += hdr['FWHM']
        verr += np.float(hdr['VERR'])**2
        spec = np.genfromtxt(fname_txt, names='wave, flux', dtype='f4, f4')
        err = np.genfromtxt(fname_txt.replace('spec','err'), 
                names='wave, flux', dtype='f4, f4')
        spec = spec['flux'].view(np.ma.masked_array)
        err = err['flux'].view(np.ma.masked_array)
        # reset bad error values to 1
        #err['flux'] = np.where(err['flux'] <= 0, 1, err['flux']) 
        bad_err = err <= 0
        spec[bad_err] = np.ma.masked
        err[bad_err] = np.ma.masked

        spectra[:, i+1] = spec
        spectra_err[:, i+1] = err
        if scale_spectra:
            if use_ratios:
                # use the specified region to determine te ratio of spectra
                good = np.where((spec > ratio_range[0]) & 
                        (spec < ratio_range[1]))
                ratio.append(np.median(spec_ref[good]/spec[good]))
            else:
                spec_good_err = err > 0
                # identify overlap between sides
                wgd = (err_ref > 0) & (err > 0)

                ratio.append(match_spectra_leastsq(spec[wgd], 
                        spec_ref[wgd], err[wgd], 
                        err_ref[wgd]))

            

    spec_avg, sum_weights = np.average(spectra*ratio, weights=1./(spectra_err*ratio)**2, axis=1, returned=True)
    spec_err = 1./np.sqrt(sum_weights)
    # output coadded spectra and uncertainties
    f = open('%s.spec.txt' % out_name, 'w')
    g = open('%s.err.txt' % out_name, 'w')
    h = open('%s.snr.txt' % out_name, 'w')
    # add some header keywords
    for key in hdr_save.keys():
        f.write('# %s = %s\n' % (key, hdr_save[key]))
    if one_side:
        # exposure time and velocity error are only well-defined for
        # data combined from a single side
        f.write('# FWHM = %.2f\n' % float(seeing/len(spec_list_fits)))
        f.write('# VERR = %.2f\n' % np.sqrt(verr))
        f.write('# MJD = %.6f\n' % (mjd + exptime/(2.*60.*60.*24.)))
    else:
        # when combining sides, use the MJD and EXPTIME from the combined blue side
        f.write('# EXPTIME = %.0f\n' % exptime_blue)
        f.write('# MJD = %.6f\n' % mjd_blue)

    for x, y, z in zip(wave, spec_avg, spec_err):
        f.write('%.3f %.5g\n' % (x, y))
        g.write('%.3f %.5g\n' % (x, z))
        h.write('%.3f %.5g\n' % (x, y/z))
    f.close()
    g.close()
    h.close()
    # save as 1D IRAF FITS files
    iraf.delete('%s.spec.fits' % out_name, verify="no")
    iraf.delete('%s.err.fits' % out_name, verify="no")
    iraf.delete('%s.snr.fits' % out_name, verify="no")
    iraf.rspectext('%s.spec.txt' % out_name, '%s.spec.fits' % out_name, 
            crval1 = hdr['CRVAL1'], cdelt1 = hdr['CDELT1'])
    iraf.rspectext('%s.err.txt' % out_name, '%s.err.fits' % out_name, 
            crval1 = hdr['CRVAL1'], cdelt1 = hdr['CDELT1'])
    iraf.rspectext('%s.snr.txt' % out_name, '%s.snr.fits' % out_name, 
            crval1 = hdr['CRVAL1'], cdelt1 = hdr['CDELT1'])
    # add keywords
    f = pyfits.open('%s.spec.fits' % out_name)
    for key in hdr_save.keys():
        f[0].header.update(key, hdr_save[key])
    f[0].header.update('DATE-OBS', date_obs)
    f[0].header.update('OBSERVAT', observat)
    f[0].header.update('EPOCH', epoch)
    if one_side:
        # exposure time and velocity error are only well-defined for
        # data combined from a single side
        f[0].header.update('EXPTIME', exptime)
        f[0].header.update('FWHM', seeing/len(spec_list_fits))
        f[0].header.update('VERR', '%.2f' % np.sqrt(verr), 'Uncertainty in km/s')
        mjd += exptime/(2.*60.*60.*24.)
    else:
        # when combining sides, use the EXPTIME from the combined blue side
        f[0].header.update('EXPTIME', exptime_blue)
        del f[0].header['VERR']
    f[0].header.update('MJD', np.round(mjd, decimals=6))

    f.writeto('%s.spec.fits' % out_name, clobber=True)
    f.close()

def normalize_to_continuum(imgID, side='blue'):
    rootname = '%s%04d' % (side,imgID)
    hdr = pyfits.getheader('%s.spec.fits' % rootname)
    iraf.unlearn('continuum')
    iraf.continuum.order = 25
    iraf.continuum.high_rej = 5
    iraf.continuum.low_rej = 2
    iraf.continuum.niterate = 10
    iraf.continuum.type = "ratio"
    iraf.continuum('%s.spec.fits' % rootname, 'norm_%s.fits' % rootname, 
            interactive="no")
    iraf.unlearn('hedit')
    iraf.hedit('norm_%s.fits' % rootname, 'np2', hdr['NAXIS1'], 
            add="yes", update="yes", verify="no")
    #iraf.imcopy('norm_%s.fits' % rootname, 'norm_%s.imh' % rootname)

def stack_plot(spec_list, offset = False, alpha=1.):
    import matplotlib.pyplot as plt

    offset_val = 0.
    for spec in spec_list:
        dat = np.genfromtxt(spec, names='wave, flux', 
            dtype='f4, f4')
        plt.plot(dat['wave'],dat['flux']+offset_val,label = spec,alpha=alpha)
        if offset:
            offset_val -= np.median(dat['flux'])
        print spec
    plt.legend()
    plt.show()

def find_gain_readnoise(side='blue', aperture='1.0'):
    """Utility for measuring the gain and readnoise from raw images.
    
    Adaptation of iraf's findgain utility."""
    
    # high statistics section of the flat
    if side == 'blue':
        section="[80:350,600:800]"
    # TODO: confirm these sections
    elif side == 'red':
        if NEW_RED_SIDE:
            section="[1300:3000,50:400]"    
        else:
            section="[215:1000,50:265]"    
    flats = find_flats(aperture, side=side)
    biases = iraf.hselect('{}????.fits'.format(side), '$I', 
        'IMGTYPE == "bias" & EXPTIME == 0', Stdout=1)
    nmin = np.min([len(flats),len(biases)])
    gain_all = []
    rdnoise_all = []
    for i in range(nmin-1):
        gn, rn = compute_gain_readnoise(flats[i], flats[i+1],
            biases[i], biases[i+1], section = section)
        gain_all.append(gn)
        rdnoise_all.append(rn)
    #print gain_all, rdnoise_all    
    print 'Gain: %.3f %.3f' % (np.average(gain_all), np.std(gain_all, ddof=1))
    print 'Readnoise: %.2f %.2f' % (np.average(rdnoise_all), np.std(rdnoise_all, ddof=1))


def compute_gain_readnoise(flat1, flat2, zero1, zero2, section="[*,*]"):
    iraf.noao(_doprint=0)
    iraf.obsutil(_doprint=0)
    iraf.imarith(flat1, '-', flat2, 'flatdif')
    iraf.imarith(zero1, '-', zero2, 'zerodif')
    s = iraf.imstat('%s%s' % (flat1, section), fields="mean", nclip=20, 
        Stdout=1, format="no")
    mean_flat1 = np.float(s[0])
    s = iraf.imstat('%s%s' % (flat2, section), fields="mean", nclip=20, 
        Stdout=1, format="no")
    mean_flat2 = np.float(s[0])
    s = iraf.imstat('%s%s' % (zero1, section), fields="mean", nclip=20, 
        Stdout=1, format="no")
    mean_zero1 = np.float(s[0])
    s = iraf.imstat('%s%s' % (zero2, section), fields="mean", nclip=20, 
        Stdout=1, format="no")
    mean_zero2 = np.float(s[0])
    s = iraf.imstat('%s%s' % ('flatdif', section), fields="stddev", nclip=20, 
        Stdout=1, format="no")
    sigma_flatdif = np.float(s[0])
    s = iraf.imstat('%s%s' % ('zerodif', section), fields="stddev", nclip=20, 
        Stdout=1, format="no")
    sigma_zerodif = np.float(s[0])
    gain = (((mean_flat1 + mean_flat2) - (mean_zero1 + mean_zero2)) / 
        ((sigma_flatdif)**2 - (sigma_zerodif)**2))
    readnoise = gain * sigma_zerodif / np.sqrt(2)
    iraf.delete('flatdif.fits', verify="no")
    iraf.delete('zerodif.fits', verify="no")
    return gain, readnoise

def combine_sides_scombine(imgID_list_blue, imgID_list_red, output=None, splot='yes'):
    """imgID_lists are lists of numbers of extracted spectra:
    eg, [41] for red0041_flux.spec.fits"""

            
    blue_files = ['blue{:04d}_flux.spec.fits'.format(imgID) for imgID in imgID_list_blue]
    red_files = ['red{:04d}_flux.spec.fits'.format(imgID) for imgID in imgID_list_red]

    input_blue_list = ','.join(blue_files)

    if output is None:
        hdr = pyfits.getheader(blue_files[0])
        obj = hdr['OBJECT'].replace(' ','_')
        # create a unique name based on the input

        def ids_to_string(idlist):
            if len(idlist) == 1:
                return "{:d}".format(idlist[0])
            else:
                return "-".join(["{:d}".format(id) for id in idlist])

        output = obj + '_' + \
            ids_to_string(imgID_list_blue) + '+' + \
            ids_to_string(imgID_list_red) + '.fits'

    # clobber the old output file
    iraf.delete(output,verify='no')
    iraf.delete(output.replace('fits','txt'),verify='no')
    
    # determine dispersion: downsample to lower-resolution spectrum
    hdr = pyfits.getheader(blue_files[0])
    dw_blue = hdr['CDELT1']
    hdr = pyfits.getheader(red_files[0])
    dw_red = hdr['CDELT1']
    dw = np.max([dw_blue,dw_red])

    # cut off blue side redder than 5500
    trim_files = [f.replace('spec','trim') for f in blue_files]
    output_trim_list = ','.join(trim_files)
    iraf.unlearn('dispcor')
    iraf.dispcor.input = input_blue_list
    iraf.dispcor.output = output_trim_list
    # wavelength solution
    iraf.dispcor.w1 = 3800
    iraf.dispcor.w2 = 5500
    iraf.dispcor.dw = dw
    iraf.dispcor.flux = 'no'
    iraf.dispcor()

    all_files = trim_files + red_files
    input_file_list = ','.join(all_files)

    iraf.unlearn('scombine')
    iraf.scombine.input = input_file_list
    iraf.scombine.output = output
    iraf.scombine.combine = 'average'
    iraf.scombine.group = 'all'
    iraf.scombine.first = 'no'
    iraf.scombine.dw = dw
    iraf.scombine.scale = 'none'
    # attempt to join sides smoothly.
    #iraf.scombine.zero = 'median'
    #iraf.scombine.sample = ''
    iraf.scombine()

    iraf.wspectext(output, output.replace('fits','txt'), header="no")

    # clean up
    iraf.delete('blue*.trim.fits',verify='no')
    
    if splot == 'yes':
        iraf.splot(output)
