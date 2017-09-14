"""
dbsp.py is a pyraf-based reduction pipeline for spectra taken with the
Palomar 200-inch Double Spectrograph.
"""

# Authors:  Eric Bellm <ebellm@caltech.edu>,
#           Branimir Sesar <bsesar@astro.caltech.edu>
# License: MIT.

from pyraf import iraf
import numpy as np
import os
import subprocess
import shutil
import inspect
import astropy.io.fits as pyfits
from scipy.optimize.minpack import leastsq
import copy
from glob import glob
import cosmics

# directory where the reduction code is stored
BASE_DIR = os.path.dirname(os.path.abspath(inspect.getfile(
                inspect.currentframe())))

def is_new_red_camera():
    """Utility for determining red camera version from image size."""
    ids = range(15)
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

    # raise ValueError('Could not locate red side files')
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

# use standard IRAF sproc, as it causes bugs for others: 
# use my modified sproc to avoid annoying doslit prompts: 
# export USEDBSPDOSLIT=1

if os.getenv('USEDBSPDOSLIT') is not None:
    #PkgName = iraf.curpack(); PkgBinary = iraf.curPkgbinary()
    iraf.reset(doslit = BASE_DIR+"/cl/doslit/")
    #iraf.task(doslitDOTpkg = 'doslit$doslit.cl', PkgName=PkgName,PkgBinary=PkgBinary)

# defaults
# (blue trace is usually around 250-260)
det_pars = {'blue': {'gain': 0.8, 'readnoise': 2.7, 'trace': 253,
                    'crval': 4345, 'cdelt': -1.072, 'arc': 'FeAr_0.5.fits',
                    'fwhm_arc': 2.8, 'pixel_size': 15.}}

if NEW_RED_SIDE:
    det_pars['red'] = {'gain': 2.8, 'readnoise': 8.5, 'trace': 166,
                       'crval': 7502, 'cdelt': 1.530, 'arc': 'HeNeAr_0.5.fits',
                       'biassec': '[4128:4141,1:440]', 'fwhm_arc': 2.4,
                       'pixel_size': 15.}
else:
    # old CCD
    det_pars['red'] = {'gain': 2.05, 'readnoise': 7.8, 'trace': 130,
                       'crval': 6600, 'cdelt': 2.46, 'arc': 'HeNeAr_0.5.fits',
                       'biassec': '[1080:1124,1:300]', 'fwhm_arc': 1.6,
                       'pixel_size': 24.}
                    # crval is in Angstrom, cdelt is Angstrom/pixel
                    # pixel size is in micron

def check_gratings_angles(ids=None):
    """Check header values for grating and angles and update dispersion values
    accordingly.

    Parameters
    ----------
    ids : list or array
        image numbers to check for header parameters
    """

    if ids is None:
        ids = range(15)
    for side in ['red','blue']:
        for id in ids:
            name = '{}{:04d}.fits'.format(side,id)
            if os.path.exists(name):
                hdr = pyfits.getheader(name)
                grating = np.int(hdr['GRATING'].split('/')[0])
                assert(grating in [158, 300, 316, 600, 1200])
                angle_toks = hdr['ANGLE'].split()
                deg = angle_toks[0]
                if len(angle_toks) > 2:
                    # ANGLE   = '27 deg 17 min'      / Grating Angle
                    min = angle_toks[2]
                else:
                    # ANGLE   = '39.2 deg'      / Grating Angle
                    min = 0.
                angle = np.float(deg) + np.float(min)/60.

                central_wavelength, dispersion = calculate_dispersion(
                    grating, angle, side=side)

                crval = np.round(central_wavelength)
                cdelt = dispersion / 1.E3 * det_pars[side]['pixel_size']

                if side == 'blue':
                    cdelt *= -1.

                det_pars[side]['crval'] = crval
                det_pars[side]['cdelt'] = cdelt

            else:
                continue

def calculate_dispersion(grating, angle, side='blue', order=1):
    """Calculate parameters needed to initialize dispersion solution.

    Modifies global variable det_pars with the computed values.
    See http://www.astro.caltech.edu/palomar/200inch/dbl_spec/dbspoverview.html
    for details.

    Parameters
    ----------
    grating : integer
        grating lines/mm
    angle : float
        grating angle in degrees
    side : {'blue' (default), 'red'}
        'blue' or 'red' to indicate the arm of the spectrograph
    order : int (default = 1)
        Spectral order

    Returns
    -------
    central_wavelength : float
        central wavelength in Angstroms
    dispersion : float
        dispersion in Angstroms / mm
    """

    MM_TO_ANGSTROM = 10000000
    INCHES_TO_MM = 25.4
    diffracted_angle = {'red':35.00, 'blue':38.50} # degrees
    focal_length = {'red':12.*INCHES_TO_MM, 'blue':9.*INCHES_TO_MM}

    line_spacing = 1./grating # mm

    theta_m = diffracted_angle[side] - angle
    central_wavelength = order * np.abs((np.sin(np.radians(theta_m)) - 
        np.sin(np.radians(angle))) * line_spacing * MM_TO_ANGSTROM) # Angstrom

    dispersion = (line_spacing * np.cos(np.radians(theta_m)) / 
        (focal_length[side] * order)) * MM_TO_ANGSTROM 

    return central_wavelength, dispersion

def mark_bad(imgID_list, side='blue'):
    """Utility for excluding specific files from further analysis.

    Saturated or mis-configured exposures are suffixed .bad so file searches
    do not find them.
    
    Parameters
    ----------
    imgID_list : list of ints or int
        image id(s) to be marked as bad.
    side : {'blue' (default), 'red'}
        'blue' or 'red' to indicate the arm of the spectrograph

    """

    assert (side in ['blue', 'red'])
    try:
        for num in imgID_list:
            name = '{:s}{:04d}.fits'.format(side, num)
            os.rename(name, name+'.bad')
    except TypeError:
        # single number
        num = imgID_list
        name = '{:s}{:04d}.fits'.format(side, num)
        os.rename(name, name+'.bad')

# run automatically
check_gratings_angles()
    

def create_arc_dome(side='both', trace=None, arcslit='0.5', overwrite=True):
    """Convenience function which subtracts bias, creates dome flats, and
    creates arc frames.
    
    Parameters
    ----------
    side : {'both' (default), 'blue', 'red'}
        'blue' or 'red' to indicate the arm of the spectrograph;
        'both' to reduce both
    trace : int
        row or column of the spectral trace, if different from default.
    arcslit : {'0.5','1.0', '1.5', '2.0'}
        string indicating the preferred arc slit width
    overwrite : boolean (default True)
        If True, overwrites any existing files; if False, skips processing.
    """

    assert ((side == 'both') or (side == 'blue') or (side == 'red'))

    if side == 'both':
        if len(glob('blue????.fits')) > 0:
            create_arc_dome(side='blue', trace=trace, arcslit=arcslit, 
                overwrite=overwrite)
        else:
            print 'No blue side files found.'
        if len(glob('red????.fits')) > 0:
            create_arc_dome(side='red', trace=trace, arcslit=arcslit, 
                overwrite=overwrite)
        else:
            print 'No red side files found.'
        return

    if trace is None:
        trace = det_pars[side]['trace']
    
    bias_subtract(side=side, trace=trace)

    if side == 'blue':
        fix_bad_column_blue()

    make_flats(side=side, overwrite=overwrite)

    if side == 'blue':
        make_arcs_blue(slit=arcslit, overwrite=overwrite)
    else:
        make_arcs_red(slit=arcslit, overwrite=overwrite)

def bias_subtract(side='blue', trace=None):
    """Use iraf.ccdproc to subtract bias from all frames.

    Parameters
    ----------
    side : {'blue' (default), 'red'}
        'blue' or 'red' to indicate the arm of the spectrograph
    trace : int
        row or column of the spectral trace, if different from default.
    """

    # update the headers
    iraf.asthedit('%s????.fits' % side, BASE_DIR + '/cal/DBSP.hdr')
    if side == 'blue':
        iraf.hedit('blue*.fits', 'DISPAXIS', 2, update="yes", 
            verify="no", add="yes", show="no")
    else:
        iraf.hedit('red*.fits', 'DISPAXIS', 1, update="yes", 
            verify="no", add="yes", show="no")

    # Need to define an instrument translation file in iraf 2.16.1
    iraf.unlearn('setinst')
    iraf.setinst.instrument = 'kpnoheaders'
    iraf.setinst.review = 'no'
    iraf.setinst.mode = 'h'
    iraf.setinst()

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
    domeflats = iraf.hselect('%s????.fits' % side, '$I', 'TURRET == "APERTURE" & APERTURE == "%s" & LAMPS == "0000000" & AIRMASS < 1.01 & IMGTYPE == "flat"' % aperture, Stdout=1)
    # find internal flat (incandescent lamp) images
    intflats = iraf.hselect('%s????.fits' % side, '$I', 'TURRET == "LAMPS" & APERTURE == "%s" & LAMPS == "0000001" & AIRMASS < 1.01' % aperture, Stdout=1)
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
    for aperture in ['0.5', '1.0', '1.5', '2.0']:
        flats = find_flats(aperture, side=side)
        if len(flats) > 0:
            if overwrite:
                iraf.delete('flat_%s_%s.fits' % (side, aperture), verify='no')
                if len(flats) < 3:
                    iraf.flatcombine(','.join(flats), output='temp', reject='pclip')
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
                    iraf.twodspec.longslit.dispaxis = 2
                else:
                    iraf.twodspec.longslit.dispaxis = 1
                iraf.response('temp[0]', 'temp[0]', 
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


def make_arcs_blue(slit='0.5', overwrite=False, zenith_only = True):
    """Creates the master FeAr arc with iraf.imcombine.

    Stores arc as FeAr_{slit}.fits.

    Parameters
    ----------
    slit : {'0.5' (default),'1.0', '1.5', '2.0'}
        string indicating the slit width
    overwrite : boolean
        If True, overwrites existing arc; if False, skips processing.
    zenith_only : boolean, default True
        If True, only use arcs taken at zenith.
        (Combining arcs taken at different elevations is not recommended
        due to flexure.)
    """

    aperture = slit

    iraf.unlearn('imcombine')
    iraf.imcombine.rdnoise = det_pars['blue']['readnoise']
    iraf.imcombine.gain = det_pars['blue']['gain']
    arcs = iraf.hselect('blue????.fits', '$I', 'TURRET == "LAMPS" & APERTURE == "{aperture}" & LAMPS == "0100000"'.format(aperture=aperture), Stdout=1)
    if zenith_only:
        try:
            arcs = iraf.hselect(','.join(arcs), '$I', 'TURRET == "LAMPS" & APERTURE == "{aperture}" & LAMPS == "0100000" & AIRMASS < 1.01'.format(aperture=aperture), Stdout=1)
        except:
            pass
    if overwrite:
        iraf.delete('FeAr_{aperture}.fits'.format(aperture=aperture), 
            verify='no')
    iraf.imcombine(','.join(arcs), 'FeAr_{}'.format(aperture), reject="none")

def make_arcs_red(slit='0.5', overwrite=False, zenith_only = True):
    """Creates the master HeNeAr arc with iraf.imcombine.

    Stores arc as HeNeAr_{slit}.fits.

    Parameters
    ----------
    slit : {'0.5' (default), '1.0', '1.5', '2.0'}
        string indicating the slit width
    overwrite : boolean
        If True, overwrites existing arc; if False, skips processing.
    zenith_only : boolean, default True
        If True, only use arcs taken at zenith.
        (Combining arcs taken at different elevations is not recommended
        due to flexure.)
    """

    aperture = slit

    iraf.unlearn('imcombine')
    iraf.imcombine.rdnoise = det_pars['red']['readnoise']
    iraf.imcombine.gain = det_pars['red']['gain']
    arcs = iraf.hselect('red????.fits', '$I', 'TURRET == "LAMPS" & APERTURE == "{aperture}" & LAMPS == "0001110"'.format(aperture=aperture), Stdout=1)
    if zenith_only:
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

    assert(flatcor in ['yes', 'no'])

    if trace is None:
        trace = det_pars[side]['trace']

    # Need to define an instrument translation file in iraf 2.16.1
    iraf.unlearn('setinst')
    iraf.setinst.instrument = 'kpnoheaders'
    iraf.setinst.review = 'no'
    iraf.setinst.mode = 'h'
    iraf.setinst()

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
        flat="flat_%s_%s" % (side, hdr['APERTURE']))

    if (side == 'blue') and ('FIXPIX' not in hdr):
        iraf.fixpix('blue????.fits', "bluebpm")

    if 'OBSERVAT' not in hdr:
        # update the headers
        iraf.asthedit(filename, BASE_DIR + '/cal/DBSP.hdr')

    # remove cosmic rays with LA Cosmic
    if remove_cosmics and ('COSMIC' not in hdr) and (hdr['EXPTIME'] > 60) and \
            (hdr['TURRET'] == 'APERTURE'):
        array, header = pyfits.getdata(filename, header=True)
        c = cosmics.cosmicsimage(array, gain=det_pars[side]['gain'], 
        readnoise=det_pars[side]['readnoise'], 
        sigclip = 4.5, sigfrac = 0.5, objlim = 2.0, satlevel=60000,
        skyOrder = 0, objectOrder = 0)
        c.run(maxiter = 3)
        header.set('COSMIC', 1, '1 if we ran LA Cosmic')
        pyfits.writeto(filename, c.cleanarray, header, clobber=True)

def store_standards(imgID_list, side='blue', trace=None, 
    arc=None, splot='no', redo='no', resize='yes', 
    caldir = "onedstds$iidscal/",
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
    arc : string
        Arc image to use if different from the default, e.g. 'HeNeAr_0.5.fits'
    splot : {'yes', 'no' (default)}
        Plot extracted spectra with iraf.splot?
    redo : {'yes', 'no' (default)}
        Redo the spectral extraction from scratch?  Passed to iraf.doslit.
        **Warning**--discards previous wavelength solution!  
        Use extract=False if you simply want to redefine calibration bandpasses
        and refit the sensitivity function.
    resize : {'yes' (default), 'no'}
        Resize the extraction aperture?  Passed to iraf.doslit.    
    caldir : string, default "onedstds$iidscal/"
        Directory to search for calibration standards in.
    crval : int or None (default)
        Spectrum reference dispersion coordinate, if different from default
        [Angstroms of central pixel]
    cdelt : int or None (default)
        Spectral dispersion, if different from default [Angstroms per pixel]
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
                telluric_cal_id = telluric_cal_id, quicklook = 'no')

    iraf.delete('std-{}'.format(side), verify='no')
    iraf.delete('sens-{}'.format(side), verify='no')

    iraf.unlearn('standard')
    iraf.standard.caldir = caldir
    iraf.standard.output = 'std-{}'.format(side)
    # use the tabulated bandpasses for the standards
    iraf.standard.bandwidth = "INDEF"
    iraf.standard.bandsep = "INDEF"
    # try these one at a time
    for imgID in imgID_list:
        # use the extracted spectrum!
        iraf.standard('%s%04d.spec.fits' % (side, imgID))

    iraf.unlearn('sensfunc')
    iraf.sensfunc.standards = 'std-{}'.format(side)
    iraf.sensfunc.sensitivity = 'sens-{}'.format(side)
    if side == 'blue':
        iraf.sensfunc.order = 3
    else:
        iraf.sensfunc.order = 6
    # varun says to use ignoreaps, but it's causing me problems downstream
    # iraf.sensfunc.ignoreaps = 'yes'
    iraf.sensfunc()

def estimateFWHM(imgID, side='blue'):
    """Use IRAF's imexam to measure the FWHM of the trace.
    
    Parameters
    ----------
    imgID : int
        File number of image to process, e.g., red0011.fits -> 11
    side : {'blue' (default), 'red'}
        'blue' or 'red' to indicate the arm of the spectrograph
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
        defkey = 'j'
    else:
        defkey = 'k'
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
    f[0].header.set('FWHM', np.round(fwhm, 2), 'FWHM estimate of the trace [pix]')
    f.writeto('%s%04d.spec.fits' % (side, imgID), clobber=True)
    f.close()
    if os.access('%s%04d_flux.spec.fits' % (side, imgID), os.F_OK):
        f = pyfits.open('%s%04d_flux.spec.fits' % (side, imgID))
        f[0].header.set('FWHM', np.round(fwhm, 2), 'FWHM estimate of the trace [pix]')
        f.writeto('%s%04d_flux.spec.fits' % (side, imgID), clobber=True)
        f.close()

def extract1D(imgID, side='blue', trace=None, arc=None, splot='no', 
        resize='yes', flux=True, telluric_cal_id=None, reextract=False, 
        sky_shift=True, redo='no', crval=None, cdelt=None, quicklook='no'):
    """Extract spectra for science objects and apply flux and telluric 
    corrections if requested.
    
    Wraps preprocess_image, iraf.doslit, iraf.calibrate, iraf.dispcor, 
    iraf.telluric, and some others.  

    Uses sky lines to correct the wavelength solution and estimate its 
    uncertainty.

    Generates text and fits spectra of the form {side}####.spec.{fits/txt}
    and uncertainties of the form {side}####.err.{fits/txt}.  Fluxed
    spectra and errors are stored as {side}####_flux.spec.{fits/txt}
    and {side}####_flux.err.{fits/txt} in units of erg/cm2/sec/Ang.

    Parameters
    ----------
    imgID : int
        File number of image to process, e.g., red0011.fits -> 11
    side : {'blue' (default), 'red'}
        'blue' or 'red' to indicate the arm of the spectrograph
    trace : int
        Row or column of the spectral trace, if different from default.
    arc : string
        Arc image to use if different from the default, e.g. 'HeNeAr_0.5.fits'
    splot : {'yes', 'no' (default)}
        Plot extracted spectra with iraf.splot?
    resize : {'yes' (default), 'no'}
        Resize the extraction aperture?  Passed to iraf.doslit.    
    telluric_cal_id : int or None (default)
        If defined, use specified image id to perform telluric correction.
    reextract : boolean (default False)
        Re-extract spectra?  If False, use existing aperture definitions
        if they exist.  If True, 
    sky_shift : boolean (default True)
        if True, use fit sky lines to adjust wavelength solution
    redo : {'yes', 'no' (default)}
        Redo the spectral extraction from scratch?  Passed to iraf.doslit.
        **Warning**--discards previous wavelength solution!  This is probably 
        not what you want!
        Use reextract=True if you simply want to redefine the apertures.
    crval : int or None (default)
        Spectrum reference dispersion coordinate, if different from default
        [Angstroms of central pixel]
    cdelt : int or None (default)
        Spectral dispersion, if different from default [Angstroms per pixel]
    quicklook : {'yes', 'no' (default)}
        Non-interactive aperture selection, tracing, and dispersion?  Passed
        to iraf.doslit.
    """

    assert (side in ['blue', 'red'])
    assert (splot in ['yes', 'no'])
    assert (redo in ['yes', 'no'])
    assert (resize in ['yes', 'no'])
    assert (quicklook in ['yes', 'no'])

    rootname = '%s%04d' % (side, imgID)

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

    # preprocess the arc image
    preprocess_image(arc, side=side, trace=trace, flatcor='no', 
        remove_cosmics=False)

    # set up doslit
    fwhm = 5.6
    iraf.unlearn('doslit')
    iraf.unlearn('sparams')
    iraf.unlearn('aidpars')

    iraf.doslit.quicklook = quicklook
    iraf.doslit.readnoise = "RON"
    iraf.doslit.gain = "GAIN"
    iraf.doslit.width = 1.4*fwhm
    iraf.doslit.crval = crval
    iraf.doslit.cdelt = cdelt
    iraf.sparams.nsum = 50
    iraf.doslit.clean = "yes"
    if side == 'blue':
        iraf.doslit.dispaxis = 2
    else:
        iraf.doslit.dispaxis = 1
    iraf.sparams.extras = "yes"
    iraf.sparams.lower = -1*int(round(iraf.doslit.width/2.))
    iraf.sparams.upper = int(round(iraf.doslit.width/2.))
    iraf.sparams.t_function = "legendre"
    iraf.sparams.t_niter = 3
    iraf.sparams.t_order = 4
    iraf.sparams.t_high = 2
    iraf.sparams.t_low = 2
    iraf.sparams.weights = "variance"
    iraf.sparams.b_order = 3
    iraf.sparams.b_niterate = 1
    iraf.sparams.select = "average"
    anullus_start = fwhm*2
    xL = np.floor(np.linspace(-80, -1*anullus_start, 10))
    xR = np.floor(np.linspace(anullus_start, 80, 10))
    background_range = ''
    for i in np.arange(xL.size-1):
        background_range += '%d:%d,' % (np.int(xL[i]), np.int(xL[i+1]-1))
    for i in np.arange(xR.size-1):
        background_range += '%d:%d,' % (np.int(xR[i]+1), np.int(xR[i+1]))
    iraf.sparams.b_sample = background_range[:-1]
    iraf.sparams.i_function = "legendre"
    iraf.sparams.i_order = 4
    if side == 'blue':
        iraf.sparams.coordlist = BASE_DIR + '/cal/FeAr_dbsp.dat'
        fwhm_arc = det_pars['blue']['fwhm_arc'] # input FWHM of arc lines here (in pixels)
    else:
        iraf.sparams.coordlist = BASE_DIR + '/cal/HeNeAr_dbsp.dat'
        fwhm_arc = det_pars['red']['fwhm_arc'] # input FWHM of arc lines here (in pixels)
    iraf.sparams.fwidth = fwhm_arc
    iraf.sparams.match = 10. # positive number is angstrom, negative is pix
    iraf.sparams.i_niterate = 5
    iraf.sparams.addfeatures = 'no'
    iraf.sparams.linearize = "yes"

    # extract 1D spectrum
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
    # normalize the response with a low-order spline
    iraf.unlearn('continuum')
    iraf.continuum.order = 5
    iraf.continuum.high_rej = 5
    iraf.continuum.low_rej = 2
    iraf.continuum.niterate = 10
    iraf.continuum.type = "ratio"
    iraf.continuum('d_trace_flat.fits', 'norm_d_trace_flat.fits',
                   interactive="no")

    # correct tellurics, if requested
    if telluric_cal_id is not None and side == 'red':
        tell_rootname = '%s%04d' % (side, telluric_cal_id)
        if not os.path.exists('norm_' + tell_rootname + '.fits'):
            normalize_to_continuum(telluric_cal_id, side=side)
        iraf.unlearn('telluric')
        iraf.telluric.input = rootname + '.ms.fits'
        iraf.telluric.output = ""
        iraf.telluric.sample = "6277:6288,6860:7000,7584:7678,9252:9842"
        iraf.telluric.interactive = "no"
        iraf.telluric.cal = 'norm_%s.fits' % tell_rootname
        iraf.telluric.ignoreaps = 'yes'
        iraf.telluric.xcorr = 'yes'
        iraf.telluric.tweakrms = 'yes'
        iraf.telluric.threshold = 0.01
        iraf.telluric()

    # measure shift with sky lines *before* fluxing to avoid floating point errors
    # measure the position and width of sky lines (do this only for exposures longer than 3 min)
    hdr = pyfits.getheader(rootname + '.fits')
    midpoint_loc = {'blue':4750,'red':7400}
    if hdr['EXPTIME'] > 120 and sky_shift:
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
            {'wavelength':[6923.21, 7340.885, 7821.51, 8430.174, 8885.83],
            'regs':['6917 6930', '7333 7348', '7813 7830', '8422 8439', '8875 8895']}}
        
        offsets = []
        for i in range(len(sky_lines[side]['wavelength'])):
            iraf.fitprofs(rootname + '.2001.fits',
                reg=sky_lines[side]['regs'][i], 
                logfile='skyfit_{:s}_{:1d}.dat'.format(side, i), 
                pos=BASE_DIR + '/cal/skyline_{:s}_{:1d}.dat'.format(side, i), 
                verbose='no')

        # dump useful data from skyfit?.dat (center width err_center err_width)
            os.system('fgrep -v "#" skyfit_{side:s}_{num:1d}.dat |perl -pe "s/0\.\n/0\./g;s/^ +//;s/\(/ /g;s/\)/ /g;s/ +/ /g;" |cut -d" " -f1,6,8,13 > wavelength_offset_{side:s}_{num:1d}.dat'.format(side=side, num=i))
            try:
                dat = np.genfromtxt('wavelength_offset_{side:s}_{num:1d}.dat'.format(side=side, num=i), usecols=(0, 2), names="center, error")
            except:
                # keep going if there's a bad profile fit
                print "Warning: bad fit for wavelength_offset_{side:s}_{num:1d}.dat".format(side=side, num=i)
                continue
            assert (dat['center'].size == 1)
            offsets.append(dat['center'] - sky_lines[side]['wavelength'][i])
        print offsets

        offset_final = np.mean(offsets)
        error_at_mid = np.std(offsets, ddof=1)/np.sqrt(len(offsets)) / \
            midpoint_loc[side]*299792.458 # uncertainty in km/s at 4750 A

    # add wavelength shifts/ uncertainty in km/s to headers
    # (CRVAL1 doesn't seem to apply the shift correctly?)
    f = pyfits.open(rootname + '.ms.fits')
    hdr = f[0].header
    if hdr['EXPTIME'] > 120 and sky_shift:
        f[0].header.set('WOFF', '%.2f' % offset_final, 'Wavelength offset from sky lines in A at {} A'.format(midpoint_loc[side]))
        f[0].header.set('VERR', '%.2f' % error_at_mid, 'Uncertainty in km/s at {} A'.format(midpoint_loc[side]))
    else:
        f[0].header.set('WOFF', '-99.99', 'Wavelength offset from sky lines in A at {} A'.format(midpoint_loc[side]))
        f[0].header.set('VERR', '-99.99', 'Uncertainty in km/s at {} A'.format(midpoint_loc[side]))
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
    f = pyfits.open(rootname + '.0001.fits')
    g = pyfits.open(rootname + '.3001.fits')
    hdr = f[0].header
    if hdr['EXPTIME'] > 120 and sky_shift:
        f[0].header.set('WOFF', '%.2f' % offset_final, 'Wavelength offset from sky lines in A at {} A'.format(midpoint_loc[side]))
        f[0].header.set('VERR', '%.2f' % error_at_mid, 'Uncertainty in km/s at {} A'.format(midpoint_loc[side]))
        g[0].header.set('WOFF', '%.2f' % offset_final, 'Wavelength offset from sky lines in A at {} A'.format(midpoint_loc[side]))
        g[0].header.set('VERR', '%.2f' % error_at_mid, 'Uncertainty in km/s at {} A'.format(midpoint_loc[side]))
        iraf.specshift(rootname + '.0001', '%.3f' % (-offset_final))
        iraf.specshift(rootname + '.3001', '%.3f' % (-offset_final))
    else:
        f[0].header.set('WOFF', '-99.99', 'Wavelength offset from sky lines in A at {} A'.format(midpoint_loc[side]))
        f[0].header.set('VERR', '-99.99', 'Uncertainty in km/s at {} A'.format(midpoint_loc[side]))
        g[0].header.set('WOFF', '-99.99', 'Wavelength offset from sky lines in A at {} A'.format(midpoint_loc[side]))
        g[0].header.set('VERR', '-99.99', 'Uncertainty in km/s at {} A'.format(midpoint_loc[side]))
    f.writeto(rootname + '.0001.fits', clobber=True)
    g.writeto(rootname + '.3001.fits', clobber=True)
    f.close()
    g.close()

    # normalize the spectrum and the uncertainty using the response function
    iraf.imarith(rootname + '.0001', '/', 'norm_d_trace_flat.fits', rootname + '.0001')
    iraf.imarith(rootname + '.3001', '/', 'norm_d_trace_flat.fits', rootname + '.3001')
    iraf.delete('*trace_flat.fits', verify="no")


    # flux, if requested
    if flux:
        iraf.unlearn('calibrate')
        # mode switch to make the input noninteractive
        iraf.calibrate.mode = 'h'
        iraf.calibrate.input = rootname+'.0001,'+rootname+'.3001'
        iraf.calibrate.output = rootname+'_flux.0001,'+rootname+'_flux.3001'
        # TODO: run setairmass; ensure extinction is set up correctly
        iraf.calibrate.extinction = ''
        # I'm not sure yet why this gets moved to .0001...
        iraf.calibrate.sensitivity = 'sens-{}.0001'.format(side)
        iraf.calibrate.ignoreaps = 'yes'
        iraf.calibrate()

    # output to text files
    def text_output(rootname,hdr_arc,flux=False):
        """Generate ASCII text spectra of the form rootname{_flux}.spec.txt and
        rootname{_flux}.err.txt.

        Parameters
        ----------
        rootname : string
            Base name of the input image
        hdr_arc : dict
            PyFITS-generated dictionary of the arc header
        flux : boolean (default False)
            Is the spectrum flux-corrected?
        """
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
    text_output(rootname, hdr_arc, flux=False)
    if flux:
        text_output(rootname, hdr_arc, flux=True)

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
    if hdr['EXPTIME'] > 120 and sky_shift:
        print "Wavelengths are offset by %.3f A, zero-point uncertainty is %.2f km/s at %.0f A." % (offset_final, error_at_mid, midpoint_loc[side])
    snr_loc = {'blue':4000,'red':7000}
    wave1 = np.int(np.floor((snr_loc[side]-10 - hdr_arc['CRVAL1'])/hdr_arc['CDELT1']))
    wave2 = np.int(np.floor((snr_loc[side]+10 - hdr_arc['CRVAL1'])/hdr_arc['CDELT1']))
    try:
        s = iraf.imstat(rootname + '.snr.fits[%d:%d]' % (wave1, wave2), 
                fields='mean', nclip=20, Stdout=1, format="no")
        print "SNR = %.1f at %d A" % (np.float(s[0]), snr_loc[side])
    except iraf.IrafError:
        print "Warning: could not imstat SNR"

    # measure FWHM of the trace
    estimateFWHM(imgID, side=side)

def combine_sides(imgID_list_blue, imgID_list_red, output=None, 
    save_sides=False, splot='yes'):
    """Downsample extracted blue and red spectra onto a common wavelength grid 
    and coadd, weighting by uncertainties.

    Takes output of extract1D and wraps iraf.dispcor and coadd_spectra.
    
    Parameters
    ----------
    imgID_list_blue : list of ints
        List of file numbers for blue spectra to combine
        e.g., blue0011.spec.fits and blue0015.spec.fits -> [11,15]
    imgID_list_red : list of ints
        List of file numbers for red spectra to combine
        e.g., red0011.spec.fits and red0015.spec.fits -> [11,15]
    output : string or None (default)
        File name of combined spectra.  If None, defaults to the OBJECT name
        of the first blue spectrum and the blue and red imgIDs.
    save_sides : boolean (default False)
        if True, save coadded sides with same formated names
    splot : {'yes' (default), 'no'}
        Plot the combined spectrum?
    """
            
    blue_files = ['blue{:04d}_flux.spec.fits'.format(imgID) for imgID in imgID_list_blue]
    red_files = ['red{:04d}_flux.spec.fits'.format(imgID) for imgID in imgID_list_red]
    blue_err_files = ['blue{:04d}_flux.err.fits'.format(imgID) for imgID in imgID_list_blue]
    red_err_files = ['red{:04d}_flux.err.fits'.format(imgID) for imgID in imgID_list_red]


    if output is None:
        hdr = pyfits.getheader(blue_files[0])
        obj = hdr['OBJECT'].replace(' ', '_')
        # create a unique name based on the input

        def ids_to_string(idlist):
            """Create hyphen-delimited string from list of ints"""
            if len(idlist) == 1:
                return "{:d}".format(idlist[0])
            else:
                return "-".join(["{:d}".format(id) for id in idlist])

        output = obj + '_' + \
            ids_to_string(imgID_list_blue) + '+' + \
            ids_to_string(imgID_list_red) 
        if save_sides:
            output_blue = obj + '_blue_' + \
                ids_to_string(imgID_list_blue) 
            output_red = obj + '_red_' + \
                ids_to_string(imgID_list_blue) 

    # clobber the old output files if they exist
    iraf.delete(output+'.*.fits', verify='no')
    iraf.delete(output+'.*.txt', verify='no')
    
    # determine dispersion: downsample to lower-resolution spectrum
    hdr = pyfits.getheader(blue_files[0])
    dw_blue = hdr['CDELT1']
    hdr = pyfits.getheader(red_files[0])
    dw_red = hdr['CDELT1']
    dw = np.max([dw_blue, dw_red])

    # find wavelength ranges
    def wavelength_range(fits_list):
        """Given an input list of fits spectra from extract1D, 
        return the absolute minimum and maximum wavelength ranges"""
        mins = []
        maxes = []
        for fname in fits_list:
            spec = np.genfromtxt(fname.replace('fits', 'txt'), names='wave, flux', 
                    dtype='f4, f4')
            mins.append(spec['wave'].min())
            maxes.append(spec['wave'].max())
        return [np.array(mins).min(), np.array(maxes).max()]

    blue_range = wavelength_range(blue_files)
    red_range = wavelength_range(red_files)

    # find overlap region
    if red_range[0] >= blue_range[1]:
        if save_sides:
            coadd_spectra(blue_files, output_blue)
            coadd_spectra(red_files, output_red)
            print('No overlap in wavelength solution between sides, exiting!')
            return

    # specify total spectral range
    w1 = blue_range[0]
    w2 = red_range[1]

    # re-disperse to common wavelength solution
    def redisperse_list(files,dw,w1,w2,key='spec'):
        """Re-disperse a list of spectra.

        Wraps iraf.dispcor.

        Parameters
        ----------
        files : list of strings
            Files to disperse
        dw : int
            Spectral dispersion [Angstroms per pixel]
        w1 : int
            Minimum wavelength [Angstrom]
        w2 : int
            Maximum wavelength [Angstrom]
        key : {'spec' (default) or 'err'}
            Whether the files are spectra or uncertainties.
        """
        input_list = ','.join(files)
        disp_files = [f.replace(key, key+'-disp') for f in files]
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
            iraf.wspectext(output, output.replace('fits', 'txt'), header="no")

        return disp_files

    # delete any lingering files
    iraf.delete('*-disp.fits', verify='no')
    iraf.delete('*-disp.txt', verify='no')

    blue_files_redisp = redisperse_list(blue_files, dw, w1, w2)
    red_files_redisp = redisperse_list(red_files, dw, w1, w2)
    blue_err_files_redisp = redisperse_list(blue_err_files, dw, w1, w2, key='err')
    red_err_files_redisp = redisperse_list(red_err_files, dw, w1, w2, key='err')

    # combine individual sides
    coadd_spectra(blue_files_redisp, 'tmp-blue')
    coadd_spectra(red_files_redisp, 'tmp-red')

    # for simplicity, just do this again if we're saving the sides
    if save_sides:
        coadd_spectra(blue_files_redisp, output_blue)
        coadd_spectra(red_files_redisp, output_red)

    # find optimum weighting between sides

    # combine sides, weighted by uncertainties
    coadd_spectra(['tmp-blue.spec.fits', 'tmp-red.spec.fits'], output,
        one_side=False)

    # clean up
    iraf.delete('*-disp.fits', verify='no')
    iraf.delete('*-disp.txt', verify='no')
    iraf.delete('tmp-*.fits', verify='no')
    iraf.delete('tmp-*.txt', verify='no')
    
    if splot == 'yes':
        iraf.splot(output+'.spec')

def match_spectra_leastsq(y, yref, yerr, yreferr):
    """Determine the best-fit ratio between two spectra.

    Parameters
    ----------
    y : array
        Spectrum
    yref : array
        Reference spectrum
    yerr : array
        Uncertainty in spectrum
    yreferr : array
        Uncertainty in reference spectrum
    """
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

def simple_coadd(imgID_list, side='blue'):
    """Utility for coadding spectra from a single side.
    
    Parameters
    ----------
    imgID_list : list of ints or int
        image id(s) to be coadded.
    side : {'blue' (default), 'red', 'both'}
        'blue' or 'red' to indicate the arm of the spectrograph

    """
    assert side in ('blue', 'red', 'both')
    if side == 'both':
        simple_coadd(imgID_list, side='blue')
        simple_coadd(imgID_list, side='red')
        return
    spec_list_fits = ['{}{:04d}_flux.spec.fits'.format(side,id) for id in imgID_list]
    out_name = side + '+'.join(['{:03d}'.format(i) for i in imgID_list]) + '_flux'
    coadd_spectra(spec_list_fits, out_name, scale_spectra=False)


def coadd_spectra(spec_list_fits, out_name, scale_spectra=True,
    use_ratios=False, ratio_range=[4200, 4300], 
    one_side=True):
    """Scales input 1D spectra onto the same scale multiplicatively
    and then combines spectra using an uncertainty-weighted mean.
       
    Parameters
    ----------
    spec_list_fits : list of strings
        List of spectra extracted by extract1D from a single side
    out_name : string
        Name of output file
    scale_spectra : boolean (default True)
        If True, fit and apply a multiplicative scale factor to match 
        the spectra before coadding
    use_ratios : boolean (default False)
        If True, use the median ratio of the spectra in the ratio_range
        to determine the scaling factor; otherwise match with least squares
        in overlap region
    ratio_range : list of ints 
        Two-element array of wavelengths (in Angstroms) where the ratio
        of the spectra are computed if scale_spectra and use_ratios are
        both True.
    one_side : boolean (default True)
        Are the data from a single side of the spectrograph?  If so,
        compute statistics for the coadd.
    """

    spec_list_txt = [f.replace('fits', 'txt') for f in spec_list_fits]

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
    err_ref = np.genfromtxt(spec_list_txt[0].replace('spec', 'err'), 
            names='wave, flux', dtype='f4, f4')
    wave = spec_ref['wave']
    spec_ref = spec_ref['flux'].view(np.ma.masked_array)
    err_ref = err_ref['flux'].view(np.ma.masked_array)


    # err_ref['flux'] = np.where(err_ref['flux'] <= 0, 1, err_ref['flux']) # reset bad error values to 1
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
        err = np.genfromtxt(fname_txt.replace('spec', 'err'), 
                names='wave, flux', dtype='f4, f4')
        spec = spec['flux'].view(np.ma.masked_array)
        err = err['flux'].view(np.ma.masked_array)
        # reset bad error values to 1
        # err['flux'] = np.where(err['flux'] <= 0, 1, err['flux']) 
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
        f[0].header.set(key, hdr_save[key])
    f[0].header.set('DATE-OBS', date_obs)
    f[0].header.set('OBSERVAT', observat)
    f[0].header.set('EPOCH', epoch)
    if one_side:
        # exposure time and velocity error are only well-defined for
        # data combined from a single side
        f[0].header.set('EXPTIME', exptime)
        f[0].header.set('FWHM', seeing/len(spec_list_fits))
        f[0].header.set('VERR', '%.2f' % np.sqrt(verr), 'Uncertainty in km/s')
        mjd += exptime/(2.*60.*60.*24.)
    else:
        # when combining sides, use the EXPTIME from the combined blue side
        f[0].header.set('EXPTIME', exptime_blue)
        try:
            del f[0].header['VERR']
        except KeyError:
            pass
    f[0].header.set('MJD', np.round(mjd, decimals=6))

    f.writeto('%s.spec.fits' % out_name, clobber=True)
    f.close()

def normalize_to_continuum(imgID, side='blue'):
    """Normalize a spectrum to the continuum.

    Wraps iraf.continuum.

    Parameters
    ----------
    imgID : int
        File number of image to process, e.g., red0011.fits -> 11
    side : {'blue' (default), 'red'}
        'blue' or 'red' to indicate the arm of the spectrograph
    """

    rootname = '%s%04d' % (side, imgID)
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
    # iraf.imcopy('norm_%s.fits' % rootname, 'norm_%s.imh' % rootname)

def read_spectrum(filename):

    def _load_fits_spectrum(filename):
        hdulist = pyfits.open(filename)
        hdr = hdulist[0].header
        val = hdulist[0].data
        crpix1 = hdr['CRPIX1']
        crval1 = hdr['CRVAL1']
        cd1_1 = hdr['CDELT1']
        spec_length = len(val)
        wave = cd1_1 * (np.arange(spec_length) - (crpix1-1)) + crval1
        return wave, val

    errfile = filename.replace('spec','err')

    if filename.endswith('txt'):
        dat = np.genfromtxt(filename, names='wave, flux', 
            dtype='f4, f4')
        if os.path.exists(errfile):
            errdat = np.genfromtxt(errfile, names='wave, err', 
                dtype='f4, f4')
        else:
            errdat = None

    elif filename.endswith('fits'):
        wave, flux = _load_fits_spectrum(filename)
        dat = {'wave':wave,'flux':flux}
        if os.path.exists(errfile):
            ewave, err = _load_fits_spectrum(errfile)
            errdat = {'wave':ewave,'err':err}
        else:
            errdat = None

    return dat, errdat




def stack_plot(spec_list, offset = False, alpha=1., normalize=False, 
    legend=True):
    """Plot several spectra on top of each other with matplotlib.
    Consider also iraf.specplot('spec1,spec2,spec3').

    Parameters
    ----------
    spec_list : list of strings
        List of text or fits spectra extracted by extract1D 
        (e.g., ['red0001_flux.spec.txt',red0002_flux.spec.fits'])
    offset : boolean (default False)
        Plot with a vertical offset between spectra?
    alpha : float, default 1.0 
        Opacity of the plot lines
    normalize : boolean (default False)
        Use a robust normalization for plotting multiple spectra on top of 
        each other when flux varies.
    legend : boolean (default True)
        Include a legend with filenames?
    """

    import matplotlib.pyplot as plt

    offset_val = 0.
    for spec in spec_list:
        dat, errdat = read_spectrum(spec)
        if normalize:
            pcts = np.percentile(dat['flux'],[5,95])
            dat['flux'] /= ((pcts[1]-pcts[0])/0.9)

        plt.plot(dat['wave'], dat['flux']+offset_val, label = spec, alpha=alpha)
        if offset:
            offset_val -= np.median(dat['flux'])
        print spec
    if legend:
        plt.legend()
    plt.show()

def find_gain_readnoise(side='blue', aperture='1.0'):
    """Utility for measuring the gain and readnoise from raw images.
    
    Adaptation of iraf's findgain utility.

    Parameters
    ----------
    side : {'blue' (default), 'red'}
        'blue' or 'red' to indicate the arm of the spectrograph
    aperture : {'0.5','1.0', '1.5', '2.0'}
        string indicating the slit width
    """
    
    # high statistics section of the flat
    if side == 'blue':
        section = "[81:350,600:800]"
    # TODO: confirm these sections
    elif side == 'red':
        if NEW_RED_SIDE:
            section = "[1300:3000,50:400]"    
        else:
            section = "[215:1000,50:265]"    

    flats = find_flats(aperture, side=side)
    biases = iraf.hselect('{}????.fits'.format(side), '$I', 
        'IMGTYPE == "bias" & EXPTIME == 0', Stdout=1)
    nmin = np.min([len(flats), len(biases)])
    gain_all = []
    rdnoise_all = []
    for i in range(nmin-1):
        gn, rn = compute_gain_readnoise(flats[i], flats[i+1],
            biases[i], biases[i+1], section = section)
        gain_all.append(gn)
        rdnoise_all.append(rn)
    # print gain_all, rdnoise_all    
    print 'Gain: %.3f %.3f' % (np.average(gain_all), np.std(gain_all, ddof=1))
    print 'Readnoise: %.2f %.2f' % (np.average(rdnoise_all), np.std(rdnoise_all, ddof=1))


def compute_gain_readnoise(flat1, flat2, zero1, zero2, section="[*,*]"):
    """DEPRECATED--use find_gain_readnoise instead.

    Utility for measuring the gain and readnoise from raw images.
    
    Adaptation of iraf's findgain utility.

    Parameters
    ----------
    flat1, flat2 : string
        Filenames of flat images
    zero1, zero2 : string
        Filenames of bias images
    section : string 
        Section of the files to compute stats on (e.g., "[1300:3000,50:400]")
    """

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
    """DEPRECATED--use combine_sides instead. 
    Downsample extracted blue and red spectra onto a common wavelength grid 
    and coadd using iraf.scombine.

    Takes output of extract1D and wraps iraf.dispcor and iraf.scombine.
    
    Parameters
    ----------
    imgID_list_blue : list of ints
        List of file numbers for blue spectra to combine
        e.g., blue0011.spec.fits and blue0015.spec.fits -> [11,15]
    imgID_list_red : list of ints
        List of file numbers for red spectra to combine
        e.g., red0011.spec.fits and red0015.spec.fits -> [11,15]
    output : string or None (default)
        File name of combined spectra.  If None, defaults to the OBJECT name
        of the first blue spectrum and the blue and red imgIDs.
    splot : {'yes' (default), 'no'}
        Plot the combined spectrum?
    """

            
    blue_files = ['blue{:04d}_flux.spec.fits'.format(imgID) for imgID in imgID_list_blue]
    red_files = ['red{:04d}_flux.spec.fits'.format(imgID) for imgID in imgID_list_red]

    input_blue_list = ','.join(blue_files)

    if output is None:
        hdr = pyfits.getheader(blue_files[0])
        obj = hdr['OBJECT'].replace(' ', '_')
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
    iraf.delete(output, verify='no')
    iraf.delete(output.replace('fits', 'txt'), verify='no')
    
    # determine dispersion: downsample to lower-resolution spectrum
    hdr = pyfits.getheader(blue_files[0])
    dw_blue = hdr['CDELT1']
    hdr = pyfits.getheader(red_files[0])
    dw_red = hdr['CDELT1']
    dw = np.max([dw_blue, dw_red])

    # cut off blue side redder than 5500
    trim_files = [f.replace('spec', 'trim') for f in blue_files]
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
    # iraf.scombine.zero = 'median'
    # iraf.scombine.sample = ''
    iraf.scombine()

    iraf.wspectext(output, output.replace('fits', 'txt'), header="no")

    # clean up
    iraf.delete('blue*.trim.fits', verify='no')
    
    if splot == 'yes':
        iraf.splot(output)

def batch_process(minID, maxID, side='blue', **kwargs):
    """Convenience function for reducing large numbers of consecutive spectra.

    Skips any missing files.

    Parameters
    ----------
    side : {'blue' (default), 'red', 'both'}
        'blue' or 'red' to indicate the arm of the spectrograph; 'both'
        to loop through both
    minID : int
        Minimum file number of image to process, e.g., red0011.fits -> 11
    maxID : int
        Maximum file number of image to process, e.g., red0020.fits -> 20
    Other keyword arguments (e.g., quicklook, flux) are passed to extract1D.
    """

    if side == 'both':
        sides = ['blue','red']
    else:
        sides = [side]
    for side in sides:
        for i in range(minID, maxID+1, 1):
            filename = '%s%04d.fits' % (side, i)
            if os.path.exists(filename):
                try:
                    extract1D(i, side=side, **kwargs)
                except iraf.IrafError:
                    # some errors just require you to try again...
                    print 'Hit error, retrying...'
                    extract1D(i, side=side, **kwargs)

def auto_join(blue_range, red_range, coadd_consecutive=True, do_joins=False):
    """Convenience function for joining large numbers of consecutive spectra.

    Skips any missing files.

    Warning! Not currently well-tested.

    Note that the state machine doesn't currently cover all possible cases;
    some pathological cases may cause failures, but only matched data will 
    ever be coadded.

    Parameters
    ----------
    blue_range: list
        Range of blue side file numbers to process, e.g., 
        blue0020.fits-blue0050 -> [20,50]
    red_range: list
        Range of red side file numbers to process, e.g., 
        red0020.fits-red0050 -> [20,50]
    coadd_consecutive: boolean, default True
        Should consecutive observations of the same source be joined?
    do_joins: boolean, default False
        Should the coadd_map be passed to combine_sides?
    """

    coadd_map = []
    current_obj_ra = ''
    current_blue = []
    current_red = []
    objects = []
    b = blue_range[0]
    r = red_range[0]
    while (b <= blue_range[1]) or (r <= red_range[1]):
        print b, r

        def load_hdrpars(i,side='blue'):
            filename = '%s%04d.fits' % (side, i)
            file_exists = os.path.exists(filename)
            if file_exists:
                hdr = pyfits.getheader(filename)
                # just use the nearest arcsec--small telescope drifts otherwise
                return True, hdr['OBJECT'],hdr['RA'][:8],hdr['DEC'][:9]
            else:
                return False, None, None, None

        bfileexists, bobj, bra, bdec = load_hdrpars(b,side='blue')
        rfileexists, robj, rra, rdec = load_hdrpars(r,side='red')

        if bfileexists and rfileexists and (bra == rra) and (bdec == rdec):
            # both sides observe same object
            if (rra == current_obj_ra) and coadd_consecutive:
                # which matches the previous object
                current_blue.append(b)
                current_red.append(r)
                current_obj = robj
            else:
                # both sides observe a new object
                if current_obj_ra != '': # starting the list
                    coadd_map.append((current_blue, current_red))
                current_obj = robj
                objects.append(current_obj)
                current_blue = [b]
                current_red = [r]
                current_obj_ra = rra
            b+=1
            r+=1
        else:
            # both sides observe different objects (or one side is missing)
            if rfileexists and (rra == current_obj_ra) and coadd_consecutive:
                current_red.append(r)
                r+=1
            elif bfileexists and (bra == current_obj_ra) and coadd_consecutive:
                current_blue.append(b)
                b+=1
            else:
                # some other state. save last object
                coadd_map.append((current_blue, current_red))
                objects.append(current_obj)

                # peek ahead
                _, nbobj, nbra, nbdec = load_hdrpars(b+1,side='blue')
                _, nrobj, nrra, nrdec = load_hdrpars(r+1,side='red')

                # does current blue match either of next objects?
                if bfileexists:
                    if (bra != nbra) and (bra != nrra):
                        # no--write it out by itself
                        coadd_map.append(([b],[]))
                        current_blue = []
                        objects.append(bobj)
                    else:
                        # save and continue
                        current_blue = [b]
                        current_obj = bobj
                b+=1

                # does current red match either of next objects?
                if rfileexists:
                    if (rra != nbra) and (rra != nrra):
                        # no--write it out by itself
                        coadd_map.append(([],[r]))
                        current_red = []
                        objects.append(robj)
                    else:
                        # save and continue
                        current_red = [r]
                        current_obj = robj
                        current_ra = rra
                r+=1

    # save final object
    coadd_map.append((current_blue, current_red))

    for x in zip(objects, coadd_map):
        print x
    if do_joins:
        for lists in coadd_map:
            combine_sides(lists[0], lists[1],splot='no')

    return coadd_map, objects

                
                



def sync(raw='./raw'):
    """Convenience routine for on-the-fly reduction that copies new files 
    from the raw directory into the current working directory without
    overwriting existing files.

    Because the pipeline modifies images in place, we need the --ignore-existing
    argument to rsync to avoid reverting our processed images to the originals.
    The -d argument allows us to decend into the specified subdirectory.
    """

    subprocess.call(['rsync','-d','--ignore-existing',raw+'/','.'])
