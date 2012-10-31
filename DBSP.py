#!/usr/bin/env python
from pyraf import iraf
import numpy as np
import os
import pyfits
import cosmics
from mpfit import mpfit
import copy
from glob import glob

# directory where the reduction code is stored
BASE_DIR = '/home/ebellm/observing/reduction/dbsp/'

# load IRAF packages
iraf.noao(_doprint=0)
iraf.imred(_doprint=0)
iraf.ccdred(_doprint=0)
iraf.kpnoslit(_doprint=0)
iraf.astutil(_doprint=0)
iraf.onedspec(_doprint=0)

# defaults
# (blue trace is usually around 250-260)
det_pars = {'blue':{'gain':0.72,'readnoise':2.5,'trace':253,
                    'crval':4345, 'cdelt':-1.072, 'arc':'FeAr_0.5.fits'}, 
            'red':{'gain':2.8,'readnoise':8,'trace':166,
                    'crval':7502, 'cdelt':1.530, 'arc':'HeNeAr_0.5.fits'}}

def mark_bad(side,numbers):
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

def matchSpectra(p, spectra=None, spectraErr=None, fjac=None):
    # Parameter values are passed in "p"
    # If fjac==None then partial derivatives should not be
    # computed.  It will always be None if MPFIT is called with default
    # flag.
    # Non-negative status value means MPFIT should continue, negative means
    # stop the calculation.
    status = 0
    res = np.sum(p*spectra[:, 1:], axis=1) - spectra[:, 0]
    totErr = np.sqrt(np.sum(spectraErr**2, axis=1))
    return [status, res/totErr]

def createArcDome(side='blue',trace=None,arcslit=0.5,overwrite=False):

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
    # update the headers
    iraf.asthedit('%s????.fits' % side, '/home/bsesar/opt/python/DBSP.hdr')

    # bias subtraction using the overscan
    filenames = glob("%s????.fits" % side)
    hdr = pyfits.getheader(filenames[0])
    iraf.unlearn('ccdproc')
    iraf.ccdproc.zerocor = "no"
    iraf.ccdproc.flatcor = "no"
    iraf.ccdproc.fixpix = "no"
#    iraf.ccdproc.fixfile = "../bluebpm"
    iraf.ccdproc.biassec = hdr['BSEC1']
    if side == 'blue':
        iraf.ccdproc.trimsec = "[%d:%d,*]" % (trace-100, trace+100)
    else:
        iraf.ccdproc.trimsec = "[*,%d:%d]" % (trace-100, trace+100)
    iraf.ccdproc.ccdtype = ""
    iraf.ccdproc.darkcor = "no"
    iraf.ccdproc.function = "spline3"
    iraf.ccdproc.order = 3
    iraf.ccdproc.niterate = 3
    iraf.ccdproc('%s????.fits' % side)


def fix_bad_column_blue():
    # find the bad column using a science exposure
    science = iraf.hselect('blue0*.fits', '$I', 'TURRET == "APERTURE" & LAMPS == "0000000"', Stdout=1)
    science = iraf.hselect(','.join(science), '$I', 'TURRET == "APERTURE" & LAMPS == "0000000" & AIRMASS != "1.000"', Stdout=1)
    f = pyfits.open(science[0])
    bad_column = f[0].data[1062,:].argmin() + 1
    f.close()
    f = open('bluebpm', 'w')
    f.write('%d %d 1 2835\n' % (bad_column, bad_column))
    f.close()
    iraf.fixpix('blue????.fits', "bluebpm")

def make_flats(side='blue',overwrite=False):
    # create dome flat images
    iraf.unlearn('flatcombine')
    iraf.flatcombine.ccdtype = ""
    iraf.flatcombine.process = "no"
    iraf.flatcombine.subsets = "no"
    iraf.flatcombine.rdnoise = "RON"
    iraf.flatcombine.gain = "GAIN"
    for aperture in ['0.5','1.0', '1.5', '2.0']:
        # find dome flat images
        domeflats = iraf.hselect('%s0*.fits' % side, '$I', 'TURRET == "APERTURE" & APERTURE == "%s" & LAMPS == "0000000"' % aperture, Stdout=1)
        domeflats = iraf.hselect(','.join(domeflats), '$I', 'TURRET == "APERTURE" & APERTURE == "%s" & LAMPS == "0000000" & AIRMASS == "1.000"' % aperture, Stdout=1)
        # find internal flat (incandescent lamp) images
        intflats = iraf.hselect('%s0*.fits' % side, '$I', 'TURRET == "LAMPS" & APERTURE == "%s" & LAMPS == "0000001"' % aperture, Stdout=1)
        intflats = iraf.hselect(','.join(intflats), '$I', 'TURRET == "LAMPS" & APERTURE == "%s" & LAMPS == "0000001" & AIRMASS == "1.000"' % aperture, Stdout=1)
        # dome flats are prefered over internal flats
        flats = []
        if (len(intflats) > 0) & (len(domeflats) == 0):
            flats = intflats
            print "Using %d internal flats for the %s arcsec slit." % (len(intflats), aperture)
        if len(domeflats) > 3:
            flats = domeflats
            print "Using %d dome flats for the %s arcsec slit." % (len(domeflats), aperture)
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
                iraf.response.niterate = 0
                iraf.response('temp', 'temp', 
                    'flat_%s_%s.fits' % (side, aperture), interactive="yes")
                iraf.delete('temp.fits', verify="no")

                # measure flat-field error from sigma images
                iraf.unlearn('imcombine')
                iraf.imcombine.reject = 'avsigclip'
                iraf.imcombine(','.join(flats), output='flat', sigma='sigma', scale='mode')
                iraf.imarith('sigma', '/', 'flat', 'frac')
                s = iraf.imstat('frac.fits', fields="mean", nclip=20, Stdout=1, format="no")
                print np.float(s[0])
                iraf.delete('flat.fits', verify="no")
                iraf.delete('sigma.fits', verify="no")
                iraf.delete('frac.fits', verify="no")
        else:
            print "No dome or internal flats for the %s arcsec slit." % aperture


def make_arcs_blue(slit=0.5, overwrite=False):
    # create the master arc with FeAr lamps
    aperture = "{:3.1f}".format(slit)

    iraf.unlearn('imcombine')
    iraf.imcombine.rdnoise = det_pars['blue']['readnoise']
    iraf.imcombine.gain = det_pars['blue']['gain']
    arcs = iraf.hselect('blue0*.fits', '$I', 'TURRET == "LAMPS" & APERTURE == "{aperture}" & LAMPS == "0100000"'.format(aperture=aperture), Stdout=1)
    try:
        arcs = iraf.hselect(','.join(arcs), '$I', 'TURRET == "LAMPS" & APERTURE == "{aperture}" & LAMPS == "0100000" & AIRMASS == "1.000"'.format(aperture=aperture), Stdout=1)
    except:
        pass
    if overwrite:
        iraf.delete('FeAr_{aperture}.fits'.format(aperture=aperture), 
            verify='no')
    iraf.imcombine(','.join(arcs), 'FeAr_{}'.format(aperture), reject="none")

def make_arcs_red(slit=0.5, overwrite=False):
    # create the master arc with HeNeAr lamps
    aperture = "{:3.1f}".format(slit)

    iraf.unlearn('imcombine')
    iraf.imcombine.rdnoise = det_pars['red']['readnoise']
    iraf.imcombine.gain = det_pars['red']['gain']
    arcs = iraf.hselect('red0*.fits', '$I', 'TURRET == "LAMPS" & APERTURE == "{aperture}" & LAMPS == "0001110"'.format(aperture=aperture), Stdout=1)
    try:
        arcs = iraf.hselect(','.join(arcs), '$I', 'TURRET == "LAMPS" & APERTURE == "{aperture}" & LAMPS == "0001110" & AIRMASS == "1.000"'.format(aperture=aperture), Stdout=1)
    except:
        pass
    if overwrite:
        iraf.delete('HeNeAr_{aperture}.fits'.format(aperture=aperture), 
            verify='no')
    iraf.imcombine(','.join(arcs), 'HeNeAr_{}'.format(aperture), reject="none")


def extract1D(imgID, side='blue', trace=None, arc=None, splot='no', redo='no', resize='yes', crval=None, cdelt=None):

    assert (side in ['blue','red'])
    assert (splot in ['yes','no'])
    assert (redo in ['yes','no'])
    assert (resize in ['yes','no'])

    if trace is None:
        trace = det_pars[side]['trace']

    if arc is None:
        arc = det_pars[side]['arc']

    if crval is None:
        crval = det_pars[side]['crval']

    if cdelt is None:
        cdelt = det_pars[side]['cdelt']

    # bias subtraction using the overscan
    hdr = pyfits.getheader('%s%04d.fits' % (side,imgID))
    iraf.unlearn('ccdproc')
    iraf.ccdproc.zerocor = "no"
    iraf.ccdproc.flatcor = "no"
    iraf.ccdproc.fixpix = "no"
    iraf.ccdproc.biassec = hdr['BSEC1']
    if side == 'blue':
        iraf.ccdproc.trimsec = "[%d:%d,*]" % (trace-100, trace+100)
    else:
        iraf.ccdproc.trimsec = "[*,%d:%d]" % (trace-100, trace+100)
    iraf.ccdproc.ccdtype = ""
    iraf.ccdproc.darkcor = "no"
    iraf.ccdproc.function = "spline3"
    iraf.ccdproc.order = 3
    iraf.ccdproc.niterate = 3
    iraf.ccdproc('%s%04d.fits' % (side,imgID), 
        flat="flat_%s_%s" % (side,hdr['APERTURE']), flatcor="yes")

    if (side == 'blue') and ('FIXPIX' not in hdr):
        iraf.fixpix('blue????.fits', "bluebpm")

    if 'OBSERVAT' not in hdr:
        # update the headers
        iraf.asthedit('%s%04d' % (side,imgID), '/home/bsesar/opt/python/DBSP.hdr')

    # remove cosmic rays with LA Cosmic
    if 'COSMIC' not in hdr and hdr['EXPTIME'] > 60:
        array, header = pyfits.getdata('%s%04d.fits' % (side,imgID), header=True)
        c = cosmics.cosmicsimage(array, gain=det_pars[side]['gain'], 
        readnoise=det_pars[side]['readnoise'], 
        sigclip = 4.5, sigfrac = 0.5, objlim = 2.0, satlevel=60000)
        c.run(maxiter = 3)
        header.update('COSMIC', 1, '1 if we ran LA Cosmic')
        pyfits.writeto('%s%04d.fits' % (side,imgID), c.cleanarray, header, clobber=True)

    # process the arc image
    if arc != 'FeAr_0.5.fits':
        iraf.ccdproc(arc)
        hdr_arc = pyfits.getheader(arc)

        if 'FIXPIX' not in hdr_arc:
            iraf.fixpix(arc, "%sbpm")

        if 'OBSERVAT' not in hdr_arc:
            # update the headers
            iraf.asthedit(arc, '/home/bsesar/opt/python/DBSP.hdr')

        # remove cosmic rays with LA Cosmic
        if 'COSMIC' not in hdr_arc:
            array, header = pyfits.getdata(arc, header=True)
            c = cosmics.cosmicsimage(array, gain=det_pars[side]['gain'], 
            readnoise=det_pars[side]['readnoise'], 
            sigclip = 4.5, sigfrac = 0.5, objlim = 2.0, satlevel=60000)
            c.run(maxiter = 3)
            header.update('COSMIC', 1, '1 if we ran LA Cosmic')
            pyfits.writeto(arc, c.cleanarray, header, clobber=True)

    # set up doslit
    fwhm = 4.6
    iraf.unlearn('doslit')
    iraf.unlearn('apslitproc')
    iraf.doslit.readnoise = "RON"
    iraf.doslit.gain = "GAIN"
    iraf.doslit.width = 3*fwhm
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
    iraf.doslit.b_order = 2
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
    fwhm_arc = 2.8 # input FWHM of arc lines here
    iraf.doslit.i_function = "legendre"
    iraf.doslit.i_order = 4
    if side == 'blue':
        #iraf.doslit.coordlist = "/home/bsesar/opt/python/brani_DBSP.lst"
        iraf.doslit.coordlist = BASE_DIR + 'dbsp_cal/brani_FeAr_dbsp.dat'
    else:
        #iraf.doslit.coordlist = "/home/bsesar/opt/python/henear.dat"
        iraf.doslit.coordlist = BASE_DIR + 'dbsp_cal/brani_HeNeAr_dbsp.dat'
    iraf.doslit.fwidth = fwhm_arc
    iraf.doslit.match = 10
    iraf.doslit.i_niterate = 3
    iraf.doslit.addfeatures = 'no'
    iraf.doslit.linearize = "yes"

    # extract 1D spectrum
    print arc, iraf.doslit.crval, iraf.doslit.cdelt
    iraf.doslit('%s%04d.fits' % (side,imgID), arcs=arc, splot=splot, redo=redo, resize=resize)

    # measure the position and width of sky lines (do this only for exposures longer than 3 min)
    if hdr['EXPTIME'] > 180:
        iraf.unlearn('scopy')
        iraf.delete('%s%04d.2001.fits' % (side,imgID), verify='no')
        iraf.scopy('%s%04d.ms.fits' % (side,imgID), '%s%04d' % (side,imgID), band=3, format="onedspec")
        iraf.unlearn('fitprofs')
        iraf.fitprofs.gfwhm = fwhm_arc
        iraf.fitprofs.nerrsample = 100
        iraf.fitprofs.sigma0 = 2.5
        iraf.fitprofs.invgain = 1.389
        iraf.delete('skyfit?.dat', verify='no')

        sky_lines = {'blue':
            {'wavelength':[4046.565, 4358.335, 5460.750, 5577.340],
            'regs':['4025 4070', '4340 4380', '5440 5480', '5560 5590']},
            'red':
            {'wavelength':[6300.304,6863.955,7340.885,7821.503,8430.174,8827.096],
            'regs':['6270 6320', '6840 6880', '7330 7355', '7810 7835','8420 8450','8800 8835']}}
        
        offsets = []
        for i in range(len(sky_lines[side]['wavelength'])):
            iraf.fitprofs( '%s%04d.2001.fits' % (side,imgID), 
                reg=sky_lines[side]['regs'][i], 
                logfile='skyfit_{:s}_{:1d}.dat'.format(side,i), 
                pos=BASE_DIR + 'dbsp_cal/skyline_{:s}_{:1d}.dat'.format(side,i), 
                verbose='no')
        #iraf.fitprofs( '%s%04d.2001.fits' % (side,imgID), reg='4025 4070', logfile='skyfit1.dat', pos='../skyline2.dat', verbose='no')
        #iraf.fitprofs('%s%04d.2001.fits' % (side,imgID), reg='4340 4380', logfile='skyfit2.dat', pos='../skyline4.dat', verbose='no')
        #iraf.fitprofs( '%s%04d.2001.fits' % (side,imgID), reg='5440 5480', logfile='skyfit3.dat', pos='../skyline3.dat', verbose='no')
        #iraf.fitprofs( '%s%04d.2001.fits' % (side,imgID), reg='5560 5590', logfile='skyfit4.dat', pos='../skyline.dat', verbose='no')

        # dump useful data from skyfit?.dat (center width err_center err_width)
        #for i in [1,2,3,4]:
            os.system('fgrep -v "#" skyfit_{side:s}_{num:1d}.dat |perl -pe "s/0\.\n/0\./g;s/^ +//;s/\(/ /g;s/\)/ /g;s/ +/ /g;" |cut -d" " -f1,6,8,13 > wavelength_offset_{side:s}_{num:1d}.dat'.format(side=side,num=i))
            dat = np.genfromtxt('wavelength_offset_{side:s}_{num:1d}.dat'.format(side=side,num=i) , usecols=(0,2), names="center, error")
            offsets.append(dat['center'] - sky_lines[side]['wavelength'][i])
            print offsets

        offset_final = np.mean(offsets)
        error_at_4750 = np.std(offsets, ddof=1)/np.sqrt(len(offsets))/4750*299792.458 # uncertainty in km/s at 4750 A
        
        # correct CRVAL1 value and add uncertainty in wavelength zero-point
        #lambda0 = np.array([4046.565, 4358.335, 5460.750, 5577.340])
        #dat = np.genfromtxt('wavelength_offset1.dat', usecols=(0,2), names="center, error")
        #offsets = np.zeros((dat.size,4))
        #for i in np.arange(4):
        #    dat = np.genfromtxt('wavelength_offset%d.dat' % int(i+1), usecols=(0,2), names="center, error")
        #    offsets[:, i] = dat["center"] - lambda0[i] # subtract this from CRVAL1

        #offset_final = np.zeros(offsets[:,0].size)
        #error_at_4750 = np.zeros(offsets[:,0].size)
        #for i in np.arange(offsets[:,0].size):
        #    offset_final[i] = np.average(offsets[i])
        #    error_at_4750[i] = np.std(offsets[i], ddof=1)/np.sqrt(4.)/4750*299792.458 # uncertainty in km/s at 4750 A

    # correct CRVAL1 and add uncertainty in km/s at 4750
    f = pyfits.open('%s%04d.ms.fits' % (side,imgID))
    hdr = f[0].header
    if hdr['EXPTIME'] > 180:
        f[0].header.update('VERR', '%.2f' % error_at_4750, 'Uncertainty in km/s at 4750 A')
    else:
        f[0].header.update('VERR', '-99.99', 'Uncertainty in km/s at 4750 A')
    f.writeto('%s%04d.ms.fits' % (side,imgID), clobber=True)
    f.close()

    # extract counts and uncertainties
    iraf.unlearn('scopy')
    iraf.delete('%s%04d.0001.fits' % (side,imgID), verify='no')
    iraf.delete('%s%04d.2001.fits' % (side,imgID), verify='no')
    iraf.delete('%s%04d.3001.fits' % (side,imgID), verify='no')
    iraf.scopy('%s%04d.ms.fits' % (side,imgID), '%s%04d' % (side,imgID), band=1, format="onedspec")
    iraf.scopy('%s%04d.ms.fits' % (side,imgID), '%s%04d' % (side,imgID), band=4, format="onedspec")
    # correct wavelength calibration using sky lines
    if hdr['EXPTIME'] > 180:
        iraf.specshift('%s%04d.0001' % (side,imgID), '%.3f' % (-offset_final))
        iraf.specshift('%s%04d.3001' % (side,imgID), '%.3f' % (-offset_final))

    # output to text files
    hdr_arc = pyfits.getheader('%s.ms.fits' % os.path.splitext(arc)[0])
    iraf.delete('%s%04d.spec.txt' % (side,imgID), verify='no')
    iraf.delete('%s%04d.err.txt' % (side,imgID), verify='no')
    iraf.delete('%s%04d.spec.fits' % (side,imgID), verify='no')
    iraf.delete('%s%04d.err.fits' % (side,imgID), verify='no')
    iraf.dispcor('%s%04d.0001' % (side,imgID), '%s%04d.spec' % (side,imgID), w1=hdr_arc['CRVAL1'], dw=hdr_arc['CDELT1'], nw=hdr_arc['NAXIS1'])
    iraf.wspectext('%s%04d.spec.fits' % (side,imgID), '%s%04d.spec.txt' % (side,imgID), header="no")
    iraf.dispcor('%s%04d.3001' % (side,imgID), '%s%04d.err' % (side,imgID), w1=hdr_arc['CRVAL1'], dw=hdr_arc['CDELT1'], nw=hdr_arc['NAXIS1'], blank=1.0)
    iraf.wspectext('%s%04d.err.fits' % (side,imgID), '%s%04d.err.txt' % (side,imgID), header="no")

    # calculate SNR
    iraf.delete('%s%04d.snr.fits' % (side,imgID), verify='no')
    iraf.imarith('%s%04d.spec.fits' % (side,imgID), '/', '%s%04d.err.fits' % (side,imgID), '%s%04d.snr.fits' % (side,imgID))

    # cleanup
    iraf.delete('skyfit?.dat', verify='no')
    iraf.delete('wavelength_offset?.dat', verify='no')
    iraf.delete('%s%04d.ms.fits' % (side,imgID), verify="no")
    iraf.delete('%s%04d.0001.fits' % (side,imgID), verify="no")
    iraf.delete('%s%04d.3001.fits' % (side,imgID), verify="no")

    # statistics
    hdr = pyfits.getheader('%s%04d.spec.fits' % (side,imgID))
    if hdr['EXPTIME'] > 180:
        print "Wavelengths are offset by %.3f A, zero-point uncertainty is %.2f km/s at 4750 A." % (offset_final, error_at_4750)
    wave1 = np.int(np.floor((3990 - hdr_arc['CRVAL1'])/hdr_arc['CDELT1']))
    wave2 = np.int(np.floor((4010 - hdr_arc['CRVAL1'])/hdr_arc['CDELT1']))
    s = iraf.imstat('%s%04d.snr.fits[%d:%d]' % (side,imgID, wave1, wave2), fields='mean', nclip=20, Stdout=1, format="no")
    print "SNR = %.1f at 4000 A" % np.float(s[0])

def combineSpectra(specList, outSpecFname, use_ratios=False):
    """Scales input 1D spectra onto the same scale
       and then combines spectra using a weighted mean.
    """

    # first spectrum in the list is always the reference spectrum
    hdr = pyfits.getheader('%s.spec.fits' % specList[0])
    mjd = hdr['MJD']
    date_obs = hdr['DATE-OBS']
    epoch = hdr['EPOCH']
    observat = hdr['OBSERVAT']
    exptime = hdr['EXPTIME']
    verr = np.float(hdr['VERR'])**2
    spec_ref = np.genfromtxt('%s.spec.txt' % specList[0], names='wave, flux', dtype='f4, f4')
    err_ref = np.genfromtxt('%s.err.txt' % specList[0], names='wave, flux', dtype='f4, f4')
    err_ref['flux'] = np.where(err_ref['flux'] <= 0, 1, err_ref['flux']) # reset bad error values to 1
    wave = spec_ref['wave']

    # spectra and their errors will be stored here
    spectra = np.zeros((spec_ref.size, len(specList)), dtype='f4')
    spectraErr = np.zeros((spec_ref.size, len(specList)), dtype='f4')

    spectra[:, 0] = spec_ref['flux']
    spectraErr[:, 0] = err_ref['flux']

    if use_ratios:
        ratio = [1]

    for i, fname in enumerate(specList[1:]):
        hdr = pyfits.getheader('%s.spec.fits' % fname)
        exptime += hdr['EXPTIME']
        verr += np.float(hdr['VERR'])**2
        spec = np.genfromtxt('%s.spec.txt' % fname, names='wave, flux', dtype='f4, f4')
        err = np.genfromtxt('%s.err.txt' % fname, names='wave, flux', dtype='f4, f4')
        err['flux'] = np.where(err['flux'] <= 0, 1, err['flux']) # reset bad error values to 1
        spectra[:, i+1] = spec['flux']
        spectraErr[:, i+1] = err['flux']
        if use_ratios:
            # use the 4200-4300 A region to determine te ratio of spectra
            good = np.where((spec['wave'] > 4200) & (spec['wave'] < 4300))
            ratio.append(np.median(spec_ref['flux'][good]/spec['flux'][good]))

    if use_ratios:
        ratio = np.array(ratio)
        print ratio
    else:
        # match spectra using least-squares
        fa = {'spectra':spectra, 'spectraErr':spectraErr}
        p0 = np.ones(len(specList)-1, dtype='f8')
        parbase={'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}
        parinfo=[]
        for i in range(len(p0)):
            parinfo.append(copy.deepcopy(parbase))
        m = mpfit(matchSpectra, p0,functkw=fa, parinfo=parinfo, quiet=True)
        if (m.status <= 0):
            print 'error message = ', m.errmsg
        ratio = np.append(1, m.params)
        print ratio

    spec_avg, sum_weights = np.average(spectra*ratio, weights=1./(spectraErr*ratio)**2, axis=1, returned=True)
    spec_err = 1./np.sqrt(sum_weights)
    # output coadded spectra and uncertainties
    f = open('coadd.txt', 'w')
    g = open('err.txt', 'w')
    h = open('snr.txt', 'w')
    for x, y, z in zip(wave, spec_avg, spec_err):
        f.write('%.3f %.3f\n' % (x, y))
        g.write('%.3f %.3f\n' % (x, z))
        h.write('%.3f %.3f\n' % (x, y/z))
    f.close()
    g.close()
    h.close()
    # save as 1D IRAF FITS files
    iraf.delete('%s.spec.fits' % outSpecFname, verify="no")
    iraf.delete('%s.err.fits' % outSpecFname, verify="no")
    iraf.delete('%s.snr.fits' % outSpecFname, verify="no")
    iraf.rspectext('coadd.txt', '%s.spec.fits' % outSpecFname, crval1 = hdr['CRVAL1'], cdelt1 = hdr['CDELT1'])
    iraf.rspectext('err.txt', '%s.err.fits' % outSpecFname, crval1 = hdr['CRVAL1'], cdelt1 = hdr['CDELT1'])
    iraf.rspectext('snr.txt', '%s.snr.fits' % outSpecFname, crval1 = hdr['CRVAL1'], cdelt1 = hdr['CDELT1'])
    # add EXPTIME and MJD keywords
    mjd += exptime/(2.*60.*60.*24.)
    f = pyfits.open('%s.spec.fits' % outSpecFname)
    hdr = f[0].header
    f[0].header.update('DATE-OBS', date_obs)
    f[0].header.update('MJD', np.round(mjd, decimals=6))
    f[0].header.update('EXPTIME', exptime)
    f[0].header.update('OBSERVAT', observat)
    f[0].header.update('EPOCH', epoch)
    f[0].header.update('VERR', '%.2f' % np.sqrt(verr), 'Uncertainty in km/s at 4750 A')
    f.writeto('%s.spec.fits' % outSpecFname, clobber=True)
    f.close()
    iraf.delete('snr.txt', verify="no")
    iraf.delete('coadd.txt', verify="no")
    iraf.delete('err.txt', verify="no")

def addRaDec(specFname):
    radec = np.genfromtxt('radec_targets', names='name, ra, dec, rad, decd', dtype='|S90, |S90, |S90, f8, f8')
    good = np.where(radec['name'] == specFname)[0]
    radec = radec[good]
    f = pyfits.open('%s.spec.fits' % specFname)
    f[0].header.update('RA', str(radec['ra'][0]))
    f[0].header.update('DEC', str(radec['dec'][0]))
    f[0].header.update('RAD', radec['rad'][0])
    f[0].header.update('DECD', radec['decd'][0])
    f[0].header.update('OBJECT', specFname)
    f.writeto('%s.spec.fits' % specFname, clobber=True)
    f.close()

def normalizeToContinuum(specFname):
    hdr = pyfits.getheader('%s.spec.fits' % specFname)
    iraf.unlearn('continuum')
    iraf.continuum.order = 25
    iraf.continuum.high_rej = 5
    iraf.continuum.low_rej = 2
    iraf.continuum.niterate = 10
    iraf.continuum.type = "ratio"
    iraf.continuum('%s.spec.fits' % specFname, 'n_%s.fits' % specFname, interactive="no")
    iraf.unlearn('hedit')
    iraf.hedit('n_%s.fits' % specFname, 'np2', hdr['NAXIS1'], add="yes", update="yes", verify="no")
    iraf.imcopy('n_%s.fits' % specFname, 'n_%s.imh' % specFname)

def velErrEstimate(specFname, N_samples=100):

    # load a spectrum and its uncertainty
    spec, hdr = pyfits.getdata('%s.spec.fits' % specFname, header=True)
    err = pyfits.getdata('%s.err.fits' % specFname)
    wave = hdr['CRVAL1'] + hdr['CDELT1']*np.arange(spec.size)

    # reset bad values
    bad = np.where((spec <= 0) | (err <= 0))
    spec[bad] = 0.0
    err[bad] = 1.0

    # create a mock sample of spectra
    if not os.access("samples", os.F_OK):
        os.mkdir("samples")
    g = open('mock.lst', 'w')
    for i in np.arange(N_samples):
#        mock = spec + np.abs(err*np.random.randn(err.size))
        mock = spec + np.random.poisson(err**2)
        f = open('temp.txt', 'w')
        for x, y in zip(wave, mock):
            f.write('%.2f %.2f\n' % (x,y))
        f.close()
        iraf.rspectext('temp.txt', 'samples/mock_%d.fits' % i, crval1 = hdr['CRVAL1'], cdelt1 = hdr['CDELT1'])
        g.write('samples/mock_%d.fits\n' % i)
    g.close()
    os.unlink('temp.txt')

    ## run FXCOR
    # load packages
    iraf.noao()
    iraf.rv()
    iraf.imutil()
    iraf.astutil()
    iraf.onedspec()

    # set up observatory
    iraf.unlearn('observatory')

    iraf.unlearn('keywpars')

    iraf.unlearn('continpars')
    iraf.continpars.order = 25
    iraf.continpars.high_rej = 5
    iraf.continpars.low_rej = 2
    iraf.continpars.niterate = 10

    iraf.unlearn('filtpars')
    #iraf.filtpars.cuton = 6
    #iraf.filtpars.fullon = 11
    #iraf.filtpars.cutoff = 568
    #iraf.filtpars.fulloff = 1138

    iraf.unlearn('fxcor')
    #iraf.fxcor.osample = '3920-3950,3960-3990,4090-4110,4330-4350,4850-4870'
    #iraf.fxcor.rsample = '3920-3950,3960-3990,4090-4110,4330-4350,4850-4870'
    #iraf.fxcor.osample = '4090-4110,4330-4350,4850-4870'
    #iraf.fxcor.rsample = '4090-4110,4330-4350,4850-4870'
    iraf.fxcor.peak = "yes"
    iraf.fxcor.height = 0.75
    iraf.fxcor.maxwidth = 2000.
    #iraf.fxcor.filter = "both"

    ## measure radial velocities using various templates
    iraf.fxcor.osample = '4092-4112,4331-4351,4851-4871'
    iraf.fxcor.rsample = '4092-4112,4331-4351,4851-4871'
    #iraf.fxcor("@objects.lst", "@/home/bsesar/projects/DBSP/F2_F9_templates_DBSP_1.0.lst", interactive="no", output = "largest_fxcor", rebin="largest")
    iraf.fxcor("@mock.lst", "@/home/bsesar/projects/DBSP/kathy_RR_elodie_templates.lst", interactive="no", output = "mock_fxcor", rebin="largest")

    dat = np.genfromtxt('mock_fxcor.txt', names='name, height, tdr, rv, rvErr', usecols=(0,6,8,9,12), dtype='|S90,f4,f4,f4,f4', comments='#')
    good = np.isfinite(dat['tdr'])
    dat = dat[good]
    sI = np.argsort(dat['tdr'])
    dat = dat[sI[::-1]]
    os.system('rm mock_fxcor.*')

    # loop over objects and calculate velocities
    n_templates = 50 # number of templates to average
    rv_stats = np.zeros(np.unique(dat['name']).size, dtype=[('name', '|S90'), ('rv_mean', 'f4'), ('rv_meanErr', 'f4'), ('rv_best', 'f4'), ('rv_bestErr', 'f4'), ('GSR_corr', 'f4')])
    for i,name in enumerate(np.unique(dat['name'])):
        hdr = pyfits.getheader(name, 0)
        good = np.where(dat['name'] == name)
        rv = (dat['rv'][good])[0:n_templates]
        rvErr = (dat['rvErr'][good])[0:n_templates]
        med = np.median(rv)
        rms = 0.741*(np.percentile(rv, 75) - np.percentile(rv, 25))
        good = np.where(np.abs(rv-med) < 3*rms)
        rv = rv[good]
        rvErr = rvErr[good]
        rv_mean, sum_weights = np.average(rv, weights=1./rvErr**2, returned=True)
        rv_meanErr = 1./np.sqrt(sum_weights)
        rv_stats['name'][i] = name
        rv_stats['rv_mean'][i] = rv_mean
        #rv_stats['rv_meanErr'][i] = np.average(rvErr)
        rv_stats['rv_meanErr'][i] = rv_meanErr
        rv_stats['rv_best'][i] = rv[0]
        rv_stats['rv_bestErr'][i] = rvErr[0]

    os.system("rm -rf samples")
    os.unlink("mock.lst")

    return np.std(rv_stats['rv_mean'], ddof=1)

if __name__ == '__main__':

    pass
    #import DBSP

    # create dome flats and the master arc
    #DBSP.createArcDome(overwrite=True)

    # process one science exposure
    #DBSP.extract1D_blue(61, redo="yes")
