### Known Problems & Limitations
* expects FeAr and (simulataneous) HeNeAr arcs in the 0.5 arcsec slit
* any absorption features in the "featureless" telluric calibration will become emission features in the corrected spectrum
* fluxing is not correcting for extinction?
* only one object can be extracted from a given image
* some fluxing weirdness at the shortest blue wavelengths
* joining red and blue sides with combine_sides often introduces artifacts

### Todos  
* fluxing should correct for extinction
* improve telluric correction defaults
* improve robustness of joining red and blue sides (bad fluxing kills it)
* wrapper for single extract-combine-plot run
* script to "undo" various parts of the analysis?  eg, start from scratch w/ standards
* Brani suggests only using arcs taken in a single batch--need to adjust code logic
* automate wavelength coordinate assignment for autoidentify?
* calculate gain & readnoise from cal files
* set fwhm_arc correctly (currently only uses a default)
* write auto-joiner to speed joining and coaddition
* write a log file to record processing steps taken
* add pipeline version and reducer information to output headers

### Revisions

#### 0.2.1 dev
(ongoing)

Bug fixes:

* allow check_gratings_angles to handle decimal ANGLE keywords
* handle bug where flatcombine expected multiextension FITS 
* Provide smoother handling if create_arc_dome(side='both') is called when data from only one side are present


Enhancements:

* harmonize aperture formats across routines

#### 0.2.0 
April 18, 2014

Enhancements:

* detect and attempt to automatically estimate dispersion parameters for any grating and angle
* add quicklook reduction mode
* add batch_process interface 
* reduce repetitive prompting in `doslit` and `calibrate`
* add license
* expand documentation

Bug fixes:

* avoid instrument translation file problem in IRAF 2.16.1



#### 0.1.0 
July 3, 2013 

Initial public release.

