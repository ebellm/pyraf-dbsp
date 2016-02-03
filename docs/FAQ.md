### FAQ

*Why isn't autoidentify putting the arc lines in the right place?  Why is my wavelength solution bad?*

Are the crval and cdelt values appropriate for the CCD, grating, and angle you're using?   The `check_gratings_angles()` function attempts to automatically determine the correct values, but it is not foolproof.

You can use [this calculator](http://www.astro.caltech.edu/cgi-bin/grangle3.cgi) to determine approximate `crval` (center wavelength in angstroms) and `cdelt` values (disperson in Angstroms/pixel) for your grating. Input the grating and side you are using and a guess for the center wavelength; click calculate and note the reported grating angle. Edit your center wavelength guess until you find your grating angle.  The calculator then gives the center wavelength and dispersion you want.

reset to new values with:

	det_pars['red']['crval'] = 7330
	det_pars['red']['cdelt'] = 3.022

*I don't have any 0.5" slit arcs--what can I do?*

Choose another slit width and set it as the default:

    # use 1.0 arcsecond slit for arcs
    create_arc_dome(arcslit='1.0')

    # set default arc files
    det_pars['blue']['arc']  = 'FeAr_1.0.fits'
    det_pars['red']['arc']  = 'HeNeAr_1.0.fits



*I'm getting `"ERROR (1, "image keyword AIRMASS not found")"` when I run `create_arc_dome()`.*

One of your images is missing a header keyword--find it and mark it bad.

*I'm not happy with the fluxing--what can I do?*

Run `iraf.sensfunc()` and tweak to test, then repeat your `extract1D()` call.

*I'm getting an error:*

	ERROR: An unexpected error occurred while tokenizing input
	The following traceback may be corrupted or invalid
	The error message is: ('EOF in multi-line statement', (11, 0))

 	File "<tokenize>", line 2
    	if (Vars.dispcor and Vars.fluxcal1):
	IndentationError: unindent does not match any outer indentation level
	
This error is mysterious; it pops up irregularly, particularly if a routine has aborted unnaturally.  Try restarting `ipython` and running your command again, but be aware that in some cases it may be necessary to wipe the directory and restart from the raw images.

