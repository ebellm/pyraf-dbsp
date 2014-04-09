
Quick Start/Command Summary

Complete night reduction:

cd to your data directory
mdkir raw
cp *.fits raw
	(This code will overwrite your files--you must save the originals
	first!!!)

ipython

%run /path/to/dbsp.py
(Caltech astro users can do:
%run /home/ebellm/observing/reduction/dbsp/dbsp.py

mark_bad([47,49,50],side='blue')
mark_bad([35],side='red')

create_arc_dome()

# for telluric correction, if desired
extract1D(77,side='red')

store_standards([41,42,43], side='blue')
store_standards([41,42,43], side='red', telluric_cal_id = 77)

extract1D(61,side='blue',flux=True)
extract1D(63,side='red',flux=True, telluric_cal_id = 77)

combine_sides([61],[63])

Tips for on the fly reduction:

Start reducing data after you have taken your first standard star exposure

Copy (or rsync) new data into your raw subdirectory; then call sync() to 
	bring the new files into your existing directory without overwriting
	those you've already processed.



cd to your data directory
mdkir raw
cp *.fits raw
	(This code will overwrite your files--you must save the originals
	first!!!)

ipython

%run /home/ebellm/observing/reduction/dbsp/dbsp.py


change the names of any files you don't want to process (use your log):
	especially focus/test exposures from the beginning of the night
mark_bad([41,43,50],side='blue')

create_arc_dome(side = 'blue')

	fit normalization spectrum for temp interactively?  yes
		
		j to get residuals
		d to delete any weird spikes (symmetric ringing okay)
		f to refit
		q to quit

		> click in graphics window and hit ? for a list of commands
		(see response section of
		http://iraf.noao.edu/tutorials/doslit/doslit.html)
		This function fitting step uses a standard IRAF tool called icfit for
		interactive curve fitting. Within this mode you can type '?' to get a
		list of the commands. Typically one would use only :function and :order
		to change the function and order of the fit, and 'f' to cause the fit
		to be redone after changing any fitting parameters. Feel free to
		experiment and when you are done and have an acceptable fit exit with
		'q'. 

		(a spline fit is standard, so typing 'no' in the terminal is probably
		okay) (note that commands in the graphics window typically require you
		to click to focus; then you may have to click back to type responses in
		the terminal.) 


create_arc_dome(side = 'red')

# look at your flats in ds9 before continuing to make sure they don't have any
# weird features!

# for telluric correction, if desired
extract1D(77,side='red')

store_standards([41,42,43], side='blue')
store_standards([41,42,43], side='red',telluric_cal_id = 77)


extract1D(61,side='blue',flux=True)

	options: 
		redo: 
			generates new wavelength solution.  Otherwise, use stored
			version for specified arc file
		arc:
			specifies arc file to use for wavelength solution.  
			Default is 0.5" FeAr (blue) and HeNeAr (red)
		crval, cdelt:
			specifies central wavelengths and dispersions.  Loads default
			values for the "standard" PTF setup [TODO: describe], but 
			it is crucial to specify these correctly if using a different
			setup.
			
	edit apertures for (file)? 
		'yes', in terminal
		again, see http://iraf.noao.edu/tutorials/doslit/doslit.html
		(http://stsdas.stsci.edu/gethelp/HelpSys.html better for parameter name
		search)
		(http://www.twilightlandscapes.com/IRAFtutorial/IRAFintro_06.html)

		d to delete trace, m to set it
		b to enter background editing
			z to delete background intervals
			s s (with cursor positions) to mark new fit regions
			f to fit
			q to quit

	fit traced positions interactively?
		'yes' in window--brings up icfit again
		? for help
		:order # to change fit order.  default of 4 is probably fine; don't
		worry about excursions on the faint ends of the trace
		d to delete any points biasing the fit
		s (twice) to change the sampling region (eg to exclude faint ends)
		f to refit
		(aiming for RMS < 0.07 or so)

	autoidentify, reidentify
		get a plot with labelled lines
		f to fit
		j to switch to residuals plot
		d to delete outliers
			blue side can be bad below 4000:  delete all the points with d,
				hit q, f, l (to reidentify lines)
			good RMS is < 0.07 or so
		q to return to spectrum


		w x to zoom in
		m to specify a line (with cursor hovering over it): enter the wavelength in angstroms
		w a to zoom out
		do it for a second line (may already be identified)
		f to fit
		d to delete outliers
		l or y to ask for more lines
		f to fit again
		q twice to save solution

		when reidentifying, bad RMS may often be fixed just by hitting l then q

	change wavelength coordinate assignments?
		sets wavelength range & binning
			sensible defaults: 5500-10000, 1.525 red (new)
							   5500-7800, 2.47 red (old)
							   3700-5700, 1.07 blue
	
	standards:
		? for list
			or iraf.page('onedstds$README')
		all of the commonly-used ones are in onedstds#iidscal:
			g191b2b
			feige34
			bd332642
			bd284211
		edit bandpasses:
			(choose smooth regions)
			a a (with mouse pointer at two positions) to place new bands
			d to delete them
			q to quit and save

		sensfunc is the fitting function
			? for help
			s over graphs to eliminate mean shifts due to non-photometric 
				conditions (toggles)
			d to delete bad points
			make sure the fitted function doesn't go up after the last
			points--it will blow up the noise.  Also consider decreasing
			the order of the fit (:order 4) to avoid spline artifacts

combine_sides([64],[71,72])

		


colon comands: click in grey bar at bottom of graphics window
help: click graphics window, hit ?
window commands: click graphics window, type w ?

iraf.help('doslit')
iraf.dir('onedstds$')
iraf.type('onedstds$README') (or iraf.page)

FAQ:

Why isn't autoidentify putting the arc lines in the right place?  Why is my
wavelength solution bad?

Are the crval and cdelt values appropriate for the CCD, grating, and angle
you're using?  Defaults are:

You can use http://www.astro.caltech.edu/cgi-bin/grangle3.cgi to determine
approximate crval and cdelt values for your grating. 
Input the grating and side you are using and a guess for the center
wavelength; click calculate and note the reported grating angle.  
Edit your center wavelength guess until you find your grating angle.  
The calculator then gives the center wavelength and dispersion you want.


I'm getting "ERROR (1, "image keyword `AIRMASS' not found")" when I run
create_arc_dome.

One of your images is missing a header keyword--find it and mark it bad.

I'm not happy with the fluxing--what can I do?

run iraf.sensfunc() and tweak to test, then repeat your extract1D call.



--

KNOWN PROBLEMS/LIMITATIONS:
	absorption features in the "featureless" telluric calibration will 
	become emission features in the corrected spectrum

	wavelength solution will fail if gratings/angles other than "PTF standard"
	are used--should at least warn the user if that is the case

	fluxing is not correcting for extinction

	only one object can be extracted from a given image

	some fluxing weirdness at the shortest blue wavelengths

todos:  
add a license
fluxing should correct for extinction
provide more versatile code for using other gratings/angles/, or at least
	specifying crval/cdelt [could just modify the dict directly]
figure out how to turn off extra "Enters" where possible
improve telluric correction defaults
improve robustness of joining red and blue sides (bad fluxing kills it)
wrapper for single extract-combine-plot run
script to "undo" various parts of the analysis?  eg, start from scratch w/
	standards
brani suggests only using arcs taken in a single batch--need to adjust code
	logic

(calculate gain & readnoise from cal files)

---
