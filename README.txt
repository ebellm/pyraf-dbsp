todos:  
fit HeNeAr line list (add reddest lines; use c to find bad lines & comment
	them out?)
generate all arcs at once
flux spectra ("standards" in doslit)
combine red & blue sides
fix stats code for red side (quote at midpoint, 7500 for red; check if
	dispersion is included)

(calculate gain & readnoise from cal files)

---

cd to your data directory
mdkir raw
cp *.fits raw

ipython

%run ~/observing/reduction/dbsp/DBSP.py


change the names of any files you don't want to process (use your log):
mark_bad('blue',[41,43,50])

createArcDome(side = 'blue')

	fit normalization spectrum for temp interactively?
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

createArcDome(side = 'red')

store_standards([41,42,43], side='blue',redo='yes')
store_standards([41,42,43], side='red',redo='yes')
	for now, only use a single standard

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

	fit traced positions interactively?
		'yes' in window--brings up icfit again
		? for help
		:order # to change fit order.  default of 4 is probably fine; don't
		worry about excursions on the faint ends of the trace
		d to delete any points biasing the fit
		f to refit
		(aiming for RMS < 0.07 or so)

	autoidentify
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

	change wavelength coordinate assignments?
		sets wavelength range & binning
			sensible defaults: 5500-10000, 1.525 red (new)
							   3800-5700, 1.07 blue
	
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

combine_sides([64],[71,72],output='mytarget')

		


colon comands: click in grey bar at bottom of graphics window
help: click graphics window, hit ?
window commands: click graphics window, type w ?

iraf.help('doslit')
iraf.dir('onedstds$')
iraf.type('onedstds$README') (or iraf.page)
