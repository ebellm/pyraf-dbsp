### Quick Start/Command Summary

Complete night reduction:

`cd` to your data directory.

Save a copy of your files in a "raw" directory--this code overwrites the originals!
	
	mdkir raw
	cp *.fits raw


Start `ipython` and load the script:

	%run /path/to/dbsp.py

For users on the Caltech astro network, log in to soroban and execute:

	export PATH="/scr/ebellm/anaconda/bin:$PATH"
	source activate iraf27
	mkiraf  # choose xgterm
	ipython 
	%run /home/ebellm/observing/reduction/dbsp/dbsp.py

Exclude any images that you don't want to analyze (use your log; especially focus/test exposures from the beginning of the night):

	mark_bad([47,49,50],side='blue')
	mark_bad([35],side='red')

Create arcs and dome flats:

	create_arc_dome()

Process flux calibration standards, if desired.  If not, skip this step and set `flux=False` in `extract1D()`.

	store_standards([41,42,43], side='blue')

Extract data:

	extract1D(61,side='blue')

For basic telluric correction on the red side, first extract an appropriate telluric calibrator, then pass it to `store_standards` and `extract1D`:

	extract1D(77,side='red', Flux=False)

	store_standards([41,42,43], side='red', telluric_cal_id = 77)

	extract1D(63,side='red',flux=True, telluric_cal_id = 77)

To process a large number of science spectra in a row:
	
	batch_process(20, 45, side='blue', quicklook='no')		
Finally, join spectra from the blue and red sides, identify pairs (or more) of images:

	combine_sides([61],[63,64])
