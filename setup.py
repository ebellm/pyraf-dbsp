from distutils.core import setup

setup(
    name='DBSP',
    version='0.2.1',
    author='Eric Bellm'
    author_email='ebellm@caltech.edu',
    packages=['dbsp'],
    url='http://pypi.python.org/pypi/DBSP/',
    license='LICENSE.txt',
    description='Pyraf-based spectroscopic reduction pipeline for the Double Spectrograph on the Palomar 200-inch telescope.',
    long_description="""A PyRAF-based reduction pipeline for spectra taken with the Palomar 200-inch Double Spectrograph.

	This pipeline provides a simplified interface for basic reduction of single-object spectra with minimal overhead.  It is suitable for quicklook classification of transients as well as moderate-precision (few km/s) radial velocity work.""",
    install_requires=[
        "pyraf >= 2.0",
        "numpy >= 1.5",
		"astropy >= 1.0",
		"scipy >= 0.10.0"
    ],
)
