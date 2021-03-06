Description
===========

Code to compute polyspectrum from a 3D grid following the method described in `Watkinson et al. (2017) <http://adsabs.harvard.edu/abs/2017MNRAS.472.2436W>`__.

**IMPORTANT:**: The current input file only supports powerspectra and bispectra. Test results analogue to Fig. 7 in `Watkinson et al. (2017) <http://adsabs.harvard.edu/abs/2017MNRAS.472.2436W>`__ can be found `here <https://github.com/annehutter/polyspectrum/blob/master/bispectra_tests.pdf>`__. Suggestions for  improvements (particularly run time) are very welcome and can be sent to `a.k.hutter@rug.nl <a.k.hutter@rug.nl>`__.

Installation
============

Pre-requisities
---------------

Serial run
``````````

1. fftw3 library: ``fftw3 >= 3.3.3``

Parallel run
````````````

1. MPI library
2. fftw3 & fftw3-mpi library: ``fftw3 >= 3.3.3``

FFTW3
'''''

Go to the `FFTW webpage <http://www.fftw.org/download.html>`__ to install fftw3. Ensure to compile the library with the ``enable-mpi`` flag for parallel runs
::
    
    $ ./configure --enable-mpi
    $ make
    $ make install
    
Note: To create the dynamic libraries, run configure with the ``--enable-shared`` flag. 


Download & Build
----------------

::

    $ git clone https://github.com/annehutter/polyspectrum.git
    $ make

This will download the code and first test case from the github directory and compile the source code.

Execution
---------

The first test case can then be run by
::

    $ ./polyspectrum iniFile.ini

``iniFile.ini`` contains all input parameters that are needed for any runs. For a different simulation the code does not need to be recompiled but just this parameter file iniFile.ini to be adapted.


Parameter file
``````````````

The following parameters are specified in ``iniFile.ini``:

**Grid**
''''''''

- ``gridsize``: size of the 3D grid along one axis
- ``boxsize``: comoving boxsize in Mpc/h

- ``gasInputsInDoublePrecision``: set to 0 for single, 1 for double precision for density fields being read
- ``ionInputsInDoublePrecision``: set to 0 for single, 1 for double precision for ionization fields being read
- ``densityFile``: path to 3D density grid
- ``ionFile``: path to 3D ionization grid

**Cosmology**
'''''''''''''

- ``hubble_h``: H = 100*hubble_h km/s/Mpc
- ``omega_b``: baryon density parameter
- ``omega_m``: matter density parameter
- ``omega_l``: lambda density parameter
- ``Y``: mass fraction of Helium in the primordial gas (assumed to consist of H and He)

**Polyspectrum**
''''''''''''''''

- ``whichField``: field from which the polyspectrum is calculated; options are: ``DENS`` for density, ``XHII`` for ionization fraction and ``XHI_DENS`` for neutral gas density
- ``useOnlyPhase``: ``1`` yields the phase-only bispectrum 
- ``n``: number of vectors, i.e. ``n=2`` yields powerspectrum, ``n=3`` yields bispectrum
- ``equilateral``: ``1`` yields bispectrum for equilateral triangles with ``numValues`` providing the number of k-values depending on the box and gridsize of the simulation box, for all values not equal to ``1`` non-equilateral triangles are computed
- ``k1``: length of first vector in Mpc/h
- ``k2``: length of second vector in Mpc/h (if n>2)
- ``numValues``: set to 1, if only single value for polyspectrum should be calculated; for >1 and ``n>2`` number of k bins along ``theta = 0 to 180``
- ``theta``: angle between ``k1`` and ``k2`` in rad
- ``kbinwidth``: bin width in k-space in units of grid cells for computing the bispectra
- ``kbinningCase``: ``DEFAULT`` assumes for all k vectors a binwidth in k-space of ``kbinwidth``, ``GIVEN_BINNING`` assumes for ``k1`` and ``k2`` a binwidth in k-space of ``kbinwidth`` and for ``k3`` a width that corresponds to ``cosThetaBinwidth``, ``DERIVED_BINNING`` assumes for ``k1`` and ``k2`` a binwidth in k-space of ``kbinwidth`` and for ``k3`` it derives the binwidth from the allowed ranges given by the uncertainties of ``k1`` and ``k2``
- ``cosThetaBinwidth``: binwidth in cos(theta)

**Output**
''''''''''

- ``output_dir``: directory where output should be written
- ``output_basename``: output name of the run
- ``write_numpolygons``: ``1`` if output files should contain number of polygons, otherwise they are not included in the output
