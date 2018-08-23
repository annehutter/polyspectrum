Description
===========

Code to compute polyspectrum from a 3D grid following the method described in Watkinson et al. (2017).

IMPORTANT
=========

This is work in progress and the code has not been finalized. Suggestions can be sent to `a.k.hutter@rug.nl`.

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
