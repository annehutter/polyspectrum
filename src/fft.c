#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "fft.h"

/* THE BELOW DEFINITIONS ARE PURELY FFTW, FOWARD IS IN AGREEMENT WITH NUMPY'S FFTN, WHILE BACKWARD NEEDS TO BE MULTIPLIED WITH 1/N^3 */

void fft_real_to_kspace(int nbins, fftw_complex *input, fftw_complex *output)
{
#ifdef __MPI
    ptrdiff_t local_n0, local_0_start;
#else
    ptrdiff_t local_n0;
#endif
    fftw_plan plan_input;
        
#ifdef __MPI    
    fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
    plan_input = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, input, output, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE); 
#else 
    plan_input = fftw_plan_dft_3d(nbins, nbins, nbins, input, output, FFTW_FORWARD, FFTW_ESTIMATE);
#endif
    
    fftw_execute(plan_input);
    fftw_destroy_plan(plan_input);
}

void fft_k_to_realspace(int nbins, fftw_complex *input, fftw_complex *output)
{
#ifdef __MPI
    ptrdiff_t local_n0, local_0_start;
#else
    ptrdiff_t local_n0;
#endif
    fftw_plan plan_input;
        
#ifdef __MPI    
    fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);    
    plan_input = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, input, output, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE); 
#else 
    plan_input = fftw_plan_dft_3d(nbins, nbins, nbins, input, output, FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
    
    fftw_execute(plan_input);
    fftw_destroy_plan(plan_input);
}
