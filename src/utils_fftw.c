#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "utils_fftw.h"

fftw_complex *allocate_3D_array_fftw_complex(int nbins)
{
    fftw_complex *newArray = NULL;
#ifdef __MPI    
    ptrdiff_t alloc_local, local_n0, local_0_start;
    alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
    newArray = fftw_alloc_complex(alloc_local);
#else
    newArray = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
#endif
    
    return newArray;
}

void product_3D_fftw_arrays(int nbins, int local_n0, fftw_complex *array1, fftw_complex *array2, fftw_complex *output)
{
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int l=0; l<nbins; l++)
            {
                output[i*nbins*nbins + j*nbins + l] = array1[i*nbins*nbins + j*nbins + l] * array2[i*nbins*nbins + j*nbins + l];
            }
        }
    }
}

void product3_3D_fftw_arrays(int nbins, int local_n0, fftw_complex *array1, fftw_complex *array2, fftw_complex *array3, fftw_complex *output)
{
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int l=0; l<nbins; l++)
            {
                output[i*nbins*nbins + j*nbins + l] = array1[i*nbins*nbins + j*nbins + l] * array2[i*nbins*nbins + j*nbins + l] * array3[i*nbins*nbins + j*nbins +l];
            }
        }
    }
}

fftw_complex sum_3D_fftw_array(int nbins, int local_n0, fftw_complex *array)
{
    fftw_complex sum = 0.+0.*I;
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int l=0; l<nbins; l++)
            {
                sum += array[i*nbins*nbins + j*nbins + l];
            }
        }
    }
    
    return sum;
}

void calc_phase_3D_fftw_array(int nbins, int local_n0, fftw_complex *array)
{    
    fftw_complex value = 0. + 0.*I;
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int l=0; l<nbins; l++)
            {
                value = array[i*nbins*nbins + j*nbins + l];
                array[i*nbins*nbins + j*nbins + l] = creal(value)/cabs(value);
            }
        }
    }
}