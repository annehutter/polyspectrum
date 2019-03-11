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

#include "confObj.h"
#include "grid.h"

#include "fft.h"
#include "utils_fftw.h"
#include "filter.h"
#include "polyspectrum.h"

fftw_complex *generate_values_polygons(int nbins, int local_n0, fftw_complex *kfilter, fftw_complex *fft_array)
{
    fftw_complex *thisArray = allocate_3D_array_fftw_complex(nbins);
    fftw_complex *product = allocate_3D_array_fftw_complex(nbins);
    
    // multiply kfilter and fft_array
    product_3D_fftw_arrays(nbins, local_n0, kfilter, fft_array, product);

    // FFT kfilter to realspace
    fft_k_to_realspace(nbins, product, thisArray);
    
    fftw_free(product);
    
    return thisArray;
}


fftw_complex *generate_num_polygons(int nbins, fftw_complex *kfilter)
{
    fftw_complex *thisArray = allocate_3D_array_fftw_complex(nbins);

    // FFT kfilter to realspace
    fft_k_to_realspace(nbins, kfilter, thisArray);
    
    return thisArray;
}

double polyspectrum(int nbins, int local_n0, int local_n0_start, fftw_complex *fft_array, int n, double *k, double kbinwidth, double boxsize)
{
    fftw_complex *kfilter = NULL;
    fftw_complex *numPolygons = NULL;
    fftw_complex *valuesPolygons = NULL;
    
    fftw_complex *productNumPolygons = allocate_3D_array_fftw_complex(nbins);
    fftw_complex *productValuesPolygons = allocate_3D_array_fftw_complex(nbins);
    
    double sumNumPolygons = 0.;
    double sumValuesPolygons = 0.;
    
    double volume = boxsize * boxsize * boxsize;
    double result = 0.;

    for(int i=0; i<local_n0*nbins*nbins; i++)
    {
        productNumPolygons[i] = 1.+0.*I;
        productValuesPolygons[i] = 1.+0.*I;
    }
    
    for(int i=0; i<n; i++)
    {
        /* ------------------------------------ */
        /* Generating delta(n, k_i) & I(n, k_i) */
        /* ------------------------------------ */

        kfilter = generate_kfilter(nbins, local_n0, local_n0_start, k[i], kbinwidth, boxsize);
        
        /* I(n, k_i) */
        numPolygons = generate_num_polygons(nbins, kfilter);
        
        /* delta(n, k_i) */
        valuesPolygons = generate_values_polygons(nbins, local_n0, kfilter, fft_array);

        /* product */
        product_3D_fftw_arrays(nbins, local_n0, productNumPolygons, numPolygons, productNumPolygons);
        product_3D_fftw_arrays(nbins, local_n0, productValuesPolygons, valuesPolygons, productValuesPolygons);

        fftw_free(kfilter);
        fftw_free(numPolygons);
        fftw_free(valuesPolygons);
    }

    /* ---------------------------------------- */
    /* Summing delta(n, k_i) & I(n, k_i) over n */
    /* ---------------------------------------- */
    
    sumNumPolygons = creal(sum_3D_fftw_array(nbins, local_n0, productNumPolygons));
    sumValuesPolygons = creal(sum_3D_fftw_array(nbins, local_n0, productValuesPolygons));
    
#ifdef __MPI
    double recvSumNumPolygons = 0.;
    double recvSumValuesPolygons = 0.;

    MPI_Allreduce(&sumNumPolygons, &recvSumNumPolygons, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&sumValuesPolygons, &recvSumValuesPolygons, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    sumNumPolygons = recvSumNumPolygons;
    sumValuesPolygons = recvSumValuesPolygons;
#endif
    
    result = pow(volume, n-1) * pow(nbins, -n*3) * creal(sumValuesPolygons) / creal(sumNumPolygons);
        
    fftw_free(productNumPolygons);
    fftw_free(productValuesPolygons);
    
    return result;
}

double num_polygons(int nbins, int local_n0, int local_n0_start, int n, double *k, double kbinwidth, double boxsize)
{
    fftw_complex *kfilter = NULL;
    fftw_complex *numPolygons = NULL;
    
    fftw_complex *productNumPolygons = allocate_3D_array_fftw_complex(nbins);
    
    double sumNumPolygons = 0.;
    
    double result = 0.;

    for(int i=0; i<local_n0*nbins*nbins; i++)
    {
        productNumPolygons[i] = 1.+0.*I;
    }
    
    for(int i=0; i<n; i++)
    {
        /* ------------------------------------ */
        /* Generating I(n, k_i)                 */
        /* ------------------------------------ */

        kfilter = generate_kfilter(nbins, local_n0, local_n0_start, k[i], kbinwidth, boxsize);
        
        /* I(n, k_i) */
        numPolygons = generate_num_polygons(nbins, kfilter);
        
        /* product */
        product_3D_fftw_arrays(nbins, local_n0, productNumPolygons, numPolygons, productNumPolygons);

        fftw_free(kfilter);
        fftw_free(numPolygons);
    }

    /* ---------------------------------------- */
    /* Summing I(n, k_i) over n                 */
    /* ---------------------------------------- */
    
    sumNumPolygons = creal(sum_3D_fftw_array(nbins, local_n0, productNumPolygons));
    
#ifdef __MPI
    double recvSumNumPolygons = 0.;

    MPI_Allreduce(&sumNumPolygons, &recvSumNumPolygons, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    sumNumPolygons = recvSumNumPolygons;
#endif
    
    result = creal(sumNumPolygons);
        
    fftw_free(productNumPolygons);
    
    return result;
}
