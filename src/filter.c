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
#include "filter.h"

#define PI acos(-1.0)


fftw_complex *generate_kfilter(int nbins, int local_n0, int local_n0_start, double k, double boxsize)
{
    fftw_complex *kfilter = allocate_3D_array_fftw_complex(nbins);
    
    construct_kfilter(nbins, local_n0, local_n0_start, kfilter, k, boxsize);
    
    return kfilter;
}

void construct_kfilter(int nbins, int local_n0, int local_n0_start, fftw_complex *array, double k, double boxsize)
{
    double k_in_kf = k / (2.*PI) * boxsize;        // in units of k_F = 2*pi/L
    int half_nbins = nbins/2;
    int mx = 0, my = 0, mz = 0;
    double m = 0.;
    double expr = 0.;
    
    for(int i=0; i<local_n0; i++)
    {
        if(i+local_n0_start<half_nbins) mx = i+local_n0_start;
        else mx = -nbins + i + local_n0_start;
        
        for(int j=0; j<nbins; j++)
        {
            if(j<half_nbins) my = j;
            else my = -nbins + j;
            
            for(int l=0; l<nbins; l++)
            {
                if(l<half_nbins) mz = l;
                else mz = -nbins + l;
                
                m = sqrt((double)(mx*mx + my*my + mz*mz));
                expr = fabs(k_in_kf - m);
                array[i*nbins*nbins + j*nbins + l] = (expr<=0.5) ? 1.+0.*I : 0.+0.*I;
            }
        }
    }
}
