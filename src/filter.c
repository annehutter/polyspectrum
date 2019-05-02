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

/*------------------------------------------------------------------------*/
/* FILTER WITH MINIMAL KBINWIDTH */
/*------------------------------------------------------------------------*/

fftw_complex *generate_kfilter(int nbins, int local_n0, int local_n0_start, double k, double binwidth, double boxsize)
{
    fftw_complex *kfilter = allocate_3D_array_fftw_complex(nbins);
    
    construct_kfilter(nbins, local_n0, local_n0_start, kfilter, k, binwidth, boxsize);
    
    return kfilter;
}

void construct_kfilter(int nbins, int local_n0, int local_n0_start, fftw_complex *array, double k, double binwidth, double boxsize)
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
                array[i*nbins*nbins + j*nbins + l] = (expr<=0.5*binwidth) ? 1.+0.*I : 0.+0.*I;
            }
        }
    }
}

/*------------------------------------------------------------------------*/
/* FILTER WITH DETERMINING MINIMUM & MAXIMUM KBINWIDTH */
/*------------------------------------------------------------------------*/

fftw_complex *generate_kfilter_with_determining_kbinwidth(int nbins, int local_n0, int local_n0_start, double *kmin, double *kmax, double k, double binwidth, double boxsize)
{
    fftw_complex *kfilter = allocate_3D_array_fftw_complex(nbins);
    
    construct_kfilter_with_determining_kbinwidth(nbins, local_n0, local_n0_start, kfilter, kmin, kmax, k, binwidth, boxsize);
    
    return kfilter;
}

void construct_kfilter_with_determining_kbinwidth(int nbins, int local_n0, int local_n0_start, fftw_complex *array, double *kmin, double *kmax, double k, double binwidth, double boxsize)
{
    double k_in_kf = k / (2.*PI) * boxsize;        // in units of k_F = 2*pi/L
    int half_nbins = nbins/2;
    int mx = 0, my = 0, mz = 0;
    double m = 0.;
    double expr = 0.;
    double mmin = 1.e30;
    double mmax = 0.;
    
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
                array[i*nbins*nbins + j*nbins + l] = (expr<=0.5*binwidth) ? 1.+0.*I : 0.+0.*I;
                if(expr<=0.5*binwidth)
                    get_corners(mx, my, mz, binwidth, &mmin, &mmax);
            }
        }
    }
    
    *kmin = mmin * 2.*PI / boxsize;
    *kmax = mmax * 2.*PI / boxsize;
}


/*------------------------------------------------------------------------*/
/* FILTER WITH IMPOSING MINIMUM & MAXIMUM K-VALUES */
/*------------------------------------------------------------------------*/

fftw_complex *generate_kfilter_with_kbinwidth(int nbins, int local_n0, int local_n0_start, double kmin, double kmax, double binwidth, double boxsize)
{
    fftw_complex *kfilter = allocate_3D_array_fftw_complex(nbins);
    
    construct_kfilter_with_kbinwidth(nbins, local_n0, local_n0_start, kfilter, kmin, kmax, binwidth, boxsize);
    
    return kfilter;
}

void construct_kfilter_with_kbinwidth(int nbins, int local_n0, int local_n0_start, fftw_complex *array, double kmin, double kmax, double binwidth, double boxsize)
{
    int half_nbins = nbins/2;
    int mx = 0, my = 0, mz = 0;
    double m = 0.;
    double expr_min = 0.;
    double expr_max = 0.;
    double mmin = kmin / (2.*PI) * boxsize;        // in units of k_F = 2*pi/L
    double mmax = kmax / (2.*PI) * boxsize;        // in units of k_F = 2*pi/L
    
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
                expr_min = mmin - m;
                expr_max = m - mmax;
                array[i*nbins*nbins + j*nbins + l] = (expr_min<=0.5*binwidth && expr_max<=0.5*binwidth) ? 1.+0.*I : 0.+0.*I;
            }
        }
    }
}


/*------------------------------------------------------------------------*/
/* FILTERS AS USED IN WATKINSON ET AL. 2017, 2018 */
/*------------------------------------------------------------------------*/

fftw_complex *generate_kfilter_Watkinson(int nbins, int local_n0, int local_n0_start, double k, double binwidth, double boxsize)
{
    fftw_complex *kfilter = allocate_3D_array_fftw_complex(nbins);
    
    construct_kfilter_Watkinson(nbins, local_n0, local_n0_start, kfilter, k, binwidth, boxsize);
    
    return kfilter;
}

void construct_kfilter_Watkinson(int nbins, int local_n0, int local_n0_start, fftw_complex *array, double k, double binwidth, double boxsize)
{
    double k_in_kf = k / (2.*PI) * boxsize;        // in units of k_F = 2*pi/L
    int half_nbins = nbins/2;
    int mx = 0, my = 0, mz = 0;
    double pmin = 1.e30;
    double pmax = 0.;
    
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
                
                pmin = 1.e30;
                pmax = 0.;
                get_corners(mx, my, mz, binwidth, &pmin, &pmax);
                array[i*nbins*nbins + j*nbins + l] = (pmin<=k_in_kf && k_in_kf<=pmax) ? 1.+0.*I : 0.+0.*I;                   
            }
        }
    }
}

fftw_complex *generate_kfilter_Watkinson_kn(int nbins, int local_n0, int local_n0_start, double kmin, double kmax, double binwidth, double boxsize)
{
    fftw_complex *kfilter = allocate_3D_array_fftw_complex(nbins);
    
    construct_kfilter_Watkinson_kn(nbins, local_n0, local_n0_start, kfilter, kmin, kmax, binwidth, boxsize);
    
    return kfilter;
}

void construct_kfilter_Watkinson_kn(int nbins, int local_n0, int local_n0_start, fftw_complex *array, double kmin, double kmax, double binwidth, double boxsize)
{
    int half_nbins = nbins/2;
    int mx = 0, my = 0, mz = 0;
    double mmin = kmin / (2.*PI) * boxsize;        // in units of k_F = 2*pi/L
    double mmax = kmax / (2.*PI) * boxsize;        // in units of k_F = 2*pi/L
    double pmin = 1.e30;
    double pmax = 0.;
    
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
                
                pmin = 1.e30;
                pmax = 0.;
                get_corners(mx, my, mz, binwidth, &pmin, &pmax);
                array[i*nbins*nbins + j*nbins + l] = ((pmax>=mmax && pmin<=mmax) || (pmax<=mmax && pmin>=mmin) || (pmax>=mmin && pmin<=mmin)) ? 1.+0.*I : 0.+0.*I;                   
            }
        }
    }
}


/*------------------------------------------------------------------------*/
/* AUXILIARY FUNCTIONS */
/*------------------------------------------------------------------------*/

void get_corners(int mx, int my, int mz, double binwidth, double *mmin, double *mmax)
{
    double x = (double)mx;
    double y = (double)my;
    double z = (double)mz;
    double effBinwidth = 0.5*binwidth;
    
    int numCorners = 8;
    double corner[numCorners];
    
    double minValue = (*mmin) * (*mmin);
    double maxValue = (*mmax) * (*mmax);
    
    corner[0] = (x-effBinwidth)*(x-effBinwidth) + (y-effBinwidth)*(y-effBinwidth) + (z-effBinwidth)*(z-effBinwidth);
    corner[1] = (x+effBinwidth)*(x+effBinwidth) + (y-effBinwidth)*(y-effBinwidth) + (z-effBinwidth)*(z-effBinwidth);
    corner[2] = (x-effBinwidth)*(x-effBinwidth) + (y+effBinwidth)*(y+effBinwidth) + (z-effBinwidth)*(z-effBinwidth);
    corner[3] = (x-effBinwidth)*(x-effBinwidth) + (y-effBinwidth)*(y-effBinwidth) + (z+effBinwidth)*(z+effBinwidth);
    corner[4] = (x+effBinwidth)*(x+effBinwidth) + (y+effBinwidth)*(y+effBinwidth) + (z-effBinwidth)*(z-effBinwidth);
    corner[5] = (x+effBinwidth)*(x+effBinwidth) + (y-effBinwidth)*(y-effBinwidth) + (z+effBinwidth)*(z+effBinwidth);
    corner[6] = (x-effBinwidth)*(x-effBinwidth) + (y+effBinwidth)*(y+effBinwidth) + (z+effBinwidth)*(z+effBinwidth);
    corner[7] = (x+effBinwidth)*(x+effBinwidth) + (y+effBinwidth)*(y+effBinwidth) + (z+effBinwidth)*(z+effBinwidth);

    for(int i=0; i<numCorners; i++)
    {
        if(corner[i] > maxValue)
            maxValue = corner[i];
        if(corner[i] < minValue)
            minValue = corner[i];
    }
    
    *mmin = sqrt(minValue);
    *mmax = sqrt(maxValue);
}

void calc_k3_min_and_max(double *k, double **kmin, double **kmax)
{
    /* get cos(theta) value*/
    double cosTheta = (k[0]*k[0] + k[1]*k[1] - k[2]*k[2]) / (2.*k[0]*k[1]);
    double *kminLoc = *kmin;
    double *kmaxLoc = *kmax; 
    
    double temp[8];
    double tempMin = 1.e30;
    double tempMax = 0.;
    temp[0] = (*kmin)[0]*(*kmin)[0] + (*kmin)[1]*(*kmin)[1] - 2.*(*kmin)[0]*(*kmin)[1]*cosTheta;
    temp[1] = (*kmin)[0]*(*kmin)[0] + (*kmax)[1]*(*kmax)[1] - 2.*(*kmin)[0]*(*kmax)[1]*cosTheta;
    temp[2] = (*kmax)[0]*(*kmax)[0] + (*kmin)[1]*(*kmin)[1] - 2.*(*kmax)[0]*(*kmin)[1]*cosTheta;
    temp[3] = (*kmax)[0]*(*kmax)[0] + (*kmax)[1]*(*kmax)[1] - 2.*(*kmax)[0]*(*kmax)[1]*cosTheta;
    temp[4] = k[0]*k[0] + (*kmax)[1]*(*kmax)[1] - 2.*k[0]*(*kmax)[1]*cosTheta;
    temp[5] = (*kmax)[0]*(*kmax)[0] + k[1]*k[1] - 2.*(*kmax)[0]*k[1]*cosTheta;
    temp[6] = k[0]*k[0] + (*kmin)[1]*(*kmin)[1] - 2.*k[0]*(*kmin)[1]*cosTheta;
    temp[7] = (*kmin)[0]*(*kmin)[0] + k[1]*k[1] - 2.*(*kmin)[0]*k[1]*cosTheta;
    
    for(int i=0; i<8; i++)
    {
        if(temp[i] < tempMin)
            tempMin = temp[i];
        if(temp[i] > tempMax)
            tempMax = temp[i];
    }
    kminLoc[2] = sqrt(tempMin);
    kmaxLoc[2] = sqrt(tempMax);
}

void calc_kn_min_and_max(int n, double *k, double **kmin, double **kmax)
{
    (*kmin[n-1]) = k[n-1];
    (*kmax[n-1]) = k[n-1];
}