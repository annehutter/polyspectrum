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
#include "utils_fftw.h"
#include "fft.h"
#include "input_field.h"

void generate_XHIdens_field(grid_t *thisGrid)
{
    int nbins = thisGrid->nbins;
    int local_n0 = thisGrid->local_n0;
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int l=0; l<nbins; l++)
            {
                thisGrid->signal21cm[i*nbins*nbins+j*nbins+l] = creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+l]) * (1. - creal(i*nbins*nbins+j*nbins+l));
            }
        }
    }
}


fftw_complex *get_FT_field(grid_t *thisGrid, confObj_t simParam)
{
    fftw_complex *thisField = NULL;

    //read in fields (option 1: ionization field, option 2: density & ionization field ->21cm field)
    read_boxsize(thisGrid, simParam->box_size);
    printf("%s++\n", simParam->which_field);
    
    if(strcmp(simParam->which_field, "DENS") == 0)
    {
        read_array(thisGrid->igm_density, thisGrid, simParam->density_file, simParam->gas_inputs_in_dp);
        thisField = thisGrid->igm_density;
    }
    else if(strcmp(simParam->which_field, "XHII") == 0)
    {
        read_array(thisGrid->XHII, thisGrid, simParam->ion_file, simParam->ion_inputs_in_dp);
        thisField = thisGrid->XHII;
    }
    else if(strcmp(simParam->which_field, "XHI_DENS") == 0)
    {
        read_array(thisGrid->igm_density, thisGrid, simParam->density_file, simParam->gas_inputs_in_dp);
        read_array(thisGrid->XHII, thisGrid, simParam->ion_file, simParam->ion_inputs_in_dp);
        generate_XHIdens_field(thisGrid);
        thisField = thisGrid->signal21cm;
    }
    else
    {
        fprintf(stderr, "Not supported option. Please use 'DENS', 'XHII' or 'XHI_DENS'\n");
        exit(EXIT_FAILURE);
    }
    
    if(thisGrid->local_0_start == 0) printf("Chosen field is %s. Mean value from FT is %e\n", simParam->which_field, creal(thisField[0]));
    
    /* allocating Fourier transformed grid */
    fftw_complex *thisFTfield = allocate_3D_array_fftw_complex(thisGrid->nbins);
        
    //Fourier transforming field to k-space
    fft_real_to_kspace(thisGrid->nbins, thisField, thisFTfield);
    printf("done FFT\n");
    
    return thisFTfield;
}
