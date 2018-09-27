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

#include "utils.h"
#include "confObj.h"
#include "grid.h"

#include "fft.h"
#include "utils_fftw.h"
#include "filter.h"
#include "kvectors.h"
#include "polyspectrum.h"

#include "input_field.h"
#include "run.h"

void run(confObj_t simParam)
{
    kvectors_t *theseKvectors = read_params_to_kvectors(simParam);
    grid_t *thisGrid = initGrid_with_values(simParam->grid_size);

    /* READING NECESSARY GRIDS & DERIVE THEIR FOURIER TRANSFORMATION */
    fftw_complex *thisFTfield = get_FT_field(thisGrid, simParam);

    double *polyspec = allocate_array_double(theseKvectors->numValues, "polyspectrum");
    
    for(int i=0; i<theseKvectors->numValues; i++)
    {
        theseKvectors->kpolygon[theseKvectors->n-1] = theseKvectors->k[i];
        if(theseKvectors->n == 2) theseKvectors->kpolygon[0] = theseKvectors->k[i];
        if(theseKvectors->n == 3)
        {
            if(theseKvectors->kpolygon[0] == theseKvectors->kpolygon[1])
            {
                for(int j=0; j<3; j++) theseKvectors->kpolygon[j] = theseKvectors->k[i];
            }
        }
        
        polyspec[i] = polyspectrum(thisGrid->nbins, thisGrid->local_n0, thisGrid->local_0_start, thisFTfield, theseKvectors->n, theseKvectors->kpolygon, thisGrid->box_size);
        
        if(theseKvectors->n == 2)
        {
            printf("k = %e\t P(k) = %e\n", theseKvectors->k[i], polyspec[i]);
        }
        if(theseKvectors->n == 3)
        {
            printf("k1 = %e\t k2 = %e\n", theseKvectors->kpolygon[0], theseKvectors->kpolygon[1]);
            printf("k3 = %e\t B(k) = %e\n", theseKvectors->k[i], polyspec[i]);
        }
    }
    
    // only rank 0 should output
    if(thisGrid->local_0_start == 0) 
        save_polyspectrum(simParam, theseKvectors->numValues, theseKvectors->theta, theseKvectors->k, polyspec);
    
    free(polyspec);
    fftw_free(thisFTfield);
    
    deallocate_grid(thisGrid);
    deallocate_kvectors(theseKvectors);
}


void save_polyspectrum(confObj_t simParam, int num, double *theta, double *k, double *polyspectrum)
{
    FILE *f;
    char ending[MAXLENGTH];
    char *filename = NULL;
    
    if(simParam->n == 2)
        sprintf(ending, "_ps.dat");
    else if (simParam->n == 3)
    {
        if(simParam->num_values == 1)
            sprintf(ending, "_bs_k_%4.2e_%4.2e_theta_%4.2e.dat", simParam->k1, simParam->k2, simParam->theta);
        else sprintf(ending, "_bs_k_%4.2e_%4.2e.dat", simParam->k1, simParam->k2);
    }

    filename = concat3(simParam->output_dir, "/", simParam->output_basename, ending);

    f = fopen(filename, "wb");

    if(simParam->n == 2 || theta == NULL)
    {
        fprintf(f, "# k [h^-1 Mpc]\t Polyspectrum_%d(k)\n", simParam->n);
        for(int i=0; i<num; i++)
        {
            fprintf(f, "%e\t%e\n", k[i], polyspectrum[i]);
        }
    }
    else if(simParam->n ==3)
    {
        fprintf(f, "# theta [rad]\t k [h^-1 Mpc]\t Polyspectrum_%d(k)\n", simParam->n);
        for(int i=0; i<num; i++)
        {
            fprintf(f, "%e\t%e\t%e\n", theta[i], k[i], polyspectrum[i]);
        } 
    }
    
    fclose(f);
    
    free(filename);
}
