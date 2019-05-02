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

void run(confObj_t simParam, int size, int thisRank)
{
    /* READING INPUT PARAMETERS */
    kvectors_t *theseKvectors = read_params_to_kvectors(simParam);
    
    /* INITIALIZING GRIDS */
    grid_t *thisGrid = initGrid_with_values(simParam->grid_size);
    if(directory_exist(simParam->output_dir) == 0)
    {
        fprintf(stderr, "Directory for output does not exist!\n");
        exit(EXIT_FAILURE);
    }

    /* READING NECESSARY GRIDS & DERIVE THEIR FOURIER TRANSFORMATION */
    fftw_complex *thisFTfield = get_FT_field(thisGrid, simParam);

    /* ALLOCATION */
    double *polyspec = allocate_array_double(theseKvectors->numValues, "polyspectrum");
    double *numPolygons = allocate_array_double(theseKvectors->numValues, "polyspectrum");

    /* COMPUTING POLYSPECTRUM */
    for(int i=0; i<theseKvectors->numValues; i++)
    {
        /* Setting kpolygon to k-values for respective polyspectrum */
        theseKvectors->kpolygon[theseKvectors->n-1] = theseKvectors->k[i];
        if(theseKvectors->n == 2) theseKvectors->kpolygon[0] = theseKvectors->k[i];
        if(theseKvectors->n == 3)
        {
            if(simParam->equilateral == 1)
            {
                for(int j=0; j<theseKvectors->n; j++) theseKvectors->kpolygon[j] = theseKvectors->k[i];
                
                /* Setting minimum and maximum k-values if required */
                if(theseKvectors->kbinningCase >= 2)
                {
                    theseKvectors->kpolygonMin[theseKvectors->n-1] = calc_k3min(theseKvectors, 0);
                    theseKvectors->kpolygonMax[theseKvectors->n-1] = calc_k3max(theseKvectors, 0);
                }
            }
            else
            {
                /* Setting minimum and maximum k-values if required */
                if(theseKvectors->kbinningCase >= 2)
                {
                    theseKvectors->kpolygonMin[theseKvectors->n-1] = calc_k3min(theseKvectors, i);
                    theseKvectors->kpolygonMax[theseKvectors->n-1] = calc_k3max(theseKvectors, i);
                }
            }
        }

        
        /* Computing polyspectrum */
        polyspec[i] = polyspectrum(thisGrid->nbins, thisGrid->local_n0, thisGrid->local_0_start, thisFTfield, theseKvectors->n, theseKvectors->kpolygon, theseKvectors->kpolygonMin, theseKvectors->kpolygonMax, simParam->kbinwidth, theseKvectors->kbinningCase, thisGrid->box_size);
        
        /* Computing number of polygons */
        if(simParam->write_numpolygons == 1)
            numPolygons[i] = num_polygons(thisGrid->nbins, thisGrid->local_n0, thisGrid->local_0_start, theseKvectors->n, theseKvectors->kpolygon, theseKvectors->kpolygonMin, theseKvectors->kpolygonMax, simParam->kbinwidth, theseKvectors->kbinningCase, thisGrid->box_size);
        
        /* some output */
        if(thisRank == 0)
        {
            if(theseKvectors->n == 2)
                printf("k = %e\t P(k) = %e\t", theseKvectors->k[i], polyspec[i]);
            
            if(theseKvectors->n == 3)
                printf("k1 = %e\t k2 = %e\tk3 = %e\t B(k) = %e\t", theseKvectors->kpolygon[0], theseKvectors->kpolygon[1], theseKvectors->k[i], polyspec[i]);
            
            if(simParam->write_numpolygons == 1)
                printf("Npolygons = %e\n", numPolygons[i]);
            
            printf("\n");
        }
    }
    
    /* WRITING TO FILE (only rank 0) */
    if(thisGrid->local_0_start == 0 && thisRank == 0)
    {
        if(simParam->write_numpolygons == 1)
            save_polyspectrum_numpolygons(simParam, theseKvectors->numValues, theseKvectors->theta, theseKvectors->k, polyspec, numPolygons);
        else
            save_polyspectrum(simParam, theseKvectors->numValues, theseKvectors->theta, theseKvectors->k, polyspec);
    }
    
    /* DEALLOCAION */
    free(polyspec);
    free(numPolygons);
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
        if(simParam->equilateral == 1)
            sprintf(ending, "_bs_equilateral.dat");
        else
        {
            if(simParam->num_values == 1)
                sprintf(ending, "_bs_k_%4.2e_%4.2e_theta_%4.2e.dat", simParam->k1, simParam->k2, simParam->theta);
            else sprintf(ending, "_bs_k_%4.2e_%4.2e.dat", simParam->k1, simParam->k2);
        }
    }

    filename = concat3(simParam->output_dir, "/", simParam->output_basename, ending);

    printf("FILENAME : %s\n", filename);
    
    f = fopen(filename, "wb");

    if(simParam->n == 2 || simParam->equilateral == 1)
    {
        fprintf(f, "# k [h^-1 Mpc]\t Polyspectrum_%d(k)\n", simParam->n);
        for(int i=0; i<num; i++)
        {
            fprintf(f, "%e\t%e\n", k[i], polyspectrum[i]);
        }
    }
    else if(simParam->n == 3)
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

void save_polyspectrum_numpolygons(confObj_t simParam, int num, double *theta, double *k, double *polyspectrum, double *numpolygons)
{
    FILE *f;
    char ending[MAXLENGTH];
    char *filename = NULL;
    
    if(simParam->n == 2)
        sprintf(ending, "_ps.dat");
    else if (simParam->n == 3)
    {
        if(simParam->equilateral == 1)
            sprintf(ending, "_bs_equilateral.dat");
        else
        {
            if(simParam->num_values == 1)
                sprintf(ending, "_bs_k_%4.2e_%4.2e_theta_%4.2e.dat", simParam->k1, simParam->k2, simParam->theta);
            else sprintf(ending, "_bs_k_%4.2e_%4.2e.dat", simParam->k1, simParam->k2);
        }
    }

    filename = concat3(simParam->output_dir, "/", simParam->output_basename, ending);

    printf("FILENAME : %s\n", filename);
    
    f = fopen(filename, "wb");

    if(simParam->n == 2 || simParam->equilateral == 1)
    {
        fprintf(f, "# k [h^-1 Mpc]\t Polyspectrum_%d(k)\n", simParam->n);
        for(int i=0; i<num; i++)
        {
            fprintf(f, "%e\t%e\t%e\n", k[i], polyspectrum[i], numpolygons[i]);
        }
    }
    else if(simParam->n == 3)
    {
        fprintf(f, "# theta [rad]\t k [h^-1 Mpc]\t Polyspectrum_%d(k)\n", simParam->n);
        for(int i=0; i<num; i++)
        {
            fprintf(f, "%e\t%e\t%e\t%e\n", theta[i], k[i], polyspectrum[i], numpolygons[i]);
        } 
    }
    
    fclose(f);
    
    free(filename);
}