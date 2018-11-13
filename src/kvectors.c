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
#include "kvectors.h"

#define PI acos(-1.0)

kvectors_t *initKvectors()
{
    kvectors_t *newKvectors = malloc(sizeof(kvectors_t));
    if(newKvectors == NULL)
    {
        fprintf(stderr, "ERROR: initKvectors: Not enough memory to allocate kvectors.\n");
        exit(EXIT_FAILURE);
    }
    
    newKvectors->n = 0;
    newKvectors->kpolygon = NULL;
    
    newKvectors->numValues = 0;
    newKvectors->theta = NULL;
    newKvectors->k = NULL;
    
    return newKvectors;
}

kvectors_t *read_params_to_kvectors(confObj_t simParam)
{
    kvectors_t *theseKvectors = initKvectors();
    
    int n = simParam->n;
    double k1 = simParam->k1;
    double k2 = simParam->k2;
    double theta = simParam->theta;
    int numValues = simParam->num_values;
    int grid_size = simParam->grid_size;
    double box_size = simParam->box_size;
    
    theseKvectors->n = n;
    theseKvectors->kpolygon = allocate_array_double(n, "kpolygon");
    theseKvectors->numValues = numValues;

    if(n<2)
    {
        fprintf(stderr, "No valid input for n being the number of k-vectors, as their sum is unequal to zero: n>=2\n");
        exit(EXIT_FAILURE);
    }
    else if(n == 2)
    {
        if(simParam->num_values >1)
        {
            theseKvectors->numValues = grid_size/2;
            theseKvectors->k = generate_k_values_powerspectrum(grid_size, box_size);
        }
        else
        {
            assert(theseKvectors->numValues == 1);
            numValues = 1;
            theseKvectors->k = allocate_array_double(numValues, "k");
            
            theseKvectors->k[0] = k1;
            theseKvectors->k[1] = k1;
        }
        
        for(int i=0; i<n; i++) theseKvectors->kpolygon[i] = theseKvectors->k[i];
    }
    else if(n == 3)
    {
        if(simParam->equilateral == 1)
        {
            theseKvectors->kpolygon[0] = k1;
            theseKvectors->kpolygon[1] = k2;
            
            assert(theseKvectors->numValues > 0);
            theseKvectors->k = generate_k_values_num(numValues, grid_size, box_size);
        }
        else
        {   
            theseKvectors->kpolygon[0] = k1;
            theseKvectors->kpolygon[1] = k2;
            if(simParam->num_values > 1)
            {
                /* multiple values */
                theseKvectors->theta = generate_theta_values(numValues);
                theseKvectors->k = generate_k_values_bispectrum(numValues, k1, k2);
            }
            else
            {
                /* single values */
                assert(theseKvectors->numValues == 1);
                numValues = 1;
                theseKvectors->theta = allocate_array_double(numValues, "theta");
                theseKvectors->k = allocate_array_double(numValues, "k");
                
                theseKvectors->theta[0] = theta;
                theseKvectors->k[0] = k1*k1 + k2*k2 - 2.*k1*k2*cos(theta);
            }
        }
    }
    else
    {
        printf("Not implemented yet\n");
        exit(EXIT_FAILURE);
    }
    
    return theseKvectors;
}

void deallocate_kvectors(kvectors_t *theseKvectors)
{
    if(theseKvectors->kpolygon != NULL) free(theseKvectors->kpolygon);
    if(theseKvectors->theta != NULL) free(theseKvectors->theta);
    if(theseKvectors->k != NULL) free(theseKvectors->k);
    
    free(theseKvectors);
}

double *generate_cosTheta_values(int numValues)
{
    double *cosTheta = allocate_array_double(numValues, "theta");

    for(int i=0; i<numValues; i++)
    {
        cosTheta[i] = 1. - 2. * (double)i / (double)numValues; /* from 1 to -1 */ 
    }

    return cosTheta;
}

double *generate_theta_values(int numValues)
{
    double *theta = allocate_array_double(numValues, "theta");

    for(int i=0; i<numValues; i++)
    {
        theta[i] = PI * (double)i / (double)numValues; /* from 1 to -1 */ 
    }

    return theta;
}

double *generate_k_values_bispectrum(int numValues, double k1, double k2)
{
    double *k = allocate_array_double(numValues, "k");
    double *theta = generate_theta_values(numValues);
    
    for(int i=0; i<numValues; i++)
    {
        k[i] = sqrt(k1*k1 + k2*k2 - 2.*k1*k2*cos(theta[i]));
        printf("%d: theta = %e\t cos(theta) = %e\t k = %e\n", i, theta[i], cos(theta[i]), k[i]);
    }
    
    free(theta);
    
    return k;
}

double *generate_k_values_powerspectrum(int grid_size, double box_size)
{
    int mid = grid_size/2;
    double *k = allocate_array_double(mid, "k");
    
    for(int i=0; i<mid; i++)
    {
        k[i] = 2.*PI/box_size * (double)(i+1);
    }
    
    return k;
}

double *generate_k_values_num(int numValues, int grid_size, double box_size)
{
    double *k = allocate_array_double(numValues, "k");
    double kL = 2.*PI/box_size;
    
    for(int i=0; i<numValues; i++)
    {
        k[i] = kL * pow(grid_size/2,(double)(i+1)/(double)numValues);
    }
    
    return k;
}
