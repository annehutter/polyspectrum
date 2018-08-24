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

void run(confObj_t simParam)
{
    kvectors_t *theseKvectors = read_params_to_kvectors(simParam);
    grid_t *thisGrid = initGrid_with_values(simParam->grid_size);

    /* READING NECESSARY GRIDS & DERIVE THEIR FOURIER TRANSFORMATION */
    fftw_complex *thisFTfield = get_FT_field(thisGrid, simParam);

    for(int i=0; i<theseKvectors->numValues; i++)
    {
        theseKvectors->kpolygon[theseKvectors->n-1] = theseKvectors->k[i];
        if(theseKvectors->n == 2) 
        {
            theseKvectors->kpolygon[0] = theseKvectors->k[i];
            printf("k = %e\t P(k) = %e\n", theseKvectors->k[i], polyspectrum(thisGrid->nbins, thisGrid->local_n0, thisGrid->local_0_start, thisFTfield, theseKvectors->n, theseKvectors->kpolygon, thisGrid->box_size));
        }
        if(theseKvectors->n == 3)
            printf("k1 = %e\t k2 = %e\n", theseKvectors->kpolygon[0], theseKvectors->kpolygon[1]);
            printf("k3 = %e\t B(k) = %e\n", theseKvectors->k[i], polyspectrum(thisGrid->nbins, thisGrid->local_n0, thisGrid->local_0_start, thisFTfield, theseKvectors->n, theseKvectors->kpolygon, thisGrid->box_size));
    }
    
    fftw_free(thisFTfield);
    
    deallocate_grid(thisGrid);
    deallocate_kvectors(theseKvectors);
}
