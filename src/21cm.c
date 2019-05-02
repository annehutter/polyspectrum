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

#include "phys_const.h"
#include "confObj.h"
#include "grid.h"

#define PI acos(-1.)
#define CUB(X) ((X)*(X)*(X))

void generate_21cm_field(grid_t *thisGrid, confObj_t simParam)
{
    int nbins = thisGrid->nbins;
    int local_n0 = thisGrid->local_n0;
    
    const double factor = 3.*CUB(lambda_21)*A10*T21 / (8.*PI*G*mp_g) * simParam->h*1.e7/Mpc_cm * 
                        simParam->omega_b * sqrt(simParam->omega_m) * (1.-simParam->Y) * sqrt(1. + simParam->redshift);
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                thisGrid->signal21cm[i*nbins*nbins + j*nbins + k] = factor * creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k]) * (1. - creal(thisGrid->XHII[i*nbins*nbins+j*nbins+k]))+ 0.*I;
            }
        }
    }
}
