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

#define PI acos(-1.0)


int main (int argc, /*const*/ char * argv[]) { 
    int size = 1;
    int myRank = 0;

    char iniFile[MAXLENGTH];
    confObj_t simParam;
    kvectors_t *kVectors;
    
    double t1, t2;

#ifdef __MPI
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); 
    
    t1 = MPI_Wtime();
    
    fftw_mpi_init();
#else
    t1 = time(NULL);
#endif

    //parse command line arguments and be nice to user
    if (argc != 2) {
        printf("cifog: (C)  - Use at own risk...\n");
        printf("USAGE:\n");
        printf("cifog iniFile\n");
        
        exit(EXIT_FAILURE);
    } else {
        strcpy(iniFile, argv[1]);
    }
        
    //read parameter file
    simParam = readConfObj(iniFile);
    kVectors = read_params_to_kvectors(simParam);
    
    //grid allocation
    grid_t *grid = initGrid();
    read_files_to_grid(grid, simParam);
    
    //read in fields (option 1: ionization field, option 2: density & ionization field ->21cm field)
    read_array(grid->igm_density, grid, simParam->density_file, simParam->gas_inputs_in_dp);
    if(myRank == 0) printf("mean = %e\n", creal(grid->igm_density[0]));
    fftw_complex *output = allocate_3D_array_fftw_complex(grid->nbins);
        
    //FFT field to k-space
    fft_real_to_kspace(grid->nbins, grid->igm_density, output);
    save_to_file(grid->igm_density, grid, "test_density.dat");
    printf("done FFT\n");

    for(int i=0; i<kVectors->numValues; i++)
    {
        kVectors->kpolygon[kVectors->n-1] = kVectors->k[i];
        if(kVectors->n == 2) 
        {
            kVectors->kpolygon[0] = kVectors->k[i];
            printf("k = %e\t P(k) = %e\n", kVectors->k[i], polyspectrum(grid->nbins, grid->local_n0, grid->local_0_start, output, kVectors->n, kVectors->kpolygon, grid->box_size));
        }
        if(kVectors->n == 3)
            printf("k3 = %e\t B(k) = %e\n", kVectors->k[i], polyspectrum(grid->nbins, grid->local_n0, grid->local_0_start, output, kVectors->n, kVectors->kpolygon, grid->box_size));
    }
    
    //deallocation
    fftw_free(output);    
    deallocate_grid(grid);
    deallocate_kvectors(kVectors);
    confObj_del(&simParam);

#ifdef __MPI
    fftw_mpi_cleanup();
        
    t2 = MPI_Wtime();
    printf("Execution took %f s\n", t2-t1);
    MPI_Finalize();
#else
    fftw_cleanup();
    
    t2 = time(NULL);
    printf("Execution took %f s\n", t2-t1);
#endif
    
    return 0;
}

