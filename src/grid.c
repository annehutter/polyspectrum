#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "utils.h"
#include "grid.h"

/* functions for grid -----------------------------------------------------------------------*/

grid_t *initGrid()
{
    grid_t *newGrid;
    newGrid = malloc(sizeof(grid_t));
    if(newGrid == NULL)
    {
        fprintf(stderr, "ERROR: initGrid Not enought memory to allocate Grid.\n");
                exit(EXIT_FAILURE);
    }
    
    newGrid->nbins = 0;
    newGrid->box_size =0.;

    newGrid->xmin = 0.;
    newGrid->ymin = 0.;
    newGrid->zmin = 0.;
    
    newGrid->igm_density = NULL;
    
    //hydrogen
    newGrid->XHII = NULL;    
    newGrid->signal21cm = NULL;
        
    //domain decomposition
    newGrid->local_n0 = 0;
    newGrid->local_0_start = 0;
    
    return newGrid;
}

grid_t *initGrid_with_values(int nbins)
{
    grid_t *thisGrid = initGrid();
    
#ifdef __MPI
    ptrdiff_t alloc_local, local_n0, local_0_start;
#else
    ptrdiff_t local_n0;
#endif
    
    thisGrid->nbins = nbins;    
    thisGrid->local_n0 = nbins;
    thisGrid->local_0_start = 0;
    
#ifdef __MPI    
    alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
    
    thisGrid->local_n0 = local_n0;
    thisGrid->local_0_start = local_0_start;
    
    thisGrid->igm_density = fftw_alloc_complex(alloc_local);
    thisGrid->XHII = fftw_alloc_complex(alloc_local);
    thisGrid->signal21cm = fftw_alloc_complex(alloc_local);
#else
    thisGrid->igm_density = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);    
    thisGrid->XHII = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    thisGrid->signal21cm = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
#endif

    initialize_grid(thisGrid->igm_density, nbins, local_n0, 1.);
    initialize_grid(thisGrid->XHII, nbins, local_n0, 0.);
    initialize_grid(thisGrid->signal21cm, nbins, local_n0, 0.);
    
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    return thisGrid;
}

void read_boxsize(grid_t *thisGrid, double box_size)
{
    thisGrid->box_size = box_size;
}

void read_array(fftw_complex *toThisArray, grid_t *thisGrid, char *filename, int double_precision)
{
#ifdef __MPI
    ptrdiff_t local_n0, local_0_start;
#else
    ptrdiff_t local_n0;
#endif
    int nbins;
    
    nbins = thisGrid->nbins;
    local_n0 = thisGrid->local_n0;
    
    if(double_precision == 1)
    {
#ifdef __MPI
    local_0_start = thisGrid->local_0_start;
    read_grid_doubleprecision(toThisArray, nbins, local_n0, local_0_start, filename);
#else
    read_grid_doubleprecision(toThisArray, nbins, local_n0, filename);
#endif
    }else{
#ifdef __MPI
    local_0_start = thisGrid->local_0_start;
    read_grid(toThisArray, nbins, local_n0, local_0_start, filename);
#else
    read_grid(toThisArray, nbins, local_n0, filename);
#endif
    }
}

#ifdef __MPI
void read_grid(fftw_complex *toThisArray, int nbins, int local_n0, int local_0_start, char *filename)
#else
void read_grid(fftw_complex *toThisArray, int nbins, int local_n0, char *filename)
#endif
{
    float *tmparray;
        
    tmparray = (float*)malloc(sizeof(float)*local_n0*nbins*nbins);
    if(tmparray == NULL)
    {
        fprintf(stderr, "tmparray in read_grid (grid.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
#ifdef __MPI
    int success;
    int resultlen;
    char msg[MPI_MAX_ERROR_STRING];

    MPI_File mpifile;
    MPI_Offset offset;
    MPI_Status status;
    
    offset = (local_0_start*nbins*nbins*sizeof(float));
    
    success = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,MPI_INFO_NULL, &mpifile);
    if(success != MPI_SUCCESS)
    {
        MPI_Error_string(success, msg, &resultlen);
        fprintf(stderr, "MPI_File_open(): %s\n", msg);
        exit(-1);
    }
    MPI_File_read_at_all(mpifile,offset,tmparray, local_n0*nbins*nbins,MPI_FLOAT,&status);
    MPI_File_close(&mpifile);
#else
    FILE *fp;
    fp = fopen(filename, "rb");
    fread(tmparray, sizeof(float), nbins*nbins*nbins, fp);
    fclose(fp);
#endif
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                toThisArray[i*nbins*nbins+j*nbins+k] = (double)tmparray[i*nbins*nbins+j*nbins+k]+0.*I;
            }
        }
    }
    free(tmparray);
}


#ifdef __MPI
void read_grid_doubleprecision(fftw_complex *toThisArray, int nbins, int local_n0, int local_0_start, char *filename)
#else
void read_grid_doubleprecision(fftw_complex *toThisArray, int nbins, int local_n0, char *filename)
#endif
{
    double *tmparray;
    
    tmparray = (double*)malloc(sizeof(double)*local_n0*nbins*nbins);
    if(tmparray == NULL)
    {
        fprintf(stderr, "tmparray in read_grid (grid.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
#ifdef __MPI
    int success;
    int resultlen;
    char msg[MPI_MAX_ERROR_STRING];

    MPI_File mpifile;
    MPI_Offset offset;
    MPI_Status status;
    
    offset = (local_0_start*nbins*nbins*sizeof(double));
    
    success = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,MPI_INFO_NULL, &mpifile);
    if(success != MPI_SUCCESS)
    {
        MPI_Error_string(success, msg, &resultlen);
        fprintf(stderr, "MPI_File_open(): %s\n", msg);
        exit(-1);
    }
    MPI_File_read_at_all(mpifile,offset,tmparray, local_n0*nbins*nbins,MPI_DOUBLE,&status);
    MPI_File_close(&mpifile);
#else
    FILE *fp;
    fp = fopen(filename, "rb");
    fread(tmparray, sizeof(double), nbins*nbins*nbins, fp);
    fclose(fp);
#endif
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                toThisArray[i*nbins*nbins+j*nbins+k] = tmparray[i*nbins*nbins+j*nbins+k]+0.*I;
            }
        }
    }
    free(tmparray);
}

void initialize_grid(fftw_complex *thisArray, int nbins, int local_n0, double value)
{
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                thisArray[i*nbins*nbins+j*nbins+k] = value +0.*I;
            }
        }
    }
}

void deallocate_grid(grid_t *thisGrid)
{    
    if(thisGrid->igm_density != NULL) fftw_free(thisGrid->igm_density);
    if(thisGrid->XHII != NULL) fftw_free(thisGrid->XHII);
    if(thisGrid->signal21cm != NULL) fftw_free(thisGrid->signal21cm);
    
    free(thisGrid);
}

#ifdef __MPI
void write_grid_to_file_float(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename)
#else
void write_grid_to_file_float(fftw_complex *thisArray, int nbins, int local_n0, char *filename)
#endif
{
    float *tmparray;
    
    tmparray = (float*)malloc(sizeof(float)*local_n0*nbins*nbins);
    if(tmparray == NULL)
    {
        fprintf(stderr, "tmparray in write_grid_to_file_float (grid.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                tmparray[i*nbins*nbins+j*nbins+k] = (float)creal(thisArray[i*nbins*nbins+j*nbins+k]);
            }
        }
    }
        
#ifdef __MPI
    MPI_File mpifile;
    MPI_Offset offset;
    MPI_Status status;
    
    offset = (local_0_start*nbins*nbins*sizeof(float));

    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifile);
    if (mpifile == MPI_FILE_NULL)
    {
        fprintf(stderr, "Could not open file %s\n", filename);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);     
    }
    MPI_File_write_at_all(mpifile, offset, tmparray, local_n0*nbins*nbins, MPI_FLOAT, &status);
    MPI_File_close(&mpifile);
#else
    FILE * fp;
    
    fp = fopen(filename, "wb");
    if (fp == NULL)
    {
        fprintf(stderr, "Could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    } 
    fwrite(tmparray, sizeof(float), nbins*nbins*nbins, fp);
    fclose(fp);
#endif
    
    free(tmparray);
}

#ifdef __MPI
void write_grid_to_file_double(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename)
#else
void write_grid_to_file_double(fftw_complex *thisArray, int nbins, int local_n0, char *filename)
#endif
{
    double *tmparray;
    
    tmparray = (double*)malloc(sizeof(double)*local_n0*nbins*nbins);
    if(tmparray == NULL)
    {
        fprintf(stderr, "tmparray in write_grid_to_file_double (grid.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                tmparray[i*nbins*nbins+j*nbins+k] = (double)creal(thisArray[i*nbins*nbins+j*nbins+k]);
            }
        }
    }
    
#ifdef __MPI
    MPI_File mpifile;
    MPI_Offset offset;
    MPI_Status status;
    
    offset = (local_0_start*nbins*nbins*sizeof(double));
    
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifile);
    if (mpifile == MPI_FILE_NULL)
    {
        fprintf(stderr, "Could not open file %s\n", filename);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); 
    }
    MPI_File_write_at_all(mpifile, offset, tmparray, local_n0*nbins*nbins, MPI_DOUBLE, &status);
    MPI_File_close(&mpifile);
#else
    FILE * fp;
    
    fp = fopen(filename, "wb");
    if (fp == NULL)
    {
        fprintf(stderr, "Could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    } 
    fwrite(tmparray, sizeof(double), nbins*nbins*nbins, fp);
    fclose(fp);
#endif
    
    free(tmparray);
}

void save_to_file(fftw_complex *thisArray, grid_t *thisGrid, char *filename)
{
#ifdef __MPI
    write_grid_to_file_double(thisArray, thisGrid->nbins, thisGrid->local_n0, thisGrid->local_0_start, filename);
#else
    write_grid_to_file_double(thisArray, thisGrid->nbins, thisGrid->local_n0, filename);
#endif
}

double get_mean_grid(fftw_complex *thisArray, int nbins, int local_n0)
{
    double sum = 0.;
#ifdef __MPI
    double sum_all = 0.;
#endif
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                sum += creal(thisArray[i*nbins*nbins+j*nbins+k]);
            }
        }
    }
    
#ifdef __MPI
    MPI_Allreduce(&sum, &sum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sum = sum_all / (double)(nbins*nbins*nbins);
#else
    sum = sum / (double)(nbins*nbins*nbins);
#endif
    
    return sum;
}
