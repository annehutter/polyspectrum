#ifndef GRID_H
#define GRID_H

/* structure for grid */

typedef struct
{
    int nbins;
    float box_size;
    
    float xmin, ymin, zmin;
    
    fftw_complex *igm_density;
    
    // hydrogen
    fftw_complex *XHII;    
    fftw_complex *signal21cm;
    
    //domain decomposition
    int local_n0;
    int local_0_start;
} grid_t;

/* functions */

grid_t *initGrid();
grid_t *initGrid_with_values(int nbins);
void read_boxsize(grid_t *thisGrid, double box_size);
void read_array(fftw_complex *toThisArray, grid_t *thisGrid, char *filename, int double_precision);

#ifdef __MPI
void read_grid(fftw_complex *toThisArray, int nbins, int local_n0, int local_0_start, char *filename);
void read_grid_doubleprecision(fftw_complex *toThisArray, int nbins, int local_n0, int local_0_start, char *filename);
#else
void read_grid(fftw_complex *toThisArray, int nbins, int local_n0, char *filename);
void read_grid_doubleprecision(fftw_complex *toThisArray, int nbins, int local_n0, char *filename);
#endif

void initialize_grid(fftw_complex *thisArray, int nbins, int local_n0, double value);
void deallocate_grid(grid_t *thisGrid);

#ifdef __MPI
void write_grid_to_file_float(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename);
void write_grid_to_file_double(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename);
#else
void write_grid_to_file_float(fftw_complex *thisArray, int nbins, int local_n0, char *filename);
void write_grid_to_file_double(fftw_complex *thisArray, int nbins, int local_n0, char *filename);
#endif

void save_to_file(fftw_complex *thisArray, grid_t *thisGrid, char *filename);

double get_mean_grid(fftw_complex *thisArray, int nbins, int local_n0);

#endif
